"""
Build a REAL species tree for the pilot from the NCBI Taxonomy database, using
the Entrez E-utilities HTTP API directly (pure standard library -- no Biopython,
no Selenium, no UI).

For each taxid we fetch its full lineage (efetch db=taxonomy), assemble the
shared-ancestry topology as a trie, collapse unary internal nodes, assign unit
branch lengths, and emit a Newick tree whose tips are labelled to match
data/usage.tsv. The result is a citable NCBI-taxonomy topology that
phylo.py / model_a_usage.py --tree consume directly.

Branch lengths are 1 per surviving internal edge (a cladogram); PGLS under
Brownian motion is well-defined on this (shared path length = shared lineage
depth). Swap in dated branch lengths later if desired.

Usage:
    python build_tree.py                      # -> data/tree_ncbi.nwk
    python build_tree.py --email you@x.org    # polite Entrez identification
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Tuple

from build_pilot import CURATED
from fetch_data import DATA_DIR, ensure_data_dir

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def species_in_usage(usage_path: str) -> List[str]:
    if not os.path.exists(usage_path):
        return [name for _tid, name in CURATED]
    seen = []
    with open(usage_path, "r", encoding="utf-8") as fh:
        next(fh, None)  # header
        for line in fh:
            sp = line.split("\t", 1)[0]
            if sp and sp not in seen:
                seen.append(sp)
    return seen


def fetch_taxonomy_xml(taxids: List[int], email: Optional[str]) -> str:
    params = {
        "db": "taxonomy",
        "id": ",".join(str(t) for t in taxids),
        "retmode": "xml",
        "tool": "comparative-robustness",
    }
    if email:
        params["email"] = email
    url = EFETCH + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "comparative-robustness/1.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.read().decode("utf-8", errors="replace")


def parse_lineages(xml_text: str) -> Dict[int, List[Tuple[int, str]]]:
    """taxid -> [(ancestor_taxid, name), ... , (self_taxid, self_name)] root-first."""
    root = ET.fromstring(xml_text)
    out: Dict[int, List[Tuple[int, str]]] = {}
    for taxon in root.findall("Taxon"):
        tid = int(taxon.findtext("TaxId"))
        name = taxon.findtext("ScientificName") or str(tid)
        path: List[Tuple[int, str]] = []
        lex = taxon.find("LineageEx")
        if lex is not None:
            for anc in lex.findall("Taxon"):
                atid = int(anc.findtext("TaxId"))
                aname = anc.findtext("ScientificName") or str(atid)
                path.append((atid, aname))
        path.append((tid, name))
        out[tid] = path
    return out


# --- trie of lineages -> Newick ------------------------------------------------

class Node:
    __slots__ = ("taxid", "name", "children", "leaf_label")

    def __init__(self, taxid: int, name: str):
        self.taxid = taxid
        self.name = name
        self.children: Dict[int, "Node"] = {}
        self.leaf_label: Optional[str] = None


def build_trie(label_to_lineage: Dict[str, List[Tuple[int, str]]]) -> Node:
    root = Node(1, "root")
    for label, lineage in label_to_lineage.items():
        node = root
        for tid, name in lineage:
            if tid not in node.children:
                node.children[tid] = Node(tid, name)
            node = node.children[tid]
        node.leaf_label = label  # mark the tip
    return root


def collapse_unary(node: Node) -> Node:
    """Collapse chains of single-child internal nodes (keep branch structure)."""
    for k in list(node.children):
        node.children[k] = collapse_unary(node.children[k])
    # collapse: if exactly one child and this node is not a labelled tip, splice
    while len(node.children) == 1 and node.leaf_label is None:
        (only,) = node.children.values()
        node.taxid, node.name = only.taxid, only.name
        node.leaf_label = only.leaf_label
        node.children = only.children
    return node


def to_newick(node: Node) -> str:
    if not node.children:
        return node.leaf_label or node.name
    inner = ",".join(f"{to_newick(c)}:1.0" for c in node.children.values())
    label = node.leaf_label or ""
    return f"({inner}){label}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", type=str, default=None)
    ap.add_argument("--usage", type=str, default=os.path.join(DATA_DIR, "usage.tsv"))
    ap.add_argument("--out", type=str, default=os.path.join(DATA_DIR, "tree_ncbi.nwk"))
    args = ap.parse_args()

    wanted = species_in_usage(args.usage)
    label_to_taxid = {name: tid for tid, name in CURATED}
    pairs = [(label, label_to_taxid[label]) for label in wanted if label in label_to_taxid]
    missing = [w for w in wanted if w not in label_to_taxid]
    if missing:
        print(f"  (no taxid mapping for: {missing}) -- skipped")
    taxids = [tid for _l, tid in pairs]
    print(f"Fetching NCBI taxonomy for {len(taxids)} species ...")

    xml_text = fetch_taxonomy_xml(taxids, args.email)
    time.sleep(0.4)
    lineages_by_tid = parse_lineages(xml_text)

    label_to_lineage = {label: lineages_by_tid[tid]
                        for label, tid in pairs if tid in lineages_by_tid}
    got = set(label_to_lineage)
    for label, tid in pairs:
        if tid not in lineages_by_tid:
            print(f"  WARNING: no lineage returned for {label} (taxid {tid})")

    root = collapse_unary(build_trie(label_to_lineage))
    newick = to_newick(root) + ";"
    ensure_data_dir()
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(newick + "\n")
    print(f"Wrote {len(got)} tips -> {args.out}")
    print("Tips:", " ".join(sorted(got)))


if __name__ == "__main__":
    main()
