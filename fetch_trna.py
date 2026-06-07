"""
Fully automated tRNA gene-copy puller from GtRNAdb (pure standard library, no
Selenium, no UI).

Pipeline per species:
  1. scrape the GtRNAdb genome index (browse.html) for all genome summary pages
     and their organism names;
  2. match each pilot species (genus+species) to a genome entry;
  3. on the matched summary page, find the static tRNAscan-SE results tarball
     ('<assembly>-tRNAs.tar.gz');
  4. download it, read the FASTA member, count tRNA genes per anticodon from the
     headers (e.g. '...tRNA-Ala-AGC-1-1');
  5. map anticodon -> cognate codon and write data/trna.tsv
     (columns: species, codon, trna_copy) for model_a_usage.py --trna.

Best-effort: species that cannot be matched (GtRNAdb naming/coverage gaps,
esp. some prokaryotes) are reported and skipped. tRNA feeds only Model A.

Usage:
    python fetch_trna.py                      # all pilot species in data/usage.tsv
    python fetch_trna.py --only Homo_sapiens  # one species (validate the path)
"""

from __future__ import annotations

import argparse
import gzip
import io
import os
import re
import sys
import tarfile
import time
import urllib.request
from collections import Counter
from typing import Dict, List, Optional, Tuple

from build_pilot import CURATED
from fetch_data import DATA_DIR, anticodon_to_codon, ensure_data_dir, write_long_tsv

BROWSE_URL = "https://gtrnadb.ucsc.edu/browse.html"
BASE = "https://gtrnadb.ucsc.edu"

# Known-good overrides (label -> summary URL) so the path is verifiable even if
# index naming drifts.
OVERRIDES = {
    "Homo_sapiens": "https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-summary.html",
}

# browse.html anchors look like:
#   <a href="genomes/eukaryota/Hsapi38/">Homo sapiens (GRCh38/hg38)</a>
_GENOME_RE = re.compile(
    r'href="(?:https?://[^"]+)?/?genomes/(\w+)/([^"/]+)/"[^>]*>([^<]+)</a>',
    re.IGNORECASE)
_TARBALL_RE = re.compile(r'href="([^"]*?-tRNAs\.tar\.gz)"', re.IGNORECASE)
_ANTICODON_RE = re.compile(r'tRNA-[A-Za-z]+-([ACGT]{3})', re.IGNORECASE)


def _get(url: str, timeout: float = 60.0) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "comparative-robustness/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def scrape_index() -> List[Tuple[str, str]]:
    """Return [(organism_name_lower, summary_url), ...] from browse.html."""
    html = _get(BROWSE_URL).decode("utf-8", errors="replace")
    out = []
    for m in _GENOME_RE.finditer(html):
        domain, code, name = m.group(1), m.group(2), m.group(3)
        summary_url = f"{BASE}/genomes/{domain}/{code}/{code}-summary.html"
        out.append((name.strip().lower(), summary_url))
    return out


def match_species(label: str, index: List[Tuple[str, str]]) -> Optional[str]:
    """Match 'Genus_species[_strain]' to an index entry by genus+species tokens."""
    toks = label.lower().split("_")
    if len(toks) < 2:
        return None
    genus, species = toks[0], toks[1]
    # exact genus+species containment; prefer the shortest matching name
    cands = [url for name, url in index if genus in name and species in name]
    if not cands:
        # genus-only fallback (e.g. strain naming differences)
        cands = [url for name, url in index if genus in name]
    return min(cands, key=len) if cands else None


def tarball_url_from_summary(summary_url: str) -> Optional[str]:
    html = _get(summary_url).decode("utf-8", errors="replace")
    m = _TARBALL_RE.search(html)
    if not m:
        return None
    url = m.group(1)
    if url.startswith("/"):
        url = BASE + url
    elif not url.startswith("http"):
        url = summary_url.rsplit("/", 1)[0] + "/" + url
    return url


def anticodon_counts_from_tarball(tar_url: str) -> Counter:
    raw = _get(tar_url)
    counts: Counter = Counter()
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as tf:
        # prefer the mature/primary FASTA member
        fa_members = [m for m in tf.getmembers() if m.name.lower().endswith(".fa")]
        # rank: plain '-tRNAs.fa' before 'mature'/'pre' variants
        fa_members.sort(key=lambda m: (("mature" in m.name.lower()),
                                       ("pre" in m.name.lower()), len(m.name)))
        if not fa_members:
            return counts
        f = tf.extractfile(fa_members[0])
        if f is None:
            return counts
        text = f.read().decode("utf-8", errors="replace")
    for line in text.splitlines():
        if line.startswith(">"):
            m = _ANTICODON_RE.search(line)
            if m:
                counts[m.group(1).upper()] += 1
    return counts


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", type=str, default=None, help="restrict to one species label")
    ap.add_argument("--usage", type=str, default=os.path.join(DATA_DIR, "usage.tsv"))
    ap.add_argument("--out", type=str, default=os.path.join(DATA_DIR, "trna.tsv"))
    args = ap.parse_args()

    # species set
    labels = [name for _tid, name in CURATED]
    if os.path.exists(args.usage):
        present = set()
        with open(args.usage, "r", encoding="utf-8") as fh:
            next(fh, None)
            for line in fh:
                present.add(line.split("\t", 1)[0])
        labels = [l for l in labels if l in present]
    if args.only:
        labels = [args.only]

    print("Scraping GtRNAdb genome index ...")
    try:
        index = scrape_index()
        print(f"  index entries: {len(index)}")
    except Exception as e:
        print(f"  index scrape failed ({e!r}); relying on OVERRIDES only")
        index = []

    rows: List[Tuple[str, str, int]] = []
    kept, skipped = [], []
    for label in labels:
        summary = OVERRIDES.get(label) or match_species(label, index)
        if not summary:
            skipped.append((label, "no index match"))
            continue
        try:
            tar = tarball_url_from_summary(summary)
            if not tar:
                skipped.append((label, "no tarball link"))
                continue
            counts = anticodon_counts_from_tarball(tar)
            if not counts:
                skipped.append((label, "no anticodons parsed"))
                continue
        except Exception as e:
            skipped.append((label, repr(e)[:70]))
            continue
        for ac, n in counts.items():
            rows.append((label, anticodon_to_codon(ac), n))
        kept.append(label)
        print(f"  ok   {label:<32} anticodons={len(counts)} genes={sum(counts.values())}")
        time.sleep(0.5)

    if rows:
        ensure_data_dir()
        n = write_long_tsv(rows, args.out, value_col="trna_copy")
        print(f"\nKept {len(kept)} species, {n} rows -> {args.out}")
    else:
        print("\nNo tRNA data fetched.")
    if skipped:
        print("Skipped:")
        for label, why in skipped:
            print(f"  - {label}: {why}")


if __name__ == "__main__":
    main()
