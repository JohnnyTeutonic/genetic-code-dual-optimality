"""
DIRECTION #2, genome-scale empirical leg.

Test the prediction of the dual-optimality theorem (Remark 1): genomes sit ABOVE
the theoretical minimum decoding-set size by a redundancy factor that SHRINKS in
reduced genomes. We pull a broad GtRNAdb sample across the three domains, compute
each genome's decoding repertoire (distinct anticodons), fetch genome size from
NCBI, and regress the redundancy factor on log genome size and domain.

Pipeline (pure stdlib + numpy/statsmodels), resumable: fetched rows are cached in
data/trna_scale.tsv; rerun analyses without refetching unless --refetch.

Usage:
    python scale_trna.py --per-domain 25     # fetch + analyse
    python scale_trna.py                      # analyse cached tsv
"""

from __future__ import annotations

import argparse
import json
import os
import re
import time
import urllib.parse
import urllib.request
from collections import OrderedDict
from typing import Dict, List, Optional, Tuple

import fetch_trna
from wobble_decoding import decoding_set_size
import codes

OUT = os.path.join("data", "trna_scale.tsv")
MIN_SUPER = decoding_set_size(codes.STANDARD_CODE, "superwobble")  # 23
ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def domain_of(summary_url: str) -> str:
    m = re.search(r"/genomes/(\w+)/", summary_url)
    return m.group(1) if m else "?"


def clean_name(raw: str) -> str:
    return re.sub(r"\s*\(.*?\)\s*", "", raw).strip()


def sample_genomes(index, per_domain: int):
    """One genome per organism, up to per_domain per GtRNAdb domain, STRIDE-sampled
    evenly across each domain's (alphabetically ordered) genome list so the sample
    spreads taxonomically rather than clustering on the first few clades."""
    # collect unique organisms per domain, preserving index order
    by_domain: Dict[str, "OrderedDict[str,str]"] = {}
    for name, url in index:
        dom = domain_of(url)
        org = clean_name(name)
        d = by_domain.setdefault(dom, OrderedDict())
        if org not in d:
            d[org] = url
    out = []
    for dom, d in by_domain.items():
        items = list(d.items())
        if per_domain >= len(items):
            chosen = items
        else:
            step = len(items) / per_domain
            idx = sorted({int(i * step) for i in range(per_domain)})
            chosen = [items[i] for i in idx]
        for org, url in chosen:
            out.append((org, dom, url))
    return out


def ncbi_genome_size(organism: str) -> Optional[int]:
    """Best-effort total assembly length (bp) via NCBI assembly esummary."""
    try:
        q = urllib.parse.urlencode({"db": "assembly",
                                    "term": f"{organism}[Organism]",
                                    "retmax": "1", "retmode": "json",
                                    "tool": "comparative-robustness"})
        ids = json.loads(fetch_trna._get(ESEARCH + "?" + q).decode("utf-8", "replace"))
        idlist = ids.get("esearchresult", {}).get("idlist", [])
        if not idlist:
            return None
        q2 = urllib.parse.urlencode({"db": "assembly", "id": idlist[0],
                                     "retmode": "xml", "tool": "comparative-robustness"})
        summ = fetch_trna._get(ESUMMARY + "?" + q2).decode("utf-8", "replace")
        # stats live in the <Meta> XML block; first total_length is sequence_tag="all"
        m = re.search(r'category="total_length"[^>]*>(\d+)<', summ)
        return int(m.group(1)) if m else None
    except Exception:
        return None


def fetch(per_domain: int) -> None:
    print("Scraping GtRNAdb index ...")
    index = fetch_trna.scrape_index()
    sample = sample_genomes(index, per_domain)
    print(f"  sampled {len(sample)} genomes across domains")
    fetch_trna.ensure = None
    rows = []
    for org, dom, summary in sample:
        try:
            tar = fetch_trna.tarball_url_from_summary(summary)
            if not tar:
                continue
            counts = fetch_trna.anticodon_counts_from_tarball(tar)
            if not counts:
                continue
            size = ncbi_genome_size(org)
            rows.append((org, dom, len(counts), sum(counts.values()), size or ""))
            print(f"  ok {org[:34]:<34} {dom:<9} rep={len(counts):>3} "
                  f"genes={sum(counts.values()):>4} size={size}")
            time.sleep(0.4)
        except Exception as e:
            print(f"  skip {org[:40]}: {repr(e)[:50]}")
    os.makedirs("data", exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("organism\tdomain\trepertoire\ttotal_genes\tgenome_size_bp\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"\nWrote {len(rows)} rows -> {OUT}")


def analyse() -> None:
    if not os.path.exists(OUT):
        print(f"No {OUT}; run with --per-domain N first.")
        return
    import numpy as np
    import pandas as pd
    import statsmodels.formula.api as smf

    df = pd.read_csv(OUT, sep="\t")
    df["redundancy"] = df["repertoire"] / MIN_SUPER
    print(f"theoretical min superwobble decoding set = {MIN_SUPER}\n")
    print("Repertoire & redundancy by domain:")
    print(df.groupby("domain")[["repertoire", "redundancy", "total_genes"]]
          .agg(["count", "mean"]).round(2))

    sized = df.dropna(subset=["genome_size_bp"]).copy()
    sized = sized[sized["genome_size_bp"] != ""]
    sized["genome_size_bp"] = sized["genome_size_bp"].astype(float)
    sized = sized[sized["genome_size_bp"] > 0]
    sized["log_size"] = np.log10(sized["genome_size_bp"])
    print(f"\nGenome size resolved for {len(sized)}/{len(df)} genomes.")
    if len(sized) >= 8:
        m = smf.ols("redundancy ~ log_size", sized).fit()
        print(f"\n[pooled] redundancy ~ log10(genome_size):")
        print(f"  slope = {m.params['log_size']:.3f}  "
              f"p = {m.pvalues['log_size']:.3g}  R2 = {m.rsquared:.3f}  n={int(m.nobs)}")
        rho = sized[["log_size", "redundancy"]].corr(method="spearman").iloc[0, 1]
        print(f"  Spearman(log_size, redundancy) = {rho:+.3f}")

        # CONFOUND CONTROL 1: genome size as covariate alongside domain
        try:
            md = smf.ols("redundancy ~ log_size + C(domain)", sized).fit()
            print(f"\n[domain-controlled] redundancy ~ log_size + C(domain):")
            print(f"  log_size slope = {md.params['log_size']:.3f}  "
                  f"p = {md.pvalues['log_size']:.3g}  (n={int(md.nobs)}, R2={md.rsquared:.3f})")
            print("  => if log_size stays significant here, the size effect is NOT")
            print("     merely the prokaryote/eukaryote contrast.")
        except Exception as e:
            print(f"  (domain-controlled model failed: {e})")

        # CONFOUND CONTROL 2: within-domain regressions (fully confound-free)
        print("\n[within-domain] redundancy ~ log_size, per domain:")
        for dom, sub in sized.groupby("domain"):
            if len(sub) >= 8 and sub["log_size"].std() > 0:
                mm = smf.ols("redundancy ~ log_size", sub).fit()
                rr = sub[["log_size", "redundancy"]].corr(method="spearman").iloc[0, 1]
                print(f"  {dom:<10} n={len(sub):>3}  slope={mm.params['log_size']:+.3f}  "
                      f"p={mm.pvalues['log_size']:.3g}  Spearman={rr:+.3f}  "
                      f"size_range=[{sub['genome_size_bp'].min()/1e6:.1f},"
                      f"{sub['genome_size_bp'].max()/1e6:.0f}]Mb")
            else:
                print(f"  {dom:<10} n={len(sub):>3}  (insufficient size spread)")
        print("\n  Prediction: positive within-domain slopes => larger genomes carry")
        print("  more decoding redundancy independent of the domain contrast, as the")
        print("  superwobble dual-optimality (Remark 1) predicts.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-domain", type=int, default=0,
                    help="genomes per domain to fetch (0 = analyse cached only)")
    ap.add_argument("--refetch", action="store_true")
    args = ap.parse_args()
    if args.per_domain and (args.refetch or not os.path.exists(OUT)):
        fetch(args.per_domain)
    analyse()


if __name__ == "__main__":
    main()
