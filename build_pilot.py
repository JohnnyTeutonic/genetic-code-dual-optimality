"""
Build a REAL pilot codon-usage dataset by downloading per-species tables from
the Kazusa Codon Usage Database (GenBank-derived counts).

This fetches genuine data at runtime -- nothing here is synthesised. Species
whose Kazusa table is missing/empty are skipped and reported. Output:
data/usage.tsv with columns species, codon, count (species labelled by binomial
so it matches data/tree.nwk).

Endpoint (numeric counts incl. amino-acid annotation):
    https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?aa=1&style=N&species=<taxid>

Usage:
    python build_pilot.py            # downloads the default curated set
"""

from __future__ import annotations

import os
import sys
import time
import urllib.request

from fetch_data import DATA_DIR, parse_kazusa_text, write_long_tsv

KAZUSA_URL = "https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?aa=1&style=N&species={taxid}"

# Curated, phylogenetically spread set: (Kazusa GenBank taxid, binomial label).
# Strain-level taxids are used where the species-level node has no table.
CURATED = [
    (9606,  "Homo_sapiens"),
    (10090, "Mus_musculus"),
    (9913,  "Bos_taurus"),
    (9031,  "Gallus_gallus"),
    (8355,  "Xenopus_laevis"),
    (7955,  "Danio_rerio"),
    (7227,  "Drosophila_melanogaster"),
    (6239,  "Caenorhabditis_elegans"),
    (7460,  "Apis_mellifera"),
    (4932,  "Saccharomyces_cerevisiae"),
    (4896,  "Schizosaccharomyces_pombe"),
    (5141,  "Neurospora_crassa"),
    (3702,  "Arabidopsis_thaliana"),
    (4530,  "Oryza_sativa"),
    (3055,  "Chlamydomonas_reinhardtii"),
    (83333, "Escherichia_coli_K12"),
    (1423,  "Bacillus_subtilis"),
    (1773,  "Mycobacterium_tuberculosis"),
    (1148,  "Synechocystis_sp_PCC6803"),
    (36329, "Plasmodium_falciparum_3D7"),
    (44689, "Dictyostelium_discoideum"),
    (2287,  "Sulfolobus_solfataricus"),
    (2190,  "Methanocaldococcus_jannaschii"),
]


def fetch_one(taxid: int, timeout: float = 30.0) -> str:
    url = KAZUSA_URL.format(taxid=taxid)
    req = urllib.request.Request(url, headers={"User-Agent": "comparative-robustness-pilot/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def main() -> None:
    rows = []
    kept, skipped = [], []
    for taxid, name in CURATED:
        try:
            text = fetch_one(taxid)
            recs = parse_kazusa_text(text, name)
        except Exception as e:  # network error or no table
            skipped.append((name, taxid, repr(e)[:80]))
            continue
        if len(recs) < 60:  # Kazusa "Not found" yields ~0 codon matches
            skipped.append((name, taxid, f"only {len(recs)} codons parsed"))
            continue
        rows.extend(recs)
        kept.append(name)
        print(f"  ok   {name:<32} taxid={taxid:<7} codons={len(recs)}")
        time.sleep(0.5)  # be polite to Kazusa

    if not rows:
        sys.exit("No species fetched -- check network/endpoint.")

    out = os.path.join(DATA_DIR, "usage.tsv")
    n = write_long_tsv(rows, out)
    print(f"\nKept {len(kept)} species, {n} rows -> {out}")
    if skipped:
        print("Skipped:")
        for name, taxid, why in skipped:
            print(f"  - {name} (taxid {taxid}): {why}")
    print("\nSpecies labels (use these in data/tree.nwk):")
    print("  " + " ".join(kept))


if __name__ == "__main__":
    main()
