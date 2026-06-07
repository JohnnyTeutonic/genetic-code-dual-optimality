"""
Data fetchers for the real-data run of Model A (and phylogenetic correction).

These hit public resources and so are the only networked part of the project.
Everything in metric.py / model_b_reassign.py runs offline.

Targets:
  * Codon usage      -- CoCoPUTs / HIVE-CUTs (per-species codon counts).
                        https://dnahive.fda.gov/dna.cgi?cmd=cuts_main
                        Kazusa mirror format also supported.
  * tRNA gene copies -- GtRNAdb (per-genome tRNA gene counts by anticodon).
                        http://gtrnadb.ucsc.edu/
  * Phylogeny        -- NCBI Taxonomy / TimeTree, for PGLS covariance (TODO).

Bulk downloads should be cached locally under ./data/. The two parsers below
(Kazusa text blocks; CoCoPUTs/HIVE-CUTs wide TSV) convert whichever export you
pull into the long TSV that model_a_usage.py consumes via --input:

    species    codon    count

Run `python fetch_data.py --self-test` to validate the parsers offline.
Phylogeny (NCBI Taxonomy / TimeTree) for PGLS covariance remains a TODO hook.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from typing import Dict, List, Tuple

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

Row = Tuple[str, str, int]  # (species, codon, count)

# Matches a Kazusa codon entry. Handles both the simple format
#   "UUU 17.6( 174557)"
# and the aa=1&style=N format
#   "UUU F 0.46 17.6 (714298)"
# by capturing the codon and the integer count in the following parentheses.
_KAZUSA_RE = re.compile(r"([ACGTU]{3})[^()]*?\(\s*(\d+)\s*\)", re.IGNORECASE)


def ensure_data_dir() -> str:
    os.makedirs(DATA_DIR, exist_ok=True)
    return DATA_DIR


def _to_dna(codon: str) -> str:
    return codon.upper().replace("U", "T")


_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def anticodon_to_codon(anticodon: str) -> str:
    """Cognate codon for a tRNA anticodon (reverse complement, DNA)."""
    ac = _to_dna(anticodon)
    return ac.translate(_COMPLEMENT)[::-1]


# ---------------------------------------------------------------------------
# Kazusa text format (one species block of 64 "codon freq(count)" entries)
# ---------------------------------------------------------------------------

def parse_kazusa_text(text: str, species: str) -> List[Row]:
    """Parse a Kazusa codon-usage text block into (species, codon, count) rows.

    Kazusa entries look like 'UUU 17.6(174557)'; we keep the integer count and
    convert RNA->DNA. Returns one row per codon found (expects 64).
    """
    rows: List[Row] = []
    for codon, count in _KAZUSA_RE.findall(text):
        rows.append((species, _to_dna(codon), int(count)))
    if not rows:
        raise ValueError("No Kazusa-format codon entries found in input text.")
    return rows


def parse_kazusa_file(path: str, species: str) -> List[Row]:
    with open(path, "r", encoding="utf-8") as fh:
        return parse_kazusa_text(fh.read(), species)


# ---------------------------------------------------------------------------
# CoCoPUTs / HIVE-CUTs wide TSV (one row per species; 64 codon columns)
# ---------------------------------------------------------------------------

def parse_cocoputs_tsv(path: str,
                       species_col_candidates=("Species", "Organism", "Name", "Taxid")) -> List[Row]:
    """Melt a CoCoPUTs-style wide table (codon columns of counts) to long rows.

    Detects which columns are codons (3 chars over A/C/G/T/U) and emits
    (species, codon, count) for each. Non-integer / blank cells are skipped.
    """
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Empty CoCoPUTs file.")
        species_col = next((c for c in species_col_candidates if c in reader.fieldnames), None)
        if species_col is None:
            raise ValueError(f"No species column among {species_col_candidates}; "
                             f"got {reader.fieldnames}")
        codon_cols = [c for c in reader.fieldnames
                      if len(c) == 3 and set(c.upper()) <= set("ACGTU")]
        if not codon_cols:
            raise ValueError("No codon columns detected in CoCoPUTs header.")
        rows: List[Row] = []
        for rec in reader:
            sp = str(rec[species_col]).strip()
            if not sp:
                continue
            for col in codon_cols:
                val = (rec.get(col) or "").strip().replace(",", "")
                if not val:
                    continue
                try:
                    count = int(float(val))
                except ValueError:
                    continue
                rows.append((sp, _to_dna(col), count))
        if not rows:
            raise ValueError("CoCoPUTs file parsed but produced no rows.")
        return rows


# ---------------------------------------------------------------------------
# GtRNAdb tRNA gene copy numbers (translational-selection covariate)
# ---------------------------------------------------------------------------

def parse_gtrnadb_tsv(path: str,
                      species_col_candidates=("Species", "Organism", "Genome", "Name"),
                      anticodon_col_candidates=("Anticodon", "AntiCodon", "anticodon"),
                      count_col_candidates=("Count", "GeneCount", "Genes", "count")) -> List[Row]:
    """Parse a long GtRNAdb-style table to (species, codon, trna_gene_count) rows.

    Expects columns identifying species, anticodon, and tRNA gene count. Each
    anticodon is mapped to its cognate codon (reverse complement). The result is
    the per-codon cognate tRNA gene-copy number, a standard translational-
    selection covariate for Model A.
    """
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Empty GtRNAdb file.")

        def pick(cands):
            return next((c for c in cands if c in reader.fieldnames), None)

        sp_c = pick(species_col_candidates)
        ac_c = pick(anticodon_col_candidates)
        ct_c = pick(count_col_candidates)
        if not (sp_c and ac_c and ct_c):
            raise ValueError(f"GtRNAdb needs species/anticodon/count columns; got {reader.fieldnames}")
        rows: List[Row] = []
        for rec in reader:
            sp = str(rec[sp_c]).strip()
            ac = str(rec[ac_c]).strip()
            val = (rec.get(ct_c) or "").strip().replace(",", "")
            if not sp or not ac or not val:
                continue
            if len(ac) != 3 or not (set(ac.upper()) <= set("ACGTU")):
                continue
            try:
                count = int(float(val))
            except ValueError:
                continue
            rows.append((sp, anticodon_to_codon(ac), count))
        if not rows:
            raise ValueError("GtRNAdb file parsed but produced no rows.")
        return rows


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_long_tsv(rows: List[Row], out_path: str, value_col: str = "count") -> int:
    ensure_data_dir()
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species", "codon", value_col])
        w.writerows(rows)
    return len(rows)


# ---------------------------------------------------------------------------
# Offline self-test (no network): validates the Kazusa + CoCoPUTs parsers
# ---------------------------------------------------------------------------

_KAZUSA_SAMPLE = """
UUU 17.6( 174557)  UCU 15.2( 150748)  UAU 12.2( 120882)  UGU 10.6( 105014)
UUC 20.3( 201066)  UCC 17.7( 175477)  UAC 15.3( 151651)  UGC 12.6( 124403)
UUA  7.7(  76470)  UCA 12.2( 120765)  UAA  1.0(   9568)  UGA  1.6(  15710)
UUG 12.9( 128109)  UCG  4.4(  43484)  UAG  0.8(   7886)  UGG 13.2( 130780)
CUU 13.2( 130570)  CCU 17.5( 173522)  CAU 10.9( 107694)  CGU  4.5(  44551)
CUC 19.6( 194021)  CCC 19.8( 196481)  CAC 15.1( 149665)  CGC 10.4( 102931)
CUA  7.2(  71375)  CCA 16.9( 167308)  CAA 12.3( 121853)  CGA  6.2(  61592)
CUG 39.6( 392515)  CCG  6.9(  68328)  CAG 34.2( 339101)  CGG 11.4( 113026)
AUU 16.0( 158774)  ACU 13.1( 129571)  AAU 17.0( 168258)  AGU 12.1( 119714)
AUC 20.8( 206332)  ACC 18.9( 187432)  AAC 19.1( 189185)  AGC 19.5( 193458)
AUA  7.5(  74299)  ACA 15.1( 149581)  AAA 24.4( 241729)  AGA 12.2( 120897)
AUG 22.0( 218014)  ACG  6.1(  60964)  AAG 31.9( 316297)  AGG 12.0( 119007)
GUU 11.0( 109003)  GCU 18.4( 182496)  GAU 21.8( 215929)  GGU 10.8( 107150)
GUC 14.5( 143716)  GCC 27.7( 274876)  GAC 25.1( 248856)  GGC 22.2( 220164)
GUA  7.1(  70337)  GCA 15.8( 156244)  GAA 29.0( 287644)  GGA 16.5( 163524)
GUG 28.1( 278605)  GCG  7.4(  73296)  GAG 39.6( 392381)  GGG 16.4( 162435)
"""

_COCOPUTS_SAMPLE_HEADER = "Species\tTaxid\tTTT\tTTC\tGGG\n"
_COCOPUTS_SAMPLE_ROWS = "Homo sapiens\t9606\t174557\t201066\t162435\nEscherichia coli\t562\t90000\t80000\t70000\n"


def _self_test() -> None:
    k = parse_kazusa_text(_KAZUSA_SAMPLE, "Homo sapiens")
    assert len(k) == 64, f"Kazusa: expected 64 codons, got {len(k)}"
    assert ("Homo sapiens", "TTT", 174557) in k, "Kazusa: TTT count/U->T conversion failed"
    assert all(set(c) <= set("ACGT") for _s, c, _n in k), "Kazusa: RNA not converted"

    tmp = os.path.join(ensure_data_dir(), "_cocoputs_sample.tsv")
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(_COCOPUTS_SAMPLE_HEADER + _COCOPUTS_SAMPLE_ROWS)
    c = parse_cocoputs_tsv(tmp)
    os.remove(tmp)
    assert len(c) == 6, f"CoCoPUTs: expected 6 rows (2 species x 3 codons), got {len(c)}"
    assert ("Homo sapiens", "GGG", 162435) in c, "CoCoPUTs: melt failed"

    # anticodon->codon: tRNA-Ala anticodon AGC reads codon GCT (revcomp).
    assert anticodon_to_codon("AGC") == "GCT", anticodon_to_codon("AGC")
    assert anticodon_to_codon("UGC") == "GCA", anticodon_to_codon("UGC")
    gt = os.path.join(ensure_data_dir(), "_gtrnadb_sample.tsv")
    with open(gt, "w", encoding="utf-8") as fh:
        fh.write("Species\tAnticodon\tCount\n"
                 "Homo sapiens\tAGC\t29\n"
                 "Homo sapiens\tUGC\t5\n")
    g = parse_gtrnadb_tsv(gt)
    os.remove(gt)
    assert ("Homo sapiens", "GCT", 29) in g, f"GtRNAdb mapping failed: {g}"
    print("self-test OK: Kazusa (64 codons), CoCoPUTs (melt), GtRNAdb (anticodon->codon) pass.")


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert codon-usage dumps to species/codon/count TSV.")
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--kazusa", type=str, help="path to a Kazusa codon-usage text block")
    src.add_argument("--cocoputs", type=str, help="path to a CoCoPUTs/HIVE-CUTs wide TSV")
    src.add_argument("--gtrnadb", type=str, help="path to a GtRNAdb species/anticodon/count TSV")
    src.add_argument("--self-test", action="store_true", help="run offline parser self-test")
    ap.add_argument("--species", type=str, default="unknown", help="species label for --kazusa")
    ap.add_argument("--out", type=str, default=os.path.join(DATA_DIR, "usage.tsv"))
    args = ap.parse_args()

    if args.self_test:
        _self_test()
        return
    if args.kazusa:
        rows = parse_kazusa_file(args.kazusa, args.species)
        n = write_long_tsv(rows, args.out)
    elif args.cocoputs:
        rows = parse_cocoputs_tsv(args.cocoputs)
        n = write_long_tsv(rows, args.out)
    elif args.gtrnadb:
        rows = parse_gtrnadb_tsv(args.gtrnadb)
        out = args.out if args.out != os.path.join(DATA_DIR, "usage.tsv") else os.path.join(DATA_DIR, "trna.tsv")
        n = write_long_tsv(rows, out, value_col="trna_copy")
        print(f"Wrote {n} rows to {out}")
        return
    else:
        print(__doc__)
        return
    print(f"Wrote {n} rows to {args.out}")


if __name__ == "__main__":
    main()
