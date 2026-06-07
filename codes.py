"""
Genetic codes used by the comparative analysis.

STANDARD_CODE is NCBI translation table 1. Variant codes are encoded as
overrides on the standard code, taken from the NCBI genetic-code tables
(https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). Each variant is a
documented natural reassignment and constitutes one "natural experiment" in
which block sizes differ from the standard code.

Amino acids use three-letter codes; "Stop" denotes a termination signal,
consistent with metric.py.
"""

from __future__ import annotations

from typing import Dict

Code = Dict[str, str]

# NCBI table 1 -- the Standard Code (DNA, T/C/A/G).
STANDARD_CODE: Code = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "Stop", "TAG": "Stop",
    "TGT": "Cys", "TGC": "Cys", "TGA": "Stop", "TGG": "Trp",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

# Variant codes as {NCBI table id: (name, overrides)}.
VARIANT_OVERRIDES = {
    2:  ("Vertebrate Mitochondrial",
         {"AGA": "Stop", "AGG": "Stop", "ATA": "Met", "TGA": "Trp"}),
    3:  ("Yeast Mitochondrial",
         {"ATA": "Met", "CTT": "Thr", "CTC": "Thr", "CTA": "Thr", "CTG": "Thr", "TGA": "Trp"}),
    4:  ("Mold/Protozoan/Coelenterate Mito & Mycoplasma/Spiroplasma",
         {"TGA": "Trp"}),
    5:  ("Invertebrate Mitochondrial",
         {"AGA": "Ser", "AGG": "Ser", "ATA": "Met", "TGA": "Trp"}),
    6:  ("Ciliate/Dasycladacean/Hexamita Nuclear",
         {"TAA": "Gln", "TAG": "Gln"}),
    9:  ("Echinoderm/Flatworm Mitochondrial",
         {"AAA": "Asn", "AGA": "Ser", "AGG": "Ser", "TGA": "Trp"}),
    10: ("Euplotid Nuclear",
         {"TGA": "Cys"}),
    12: ("Alternative Yeast Nuclear (Candida)",
         {"CTG": "Ser"}),
    13: ("Ascidian Mitochondrial",
         {"AGA": "Gly", "AGG": "Gly", "ATA": "Met", "TGA": "Trp"}),
    14: ("Alternative Flatworm Mitochondrial",
         {"AAA": "Asn", "AGA": "Ser", "AGG": "Ser", "TAA": "Tyr", "TGA": "Trp"}),
    16: ("Chlorophycean Mitochondrial",
         {"TAG": "Leu"}),
    21: ("Trematode Mitochondrial",
         {"TGA": "Trp", "ATA": "Met", "AGA": "Ser", "AGG": "Ser", "AAA": "Asn"}),
    22: ("Scenedesmus obliquus Mitochondrial",
         {"TCA": "Stop", "TAG": "Leu"}),
    24: ("Rhabdopleuridae Mitochondrial",
         {"AGA": "Ser", "AGG": "Lys", "TGA": "Trp"}),
    25: ("Candidate Division SR1 / Gracilibacteria",
         {"TGA": "Gly"}),
    33: ("Cephalodiscidae Mitochondrial",
         {"TAA": "Tyr", "AGA": "Ser", "AGG": "Lys", "TGA": "Trp"}),
}


# Broad lineage group for each variant table, used to count INDEPENDENT origins
# of a reassignment (a reassignment shared across distinct groups = convergent =
# independent evidence). Coarse, deliberately conservative; documented for the
# write-up. Mitochondrial codes are grouped by host lineage.
LINEAGE_GROUP = {
    2:  "vertebrate_mito",
    3:  "fungal_mito",
    4:  "protist_mito_and_mollicutes",
    5:  "invertebrate_mito",
    6:  "ciliate_nuclear",
    9:  "echinoderm_flatworm_mito",
    10: "euplotid_nuclear",
    12: "candida_nuclear",
    13: "ascidian_mito",
    14: "flatworm_mito",
    16: "chlorophyte_mito",
    21: "trematode_mito",
    22: "chlorophyte_mito",      # Scenedesmus = chlorophyte (same group as 16)
    24: "pterobranch_mito",
    25: "gracilibacteria",
    33: "pterobranch_mito",      # Cephalodiscidae = pterobranch (same group as 24)
}


def variant_code(table_id: int) -> Code:
    """Materialise a variant code as standard + documented overrides."""
    name, overrides = VARIANT_OVERRIDES[table_id]
    code = dict(STANDARD_CODE)
    code.update(overrides)
    return code


def all_variant_codes() -> Dict[int, "tuple[str, Code]"]:
    return {tid: (name, variant_code(tid)) for tid, (name, _) in VARIANT_OVERRIDES.items()}
