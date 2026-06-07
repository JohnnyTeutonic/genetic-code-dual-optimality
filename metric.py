"""
Faithful re-implementation of the submitted wobble paper's error-resistance
functional, plus the per-codon / reassignment marginal decompositions that the
comparative regressions (Models A and B) require.

The functional is the SILENT-EDGE COUNT of a genetic code, exactly as defined in
the paper's Lean development and `robustness_scoring.py`:

    silentEdge(c1, c2)  iff  c1, c2 differ at exactly one position AND code[c1] == code[c2]

    S(code)             = number of undirected silent single-substitution edges
    S_pos(code, p)      = silent edges restricted to mutations at position p
    S_ti(code)          = silent edges that are transitions (purine<->purine, pyr<->pyr)

The Lean theorems prove that third-position (wobble) silent edges are maximised
(= 96) exactly by wobble-fibre-constancy, and silent transitions are maximised
(= 32) under the same condition. This module reproduces S, S_pos, S_ti and adds:

    silent_degree(code, codon)      per-codon contribution (Model A predictor)
    delta_S(code, codon, new_aa)    marginal cost of a reassignment (Model B predictor)

Self-contained: standard library only.
"""

from __future__ import annotations

from collections import defaultdict
from itertools import product
from typing import Dict, Iterable, List, Optional, Tuple

Code = Dict[str, str]

BASES: Tuple[str, ...] = ("T", "C", "A", "G")

# Transitions: pyrimidine<->pyrimidine (T<->C), purine<->purine (A<->G).
TRANSITION_PAIRS = {("T", "C"), ("C", "T"), ("A", "G"), ("G", "A")}


# ---------------------------------------------------------------------------
# Codon-space geometry
# ---------------------------------------------------------------------------

def all_codons(alphabet: Iterable[str] = BASES, length: int = 3) -> List[str]:
    return ["".join(c) for c in product(tuple(alphabet), repeat=length)]


def neighbours_at_position(codon: str, pos: int, alphabet: Iterable[str] = BASES) -> List[str]:
    """All codons differing from `codon` at exactly position `pos`."""
    return [codon[:pos] + b + codon[pos + 1:] for b in alphabet if b != codon[pos]]


def is_transition(c1: str, c2: str, pos: int) -> bool:
    return (c1[pos], c2[pos]) in TRANSITION_PAIRS


# ---------------------------------------------------------------------------
# The functional S and its restrictions
# ---------------------------------------------------------------------------

def silent_edges(code: Code, position: Optional[int] = None,
                 alphabet: Iterable[str] = BASES, codon_length: int = 3) -> int:
    """Undirected count of silent single-substitution edges.

    position=None -> all positions; otherwise restrict to that position.
    Each undirected edge counted once (lexicographic c < nb guard).
    """
    alphabet = tuple(alphabet)
    positions = range(codon_length) if position is None else [position]
    total = 0
    for pos in positions:
        for c in code:
            for nb in neighbours_at_position(c, pos, alphabet):
                if nb in code and nb > c and code[c] == code[nb]:
                    total += 1
    return total


def silent_transitions(code: Code, position: Optional[int] = None,
                       alphabet: Iterable[str] = BASES, codon_length: int = 3) -> int:
    """Silent edges that are transitions."""
    alphabet = tuple(alphabet)
    positions = range(codon_length) if position is None else [position]
    total = 0
    for pos in positions:
        for c in code:
            for nb in neighbours_at_position(c, pos, alphabet):
                if nb in code and nb > c and code[c] == code[nb] and is_transition(c, nb, pos):
                    total += 1
    return total


def functional_report(code: Code, alphabet: Iterable[str] = BASES, codon_length: int = 3) -> Dict:
    """Full S report, mirroring the paper's headline numbers."""
    by_pos = {p: silent_edges(code, p, alphabet, codon_length) for p in range(codon_length)}
    ti_by_pos = {p: silent_transitions(code, p, alphabet, codon_length) for p in range(codon_length)}
    return {
        "total_silent_edges": sum(by_pos.values()),
        "silent_edges_by_position": by_pos,
        "third_position_silent_edges": by_pos.get(codon_length - 1, 0),
        "total_silent_transitions": sum(ti_by_pos.values()),
        "third_position_silent_transitions": ti_by_pos.get(codon_length - 1, 0),
        "distinct_signals": len(set(code.values())),
    }


# ---------------------------------------------------------------------------
# Per-codon decomposition (Model A predictor)
# ---------------------------------------------------------------------------

def silent_degree(code: Code, codon: str, alphabet: Iterable[str] = BASES,
                  codon_length: int = 3) -> int:
    """Number of synonymous Hamming-1 neighbours of `codon` (directed count).

    This is the per-codon contribution to 2*S: sum over all codons of
    silent_degree == 2 * silent_edges(code). It is the natural per-codon
    error-robustness predictor for the within-block usage regression.
    """
    alphabet = tuple(alphabet)
    aa = code.get(codon)
    if aa is None:
        return 0
    deg = 0
    for pos in range(codon_length):
        for nb in neighbours_at_position(codon, pos, alphabet):
            if code.get(nb) == aa:
                deg += 1
    return deg


def silent_degree_by_position(code: Code, codon: str, pos: int,
                              alphabet: Iterable[str] = BASES) -> int:
    """Silent degree of `codon` restricted to substitutions at `pos`
    (the third-position component isolates wobble robustness)."""
    alphabet = tuple(alphabet)
    aa = code.get(codon)
    if aa is None:
        return 0
    return sum(1 for nb in neighbours_at_position(codon, pos, alphabet) if code.get(nb) == aa)


def per_codon_table(code: Code, alphabet: Iterable[str] = BASES, codon_length: int = 3) -> List[Dict]:
    """Per-codon predictor table for Model A."""
    rows = []
    for codon in sorted(code):
        aa = code[codon]
        rows.append({
            "codon": codon,
            "amino_acid": aa,
            "block_size": sum(1 for v in code.values() if v == aa),
            "silent_degree": silent_degree(code, codon, alphabet, codon_length),
            "silent_degree_wobble": silent_degree_by_position(code, codon, codon_length - 1, alphabet),
            "gc3": 1 if codon[-1] in ("G", "C") else 0,
            "third_base": codon[-1],
        })
    return rows


# ---------------------------------------------------------------------------
# Reassignment marginal (Model B predictor)
# ---------------------------------------------------------------------------

def delta_S(code: Code, codon: str, new_aa: str,
            alphabet: Iterable[str] = BASES, codon_length: int = 3) -> int:
    """Change in total silent edges when `codon` is reassigned to `new_aa`.

    Only edges incident to `codon` change, so this is computed locally rather
    than by recomputing S over the whole code. Negative => robustness cost.
    """
    alphabet = tuple(alphabet)
    old_aa = code.get(codon)
    if old_aa is None or old_aa == new_aa:
        return 0
    before = after = 0
    for pos in range(codon_length):
        for nb in neighbours_at_position(codon, pos, alphabet):
            nb_aa = code.get(nb)
            if nb_aa is None:
                continue
            before += 1 if nb_aa == old_aa else 0
            after += 1 if nb_aa == new_aa else 0
    return after - before


def reassignment_events(reference: Code, variant: Code) -> List[Dict]:
    """Codon-level differences between a reference (standard) and a variant code.

    For each reassigned codon, reports the realised ΔS (cost/benefit under the
    paper's functional) of moving from the reference assignment to the variant
    assignment, evaluated on the reference background.
    """
    events = []
    for codon, ref_aa in reference.items():
        var_aa = variant.get(codon)
        if var_aa is not None and var_aa != ref_aa:
            events.append({
                "codon": codon,
                "from_aa": ref_aa,
                "to_aa": var_aa,
                "delta_S": delta_S(reference, codon, var_aa),
                "gc3": 1 if codon[-1] in ("G", "C") else 0,
            })
    return events


if __name__ == "__main__":
    from codes import STANDARD_CODE

    rep = functional_report(STANDARD_CODE)
    print("Standard code S-functional:")
    for k, v in rep.items():
        print(f"  {k}: {v}")
    # Sanity: per-codon degrees sum to 2*S
    deg_sum = sum(silent_degree(STANDARD_CODE, c) for c in STANDARD_CODE)
    print(f"  sum(silent_degree) = {deg_sum}  (should equal 2*total_silent_edges = {2*rep['total_silent_edges']})")
