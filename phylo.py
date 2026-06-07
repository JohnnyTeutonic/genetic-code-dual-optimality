"""
Phylogenetic correction for Model A.

Cross-organism regressions violate independence: species share ancestry. This
module provides the standard fix -- Phylogenetic Generalised Least Squares
(PGLS) under a Brownian-motion model -- without external phylo libraries.

Pipeline:
  parse_newick(s)        -> tree (nested dict with branch lengths)
  vcv_matrix(tree, taxa) -> Brownian-motion covariance C_ij = shared root->MRCA
                            path length (C_ii = root->tip distance)
  pgls(y, X, C)          -> GLS fit y = X beta + e, Cov(e) ∝ C

Use at the SPECIES level: reduce the multilevel usage data to one scalar per
species (the within-species "robustness tilt", i.e. the OLS slope of RSCU on the
per-codon wobble silent-degree), then test whether the mean tilt is positive
after accounting for phylogeny.

Requires numpy only.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

Tree = dict  # {"name": str|None, "length": float, "children": [Tree, ...]}


# ---------------------------------------------------------------------------
# Newick parsing (supports names, branch lengths, nested clades)
# ---------------------------------------------------------------------------

def parse_newick(s: str) -> Tree:
    s = s.strip()
    if s.endswith(";"):
        s = s[:-1]
    pos = 0

    def parse_clade() -> Tree:
        nonlocal pos
        node: Tree = {"name": None, "length": 0.0, "children": []}
        if s[pos] == "(":
            pos += 1  # consume '('
            while True:
                node["children"].append(parse_clade())
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
        # read optional name
        start = pos
        while pos < len(s) and s[pos] not in ",():":
            pos += 1
        name = s[start:pos]
        if name:
            node["name"] = name
        # read optional branch length
        if pos < len(s) and s[pos] == ":":
            pos += 1
            start = pos
            while pos < len(s) and s[pos] not in ",()":
                pos += 1
            node["length"] = float(s[start:pos])
        return node

    return parse_clade()


# ---------------------------------------------------------------------------
# Brownian-motion variance-covariance matrix
# ---------------------------------------------------------------------------

def _tip_depths(tree: Tree) -> Dict[str, float]:
    """Root-to-tip distance for each named tip."""
    depths: Dict[str, float] = {}

    def walk(node: Tree, acc: float) -> None:
        d = acc + node["length"]
        if not node["children"]:
            if node["name"] is not None:
                depths[node["name"]] = d
        else:
            for ch in node["children"]:
                walk(ch, d)

    # root branch length excluded from depth baseline
    for ch in tree["children"] if tree["children"] else []:
        walk(ch, 0.0)
    if not tree["children"] and tree["name"] is not None:
        depths[tree["name"]] = tree["length"]
    return depths


def _shared_depth(tree: Tree, a: str, b: str) -> float:
    """Root-to-MRCA path length shared by tips a and b (off-diagonal C_ab)."""
    def tips_of(node: Tree) -> set:
        if not node["children"]:
            return {node["name"]}
        out = set()
        for ch in node["children"]:
            out |= tips_of(ch)
        return out

    def walk(node: Tree, acc: float) -> float:
        # acc = depth at the START of this node's branch (above it)
        here = acc + node["length"]
        for ch in node["children"]:
            t = tips_of(ch)
            if a in t and b in t:
                return walk(ch, here)  # both descend through this child -> go deeper
        return here  # this node is the MRCA of a and b

    # start below the root (root length not counted)
    return walk(tree, 0.0) if tree["children"] else 0.0


def vcv_matrix(tree: Tree, taxa: List[str]) -> "np.ndarray":
    depths = _tip_depths(tree)
    n = len(taxa)
    C = np.zeros((n, n))
    for i, a in enumerate(taxa):
        C[i, i] = depths[a]
        for j in range(i + 1, n):
            b = taxa[j]
            shared = _shared_depth(tree, a, b)
            C[i, j] = C[j, i] = shared
    return C


# ---------------------------------------------------------------------------
# PGLS (generalised least squares with covariance C)
# ---------------------------------------------------------------------------

def pgls(y: "np.ndarray", X: "np.ndarray", C: "np.ndarray") -> Dict:
    """GLS: y = X beta + e, Cov(e) = sigma^2 C. Returns coefs, SEs, t, p."""
    from scipy import stats

    Cinv = np.linalg.pinv(C)
    XtCi = X.T @ Cinv
    beta = np.linalg.solve(XtCi @ X, XtCi @ y)
    resid = y - X @ beta
    n, k = X.shape
    dof = n - k
    sigma2 = float(resid.T @ Cinv @ resid) / dof
    cov_beta = sigma2 * np.linalg.pinv(XtCi @ X)
    se = np.sqrt(np.diag(cov_beta))
    tvals = beta / se
    pvals = 2 * stats.t.sf(np.abs(tvals), dof)
    return {"beta": beta, "se": se, "t": tvals, "p": pvals, "dof": dof, "sigma2": sigma2}


# ---------------------------------------------------------------------------
# Grafen (1989) branch lengths: replace arbitrary unit edges with lengths
# derived from clade sizes (node height = #descendant tips - 1). Yields an
# ultrametric tree whose VCV reflects nested clade structure rather than the
# raw count of taxonomic ranks.
# ---------------------------------------------------------------------------

def _count_tips(node: Tree) -> int:
    if not node["children"]:
        return 1
    return sum(_count_tips(c) for c in node["children"])


def assign_grafen_lengths(tree: Tree) -> Tree:
    """Set each branch length to (parent_height - node_height) with node height
    = (#descendant tips - 1); tips get height 0. Mutates and returns the tree."""
    def heights(node: Tree) -> int:
        if not node["children"]:
            node["_h"] = 0.0
            return 1
        sz = sum(heights(c) for c in node["children"])
        node["_h"] = float(sz - 1)
        return sz

    heights(tree)

    def setlen(node: Tree, parent_h: float) -> None:
        node["length"] = parent_h - node["_h"]
        for ch in node["children"]:
            setlen(ch, node["_h"])

    setlen(tree, tree["_h"])  # root branch length -> 0
    return tree


# ---------------------------------------------------------------------------
# Pagel's lambda: ML estimate of phylogenetic signal, with the slope re-fit at
# the ML lambda and a likelihood-ratio test of lambda=0 (no signal).
# ---------------------------------------------------------------------------

def _lambda_transform(C: "np.ndarray", lam: float) -> "np.ndarray":
    """Scale off-diagonal covariances by lam, preserving the diagonal."""
    d = np.diag(np.diag(C))
    return lam * C + (1.0 - lam) * d


def _gls_loglik(y: "np.ndarray", X: "np.ndarray", C: "np.ndarray") -> float:
    """ML log-likelihood of the GLS model y = X beta + e, Cov(e) = sigma^2 C."""
    n = X.shape[0]
    sign, logdet = np.linalg.slogdet(C)
    if sign <= 0:
        return -np.inf
    Cinv = np.linalg.pinv(C)
    XtCi = X.T @ Cinv
    beta = np.linalg.solve(XtCi @ X, XtCi @ y)
    resid = y - X @ beta
    sigma2 = float(resid.T @ Cinv @ resid) / n
    if sigma2 <= 0:
        return -np.inf
    return -0.5 * (n * np.log(2.0 * np.pi * sigma2) + logdet + n)


def pgls_lambda(y: "np.ndarray", X: "np.ndarray", C: "np.ndarray") -> Dict:
    """Profile Pagel's lambda over [0,1] by ML, re-fit the slope at the optimum,
    and LR-test lambda against 0 (1 dof). Returns the pgls() dict augmented with
    'lambda', 'loglik', and 'lrt_p' (signal-vs-none)."""
    from scipy import stats

    def negll(lam: float) -> float:
        lam = min(max(lam, 0.0), 1.0)
        return -_gls_loglik(y, X, _lambda_transform(C, lam))

    grid = np.linspace(0.0, 1.0, 41)
    best = min(grid, key=negll)
    lo, hi = max(0.0, best - 0.025), min(1.0, best + 0.025)
    gr = (np.sqrt(5) - 1) / 2
    a, b = lo, hi
    c1, c2 = b - gr * (b - a), a + gr * (b - a)
    for _ in range(40):
        if negll(c1) < negll(c2):
            b, c2 = c2, c1
            c1 = b - gr * (b - a)
        else:
            a, c1 = c1, c2
            c2 = a + gr * (b - a)
    lam = min(max((a + b) / 2, 0.0), 1.0)

    res = pgls(y, X, _lambda_transform(C, lam))
    ll_lam = _gls_loglik(y, X, _lambda_transform(C, lam))
    ll_0 = _gls_loglik(y, X, _lambda_transform(C, 0.0))
    lr = 2.0 * (ll_lam - ll_0)
    res["lambda"] = lam
    res["loglik"] = ll_lam
    res["lrt_p"] = float(stats.chi2.sf(max(lr, 0.0), 1))
    return res


# ---------------------------------------------------------------------------
# Offline self-test
# ---------------------------------------------------------------------------

def _self_test() -> None:
    # ((A:1,B:1):1,(C:1,D:1):1);  ultrametric, root-tip depth = 2 for all tips
    t = parse_newick("((A:1,B:1):1,(C:1,D:1):1);")
    taxa = ["A", "B", "C", "D"]
    C = vcv_matrix(t, taxa)
    # Diagonals = 2 (root->tip); A,B share 1 (their MRCA); A,C share 0.
    assert np.allclose(np.diag(C), [2, 2, 2, 2]), f"diag wrong: {np.diag(C)}"
    assert np.isclose(C[0, 1], 1.0), f"C[A,B]={C[0,1]} expected 1"
    assert np.isclose(C[0, 2], 0.0), f"C[A,C]={C[0,2]} expected 0"
    assert np.isclose(C[2, 3], 1.0), f"C[C,D]={C[2,3]} expected 1"
    # PGLS sanity: constant predictor recovers the (phylo-weighted) mean.
    y = np.array([1.0, 1.2, 0.9, 1.1])
    X = np.ones((4, 1))
    res = pgls(y, X, C)
    assert res["beta"].shape == (1,)

    # Grafen branch lengths: cladogram of 4 tips -> ultrametric, root height 3.
    tg = assign_grafen_lengths(parse_newick("((A,B),(C,D));"))
    Cg = vcv_matrix(tg, taxa)
    assert np.allclose(np.diag(Cg), 3.0), f"grafen diag {np.diag(Cg)} != 3"
    # A,B share root->their-MRCA = 3 - 1 = 2; A,C share root only = 0.
    assert np.isclose(Cg[0, 1], 2.0), f"grafen C[A,B]={Cg[0,1]} expected 2"
    assert np.isclose(Cg[0, 2], 0.0), f"grafen C[A,C]={Cg[0,2]} expected 0"

    # Pagel lambda: data generated with strong covariance -> lambda near 1 and
    # significant LRT; independent noise -> lambda near 0.
    rng = np.random.default_rng(0)
    big = parse_newick("((A:0.05,B:0.05):1,(C:0.05,D:0.05):1);")  # strong within-pair cov
    Cb = vcv_matrix(big, taxa)
    L = np.linalg.cholesky(Cb + 1e-9 * np.eye(4))
    sig = np.array([(L @ rng.standard_normal(4)) for _ in range(400)])
    lam_hi = np.mean([pgls_lambda(s, X, Cb)["lambda"] for s in sig[:40]])
    assert lam_hi > 0.5, f"expected high lambda, got {lam_hi}"
    print("phylo self-test OK: Newick, BM vcv, PGLS, Grafen lengths, Pagel lambda.")


if __name__ == "__main__":
    _self_test()
