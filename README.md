# Two faces of one optimum

**Wobble error-resistance and minimal tRNA decoding coincide in the genetic code.**

Code and data accompanying the paper by **Jonathan Reich** (Independent Researcher).

The standard genetic code has been called "optimal" in two separate literatures:
it is highly error-resistant (synonymous codons cluster under single mutations),
and it is decoded by a remarkably small set of tRNAs. This work proves these are
the **same optimum**, and confirms the resulting prediction across 390 genomes.

## The result

**Theorem (dual optimality).** Over all codon-to-symbol maps, the third-position
silent-edge count `S_3` satisfies `S_3 <= 96`, and the superwobble decoding number
`D_sw` satisfies `D_sw >= 16`, with equality in *both* if and only if the code is
**wobble-fibre constant** (every two-base codon prefix is a single four-codon
family box). Hence `argmin D_sw = argmax S_3`: the minimal-tRNA codes are exactly
the maximally wobble-robust codes. (Per-box lemma: merging two amino-acid groups
of sizes `a, b` raises that box's `S_3` by `ab` and lowers its `D_sw` by 1.)

The decoding model reproduces the classical **31-tRNA Crick minimum** and the
**22-tRNA vertebrate-mitochondrial set**; the standard code scores `S_3 = 64`.

**Prediction and test.** The coupling is specific to *superwobble* (the documented
four-way decoding used by organelles and reduced genomes), so decoding redundancy
above the theoretical minimum should grow with genome size. Across 390 genomes
(GtRNAdb tRNA gene sets; NCBI assembly sizes) it does:

| Analysis | slope / decade | p |
|---|---|---|
| OLS, pooled | +0.15 | 2e-18 |
| + domain covariate | +0.21 | 9e-15 |
| within bacteria / eukaryotes / archaea | positive | all significant |
| PGLS, Brownian (unit edges) | +0.26 | 6e-10 |
| PGLS, Grafen branch lengths | +0.29 | 8e-9 |
| PGLS, Pagel's lambda (ML, lambda = 0.83) | +0.25 | 6e-9 |

The size effect is positive and significant under every phylogenetic branch-length
model; the ML phylogenetic signal `lambda = 0.83` (LR test vs 0: `p = 5e-25`)
shows the effect survives strong phylogenetic non-independence.

## Layout

```
wobble_decoding.py   S_3 / D_sw engine + theorem & decoding-model validation
metric.py            silent-edge functional (stdlib only)
codes.py             standard and variant genetic-code tables
phylo.py             Newick / Brownian VCV / Grafen lengths / Pagel's lambda (+self-test)
scale_trna.py        genome-size decoding-redundancy regression (OLS, domain, within-domain)
scale_pgls.py        phylogenetically controlled regression (unit / Grafen / Pagel lambda)
build_tree.py        NCBI-taxonomy tree builder
fetch_trna.py        GtRNAdb tRNA-set fetcher
build_pilot.py       curated genome list (fetch helper)
fetch_data.py        data-IO helpers
data/                cached inputs so analyses run offline
paper/               manuscript (PDF + LaTeX source)
```

## Requirements

Python 3.9+ and (see `requirements.txt`):

```
pip install -r requirements.txt   # numpy, scipy, pandas, statsmodels
```

## Reproduce

```bash
# 1. Theorem engine + decoding-model validation (Crick 31, mito 22, standard S_3 = 64)
python wobble_decoding.py

# 2. Phylogenetics self-test (Newick, Brownian VCV, Grafen lengths, Pagel's lambda)
python phylo.py

# 3. Genome-size redundancy regression, from cached data, offline
#    (pooled OLS, domain-controlled, and within each domain)
python scale_trna.py

# 4. Phylogenetically controlled PGLS under three branch-length models
python scale_pgls.py
```

### Regenerating the data from scratch (optional; requires network)

The `data/` TSVs are cached so the analyses above run offline. To rebuild from
the primary sources:

```bash
python scale_trna.py --per-domain 150 --refetch   # GtRNAdb + NCBI -> data/trna_scale.tsv
python scale_pgls.py                                # resolves taxids, rebuilds the tree
```

## Data provenance

- **tRNA gene sets:** GtRNAdb (tRNAscan-SE results), https://gtrnadb.ucsc.edu
- **Genome sizes and taxonomy:** NCBI Assembly / Taxonomy (E-utilities)

## Citation

> Reich, J. (2026). *Two faces of one optimum: wobble error-resistance and minimal
> tRNA decoding coincide in the genetic code.* Preprint (bioRxiv).

## License

Code is released under the MIT License (see `LICENSE`). Derived data are credited
to GtRNAdb and NCBI under their respective terms.
