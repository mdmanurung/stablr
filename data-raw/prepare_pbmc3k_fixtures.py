"""Generate PBMC 3k fixtures for the stablr cebmf/metacells vignettes.

Provenance: 10x Genomics PBMC 3k filtered gene-barcode matrix (hg19),
https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

Steps: gene filter -> CP10k + log1p -> 1,200 HVGs -> PCA-50 -> kNN(k=10)
-> Louvain clusters -> cEBMF (cebmf_torch, K=10, CPU, seeded).

Outputs (in /workspace/pbmc3k/out/): counts_hvg.mtx (cells x genes),
genes_hvg.txt, barcodes.txt, clusters.csv, cebmf_L.csv, cebmf_F.csv,
run_metadata.json.
"""
import json
import subprocess
import tarfile
from pathlib import Path

import igraph as ig
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import torch

SEED = 20261020
BASE = Path("/workspace/pbmc3k")
OUT = BASE / "out"
OUT.mkdir(exist_ok=True)

# ── 1. extract + load ────────────────────────────────────────────────────────
tgz = BASE / "pbmc3k.tar.gz"
if not (BASE / "filtered_gene_bc_matrices").exists():
    with tarfile.open(tgz) as tf:
        tf.extractall(BASE)
mat_dir = BASE / "filtered_gene_bc_matrices" / "hg19"
X = scipy.io.mmread(mat_dir / "matrix.mtx").tocsr().T  # cells x genes
genes = pd.read_csv(mat_dir / "genes.tsv", sep="\t", header=None,
                    names=["ensembl", "symbol"])
barcodes = pd.read_csv(mat_dir / "barcodes.tsv", header=None, names=["barcode"])
print("loaded:", X.shape, "cells x genes")

# ── 2. gene filter + normalize + HVG ─────────────────────────────────────────
min_cells = 10
keep = np.asarray((X > 0).sum(axis=0)).ravel() >= min_cells
X = X[:, keep]
genes = genes.loc[keep].reset_index(drop=True)
print("after filter:", X.shape)

lib = np.asarray(X.sum(axis=1)).ravel()
Xn = X.multiply(1e4 / lib[:, None]).tocsr()
Xn.data = np.log1p(Xn.data)                      # log1p CP10k, sparse

rng = np.random.default_rng(SEED)
# HVG by variance of log-normalized values (dense for 2.6k x 3.3k is fine)
Xd = Xn.toarray()
var = Xd.var(axis=0)
n_hvg = 1200
hvg_idx = np.argsort(var)[::-1][:n_hvg]
hvg_idx.sort()
Xh = Xd[:, hvg_idx]
genes_hvg = genes.loc[hvg_idx].reset_index(drop=True)
# unique-ify symbols for R dimnames
symbols = genes_hvg["symbol"].astype(str)
dup = symbols.duplicated()
symbols[dup] = symbols[dup] + "." + genes_hvg.loc[dup, "ensembl"].astype(str)
print("HVG:", Xh.shape)

# ── 3. PCA + kNN + Louvain clusters ──────────────────────────────────────────
Xc = Xh - Xh.mean(axis=0)
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
pcs = U[:, :50] * S[:50]
# exact kNN (k=10) on PCs
d2 = ((pcs[:, None, :] - pcs[None, :, :]) ** 2).sum(-1)
np.fill_diagonal(d2, np.inf)
knn = np.argsort(d2, axis=1)[:, :10]
edges = [(i, j) for i in range(len(pcs)) for j in knn[i]]
g = ig.Graph(edges=edges, n=len(pcs), directed=False)
g.simplify()
# deterministic Louvain: route igraph's C RNG through Python's seeded RNG
import random
random.seed(SEED)
ig.set_random_number_generator(random)
cl = g.community_multilevel(weights=None, resolution=1.0)
labels = np.array(cl.membership) + 1  # 1-based
print("clusters:", pd.Series(labels).value_counts().sort_index().to_dict())

# ── 4. cEBMF ─────────────────────────────────────────────────────────────────
from cebmf_torch import cEBMF
torch.manual_seed(SEED)
Y = torch.tensor(Xh, dtype=torch.float32)
model = cEBMF(Y, K=10, prior_L="norm", prior_F="norm")
fit = model.fit(maxit=50)
L = fit.L.detach().cpu().numpy()
F = fit.F.detach().cpu().numpy()
print("cEBMF: L", L.shape, "F", F.shape, "tau", float(fit.tau))

# ── 5. write outputs ─────────────────────────────────────────────────────────
# Ship RAW counts for the HVG set (the metacells vignette normalizes and
# aggregates live; raw integers are also smaller on disk than normalized
# doubles).
scipy.io.mmwrite(OUT / "counts_hvg_raw.mtx", sp.csr_matrix(X[:, hvg_idx]))
barcodes.to_csv(OUT / "barcodes.txt", index=False, header=False)
genes_hvg_out = genes_hvg.copy()
genes_hvg_out["symbol_unique"] = symbols.values
genes_hvg_out.to_csv(OUT / "genes_hvg.txt", index=False, header=False)
pd.DataFrame({"barcode": barcodes["barcode"], "cluster": labels}).to_csv(
    OUT / "clusters.csv", index=False)
pd.DataFrame(L, columns=[f"cEBMF{k}" for k in range(1, 11)]).to_csv(
    OUT / "cebmf_L.csv", index=False)
Fdf = pd.DataFrame(F, columns=[f"cEBMF{k}" for k in range(1, 11)])
Fdf.insert(0, "symbol", symbols.values)
Fdf.to_csv(OUT / "cebmf_F.csv", index=False)

meta = {
    "seed": SEED,
    "source_url": "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
    "source_md5": subprocess.check_output(["md5sum", str(tgz)], text=True).split()[0],
    "n_cells": int(X.shape[0]), "n_genes_hvg": int(Xh.shape[1]),
    "min_cells_gene_filter": min_cells,
    "normalization": "CP10k + log1p",
    "pca_components": 50, "knn_k": 10, "cluster_method": "Louvain (igraph community_multilevel, resolution 1.0)",
    "cebmf": {"package": "cebmf_torch", "K": 10, "prior_L": "norm",
              "prior_F": "norm", "maxit": 50, "device": "cpu"},
    "versions": {"numpy": np.__version__, "torch": torch.__version__,
                 "igraph": ig.__version__, "pandas": pd.__version__},
}
(OUT / "run_metadata.json").write_text(json.dumps(meta, indent=2))
print("fixtures written to", OUT)
