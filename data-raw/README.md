# Bundled Onset of Labor data provenance

The bounded examples in `inst/extdata` derive from the Onset of Labor files in
the Dryad dataset “Discovery of sparse, reliable omic biomarkers with Stabl”
(<https://doi.org/10.5061/dryad.stqjq2c7d>), released under CC0. The Onset of
Labor source study is <https://doi.org/10.1126/scitranslmed.abd9898>.

Dryad publishes these files in `data.zip` (SHA-256
`302dfd9919d416f0afb394d1560dad44180894210389bc9bd91d47f501e34ddf`).
The identical archive is mirrored in `Sample Data/data.zip` at the pinned
public reference commit
`gregbellan/Stabl@1d07f85a13cfbecb4f08ce21075bf4fbb8e34678`.

Only the sample identifier, DOS outcome, and the first 100 CyTOF and proteomic
feature columns (in source order) are retained. The six source tables already
encode the study's training and validation split; no new split is generated.
Values are the study-exported processed measurements; the
package applies no additional normalization or imputation.

Run `prepare_ool_data.R` on the six exported Dryad tables to rebuild the
fixtures. The script verifies these source SHA-256 checksums before writing:

- Training CyTOF: `c17d4d2ad44a2a8b1da887a58b0cf190c1f74c23755f0f80c6231945f97f47f6`
- Training proteomics: `5037b8e2ff3c28c4d562bb32e6d04977c1532f73d659bdccaefa18271785dfba`
- Training DOS: `11889faf2d1c3846025f57fa89e3a92f036aeff7162596a4d69f7eccf72695ec`
- Validation CyTOF: `912053964c201c215ce557ab9bc43488a74da3f9c1b2b331f5ec5b85de40ccfd`
- Validation proteomics: `55fb10bd1d2093d8bc593e7d4509fb51389f8882b842db1e8bf8139ff6584c71`
- Validation DOS: `f1bf4ed0dc1f577fb59e44b4a9014bcc4b7516677705c524790cadfc226c8ab6`

The committed output-file SHA-256 checksums are:

- `ool_cytof_train.csv.gz`: `2fffa0d6bdfeea18886e7b6036fb2c57ee6edd21220112c673182c0f44db4b97`
- `ool_cytof_valid.csv.gz`: `e8d1ffd721f9f41994bd98b78cca9431f60ab25ca59e5f2b5fc8577f9b313fa3`
- `ool_dos_train.csv`: `fcdf4f931cbb40ab084da86da5b48f454e43cca945f070e301ee101c438dbe70`
- `ool_dos_valid.csv`: `90223bfdde2233108ea7b046a0c8c9a8228afbcb56005e8a32dc0640472ce530`
- `ool_proteomics_train.csv.gz`: `d391b512f984e4d3c0141e3b53fcaefc82a3b006a59549c09114a25752e34c4c`
- `ool_proteomics_valid.csv.gz`: `8072550a5e0e415bb235e9e705ade8f3033249a3a62fb550af0fb3e256ff8a41`

# Bundled PBMC 3k single-cell data provenance

The fixtures in `inst/extdata/pbmc3k` derive from the 10x Genomics PBMC 3k
filtered gene-barcode matrix (hg19), downloaded from
<https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz>
(md5 of the tarball recorded in `run_metadata.json` as
`f76d73bf7a582aaf32b798bcae9574d3`). The dataset is the classic ~2,700 PBMC
sample from a healthy donor, distributed by 10x Genomics and originally
profiled in Zheng et al. (2017) <https://doi.org/10.1038/ncomms14049>.

Processing (all seeded and deterministic; see `run_metadata.json` for the full
parameter and version record): genes detected in fewer than 10 cells are
removed; counts are CP10k-normalized and log1p-transformed; the 1,200 most
variable genes are retained; cells are embedded in 50 principal components;
an exact kNN graph (k = 10) feeds Louvain clustering
(`igraph::community_multilevel`, resolution 1.0, seed 20261020), giving the
nine `c1`-`c9` labels; the log-normalized HVG matrix is factorized with the
Python package `cebmf_torch` (`cEBMF`, K = 10, normal priors on both factors,
50 iterations, CPU), yielding the cell loading matrix `L` (2,700 x 10) and
gene factor matrix `F` (1,200 x 10).

Three fixtures are shipped: `pbmc3k_counts_hvg.rds` (raw counts of the 1,200
HVGs as a dgCMatrix; normalization is repeated live in the vignettes),
`pbmc3k_clusters.rds` (named factor of Louvain labels), and
`pbmc3k_cebmf.rds` (the cEBMF `L` and `F` matrices). The cEBMF fit is
precomputed so the vignettes and the R CMD check gate stay hermetic (no
Python/torch dependency at build time); the vignette shows the Python code in
unevaluated chunks.

Run `prepare_pbmc3k_fixtures.py` (downloads and factors; requires `torch`,
`igraph`, `cebmf_torch`) followed by `prepare_pbmc3k_fixtures.R` (assembles
the .rds fixtures and the SHA-256 manifest) to rebuild the fixtures. The
committed output-file SHA-256 checksums (also in
`inst/extdata/pbmc3k/SHA256SUMS`) are:

- `pbmc3k_counts_hvg.rds`: `cb2a140940408712d8e23f0054c89d7f011f8360fb0164452fb3e9a7da659653`
- `pbmc3k_clusters.rds`: `a79b3639d9ab35a12a0cb52bf7e0cad1317b6ce06b1d18f556899fffdd4db255`
- `pbmc3k_cebmf.rds`: `6d168823e9e1689dc5adf0ab46aa137e5fa5d7e24f8cc80e2e28b10059ea0e76`
- `run_metadata.json`: `a63e5d7cceee4b13dfeb57777cda6d90ebe3f4f99e7feaf40b49789429924edf`
