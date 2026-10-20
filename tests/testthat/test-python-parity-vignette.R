# Contract tests for the 1.0.0 vignette bundle (Phase 13 PR-B):
# inst/extdata/python_parity_vignette/ + the stablr-python-parity and
# stablr-reading-plots articles.
#
# Gates encoded here (bundle INTEGRATION.md):
#   1. every bundled asset is installed and readable;
#   2. each fixture matches its SHA-256 manifest entry;
#   3. each fixture is byte-identical to the canonical frozen fixture in
#      tests/testthat/fixtures/python_parity/ (copies, not substitutes);
#   4. both vignettes render without Python, using only installed-package
#      system.file() resources.

.parity_vignette_dir <- function() {
  d <- system.file("extdata", "python_parity_vignette", package = "stablr")
  if (!nzchar(d)) {
    # Source-tree fallback for devtools::test() without pkgload shims.
    d <- file.path("..", "..", "inst", "extdata", "python_parity_vignette")
  }
  normalizePath(d, mustWork = FALSE)
}

.parity_vignette_sha256 <- function(path) {
  sha256sum <- Sys.which("sha256sum")
  if (nzchar(sha256sum)) {
    output <- system2(sha256sum, path, stdout = TRUE, stderr = TRUE)
    return(tolower(strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]))
  }
  shasum <- Sys.which("shasum")
  if (nzchar(shasum)) {
    output <- system2(shasum, c("-a", "256", path), stdout = TRUE, stderr = TRUE)
    return(tolower(strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]))
  }
  NA_character_
}

.parity_vignette_files <- c(
  "fdp_plus_by_lambda.csv",
  "fdp_plus_curve.csv",
  "provenance.json",
  "solver_cases.csv",
  "stability_scores.csv",
  "stabl_fit_reference.csv"
)

test_that("python_parity_vignette assets are installed and readable", {
  d <- .parity_vignette_dir()
  expect_true(dir.exists(d))
  for (f in c(.parity_vignette_files, "manifest.sha256", "README.md")) {
    p <- file.path(d, f)
    expect_true(file.exists(p), info = f)
    expect_gt(file.size(p), 0L)
  }
})

test_that("python_parity_vignette fixtures match their SHA-256 manifest", {
  d <- .parity_vignette_dir()
  manifest <- readLines(file.path(d, "manifest.sha256"), warn = FALSE)
  manifest <- manifest[nzchar(trimws(manifest))]
  expect_length(manifest, length(.parity_vignette_files))
  for (line in manifest) {
    parts <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
    expect_gte(length(parts), 2L)
    expected <- tolower(parts[[1L]])
    fname <- sub("^\\*", "", parts[[2L]])
    observed <- .parity_vignette_sha256(file.path(d, fname))
    skip_if(is.na(observed), "no SHA-256 tool available")
    expect_identical(observed, expected, info = fname)
  }
})

test_that("python_parity_vignette fixtures are byte-identical to the canonical frozen fixtures", {
  d <- .parity_vignette_dir()
  canonical <- test_path("fixtures", "python_parity")
  skip_if_not(dir.exists(canonical), "canonical fixture tree not available")
  for (f in .parity_vignette_files) {
    a <- readBin(file.path(d, f), "raw", n = file.size(file.path(d, f)))
    b <- readBin(file.path(canonical, f), "raw",
                 n = file.size(file.path(canonical, f)))
    expect_identical(a, b, info = f)
  }
})

test_that("both bundle vignettes render without Python from installed resources", {
  skip_on_cran()
  skip_if_not_installed("knitr")
  skip_if_not_installed("rmarkdown")
  skip_if_not(dir.exists(.parity_vignette_dir()),
              "python_parity_vignette extdata not installed")

  vignette_src <- c(
    system.file("doc", "stablr-python-parity.Rmd", package = "stablr"),
    file.path("..", "..", "vignettes", "stablr-python-parity.Rmd")
  )
  vignette_src <- vignette_src[file.exists(vignette_src)]
  plots_src <- c(
    system.file("doc", "stablr-reading-plots.Rmd", package = "stablr"),
    file.path("..", "..", "vignettes", "stablr-reading-plots.Rmd")
  )
  plots_src <- plots_src[file.exists(plots_src)]
  expect_gt(length(vignette_src), 0L)
  expect_gt(length(plots_src), 0L)

  # Python must be absent and unnecessary: refuse to render if reticulate
  # is loaded, and render in a clean temp library path context.
  expect_false("reticulate" %in% loadedNamespaces(),
               info = "reticulate must not be needed for these articles")

  for (src in c(vignette_src[[1L]], plots_src[[1L]])) {
    out <- file.path(tempdir(), paste0(tools::file_path_sans_ext(basename(src)),
                                       ".html"))
    render_err <- tryCatch(
      {
        suppressWarnings(rmarkdown::render(
          src,
          output_format = "html_document",
          output_file = out,
          quiet = TRUE,
          envir = new.env(parent = globalenv())
        ))
        NULL
      },
      error = function(e) conditionMessage(e)
    )
    expect_null(render_err, info = basename(src))
    expect_true(file.exists(out))
    expect_gt(file.size(out), 10000L)
  }
})

# ── PBMC 3k single-cell vignette fixtures (stablr-cebmf / stablr-metacells) ──
# Same contract style as the python_parity_vignette bundle: assets installed,
# SHA-256 manifest honoured, structure consistent, and the vignettes' seeded
# analyses reproduce their canonical marker genes.

.pbmc3k_dir <- function() {
  d <- system.file("extdata", "pbmc3k", package = "stablr")
  if (!nzchar(d)) {
    # Source-tree fallback for devtools::test() without pkgload shims.
    d <- file.path("..", "..", "inst", "extdata", "pbmc3k")
  }
  normalizePath(d, mustWork = FALSE)
}

.pbmc3k_files <- c(
  "pbmc3k_counts_hvg.rds",
  "pbmc3k_clusters.rds",
  "pbmc3k_cebmf.rds",
  "run_metadata.json"
)

test_that("pbmc3k assets are installed and readable", {
  d <- .pbmc3k_dir()
  expect_true(dir.exists(d))
  for (f in c(.pbmc3k_files, "SHA256SUMS")) {
    p <- file.path(d, f)
    expect_true(file.exists(p), info = f)
    expect_gt(file.size(p), 0L)
  }
})

test_that("pbmc3k fixtures match their SHA-256 manifest", {
  d <- .pbmc3k_dir()
  manifest <- readLines(file.path(d, "SHA256SUMS"), warn = FALSE)
  manifest <- manifest[nzchar(trimws(manifest))]
  expect_length(manifest, length(.pbmc3k_files))
  for (line in manifest) {
    parts <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
    expect_gte(length(parts), 2L)
    expected <- tolower(parts[[1L]])
    fname <- sub("^\\*", "", parts[[2L]])
    observed <- .parity_vignette_sha256(file.path(d, fname))
    skip_if(is.na(observed), "no SHA-256 tool available")
    expect_identical(observed, expected, info = fname)
  }
})

test_that("pbmc3k fixtures are structurally consistent", {
  d <- .pbmc3k_dir()
  counts   <- readRDS(file.path(d, "pbmc3k_counts_hvg.rds"))
  clusters <- readRDS(file.path(d, "pbmc3k_clusters.rds"))
  cebmf    <- readRDS(file.path(d, "pbmc3k_cebmf.rds"))

  expect_s4_class(counts, "dgCMatrix")
  expect_equal(dim(counts), c(2700L, 1200L))
  expect_true(all(counts@x >= 0))                 # raw counts
  expect_true(all(counts@x == floor(counts@x)))   # integer-valued

  expect_s3_class(clusters, "factor")
  expect_length(clusters, 2700L)
  expect_identical(names(clusters), rownames(counts))

  expect_equal(dim(cebmf$L), c(2700L, 10L))
  expect_identical(rownames(cebmf$L), rownames(counts))
  expect_equal(dim(cebmf$F), c(1200L, 10L))
  expect_identical(rownames(cebmf$F), colnames(counts))
  expect_false(anyNA(cebmf$L))
  expect_false(anyNA(cebmf$F))
})

test_that("the cebmf vignette selection recovers the B-cell receptor module", {
  skip_on_cran()
  d <- .pbmc3k_dir()
  counts   <- readRDS(file.path(d, "pbmc3k_counts_hvg.rds"))
  clusters <- readRDS(file.path(d, "pbmc3k_clusters.rds"))
  cebmf    <- readRDS(file.path(d, "pbmc3k_cebmf.rds"))

  # Exact vignette pipeline (stablr-cebmf.Rmd sections 3-5).
  keep <- clusters != "c9"
  clusters <- droplevels(clusters[keep])
  counts   <- counts[keep, ]
  L        <- cebmf$L[keep, ]
  y <- factor(ifelse(clusters == "c2", "B", "other"))
  names(y) <- rownames(counts)
  lib <- Matrix::rowSums(counts)
  xn  <- log1p(1e4 * counts / lib)
  v   <- Matrix::colSums(xn^2) / nrow(xn) - (Matrix::colSums(xn) / nrow(xn))^2
  top <- order(v, decreasing = TRUE)[1:500]
  denoised <- L %*% t(cebmf$F)
  colnames(denoised) <- rownames(cebmf$F)

  fit <- stabl_fit(
    x = denoised[, top], y = y, family = "binomial",
    lambda_grid = "auto", n_bootstraps = 100L, n_lambda = 20L,
    random_state = 42L
  )
  sel <- get_feature_names(fit)
  expect_gt(length(sel), 0L)
  expect_lte(length(sel), 10L)
  expect_true(all(c("CD79A", "CD79B") %in% sel))
})

test_that("the metacells vignette selection recovers monocyte markers", {
  skip_on_cran()
  d <- .pbmc3k_dir()
  counts   <- readRDS(file.path(d, "pbmc3k_counts_hvg.rds"))
  clusters <- readRDS(file.path(d, "pbmc3k_clusters.rds"))

  # Exact vignette pipeline (stablr-metacells.Rmd sections 2-5).
  keep <- clusters != "c9"
  counts   <- counts[keep, ]
  clusters <- droplevels(clusters[keep])
  lib <- Matrix::rowSums(counts)
  xn  <- log1p(1e4 * counts / lib)
  pc <- prcomp(as.matrix(xn), rank. = 30)$x
  dd <- as.matrix(dist(pc))
  n  <- nrow(pc)
  k  <- 25L
  assigned <- rep(NA_integer_, n)
  mc <- 0L
  for (i in order(pc[, 1])) {
    if (!is.na(assigned[i])) next
    pool <- which(is.na(assigned))
    if (length(pool) < 2L * k) break
    mc <- mc + 1L
    assigned[pool[order(dd[i, pool])[1:k]]] <- mc
  }
  left <- which(is.na(assigned))
  for (i in left) {
    idx <- which(!is.na(assigned))
    assigned[i] <- assigned[idx[which.min(dd[i, idx])]]
  }
  expect_false(anyNA(assigned))
  expect_equal(length(unique(assigned)), 106L)   # sizes 25-32, all 2689 cells assigned
  mc_counts <- as.matrix(Matrix::fac2sparse(factor(assigned)) %*% counts)
  mc_xn     <- log1p(1e4 * mc_counts / rowSums(mc_counts))
  maj <- tapply(clusters, assigned, function(v) {
    tb <- table(v)
    c(names(which.max(tb)), max(tb) / length(v))
  })
  mc_lab <- factor(vapply(maj, `[`, character(1L), 1L))
  purity <- as.numeric(vapply(maj, `[`, character(1L), 2L))
  keepmc <- purity >= 0.6
  expect_equal(sum(!keepmc), 15L)                # purity QC count
  y <- factor(ifelse(mc_lab[keepmc] == "c4", "monocyte", "other"))
  names(y) <- rownames(mc_xn)[keepmc]

  fit <- stabl_fit(
    x = mc_xn[keepmc, ], y = y, family = "binomial",
    lambda_grid = "auto", n_bootstraps = 200L, n_lambda = 20L,
    stratify_bootstrap = TRUE, bootstrap_strata = y,
    random_state = 42L
  )
  sel <- get_feature_names(fit)
  expect_gt(length(sel), 0L)
  expect_true(all(c("S100A8", "CD14") %in% sel))
})
