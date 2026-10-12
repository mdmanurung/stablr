#!/usr/bin/env Rscript
# =====================================================================
# Generator for the multiview-parity fixtures (improvement item F4).
#
# Freezes the vendored cooperative backend (R/cooperative-*.R, a fork of
# CRAN multiview 1.0) against reference values produced by CRAN
# multiview 1.0 itself on fixed input data, and pins the three deliberate
# divergences documented in inst/PATCHES.cooperative.md:
#   D4/B1  offset subsetting in build_predmat
#   D5/B2  halved-step warm-start residual weighting
#   D6/B3  lower.limits plumbing to elnet.fit
#
# Generation-time dependencies (NOT needed to run the package tests):
#   * CRAN multiview 1.0  -- reference engine
#   * pkgload             -- loads the vendored backend incl. compiled C++
#   * jsonlite            -- provenance record
# Tests consume only the frozen files listed in manifest.sha256.
#
# Run from the package root:
#   R_LIBS=<lib> Rscript tests/testthat/fixtures/multiview_parity/generate.R
# =====================================================================

pkg_root <- getwd()
stopifnot(
  file.exists(file.path(pkg_root, "DESCRIPTION")),
  dir.exists(file.path(pkg_root, "tests", "testthat"))
)

suppressPackageStartupMessages({
  library(glmnet)
  library(multiview)
})
stopifnot(as.character(packageVersion("multiview")) == "1.0")
stopifnot(requireNamespace("pkgload", quietly = TRUE),
          requireNamespace("jsonlite", quietly = TRUE))

out_dir <- file.path(pkg_root, "tests", "testthat", "fixtures", "multiview_parity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load the vendored backend with its compiled C++ routines (wls_exp etc.);
# sourcing the R files alone is not enough for the weighted/offset paths.
pkgload::load_all(pkg_root, quiet = TRUE)

## ------------------------------------------------------------------
## 1. Fixed inputs (seed 2026, n = 80, two views of p = 6 features)
## ------------------------------------------------------------------
set.seed(2026)
n <- 80L
p <- 6L
view_a <- matrix(rnorm(n * p), nrow = n, dimnames = list(NULL, paste0("a", seq_len(p))))
view_b <- matrix(rnorm(n * p), nrow = n, dimnames = list(NULL, paste0("b", seq_len(p))))
y_gaussian <- drop(1.2 * view_a[, 1] - 0.9 * view_b[, 2] + 0.5 * view_a[, 3] + rnorm(n, sd = 0.5))
eta_binomial <- 0.8 * view_a[, 1] - 0.6 * view_b[, 2] + 0.4 * view_a[, 5]
y_binomial <- rbinom(n, 1, stats::plogis(eta_binomial))
weights_unit <- rep(1, n)
weights_weighted <- rep(1, n)
weights_weighted[1:10] <- 20
weights_weighted[71:80] <- 0.25
offset_zero <- rep(0, n)
offset_nonzero <- round(0.3 * cos(seq_len(n) / 8), 6)
# aggressive offset aligned with the outcome: the only regime where the
# D5/B2 halved-step residual divergence can bite (needs halving to trigger)
offset_separating <- round(4 * (2 * y_binomial - 1), 6)
foldid <- sample(rep(1:5, length.out = n))
lambda_gaussian <- c(0.08, 0.04, 0.02, 0.01)
lambda_binomial <- c(0.05, 0.02, 0.01, 0.005)
rho <- 0.3
x_list <- list(a = view_a, b = view_b)

inputs <- list(
  seed = 2026L, n = n, p = p, rho = rho,
  view_a = view_a, view_b = view_b,
  y_gaussian = y_gaussian, y_binomial = y_binomial,
  weights_unit = weights_unit, weights_weighted = weights_weighted,
  offset_zero = offset_zero, offset_nonzero = offset_nonzero,
  offset_separating = offset_separating,
  foldid = foldid,
  lambda_gaussian = lambda_gaussian, lambda_binomial = lambda_binomial
)
saveRDS(inputs, file.path(out_dir, "inputs.rds"))

## ------------------------------------------------------------------
## 2. Engine runners (vendored = stablr fork; cran = multiview 1.0)
## ------------------------------------------------------------------
.muffled <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(
    expr,
    warning = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleWarning")
    }
  )
  list(value = val, warnings = unique(msgs))
}

run_vendored <- function(y, family, weights = NULL, offset = NULL, lambda,
                         lower.limits = -Inf) {
  .muffled(stablr:::.cooperative_backend_fit(
    x_list, y, family = family, rho = rho, weights = weights,
    offset = offset, lambda = lambda, lower.limits = lower.limits
  ))
}
run_cran <- function(y, family, weights = NULL, offset = NULL, lambda,
                     lower.limits = -Inf) {
  .muffled(multiview::multiview(
    x_list, y, family = family, rho = rho, weights = weights,
    offset = offset, lambda = lambda, lower.limits = lower.limits
  ))
}

coef_frame <- function(fit, lambda, family_label, engine_label) {
  b <- as.matrix(fit$beta)
  stopifnot(nrow(b) == 2L * p, ncol(b) == length(lambda))
  out <- expand.grid(feature = rownames(b), lambda_index = seq_along(lambda),
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  out <- out[order(out$lambda_index, out$feature), ]
  fidx <- match(out$feature, rownames(b))
  out$view <- substr(out$feature, 1L, 1L)
  out$lambda <- lambda[out$lambda_index]
  out$coefficient <- b[cbind(fidx, out$lambda_index)]
  out$selected <- out$coefficient != 0
  out$family <- family_label
  out$engine <- engine_label
  out[, c("family", "engine", "view", "feature", "lambda_index", "lambda",
          "coefficient", "selected")]
}

engine_diff <- function(df) {
  v <- df[df$engine == "vendored", , drop = FALSE]
  k <- df[df$engine == "cran_1_0", , drop = FALSE]
  key_v <- paste(v$family, v$feature, v$lambda_index)
  key_k <- paste(k$family, k$feature, k$lambda_index)
  stopifnot(setequal(key_v, key_k))
  k2 <- k[match(key_v, key_k), , drop = FALSE]
  max(abs(v$coefficient - k2$coefficient))
}

parity_case <- function(y, family_key, weights, offset, lambda, family_label) {
  fam_obj <- switch(family_key, gaussian = stats::gaussian(), binomial = stats::binomial())
  v <- run_vendored(y, family_key, weights, offset, lambda)$value
  k <- run_cran(y, fam_obj, weights, offset, lambda)$value
  rbind(coef_frame(v, lambda, family_label, "vendored"),
        coef_frame(k, lambda, family_label, "cran_1_0"))
}

write_fixture_csv <- function(df, name) {
  write.csv(df, file.path(out_dir, name), row.names = FALSE)
}

## 2a. Parity: unit weights, no offset -- paths where the fork must match
##     CRAN multiview 1.0 exactly.
write_fixture_csv(
  parity_case(y_gaussian, "gaussian", NULL, NULL, lambda_gaussian, "gaussian"),
  "reference_parity_gaussian.csv")
write_fixture_csv(
  parity_case(y_binomial, "binomial", NULL, NULL, lambda_binomial, "binomial"),
  "reference_parity_binomial.csv")

## 2b. Parity under non-unit case weights, still no offset (D5/B2 must not
##     bite: with a zero offset both residual expressions agree).
write_fixture_csv(
  rbind(
    parity_case(y_gaussian, "gaussian", weights_weighted, NULL,
                lambda_gaussian, "gaussian_weighted"),
    parity_case(y_binomial, "binomial", weights_weighted, NULL,
                lambda_binomial, "binomial_weighted")),
  "reference_parity_weighted.csv")

## 2c. D5/B2 fixture: weights AND a non-zero offset.
b2 <- rbind(
  parity_case(y_gaussian, "gaussian", weights_weighted, offset_nonzero,
              lambda_gaussian, "gaussian_weighted_offset"),
  parity_case(y_binomial, "binomial", weights_weighted, offset_nonzero,
              lambda_binomial, "binomial_weighted_offset"))
b2_halved <- NULL
b2_halved_info <- list(completed = FALSE,
                       note = paste("aggressive separating-offset case did not",
                                    "complete for both engines; only the mild",
                                    "offset cases are frozen"))
agg_v <- try(run_vendored(y_binomial, "binomial", weights_weighted,
                          offset_separating, lambda_binomial), silent = TRUE)
agg_k <- try(run_cran(y_binomial, stats::binomial(), weights_weighted,
                      offset_separating, lambda_binomial), silent = TRUE)
if (!inherits(agg_v, "try-error") && !inherits(agg_k, "try-error")) {
  b2_halved <- rbind(
    coef_frame(agg_v$value, lambda_binomial, "binomial_weighted_offset_separating", "vendored"),
    coef_frame(agg_k$value, lambda_binomial, "binomial_weighted_offset_separating", "cran_1_0"))
  b2_halved_info <- list(
    completed = TRUE,
    halving_warnings_vendored = paste(agg_v$warnings, collapse = " | "),
    halving_warnings_cran = paste(agg_k$warnings, collapse = " | "))
}
b2_all <- rbind(b2, b2_halved)
write_fixture_csv(b2_all, "reference_b2_weighted_offset.csv")

## 2d. D6/B3 fixture: binding lower.limits (feature b2 would naturally sit
##     near -0.76 at the densest lambda).
lower_limits <- rep(-Inf, 2L * p)
lower_limits[p + 2L] <- -0.2  # feature "b2"
b3_v <- run_vendored(y_gaussian, "gaussian", NULL, NULL,
                     lambda_gaussian, lower.limits = lower_limits)$value
b3_k <- run_cran(y_gaussian, stats::gaussian(), NULL, NULL,
                 lambda_gaussian, lower.limits = lower_limits)$value
b3 <- rbind(
  coef_frame(b3_v, lambda_gaussian, "gaussian_lower_limits", "vendored"),
  coef_frame(b3_k, lambda_gaussian, "gaussian_lower_limits", "cran_1_0"))
write_fixture_csv(b3, "reference_b3_lower_limits.csv")

## 2e. D4/B1 fixture: build_predmat with a non-NULL offset. Upstream
##     multiview 1.0 self-references `offset_sub` before it exists
##     (R/build_predmat.R:26) and errors; the vendored fix subsets rows.
outlist <- vector("list", max(foldid))
for (k in seq_len(max(foldid))) {
  which_k <- foldid == k
  x_sub <- lapply(x_list, function(m) m[!which_k, , drop = FALSE])
  outlist[[k]] <- .muffled(stablr:::.mv_multiview(
    x_list = x_sub, y = y_gaussian[!which_k], family = stats::gaussian(),
    rho = rho, offset = offset_nonzero[!which_k], lambda = lambda_gaussian
  ))$value
}
predmat_v <- stablr:::build_predmat(outlist, lambda_gaussian, x_list,
                                    offset_nonzero, foldid, "lambda",
                                    stats::gaussian())
cran_b1_behavior <- tryCatch({
  multiview:::build_predmat(outlist, lambda_gaussian, x_list, offset_nonzero,
                            foldid, "lambda", stats::gaussian())
  "completed without error"
}, error = function(e) paste0("error: ", conditionMessage(e)))

b1 <- expand.grid(sample_index = seq_len(n), lambda_index = seq_along(lambda_gaussian),
                  KEEP.OUT.ATTRS = FALSE)
b1 <- b1[order(b1$sample_index, b1$lambda_index), ]
b1$prediction <- predmat_v[cbind(b1$sample_index, b1$lambda_index)]
b1$lambda <- lambda_gaussian[b1$lambda_index]
write_fixture_csv(b1[, c("sample_index", "lambda_index", "lambda", "prediction")],
                  "reference_b1_predmat_offset.csv")

## 2f. D4/B1 end-to-end: the cv path with a non-NULL offset (pre-fix this
##     errors inside build_predmat).
cv_v <- .muffled(stablr:::.cooperative_backend_cv(
  x_list, y_gaussian, family = "gaussian", rho = rho, foldid = foldid,
  offset = offset_nonzero, lambda = lambda_gaussian
))$value
cv_fields <- c("lambda", "cvm", "cvsd", "cvlo", "cvup")
if (!inherits(cv_v, "try-error") && all(cv_fields %in% names(cv_v))) {
  write.csv(data.frame(lambda = cv_v$lambda, cvm = cv_v$cvm, cvsd = cv_v$cvsd,
                       cvlo = cv_v$cvlo, cvup = cv_v$cvup),
            file.path(out_dir, "reference_b1_cv_offset.csv"), row.names = FALSE)
  cv_fixture_ok <- TRUE
} else {
  cv_fixture_ok <- FALSE
}

## ------------------------------------------------------------------
## 3. Provenance record (machine-readable, feeds the B4 divergence log)
## ------------------------------------------------------------------
prov <- list(
  fixture_set = "multiview_parity (improvement item F4)",
  generated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  reference_engine = list(package = "multiview",
                          version = as.character(packageVersion("multiview")),
                          source = "CRAN"),
  vendored_engine = list(package = "stablr", files = "R/cooperative-*.R",
                         patch_log = "inst/PATCHES.cooperative.md"),
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  glmnet_version = as.character(packageVersion("glmnet")),
  seed = 2026L,
  inputs = "inputs.rds",
  design = list(n = n, p = p, views = c("a", "b"), rho = rho,
                families = c("gaussian", "binomial"),
                weights = c("unit", "weighted"),
                offsets = c("zero", "nonzero", "separating")),
  parity_cases = c("reference_parity_gaussian.csv",
                   "reference_parity_binomial.csv",
                   "reference_parity_weighted.csv"),
  divergence_cases = list(
    D4_B1 = list(
      patch = "offset subsetting in build_predmat",
      fixtures = c("reference_b1_predmat_offset.csv",
                   if (cv_fixture_ok) "reference_b1_cv_offset.csv" else NULL),
      upstream_1_0_behavior = cran_b1_behavior),
    D5_B2 = list(
      patch = "halved-step warm-start residual weighting",
      fixture = "reference_b2_weighted_offset.csv",
      max_abs_coef_diff_vendored_vs_cran =
        unname(round(engine_diff(b2_all), 12)),
      aggressive_offset_case = b2_halved_info),
    D6_B3 = list(
      patch = "lower.limits plumbing to elnet.fit",
      fixture = "reference_b3_lower_limits.csv",
      max_abs_coef_diff_vendored_vs_cran =
        unname(round(engine_diff(b3), 12)))
  )
)
writeLines(jsonlite::toJSON(prov, auto_unbox = TRUE, pretty = TRUE, digits = NA),
           file.path(out_dir, "provenance.json"))

## ------------------------------------------------------------------
## 4. SHA-256 manifest (same pattern as fixtures/python_parity)
## ------------------------------------------------------------------
sha256_of <- function(path) {
  tool <- Sys.which("sha256sum")
  if (nzchar(tool)) {
    return(tolower(strsplit(system2(tool, path, stdout = TRUE)[[1L]],
                            "[[:space:]]+")[[1L]][[1L]]))
  }
  tool <- Sys.which("shasum")
  stopifnot(nzchar(tool))
  tolower(strsplit(system2(tool, c("-a", "256", path), stdout = TRUE)[[1L]],
                   "[[:space:]]+")[[1L]][[1L]])
}
files <- sort(setdiff(list.files(out_dir), "manifest.sha256"))
manifest <- vapply(files, function(f) {
  sprintf("%s  %s", sha256_of(file.path(out_dir, f)), f)
}, character(1))
writeLines(manifest, file.path(out_dir, "manifest.sha256"))

## ------------------------------------------------------------------
## 5. Diagnostics
## ------------------------------------------------------------------
cat("Wrote", length(files), "fixture files + manifest.sha256\n")
cat("parity gaussian  max|vendored-cran|:", engine_diff(read.csv(file.path(out_dir, "reference_parity_gaussian.csv"))), "\n")
cat("parity binomial  max|vendored-cran|:", engine_diff(read.csv(file.path(out_dir, "reference_parity_binomial.csv"))), "\n")
cat("parity weighted  max|vendored-cran|:", engine_diff(read.csv(file.path(out_dir, "reference_parity_weighted.csv"))), "\n")
cat("B2 max|vendored-cran|:", engine_diff(b2_all), "\n")
cat("B2 aggressive case completed:", b2_halved_info$completed, "\n")
cat("B3 max|vendored-cran|:", engine_diff(b3), "\n")
cat("B1 upstream 1.0 behavior:", cran_b1_behavior, "\n")
cat("B1 cv fixture written:", cv_fixture_ok, "\n")
