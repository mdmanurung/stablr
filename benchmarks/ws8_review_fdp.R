# WS-8 review follow-up 1: per-environment FDP+ diagnostics.
#
# Closes predesign success criterion 2 ("FDP+ diagnostics are computable for
# every per-environment fit").  The diagnostics were computed inside every
# stabl_fit at run time (field `min_fdr_`) but never extracted or reported.
# This script reads the objects saved by ws8_aging_showcase.R (no refitting)
# and writes ws8_fdp_diagnostics.csv.
#
# Usage: run from the repo root, e.g.
#   Rscript benchmarks/ws8_review_fdp.R
#
# Application case study — NOT calibration evidence.

pkgload::load_all(".", compile = FALSE, quiet = TRUE)

obj_path <- Sys.getenv("WS8_OBJECTS", "/workspace/ws8_out/ws8_objects.rds")
out_dir  <- Sys.getenv("WS8_OUT_DIR", "benchmarks")
obj  <- readRDS(obj_path)
fits <- obj$env_rf$env_fits

fdp <- do.call(rbind, lapply(names(fits), function(e) {
  f   <- fits[[e]]
  tau <- f$fdr_threshold_range
  i05 <- match(0.5, tau)
  stopifnot(!is.na(i05))
  q_real <- apply(f$stabl_scores_, 1, max)           # max-over-lambda, real
  q_art  <- apply(f$stabl_scores_artificial_, 1, max) # max-over-lambda, decoys
  data.frame(
    cohort              = e,
    n_donors            = sum(obj$env == e),
    min_fdp_plus        = f$min_fdr_,
    tau_at_min_fdp      = f$fdr_min_threshold_,
    fdp_plus_at_0.5     = unname(f$FDRs_[i05]),
    n_real_score_ge_0.5 = sum(q_real >= 0.5),
    n_art_score_ge_0.5  = sum(q_art >= 0.5),
    stringsAsFactors    = FALSE
  )
}))

write.csv(fdp, file.path(out_dir, "ws8_fdp_diagnostics.csv"), row.names = FALSE)
cat("Per-environment FDP+ diagnostics (rf + knockoff_equi, 200 bootstraps):\n")
print(fdp, row.names = FALSE)
