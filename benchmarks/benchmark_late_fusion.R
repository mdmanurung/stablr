#!/usr/bin/env Rscript

# Reprofile the same late-fusion workload in legacy and five-fold OOF modes.
# Raw profiles belong outside the repository, for example:
# Rscript bench/benchmark_late_fusion.R /tmp/stablr-v0.1.1-profile 5

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1L]] else "/tmp/stablr-v0.1.1-profile"
repetitions <- if (length(args) >= 2L) as.integer(args[[2L]]) else 5L
if (is.na(repetitions) || repetitions < 3L) stop("repetitions must be >= 3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required")
pkgload::load_all(".", quiet = TRUE)

set.seed(220722L)
n <- 60L; p <- 10L
ids <- paste0("s", seq_len(n))
x1 <- matrix(rnorm(n * p), n, p,
             dimnames = list(ids, paste0("a", seq_len(p))))
x2 <- matrix(rnorm(n * p), n, p,
             dimnames = list(ids, paste0("b", seq_len(p))))
y <- setNames(x1[, 1L] - 0.7 * x2[, 1L] + rnorm(n), ids)
common <- list(
  x_train_list = list(a = x1, b = x2), y_train = y,
  lambda_grid = data.frame(lambda = c(0.2, 0.1, 0.05)),
  artificial_type = "random_permutation", n_bootstraps = 10L,
  late_fusion = TRUE, late_fusion_nfolds = 5L, n_iter_lf = 100L,
  random_state = 220722L, workers = 1L
)

run_one <- function(mode, profile = FALSE) {
  if (profile) {
    Rprof(file.path(out_dir, paste0("late_fusion_", mode, ".Rprof")),
          interval = 0.01, memory.profiling = TRUE)
    on.exit(Rprof(NULL), add = TRUE)
    Rprofmem(file.path(out_dir, paste0("late_fusion_", mode, ".Rprofmem")))
    on.exit(Rprofmem(NULL), add = TRUE)
  }
  invisible(gc(reset = TRUE))
  elapsed <- system.time(do.call(
    stabl_multiomic_train_validate,
    c(common, list(late_fusion_training = mode))
  ))[["elapsed"]]
  memory <- gc()
  max_bytes <- memory["Ncells", "max used"] * 56 +
    memory["Vcells", "max used"] * 8
  c(elapsed_sec = unname(elapsed), gc_max_bytes = unname(max_bytes))
}

rows <- list(); k <- 1L
for (mode in c("python_legacy", "oof")) {
  run_one(mode, profile = TRUE)
  for (i in seq_len(repetitions)) {
    values <- run_one(mode)
    rows[[k]] <- data.frame(mode = mode, repetition = i,
                            elapsed_sec = values[["elapsed_sec"]],
                            gc_max_bytes = values[["gc_max_bytes"]])
    k <- k + 1L
  }
  prof <- summaryRprof(
    file.path(out_dir, paste0("late_fusion_", mode, ".Rprof")), memory = "both"
  )
  utils::write.csv(prof$by.total,
                   file.path(out_dir, paste0("late_fusion_", mode, "_calls.csv")))
}
results <- do.call(rbind, rows)
summaries <- do.call(rbind, lapply(split(results, results$mode), function(x) {
  data.frame(
    mode = x$mode[[1L]], repetitions = nrow(x),
    median_elapsed_sec = stats::median(x$elapsed_sec),
    elapsed_q1_sec = unname(stats::quantile(x$elapsed_sec, 0.25)),
    elapsed_q3_sec = unname(stats::quantile(x$elapsed_sec, 0.75)),
    median_gc_max_bytes = stats::median(x$gc_max_bytes),
    gc_max_q1_bytes = unname(stats::quantile(x$gc_max_bytes, 0.25)),
    gc_max_q3_bytes = unname(stats::quantile(x$gc_max_bytes, 0.75))
  )
}))
ratio <- summaries$median_elapsed_sec[summaries$mode == "oof"] /
  summaries$median_elapsed_sec[summaries$mode == "python_legacy"]
summaries$oof_to_legacy_elapsed_ratio <- ratio
utils::write.csv(results, file.path(out_dir, "late_fusion_benchmark_replicates.csv"),
                 row.names = FALSE)
utils::write.csv(summaries, file.path(out_dir, "late_fusion_benchmark_summary.csv"),
                 row.names = FALSE)
writeLines(c(
  paste("R", R.version.string),
  paste("stablr", as.character(utils::packageVersion("stablr"))),
  paste("repetitions", repetitions),
  "Identical workload for both modes: n=60, two p=10 omics, 10 bootstraps.",
  "OOF uses five stacking folds plus the full refit; legacy uses one full fit.",
  "gc_max_bytes is the base-R GC maximum after reset, not process peak RSS.",
  "Rprof and Rprofmem files are raw call-stack/allocation evidence."
), file.path(out_dir, "late_fusion_benchmark_manifest.txt"))
