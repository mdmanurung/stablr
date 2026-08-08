#!/usr/bin/env Rscript

# Paired, independent-test validation of OOF versus historical late fusion.
`%||%` <- function(x, y) if (is.null(x)) y else x
.find_source_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script <- grep("^--file=", args, value = TRUE)
  starts <- c(getwd(), if (length(script)) dirname(sub("^--file=", "", script[[1L]])))
  for (start in starts) {
    path <- normalizePath(start, mustWork = TRUE)
    repeat {
      if (file.exists(file.path(path, "DESCRIPTION")) &&
          file.exists(file.path(path, "R", "stabl_fit.R"))) return(path)
      parent <- dirname(path)
      if (identical(parent, path)) break
      path <- parent
    }
  }
  NA_character_
}
.root <- .find_source_root()

.git_provenance <- function(root) {
  empty <- list(commit = NA_character_, tree = NA_character_, dirty = NA)
  if (is.na(root) || !file.exists(file.path(root, ".git")) ||
      !nzchar(Sys.which("git"))) return(empty)
  git <- function(...) {
    out <- suppressWarnings(system2(
      "git", c("-C", shQuote(root), ...), stdout = TRUE, stderr = TRUE
    ))
    if (!is.null(attr(out, "status")) && attr(out, "status") != 0L) {
      return(NA_character_)
    }
    paste(out, collapse = "\n")
  }
  status <- git("status", "--porcelain", "--untracked-files=no")
  list(
    commit = git("rev-parse", "HEAD"),
    tree = git("rev-parse", "HEAD^{tree}"),
    dirty = if (is.na(status)) NA else nzchar(status)
  )
}

.git_provenance_is_clean <- function(provenance) {
  is.character(provenance$commit) && length(provenance$commit) == 1L &&
    !is.na(provenance$commit) && nzchar(provenance$commit) &&
    is.character(provenance$tree) && length(provenance$tree) == 1L &&
    !is.na(provenance$tree) && nzchar(provenance$tree) &&
    identical(provenance$dirty, FALSE)
}

.git_provenance_is_available <- function(provenance) {
  is.character(provenance$commit) && length(provenance$commit) == 1L &&
    !is.na(provenance$commit) && nzchar(provenance$commit) &&
    is.character(provenance$tree) && length(provenance$tree) == 1L &&
    !is.na(provenance$tree) && nzchar(provenance$tree)
}

.git_provenance_is_stable <- function(start, end) {
  identical(start$commit, end$commit) &&
    identical(start$tree, end$tree) &&
    identical(start$dirty, end$dirty)
}

.late_validation_runtime <- new.env(parent = emptyenv())

.load_late_validation_package <- function(candidate_only = FALSE,
                                          candidate_path = NULL) {
  if (isTRUE(.late_validation_runtime$loaded)) {
    if (isTRUE(candidate_only)) {
      if (!is.character(candidate_path) || length(candidate_path) != 1L ||
          is.na(candidate_path) || !dir.exists(candidate_path)) {
        stop("Candidate validation requires an existing candidate path.",
             call. = FALSE)
      }
      expected_mode <- paste0(
        "installed_candidate:",
        normalizePath(candidate_path, winslash = "/", mustWork = TRUE)
      )
      if (!identical(.late_validation_runtime$mode, expected_mode)) {
        stop("Candidate validation cannot reuse a foreign package.",
             call. = FALSE)
      }
    }
    return(.late_validation_runtime$mode)
  }
  if (isTRUE(candidate_only)) {
    if (!is.character(candidate_path) || length(candidate_path) != 1L ||
        is.na(candidate_path) || !grepl("^/", candidate_path) ||
        !dir.exists(candidate_path)) {
      stop("Candidate validation requires an existing absolute candidate path.",
           call. = FALSE)
    }
    if (!requireNamespace("stablr", quietly = TRUE)) {
      stop("The isolated stablr candidate is unavailable.", call. = FALSE)
    }
    candidate_path <- normalizePath(
      candidate_path, winslash = "/", mustWork = TRUE
    )
    loaded_path <- normalizePath(
      getNamespaceInfo(asNamespace("stablr"), "path"),
      winslash = "/", mustWork = TRUE
    )
    if (!identical(loaded_path, candidate_path)) {
      stop("The loaded stablr namespace is not the declared candidate.",
           call. = FALSE)
    }
    mode <- paste0(
      "installed_candidate:", loaded_path
    )
  } else if ("stablr" %in% loadedNamespaces()) {
    mode <- paste0(
      "loaded:", getNamespaceInfo(asNamespace("stablr"), "path")
    )
  } else if (!is.na(.root) && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(.root, quiet = TRUE)
    mode <- paste0("source:", .root)
  } else {
    if (!requireNamespace("stablr", quietly = TRUE)) {
      stop("stablr must be installed.", call. = FALSE)
    }
    mode <- paste0("installed:", system.file(package = "stablr"))
  }
  .late_validation_runtime$loaded <- TRUE
  .late_validation_runtime$mode <- mode
  mode
}

.sha256_file <- function(path) {
  if (nzchar(Sys.which("sha256sum"))) {
    out <- system2("sha256sum", shQuote(path), stdout = TRUE)
    return(strsplit(out[[1L]], "[[:space:]]+")[[1L]][[1L]])
  }
  if (nzchar(Sys.which("shasum"))) {
    out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
    return(strsplit(out[[1L]], "[[:space:]]+")[[1L]][[1L]])
  }
  NA_character_
}

.args <- function(x) {
  out <- list(); i <- 1L
  while (i <= length(x)) {
    key <- sub("^--", "", x[[i]])
    if (i == length(x)) stop("Missing value for --", key)
    out[[key]] <- x[[i + 1L]]; i <- i + 2L
  }
  out
}

.simulate <- function(family, regime, seed, n_train = 90L, n_test = 120L,
                      p = 12L) {
  set.seed(seed)
  n <- n_train + n_test
  ids <- paste0("s", seq_len(n))
  z <- matrix(rnorm(n * p), n, p)
  x1 <- z + matrix(rnorm(n * p, sd = 0.15), n, p)
  x2 <- matrix(rnorm(n * p), n, p)
  colnames(x1) <- paste0("a", seq_len(p)); colnames(x2) <- paste0("b", seq_len(p))
  rownames(x1) <- rownames(x2) <- ids
  eta <- if (regime == "signal") 1.2 * x1[, 1] - x1[, 2] + 0.8 * x2[, 1] else rep(0, n)
  y <- switch(
    family,
    gaussian = eta + rnorm(n),
    binomial = rbinom(n, 1L, plogis(eta)),
    multinomial = {
      latent <- cbind(A = eta, B = -eta, C = 0.3 * x2[, 2]) +
        matrix(rnorm(n * 3L, sd = 0.7), n, 3L)
      factor(c("A", "B", "C")[max.col(latent)], levels = c("A", "B", "C"))
    }
  )
  names(y) <- ids
  tr <- seq_len(n_train); te <- seq.int(n_train + 1L, n)
  list(
    train_x = list(a = x1[tr, , drop = FALSE], b = x2[tr, , drop = FALSE]),
    test_x = list(a = x1[te, , drop = FALSE], b = x2[te, , drop = FALSE]),
    train_y = y[tr], test_y = y[te]
  )
}

.simulate_well_posed <- function(family, regime, seed, min_class_count = 10L,
                                  max_attempts = 100L) {
  for (attempt in 0:(max_attempts - 1L)) {
    data_seed <- seed + attempt * 1000003L
    out <- .simulate(family, regime, data_seed)
    if (identical(family, "gaussian") ||
        (min(table(out$train_y)) >= min_class_count &&
         min(table(out$test_y)) >= min_class_count)) {
      out$data_seed <- data_seed
      out$simulation_attempt <- attempt + 1L
      return(out)
    }
  }
  stop("Could not generate the predeclared minimum class counts.", call. = FALSE)
}

.metric <- function(family, truth, pred) {
  if (family == "gaussian") return(stablr:::.r_squared(as.numeric(truth), as.numeric(pred)))
  if (family == "binomial") return(stablr:::.r_auc(as.integer(truth), as.numeric(pred)))
  mean(as.character(truth) == as.character(pred$predicted_class))
}

.train_metric <- function(family, truth, fit) {
  if (family != "multinomial") return(fit$late_fusion$score)
  mean(as.character(truth) ==
         as.character(fit$late_fusion$train_predictions$predicted_class))
}

.fallback_rate <- function(fit) {
  reasons <- unlist(lapply(fit$late_fusion$provenance$folds, `[[`, "fallback_reasons"),
                    use.names = FALSE)
  if (!length(reasons)) return(NA_real_)
  mean(!is.na(reasons) & nzchar(reasons))
}

.capture_fit <- function(args) {
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      do.call(stabl_multiomic_train_validate, args),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  list(fit = fit, warnings = unique(warnings))
}

.late_fusion_common_gate_table <- function(summary, expected_replicates) {
  required <- c(
    "family", "regime", "replicates", "successful_replicates",
    "mean_legacy_optimism", "mean_oof_optimism", "mean_test_difference",
    "fallback_rate", "reduced_optimism", "noninferior", "fallback_ok"
  )
  if (!is.data.frame(summary) || !all(required %in% names(summary))) {
    stop("Late-fusion adapter summary does not match its gate schema.",
         call. = FALSE)
  }
  expected_cells <- 6L
  complete <- nrow(summary) == expected_cells &&
    all(summary$replicates == expected_replicates) &&
    all(summary$successful_replicates == expected_replicates)
  optimism_pass <- nrow(summary) == expected_cells &&
    all(!is.na(summary$reduced_optimism) & summary$reduced_optimism)
  noninferior_pass <- nrow(summary) == expected_cells &&
    all(!is.na(summary$noninferior) & summary$noninferior)
  signal <- summary[summary$regime == "signal", , drop = FALSE]
  fallback_pass <- nrow(signal) == 3L &&
    all(!is.na(signal$fallback_ok) & signal$fallback_ok)

  data.frame(
    gate_id = c(
      "cell_completeness", "reduced_optimism", "noninferiority",
      "fallback_rate"
    ),
    gate_version = c(
      "late-fusion-cell-completeness/v1",
      "late-fusion-reduced-optimism/v1",
      "late-fusion-noninferiority-0.02/v1",
      "late-fusion-signal-fallback-0.05/v1"
    ),
    scope = c(
      "all family-by-regime cells", "all family-by-regime cells",
      "all family-by-regime cells", "all signal cells"
    ),
    observed = c(
      paste0(
        sum(summary$successful_replicates), "/",
        expected_cells * expected_replicates, " replicates successful"
      ),
      paste0(
        sum(!is.na(summary$reduced_optimism) & summary$reduced_optimism),
        "/", expected_cells, " cells with reduced optimism"
      ),
      paste0(
        sum(!is.na(summary$noninferior) & summary$noninferior),
        "/", expected_cells, " cells noninferior"
      ),
      paste0(
        sum(!is.na(signal$fallback_ok) & signal$fallback_ok),
        "/3 signal cells below threshold"
      )
    ),
    criterion = c(
      "all six cells contain every predeclared successful replicate",
      "mean OOF optimism is below mean legacy optimism in every cell",
      "mean OOF-minus-legacy test score is >= -0.02 in every cell",
      "mean OOF fallback rate is < 0.05 in every signal cell"
    ),
    pass = c(complete, optimism_pass, noninferior_pass, fallback_pass),
    reason = c(
      if (complete) "all required replicates completed" else
        "one or more cells or replicates are missing or failed",
      if (optimism_pass) "every cell reduced optimism" else
        "one or more cells did not reduce optimism",
      if (noninferior_pass) "every cell met the noninferiority margin" else
        "one or more cells missed the noninferiority margin",
      if (fallback_pass) "every signal cell met the fallback threshold" else
        "one or more signal cells missed the fallback threshold"
    ),
    stringsAsFactors = FALSE
  )
}

run_late_fusion_validation <- function(out, replicates = 50L, n_bootstraps = 20L,
                                       n_iter = 500L, seed = 220711L,
                                       candidate_only = FALSE,
                                       fail_on_gates = TRUE,
                                       candidate_path = NULL) {
  if (!is.logical(candidate_only) || length(candidate_only) != 1L ||
      is.na(candidate_only) || !is.logical(fail_on_gates) ||
      length(fail_on_gates) != 1L || is.na(fail_on_gates)) {
    stop("`candidate_only` and `fail_on_gates` must be TRUE or FALSE.",
         call. = FALSE)
  }
  raw_settings <- list(replicates, n_bootstraps, n_iter, seed)
  if (any(lengths(raw_settings) != 1L)) {
    stop("Release settings must be scalar positive integers.", call. = FALSE)
  }
  values <- suppressWarnings(as.integer(unlist(
    raw_settings, use.names = FALSE
  )))
  numeric_values <- suppressWarnings(as.numeric(unlist(
    raw_settings, use.names = FALSE
  )))
  if (length(values) != 4L || anyNA(values) || anyNA(numeric_values) ||
      any(values < 1L) || !identical(as.numeric(values), numeric_values)) {
    stop("Release settings must be positive integers.", call. = FALSE)
  }
  replicates <- values[[1L]]
  n_bootstraps <- values[[2L]]
  n_iter <- values[[3L]]
  seed <- values[[4L]]
  if (isTRUE(candidate_only) && !identical(
    values, c(50L, 20L, 500L, 220711L)
  )) {
    stop("Candidate-only late-fusion settings do not match the contract.",
         call. = FALSE)
  }
  package_mode <- .load_late_validation_package(
    candidate_only,
    candidate_path = candidate_path
  )
  provenance_start <- if (isTRUE(candidate_only)) {
    list(commit = NA_character_, tree = NA_character_, dirty = NA)
  } else {
    .git_provenance(.root)
  }
  source_repository_present <- !is.na(.root) &&
    file.exists(file.path(.root, ".git"))
  if (!isTRUE(candidate_only) && (source_repository_present ||
       .git_provenance_is_available(provenance_start)) &&
      !.git_provenance_is_clean(provenance_start)) {
    stop(
      "Late-fusion release validation requires a clean Git source tree with an identifiable commit.",
      call. = FALSE
    )
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(out, winslash = "/", mustWork = TRUE)
  design <- expand.grid(
    family = c("gaussian", "binomial", "multinomial"),
    regime = c("null", "signal"), replicate = seq_len(replicates),
    stringsAsFactors = FALSE
  )
  rows <- vector("list", nrow(design))
  warning_rows <- list()
  for (i in seq_len(nrow(design))) {
    d <- design[i, ]; fit_seed <- seed + i * 101L
    dat <- .simulate_well_posed(d$family, d$regime, fit_seed)
    common <- list(
      x_train_list = dat$train_x, y_train = dat$train_y,
      lambda_grid = data.frame(lambda = c(0.2, 0.1, 0.05)),
      x_valid_list = dat$test_x, y_valid = dat$test_y,
      family = d$family, artificial_type = NULL, hard_threshold = 0.25,
      n_bootstraps = n_bootstraps, sample_fraction = 0.7,
      late_fusion = TRUE, late_fusion_nfolds = 5L,
      n_iter_lf = n_iter, random_state = fit_seed
    )
    if (d$family %in% c("binomial", "multinomial")) {
      common$stratify_bootstrap <- TRUE
      common$bootstrap_strata_train <- dat$train_y
    }
    legacy_result <- .capture_fit(c(common, list(late_fusion_training = "python_legacy")))
    oof_result <- .capture_fit(c(common, list(late_fusion_training = "oof")))
    for (mode in c("python_legacy", "oof")) {
      messages <- if (mode == "oof") oof_result$warnings else legacy_result$warnings
      if (length(messages)) warning_rows[[length(warning_rows) + 1L]] <- data.frame(
        family = d$family, regime = d$regime, replicate = d$replicate,
        seed = fit_seed, mode = mode, warning = messages,
        stringsAsFactors = FALSE
      )
    }
    if (inherits(legacy_result$fit, "error") || inherits(oof_result$fit, "error")) {
      rows[[i]] <- data.frame(
        family = d$family, regime = d$regime, replicate = d$replicate,
        seed = fit_seed, data_seed = dat$data_seed,
        simulation_attempt = dat$simulation_attempt, status = "error",
        legacy_train = NA_real_, oof_train = NA_real_,
        legacy_test = NA_real_, oof_test = NA_real_,
        legacy_optimism = NA_real_, oof_optimism = NA_real_,
        oof_fallback_rate = NA_real_,
        error = paste(
          if (inherits(legacy_result$fit, "error")) conditionMessage(legacy_result$fit) else "",
          if (inherits(oof_result$fit, "error")) conditionMessage(oof_result$fit) else "",
          sep = " | "
        ), stringsAsFactors = FALSE
      )
      next
    }
    legacy <- legacy_result$fit
    oof <- oof_result$fit
    legacy_test <- .metric(d$family, dat$test_y, legacy$late_fusion$valid_predictions)
    oof_test <- .metric(d$family, dat$test_y, oof$late_fusion$valid_predictions)
    legacy_train <- .train_metric(d$family, dat$train_y, legacy)
    oof_train <- .train_metric(d$family, dat$train_y, oof)
    rows[[i]] <- data.frame(
      family = d$family, regime = d$regime, replicate = d$replicate,
      seed = fit_seed, data_seed = dat$data_seed,
      simulation_attempt = dat$simulation_attempt,
      status = "ok", legacy_train = legacy_train,
      oof_train = oof_train, legacy_test = legacy_test,
      oof_test = oof_test,
      legacy_optimism = legacy_train - legacy_test,
      oof_optimism = oof_train - oof_test,
      oof_fallback_rate = .fallback_rate(oof), error = "", stringsAsFactors = FALSE
    )
  }
  results <- do.call(rbind, rows)
  groups <- split(results, interaction(results$family, results$regime, drop = TRUE))
  summary <- do.call(rbind, lapply(groups, function(x) data.frame(
    family = x$family[[1L]], regime = x$regime[[1L]], replicates = nrow(x),
    successful_replicates = sum(x$status == "ok"),
    mean_legacy_optimism = mean(x$legacy_optimism),
    mean_oof_optimism = mean(x$oof_optimism),
    mean_test_difference = mean(x$oof_test - x$legacy_test),
    fallback_rate = mean(x$oof_fallback_rate),
    reduced_optimism = all(x$status == "ok") &&
      mean(x$oof_optimism) < mean(x$legacy_optimism),
    noninferior = all(x$status == "ok") &&
      mean(x$oof_test - x$legacy_test) >= -0.02,
    fallback_ok = all(x$status == "ok") &&
      if (x$regime[[1L]] == "signal") mean(x$oof_fallback_rate) < 0.05 else TRUE,
    stringsAsFactors = FALSE
  )))
  results_path <- file.path(out, "late_fusion_replicates.csv")
  summary_path <- file.path(out, "late_fusion_summary.csv")
  warnings_path <- file.path(out, "late_fusion_warnings.csv")
  gates_path <- file.path(out, "late_fusion_gates.csv")
  utils::write.csv(results, results_path, row.names = FALSE)
  utils::write.csv(summary, summary_path, row.names = FALSE)
  warnings <- if (length(warning_rows)) do.call(rbind, warning_rows) else data.frame(
    family = character(), regime = character(), replicate = integer(),
    seed = integer(), mode = character(), warning = character()
  )
  utils::write.csv(warnings, warnings_path, row.names = FALSE)
  gates <- .late_fusion_common_gate_table(summary, replicates)
  utils::write.csv(gates, gates_path, row.names = FALSE)

  artifacts <- list(
    replicates = results_path,
    summary = summary_path,
    warnings = warnings_path,
    gates = gates_path
  )
  if (!isTRUE(candidate_only)) {
    provenance_end <- .git_provenance(.root)
    provenance_stable <- .git_provenance_is_stable(
      provenance_start,
      provenance_end
    )
    artifact_paths <- unlist(artifacts, use.names = FALSE)
    manifest <- file.path(out, "late_fusion_manifest.txt")
    writeLines(c(
      paste("R", R.version.string),
      paste("stablr", as.character(utils::packageVersion("stablr"))),
      paste("package_mode", package_mode),
      paste("source_git_commit", provenance_start$commit),
      paste("source_git_tree", provenance_start$tree),
      paste("source_tracked_dirty", provenance_start$dirty),
      paste("end_git_commit", provenance_end$commit),
      paste("end_git_tree", provenance_end$tree),
      paste("end_tracked_dirty", provenance_end$dirty),
      paste("source_git_stable", provenance_stable),
      paste("replicates_per_cell", replicates),
      paste("n_bootstraps", n_bootstraps), paste("n_iter", n_iter),
      "Families: gaussian, binomial, multinomial; regimes: null, signal.",
      paste(
        "Classification simulations require at least 10 samples per class",
        "in train and test; classification bootstraps are stratified."
      ),
      paste(
        "Independent test samples are never used to fit selectors,",
        "learners, or weights."
      ),
      paste(
        "Gates: fallback < 5% in signal, reduced directional optimism,",
        "test noninferiority margin 0.02."
      ),
      "artifact_sha256:",
      paste(
        basename(artifact_paths),
        vapply(artifact_paths, .sha256_file, character(1L))
      )
    ), manifest)
    artifacts$manifest <- manifest
    if (!isTRUE(provenance_stable)) {
      stop(
        paste(
          "Git source provenance changed during late-fusion validation;",
          "artifacts are not release evidence."
        ),
        call. = FALSE
      )
    }
  }
  if (isTRUE(fail_on_gates) && !all(gates$pass)) {
    stop("Late-fusion release gates failed; release must stop.", call. = FALSE)
  }
  invisible(artifacts)
}

if (sys.nframe() == 0L) {
  a <- .args(commandArgs(trailingOnly = TRUE))
  if (is.null(a$out)) stop("Usage: run_late_fusion_validation.R --out DIR [--replicates 50]")
  run_late_fusion_validation(
    a$out, as.integer(a$replicates %||% 50L),
    as.integer(a$`n-bootstraps` %||% 20L), as.integer(a$`n-iter` %||% 500L)
  )
}
