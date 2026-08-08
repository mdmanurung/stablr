#!/usr/bin/env Rscript

# Deterministic, bounded methodology validation experiments for stablr.
# This file is intentionally outside R/ so experiments can observe proposed
# methodology changes without changing package runtime behavior.

.parse_cli_args <- function(args) {
  if ("--help" %in% args || "-h" %in% args) {
    return(list(help = TRUE))
  }

  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key, call. = FALSE)
    }
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      stop("Missing value for ", key, call. = FALSE)
    }
    out[[sub("^-+", "", key)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

.usage <- function() {
  paste(
    "Usage:",
    "  run_methodology_validation.R --out /tmp/output-dir [options]",
    "",
    "Options:",
    "  --profile NAME             bounded (default) or locked release profile",
    "  --replicates N              Monte Carlo replicates per family/scenario/type. Default: 3",
    "  --n-bootstraps N            STABL bootstraps per fit. Default: 12",
    "  --n-lambda N                Auto lambda-grid size. Default: 6",
    "  --workers N                 Independent replicate workers. Default: 1",
    "  --families CSV              Default: gaussian (release: all advertised families)",
    "  --artificial-types CSV      Default: random_permutation,knockoff_equi",
    "  --scenarios CSV             Default: all bounded scenarios",
    "  --seed N                    Top-level seed. Default: 270627",
    "  --target-fdp X              Calibration reference threshold. Default: 0.1",
    "",
    "Artifacts:",
    "  methodology_validation_replicates.csv",
    "  methodology_validation_summary.csv",
    "  methodology_validation_warnings.csv",
    "  python_metrics_parity.csv",
    "  methodology_validation_gates.csv",
    "  methodology_validation_manifest.txt",
    sep = "\n"
  )
}

.advertised_families <- c("gaussian", "binomial", "multinomial", "cox")
.advertised_artificial_types <- c(
  "random_permutation", "knockoff", "knockoff_equi", "knockoff_mvr"
)

.release_profile_settings <- function() {
  list(
    replicates = 100L,
    n_bootstraps = 1000L,
    n_lambda = 30L,
    families = .advertised_families,
    artificial_types = .advertised_artificial_types,
    scenario_ids = "all"
  )
}

.csv_arg <- function(value, default) {
  if (is.null(value) || !nzchar(value)) return(default)
  trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
}

.int_arg <- function(value, default, name, min = 1L) {
  if (is.null(value) || !nzchar(value)) return(default)
  out <- suppressWarnings(as.integer(value))
  if (is.na(out) || out < min) {
    stop("`--", name, "` must be an integer >= ", min, ".", call. = FALSE)
  }
  out
}

.num_arg <- function(value, default, name, min = -Inf, max = Inf) {
  if (is.null(value) || !nzchar(value)) return(default)
  out <- suppressWarnings(as.numeric(value))
  if (is.na(out) || out < min || out > max) {
    stop("`--", name, "` must be numeric in [", min, ", ", max, "].", call. = FALSE)
  }
  out
}

.default_scenarios <- function() {
  data.frame(
    scenario = c(
      "null_independent",
      "null_correlated",
      "signal_independent",
      "signal_correlated",
      "null_high_dim",
      "signal_high_dim"
    ),
    regime = c("null", "null", "signal", "signal", "null", "signal"),
    profile = c("low_dim", "low_dim", "low_dim", "low_dim", "high_dim", "high_dim"),
    n = c(60L, 60L, 60L, 60L, 45L, 45L),
    p = c(30L, 30L, 30L, 30L, 75L, 75L),
    n_signal = c(0L, 0L, 5L, 5L, 0L, 5L),
    correlation = c(0.0, 0.7, 0.0, 0.7, 0.6, 0.6),
    noise_sd = c(1.0, 1.0, 0.8, 0.8, 1.0, 0.8),
    signal_strength = c(0, 0, 2, 2, 0, 2),
    stringsAsFactors = FALSE
  )
}

.select_scenarios <- function(ids) {
  scenarios <- .default_scenarios()
  if (identical(ids, "all")) return(scenarios)
  unknown <- setdiff(ids, scenarios$scenario)
  if (length(unknown)) {
    stop("Unknown scenario(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  scenarios[match(ids, scenarios$scenario), , drop = FALSE]
}

.ar1_matrix <- function(p, rho) {
  idx <- seq_len(p)
  outer(idx, idx, function(i, j) rho^abs(i - j))
}

.sample_multinomial <- function(probabilities) {
  apply(probabilities, 1L, function(p) sample.int(ncol(probabilities), 1L, prob = p))
}

.simulate_scenario <- function(scenario, family, seed) {
  set.seed(seed)
  sigma <- .ar1_matrix(scenario$p, scenario$correlation)
  z <- matrix(stats::rnorm(scenario$n * scenario$p), nrow = scenario$n)
  x <- z %*% chol(sigma)
  x <- scale(x)
  x <- matrix(
    as.numeric(x),
    nrow = scenario$n,
    dimnames = list(
      paste0("s", seq_len(scenario$n)),
      paste0("f", seq_len(scenario$p))
    )
  )

  if (scenario$n_signal > 0L) {
    beta <- seq(1.0, 0.45, length.out = scenario$n_signal)
    eta <- drop(x[, seq_len(scenario$n_signal), drop = FALSE] %*% beta)
    eta <- scenario$signal_strength * eta / stats::sd(eta)
    signal <- paste0("f", seq_len(scenario$n_signal))
  } else {
    eta <- rep.int(0, scenario$n)
    signal <- character(0L)
  }

  y <- switch(
    family,
    gaussian = eta + stats::rnorm(scenario$n, sd = scenario$noise_sd),
    binomial = stats::rbinom(scenario$n, 1L, stats::plogis(eta)),
    multinomial = {
      logits <- cbind(eta, -eta, rep.int(0, scenario$n))
      logits <- logits - apply(logits, 1L, max)
      probabilities <- exp(logits)
      probabilities <- probabilities / rowSums(probabilities)
      factor(.sample_multinomial(probabilities), levels = 1:3,
             labels = c("A", "B", "C"))
    },
    cox = {
      event_time <- stats::rexp(scenario$n, rate = 0.1 * exp(eta))
      censor_time <- stats::rexp(scenario$n, rate = 0.06)
      outcome <- survival::Surv(
        time = pmin(event_time, censor_time),
        event = event_time <= censor_time
      )
      rownames(outcome) <- rownames(x)
      outcome
    },
    stop("Unsupported family: ", family, call. = FALSE)
  )

  if (!identical(family, "cox")) names(y) <- rownames(x)

  list(
    x = x,
    y = y,
    signal = signal
  )
}

.fit_with_warning_capture <- function(data, family, artificial_type,
                                      n_bootstraps, n_lambda, seed) {
  warnings <- character()
  elapsed <- system.time({
    fit <- tryCatch(
      withCallingHandlers(
        stablr::stabl_fit(
          x = data$x,
          y = data$y,
          family = family,
          lambda_grid = "auto",
          n_lambda = n_lambda,
          n_bootstraps = n_bootstraps,
          artificial_type = artificial_type,
          random_state = seed
        ),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
  })

  list(
    fit = fit,
    warnings = warnings,
    elapsed_sec = unname(elapsed[["elapsed"]])
  )
}

.empty_replicate_row <- function() {
  data.frame(
    family = character(),
    scenario = character(),
    regime = character(),
    correlation = numeric(),
    profile = character(),
    replicate = integer(),
    artificial_type = character(),
    actual_artificial_type = character(),
    actual_artificial_types = character(),
    fallback_event_count = integer(),
    fallback_random_permutation_events = integer(),
    fallback_equi_events = integer(),
    fallback_history = character(),
    data_seed = integer(),
    fit_seed = integer(),
    status = character(),
    n = integer(),
    p = integer(),
    n_signal = integer(),
    n_bootstraps = integer(),
    n_lambda = integer(),
    fallback_random_permutation_warnings = integer(),
    fallback_equi_warnings = integer(),
    warning_count = integer(),
    elapsed_sec = numeric(),
    n_selected = integer(),
    true_positives = integer(),
    false_positives = integer(),
    empirical_fdp = numeric(),
    tpr = numeric(),
    min_fdp_plus = numeric(),
    fdp_threshold = numeric(),
    mean_max_real_score = numeric(),
    mean_max_artificial_score = numeric(),
    max_artificial_score = numeric(),
    selected_features = character(),
    error = character(),
    stringsAsFactors = FALSE
  )
}

.methodology_artificial_provenance <- function(fit, requested_type) {
  if (is.null(fit)) {
    return(list(
      actual_type = NA_character_,
      actual_types = NA_character_,
      fallback_event_count = NA_integer_,
      fallback_random_permutation_events = NA_integer_,
      fallback_equi_events = NA_integer_,
      fallback_history = NA_character_
    ))
  }

  provenance <- fit$artificial_provenance
  required <- c(
    "requested_type", "actual_type", "actual_types", "fallback_history"
  )
  if (!is.list(provenance) || !all(required %in% names(provenance))) {
    stop(
      "A successful methodology fit is missing structured artificial provenance.",
      call. = FALSE
    )
  }
  if (!identical(provenance$requested_type, requested_type)) {
    stop(
      "Artificial provenance requested type does not match the validation task.",
      call. = FALSE
    )
  }
  history <- provenance$fallback_history
  history_fields <- c(
    "event", "chunk", "step", "from_type", "to_type", "reason",
    "condition_class"
  )
  if (!is.data.frame(history) || !identical(names(history), history_fields)) {
    stop("Artificial fallback history has an invalid schema.", call. = FALSE)
  }
  if (nrow(history) && !identical(history$event, seq_len(nrow(history)))) {
    stop("Artificial fallback history is not in event order.", call. = FALSE)
  }

  serialized <- if (!nrow(history)) {
    ""
  } else {
    paste(
      sprintf(
        "%d:%d:%d:%s>%s:%s:%s",
        history$event,
        history$chunk,
        history$step,
        history$from_type,
        history$to_type,
        history$condition_class,
        history$reason
      ),
      collapse = " || "
    )
  }
  list(
    actual_type = provenance$actual_type,
    actual_types = paste(provenance$actual_types, collapse = ";"),
    fallback_event_count = nrow(history),
    fallback_random_permutation_events = sum(
      history$to_type == "random_permutation"
    ),
    fallback_equi_events = sum(history$to_type == "knockoff_equi"),
    fallback_history = serialized
  )
}

.replicate_row <- function(scenario, family, replicate, artificial_type,
                           data_seed, fit_seed, status,
                           n_bootstraps, n_lambda, warnings, elapsed_sec,
                           selected = character(), signal = character(),
                           fit = NULL, error = "") {
  fp <- length(setdiff(selected, signal))
  tp <- length(intersect(selected, signal))
  n_selected <- length(selected)
  art_scores <- if (!is.null(fit)) fit$stabl_scores_artificial_ else NULL
  real_scores <- if (!is.null(fit)) fit$stabl_scores_ else NULL
  provenance <- .methodology_artificial_provenance(fit, artificial_type)

  data.frame(
    family = family,
    scenario = scenario$scenario,
    regime = scenario$regime,
    correlation = scenario$correlation,
    profile = scenario$profile,
    replicate = replicate,
    artificial_type = artificial_type,
    actual_artificial_type = provenance$actual_type,
    actual_artificial_types = provenance$actual_types,
    fallback_event_count = provenance$fallback_event_count,
    fallback_random_permutation_events =
      provenance$fallback_random_permutation_events,
    fallback_equi_events = provenance$fallback_equi_events,
    fallback_history = provenance$fallback_history,
    data_seed = data_seed,
    fit_seed = fit_seed,
    status = status,
    n = scenario$n,
    p = scenario$p,
    n_signal = scenario$n_signal,
    n_bootstraps = n_bootstraps,
    n_lambda = n_lambda,
    fallback_random_permutation_warnings = sum(grepl("falling back to random permutation", warnings, fixed = TRUE)),
    fallback_equi_warnings = sum(grepl("using equi S", warnings, fixed = TRUE)),
    warning_count = length(warnings),
    elapsed_sec = elapsed_sec,
    n_selected = n_selected,
    true_positives = tp,
    false_positives = fp,
    empirical_fdp = fp / max(1L, n_selected),
    tpr = if (length(signal)) tp / length(signal) else NA_real_,
    min_fdp_plus = if (!is.null(fit)) fit$min_fdr_ else NA_real_,
    fdp_threshold = if (!is.null(fit)) fit$fdr_min_threshold_ else NA_real_,
    mean_max_real_score = if (!is.null(real_scores)) mean(apply(real_scores, 1L, max)) else NA_real_,
    mean_max_artificial_score = if (!is.null(art_scores)) mean(apply(art_scores, 1L, max)) else NA_real_,
    max_artificial_score = if (!is.null(art_scores)) max(apply(art_scores, 1L, max)) else NA_real_,
    selected_features = paste(selected, collapse = ";"),
    error = error,
    stringsAsFactors = FALSE
  )
}

.warning_rows <- function(scenario, family, replicate, artificial_type, warnings) {
  if (!length(warnings)) {
    return(data.frame(
      family = character(),
      scenario = character(),
      replicate = integer(),
      artificial_type = character(),
      warning_index = integer(),
      warning = character(),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    family = family,
    scenario = scenario$scenario,
    replicate = replicate,
    artificial_type = artificial_type,
    warning_index = seq_along(warnings),
    warning = warnings,
    stringsAsFactors = FALSE
  )
}

.run_one_condition <- function(scenario, family, replicate, artificial_type, data,
                               n_bootstraps, n_lambda, data_seed, fit_seed) {
  if (startsWith(artificial_type, "knockoff") &&
      !requireNamespace("knockoff", quietly = TRUE)) {
    return(list(
      replicate = .replicate_row(
        scenario = scenario,
        family = family,
        replicate = replicate,
        artificial_type = artificial_type,
        data_seed = data_seed,
        fit_seed = fit_seed,
        status = "skipped_missing_knockoff",
        n_bootstraps = n_bootstraps,
        n_lambda = n_lambda,
        warnings = character(),
        elapsed_sec = NA_real_,
        error = "The optional knockoff package is not installed."
      ),
      warnings = .warning_rows(
        scenario, family, replicate, artificial_type, character()
      )
    ))
  }

  result <- .fit_with_warning_capture(
    data = data,
    family = family,
    artificial_type = artificial_type,
    n_bootstraps = n_bootstraps,
    n_lambda = n_lambda,
    seed = fit_seed
  )

  if (inherits(result$fit, "error")) {
    row <- .replicate_row(
      scenario = scenario,
      family = family,
      replicate = replicate,
      artificial_type = artificial_type,
      data_seed = data_seed,
      fit_seed = fit_seed,
      status = "error",
      n_bootstraps = n_bootstraps,
      n_lambda = n_lambda,
      warnings = result$warnings,
      elapsed_sec = result$elapsed_sec,
      error = conditionMessage(result$fit)
    )
  } else {
    selected <- stablr::get_feature_names_out(result$fit)
    row <- .replicate_row(
      scenario = scenario,
      family = family,
      replicate = replicate,
      artificial_type = artificial_type,
      data_seed = data_seed,
      fit_seed = fit_seed,
      status = "ok",
      n_bootstraps = n_bootstraps,
      n_lambda = n_lambda,
      warnings = result$warnings,
      elapsed_sec = result$elapsed_sec,
      selected = selected,
      signal = data$signal,
      fit = result$fit
    )
  }

  list(
    replicate = row,
    warnings = .warning_rows(
      scenario = scenario,
      family = family,
      replicate = replicate,
      artificial_type = artificial_type,
      warnings = result$warnings
    )
  )
}

.mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

.sd_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x)
}

.se_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x) / sqrt(length(x))
}

.summarise_replicates <- function(rows, target_fdp) {
  if (!nrow(rows)) {
    return(data.frame(
      family = character(),
      scenario = character(),
      artificial_type = character(),
      actual_artificial_types = character(),
      generator_identity_complete = logical(),
      profile = character(),
      regime = character(),
      correlation = numeric(),
      n = integer(),
      p = integer(),
      n_signal = integer(),
      replicates = integer(),
      ok_replicates = integer(),
      mean_selected = numeric(),
      sd_selected = numeric(),
      se_selected = numeric(),
      mean_empirical_fdp = numeric(),
      sd_empirical_fdp = numeric(),
      se_empirical_fdp = numeric(),
      mean_tpr = numeric(),
      sd_tpr = numeric(),
      se_tpr = numeric(),
      empirical_fdp_exceedance_rate = numeric(),
      sd_empirical_fdp_exceedance = numeric(),
      se_empirical_fdp_exceedance_rate = numeric(),
      mean_min_fdp_plus = numeric(),
      mean_fdp_threshold = numeric(),
      fallback_random_permutation_rate = numeric(),
      se_fallback_random_permutation_rate = numeric(),
      fallback_equi_rate = numeric(),
      se_fallback_equi_rate = numeric(),
      mean_elapsed_sec = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  groups <- unique(rows[c("family", "scenario", "artificial_type")])
  summaries <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    idx <- rows$family == groups$family[[i]] &
      rows$scenario == groups$scenario[[i]] &
      rows$artificial_type == groups$artificial_type[[i]]
    d <- rows[idx, , drop = FALSE]
    ok <- d$status == "ok"
    actual_types <- sort(unique(d$actual_artificial_type[ok]))
    actual_types <- actual_types[!is.na(actual_types) & nzchar(actual_types)]
    fdp_exceeded <- as.numeric(d$empirical_fdp[ok] > target_fdp)
    fallback_rp <- as.numeric(
      d$fallback_random_permutation_events[ok] > 0L
    )
    fallback_equi <- as.numeric(d$fallback_equi_events[ok] > 0L)
    summaries[[i]] <- data.frame(
      family = d$family[[1L]],
      scenario = d$scenario[[1L]],
      artificial_type = d$artificial_type[[1L]],
      actual_artificial_types = paste(actual_types, collapse = ";"),
      generator_identity_complete =
        length(actual_types) == 1L &&
        identical(actual_types, d$artificial_type[[1L]]),
      profile = d$profile[[1L]],
      regime = d$regime[[1L]],
      correlation = d$correlation[[1L]],
      n = d$n[[1L]],
      p = d$p[[1L]],
      n_signal = d$n_signal[[1L]],
      replicates = nrow(d),
      ok_replicates = sum(ok),
      mean_selected = .mean_or_na(d$n_selected[ok]),
      sd_selected = .sd_or_na(d$n_selected[ok]),
      se_selected = .se_or_na(d$n_selected[ok]),
      mean_empirical_fdp = .mean_or_na(d$empirical_fdp[ok]),
      sd_empirical_fdp = .sd_or_na(d$empirical_fdp[ok]),
      se_empirical_fdp = .se_or_na(d$empirical_fdp[ok]),
      mean_tpr = .mean_or_na(d$tpr[ok]),
      sd_tpr = .sd_or_na(d$tpr[ok]),
      se_tpr = .se_or_na(d$tpr[ok]),
      empirical_fdp_exceedance_rate = .mean_or_na(fdp_exceeded),
      sd_empirical_fdp_exceedance = .sd_or_na(fdp_exceeded),
      se_empirical_fdp_exceedance_rate = .se_or_na(fdp_exceeded),
      mean_min_fdp_plus = .mean_or_na(d$min_fdp_plus[ok]),
      mean_fdp_threshold = .mean_or_na(d$fdp_threshold[ok]),
      fallback_random_permutation_rate = .mean_or_na(fallback_rp),
      se_fallback_random_permutation_rate = .se_or_na(fallback_rp),
      fallback_equi_rate = .mean_or_na(fallback_equi),
      se_fallback_equi_rate = .se_or_na(fallback_equi),
      mean_elapsed_sec = .mean_or_na(d$elapsed_sec[ok]),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, summaries)
}

.find_package_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    desc <- file.path(path, "DESCRIPTION")
    if (file.exists(desc)) {
      fields <- read.dcf(desc)
      if (identical(unname(fields[1L, "Package"]), "stablr")) return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) return(NA_character_)
    path <- parent
  }
}

.validation_runtime <- new.env(parent = emptyenv())

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

.git_provenance_is_stable <- function(start, end) {
  identical(start$commit, end$commit) &&
    identical(start$tree, end$tree) &&
    identical(start$dirty, end$dirty)
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

.load_validation_package <- function() {
  if (isTRUE(.validation_runtime$loaded)) return(.validation_runtime$mode)
  root <- .find_package_root()
  if ("stablr" %in% loadedNamespaces()) {
    mode <- paste0(
      "loaded:", getNamespaceInfo(asNamespace("stablr"), "path")
    )
  } else if (!is.na(root) && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(root, quiet = TRUE)
    mode <- paste0("source:", normalizePath(root, winslash = "/"))
  } else {
    if (!requireNamespace("stablr", quietly = TRUE)) {
      stop("stablr must be installed when source loading via pkgload is unavailable.",
           call. = FALSE)
    }
    mode <- paste0("installed:", system.file(package = "stablr"))
  }
  .validation_runtime$loaded <- TRUE
  .validation_runtime$mode <- mode
  mode
}

.metric_inputs <- function() {
  list(
    pair_pred = c("f1", "f2", "f5"),
    pair_true = c("f2", "f3", "f6"),
    list_of_lists = list(
      c("f1", "f2", "f3"),
      c("f2", "f3", "f4"),
      c("f1", "f4"),
      character(0)
    ),
    nb_total = 8L,
    d = 8L
  )
}

.similarity_scalar_observations <- function(inp) {
  adj_median <- stablr::adjusted_similarity_measure(inp$list_of_lists, inp$nb_total)
  adj_mean <- stablr::adjusted_similarity_measure(inp$list_of_lists, inp$nb_total, stat = "mean")
  pear_median <- stablr::pearson_similarity_measure(inp$list_of_lists, inp$d)
  pear_mean <- stablr::pearson_similarity_measure(inp$list_of_lists, inp$d, stat = "mean")

  c(
    jaccard_similarity_pair = stablr::jaccard_similarity(inp$pair_pred, inp$pair_true),
    adjusted_similarity_pair = stablr::adjusted_similarity(inp$pair_pred, inp$pair_true, inp$nb_total),
    pearson_similarity_pair = stablr::pearson_similarity(inp$pair_pred, inp$pair_true, inp$d),
    fdr_similarity_pair = stablr::fdr_similarity(inp$pair_pred, inp$pair_true),
    tpr_similarity_pair = stablr::tpr_similarity(inp$pair_pred, inp$pair_true),
    fscore_similarity_beta1_pair = stablr::fscore_similarity(inp$pair_pred, inp$pair_true, beta = 1),
    fscore_similarity_beta2_pair = stablr::fscore_similarity(inp$pair_pred, inp$pair_true, beta = 2),
    adjusted_similarity_measure_median_stat = adj_median$statistic,
    adjusted_similarity_measure_median_err_q25 = unname(adj_median$err[[1L]]),
    adjusted_similarity_measure_median_err_q75 = unname(adj_median$err[[2L]]),
    adjusted_similarity_measure_mean_stat = adj_mean$statistic,
    adjusted_similarity_measure_mean_err_sd = adj_mean$err,
    pearson_similarity_measure_median_stat = pear_median$statistic,
    pearson_similarity_measure_median_err_q25 = unname(pear_median$err[[1L]]),
    pearson_similarity_measure_median_err_q75 = unname(pear_median$err[[2L]]),
    pearson_similarity_measure_mean_stat = pear_mean$statistic,
    pearson_similarity_measure_mean_err_sd = pear_mean$err
  )
}

.similarity_vector_observations <- function(inp) {
  jaccard <- stablr::jaccard_matrix(inp$list_of_lists, remove_diag = TRUE)
  list(
    jaccard_matrix_remove_diag_rowmajor = as.vector(t(jaccard)),
    adjusted_similarity_values = stablr::adjusted_similarity_values(inp$list_of_lists, inp$nb_total),
    pearson_similarity_values = stablr::pearson_similarity_values(inp$list_of_lists, inp$d)
  )
}

.run_python_metrics_parity <- function(package_root, tolerance = 1e-12) {
  schema <- data.frame(
    metric = character(),
    index = integer(),
    reference = numeric(),
    observed = numeric(),
    abs_error = numeric(),
    status = character(),
    stringsAsFactors = FALSE
  )
  if (is.na(package_root)) {
    return(rbind(schema, data.frame(
      metric = "python_parity_fixture_dir",
      index = NA_integer_,
      reference = NA_real_,
      observed = NA_real_,
      abs_error = NA_real_,
      status = "skipped_package_root_not_found",
      stringsAsFactors = FALSE
    )))
  }

  fixture_dir <- file.path(package_root, "tests", "testthat", "fixtures", "python_parity")
  scalar_file <- file.path(fixture_dir, "metrics_scalars.csv")
  vector_file <- file.path(fixture_dir, "metrics_vectors.csv")
  if (!file.exists(scalar_file) || !file.exists(vector_file)) {
    return(rbind(schema, data.frame(
      metric = "python_parity_fixtures",
      index = NA_integer_,
      reference = NA_real_,
      observed = NA_real_,
      abs_error = NA_real_,
      status = "skipped_missing_fixture",
      stringsAsFactors = FALSE
    )))
  }

  scalars <- utils::read.csv(scalar_file, stringsAsFactors = FALSE, check.names = FALSE)
  vectors <- utils::read.csv(vector_file, stringsAsFactors = FALSE, check.names = FALSE)
  inp <- .metric_inputs()
  observed_scalars <- .similarity_scalar_observations(inp)
  observed_vectors <- .similarity_vector_observations(inp)

  scalar_rows <- lapply(seq_len(nrow(scalars)), function(i) {
    metric <- scalars$metric[[i]]
    observed <- unname(observed_scalars[[metric]])
    err <- abs(observed - scalars$value[[i]])
    data.frame(
      metric = metric,
      index = NA_integer_,
      reference = scalars$value[[i]],
      observed = observed,
      abs_error = err,
      status = if (is.finite(err) && err <= tolerance) "ok" else "mismatch",
      stringsAsFactors = FALSE
    )
  })

  vector_rows <- lapply(split(vectors, vectors$metric), function(ref) {
    metric <- ref$metric[[1L]]
    ref <- ref[order(ref$index), , drop = FALSE]
    observed <- observed_vectors[[metric]]
    n <- max(length(observed), nrow(ref))
    data.frame(
      metric = metric,
      index = seq_len(n),
      reference = c(ref$value, rep(NA_real_, n - nrow(ref))),
      observed = c(observed, rep(NA_real_, n - length(observed))),
      abs_error = abs(c(observed, rep(NA_real_, n - length(observed))) -
                        c(ref$value, rep(NA_real_, n - nrow(ref)))),
      status = ifelse(
        is.finite(abs(c(observed, rep(NA_real_, n - length(observed))) -
                        c(ref$value, rep(NA_real_, n - nrow(ref))))) &
          abs(c(observed, rep(NA_real_, n - length(observed))) -
                c(ref$value, rep(NA_real_, n - nrow(ref)))) <= tolerance,
        "ok",
        "mismatch"
      ),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, c(list(schema), scalar_rows, vector_rows))
}

.write_manifest <- function(path, settings, artifacts) {
  artifact_paths <- unlist(artifacts[c(
    "replicates", "summary", "warnings", "parity", "gates"
  )], use.names = FALSE)
  lines <- c(
    "stablr methodology validation manifest",
    paste0("created_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("seed: ", settings$seed),
    paste0("validation_profile: ", settings$profile),
    paste0("replicates: ", settings$replicates),
    paste0("n_bootstraps: ", settings$n_bootstraps),
    paste0("n_lambda: ", settings$n_lambda),
    paste0("target_fdp: ", settings$target_fdp),
    paste0("workers: ", settings$workers),
    paste0("package_mode: ", settings$package_mode),
    paste0("package_version: ", settings$package_version),
    paste0("source_git_commit: ", settings$git_start$commit),
    paste0("source_git_tree: ", settings$git_start$tree),
    paste0("source_tracked_dirty: ", settings$git_start$dirty),
    paste0("end_git_commit: ", settings$git_end$commit),
    paste0("end_git_tree: ", settings$git_end$tree),
    paste0("end_tracked_dirty: ", settings$git_end$dirty),
    paste0("source_git_stable: ", settings$git_stable),
    paste0("families: ", paste(settings$families, collapse = ",")),
    paste0("artificial_types: ", paste(settings$artificial_types, collapse = ",")),
    paste0("scenarios: ", paste(settings$scenarios$scenario, collapse = ",")),
    "",
    "hypotheses:",
    "1. FDP+ should avoid all-feature collapse under null settings, including correlated predictors.",
    "2. In signal settings, empirical FDP and TPR quantify the calibration-power tradeoff across artificial-feature strategies.",
    "3. Knockoff fallback frequency should be visible through warnings and additive fit provenance.",
    "4. Bundled Python metric fixtures should remain numerically identical to exported R metrics.",
    "5. The fitting rule selects the FDP+ minimizer; target_fdp is a reporting reference, not a fitting target or universal guarantee.",
    "",
    "artifacts:",
    paste0("replicates: ", artifacts$replicates),
    paste0("summary: ", artifacts$summary),
    paste0("warnings: ", artifacts$warnings),
    paste0("parity: ", artifacts$parity),
    paste0("gates: ", artifacts$gates),
    "",
    "artifact_sha256:",
    paste(basename(artifact_paths),
          vapply(artifact_paths, .sha256_file, character(1L)))
  )
  writeLines(lines, con = path)
}

.wilson_bound <- function(successes, trials, side = c("upper", "lower"),
                          confidence = 0.95) {
  side <- match.arg(side)
  if (trials <= 0L) return(NA_real_)
  z <- stats::qnorm(confidence)
  p <- successes / trials
  centre <- (p + z^2 / (2 * trials)) / (1 + z^2 / trials)
  half <- z * sqrt(p * (1 - p) / trials + z^2 / (4 * trials^2)) /
    (1 + z^2 / trials)
  if (side == "upper") min(1, centre + half) else max(0, centre - half)
}

.one_sided_mean_bound <- function(x, side = c("upper", "lower"),
                                  confidence = 0.95) {
  side <- match.arg(side)
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  half <- stats::qt(confidence, df = length(x) - 1L) *
    stats::sd(x) / sqrt(length(x))
  if (side == "upper") min(1, mean(x) + half) else max(0, mean(x) - half)
}

.selected_fraction_gate_version <-
  "SCI-04-null-selected-fraction-v2-draft"
.selected_fraction_gate_status <- "proposed_unaccepted"

.gate_cell_status <- function(cell, expected_replicates) {
  expected_ids <- seq_len(expected_replicates)
  if (!nrow(cell)) return("missing_cell")
  if (anyDuplicated(cell$replicate) ||
      !identical(sort(cell$replicate), expected_ids)) {
    return("incomplete_replicates")
  }
  if (any(cell$status != "ok")) return("fit_not_ok")
  "complete"
}

.release_gate_table <- function(rows, scenarios = NULL, families = NULL,
                                artificial_types = NULL,
                                expected_replicates = NULL) {
  if (!"actual_artificial_type" %in% names(rows)) {
    rows$actual_artificial_type <- rows$artificial_type
  }
  if (is.null(scenarios)) {
    scenarios <- unique(rows[c(
      "scenario", "regime", "profile", "correlation", "n", "p", "n_signal"
    )])
  }
  if (is.null(families)) families <- unique(rows$family)
  if (is.null(artificial_types)) {
    artificial_types <- unique(rows$artificial_type)
  }
  if (is.null(expected_replicates)) {
    expected_replicates <- if (nrow(rows)) max(rows$replicate) else 0L
  }

  gate_rows <- list()
  k <- 0L
  for (family in families) {
    for (scenario_i in seq_len(nrow(scenarios))) {
      scenario <- scenarios[scenario_i, , drop = FALSE]
      for (artificial_type in artificial_types) {
        cell <- rows[
          rows$family == family &
            rows$scenario == scenario$scenario &
            rows$artificial_type == artificial_type,
          , drop = FALSE
        ]
        cell_status <- .gate_cell_status(cell, expected_replicates)
        actual_generators <- sort(unique(cell$actual_artificial_type[
          cell$status == "ok"
        ]))
        actual_generators <- actual_generators[
          !is.na(actual_generators) & nzchar(actual_generators)
        ]
        generator_matches <- length(actual_generators) == 1L &&
          identical(actual_generators, artificial_type)
        if (identical(cell_status, "complete") && !generator_matches) {
          cell_status <- "generator_mismatch"
        }
        complete <- identical(cell_status, "complete")

        if (identical(scenario$regime, "null")) {
          select_any <- if (complete) {
            .wilson_bound(sum(cell$n_selected > 0L), expected_replicates, "upper")
          } else NA_real_
          # SCI-04: replicate-level selected fractions are continuous bounded
          # observations, not binomial successes. The proposed v2 decision
          # record is deliberately unaccepted, so this release gate stays
          # unavailable and fails closed until independent review is recorded.
          selected_fraction <- NA_real_
          collapses <- if (complete) {
            sum(cell$n_selected >= 0.9 * cell$p)
          } else NA_real_
          gates <- data.frame(
            gate = c("null_select_any", "null_selected_fraction",
                     "null_90pct_collapse"),
            bound = c(select_any, selected_fraction, collapses),
            criterion = c(
              "<= 0.10",
              "<= 0.10 after independent acceptance",
              "== 0"
            ),
            gate_version = c(
              "null-select-any-wilson-v1",
              .selected_fraction_gate_version,
              "null-90pct-collapse-count-v1"
            ),
            decision_status = c(
              "active",
              .selected_fraction_gate_status,
              "active"
            ),
            pass = complete & c(
              is.finite(select_any) && select_any <= 0.10,
              FALSE,
              is.finite(collapses) && collapses == 0L
            ),
            stringsAsFactors = FALSE
          )
        } else {
          fdp_upper <- if (complete) {
            .one_sided_mean_bound(cell$empirical_fdp, "upper")
          } else NA_real_
          tpr_lower <- if (complete) {
            .one_sided_mean_bound(cell$tpr, "lower")
          } else NA_real_
          gates <- data.frame(
            gate = c("signal_mean_fdp", "signal_tpr"),
            bound = c(fdp_upper, tpr_lower),
            criterion = c("<= 0.12", ">= 0.50"),
            gate_version = c(
              "signal-mean-fdp-t-v1",
              "signal-tpr-t-v1"
            ),
            decision_status = "active",
            pass = complete & c(
              is.finite(fdp_upper) && fdp_upper <= 0.12,
              is.finite(tpr_lower) && tpr_lower >= 0.50
            ),
            stringsAsFactors = FALSE
          )
        }

        k <- k + 1L
        gate_rows[[k]] <- cbind(
          data.frame(
            family = family,
            scenario = scenario$scenario,
            regime = scenario$regime,
            profile = scenario$profile,
            artificial_type = artificial_type,
            actual_artificial_type = paste(
              actual_generators,
              collapse = ";"
            ),
            expected_replicates = expected_replicates,
            observed_replicates = nrow(cell),
            ok_replicates = sum(cell$status == "ok"),
            cell_status = cell_status,
            stringsAsFactors = FALSE
          ),
          gates
        )
      }
    }
  }
  do.call(rbind, gate_rows)
}

run_methodology_validation <- function(out,
                                       profile = c("bounded", "release"),
                                       replicates = 3L,
                                       n_bootstraps = 12L,
                                       n_lambda = 6L,
                                       families = "gaussian",
                                       artificial_types = c("random_permutation", "knockoff_equi"),
                                       scenario_ids = "all",
                                       seed = 270627L,
                                       target_fdp = 0.1,
                                       workers = 1L) {
  profile <- match.arg(profile)
  package_root <- .find_package_root()
  git_start <- .git_provenance(package_root)
  if (identical(profile, "release") &&
      !.git_provenance_is_clean(git_start)) {
    stop(
      "Locked release methodology validation requires a clean Git source tree with an identifiable commit.",
      call. = FALSE
    )
  }
  package_mode <- .load_validation_package()
  if (identical(profile, "release")) {
    locked <- .release_profile_settings()
    replicates <- locked$replicates
    n_bootstraps <- locked$n_bootstraps
    n_lambda <- locked$n_lambda
    families <- locked$families
    artificial_types <- locked$artificial_types
    scenario_ids <- locked$scenario_ids
  }
  if (missing(out) || is.null(out) || !nzchar(out)) {
    stop("`out` must be a non-empty output directory.", call. = FALSE)
  }
  replicates <- as.integer(replicates)
  n_bootstraps <- as.integer(n_bootstraps)
  n_lambda <- as.integer(n_lambda)
  seed <- as.integer(seed)
  workers <- as.integer(workers)
  if (any(is.na(c(replicates, n_bootstraps, n_lambda, seed, workers))) ||
      replicates < 1L || n_bootstraps < 1L || n_lambda < 1L || workers < 1L) {
    stop("`replicates`, `n_bootstraps`, `n_lambda`, `seed`, and `workers` must be positive integers.",
         call. = FALSE)
  }
  if (workers > 1L && .Platform$OS.type == "windows") {
    stop("`workers > 1` is currently supported only on Unix-alike systems.",
         call. = FALSE)
  }
  if (!is.numeric(target_fdp) || length(target_fdp) != 1L ||
      is.na(target_fdp) || target_fdp < 0 || target_fdp > 1) {
    stop("`target_fdp` must be a numeric scalar in [0, 1].", call. = FALSE)
  }
  if (!is.character(artificial_types) || length(artificial_types) == 0L ||
      anyNA(artificial_types) || any(!nzchar(artificial_types))) {
    stop("`artificial_types` must be a non-empty character vector.", call. = FALSE)
  }
  if (anyDuplicated(artificial_types)) {
    stop(
      "`artificial_types` must not contain duplicates; duplicate same-seed runs are not independent replicates.",
      call. = FALSE
    )
  }
  if (!all(artificial_types %in% .advertised_artificial_types)) {
    stop(
      "Unknown `artificial_types`: ",
      paste(setdiff(artificial_types, .advertised_artificial_types), collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.character(families) || !length(families) || anyNA(families) ||
      any(!nzchar(families)) || anyDuplicated(families)) {
    stop("`families` must be a unique, non-empty character vector.", call. = FALSE)
  }
  if (!all(families %in% .advertised_families)) {
    stop(
      "Unknown `families`: ",
      paste(setdiff(families, .advertised_families), collapse = ", "),
      call. = FALSE
    )
  }

  scenarios <- .select_scenarios(scenario_ids)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(out, mustWork = TRUE)

  tasks <- expand.grid(
    family_i = seq_along(families),
    scenario_i = seq_len(nrow(scenarios)),
    replicate = seq_len(replicates),
    KEEP.OUT.ATTRS = FALSE
  )
  run_task <- function(task_i) {
    task <- tasks[task_i, , drop = FALSE]
    family_i <- task$family_i[[1L]]
    scenario_i <- task$scenario_i[[1L]]
    replicate <- task$replicate[[1L]]
    family <- families[[family_i]]
    scenario <- scenarios[scenario_i, , drop = FALSE]
    data_seed <- seed + family_i * 1000000L +
      scenario_i * 100000L + replicate * 1000L
    fit_seed <- data_seed + 1L
    data <- .simulate_scenario(scenario, family = family, seed = data_seed)
    conditions <- lapply(artificial_types, function(artificial_type) {
      .run_one_condition(
        scenario = scenario, family = family, replicate = replicate,
        artificial_type = artificial_type, data = data,
        n_bootstraps = n_bootstraps, n_lambda = n_lambda,
        data_seed = data_seed, fit_seed = fit_seed
      )
    })
    list(
      replicates = lapply(conditions, `[[`, "replicate"),
      warnings = lapply(conditions, `[[`, "warnings")
    )
  }
  task_results <- if (workers == 1L) {
    lapply(seq_len(nrow(tasks)), run_task)
  } else {
    parallel::mclapply(
      seq_len(nrow(tasks)), run_task,
      mc.cores = min(workers, nrow(tasks)), mc.preschedule = TRUE
    )
  }
  replicate_rows <- unlist(lapply(task_results, `[[`, "replicates"),
                           recursive = FALSE)
  warning_rows <- unlist(lapply(task_results, `[[`, "warnings"),
                         recursive = FALSE)

  replicates_df <- if (length(replicate_rows)) {
    do.call(rbind, replicate_rows)
  } else {
    .empty_replicate_row()
  }
  warnings_df <- if (length(warning_rows)) {
    do.call(rbind, warning_rows)
  } else {
    .warning_rows(
      scenarios[1L, , drop = FALSE], "", 1L, "", character()
    )
  }
  summary_df <- .summarise_replicates(replicates_df, target_fdp = target_fdp)
  parity_df <- .run_python_metrics_parity(.find_package_root())

  artifacts <- list(
    replicates = file.path(out, "methodology_validation_replicates.csv"),
    summary = file.path(out, "methodology_validation_summary.csv"),
    warnings = file.path(out, "methodology_validation_warnings.csv"),
    parity = file.path(out, "python_metrics_parity.csv"),
    manifest = file.path(out, "methodology_validation_manifest.txt"),
    gates = file.path(out, "methodology_validation_gates.csv")
  )

  utils::write.csv(replicates_df, artifacts$replicates, row.names = FALSE)
  utils::write.csv(summary_df, artifacts$summary, row.names = FALSE)
  utils::write.csv(warnings_df, artifacts$warnings, row.names = FALSE)
  utils::write.csv(parity_df, artifacts$parity, row.names = FALSE)
  gates_df <- .release_gate_table(
    replicates_df,
    scenarios = scenarios,
    families = families,
    artificial_types = artificial_types,
    expected_replicates = replicates
  )
  utils::write.csv(gates_df, artifacts$gates, row.names = FALSE)
  git_end <- .git_provenance(package_root)
  git_stable <- .git_provenance_is_stable(git_start, git_end)
  .write_manifest(
    path = artifacts$manifest,
    settings = list(
      seed = seed,
      profile = profile,
      replicates = replicates,
      n_bootstraps = n_bootstraps,
      n_lambda = n_lambda,
      target_fdp = target_fdp,
      workers = workers,
      package_mode = package_mode,
      package_version = as.character(utils::packageVersion("stablr")),
      families = families,
      artificial_types = artificial_types,
      scenarios = scenarios,
      git_start = git_start,
      git_end = git_end,
      git_stable = git_stable
    ),
    artifacts = artifacts
  )

  if (identical(profile, "release") && !isTRUE(git_stable)) {
    stop(
      "Git source provenance changed during locked release methodology validation; artifacts are not release evidence.",
      call. = FALSE
    )
  }

  artifacts
}

.main <- function(args = commandArgs(trailingOnly = TRUE)) {
  parsed <- .parse_cli_args(args)
  if (isTRUE(parsed$help)) {
    cat(.usage(), "\n")
    return(invisible(0L))
  }
  if (is.null(parsed$out) || !nzchar(parsed$out)) {
    stop(.usage(), call. = FALSE)
  }

  profile <- if (is.null(parsed$profile)) "bounded" else parsed$profile
  defaults <- if (identical(profile, "release")) {
    list(replicates = 100L, n_bootstraps = 1000L, n_lambda = 30L)
  } else {
    list(replicates = 3L, n_bootstraps = 12L, n_lambda = 6L)
  }
  artifacts <- run_methodology_validation(
    out = parsed$out,
    profile = profile,
    replicates = .int_arg(parsed$replicates, defaults$replicates, "replicates"),
    n_bootstraps = .int_arg(parsed[["n-bootstraps"]], defaults$n_bootstraps, "n-bootstraps"),
    n_lambda = .int_arg(parsed[["n-lambda"]], defaults$n_lambda, "n-lambda"),
    families = .csv_arg(parsed$families, "gaussian"),
    artificial_types = .csv_arg(parsed[["artificial-types"]],
                                c("random_permutation", "knockoff_equi")),
    scenario_ids = .csv_arg(parsed$scenarios, "all"),
    seed = .int_arg(parsed$seed, 270627L, "seed"),
    target_fdp = .num_arg(parsed[["target-fdp"]], 0.1, "target-fdp", min = 0, max = 1),
    workers = .int_arg(parsed$workers, 1L, "workers")
  )

  message("Wrote methodology validation artifacts:")
  for (name in names(artifacts)) {
    message("  ", name, ": ", artifacts[[name]])
  }
  if (identical(profile, "release")) {
    gates <- utils::read.csv(artifacts$gates, stringsAsFactors = FALSE)
    if (!all(gates$pass)) {
      stop("Locked release methodology gates failed; release must stop.",
           call. = FALSE)
    }
  }
  invisible(artifacts)
}

if (sys.nframe() == 0L && !interactive()) {
  .main()
}
