#!/usr/bin/env Rscript

.release_cli_script_path <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command, value = TRUE)
  if (length(file_arg) != 1L) {
    stop("Could not identify the release-evidence CLI script.", call. = FALSE)
  }
  normalizePath(
    sub("^--file=", "", file_arg[[1L]]),
    winslash = "/",
    mustWork = TRUE
  )
}

.release_cli_path <- .release_cli_script_path()
.release_cli_module_path <- file.path(
  dirname(.release_cli_path), "release_evidence.R"
)
if (!file.exists(.release_cli_module_path)) {
  stop("The release-evidence module is missing beside its CLI.", call. = FALSE)
}
source(.release_cli_module_path, local = globalenv())

.release_cli_usage <- function() {
  paste(
    "Usage:",
    "  release_evidence_cli.R prepare --candidate FILE --runtime FILE --store DIR --out-ref FILE",
    "  release_evidence_cli.R run --context-ref FILE --validation-id ID --out-ref FILE",
    paste(
      "  release_evidence_cli.R assemble --context-ref FILE",
      "--packet-ref FILE --packet-ref FILE --out-ref FILE"
    ),
    sep = "\n"
  )
}

.release_cli_parse <- function(args) {
  if (!length(args) || args[[1L]] %in% c("--help", "-h")) {
    return(list(command = "help", options = list()))
  }
  command <- args[[1L]]
  if (!command %in% c("prepare", "run", "assemble")) {
    stop("Unknown command `", command, "`.", call. = FALSE)
  }
  args <- args[-1L]
  options <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args) ||
        startsWith(args[[i + 1L]], "--")) {
      stop("Every CLI option requires one value.", call. = FALSE)
    }
    key <- sub("^--", "", key)
    value <- args[[i + 1L]]
    if (identical(key, "packet-ref")) {
      options[[key]] <- c(options[[key]], value)
    } else {
      if (!is.null(options[[key]])) {
        stop("Duplicate CLI option `--", key, "`.", call. = FALSE)
      }
      options[[key]] <- value
    }
    i <- i + 2L
  }
  list(command = command, options = options)
}

.release_cli_require_options <- function(options, required) {
  if (!setequal(names(options), required) || anyDuplicated(names(options))) {
    stop(
      "Command requires exactly these options: ",
      paste(paste0("--", required), collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.release_cli_assert_candidate_code <- function(candidate) {
  candidate <- .release_normalize_absolute(candidate, "candidate", "file")
  inspection_dir <- tempfile("stablr-cli-candidate-")
  on.exit({
    if (dir.exists(inspection_dir)) .release_discard_attempt(inspection_dir)
  }, add = TRUE)
  inspection <- .release_inspect_candidate_archive(candidate, inspection_dir)
  candidate_module <- file.path(
    inspection$root, "inst", "analysis", "release_evidence.R"
  )
  candidate_cli <- file.path(
    inspection$root, "inst", "analysis", "release_evidence_cli.R"
  )
  if (!file.exists(candidate_module) || !file.exists(candidate_cli) ||
      !identical(
        .release_hash_file(candidate_module),
        .release_hash_file(.release_cli_module_path)
      ) || !identical(
        .release_hash_file(candidate_cli), .release_hash_file(.release_cli_path)
      )) {
    .release_identity_error(
      "Prepare CLI and module do not match the Candidate Source Artifact."
    )
  }
  invisible(candidate)
}

.release_cli_assert_context_code <- function(context) {
  verified <- .release_verify_context(context)
  candidate_analysis <- file.path(
    verified$record$library$path, "stablr", "analysis"
  )
  candidate_module <- file.path(candidate_analysis, "release_evidence.R")
  candidate_cli <- file.path(candidate_analysis, "release_evidence_cli.R")
  if (!identical(
    .release_hash_file(candidate_module),
    .release_hash_file(.release_cli_module_path)
  ) || !identical(
    .release_hash_file(candidate_cli), .release_hash_file(.release_cli_path)
  )) {
    .release_identity_error(
      "Run CLI and module do not match the installed candidate."
    )
  }
  verified
}

.release_cli_main <- function(parsed) {
  command <- parsed$command
  options <- parsed$options
  if (identical(command, "help")) {
    cat(.release_cli_usage(), "\n")
    return(invisible(NULL))
  }
  if (identical(command, "prepare")) {
    .release_cli_require_options(
      options, c("candidate", "runtime", "store", "out-ref")
    )
    candidate <- .release_cli_assert_candidate_code(options$candidate)
    context <- prepare_release_context(
      candidate_tarball = candidate,
      runtime = options$runtime,
      store = options$store
    )
    write_release_reference(context, options[["out-ref"]])
    return(invisible(context))
  }
  .release_cli_require_options(
    options,
    if (identical(command, "run")) {
      c("context-ref", "validation-id", "out-ref")
    } else {
      c("context-ref", "packet-ref", "out-ref")
    }
  )
  context <- read_release_reference(options[["context-ref"]])
  .release_cli_assert_context_code(context)
  if (identical(command, "run")) {
    packet <- tryCatch(
      run_required_validation(context, options[["validation-id"]]),
      stablr_validation_failed = function(e) {
        write_release_reference(e$packet_reference, options[["out-ref"]])
        stop(e)
      }
    )
    write_release_reference(packet, options[["out-ref"]])
    return(invisible(packet))
  }
  packet_refs <- lapply(
    options[["packet-ref"]], read_release_reference
  )
  release <- assemble_release_evidence(context, packet_refs)
  write_release_reference(release, options[["out-ref"]])
  invisible(release)
}

.release_cli_exit_status <- function(condition) {
  if (inherits(condition, "stablr_release_contract_error")) return(20L)
  if (inherits(condition, "stablr_release_identity_error")) return(21L)
  if (inherits(condition, "stablr_release_runtime_error")) return(22L)
  if (inherits(condition, "stablr_validation_incomplete")) return(23L)
  if (inherits(condition, "stablr_validation_failed")) return(24L)
  if (inherits(condition, "stablr_release_assembly_error")) return(25L)
  70L
}

tryCatch(
  .release_cli_main(.release_cli_parse(commandArgs(trailingOnly = TRUE))),
  error = function(e) {
    message("release evidence error: ", conditionMessage(e))
    quit(save = "no", status = .release_cli_exit_status(e), runLast = FALSE)
  }
)
