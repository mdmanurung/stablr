# Artifact-first release-evidence lifecycle for stablr.
#
# This module is intentionally internal. The only caller-facing operations are
# prepare_release_context(), run_required_validation(), and
# assemble_release_evidence(). The command-line wrapper sources this file from
# a freshly extracted candidate artifact.

`%||%` <- function(x, y) if (is.null(x)) y else x

.release_new_condition <- function(subclass, message, call = NULL, ...) {
  structure(
    c(list(message = message, call = call), list(...)),
    class = c(subclass, "stablr_release_error", "error", "condition")
  )
}

.release_abort <- function(subclass, message, call = NULL, ...) {
  stop(.release_new_condition(subclass, message, call = call, ...))
}

.release_contract_error <- function(message, ...) {
  .release_abort("stablr_release_contract_error", message, ...)
}

.release_identity_error <- function(message, ...) {
  .release_abort("stablr_release_identity_error", message, ...)
}

.release_runtime_error <- function(message, ...) {
  .release_abort("stablr_release_runtime_error", message, ...)
}

.release_assembly_error <- function(message, ...) {
  .release_abort("stablr_release_assembly_error", message, ...)
}

.release_require_jsonlite <- function() {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    .release_runtime_error(
      "Release tooling requires the contract-declared `jsonlite` package."
    )
  }
  invisible(TRUE)
}

.release_is_scalar_character <- function(x, nonempty = TRUE) {
  is.character(x) && length(x) == 1L && !is.na(x) &&
    (!nonempty || nzchar(x))
}

.release_is_scalar_logical <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}

.release_is_scalar_number <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x)
}

.release_is_sha256 <- function(x) {
  .release_is_scalar_character(x) && grepl("^[0-9a-f]{64}$", x)
}

.release_assert_exact_fields <- function(x, expected, label,
                                         condition = "contract") {
  actual <- names(x)
  valid <- is.list(x) && !is.null(actual) &&
    setequal(actual, expected) && !anyDuplicated(actual)
  if (isTRUE(valid)) return(invisible(TRUE))
  message <- paste0(
    label, " must contain exactly these fields: ",
    paste(expected, collapse = ", "), "."
  )
  if (identical(condition, "contract")) {
    .release_contract_error(message, record = label)
  }
  .release_identity_error(message, record = label)
}

.release_json_array <- function(x, label) {
  if (!is.list(x) || is.object(x)) {
    .release_contract_error(label, " must be a JSON array.")
  }
  x
}

.release_character_array <- function(x, label, unique = TRUE) {
  values <- unlist(.release_json_array(x, label), use.names = FALSE)
  if (!is.character(values) || !length(values) || anyNA(values) ||
      any(!nzchar(values)) || (unique && anyDuplicated(values))) {
    .release_contract_error(
      label, " must be a non-empty array of unique non-empty strings."
    )
  }
  values
}

.release_sort_json_object <- function(x) {
  if (!is.list(x) || is.data.frame(x)) return(x)
  if (is.null(names(x))) return(lapply(x, .release_sort_json_object))
  ordered <- x[order(names(x), method = "radix")]
  lapply(ordered, .release_sort_json_object)
}

.release_canonical_json <- function(x) {
  .release_require_jsonlite()
  jsonlite::toJSON(
    .release_sort_json_object(x),
    auto_unbox = TRUE,
    null = "null",
    na = "null",
    digits = NA,
    pretty = FALSE,
    force = TRUE
  )
}

.release_read_json <- function(path, label = "JSON record") {
  .release_require_jsonlite()
  if (!.release_is_scalar_character(path) || !file.exists(path) ||
      dir.exists(path)) {
    .release_identity_error(label, " is not an existing regular file.")
  }
  tryCatch(
    jsonlite::fromJSON(path, simplifyVector = FALSE),
    error = function(e) {
      .release_identity_error(
        label, " is not valid JSON: ", conditionMessage(e), parent = e
      )
    }
  )
}

.release_hash_file <- function(path) {
  if (!.release_is_scalar_character(path) || !file.exists(path) ||
      dir.exists(path) || grepl("[\r\n]", path)) {
    .release_identity_error("Cannot hash a missing or invalid regular file.")
  }
  tool <- Sys.which("sha256sum")
  if (!nzchar(tool)) {
    .release_runtime_error("The required `sha256sum` executable is unavailable.")
  }
  output <- suppressWarnings(system2(
    unname(tool), shQuote(path), stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status") %||% 0L
  hash <- if (length(output)) {
    strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
  } else {
    NA_character_
  }
  if (status != 0L || !.release_is_sha256(hash)) {
    .release_runtime_error(
      "SHA-256 hashing failed for `", path, "`.",
      hash_output = output,
      hash_status = status
    )
  }
  hash
}

.release_write_canonical_json <- function(record, path) {
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    .release_identity_error("JSON destination parent does not exist: ", parent)
  }
  if (file.exists(path)) {
    .release_identity_error("Refusing to replace existing JSON record: ", path)
  }
  temporary <- file.path(
    parent,
    paste0(".", basename(path), ".", .release_random_id(), ".tmp")
  )
  connection <- file(temporary, open = "wb")
  tryCatch(
    writeChar(
      as.character(.release_canonical_json(record)),
      connection,
      eos = NULL,
      useBytes = TRUE
    ),
    finally = close(connection)
  )
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    .release_identity_error("Atomic JSON publication failed: ", path)
  }
  .release_hash_file(path)
}

.release_random_id <- function(bytes = 16L) {
  source <- "/dev/urandom"
  if (.Platform$OS.type != "unix" || !file.exists(source)) {
    .release_runtime_error(
      "Release evidence v1 requires a POSIX runtime with `/dev/urandom`."
    )
  }
  connection <- file(source, open = "rb", raw = TRUE)
  on.exit(close(connection), add = TRUE)
  value <- readBin(connection, what = "raw", n = bytes)
  if (length(value) != bytes) {
    .release_runtime_error("Could not allocate a random immutable record ID.")
  }
  paste(sprintf("%02x", as.integer(value)), collapse = "")
}

.release_expected_contract <- function() {
  list(
    schema_version = "stablr.release-contract/v1",
    contract_id = "stablr-0.1.1-release-v1",
    candidate = c(package = "stablr", version = "0.1.1"),
    hash_algorithm = "sha256",
    promotable = FALSE,
    blockers = "SCI-04-null-selected-fraction-v2-pending-independent-acceptance",
    dependency_roots = c(
      "jsonlite", "glmnet", "Matrix", "matrixStats", "RColorBrewer",
      "Rcpp", "RcppEigen", "survival", "knockoff"
    ),
    environment = c(
      TZ = "UTC",
      LC_ALL = "C",
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      R_PROFILE_USER = "/dev/null",
      R_ENVIRON_USER = "/dev/null"
    ),
    validation_ids = c("methodology", "late_fusion"),
    adapters = c("methodology_validation", "late_fusion_validation"),
    settings = list(
      methodology = list(
        profile = "release",
        replicates = 100L,
        n_bootstraps = 1000L,
        n_lambda = 30L,
        families = c("gaussian", "binomial", "multinomial", "cox"),
        artificial_types = c(
          "random_permutation", "knockoff", "knockoff_equi", "knockoff_mvr"
        ),
        scenario_ids = c(
          "null_independent", "null_correlated", "signal_independent",
          "signal_correlated", "null_high_dim", "signal_high_dim"
        ),
        seed = 270627L,
        target_fdp = 0.1,
        workers = 32L
      ),
      late_fusion = list(
        replicates = 50L,
        n_bootstraps = 20L,
        n_iter = 500L,
        seed = 220711L
      )
    ),
    artifact_names = list(
      methodology = c("replicates", "summary", "warnings", "parity", "gates"),
      late_fusion = c("replicates", "summary", "warnings", "gates")
    ),
    gate_ids = list(
      methodology = c(
        "cell_completeness", "python_metrics_parity", "null_select_any",
        "null_selected_fraction", "null_select_90pct", "signal_mean_fdp",
        "signal_tpr"
      ),
      late_fusion = c(
        "cell_completeness", "reduced_optimism", "noninferiority",
        "fallback_rate"
      )
    )
  )
}

.release_validate_settings <- function(settings, expected, validation_id) {
  .release_assert_exact_fields(
    settings, names(expected),
    paste0("settings for required validation `", validation_id, "`")
  )
  for (name in names(expected)) {
    observed <- settings[[name]]
    wanted <- expected[[name]]
    if (is.list(observed) && is.character(wanted)) {
      observed <- unlist(observed, use.names = FALSE)
    }
    if (is.integer(wanted)) observed <- suppressWarnings(as.integer(observed))
    if (is.double(wanted)) observed <- suppressWarnings(as.numeric(observed))
    if (!identical(observed, wanted)) {
      .release_contract_error(
        "Contract setting `", validation_id, ".", name,
        "` does not match the locked v1 value."
      )
    }
  }
  invisible(TRUE)
}

.release_validate_artifacts <- function(artifacts, validation_id,
                                        expected_names) {
  artifacts <- .release_json_array(
    artifacts, paste0("artifacts for `", validation_id, "`")
  )
  if (length(artifacts) != length(expected_names)) {
    .release_contract_error(
      "Required validation `", validation_id,
      "` has a missing or extra artifact declaration."
    )
  }
  observed_names <- character(length(artifacts))
  for (i in seq_along(artifacts)) {
    artifact <- artifacts[[i]]
    .release_assert_exact_fields(
      artifact, c("name", "path", "schema"),
      paste0("artifact declaration ", i, " for `", validation_id, "`")
    )
    if (!all(vapply(
      artifact, .release_is_scalar_character, logical(1L)
    ))) {
      .release_contract_error("Artifact declarations must contain strings.")
    }
    path <- artifact$path
    components <- strsplit(path, "/", fixed = TRUE)[[1L]]
    if (grepl("^/", path) || any(components %in% c("", ".", "..")) ||
        !identical(basename(path), path)) {
      .release_contract_error(
        "Contract artifact paths must be plain relative file names."
      )
    }
    observed_names[[i]] <- artifact$name
  }
  if (!identical(observed_names, expected_names) ||
      anyDuplicated(observed_names)) {
    .release_contract_error(
      "Required validation `", validation_id,
      "` artifact IDs do not match the locked v1 contract."
    )
  }
  invisible(TRUE)
}

.release_validate_gates <- function(gates, validation_id, expected_ids) {
  gates <- .release_json_array(
    gates, paste0("scientific_gates for `", validation_id, "`")
  )
  if (length(gates) != length(expected_ids)) {
    .release_contract_error(
      "Required validation `", validation_id,
      "` has a missing or extra Scientific Gate declaration."
    )
  }
  observed <- character(length(gates))
  for (i in seq_along(gates)) {
    gate <- gates[[i]]
    .release_assert_exact_fields(
      gate, c("gate_id", "gate_version", "decision_status"),
      paste0("Scientific Gate declaration ", i, " for `", validation_id, "`")
    )
    if (!all(vapply(gate, .release_is_scalar_character, logical(1L)))) {
      .release_contract_error("Scientific Gate declarations must contain strings.")
    }
    if (!gate$decision_status %in% c("accepted", "proposed_unaccepted")) {
      .release_contract_error("Unknown Scientific Gate decision status.")
    }
    observed[[i]] <- gate$gate_id
  }
  if (!identical(observed, expected_ids) || anyDuplicated(observed)) {
    .release_contract_error(
      "Required validation `", validation_id,
      "` Scientific Gate IDs do not match the locked v1 contract."
    )
  }
  selected_fraction <- gates[[match("null_selected_fraction", observed)]]
  if (identical(validation_id, "methodology") &&
      (!identical(
        selected_fraction$gate_version,
        "SCI-04-null-selected-fraction-v2-draft"
      ) || !identical(
        selected_fraction$decision_status,
        "proposed_unaccepted"
      ))) {
    .release_contract_error(
      "SCI-04 must remain proposed and unaccepted until independent review."
    )
  }
  invisible(TRUE)
}

.release_validate_contract <- function(contract) {
  expected <- .release_expected_contract()
  .release_assert_exact_fields(
    contract,
    c(
      "schema_version", "contract_id", "candidate", "hash_algorithm",
      "promotable", "blockers", "runtime", "required_validations"
    ),
    "Release Contract"
  )
  if (!identical(contract$schema_version, expected$schema_version)) {
    .release_contract_error("Unknown Release Contract schema version.")
  }
  if (!identical(contract$contract_id, expected$contract_id) ||
      !identical(contract$hash_algorithm, expected$hash_algorithm)) {
    .release_contract_error("Unsupported Release Contract identity or hash algorithm.")
  }
  .release_assert_exact_fields(
    contract$candidate, c("package", "version"), "candidate constraint"
  )
  candidate <- unlist(contract$candidate, use.names = TRUE)
  if (!identical(candidate, expected$candidate)) {
    .release_contract_error("Candidate package/version constraints were weakened.")
  }
  if (!identical(contract$promotable, expected$promotable)) {
    .release_contract_error(
      "The v1 contract must remain non-promotable while SCI-04 is unaccepted."
    )
  }
  blockers <- .release_character_array(contract$blockers, "blockers")
  if (!identical(blockers, expected$blockers)) {
    .release_contract_error("Release Contract blockers do not match v1.")
  }
  .release_assert_exact_fields(
    contract$runtime,
    c("platform", "dependency_roots", "environment"),
    "runtime requirements"
  )
  if (!identical(contract$runtime$platform, "linux-posix")) {
    .release_contract_error("Unsupported Release Contract runtime platform.")
  }
  dependencies <- .release_character_array(
    contract$runtime$dependency_roots, "runtime dependency_roots"
  )
  if (!identical(dependencies, expected$dependency_roots)) {
    .release_contract_error("Runtime dependency roots do not match v1.")
  }
  environment <- unlist(contract$runtime$environment, use.names = TRUE)
  if (!identical(environment, expected$environment)) {
    .release_contract_error("Sanitized runtime environment does not match v1.")
  }

  validations <- .release_json_array(
    contract$required_validations, "required_validations"
  )
  if (length(validations) != length(expected$validation_ids)) {
    .release_contract_error("Release Contract has missing or extra validations.")
  }
  for (i in seq_along(validations)) {
    .release_assert_exact_fields(
      validations[[i]],
      c(
        "validation_id", "adapter", "settings", "artifacts",
        "scientific_gates"
      ),
      paste0("Required Validation declaration ", i)
    )
  }
  ids <- vapply(validations, `[[`, character(1L), "validation_id")
  adapters <- vapply(validations, `[[`, character(1L), "adapter")
  if (!identical(ids, expected$validation_ids) || anyDuplicated(ids)) {
    .release_contract_error("Required Validation IDs do not match v1.")
  }
  if (!identical(adapters, expected$adapters)) {
    .release_contract_error("Unknown or replaced Required Validation Adapter.")
  }
  for (i in seq_along(validations)) {
    validation <- validations[[i]]
    id <- ids[[i]]
    .release_validate_settings(validation$settings, expected$settings[[id]], id)
    .release_validate_artifacts(
      validation$artifacts, id, expected$artifact_names[[id]]
    )
    .release_validate_gates(
      validation$scientific_gates, id, expected$gate_ids[[id]]
    )
  }
  invisible(contract)
}

.release_read_contract <- function(path) {
  contract <- .release_read_json(path, "Release Contract")
  tryCatch(
    .release_validate_contract(contract),
    stablr_release_contract_error = function(e) stop(e),
    error = function(e) {
      .release_contract_error(
        "Release Contract validation failed: ", conditionMessage(e), parent = e
      )
    }
  )
  attr(contract, "path") <- normalizePath(path, winslash = "/", mustWork = TRUE)
  attr(contract, "sha256") <- .release_hash_file(path)
  contract
}

.release_validate_context_ref <- function(record) {
  .release_assert_exact_fields(
    record,
    c(
      "schema_version", "kind", "context_id", "store", "context_path",
      "context_sha256"
    ),
    "Release Context Reference",
    condition = "identity"
  )
  if (!identical(record$schema_version, "stablr.release-context-ref/v1") ||
      !identical(record$kind, "release_context") ||
      !.release_is_scalar_character(record$context_id) ||
      nchar(record$context_id) < 16L ||
      !all(vapply(
        record[c("store", "context_path")],
        function(x) .release_is_scalar_character(x) && grepl("^/", x),
        logical(1L)
      )) || !.release_is_sha256(record$context_sha256)) {
    .release_identity_error("Invalid Release Context Reference.")
  }
  invisible(record)
}

.release_validate_packet_ref <- function(record) {
  .release_assert_exact_fields(
    record,
    c(
      "schema_version", "kind", "packet_id", "validation_id", "status",
      "packet_path", "packet_sha256"
    ),
    "Packet Reference",
    condition = "identity"
  )
  valid <- identical(record$schema_version, "stablr.packet-ref/v1") &&
    identical(record$kind, "validation_packet") &&
    .release_is_scalar_character(record$packet_id) &&
    nchar(record$packet_id) >= 16L &&
    record$validation_id %in% c("methodology", "late_fusion") &&
    record$status %in% c("passing", "failed") &&
    .release_is_scalar_character(record$packet_path) &&
    grepl("^/", record$packet_path) &&
    .release_is_sha256(record$packet_sha256)
  if (!isTRUE(valid)) .release_identity_error("Invalid Packet Reference.")
  invisible(record)
}

.release_new_context_ref <- function(record) {
  .release_validate_context_ref(record)
  structure(
    .release_sort_json_object(record),
    class = c("stablr_release_context", "list")
  )
}

.release_new_packet_ref <- function(record) {
  .release_validate_packet_ref(record)
  structure(
    .release_sort_json_object(record),
    class = c("stablr_packet_ref", "list")
  )
}

write_release_reference <- function(reference, path) {
  record <- unclass(reference)
  if (inherits(reference, "stablr_release_context")) {
    .release_validate_context_ref(record)
  } else if (inherits(reference, "stablr_packet_ref")) {
    .release_validate_packet_ref(record)
  } else {
    .release_identity_error("Unsupported release reference class.")
  }
  .release_write_canonical_json(record, path)
  invisible(normalizePath(path, winslash = "/", mustWork = TRUE))
}

read_release_reference <- function(path) {
  record <- .release_read_json(path, "release reference")
  if (identical(record$schema_version, "stablr.release-context-ref/v1")) {
    return(.release_new_context_ref(record))
  }
  if (identical(record$schema_version, "stablr.packet-ref/v1")) {
    return(.release_new_packet_ref(record))
  }
  .release_identity_error("Unknown release reference schema version.")
}
