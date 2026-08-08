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

.release_typed_error <- function(subclass, message, ...) {
  details <- list(...)
  detail_names <- names(details)
  if (is.null(detail_names)) detail_names <- rep("", length(details))
  message_fields <- !nzchar(detail_names)
  message <- paste0(
    c(list(message), details[message_fields]), collapse = ""
  )
  metadata <- details[!message_fields]
  do.call(
    .release_abort,
    c(list(subclass = subclass, message = message), metadata)
  )
}

.release_contract_error <- function(message, ...) {
  .release_typed_error("stablr_release_contract_error", message, ...)
}

.release_identity_error <- function(message, ...) {
  .release_typed_error("stablr_release_identity_error", message, ...)
}

.release_runtime_error <- function(message, ...) {
  .release_typed_error("stablr_release_runtime_error", message, ...)
}

.release_assembly_error <- function(message, ...) {
  .release_typed_error("stablr_release_assembly_error", message, ...)
}

.release_incomplete_error <- function(message, ...) {
  .release_typed_error("stablr_validation_incomplete", message, ...)
}

.release_validation_failed_error <- function(message, packet_reference, ...) {
  .release_typed_error(
    "stablr_validation_failed",
    message,
    packet_reference = packet_reference,
    ...
  )
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
  switch(
    condition,
    contract = .release_contract_error(message, record = label),
    identity = .release_identity_error(message, record = label),
    runtime = .release_runtime_error(message, record = label),
    assembly = .release_assembly_error(message, record = label),
    .release_identity_error(message, record = label)
  )
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

.release_json_identical <- function(x, y) {
  identical(
    as.character(.release_canonical_json(x)),
    as.character(.release_canonical_json(y))
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
    ),
    gate_versions = list(
      methodology = c(
        "methodology-cell-completeness/v1",
        "python-metrics-parity/v1",
        "null-select-any-wilson-v1",
        "SCI-04-null-selected-fraction-v2-draft",
        "null-90pct-collapse-count-v1",
        "signal-mean-fdp-t-v1",
        "signal-tpr-t-v1"
      ),
      late_fusion = c(
        "late-fusion-cell-completeness/v1",
        "late-fusion-reduced-optimism/v1",
        "late-fusion-noninferiority-0.02/v1",
        "late-fusion-signal-fallback-0.05/v1"
      )
    ),
    gate_statuses = list(
      methodology = c(
        rep("accepted", 3L), "proposed_unaccepted", rep("accepted", 3L)
      ),
      late_fusion = rep("accepted", 4L)
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

.release_validate_gates <- function(gates, validation_id, expected_ids,
                                    expected_versions, expected_statuses) {
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
  versions <- character(length(gates))
  statuses <- character(length(gates))
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
    versions[[i]] <- gate$gate_version
    statuses[[i]] <- gate$decision_status
  }
  if (!identical(observed, expected_ids) || anyDuplicated(observed) ||
      !identical(versions, expected_versions) ||
      !identical(statuses, expected_statuses)) {
    .release_contract_error(
      "Required validation `", validation_id,
      "` Scientific Gates do not match the locked v1 contract."
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
      validation$scientific_gates,
      id,
      expected$gate_ids[[id]],
      expected$gate_versions[[id]],
      expected$gate_statuses[[id]]
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

.release_validate_release_ref <- function(record) {
  .release_assert_exact_fields(
    record,
    c(
      "schema_version", "kind", "release_id", "store", "release_path",
      "release_sha256"
    ),
    "Release Evidence Reference",
    condition = "identity"
  )
  valid <- identical(
    record$schema_version, "stablr.release-evidence-ref/v1"
  ) && identical(record$kind, "release_evidence") &&
    .release_is_scalar_character(record$release_id) &&
    nchar(record$release_id) >= 16L &&
    .release_is_scalar_character(record$store) && grepl("^/", record$store) &&
    .release_is_scalar_character(record$release_path) &&
    grepl("^/", record$release_path) &&
    .release_is_sha256(record$release_sha256)
  if (!isTRUE(valid)) {
    .release_identity_error("Invalid Release Evidence Reference.")
  }
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

.release_new_release_ref <- function(record) {
  .release_validate_release_ref(record)
  structure(
    .release_sort_json_object(record),
    class = c("stablr_release_ref", "list")
  )
}

write_release_reference <- function(reference, path) {
  record <- unclass(reference)
  if (inherits(reference, "stablr_release_context")) {
    .release_validate_context_ref(record)
  } else if (inherits(reference, "stablr_packet_ref")) {
    .release_validate_packet_ref(record)
  } else if (inherits(reference, "stablr_release_ref")) {
    .release_validate_release_ref(record)
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
  if (identical(record$schema_version, "stablr.release-evidence-ref/v1")) {
    return(.release_new_release_ref(record))
  }
  .release_identity_error("Unknown release reference schema version.")
}

.release_normalize_absolute <- function(path, label, type = c("file", "dir"),
                                        must_work = TRUE) {
  type <- match.arg(type)
  if (!.release_is_scalar_character(path) || !grepl("^/", path)) {
    .release_identity_error(label, " must be an absolute path.")
  }
  normalized <- tryCatch(
    normalizePath(path, winslash = "/", mustWork = must_work),
    error = function(e) {
      .release_identity_error(
        label, " cannot be resolved: ", conditionMessage(e), parent = e
      )
    }
  )
  if (must_work) {
    valid <- if (identical(type, "file")) {
      file.exists(normalized) && !dir.exists(normalized)
    } else {
      dir.exists(normalized)
    }
    if (!isTRUE(valid)) {
      .release_identity_error(label, " has the wrong filesystem type.")
    }
  }
  normalized
}

.release_path_within <- function(path, root, must_work = TRUE) {
  path <- normalizePath(path, winslash = "/", mustWork = must_work)
  root <- normalizePath(root, winslash = "/", mustWork = must_work)
  identical(path, root) || startsWith(path, paste0(root, "/"))
}

.release_hash_record <- function(record) {
  path <- tempfile("stablr-release-record-", fileext = ".json")
  on.exit(unlink(path), add = TRUE)
  connection <- file(path, open = "wb")
  tryCatch(
    writeChar(
      as.character(.release_canonical_json(record)),
      connection,
      eos = NULL,
      useBytes = TRUE
    ),
    finally = close(connection)
  )
  .release_hash_file(path)
}

.release_tree_record <- function(path, label = "filesystem tree") {
  root <- .release_normalize_absolute(path, label, type = "dir")
  entries <- list.files(
    root,
    all.files = TRUE,
    no.. = TRUE,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = TRUE
  )
  if (!length(entries)) {
    manifest <- list()
    return(list(manifest = manifest, sha256 = .release_hash_record(manifest)))
  }
  links <- Sys.readlink(entries)
  if (any(nzchar(links))) {
    .release_identity_error(label, " contains a symbolic link.")
  }
  info <- file.info(entries)
  if (anyNA(info$isdir)) {
    .release_identity_error(label, " contains an unreadable member.")
  }
  relative <- substring(entries, nchar(root) + 2L)
  order_index <- order(relative, method = "radix")
  entries <- entries[order_index]
  relative <- relative[order_index]
  info <- info[order_index, , drop = FALSE]
  manifest <- lapply(seq_along(entries), function(i) {
    is_directory <- isTRUE(info$isdir[[i]])
    list(
      path = relative[[i]],
      type = if (is_directory) "directory" else "file",
      mode = as.character(as.octmode(info$mode[[i]])),
      size = if (is_directory) 0 else unname(info$size[[i]]),
      sha256 = if (is_directory) NULL else .release_hash_file(entries[[i]])
    )
  })
  list(manifest = manifest, sha256 = .release_hash_record(manifest))
}

.release_apply_tree_modes <- function(root, manifest, label) {
  previous_umask <- Sys.umask("0000")
  on.exit(Sys.umask(previous_umask), add = TRUE)
  for (member in manifest) {
    path <- file.path(root, member$path)
    mode <- paste0("0", member$mode)
    if (!file.exists(path) || !isTRUE(Sys.chmod(path, mode = mode))) {
      .release_runtime_error(label, " file modes could not be preserved.")
    }
  }
  invisible(TRUE)
}

.release_validate_archive_members <- function(members) {
  if (!is.character(members) || !length(members) || anyNA(members) ||
      any(!nzchar(members))) {
    .release_identity_error("Candidate archive has an invalid member list.")
  }
  normalized <- sub("/+$", "", members)
  invalid <- vapply(normalized, function(member) {
    components <- strsplit(member, "/", fixed = TRUE)[[1L]]
    grepl("^/", member) || grepl("^[A-Za-z]:", member) ||
      grepl("[\\\\\r\n]", member) ||
      any(components %in% c("", ".", ".."))
  }, logical(1L))
  if (any(invalid) || anyDuplicated(normalized)) {
    .release_identity_error(
      "Candidate archive contains an unsafe or duplicate member name."
    )
  }
  roots <- unique(vapply(
    strsplit(normalized, "/", fixed = TRUE), `[[`, character(1L), 1L
  ))
  if (length(roots) != 1L) {
    .release_identity_error("Candidate archive must contain exactly one root.")
  }
  list(root = roots[[1L]], members = normalized)
}

.release_inspect_candidate_archive <- function(tarball, destination) {
  members <- tryCatch(
    utils::untar(tarball, list = TRUE),
    error = function(e) {
      .release_identity_error(
        "Candidate archive could not be listed: ", conditionMessage(e),
        parent = e
      )
    }
  )
  inspected <- .release_validate_archive_members(members)
  tar_tool <- Sys.which("tar")
  if (!nzchar(tar_tool)) {
    .release_runtime_error("Candidate inspection requires `tar`.")
  }
  verbose <- suppressWarnings(system2(
    unname(tar_tool),
    c("-tvzf", shQuote(tarball)),
    stdout = TRUE,
    stderr = TRUE
  ))
  if ((attr(verbose, "status") %||% 0L) != 0L) {
    .release_identity_error("Candidate archive verbose inspection failed.")
  }
  if (any(!grepl("^[-d]", verbose))) {
    .release_identity_error(
      "Candidate archive contains a link or unsupported member type."
    )
  }
  if (!dir.create(destination, recursive = TRUE, showWarnings = FALSE)) {
    .release_identity_error("Could not allocate candidate inspection directory.")
  }
  status <- tryCatch(
    utils::untar(tarball, exdir = destination),
    error = function(e) e
  )
  if (inherits(status, "error") || (!is.null(status) && status != 0L)) {
    .release_identity_error("Candidate archive extraction failed.")
  }
  root <- file.path(destination, inspected$root)
  root <- .release_normalize_absolute(root, "candidate archive root", "dir")
  if (!.release_path_within(root, destination)) {
    .release_identity_error("Candidate archive root escapes its extraction.")
  }
  extracted <- list.files(
    root, all.files = TRUE, no.. = TRUE, recursive = TRUE, full.names = TRUE,
    include.dirs = TRUE
  )
  if (length(extracted) && any(nzchar(Sys.readlink(extracted)))) {
    .release_identity_error("Extracted candidate contains a symbolic link.")
  }
  description_path <- file.path(root, "DESCRIPTION")
  contract_path <- file.path(root, "inst", "release", "release-contract.json")
  if (!file.exists(description_path) || !file.exists(contract_path)) {
    .release_identity_error(
      "Candidate archive is missing DESCRIPTION or the Release Contract."
    )
  }
  description <- tryCatch(
    read.dcf(description_path),
    error = function(e) {
      .release_identity_error("Candidate DESCRIPTION cannot be parsed.", parent = e)
    }
  )
  if (!identical(unname(description[1L, "Package"]), "stablr") ||
      !identical(unname(description[1L, "Version"]), "0.1.1")) {
    .release_identity_error("Candidate package identity is not stablr 0.1.1.")
  }
  contract <- .release_read_contract(contract_path)
  list(
    root = root,
    archive_root = inspected$root,
    description = description,
    contract = contract,
    contract_path = contract_path
  )
}

.release_probe_store_atomicity <- function(store) {
  probe_id <- .release_random_id()
  staging <- file.path(store, paste0(".atomic-probe-", probe_id))
  published <- file.path(store, paste0(".atomic-probe-published-", probe_id))
  on.exit(unlink(c(staging, published), recursive = TRUE, force = TRUE), add = TRUE)
  if (!dir.create(staging, showWarnings = FALSE)) {
    .release_runtime_error("Evidence store atomicity probe could not create staging.")
  }
  marker <- file.path(staging, "EXCLUSIVE")
  if (!dir.create(marker, showWarnings = FALSE) ||
      dir.create(marker, showWarnings = FALSE)) {
    .release_runtime_error("Evidence store lacks exclusive marker semantics.")
  }
  if (!file.rename(staging, published)) {
    .release_runtime_error("Evidence store lacks same-filesystem atomic rename.")
  }
  invisible(TRUE)
}

.release_store_candidate <- function(tarball, candidate_sha256, store) {
  candidate_dir <- file.path(store, "candidates", candidate_sha256)
  candidate_path <- file.path(candidate_dir, "source.tar.gz")
  if (dir.exists(candidate_dir)) {
    if (!dir.exists(file.path(candidate_dir, "SEALED")) ||
        !file.exists(candidate_path) ||
        !identical(.release_hash_file(candidate_path), candidate_sha256)) {
      .release_identity_error("Existing store candidate is incomplete or changed.")
    }
    return(normalizePath(candidate_path, winslash = "/", mustWork = TRUE))
  }
  parent <- file.path(store, "candidates")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  staging <- file.path(parent, paste0(".candidate-", .release_random_id()))
  if (!dir.create(staging, showWarnings = FALSE)) {
    .release_identity_error("Could not allocate store candidate staging.")
  }
  completed <- FALSE
  on.exit({
    if (!completed && dir.exists(staging)) {
      unlink(staging, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)
  staged_path <- file.path(staging, "source.tar.gz")
  if (!file.copy(tarball, staged_path, copy.mode = TRUE, copy.date = TRUE) ||
      !identical(.release_hash_file(staged_path), candidate_sha256)) {
    .release_identity_error("Store-owned candidate copy failed verification.")
  }
  if (!dir.create(file.path(staging, "SEALED"), showWarnings = FALSE)) {
    .release_identity_error("Could not seal the store-owned candidate.")
  }
  if (!file.rename(staging, candidate_dir)) {
    if (!dir.exists(candidate_dir) || !file.exists(candidate_path) ||
        !identical(.release_hash_file(candidate_path), candidate_sha256)) {
      .release_identity_error("Candidate publication lost a concurrent race.")
    }
    unlink(staging, recursive = TRUE, force = TRUE)
  }
  completed <- TRUE
  Sys.chmod(candidate_path, mode = "0444")
  normalizePath(candidate_path, winslash = "/", mustWork = TRUE)
}

.release_dependency_closure <- function(roots) {
  database <- installed.packages(noCache = TRUE)
  missing <- setdiff(roots, rownames(database))
  if (length(missing)) {
    .release_runtime_error(
      "Contract-declared dependencies are unavailable: ",
      paste(missing, collapse = ", "),
      missing_packages = missing
    )
  }
  dependencies <- unique(c(
    roots,
    unlist(
      tools::package_dependencies(
        roots,
        db = database,
        which = c("Depends", "Imports", "LinkingTo"),
        recursive = TRUE
      ),
      use.names = FALSE
    )
  ))
  dependencies <- setdiff(dependencies, c("R", "base"))
  missing <- setdiff(dependencies, rownames(database))
  if (length(missing)) {
    .release_runtime_error(
      "The dependency closure is incomplete: ", paste(missing, collapse = ", "),
      missing_packages = missing
    )
  }
  priority <- database[dependencies, "Priority"]
  dependencies <- dependencies[is.na(priority) | priority != "base"]
  sort(unique(dependencies), method = "radix")
}

.release_copy_dependency_closure <- function(roots, library) {
  packages <- .release_dependency_closure(roots)
  records <- vector("list", length(packages))
  for (i in seq_along(packages)) {
    package <- packages[[i]]
    source <- tryCatch(
      find.package(package, quiet = FALSE),
      error = function(e) {
        .release_runtime_error(
          "Could not resolve dependency `", package, "`.", parent = e
        )
      }
    )
    source <- .release_normalize_absolute(
      source, paste0("dependency `", package, "`"), "dir"
    )
    source_before <- .release_tree_record(source, paste0("dependency `", package, "`"))
    destination <- file.path(library, package)
    if (file.exists(destination) ||
        !file.copy(
          source,
          library,
          recursive = TRUE,
          copy.mode = TRUE,
          copy.date = TRUE
        )) {
      .release_runtime_error("Could not copy dependency `", package, "`.")
    }
    .release_apply_tree_modes(
      destination,
      source_before$manifest,
      paste0("dependency `", package, "`")
    )
    source_after <- .release_tree_record(source, paste0("dependency `", package, "`"))
    copied <- .release_tree_record(
      destination, paste0("copied dependency `", package, "`")
    )
    if (!identical(source_before$sha256, source_after$sha256) ||
        !identical(source_before$sha256, copied$sha256)) {
      .release_runtime_error(
        "Dependency `", package, "` changed during isolated copying."
      )
    }
    records[[i]] <- list(
      package = package,
      version = as.character(utils::packageVersion(package)),
      source_path = source,
      tree_sha256 = copied$sha256
    )
  }
  records
}

.release_sanitized_environment <- function(contract, library) {
  environment <- unlist(contract$runtime$environment, use.names = TRUE)
  environment <- c(
    environment,
    R_LIBS = library,
    R_LIBS_USER = library,
    R_LIBS_SITE = "",
    R_DEFAULT_PACKAGES = "NULL"
  )
  paste0(names(environment), "=", unname(environment))
}

.release_run_process <- function(command, args, environment, label) {
  output <- suppressWarnings(system2(
    command,
    args,
    stdout = TRUE,
    stderr = TRUE,
    env = environment
  ))
  status <- attr(output, "status") %||% 0L
  if (status != 0L) {
    .release_runtime_error(
      label, " failed with exit status ", status, ".",
      process_output = output,
      process_status = status
    )
  }
  output
}

.release_install_candidate <- function(candidate_tarball, r_executable,
                                       library, environment, log_path) {
  output <- .release_run_process(
    r_executable,
    c(
      "CMD", "INSTALL", paste0("--library=", shQuote(library)),
      "--no-multiarch", "--no-test-load", "--install-tests",
      shQuote(candidate_tarball)
    ),
    environment,
    "Candidate installation"
  )
  writeLines(output, log_path, useBytes = TRUE)
  candidate_path <- file.path(library, "stablr")
  .release_normalize_absolute(candidate_path, "installed candidate", "dir")
}

.release_write_runtime_probe <- function(path) {
  lines <- c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "library <- normalizePath(args[[1L]], winslash = '/', mustWork = TRUE)",
    "output <- args[[2L]]",
    ".libPaths(c(library, .Library))",
    "if ('stablr' %in% loadedNamespaces()) quit(status = 91L)",
    "namespace <- loadNamespace('stablr')",
    "candidate_path <- normalizePath(getNamespaceInfo(namespace, 'path'), winslash = '/', mustWork = TRUE)",
    "if (!(identical(candidate_path, file.path(library, 'stablr')) || startsWith(candidate_path, paste0(library, '/')))) quit(status = 92L)",
    "session <- utils::sessionInfo()",
    "dlls <- getLoadedDLLs()",
    "dll_paths <- unique(vapply(dlls, function(x) { path <- x[['path']]; if (is.null(path)) '' else path }, character(1L)))",
    "dll_paths <- dll_paths[nzchar(dll_paths) & file.exists(dll_paths)]",
    "record <- list(r_home = normalizePath(R.home(), winslash = '/', mustWork = TRUE), r_version = R.version.string, platform = R.version$platform, locale = list(all = Sys.getlocale(), collate = Sys.getlocale('LC_COLLATE'), ctype = Sys.getlocale('LC_CTYPE'), numeric = Sys.getlocale('LC_NUMERIC'), time = Sys.getlocale('LC_TIME')), blas = session$BLAS, lapack = session$LAPACK, candidate_path = candidate_path, candidate_version = as.character(utils::packageVersion('stablr')), dll_paths = as.list(sort(dll_paths)))",
    "jsonlite::write_json(record, output, auto_unbox = TRUE, null = 'null', pretty = FALSE)"
  )
  writeLines(lines, path, useBytes = TRUE)
}

.release_ldd_paths <- function(objects) {
  tool <- Sys.which("ldd")
  if (!nzchar(tool)) {
    .release_runtime_error("Runtime Artifact discovery requires `ldd`.")
  }
  paths <- character()
  for (object in objects) {
    output <- suppressWarnings(system2(
      unname(tool), shQuote(object), stdout = TRUE, stderr = TRUE
    ))
    if ((attr(output, "status") %||% 0L) != 0L) {
      .release_runtime_error("Native-library discovery failed for `", object, "`.")
    }
    for (line in trimws(output)) {
      candidate <- if (grepl("=>", line, fixed = TRUE)) {
        trimws(strsplit(line, "=>", fixed = TRUE)[[1L]][[2L]])
      } else {
        line
      }
      candidate <- strsplit(candidate, "[[:space:](]")[[1L]][[1L]]
      if (startsWith(candidate, "/") && file.exists(candidate)) {
        paths <- c(paths, normalizePath(candidate, winslash = "/"))
      }
    }
  }
  sort(unique(paths), method = "radix")
}

.release_build_runtime_artifact <- function(runtime, r_executable, library,
                                            package_records, environment,
                                            staging, final_context) {
  probe_script <- file.path(staging, "runtime-probe.R")
  probe_output <- file.path(staging, "runtime-probe.json")
  .release_write_runtime_probe(probe_script)
  .release_run_process(
    runtime,
    c("--vanilla", shQuote(probe_script), shQuote(library), shQuote(probe_output)),
    environment,
    "Runtime Artifact probe"
  )
  probe <- .release_read_json(probe_output, "Runtime Artifact probe output")
  if (!identical(probe$candidate_version, "0.1.1")) {
    .release_runtime_error("Isolated candidate version is not 0.1.1.")
  }
  expected_candidate <- file.path(final_context, "library", "stablr")
  staged_candidate <- file.path(staging, "library", "stablr")
  if (!identical(
    normalizePath(probe$candidate_path, winslash = "/", mustWork = TRUE),
    normalizePath(staged_candidate, winslash = "/", mustWork = TRUE)
  )) {
    .release_runtime_error("The loaded stablr namespace is not the candidate.")
  }
  package_records <- lapply(package_records, function(record) {
    installed <- file.path(library, record$package)
    record$installed_path <- file.path(
      final_context, "library", record$package
    )
    record$tree_sha256 <- .release_tree_record(
      installed, paste0("isolated dependency `", record$package, "`")
    )$sha256
    record
  })
  candidate_tree <- .release_tree_record(staged_candidate, "installed candidate")
  package_records <- c(
    package_records,
    list(list(
      package = "stablr",
      version = "0.1.1",
      source_path = expected_candidate,
      installed_path = expected_candidate,
      tree_sha256 = candidate_tree$sha256
    ))
  )
  package_records <- package_records[order(vapply(
    package_records, `[[`, character(1L), "package"
  ), method = "radix")]
  library_tree <- .release_tree_record(library, "isolated candidate library")
  objects <- unique(c(
    unlist(probe$dll_paths, use.names = FALSE),
    list.files(library, pattern = "[.]so$", recursive = TRUE, full.names = TRUE)
  ))
  objects <- objects[file.exists(objects)]
  native_paths <- sort(unique(c(objects, .release_ldd_paths(objects))),
                       method = "radix")
  native_libraries <- lapply(native_paths, function(path) {
    path <- normalizePath(path, winslash = "/", mustWork = TRUE)
    recorded_path <- if (.release_path_within(path, library)) {
      relative <- substring(path, nchar(normalizePath(
        library, winslash = "/", mustWork = TRUE
      )) + 2L)
      file.path(final_context, "library", relative)
    } else {
      path
    }
    list(path = recorded_path, sha256 = .release_hash_file(path))
  })
  base_record <- list(
    schema_version = "stablr.runtime-artifact/v1",
    rscript = list(path = runtime, sha256 = .release_hash_file(runtime)),
    r_executable = list(
      path = r_executable,
      sha256 = .release_hash_file(r_executable)
    ),
    r_home = probe$r_home,
    r_version = probe$r_version,
    platform = probe$platform,
    locale = probe$locale,
    blas = probe$blas,
    lapack = probe$lapack,
    library = list(path = file.path(final_context, "library"),
                   tree_sha256 = library_tree$sha256),
    packages = package_records,
    native_libraries = native_libraries
  )
  runtime_id <- .release_hash_record(base_record)
  c(list(runtime_id = runtime_id), base_record)
}

.release_validate_runtime_artifact <- function(record) {
  .release_assert_exact_fields(
    record,
    c(
      "runtime_id", "schema_version", "rscript", "r_executable", "r_home",
      "r_version", "platform", "locale", "blas", "lapack", "library",
      "packages", "native_libraries"
    ),
    "Runtime Artifact",
    condition = "runtime"
  )
  for (name in c("rscript", "r_executable")) {
    .release_assert_exact_fields(
      record[[name]], c("path", "sha256"),
      paste0("Runtime Artifact ", name), condition = "runtime"
    )
  }
  .release_assert_exact_fields(
    record$library, c("path", "tree_sha256"),
    "Runtime Artifact library", condition = "runtime"
  )
  .release_assert_exact_fields(
    record$locale, c("all", "collate", "ctype", "numeric", "time"),
    "Runtime Artifact locale", condition = "runtime"
  )
  path_records <- c(record[c("rscript", "r_executable")], list(record$library))
  path_valid <- vapply(path_records, function(member) {
    .release_is_scalar_character(member$path) && grepl("^/", member$path)
  }, logical(1L))
  hashes <- c(
    record$rscript$sha256,
    record$r_executable$sha256,
    record$library$tree_sha256
  )
  if (!all(path_valid) ||
      !all(vapply(hashes, .release_is_sha256, logical(1L))) ||
      !all(vapply(record$locale, .release_is_scalar_character, logical(1L))) ||
      !.release_is_scalar_character(record$r_home) ||
      !grepl("^/", record$r_home) ||
      !.release_is_scalar_character(record$r_version) ||
      !.release_is_scalar_character(record$platform)) {
    .release_runtime_error("Runtime Artifact scalar identity is invalid.")
  }
  if (!is.list(record$packages) || is.object(record$packages) ||
      !length(record$packages)) {
    .release_runtime_error("Runtime Artifact package closure is invalid.")
  }
  package_names <- character(length(record$packages))
  installed_paths <- character(length(record$packages))
  for (i in seq_along(record$packages)) {
    package <- record$packages[[i]]
    .release_assert_exact_fields(
      package,
      c(
        "package", "version", "source_path", "installed_path",
        "tree_sha256"
      ),
      paste0("Runtime Artifact package record ", i),
      condition = "runtime"
    )
    scalar_fields <- package[c(
      "package", "version", "source_path", "installed_path"
    )]
    if (!all(vapply(
      scalar_fields, .release_is_scalar_character, logical(1L)
    )) || !grepl("^/", package$source_path) ||
        !grepl("^/", package$installed_path) ||
        !.release_is_sha256(package$tree_sha256)) {
      .release_runtime_error("Runtime Artifact package identity is invalid.")
    }
    package_names[[i]] <- package$package
    installed_paths[[i]] <- package$installed_path
  }
  if (anyDuplicated(package_names) || anyDuplicated(installed_paths) ||
      sum(package_names == "stablr") != 1L) {
    .release_runtime_error("Runtime Artifact package closure is ambiguous.")
  }
  if (!is.list(record$native_libraries) ||
      is.object(record$native_libraries)) {
    .release_runtime_error("Runtime Artifact native-library closure is invalid.")
  }
  native_paths <- character(length(record$native_libraries))
  for (i in seq_along(record$native_libraries)) {
    native <- record$native_libraries[[i]]
    .release_assert_exact_fields(
      native, c("path", "sha256"),
      paste0("Runtime Artifact native-library record ", i),
      condition = "runtime"
    )
    if (!.release_is_scalar_character(native$path) ||
        !grepl("^/", native$path) || !.release_is_sha256(native$sha256)) {
      .release_runtime_error("Runtime Artifact native-library identity is invalid.")
    }
    native_paths[[i]] <- native$path
  }
  if (anyDuplicated(native_paths)) {
    .release_runtime_error("Runtime Artifact native-library paths are duplicated.")
  }
  base <- record[names(record) != "runtime_id"]
  valid <- identical(record$schema_version, "stablr.runtime-artifact/v1") &&
    .release_is_sha256(record$runtime_id) &&
    identical(.release_hash_record(base), record$runtime_id) &&
    identical(record$platform, R.version$platform)
  if (!isTRUE(valid)) {
    .release_runtime_error("Runtime Artifact record is invalid or unsupported.")
  }
  invisible(record)
}

.release_validate_context_record <- function(record) {
  .release_assert_exact_fields(
    record,
    c(
      "schema_version", "context_id", "store", "candidate", "contract",
      "runtime", "library", "required_validations"
    ),
    "Release Context",
    condition = "identity"
  )
  if (!identical(record$schema_version, "stablr.release-context/v1") ||
      !.release_is_scalar_character(record$context_id) ||
      nchar(record$context_id) < 16L ||
      !.release_is_scalar_character(record$store) ||
      !grepl("^/", record$store)) {
    .release_identity_error("Release Context record is invalid.")
  }
  .release_assert_exact_fields(
    record$candidate,
    c("path", "sha256", "package", "version", "archive_root"),
    "Release Context candidate",
    condition = "identity"
  )
  .release_assert_exact_fields(
    record$contract,
    c("path", "sha256", "contract_id", "promotable", "blockers"),
    "Release Context contract",
    condition = "identity"
  )
  .release_assert_exact_fields(
    record$runtime, c("path", "sha256", "runtime_id"),
    "Release Context runtime", condition = "identity"
  )
  .release_assert_exact_fields(
    record$library, c("path", "tree_sha256"),
    "Release Context library", condition = "identity"
  )
  hashes <- c(
    record$candidate$sha256,
    record$contract$sha256,
    record$runtime$sha256,
    record$runtime$runtime_id,
    record$library$tree_sha256
  )
  path_values <- c(
    record$candidate$path,
    record$contract$path,
    record$runtime$path,
    record$library$path
  )
  if (!all(vapply(hashes, .release_is_sha256, logical(1L))) ||
      !all(vapply(
        path_values,
        function(path) .release_is_scalar_character(path) && grepl("^/", path),
        logical(1L)
      )) ||
      !identical(record$candidate$package, "stablr") ||
      !identical(record$candidate$version, "0.1.1") ||
      !.release_is_scalar_character(record$candidate$archive_root) ||
      !.release_is_scalar_character(record$contract$contract_id) ||
      !.release_is_scalar_logical(record$contract$promotable) ||
      !is.list(record$contract$blockers) ||
      !is.list(record$required_validations)) {
    .release_identity_error("Release Context identities are invalid.")
  }
  invisible(record)
}

.release_make_read_only <- function(path, include_root = TRUE) {
  entries <- c(
    list.files(
      path, all.files = TRUE, no.. = TRUE, recursive = TRUE,
      full.names = TRUE, include.dirs = TRUE
    ),
    if (isTRUE(include_root)) path else character()
  )
  info <- file.info(entries)
  files <- entries[!info$isdir]
  directories <- entries[info$isdir]
  if (length(files)) Sys.chmod(files, mode = "0444")
  if (length(directories)) Sys.chmod(directories, mode = "0555")
  invisible(TRUE)
}

.release_discard_attempt <- function(path) {
  if (!dir.exists(path)) return(invisible(TRUE))
  entries <- list.files(
    path, all.files = TRUE, no.. = TRUE, recursive = TRUE,
    full.names = TRUE, include.dirs = TRUE
  )
  if (length(entries)) Sys.chmod(entries, mode = "0700")
  Sys.chmod(path, mode = "0700")
  unlink(path, recursive = TRUE, force = TRUE)
  invisible(!dir.exists(path))
}

prepare_release_context <- function(candidate_tarball, runtime, store) {
  if ("stablr" %in% loadedNamespaces()) {
    .release_identity_error(
      "Preparation requires a fresh R process with no preloaded stablr namespace."
    )
  }
  candidate_tarball <- .release_normalize_absolute(
    candidate_tarball, "candidate_tarball", "file"
  )
  runtime <- .release_normalize_absolute(runtime, "runtime", "file")
  if (file.access(runtime, mode = 1L) != 0L) {
    .release_runtime_error("`runtime` is not executable.")
  }
  active_runtime <- normalizePath(
    file.path(R.home("bin"), "Rscript"), winslash = "/", mustWork = TRUE
  )
  if (!identical(runtime, active_runtime)) {
    .release_runtime_error(
      "Preparation must execute under the exact supplied Rscript runtime."
    )
  }
  r_executable <- .release_normalize_absolute(
    file.path(dirname(runtime), "R"), "companion R executable", "file"
  )
  if (file.access(r_executable, mode = 1L) != 0L) {
    .release_runtime_error("The companion R executable is not executable.")
  }
  if (!dir.exists(store) &&
      !dir.create(store, recursive = TRUE, showWarnings = FALSE)) {
    .release_identity_error("Could not create the evidence store.")
  }
  store <- .release_normalize_absolute(store, "store", "dir")
  .release_probe_store_atomicity(store)

  caller_hash_before <- .release_hash_file(candidate_tarball)
  stored_tarball <- .release_store_candidate(
    candidate_tarball, caller_hash_before, store
  )
  caller_hash_after <- .release_hash_file(candidate_tarball)
  if (!identical(caller_hash_before, caller_hash_after) ||
      !identical(.release_hash_file(stored_tarball), caller_hash_before)) {
    .release_identity_error("Candidate tarball changed during store adoption.")
  }

  contexts_dir <- file.path(store, "contexts")
  dir.create(contexts_dir, recursive = TRUE, showWarnings = FALSE)
  context_id <- .release_random_id()
  staging <- file.path(contexts_dir, paste0(".attempt-", context_id))
  final_context <- file.path(contexts_dir, paste0("context-", context_id))
  if (!dir.create(staging, showWarnings = FALSE)) {
    .release_identity_error("Could not allocate Release Context staging.")
  }
  published <- FALSE
  on.exit({
    if (!published && dir.exists(staging)) {
      .release_discard_attempt(staging)
    }
  }, add = TRUE)

  inspection <- .release_inspect_candidate_archive(
    stored_tarball, file.path(staging, "archive-inspection")
  )
  contract <- inspection$contract
  contract_dir <- file.path(staging, "contract")
  library <- file.path(staging, "library")
  dir.create(contract_dir, showWarnings = FALSE)
  dir.create(library, showWarnings = FALSE)
  staged_contract <- file.path(contract_dir, "release-contract.json")
  if (!file.copy(inspection$contract_path, staged_contract) ||
      !identical(
        .release_hash_file(staged_contract),
        attr(contract, "sha256")
      )) {
    .release_identity_error("Embedded Release Contract copy failed verification.")
  }
  package_records <- .release_copy_dependency_closure(
    unlist(contract$runtime$dependency_roots, use.names = FALSE),
    library
  )
  environment <- .release_sanitized_environment(contract, library)
  candidate_path <- .release_install_candidate(
    stored_tarball,
    r_executable,
    library,
    environment,
    file.path(staging, "install.log")
  )
  installed_contract <- file.path(
    candidate_path, "release", "release-contract.json"
  )
  if (!file.exists(installed_contract) ||
      !identical(
        .release_hash_file(installed_contract),
        attr(contract, "sha256")
      )) {
    .release_identity_error(
      "Installed candidate Release Contract differs from the source artifact."
    )
  }
  .release_make_read_only(library)
  runtime_artifact <- .release_build_runtime_artifact(
    runtime,
    r_executable,
    library,
    package_records,
    environment,
    staging,
    final_context
  )
  .release_validate_runtime_artifact(runtime_artifact)
  runtime_path <- file.path(staging, "runtime-artifact.json")
  runtime_sha256 <- .release_write_canonical_json(runtime_artifact, runtime_path)
  library_tree <- .release_tree_record(library, "isolated candidate library")
  context_record <- list(
    schema_version = "stablr.release-context/v1",
    context_id = context_id,
    store = store,
    candidate = list(
      path = stored_tarball,
      sha256 = caller_hash_before,
      package = "stablr",
      version = "0.1.1",
      archive_root = inspection$archive_root
    ),
    contract = list(
      path = file.path(final_context, "contract", "release-contract.json"),
      sha256 = attr(contract, "sha256"),
      contract_id = contract$contract_id,
      promotable = contract$promotable,
      blockers = contract$blockers
    ),
    runtime = list(
      path = file.path(final_context, "runtime-artifact.json"),
      sha256 = runtime_sha256,
      runtime_id = runtime_artifact$runtime_id
    ),
    library = list(
      path = file.path(final_context, "library"),
      tree_sha256 = library_tree$sha256
    ),
    required_validations = contract$required_validations
  )
  .release_validate_context_record(context_record)
  context_path <- file.path(staging, "context.json")
  context_sha256 <- .release_write_canonical_json(context_record, context_path)
  unlink(file.path(staging, "archive-inspection"), recursive = TRUE, force = TRUE)
  unlink(file.path(staging, "runtime-probe.R"), force = TRUE)
  unlink(file.path(staging, "runtime-probe.json"), force = TRUE)
  if (!file.rename(staging, final_context)) {
    .release_identity_error("Atomic Release Context publication failed.")
  }
  published <- TRUE
  reference <- .release_new_context_ref(list(
    schema_version = "stablr.release-context-ref/v1",
    kind = "release_context",
    context_id = context_id,
    store = store,
    context_path = file.path(final_context, "context.json"),
    context_sha256 = context_sha256
  ))
  .release_make_read_only(final_context, include_root = FALSE)
  .release_verify_context(reference, require_seal = FALSE)
  if (!dir.create(file.path(final_context, "SEALED"), showWarnings = FALSE)) {
    .release_identity_error("Release Context seal-last publication failed.")
  }
  .release_verify_context(reference)
  reference
}

.release_verify_context <- function(context, require_seal = TRUE) {
  if (is.character(context) && length(context) == 1L) {
    context <- read_release_reference(context)
  }
  if (!inherits(context, "stablr_release_context")) {
    .release_identity_error("`context` is not a Release Context Reference.")
  }
  reference <- unclass(context)
  .release_validate_context_ref(reference)
  store <- .release_normalize_absolute(reference$store, "evidence store", "dir")
  if (!identical(store, reference$store)) {
    .release_identity_error("Release Context Reference store is not canonical.")
  }
  context_path <- .release_normalize_absolute(
    reference$context_path, "Release Context record", "file"
  )
  context_dir <- dirname(context_path)
  expected_context_dir <- file.path(
    store, "contexts", paste0("context-", reference$context_id)
  )
  expected_context_path <- file.path(expected_context_dir, "context.json")
  seal_exists <- dir.exists(file.path(context_dir, "SEALED"))
  if (!identical(context_dir, expected_context_dir) ||
      !identical(context_path, expected_context_path) ||
      (isTRUE(require_seal) && !seal_exists) ||
      (!isTRUE(require_seal) && seal_exists) ||
      !identical(.release_hash_file(context_path), reference$context_sha256)) {
    .release_identity_error("Release Context seal or record identity failed.")
  }
  record <- .release_read_json(context_path, "Release Context")
  .release_validate_context_record(record)
  if (!identical(record$context_id, reference$context_id) ||
      !identical(record$store, store)) {
    .release_identity_error("Release Context Reference targets foreign state.")
  }
  expected_candidate <- file.path(
    store, "candidates", record$candidate$sha256, "source.tar.gz"
  )
  expected_contract <- file.path(
    context_dir, "contract", "release-contract.json"
  )
  expected_runtime <- file.path(context_dir, "runtime-artifact.json")
  expected_library <- file.path(context_dir, "library")
  if (!identical(record$candidate$path, expected_candidate) ||
      !identical(record$contract$path, expected_contract) ||
      !identical(record$runtime$path, expected_runtime) ||
      !identical(record$library$path, expected_library)) {
    .release_identity_error("Release Context contains a foreign path.")
  }
  if (!identical(
    .release_hash_file(record$candidate$path), record$candidate$sha256
  )) {
    .release_identity_error("Candidate Source Artifact identity changed.")
  }
  contract <- .release_read_contract(record$contract$path)
  if (!identical(attr(contract, "sha256"), record$contract$sha256) ||
      !identical(contract$contract_id, record$contract$contract_id) ||
      !identical(contract$promotable, record$contract$promotable) ||
      !.release_json_identical(contract$blockers, record$contract$blockers) ||
      !.release_json_identical(
        contract$required_validations, record$required_validations
      )) {
    .release_identity_error("Release Contract identity changed.")
  }
  if (!identical(
    .release_hash_file(record$runtime$path), record$runtime$sha256
  )) {
    .release_runtime_error("Runtime Artifact record changed.")
  }
  runtime <- .release_read_json(record$runtime$path, "Runtime Artifact")
  .release_validate_runtime_artifact(runtime)
  if (!identical(runtime$runtime_id, record$runtime$runtime_id) ||
      !identical(runtime$library$path, record$library$path) ||
      !identical(
        runtime$library$tree_sha256, record$library$tree_sha256
      ) || !identical(
        runtime$r_home,
        normalizePath(R.home(), winslash = "/", mustWork = TRUE)
      ) || !identical(runtime$r_version, R.version.string) ||
      !identical(runtime$platform, R.version$platform)) {
    .release_runtime_error("Runtime Artifact identity changed.")
  }
  if (!identical(
    .release_hash_file(runtime$rscript$path), runtime$rscript$sha256
  ) || !identical(
    .release_hash_file(runtime$r_executable$path),
    runtime$r_executable$sha256
  )) {
    .release_runtime_error("R executable identity changed.")
  }
  native_ok <- vapply(runtime$native_libraries, function(member) {
    file.exists(member$path) &&
      identical(.release_hash_file(member$path), member$sha256)
  }, logical(1L))
  if (length(native_ok) && !all(native_ok)) {
    .release_runtime_error("Resolved native-library identity changed.")
  }
  package_names <- vapply(
    runtime$packages, `[[`, character(1L), "package"
  )
  expected_package_paths <- file.path(record$library$path, package_names)
  package_ok <- vapply(seq_along(runtime$packages), function(i) {
    package <- runtime$packages[[i]]
    identical(package$installed_path, expected_package_paths[[i]]) &&
      dir.exists(package$installed_path) &&
      identical(
        .release_tree_record(
          package$installed_path,
          paste0("isolated package `", package$package, "`")
        )$sha256,
        package$tree_sha256
      )
  }, logical(1L))
  installed_members <- sort(
    list.dirs(record$library$path, recursive = FALSE, full.names = FALSE),
    method = "radix"
  )
  if (!all(package_ok) ||
      !identical(sort(package_names, method = "radix"), installed_members)) {
    .release_runtime_error("Isolated package closure identity changed.")
  }
  library_tree <- .release_tree_record(
    record$library$path, "isolated candidate library"
  )
  if (!identical(library_tree$sha256, record$library$tree_sha256)) {
    .release_runtime_error("Isolated candidate library identity changed.")
  }
  list(reference = context, record = record, contract = contract,
       runtime = runtime)
}

.release_validation_declaration <- function(contract, validation_id) {
  if (!.release_is_scalar_character(validation_id)) {
    .release_contract_error("`validation_id` must be a non-empty string.")
  }
  ids <- vapply(
    contract$required_validations, `[[`, character(1L), "validation_id"
  )
  index <- match(validation_id, ids)
  if (is.na(index)) {
    .release_contract_error(
      "Unknown Required Validation `", validation_id, "`."
    )
  }
  contract$required_validations[[index]]
}

.release_gate_records_from_data_frame <- function(gates) {
  required <- c(
    "gate_id", "gate_version", "scope", "observed", "criterion", "pass",
    "reason"
  )
  if (!is.data.frame(gates) || !identical(names(gates), required)) {
    stop("Adapter common gate table has an invalid schema.", call. = FALSE)
  }
  lapply(seq_len(nrow(gates)), function(i) {
    list(
      gate_id = as.character(gates$gate_id[[i]]),
      gate_version = as.character(gates$gate_version[[i]]),
      scope = as.character(gates$scope[[i]]),
      observed = as.character(gates$observed[[i]]),
      criterion = as.character(gates$criterion[[i]]),
      pass = as.logical(gates$pass[[i]]),
      reason = as.character(gates$reason[[i]])
    )
  })
}

.release_adapter_artifact_records <- function(declaration, artifacts,
                                               artifact_dir) {
  expected_names <- vapply(
    declaration$artifacts, `[[`, character(1L), "name"
  )
  if (!is.list(artifacts) || is.null(names(artifacts)) ||
      !identical(names(artifacts), expected_names)) {
    stop("Adapter returned a foreign artifact set.", call. = FALSE)
  }
  expected_files <- vapply(
    declaration$artifacts, `[[`, character(1L), "path"
  )
  actual_files <- sort(
    list.files(artifact_dir, all.files = TRUE, no.. = TRUE),
    method = "radix"
  )
  if (!identical(actual_files, sort(expected_files, method = "radix"))) {
    stop("Adapter wrote an undeclared or incomplete artifact set.",
         call. = FALSE)
  }
  lapply(seq_along(declaration$artifacts), function(i) {
    declared <- declaration$artifacts[[i]]
    actual <- normalizePath(
      artifacts[[declared$name]], winslash = "/", mustWork = TRUE
    )
    expected <- normalizePath(
      file.path(artifact_dir, declared$path),
      winslash = "/",
      mustWork = TRUE
    )
    if (!identical(actual, expected) || dir.exists(actual) ||
        nzchar(Sys.readlink(actual))) {
      stop("Adapter artifact path or type changed.", call. = FALSE)
    }
    list(
      name = declared$name,
      path = file.path("artifacts", declared$path),
      schema = declared$schema
    )
  })
}

.release_adapter_child_main <- function(context_reference_path, validation_id,
                                        attempt, result_path) {
  context_reference_path <- .release_normalize_absolute(
    context_reference_path, "child Release Context Reference", "file"
  )
  attempt <- .release_normalize_absolute(
    attempt, "Validation Attempt", "dir"
  )
  if (!.release_is_scalar_character(result_path) ||
      !identical(dirname(result_path), attempt) ||
      !identical(basename(result_path), "adapter-result.json") ||
      file.exists(result_path)) {
    .release_identity_error("Adapter result destination is invalid.")
  }
  context <- read_release_reference(context_reference_path)
  verified <- .release_verify_context(context)
  if (!identical(dirname(attempt), file.path(
    verified$record$store, "attempts"
  )) || !grepl("^[.]attempt-[0-9a-f]{32}$", basename(attempt))) {
    .release_identity_error("Validation Attempt is outside its owned store path.")
  }
  declaration <- .release_validation_declaration(
    verified$contract, validation_id
  )
  candidate_path <- file.path(verified$record$library$path, "stablr")
  loaded_path <- normalizePath(
    getNamespaceInfo(asNamespace("stablr"), "path"),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(loaded_path, candidate_path)) {
    .release_identity_error("Adapter child loaded a foreign stablr namespace.")
  }
  artifact_dir <- file.path(attempt, "artifacts")
  if (!dir.create(artifact_dir, showWarnings = FALSE)) {
    .release_identity_error("Adapter artifact directory could not be created.")
  }
  settings <- declaration$settings

  if (identical(declaration$adapter, "methodology_validation")) {
    adapter <- new.env(parent = globalenv())
    sys.source(
      file.path(candidate_path, "analysis", "run_methodology_validation.R"),
      envir = adapter
    )
    artifacts <- adapter$run_methodology_validation(
      out = artifact_dir,
      profile = settings$profile,
      replicates = as.integer(settings$replicates),
      n_bootstraps = as.integer(settings$n_bootstraps),
      n_lambda = as.integer(settings$n_lambda),
      families = unlist(settings$families, use.names = FALSE),
      artificial_types = unlist(
        settings$artificial_types, use.names = FALSE
      ),
      scenario_ids = unlist(settings$scenario_ids, use.names = FALSE),
      seed = as.integer(settings$seed),
      target_fdp = as.numeric(settings$target_fdp),
      workers = as.integer(settings$workers),
      candidate_only = TRUE,
      candidate_path = candidate_path
    )
    details <- utils::read.csv(
      artifacts$gates, stringsAsFactors = FALSE, check.names = FALSE
    )
    parity <- utils::read.csv(
      artifacts$parity, stringsAsFactors = FALSE, check.names = FALSE
    )
    gates <- adapter$.methodology_common_gate_table(details, parity)
  } else if (identical(declaration$adapter, "late_fusion_validation")) {
    adapter <- new.env(parent = globalenv())
    sys.source(
      file.path(candidate_path, "analysis", "run_late_fusion_validation.R"),
      envir = adapter
    )
    artifacts <- adapter$run_late_fusion_validation(
      out = artifact_dir,
      replicates = as.integer(settings$replicates),
      n_bootstraps = as.integer(settings$n_bootstraps),
      n_iter = as.integer(settings$n_iter),
      seed = as.integer(settings$seed),
      candidate_only = TRUE,
      fail_on_gates = FALSE,
      candidate_path = candidate_path
    )
    gates <- utils::read.csv(
      artifacts$gates, stringsAsFactors = FALSE, check.names = FALSE
    )
  } else {
    .release_contract_error("Required Validation Adapter is unsupported.")
  }

  artifact_records <- .release_adapter_artifact_records(
    declaration, artifacts, artifact_dir
  )
  gate_records <- .release_gate_records_from_data_frame(gates)
  result <- list(
    schema_version = "stablr.adapter-result/v1",
    validation_id = validation_id,
    adapter = declaration$adapter,
    bounded = FALSE,
    artifacts = artifact_records,
    gates = gate_records,
    metadata = list(
      candidate_path = candidate_path,
      candidate_version = as.character(utils::packageVersion("stablr")),
      settings_sha256 = .release_hash_record(settings)
    )
  )
  .release_write_canonical_json(result, result_path)
  invisible(result_path)
}

.release_csv_schema_columns <- function(schema) {
  switch(
    schema,
    "methodology-replicates/v2" = c(
      "family", "scenario", "regime", "correlation", "profile", "replicate",
      "artificial_type", "actual_artificial_type", "actual_artificial_types",
      "fallback_event_count", "fallback_random_permutation_events",
      "fallback_equi_events", "fallback_history", "data_seed", "fit_seed",
      "status", "n", "p", "n_signal", "n_bootstraps", "n_lambda",
      "fallback_random_permutation_warnings", "fallback_equi_warnings",
      "warning_count", "elapsed_sec", "n_selected", "true_positives",
      "false_positives", "empirical_fdp", "tpr", "min_fdp_plus",
      "fdp_threshold", "mean_max_real_score", "mean_max_artificial_score",
      "max_artificial_score", "selected_features", "error"
    ),
    "methodology-summary/v2" = c(
      "family", "scenario", "artificial_type", "actual_artificial_types",
      "generator_identity_complete", "profile", "regime", "correlation",
      "n", "p", "n_signal", "replicates", "ok_replicates",
      "mean_selected", "sd_selected", "se_selected", "mean_empirical_fdp",
      "sd_empirical_fdp", "se_empirical_fdp", "mean_tpr", "sd_tpr",
      "se_tpr", "empirical_fdp_exceedance_rate",
      "sd_empirical_fdp_exceedance", "se_empirical_fdp_exceedance_rate",
      "mean_min_fdp_plus", "mean_fdp_threshold",
      "fallback_random_permutation_rate",
      "se_fallback_random_permutation_rate", "fallback_equi_rate",
      "se_fallback_equi_rate", "mean_elapsed_sec"
    ),
    "methodology-warnings/v2" = c(
      "family", "scenario", "replicate", "artificial_type", "warning_index",
      "warning"
    ),
    "python-metrics-parity/v1" = c(
      "metric", "index", "reference", "observed", "abs_error", "status"
    ),
    "methodology-gates/v2" = c(
      "family", "scenario", "regime", "profile", "artificial_type",
      "actual_artificial_type", "expected_replicates", "observed_replicates",
      "ok_replicates", "cell_status", "gate", "bound", "criterion",
      "gate_version", "decision_status", "pass"
    ),
    "late-fusion-replicates/v1" = c(
      "family", "regime", "replicate", "seed", "data_seed",
      "simulation_attempt", "status", "legacy_train", "oof_train",
      "legacy_test", "oof_test", "legacy_optimism", "oof_optimism",
      "oof_fallback_rate", "error"
    ),
    "late-fusion-summary/v1" = c(
      "family", "regime", "replicates", "successful_replicates",
      "mean_legacy_optimism", "mean_oof_optimism", "mean_test_difference",
      "fallback_rate", "reduced_optimism", "noninferior", "fallback_ok"
    ),
    "late-fusion-warnings/v1" = c(
      "family", "regime", "replicate", "seed", "mode", "warning"
    ),
    "common-gates/v1" = c(
      "gate_id", "gate_version", "scope", "observed", "criterion", "pass",
      "reason"
    ),
    NULL
  )
}

.release_validate_csv_artifact <- function(path, schema) {
  columns <- .release_csv_schema_columns(schema)
  if (is.null(columns)) {
    return(list(valid = FALSE, reason = paste0("unknown schema `", schema, "`")))
  }
  parsed <- tryCatch(
    utils::read.csv(
      path, stringsAsFactors = FALSE, check.names = FALSE,
      na.strings = c("NA")
    ),
    error = identity
  )
  if (inherits(parsed, "error")) {
    return(list(
      valid = FALSE,
      reason = paste0("CSV parse failed: ", conditionMessage(parsed))
    ))
  }
  if (!identical(names(parsed), columns)) {
    return(list(valid = FALSE, reason = "CSV columns do not match the schema"))
  }
  allow_empty <- schema %in% c(
    "methodology-warnings/v2", "late-fusion-warnings/v1"
  )
  if (!allow_empty && nrow(parsed) == 0L) {
    return(list(valid = FALSE, reason = "required CSV has no data rows"))
  }
  if (schema %in% c("methodology-gates/v2", "common-gates/v1") &&
      (anyNA(parsed$pass) || !is.logical(parsed$pass))) {
    return(list(valid = FALSE, reason = "gate pass values are not complete logicals"))
  }
  list(valid = TRUE, reason = "schema validated")
}

.release_validate_adapter_result <- function(result, declaration,
                                             candidate_path) {
  .release_assert_exact_fields(
    result,
    c(
      "schema_version", "validation_id", "adapter", "bounded", "artifacts",
      "gates", "metadata"
    ),
    "Adapter Result",
    condition = "identity"
  )
  .release_assert_exact_fields(
    result$metadata,
    c("candidate_path", "candidate_version", "settings_sha256"),
    "Adapter Result metadata",
    condition = "identity"
  )
  if (!identical(result$schema_version, "stablr.adapter-result/v1") ||
      !identical(result$validation_id, declaration$validation_id) ||
      !identical(result$adapter, declaration$adapter) ||
      !identical(result$bounded, FALSE) ||
      !identical(result$metadata$candidate_path, candidate_path) ||
      !identical(result$metadata$candidate_version, "0.1.1") ||
      !.release_is_sha256(result$metadata$settings_sha256) ||
      !identical(
        result$metadata$settings_sha256,
        .release_hash_record(declaration$settings)
      ) || !is.list(result$artifacts) || !is.list(result$gates)) {
    .release_identity_error("Adapter Result identity is invalid.")
  }
  invisible(result)
}

.release_normalize_gate_records <- function(raw_gates, declaration) {
  expected <- declaration$scientific_gates
  expected_ids <- vapply(expected, `[[`, character(1L), "gate_id")
  raw_ids <- vapply(raw_gates, function(gate) {
    if (is.list(gate) && .release_is_scalar_character(gate$gate_id)) {
      gate$gate_id
    } else {
      ""
    }
  }, character(1L))
  table_valid <- length(raw_gates) == length(expected) &&
    !anyDuplicated(raw_ids) && identical(raw_ids, expected_ids)
  normalized <- vector("list", length(expected))

  for (i in seq_along(expected)) {
    declared <- expected[[i]]
    matches <- which(raw_ids == declared$gate_id)
    valid <- length(matches) == 1L
    gate <- if (valid) raw_gates[[matches]] else NULL
    if (valid) {
      valid <- tryCatch({
        .release_assert_exact_fields(
          gate,
          c(
            "gate_id", "gate_version", "scope", "observed", "criterion",
            "pass", "reason"
          ),
          paste0("Adapter gate `", declared$gate_id, "`"),
          condition = "identity"
        )
        strings <- gate[c(
          "gate_id", "gate_version", "scope", "observed", "criterion", "reason"
        )]
        all(vapply(strings, .release_is_scalar_character, logical(1L))) &&
          .release_is_scalar_logical(gate$pass) &&
          identical(gate$gate_version, declared$gate_version)
      }, stablr_release_error = function(e) FALSE)
    }
    if (isTRUE(valid)) {
      if (!identical(declared$decision_status, "accepted") &&
          identical(gate$pass, TRUE)) {
        gate$pass <- FALSE
        gate$reason <- paste0(
          "gate decision status is `", declared$decision_status,
          "`; a scientific pass is unavailable"
        )
      }
      normalized[[i]] <- gate
    } else {
      table_valid <- FALSE
      normalized[[i]] <- list(
        gate_id = declared$gate_id,
        gate_version = declared$gate_version,
        scope = "required validation",
        observed = "missing or invalid",
        criterion = "adapter must report exactly one schema-valid gate",
        pass = FALSE,
        reason = "required gate record is missing, duplicated, or invalid"
      )
    }
  }
  if (!isTRUE(table_valid) && length(normalized) &&
      all(vapply(normalized, `[[`, logical(1L), "pass"))) {
    normalized[[1L]]$pass <- FALSE
    normalized[[1L]]$reason <- "gate table contains an unknown or extra record"
  }
  list(gates = normalized, complete = isTRUE(table_valid))
}

.release_collect_artifact_records <- function(result_artifacts, declaration,
                                              artifact_dir) {
  expected <- declaration$artifacts
  expected_names <- vapply(expected, `[[`, character(1L), "name")
  expected_files <- vapply(expected, `[[`, character(1L), "path")
  actual_entries <- list.files(
    artifact_dir, all.files = TRUE, no.. = TRUE, full.names = FALSE
  )
  extras <- setdiff(actual_entries, expected_files)
  if (length(extras)) {
    .release_identity_error(
      "Adapter wrote undeclared filesystem entries: ",
      paste(extras, collapse = ", "), "."
    )
  }
  result_names <- vapply(result_artifacts, function(artifact) {
    if (is.list(artifact) && .release_is_scalar_character(artifact$name)) {
      artifact$name
    } else {
      ""
    }
  }, character(1L))
  declarations_complete <- length(result_artifacts) == length(expected) &&
    !anyDuplicated(result_names) && identical(result_names, expected_names)
  records <- vector("list", length(expected))

  for (i in seq_along(expected)) {
    declared <- expected[[i]]
    relative <- file.path("artifacts", declared$path)
    path <- file.path(artifact_dir, declared$path)
    matches <- which(result_names == declared$name)
    declaration_valid <- length(matches) == 1L
    if (declaration_valid) {
      reported <- result_artifacts[[matches]]
      declaration_valid <- tryCatch({
        .release_assert_exact_fields(
          reported, c("name", "path", "schema"),
          paste0("Adapter artifact `", declared$name, "`"),
          condition = "identity"
        )
        identical(reported$name, declared$name) &&
          identical(reported$path, relative) &&
          identical(reported$schema, declared$schema)
      }, stablr_release_error = function(e) FALSE)
    }
    if (!isTRUE(declaration_valid)) declarations_complete <- FALSE

    link_target <- Sys.readlink(path)
    present <- file.exists(path)
    is_link <- !is.na(link_target) && nzchar(link_target)
    if (is_link || (present && dir.exists(path))) {
      .release_identity_error("Adapter artifact is not a regular owned file.")
    }
    schema <- if (present) {
      .release_validate_csv_artifact(path, declared$schema)
    } else {
      list(valid = FALSE, reason = "required artifact is missing")
    }
    schema_valid <- isTRUE(schema$valid) && isTRUE(declaration_valid)
    reason <- if (!isTRUE(declaration_valid)) {
      paste("adapter declaration is missing or invalid;", schema$reason)
    } else {
      schema$reason
    }
    records[[i]] <- list(
      name = declared$name,
      path = relative,
      schema = declared$schema,
      present = present,
      schema_valid = schema_valid,
      sha256 = if (present) .release_hash_file(path) else "",
      size = if (present) unname(file.info(path)$size) else 0,
      reason = reason
    )
  }
  list(
    artifacts = records,
    complete = isTRUE(declarations_complete) && all(vapply(
      records,
      function(record) isTRUE(record$present) && isTRUE(record$schema_valid),
      logical(1L)
    ))
  )
}

.release_file_manifest_entry <- function(root, relative) {
  components <- strsplit(relative, "/", fixed = TRUE)[[1L]]
  if (!.release_is_scalar_character(relative) || grepl("^/", relative) ||
      any(components %in% c("", ".", ".."))) {
    .release_identity_error("Packet manifest path is unsafe.")
  }
  path <- file.path(root, relative)
  if (!file.exists(path) || dir.exists(path) || nzchar(Sys.readlink(path))) {
    .release_identity_error("Packet manifest member is not a regular file.")
  }
  list(
    path = relative,
    sha256 = .release_hash_file(path),
    size = unname(file.info(path)$size)
  )
}

.release_validate_packet_record <- function(packet, declaration) {
  .release_assert_exact_fields(
    packet,
    c(
      "schema_version", "packet_id", "validation_id", "status",
      "promotable", "candidate_sha256", "runtime_id", "contract_sha256",
      "artifacts", "gates", "manifest_sha256"
    ),
    "Validation Run Packet",
    condition = "identity"
  )
  valid <- identical(packet$schema_version, "stablr.validation-run-packet/v1") &&
    .release_is_scalar_character(packet$packet_id) &&
    nchar(packet$packet_id) >= 16L &&
    identical(packet$validation_id, declaration$validation_id) &&
    packet$status %in% c("passing", "failed") &&
    .release_is_scalar_logical(packet$promotable) &&
    .release_is_sha256(packet$candidate_sha256) &&
    .release_is_sha256(packet$runtime_id) &&
    .release_is_sha256(packet$contract_sha256) &&
    .release_is_sha256(packet$manifest_sha256) &&
    is.list(packet$artifacts) && is.list(packet$gates)
  if (!isTRUE(valid)) {
    .release_identity_error("Validation Run Packet identity is invalid.")
  }
  invisible(packet)
}

.release_validate_packet_manifest <- function(manifest, packet_id) {
  .release_assert_exact_fields(
    manifest, c("schema_version", "packet_id", "files"),
    "Packet Manifest", condition = "identity"
  )
  if (!identical(manifest$schema_version, "stablr.packet-manifest/v1") ||
      !identical(manifest$packet_id, packet_id) ||
      !is.list(manifest$files) || !length(manifest$files)) {
    .release_identity_error("Packet Manifest identity is invalid.")
  }
  paths <- character(length(manifest$files))
  for (i in seq_along(manifest$files)) {
    member <- manifest$files[[i]]
    .release_assert_exact_fields(
      member, c("path", "sha256", "size"),
      paste0("Packet Manifest member ", i), condition = "identity"
    )
    components <- if (.release_is_scalar_character(member$path)) {
      strsplit(member$path, "/", fixed = TRUE)[[1L]]
    } else {
      ""
    }
    if (!.release_is_scalar_character(member$path) ||
        grepl("^/", member$path) ||
        any(components %in% c("", ".", "..")) ||
        !.release_is_sha256(member$sha256) ||
        !.release_is_scalar_number(member$size) || member$size < 0) {
      .release_identity_error("Packet Manifest member identity is invalid.")
    }
    paths[[i]] <- member$path
  }
  if (anyDuplicated(paths) ||
      !identical(paths, sort(paths, method = "radix"))) {
    .release_identity_error("Packet Manifest paths are ambiguous or unsorted.")
  }
  invisible(manifest)
}

.release_verify_packet <- function(context, packet_reference,
                                   require_seal = TRUE) {
  verified <- .release_verify_context(context)
  if (is.character(packet_reference) && length(packet_reference) == 1L) {
    packet_reference <- read_release_reference(packet_reference)
  }
  if (!inherits(packet_reference, "stablr_packet_ref")) {
    .release_identity_error("`packet_reference` is not a Packet Reference.")
  }
  reference <- unclass(packet_reference)
  .release_validate_packet_ref(reference)
  declaration <- .release_validation_declaration(
    verified$contract, reference$validation_id
  )
  packet_path <- .release_normalize_absolute(
    reference$packet_path, "Validation Run Packet", "file"
  )
  packet_dir <- dirname(packet_path)
  expected_dir <- file.path(
    verified$record$store, "packets", paste0("packet-", reference$packet_id)
  )
  seal_exists <- dir.exists(file.path(packet_dir, "SEALED"))
  if (!identical(packet_dir, expected_dir) ||
      !identical(packet_path, file.path(expected_dir, "packet.json")) ||
      (isTRUE(require_seal) && !seal_exists) ||
      (!isTRUE(require_seal) && seal_exists) ||
      !identical(.release_hash_file(packet_path), reference$packet_sha256)) {
    .release_identity_error("Packet Reference seal or identity failed.")
  }
  packet <- .release_read_json(packet_path, "Validation Run Packet")
  .release_validate_packet_record(packet, declaration)
  if (!identical(packet$packet_id, reference$packet_id) ||
      !identical(packet$status, reference$status) ||
      !identical(
        packet$candidate_sha256, verified$record$candidate$sha256
      ) || !identical(packet$runtime_id, verified$runtime$runtime_id) ||
      !identical(packet$contract_sha256, verified$record$contract$sha256)) {
    .release_identity_error("Validation Run Packet targets foreign evidence.")
  }

  manifest_path <- file.path(packet_dir, "packet-manifest.json")
  if (!file.exists(manifest_path) ||
      !identical(.release_hash_file(manifest_path), packet$manifest_sha256)) {
    .release_identity_error("Packet Manifest identity changed.")
  }
  manifest <- .release_read_json(manifest_path, "Packet Manifest")
  .release_validate_packet_manifest(manifest, packet$packet_id)
  manifest_ok <- vapply(manifest$files, function(member) {
    path <- file.path(packet_dir, member$path)
    file.exists(path) && !dir.exists(path) && !nzchar(Sys.readlink(path)) &&
      identical(.release_hash_file(path), member$sha256) &&
      identical(
        as.numeric(unname(file.info(path)$size)),
        as.numeric(member$size)
      )
  }, logical(1L))
  if (!all(manifest_ok)) {
    .release_identity_error("Packet Manifest file closure changed.")
  }
  manifest_paths <- vapply(manifest$files, `[[`, character(1L), "path")
  actual_files <- sort(
    list.files(
      packet_dir, all.files = TRUE, no.. = TRUE, recursive = TRUE,
      full.names = FALSE, include.dirs = FALSE
    ),
    method = "radix"
  )
  expected_files <- sort(
    c(manifest_paths, "packet-manifest.json", "packet.json"),
    method = "radix"
  )
  if (!identical(actual_files, expected_files)) {
    .release_identity_error("Validation Run Packet has an undeclared file.")
  }
  directories <- list.dirs(
    packet_dir, recursive = FALSE, full.names = FALSE
  )
  directories <- sort(directories[nzchar(directories)], method = "radix")
  expected_directories <- c("artifacts", if (isTRUE(require_seal)) "SEALED")
  if (!identical(directories, sort(expected_directories, method = "radix"))) {
    .release_identity_error("Validation Run Packet directory closure changed.")
  }
  content_paths <- file.path(packet_dir, expected_files)
  writable_bits <- bitwAnd(
    as.integer(file.info(content_paths)$mode),
    as.integer(as.octmode("0222"))
  )
  if (any(writable_bits != 0L)) {
    .release_identity_error("Sealed packet content remains writable.")
  }

  result_path <- file.path(packet_dir, "adapter-result.json")
  result <- .release_read_json(result_path, "Adapter Result")
  candidate_path <- file.path(verified$record$library$path, "stablr")
  .release_validate_adapter_result(result, declaration, candidate_path)
  recollected <- .release_collect_artifact_records(
    result$artifacts,
    declaration,
    file.path(packet_dir, "artifacts")
  )
  gates <- .release_normalize_gate_records(result$gates, declaration)
  if (!.release_json_identical(packet$artifacts, recollected$artifacts) ||
      !.release_json_identical(packet$gates, gates$gates)) {
    .release_identity_error("Packet lifecycle records changed after derivation.")
  }
  gate_pass <- gates$complete && all(vapply(
    gates$gates, function(gate) identical(gate$pass, TRUE), logical(1L)
  ))
  passing <- recollected$complete && gate_pass
  expected_status <- if (passing) "passing" else "failed"
  accepted <- all(vapply(
    declaration$scientific_gates,
    function(gate) identical(gate$decision_status, "accepted"),
    logical(1L)
  ))
  expected_promotable <- passing && isTRUE(verified$contract$promotable) &&
    accepted
  if (!identical(packet$status, expected_status) ||
      !identical(packet$promotable, expected_promotable)) {
    .release_identity_error("Packet status or promotability was not derived.")
  }
  list(
    reference = packet_reference,
    packet = packet,
    manifest = manifest,
    context = verified,
    declaration = declaration
  )
}

.release_finalize_validation_attempt <- function(verified, declaration,
                                                 attempt, packet_id, result) {
  .release_verify_context(verified$reference)
  candidate_path <- file.path(verified$record$library$path, "stablr")
  .release_validate_adapter_result(result, declaration, candidate_path)
  artifact_dir <- file.path(attempt, "artifacts")
  if (!dir.exists(artifact_dir)) {
    .release_identity_error("Validation Attempt has no artifact directory.")
  }
  artifacts <- .release_collect_artifact_records(
    result$artifacts, declaration, artifact_dir
  )
  gates <- .release_normalize_gate_records(result$gates, declaration)
  gate_pass <- gates$complete && all(vapply(
    gates$gates, function(gate) identical(gate$pass, TRUE), logical(1L)
  ))
  passing <- artifacts$complete && gate_pass
  status <- if (passing) "passing" else "failed"
  accepted <- all(vapply(
    declaration$scientific_gates,
    function(gate) identical(gate$decision_status, "accepted"),
    logical(1L)
  ))
  promotable <- passing && isTRUE(verified$contract$promotable) && accepted

  control_files <- c(
    "adapter-result.json", "context-reference.json", "process.log"
  )
  present_artifacts <- vapply(
    artifacts$artifacts,
    function(record) if (isTRUE(record$present)) record$path else "",
    character(1L)
  )
  manifest_paths <- sort(
    c(control_files, present_artifacts[nzchar(present_artifacts)]),
    method = "radix"
  )
  manifest <- list(
    schema_version = "stablr.packet-manifest/v1",
    packet_id = packet_id,
    files = lapply(
      manifest_paths,
      function(path) .release_file_manifest_entry(attempt, path)
    )
  )
  .release_validate_packet_manifest(manifest, packet_id)
  manifest_path <- file.path(attempt, "packet-manifest.json")
  manifest_sha256 <- .release_write_canonical_json(manifest, manifest_path)
  packet <- list(
    schema_version = "stablr.validation-run-packet/v1",
    packet_id = packet_id,
    validation_id = declaration$validation_id,
    status = status,
    promotable = promotable,
    candidate_sha256 = verified$record$candidate$sha256,
    runtime_id = verified$runtime$runtime_id,
    contract_sha256 = verified$record$contract$sha256,
    artifacts = artifacts$artifacts,
    gates = gates$gates,
    manifest_sha256 = manifest_sha256
  )
  .release_validate_packet_record(packet, declaration)
  packet_path <- file.path(attempt, "packet.json")
  packet_sha256 <- .release_write_canonical_json(packet, packet_path)
  .release_make_read_only(attempt, include_root = FALSE)

  packets_dir <- file.path(verified$record$store, "packets")
  dir.create(packets_dir, recursive = TRUE, showWarnings = FALSE)
  final_packet <- file.path(packets_dir, paste0("packet-", packet_id))
  if (file.exists(final_packet) || dir.exists(final_packet) ||
      !file.rename(attempt, final_packet)) {
    .release_identity_error("Atomic Validation Run Packet publication failed.")
  }
  reference <- .release_new_packet_ref(list(
    schema_version = "stablr.packet-ref/v1",
    kind = "validation_packet",
    packet_id = packet_id,
    validation_id = declaration$validation_id,
    status = status,
    packet_path = file.path(final_packet, "packet.json"),
    packet_sha256 = packet_sha256
  ))
  .release_verify_packet(
    verified$reference, reference, require_seal = FALSE
  )
  if (!dir.create(file.path(final_packet, "SEALED"), showWarnings = FALSE)) {
    .release_identity_error("Validation Run Packet seal-last publication failed.")
  }
  .release_verify_packet(verified$reference, reference)
  reference
}

run_required_validation <- function(context, validation_id) {
  verified <- .release_verify_context(context)
  declaration <- .release_validation_declaration(
    verified$contract, validation_id
  )
  packet_id <- .release_random_id()
  attempts_dir <- file.path(verified$record$store, "attempts")
  dir.create(attempts_dir, recursive = TRUE, showWarnings = FALSE)
  attempt <- file.path(attempts_dir, paste0(".attempt-", packet_id))
  if (!dir.create(attempt, showWarnings = FALSE)) {
    .release_identity_error("Could not allocate a Validation Attempt.")
  }
  artifact_dir <- file.path(attempt, "artifacts")
  temporary_dir <- file.path(attempt, "tmp")
  if (!dir.create(temporary_dir, showWarnings = FALSE)) {
    .release_identity_error("Could not allocate Validation Attempt temporary state.")
  }
  context_path <- file.path(attempt, "context-reference.json")
  write_release_reference(verified$reference, context_path)
  result_path <- file.path(attempt, "adapter-result.json")
  process_log <- file.path(attempt, "process.log")
  child <- .release_normalize_absolute(
    file.path(
      verified$record$library$path,
      "stablr", "analysis", "release_validation_child.R"
    ),
    "candidate validation child",
    "file"
  )
  environment <- c(
    .release_sanitized_environment(
      verified$contract, verified$record$library$path
    ),
    paste0("TMPDIR=", temporary_dir),
    paste0("TMP=", temporary_dir),
    paste0("TEMP=", temporary_dir)
  )
  output <- suppressWarnings(system2(
    verified$runtime$rscript$path,
    c(
      "--vanilla", shQuote(child), shQuote(context_path),
      shQuote(validation_id), shQuote(attempt), shQuote(result_path)
    ),
    stdout = TRUE,
    stderr = TRUE,
    env = environment
  ))
  status <- attr(output, "status") %||% 0L
  writeLines(output, process_log, useBytes = TRUE)
  if (status != 0L) {
    .release_incomplete_error(
      "Required Validation child exited unexpectedly with status ", status,
      ". The attempt remains unsealed.",
      attempt_path = attempt,
      process_status = status,
      process_output = output
    )
  }
  if (!file.exists(result_path) || dir.exists(result_path)) {
    .release_incomplete_error(
      "Required Validation child returned without an Adapter Result. ",
      "The attempt remains unsealed.",
      attempt_path = attempt
    )
  }
  if (dir.exists(temporary_dir)) {
    .release_discard_attempt(temporary_dir)
  }
  result <- tryCatch(
    .release_read_json(result_path, "Adapter Result"),
    stablr_release_error = function(e) {
      .release_incomplete_error(
        "Adapter Result could not be read. The attempt remains unsealed.",
        attempt_path = attempt,
        parent = e
      )
    }
  )
  verified <- .release_verify_context(verified$reference)
  reference <- .release_finalize_validation_attempt(
    verified, declaration, attempt, packet_id, result
  )
  if (identical(reference$status, "failed")) {
    .release_validation_failed_error(
      "Required Validation completed but its sealed packet failed.",
      packet_reference = reference
    )
  }
  reference
}

.release_validate_explicit_packet_set <- function(verified, packets) {
  required_ids <- vapply(
    verified$contract$required_validations,
    `[[`, character(1L), "validation_id"
  )
  if (!is.list(packets) || length(packets) != length(required_ids)) {
    .release_assembly_error(
      "Assembly requires exactly one explicit packet per Required Validation."
    )
  }
  ids <- vapply(packets, function(packet) {
    if (is.list(packet) && is.list(packet$packet) &&
        .release_is_scalar_character(packet$packet$validation_id)) {
      packet$packet$validation_id
    } else {
      ""
    }
  }, character(1L))
  if (anyDuplicated(ids) || !setequal(ids, required_ids)) {
    .release_assembly_error(
      "Explicit packet set has a duplicate, missing, extra, or foreign validation."
    )
  }
  ordered <- packets[match(required_ids, ids)]
  eligible <- vapply(ordered, function(packet) {
    identical(packet$packet$status, "passing") &&
      identical(packet$packet$promotable, TRUE) &&
      inherits(packet$reference, "stablr_packet_ref")
  }, logical(1L))
  if (!all(eligible)) {
    .release_assembly_error(
      "Every explicit packet must be sealed, passing, and promotable."
    )
  }
  if (!isTRUE(verified$contract$promotable) ||
      length(verified$contract$blockers)) {
    .release_assembly_error(
      "The Release Contract is non-promotable or has unresolved blockers.",
      blockers = verified$contract$blockers
    )
  }
  ordered
}

.release_validate_release_record <- function(release, verified) {
  .release_assert_exact_fields(
    release,
    c(
      "schema_version", "release_id", "candidate_sha256", "runtime_id",
      "contract_sha256", "packet_refs", "manifest_sha256"
    ),
    "Release Evidence Packet",
    condition = "assembly"
  )
  valid <- identical(
    release$schema_version, "stablr.release-evidence-packet/v1"
  ) && .release_is_scalar_character(release$release_id) &&
    nchar(release$release_id) >= 16L &&
    identical(
      release$candidate_sha256, verified$record$candidate$sha256
    ) && identical(release$runtime_id, verified$runtime$runtime_id) &&
    identical(
      release$contract_sha256, verified$record$contract$sha256
    ) && is.list(release$packet_refs) &&
    .release_is_sha256(release$manifest_sha256)
  if (!isTRUE(valid)) {
    .release_assembly_error("Release Evidence Packet identity is invalid.")
  }
  invisible(release)
}

.release_verify_release <- function(context, release_reference,
                                    require_seal = TRUE) {
  verified <- .release_verify_context(context)
  if (is.character(release_reference) && length(release_reference) == 1L) {
    release_reference <- read_release_reference(release_reference)
  }
  if (!inherits(release_reference, "stablr_release_ref")) {
    .release_assembly_error(
      "`release_reference` is not a Release Evidence Reference."
    )
  }
  reference <- unclass(release_reference)
  .release_validate_release_ref(reference)
  release_path <- .release_normalize_absolute(
    reference$release_path, "Release Evidence Packet", "file"
  )
  release_dir <- dirname(release_path)
  expected_dir <- file.path(
    verified$record$store,
    "releases",
    paste0("release-", reference$release_id)
  )
  seal_exists <- dir.exists(file.path(release_dir, "SEALED"))
  if (!identical(release_dir, expected_dir) ||
      !identical(release_path, file.path(expected_dir, "release.json")) ||
      (isTRUE(require_seal) && !seal_exists) ||
      (!isTRUE(require_seal) && seal_exists) ||
      !identical(.release_hash_file(release_path), reference$release_sha256)) {
    .release_assembly_error("Release Evidence Reference identity failed.")
  }
  release <- .release_read_json(release_path, "Release Evidence Packet")
  .release_validate_release_record(release, verified)
  if (!identical(release$release_id, reference$release_id)) {
    .release_assembly_error("Release Evidence Reference targets foreign state.")
  }

  manifest_path <- file.path(release_dir, "release-manifest.json")
  if (!file.exists(manifest_path) ||
      !identical(.release_hash_file(manifest_path), release$manifest_sha256)) {
    .release_assembly_error("Release Manifest identity changed.")
  }
  manifest <- .release_read_json(manifest_path, "Release Manifest")
  .release_assert_exact_fields(
    manifest, c("schema_version", "release_id", "files"),
    "Release Manifest", condition = "assembly"
  )
  if (!identical(manifest$schema_version, "stablr.release-manifest/v1") ||
      !identical(manifest$release_id, release$release_id) ||
      !is.list(manifest$files)) {
    .release_assembly_error("Release Manifest record is invalid.")
  }
  manifest_paths <- character(length(manifest$files))
  manifest_ok <- vapply(seq_along(manifest$files), function(i) {
    member <- manifest$files[[i]]
    valid_fields <- tryCatch({
      .release_assert_exact_fields(
        member, c("path", "sha256", "size"),
        paste0("Release Manifest member ", i), condition = "assembly"
      )
      TRUE
    }, stablr_release_error = function(e) FALSE)
    components <- if (.release_is_scalar_character(member$path)) {
      strsplit(member$path, "/", fixed = TRUE)[[1L]]
    } else {
      ""
    }
    if (!valid_fields || !.release_is_scalar_character(member$path) ||
        grepl("^/", member$path) ||
        any(components %in% c("", ".", "..")) ||
        !.release_is_sha256(member$sha256) ||
        !.release_is_scalar_number(member$size)) {
      return(FALSE)
    }
    manifest_paths[[i]] <<- member$path
    path <- file.path(release_dir, member$path)
    file.exists(path) && !dir.exists(path) && !nzchar(Sys.readlink(path)) &&
      identical(.release_hash_file(path), member$sha256) &&
      identical(
        as.numeric(unname(file.info(path)$size)), as.numeric(member$size)
      )
  }, logical(1L))
  if (!all(manifest_ok) || anyDuplicated(manifest_paths) ||
      !identical(manifest_paths, sort(manifest_paths, method = "radix"))) {
    .release_assembly_error("Release Manifest file closure changed.")
  }
  actual_files <- sort(
    list.files(
      release_dir, all.files = TRUE, no.. = TRUE, recursive = TRUE,
      full.names = FALSE, include.dirs = FALSE
    ),
    method = "radix"
  )
  expected_files <- sort(
    c(manifest_paths, "release-manifest.json", "release.json"),
    method = "radix"
  )
  if (!identical(actual_files, expected_files)) {
    .release_assembly_error("Release Evidence Packet has an undeclared file.")
  }
  directories <- list.dirs(
    release_dir, recursive = FALSE, full.names = FALSE
  )
  directories <- sort(directories[nzchar(directories)], method = "radix")
  expected_directories <- if (isTRUE(require_seal)) "SEALED" else character()
  if (!identical(directories, expected_directories)) {
    .release_assembly_error("Release Evidence Packet directories changed.")
  }
  content_paths <- file.path(release_dir, expected_files)
  writable_bits <- bitwAnd(
    as.integer(file.info(content_paths)$mode),
    as.integer(as.octmode("0222"))
  )
  if (any(writable_bits != 0L)) {
    .release_assembly_error("Release Evidence Packet content remains writable.")
  }

  stored_context <- read_release_reference(
    file.path(release_dir, "context-reference.json")
  )
  if (!identical(stored_context, verified$reference)) {
    .release_assembly_error("Release Context Reference changed during assembly.")
  }
  required_ids <- vapply(
    verified$contract$required_validations,
    `[[`, character(1L), "validation_id"
  )
  stored_refs <- lapply(required_ids, function(id) {
    read_release_reference(file.path(
      release_dir, paste0("packet-reference-", id, ".json")
    ))
  })
  packets <- lapply(stored_refs, function(packet) {
    tryCatch(
      .release_verify_packet(verified$reference, packet),
      stablr_release_error = function(e) {
        .release_assembly_error(
          "A stored packet failed re-verification.", parent = e
        )
      }
    )
  })
  packets <- .release_validate_explicit_packet_set(verified, packets)
  normalized_refs <- lapply(packets, function(packet) {
    unclass(packet$reference)
  })
  if (!.release_json_identical(release$packet_refs, normalized_refs)) {
    .release_assembly_error("Release packet references changed after assembly.")
  }
  list(
    reference = release_reference,
    release = release,
    manifest = manifest,
    context = verified,
    packets = packets
  )
}

assemble_release_evidence <- function(context, packet_refs) {
  verified <- .release_verify_context(context)
  if (!is.list(packet_refs) || inherits(packet_refs, "stablr_packet_ref") ||
      !length(packet_refs)) {
    .release_assembly_error("`packet_refs` must be an explicit list.")
  }
  packets <- lapply(packet_refs, function(packet) {
    tryCatch(
      .release_verify_packet(verified$reference, packet),
      stablr_release_error = function(e) {
        .release_assembly_error(
          "An explicit packet is invalid or ineligible.", parent = e
        )
      }
    )
  })
  packets <- .release_validate_explicit_packet_set(verified, packets)

  release_id <- .release_random_id()
  releases_dir <- file.path(verified$record$store, "releases")
  dir.create(releases_dir, recursive = TRUE, showWarnings = FALSE)
  staging <- file.path(releases_dir, paste0(".release-", release_id))
  final_release <- file.path(releases_dir, paste0("release-", release_id))
  if (!dir.create(staging, showWarnings = FALSE)) {
    .release_assembly_error("Could not allocate Release Evidence staging.")
  }
  write_release_reference(
    verified$reference, file.path(staging, "context-reference.json")
  )
  for (packet in packets) {
    write_release_reference(
      packet$reference,
      file.path(
        staging,
        paste0("packet-reference-", packet$packet$validation_id, ".json")
      )
    )
  }
  reference_paths <- c(
    "context-reference.json",
    paste0(
      "packet-reference-",
      vapply(packets, function(packet) packet$packet$validation_id,
             character(1L)),
      ".json"
    )
  )
  reference_paths <- sort(reference_paths, method = "radix")
  manifest <- list(
    schema_version = "stablr.release-manifest/v1",
    release_id = release_id,
    files = lapply(
      reference_paths,
      function(path) .release_file_manifest_entry(staging, path)
    )
  )
  manifest_path <- file.path(staging, "release-manifest.json")
  manifest_sha256 <- .release_write_canonical_json(manifest, manifest_path)
  release <- list(
    schema_version = "stablr.release-evidence-packet/v1",
    release_id = release_id,
    candidate_sha256 = verified$record$candidate$sha256,
    runtime_id = verified$runtime$runtime_id,
    contract_sha256 = verified$record$contract$sha256,
    packet_refs = lapply(packets, function(packet) unclass(packet$reference)),
    manifest_sha256 = manifest_sha256
  )
  .release_validate_release_record(release, verified)
  release_path <- file.path(staging, "release.json")
  release_sha256 <- .release_write_canonical_json(release, release_path)
  .release_make_read_only(staging, include_root = FALSE)
  if (file.exists(final_release) || dir.exists(final_release) ||
      !file.rename(staging, final_release)) {
    .release_assembly_error("Atomic Release Evidence publication failed.")
  }
  reference <- .release_new_release_ref(list(
    schema_version = "stablr.release-evidence-ref/v1",
    kind = "release_evidence",
    release_id = release_id,
    store = verified$record$store,
    release_path = file.path(final_release, "release.json"),
    release_sha256 = release_sha256
  ))
  .release_verify_release(
    verified$reference, reference, require_seal = FALSE
  )
  if (!dir.create(file.path(final_release, "SEALED"), showWarnings = FALSE)) {
    .release_assembly_error("Release Evidence seal-last publication failed.")
  }
  .release_verify_release(verified$reference, reference)
  reference
}
