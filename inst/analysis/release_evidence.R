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
      "--no-multiarch", "--no-test-load", shQuote(candidate_tarball)
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
