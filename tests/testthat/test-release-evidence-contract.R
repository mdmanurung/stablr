test_that("RE-02 Release Contract is strict, versioned, and fail-closed", {
  env <- .release_module_env()
  contract <- env$.release_read_contract(.release_contract_source_path())

  expect_identical(contract$schema_version, "stablr.release-contract/v1")
  expect_identical(contract$contract_id, "stablr-0.1.1-release-v1")
  expect_false(contract$promotable)
  expect_identical(
    unlist(contract$blockers, use.names = FALSE),
    "SCI-04-null-selected-fraction-v2-pending-independent-acceptance"
  )
  expect_identical(
    vapply(
      contract$required_validations,
      `[[`,
      character(1L),
      "validation_id"
    ),
    c("methodology", "late_fusion")
  )
  methodology_gates <- contract$required_validations[[1L]]$scientific_gates
  expect_identical(
    vapply(methodology_gates, `[[`, character(1L), "gate_version"),
    c(
      "methodology-cell-completeness/v1",
      "python-metrics-parity/v1",
      "null-select-any-wilson-v1",
      "SCI-04-null-selected-fraction-v2-draft",
      "null-90pct-collapse-count-v1",
      "signal-mean-fdp-t-v1",
      "signal-tpr-t-v1"
    )
  )
  expect_match(attr(contract, "sha256"), "^[0-9a-f]{64}$")
})

test_that("RE-02 contract mutations raise typed contract errors", {
  env <- .release_module_env()
  original <- env$.release_read_contract(.release_contract_source_path())
  attr(original, "path") <- NULL
  attr(original, "sha256") <- NULL

  mutate_and_expect <- function(change) {
    contract <- .release_contract_copy(original)
    contract <- change(contract)
    expect_error(
      env$.release_validate_contract(contract),
      class = "stablr_release_contract_error"
    )
  }

  mutate_and_expect(function(x) {
    x$schema_version <- "stablr.release-contract/v2"
    x
  })
  mutate_and_expect(function(x) {
    x$hash_algorithm <- NULL
    x
  })
  mutate_and_expect(function(x) {
    x$required_validations <- x$required_validations[1L]
    x
  })
  mutate_and_expect(function(x) {
    x$required_validations[[2L]]$validation_id <- "methodology"
    x
  })
  mutate_and_expect(function(x) {
    x$required_validations[[1L]]$adapter <- "caller_adapter"
    x
  })
  mutate_and_expect(function(x) {
    x$required_validations[[1L]]$settings$replicates <- 99L
    x
  })
  mutate_and_expect(function(x) {
    x$required_validations[[1L]]$artifacts[[1L]]$path <- "../escape.csv"
    x
  })
  mutate_and_expect(function(x) {
    gate <- x$required_validations[[1L]]$scientific_gates[[4L]]
    gate$decision_status <- "accepted"
    x$required_validations[[1L]]$scientific_gates[[4L]] <- gate
    x
  })
  mutate_and_expect(function(x) {
    x$required_validations[[1L]]$scientific_gates[[3L]]$gate_version <-
      "caller-version"
    x
  })
  mutate_and_expect(function(x) {
    x$runtime$dependency_roots <- x$runtime$dependency_roots[-1L]
    x
  })
  mutate_and_expect(function(x) {
    x$caller_settings <- list(replicates = 1L)
    x
  })
})

test_that("RE-02 schemas are packaged JSON records", {
  schema_dir <- testthat::test_path("..", "..", "inst", "release", "schemas")
  if (!dir.exists(schema_dir)) {
    schema_dir <- system.file("release", "schemas", package = "stablr")
  }
  expected <- c(
    "packet-reference.schema.json",
    "release-context.schema.json",
    "release-contract.schema.json",
    "release-evidence-packet.schema.json",
    "release-reference.schema.json",
    "runtime-artifact.schema.json",
    "validation-run-packet.schema.json"
  )
  expect_setequal(list.files(schema_dir, pattern = "[.]json$"), expected)
  for (name in expected) {
    schema <- jsonlite::fromJSON(
      file.path(schema_dir, name), simplifyVector = FALSE
    )
    expect_identical(schema[["$schema"]],
                     "https://json-schema.org/draft/2020-12/schema")
    expect_false(is.null(schema[["$id"]]))
  }
})

test_that("RE-02 context and packet references round-trip canonically", {
  env <- .release_module_env()
  root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
  context <- env$.release_new_context_ref(list(
    schema_version = "stablr.release-context-ref/v1",
    kind = "release_context",
    context_id = paste(rep("a", 32L), collapse = ""),
    store = root,
    context_path = file.path(root, "context.json"),
    context_sha256 = paste(rep("b", 64L), collapse = "")
  ))
  packet <- env$.release_new_packet_ref(list(
    schema_version = "stablr.packet-ref/v1",
    kind = "validation_packet",
    packet_id = paste(rep("c", 32L), collapse = ""),
    validation_id = "methodology",
    status = "failed",
    packet_path = file.path(root, "packet.json"),
    packet_sha256 = paste(rep("d", 64L), collapse = "")
  ))
  release <- env$.release_new_release_ref(list(
    schema_version = "stablr.release-evidence-ref/v1",
    kind = "release_evidence",
    release_id = paste(rep("e", 32L), collapse = ""),
    store = root,
    release_path = file.path(root, "release.json"),
    release_sha256 = paste(rep("f", 64L), collapse = "")
  ))

  context_path <- tempfile("context-ref-", fileext = ".json")
  packet_path <- tempfile("packet-ref-", fileext = ".json")
  release_path <- tempfile("release-ref-", fileext = ".json")
  env$write_release_reference(context, context_path)
  env$write_release_reference(packet, packet_path)
  env$write_release_reference(release, release_path)

  expect_identical(env$read_release_reference(context_path), context)
  expect_identical(env$read_release_reference(packet_path), packet)
  expect_identical(env$read_release_reference(release_path), release)
  expect_false(grepl("[\r\n]", readChar(
    context_path, file.info(context_path)$size, useBytes = TRUE
  )))
})

test_that("RE-02 malformed references fail with identity conditions", {
  env <- .release_module_env()
  root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
  record <- list(
    schema_version = "stablr.release-context-ref/v1",
    kind = "release_context",
    context_id = paste(rep("a", 32L), collapse = ""),
    store = root,
    context_path = file.path(root, "context.json"),
    context_sha256 = paste(rep("b", 64L), collapse = ""),
    caller_override = TRUE
  )

  expect_error(
    env$.release_new_context_ref(record),
    class = "stablr_release_identity_error"
  )
  expect_error(
    env$.release_new_packet_ref(record),
    class = "stablr_release_identity_error"
  )
})

test_that("RE-02 jsonlite is declared only as release tooling Suggests", {
  description <- read.dcf(testthat::test_path("..", "..", "DESCRIPTION"))
  suggests <- trimws(strsplit(description[1L, "Suggests"], ",")[[1L]])
  imports <- trimws(strsplit(description[1L, "Imports"], ",")[[1L]])
  expect_true("jsonlite" %in% suggests)
  expect_false("jsonlite" %in% imports)
})
