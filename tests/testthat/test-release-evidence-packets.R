.release_test_gate_records <- function(declaration, pass = TRUE) {
  lapply(declaration$scientific_gates, function(gate) list(
    gate_id = gate$gate_id,
    gate_version = gate$gate_version,
    scope = "fixture scope",
    observed = "fixture observation",
    criterion = "fixture criterion",
    pass = pass,
    reason = if (pass) "fixture passed" else "fixture failed"
  ))
}

.release_test_csv <- function(env, schema, rows = 1L) {
  columns <- env$.release_csv_schema_columns(schema)
  values <- setNames(lapply(columns, function(column) {
    if (identical(column, "pass")) rep(TRUE, rows) else rep("fixture", rows)
  }), columns)
  as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
}

test_that("RE-05 validation lifecycle conditions retain typed metadata", {
  env <- .release_module_env()
  packet <- structure(list(status = "failed"), class = "fixture_packet")

  incomplete <- tryCatch(
    env$.release_incomplete_error(
      "attempt ", 7L, " stopped", attempt_path = "/attempt"
    ),
    stablr_validation_incomplete = identity
  )
  expect_s3_class(incomplete, "stablr_validation_incomplete")
  expect_identical(conditionMessage(incomplete), "attempt 7 stopped")
  expect_identical(incomplete$attempt_path, "/attempt")

  failed <- tryCatch(
    env$.release_validation_failed_error(
      "sealed failure", packet_reference = packet
    ),
    stablr_validation_failed = identity
  )
  expect_s3_class(failed, "stablr_validation_failed")
  expect_identical(failed$packet_reference, packet)
})

test_that("RE-05 common gate normalization is complete and fail-closed", {
  env <- .release_module_env()
  contract <- env$.release_read_contract(.release_contract_source_path())
  late <- env$.release_validation_declaration(contract, "late_fusion")
  raw <- .release_test_gate_records(late)

  normalized <- env$.release_normalize_gate_records(raw, late)
  expect_true(normalized$complete)
  expect_true(all(vapply(
    normalized$gates, `[[`, logical(1L), "pass"
  )))

  missing <- env$.release_normalize_gate_records(raw[-1L], late)
  expect_false(missing$complete)
  expect_false(missing$gates[[1L]]$pass)

  wrong_version <- raw
  wrong_version[[2L]]$gate_version <- "caller-version"
  invalid <- env$.release_normalize_gate_records(wrong_version, late)
  expect_false(invalid$complete)
  expect_false(invalid$gates[[2L]]$pass)

  methodology <- env$.release_validation_declaration(contract, "methodology")
  proposed <- env$.release_normalize_gate_records(
    .release_test_gate_records(methodology), methodology
  )
  selected_fraction <- vapply(
    proposed$gates, `[[`, character(1L), "gate_id"
  ) == "null_selected_fraction"
  expect_false(proposed$gates[[which(selected_fraction)]]$pass)
})

test_that("RE-05 declared CSV artifacts are schema-bound", {
  env <- .release_module_env()
  path <- tempfile("common-gates-", fileext = ".csv")
  utils::write.csv(
    .release_test_csv(env, "common-gates/v1"), path, row.names = FALSE
  )
  expect_true(env$.release_validate_csv_artifact(
    path, "common-gates/v1"
  )$valid)

  invalid <- tempfile("bad-gates-", fileext = ".csv")
  utils::write.csv(data.frame(pass = TRUE), invalid, row.names = FALSE)
  expect_false(env$.release_validate_csv_artifact(
    invalid, "common-gates/v1"
  )$valid)
  expect_false(env$.release_validate_csv_artifact(
    path, "unknown/v1"
  )$valid)
})

test_that("RE-05 artifact closure records missing files but rejects extras", {
  env <- .release_module_env()
  contract <- env$.release_read_contract(.release_contract_source_path())
  declaration <- env$.release_validation_declaration(contract, "late_fusion")
  artifact_dir <- tempfile("adapter-artifacts-")
  dir.create(artifact_dir)
  result_artifacts <- lapply(declaration$artifacts, function(artifact) {
    data <- .release_test_csv(env, artifact$schema)
    if (grepl("warnings", artifact$schema, fixed = TRUE)) data <- data[0L, ]
    utils::write.csv(
      data, file.path(artifact_dir, artifact$path), row.names = FALSE
    )
    list(
      name = artifact$name,
      path = file.path("artifacts", artifact$path),
      schema = artifact$schema
    )
  })

  complete <- env$.release_collect_artifact_records(
    result_artifacts, declaration, artifact_dir
  )
  expect_true(complete$complete)
  expect_true(all(vapply(
    complete$artifacts, `[[`, logical(1L), "present"
  )))

  missing_path <- file.path(artifact_dir, declaration$artifacts[[1L]]$path)
  expect_true(unlink(missing_path) == 0L)
  missing <- env$.release_collect_artifact_records(
    result_artifacts, declaration, artifact_dir
  )
  expect_false(missing$complete)
  expect_false(missing$artifacts[[1L]]$present)

  writeLines("extra", file.path(artifact_dir, "foreign.txt"))
  expect_error(
    env$.release_collect_artifact_records(
      result_artifacts, declaration, artifact_dir
    ),
    class = "stablr_release_identity_error"
  )
})

test_that("RE-07 explicit packet sets assemble only when complete and eligible", {
  env <- .release_module_env()
  hash <- paste(rep("a", 64L), collapse = "")
  make_packet <- function(validation_id) {
    reference <- env$.release_new_packet_ref(list(
      schema_version = "stablr.packet-ref/v1",
      kind = "validation_packet",
      packet_id = paste0(
        if (identical(validation_id, "methodology")) "b" else "c",
        paste(rep("0", 31L), collapse = "")
      ),
      validation_id = validation_id,
      status = "passing",
      packet_path = paste0("/packets/", validation_id, "/packet.json"),
      packet_sha256 = hash
    ))
    list(
      reference = reference,
      packet = list(
        validation_id = validation_id,
        status = "passing",
        promotable = TRUE
      )
    )
  }
  verified <- list(contract = list(
    required_validations = list(
      list(validation_id = "methodology"),
      list(validation_id = "late_fusion")
    ),
    promotable = TRUE,
    blockers = list()
  ))
  methodology <- make_packet("methodology")
  late <- make_packet("late_fusion")

  ordered <- env$.release_validate_explicit_packet_set(
    verified, list(late, methodology)
  )
  expect_identical(
    vapply(ordered, function(x) x$packet$validation_id, character(1L)),
    c("methodology", "late_fusion")
  )
  expect_error(
    env$.release_validate_explicit_packet_set(
      verified, list(methodology, methodology)
    ),
    class = "stablr_release_assembly_error"
  )
  late$packet$promotable <- FALSE
  expect_error(
    env$.release_validate_explicit_packet_set(
      verified, list(methodology, late)
    ),
    class = "stablr_release_assembly_error"
  )
  late$packet$promotable <- TRUE
  verified$contract$blockers <- list("SCI-04")
  expect_error(
    env$.release_validate_explicit_packet_set(
      verified, list(methodology, late)
    ),
    class = "stablr_release_assembly_error"
  )
})
