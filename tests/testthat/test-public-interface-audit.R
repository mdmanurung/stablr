.interface_source_root <- function() {
  candidates <- unique(c(
    normalizePath(testthat::test_path("..", ".."), mustWork = FALSE),
    normalizePath(getwd(), mustWork = FALSE),
    normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
  ))
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "NAMESPACE")) &&
        dir.exists(file.path(candidate, "inst", "release"))) {
      return(candidate)
    }
  }
  testthat::skip("package source interface audit is not available")
}

.read_interface_audit <- function() {
  utils::read.csv(
    file.path(
      .interface_source_root(), "inst", "release",
      "public-interface-audit.csv"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = character()
  )
}

test_that("public-interface audit exactly covers NAMESPACE", {
  root <- .interface_source_root()
  audit <- .read_interface_audit()
  expect_identical(
    names(audit),
    c(
      "interface", "kind", "generic", "class", "area", "implementation",
      "documentation", "example_status", "valid_evidence",
      "invalid_evidence", "coverage"
    )
  )
  expect_identical(anyDuplicated(audit$interface), 0L)

  namespace <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
  export_lines <- namespace[startsWith(namespace, "export(")]
  exports <- substring(export_lines, 8L, nchar(export_lines) - 1L)
  s3_lines <- namespace[startsWith(namespace, "S3method(")]
  s3_parts <- strsplit(
    substring(s3_lines, 10L, nchar(s3_lines) - 1L), ",", fixed = TRUE
  )
  s3 <- vapply(
    s3_parts, function(x) paste0(x[[1L]], ".", x[[2L]]), character(1L)
  )

  expect_length(exports, 42L)
  expect_length(s3, 19L)
  expect_setequal(audit$interface[audit$kind == "export"], exports)
  expect_setequal(audit$interface[audit$kind == "s3"], s3)
})

test_that("audited interfaces have code documentation examples and tests", {
  root <- .interface_source_root()
  audit <- .read_interface_audit()
  expect_true(all(audit$kind %in% c("export", "s3")))
  expect_true(all(audit$example_status %in% c(
    "running", "runnable_donttest", "running_workflow",
    "workflow_integration"
  )))
  expect_true(all(audit$coverage %in% c(
    "core_contract", "utility_contract", "accessor_contract",
    "vendored_parity"
  )))
  expect_true(all(file.exists(file.path(root, audit$implementation))))
  expect_true(all(file.exists(file.path(root, audit$documentation))))

  evidence_files <- unique(c(
    audit$valid_evidence,
    audit$invalid_evidence[!startsWith(audit$invalid_evidence, "not_applicable")]
  ))
  expect_true(all(file.exists(file.path(root, "tests", "testthat", evidence_files))))

  exports <- audit[audit$kind == "export", , drop = FALSE]
  for (i in seq_len(nrow(exports))) {
    rd <- readLines(file.path(root, exports$documentation[[i]]), warn = FALSE)
    text <- paste(rd, collapse = "\n")
    expect_true(
      any(rd == paste0("\\alias{", exports$interface[[i]], "}")),
      info = exports$interface[[i]]
    )
    expect_match(text, "\\\\examples\\{", info = exports$interface[[i]])
  }
})

test_that("all audited S3 methods are registered and callable", {
  audit <- .read_interface_audit()
  methods <- audit[audit$kind == "s3", , drop = FALSE]
  registered <- vapply(seq_len(nrow(methods)), function(i) {
    is.function(getS3method(
      methods$generic[[i]], methods$class[[i]], optional = TRUE
    ))
  }, logical(1L))
  expect_true(all(registered), info = paste(methods$interface[!registered], collapse = ", "))
})

test_that("editorial rubric contains every acceptance gate", {
  rubric <- paste(
    readLines(file.path(
      .interface_source_root(), "inst", "release", "editorial-acceptance.md"
    ), warn = FALSE),
    collapse = "\n"
  )
  for (term in c(
    "concrete", "consistent", "authoritative", "stock promotion",
    "Examples execute", "links resolve", "rendered page"
  )) {
    expect_match(rubric, term, fixed = TRUE)
  }
})
