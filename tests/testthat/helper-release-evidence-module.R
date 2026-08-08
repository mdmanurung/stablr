.release_evidence_source_path <- function() {
  source_path <- testthat::test_path(
    "..", "..", "inst", "analysis", "release_evidence.R"
  )
  if (file.exists(source_path)) return(source_path)
  installed <- system.file("analysis", "release_evidence.R", package = "stablr")
  if (nzchar(installed)) return(installed)
  stop("Could not locate the release-evidence module.", call. = FALSE)
}

.release_contract_source_path <- function() {
  source_path <- testthat::test_path(
    "..", "..", "inst", "release", "release-contract.json"
  )
  if (file.exists(source_path)) return(source_path)
  installed <- system.file(
    "release", "release-contract.json", package = "stablr"
  )
  if (nzchar(installed)) return(installed)
  stop("Could not locate the Release Contract.", call. = FALSE)
}

.release_module_env <- function() {
  env <- new.env(parent = globalenv())
  sys.source(.release_evidence_source_path(), envir = env)
  env
}

.release_contract_copy <- function(contract) {
  unserialize(serialize(contract, NULL))
}
