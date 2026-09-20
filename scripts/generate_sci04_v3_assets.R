#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]))
root <- normalizePath(file.path(dirname(script_path), ".."))
contract_dir <- file.path(root, "inst", "analysis", "contracts")

source(file.path(root, "inst", "analysis", "run_methodology_validation.R"),
       local = TRUE)

sha256_file <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE)
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

write_json <- function(object, path) {
  json <- jsonlite::toJSON(
    object,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    digits = NA
  )
  writeLines(json, path, useBytes = TRUE)
}

matrix_path <- file.path(contract_dir, "SCI-04-applicability-v1.csv")
reference_path <- file.path(contract_dir, "SCI-04-reference-cases-v1.json")
specification_path <- file.path(
  contract_dir,
  "SCI-04-statistical-specification-v3.json"
)

applicability <- .methodology_applicability_matrix()
applicability <- applicability[c(
  "family", "scenario", "profile", "artificial_type", "applicable",
  "reason_code"
)]
applicability <- applicability[order(
  match(applicability$family, .advertised_families),
  match(applicability$scenario, .default_scenarios()$scenario),
  match(applicability$artificial_type, .advertised_artificial_types)
), , drop = FALSE]
utils::write.csv(applicability, matrix_path, row.names = FALSE, na = "")

alpha <- 0.05
replicates <- 600L
threshold <- 0.10
radius <- sqrt(log(1 / alpha) / (2 * replicates))
boundary_mean <- threshold - radius

reference_cases <- list(
  schema_version = "stablr.sci04-reference-cases/v1",
  formula = "min(1, mean(x) + sqrt(log(1 / alpha) / (2 * n)))",
  alpha = alpha,
  expected_replicates = replicates,
  release_threshold = threshold,
  cases = list(
    list(
      id = "all_zero",
      input = list(kind = "constant", count = replicates, value = 0),
      expected_upper_bound = radius,
      expected_pass = TRUE
    ),
    list(
      id = "boundary_equality",
      input = list(
        kind = "constant",
        count = replicates,
        value = boundary_mean
      ),
      expected_upper_bound = threshold,
      expected_pass = TRUE,
      comparison_tolerance = 1e-15
    ),
    list(
      id = "above_boundary",
      input = list(
        kind = "constant",
        count = replicates,
        value = boundary_mean + 1e-10
      ),
      expected_upper_bound = threshold + 1e-10,
      expected_pass = FALSE,
      comparison_tolerance = 1e-15
    ),
    list(
      id = "upper_clamp",
      input = list(kind = "constant", count = replicates, value = 1),
      expected_upper_bound = 1,
      expected_pass = FALSE
    )
  ),
  fail_closed_cases = list(
    list(id = "missing_value", mutation = "one x value is NA"),
    list(id = "non_finite_value", mutation = "one x value is Inf"),
    list(id = "below_range", mutation = "one x value is below zero"),
    list(id = "above_range", mutation = "one x value exceeds one"),
    list(id = "incomplete_count", mutation = "599 values supplied"),
    list(id = "excess_count", mutation = "601 values supplied"),
    list(id = "duplicate_replicate", mutation = "replicate id duplicated"),
    list(id = "non_integer_count", mutation = "n_selected is non-integer"),
    list(id = "selected_exceeds_p", mutation = "n_selected exceeds p")
  )
)
write_json(reference_cases, reference_path)

matrix_hash <- sha256_file(matrix_path)
reference_hash <- sha256_file(reference_path)
specification <- list(
  schema_version = "stablr.sci04-statistical-specification/v3",
  specification_id = "SCI-04-null-selected-fraction-hoeffding-iut-v3",
  status = "implementation_authorized_pending_independent_signature",
  implementation_authority = "mdm123",
  independent_acceptance = FALSE,
  estimand = "mean_replicate_selected_feature_fraction",
  experimental_unit = "independently_seeded_monte_carlo_replicate",
  bound = list(method = "fixed_n_hoeffding", side = "upper", alpha = alpha),
  decision_structure = "intersection_union_all_applicable_cells",
  release_multiplicity_adjustment = "none",
  simultaneous_reporting_policy = "not_claimed",
  replicates_per_applicable_cell = replicates,
  release_threshold = threshold,
  boundary_equality_passes = TRUE,
  declared_cells = 96L,
  applicable_cells = 88L,
  not_applicable_cells = 8L,
  applicable_null_cells = 44L,
  malformed_data_policy = "fail_closed",
  applicability_matrix = list(
    path = "SCI-04-applicability-v1.csv",
    sha256 = matrix_hash
  ),
  reference_cases = list(
    path = "SCI-04-reference-cases-v1.json",
    sha256 = reference_hash
  )
)
write_json(specification, specification_path)

message("Wrote SCI-04 v3 assets:")
message("  ", matrix_path)
message("  ", reference_path)
message("  ", specification_path)
