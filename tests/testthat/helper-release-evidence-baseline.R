.re01_split_csv_line <- function(line) {
  chars <- strsplit(line, "", fixed = TRUE)[[1L]]
  fields <- character()
  token <- character()
  quoted <- FALSE
  i <- 1L

  while (i <= length(chars)) {
    char <- chars[[i]]
    if (identical(char, '"')) {
      token <- c(token, char)
      if (quoted && i < length(chars) && identical(chars[[i + 1L]], '"')) {
        i <- i + 1L
        token <- c(token, chars[[i]])
      } else {
        quoted <- !quoted
      }
    } else if (identical(char, ",") && !quoted) {
      fields <- c(fields, paste0(token, collapse = ""))
      token <- character()
    } else {
      token <- c(token, char)
    }
    i <- i + 1L
  }

  if (quoted) {
    stop("Multiline or unterminated quoted CSV field is not supported.", call. = FALSE)
  }
  c(fields, paste0(token, collapse = ""))
}

.re01_csv_field_name <- function(field) {
  if (startsWith(field, '"') && endsWith(field, '"')) {
    field <- substr(field, 2L, nchar(field) - 1L)
    field <- gsub('""', '"', field, fixed = TRUE)
  }
  field
}

.re01_drop_csv_column <- function(input, output, column) {
  lines <- readLines(input, warn = FALSE, encoding = "UTF-8")
  if (!length(lines)) {
    stop("Cannot normalize an empty CSV: ", input, call. = FALSE)
  }

  fields <- lapply(lines, .re01_split_csv_line)
  widths <- lengths(fields)
  if (any(widths != widths[[1L]])) {
    stop("CSV rows have inconsistent field counts: ", input, call. = FALSE)
  }

  header <- vapply(fields[[1L]], .re01_csv_field_name, character(1L))
  drop <- which(header == column)
  if (length(drop) != 1L) {
    stop(
      "Expected exactly one `", column, "` column in ", input,
      "; found ", length(drop), ".",
      call. = FALSE
    )
  }

  normalized <- vapply(
    fields,
    function(row) paste(row[-drop], collapse = ","),
    character(1L)
  )
  writeLines(normalized, output, useBytes = TRUE)
  invisible(output)
}

.re01_copy_exact <- function(input, output) {
  copied <- file.copy(input, output, overwrite = FALSE, copy.mode = TRUE)
  if (!isTRUE(copied)) {
    stop("Could not copy baseline artifact: ", input, call. = FALSE)
  }
  invisible(output)
}

.re01_normalize_colon_key <- function(lines, key) {
  prefix <- paste0(key, ": ")
  index <- which(startsWith(lines, prefix))
  if (length(index) != 1L) {
    stop("Expected exactly one manifest key `", key, "`.", call. = FALSE)
  }
  lines[[index]] <- paste0(prefix, "<normalized>")
  lines
}

.re01_normalize_space_key <- function(lines, key) {
  prefix <- paste0(key, " ")
  index <- which(startsWith(lines, prefix))
  if (length(index) != 1L) {
    stop("Expected exactly one manifest key `", key, "`.", call. = FALSE)
  }
  lines[[index]] <- paste0(prefix, "<normalized>")
  lines
}

.re01_normalize_methodology_manifest <- function(input, output) {
  lines <- readLines(input, warn = FALSE, encoding = "UTF-8")
  scalar_keys <- c(
    "created_at", "package_mode", "source_git_commit", "source_git_tree",
    "source_tracked_dirty", "end_git_commit", "end_git_tree",
    "end_tracked_dirty", "source_git_stable"
  )
  for (key in scalar_keys) {
    lines <- .re01_normalize_colon_key(lines, key)
  }

  artifact_names <- c("replicates", "summary", "warnings", "parity", "gates")
  artifact_header <- which(lines == "artifacts:")
  hash_header <- which(lines == "artifact_sha256:")
  if (length(artifact_header) != 1L || length(hash_header) != 1L) {
    stop("Methodology manifest sections are missing or duplicated.", call. = FALSE)
  }

  artifact_rows <- artifact_header + seq_along(artifact_names)
  artifact_keys <- sub(":.*$", "", lines[artifact_rows])
  if (!identical(artifact_keys, artifact_names)) {
    stop("Methodology artifact path order changed.", call. = FALSE)
  }
  lines[artifact_rows] <- paste0(artifact_names, ": <normalized>")

  hash_rows <- hash_header + seq_along(artifact_names)
  hash_names <- sub(" .*", "", lines[hash_rows])
  expected_hash_names <- c(
    "methodology_validation_replicates.csv",
    "methodology_validation_summary.csv",
    "methodology_validation_warnings.csv",
    "python_metrics_parity.csv",
    "methodology_validation_gates.csv"
  )
  if (!identical(hash_names, expected_hash_names)) {
    stop("Methodology artifact hash order changed.", call. = FALSE)
  }
  lines[hash_rows] <- paste(expected_hash_names, "<normalized>")

  writeLines(lines, output, useBytes = TRUE)
  invisible(output)
}

.re01_normalize_late_fusion_manifest <- function(input, output) {
  lines <- readLines(input, warn = FALSE, encoding = "UTF-8")
  scalar_keys <- c(
    "package_mode", "source_git_commit", "source_git_tree",
    "source_tracked_dirty", "end_git_commit", "end_git_tree",
    "end_tracked_dirty", "source_git_stable"
  )
  for (key in scalar_keys) {
    lines <- .re01_normalize_space_key(lines, key)
  }

  hash_header <- which(lines == "artifact_sha256:")
  if (length(hash_header) != 1L) {
    stop("Late-fusion hash section is missing or duplicated.", call. = FALSE)
  }
  expected_hash_names <- c(
    "late_fusion_replicates.csv",
    "late_fusion_summary.csv",
    "late_fusion_warnings.csv"
  )
  hash_rows <- hash_header + seq_along(expected_hash_names)
  hash_names <- sub(" .*", "", lines[hash_rows])
  if (!identical(hash_names, expected_hash_names)) {
    stop("Late-fusion artifact hash order changed.", call. = FALSE)
  }
  lines[hash_rows] <- paste(expected_hash_names, "<normalized>")

  writeLines(lines, output, useBytes = TRUE)
  invisible(output)
}

.normalize_re01_methodology <- function(input_dir, output_dir) {
  expected <- c(
    "methodology_validation_gates.csv",
    "methodology_validation_manifest.txt",
    "methodology_validation_replicates.csv",
    "methodology_validation_summary.csv",
    "methodology_validation_warnings.csv",
    "python_metrics_parity.csv"
  )
  actual <- sort(list.files(input_dir, all.files = FALSE, no.. = TRUE))
  if (!identical(actual, sort(expected))) {
    stop("Methodology baseline artifact set changed.", call. = FALSE)
  }
  if (dir.exists(output_dir) || file.exists(output_dir)) {
    stop("Refusing to replace normalization output: ", output_dir, call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  .re01_drop_csv_column(
    file.path(input_dir, "methodology_validation_replicates.csv"),
    file.path(output_dir, "methodology_validation_replicates.csv"),
    "elapsed_sec"
  )
  .re01_drop_csv_column(
    file.path(input_dir, "methodology_validation_summary.csv"),
    file.path(output_dir, "methodology_validation_summary.csv"),
    "mean_elapsed_sec"
  )
  for (name in c(
    "methodology_validation_warnings.csv", "python_metrics_parity.csv",
    "methodology_validation_gates.csv"
  )) {
    .re01_copy_exact(file.path(input_dir, name), file.path(output_dir, name))
  }
  .re01_normalize_methodology_manifest(
    file.path(input_dir, "methodology_validation_manifest.txt"),
    file.path(output_dir, "methodology_validation_manifest.txt")
  )
  invisible(output_dir)
}

.normalize_re01_late_fusion <- function(input_dir, output_dir) {
  expected <- c(
    "late_fusion_manifest.txt",
    "late_fusion_replicates.csv",
    "late_fusion_summary.csv",
    "late_fusion_warnings.csv"
  )
  actual <- sort(list.files(input_dir, all.files = FALSE, no.. = TRUE))
  if (!identical(actual, sort(expected))) {
    stop("Late-fusion baseline artifact set changed.", call. = FALSE)
  }
  if (dir.exists(output_dir) || file.exists(output_dir)) {
    stop("Refusing to replace normalization output: ", output_dir, call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (name in c(
    "late_fusion_replicates.csv", "late_fusion_summary.csv",
    "late_fusion_warnings.csv"
  )) {
    .re01_copy_exact(file.path(input_dir, name), file.path(output_dir, name))
  }
  .re01_normalize_late_fusion_manifest(
    file.path(input_dir, "late_fusion_manifest.txt"),
    file.path(output_dir, "late_fusion_manifest.txt")
  )
  invisible(output_dir)
}
