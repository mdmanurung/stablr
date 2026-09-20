#!/usr/bin/env Rscript

# Build sealed, tamper-evident reviewer inputs. This does not create a reviewer
# decision, trust root, signature, or release authorization.

abort <- function(...) stop(paste0(...), call. = FALSE)
sha256 <- function(path) unname(tools::sha256sum(path))

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 4L || !identical(args[c(1L, 3L)], c("--out-root", "--experiments-root"))) {
    abort("Usage: build_v011_review_bundles.R --out-root ABSOLUTE_NEW_PATH --experiments-root ABSOLUTE_PATH")
  }
  list(out_root = args[[2L]], experiments_root = args[[4L]])
}

copy_bundle <- function(id, sources, out_root) {
  bundle <- file.path(out_root, id)
  if (file.exists(bundle) || dir.exists(bundle)) abort("Review bundle already exists: ", bundle)
  if (!dir.create(bundle, recursive = FALSE, showWarnings = FALSE)) abort("Cannot create bundle: ", bundle)
  destinations <- file.path(bundle, names(sources))
  if (any(!file.exists(sources)) || any(dir.exists(sources)) ||
      any(!file.copy(sources, destinations, copy.mode = FALSE, copy.date = FALSE))) {
    abort("Could not copy complete review inputs for ", id, ".")
  }
  manifest_paths <- sort(basename(destinations), method = "radix")
  lines <- paste0(
    vapply(file.path(bundle, manifest_paths), sha256, character(1L)),
    "  ", manifest_paths
  )
  manifest_path <- file.path(bundle, "manifest.sha256")
  writeLines(lines, manifest_path, useBytes = TRUE)
  manifest_sha <- sha256(manifest_path)
  writeLines(
    c(paste("bundle_id", id), paste("manifest_sha256", manifest_sha),
      "status awaiting_independent_decision"),
    file.path(bundle, "SEALED.REVIEW"), useBytes = TRUE
  )
  observed <- vapply(file.path(bundle, manifest_paths), sha256, character(1L))
  if (!identical(unname(observed), substr(lines, 1L, 64L))) abort("Bundle verification failed: ", id)
  Sys.chmod(list.files(bundle, full.names = TRUE), "0444")
  Sys.chmod(bundle, "0555")
  list(path = normalizePath(bundle, winslash = "/", mustWork = TRUE), manifest_sha256 = manifest_sha)
}

main <- function() {
  cfg <- parse_args()
  package_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  experiments_root <- normalizePath(cfg$experiments_root, winslash = "/", mustWork = TRUE)
  if (!grepl("^/", cfg$out_root) || file.exists(cfg$out_root) || dir.exists(cfg$out_root)) {
    abort("--out-root must be a new absolute path.")
  }
  parent <- normalizePath(dirname(cfg$out_root), winslash = "/", mustWork = TRUE)
  out_root <- file.path(parent, basename(cfg$out_root))
  if (!dir.create(out_root, recursive = FALSE, showWarnings = FALSE)) abort("Cannot create output root.")

  sci_root <- file.path(package_root, "inst", "analysis", "contracts")
  sci <- copy_bundle(
    "SCI-04-review-bundle-v2",
    c(
      "SCI-04-independent-review-packet-v1.md" = file.path(sci_root, "SCI-04-independent-review-packet-v1.md"),
      "SCI-04-independent-review-packet-v2.md" = file.path(sci_root, "SCI-04-independent-review-packet-v2.md"),
      "SCI-04-null-selected-fraction-v2.md" = file.path(sci_root, "SCI-04-null-selected-fraction-v2.md"),
      "sci-decision.schema.json" = file.path(sci_root, "sci-decision.schema.json"),
      "reviewer-trust-root.schema.json" = file.path(sci_root, "reviewer-trust-root.schema.json"),
      "release-contract-v1.json" = file.path(package_root, "inst", "release", "release-contract.json")
    ),
    out_root
  )

  ool_root <- file.path(experiments_root, "analysis")
  ool <- copy_bundle(
    "OOL-parity-v3-review-bundle",
    c(
      "OOL-parity-v3-protocol.md" = file.path(ool_root, "contracts", "OOL-parity-v3-protocol.md"),
      "OOL-parity-v3-review-checklist.md" = file.path(ool_root, "contracts", "OOL-parity-v3-review-checklist.md"),
      "ool-decision.schema.json" = file.path(ool_root, "contracts", "ool-decision.schema.json"),
      "reviewer-trust-root.schema.json" = file.path(ool_root, "contracts", "reviewer-trust-root.schema.json"),
      "OOL-v2-negative-evidence-identity.json" = file.path(ool_root, "contracts", "OOL-v2-negative-evidence-identity.json"),
      "ool_python_tutorial_reference_v3.json" = file.path(ool_root, "data", "ool_python_tutorial_reference_v3.json"),
      "generate_ool_parity_v3.R" = file.path(ool_root, "generate_ool_parity_v3.R"),
      "verify_ool_parity_v3.R" = file.path(ool_root, "verify_ool_parity_v3.R"),
      "evidence_contract.R" = file.path(ool_root, "evidence_contract.R")
    ),
    out_root
  )
  Sys.chmod(out_root, "0555")
  cat(sprintf("SCI %s %s\nOOL %s %s\n", sci$path, sci$manifest_sha256, ool$path, ool$manifest_sha256))
}

if (sys.nframe() == 0L) main()
