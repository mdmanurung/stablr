#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    paste(
      "Usage: release_validation_child.R CONTEXT_REF VALIDATION_ID",
      "ATTEMPT RESULT"
    ),
    call. = FALSE
  )
}

command <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command, value = TRUE)
if (length(file_arg) != 1L) {
  stop("Could not identify the candidate adapter child script.", call. = FALSE)
}
script_path <- normalizePath(
  sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE
)
candidate_path <- dirname(dirname(script_path))
library <- dirname(candidate_path)
.libPaths(c(library, .Library))

if ("stablr" %in% loadedNamespaces()) {
  stop("The adapter child started with a preloaded stablr namespace.",
       call. = FALSE)
}
namespace <- loadNamespace("stablr")
loaded_path <- normalizePath(
  getNamespaceInfo(namespace, "path"), winslash = "/", mustWork = TRUE
)
if (!identical(loaded_path, candidate_path)) {
  stop("The adapter child did not load its installed candidate.",
       call. = FALSE)
}

module <- new.env(parent = globalenv())
sys.source(
  file.path(candidate_path, "analysis", "release_evidence.R"),
  envir = module
)
module$.release_adapter_child_main(
  context_reference_path = args[[1L]],
  validation_id = args[[2L]],
  attempt = args[[3L]],
  result_path = args[[4L]]
)
