#' stablr: Sparse and Reliable Biomarker Discovery in R
#'
#' `stablr` is an R-native implementation of STABL for sparse, stable biomarker
#' selection in high-dimensional clinical and omic datasets. The package
#' compiles C++ for its cooperative-learning backend and
#' provides a core bootstrap stability-selection engine, FDP+ threshold
#' calibration with artificial features, glmnet-family learner adapters, and
#' multi-omic workflows for per-omic, early-fusion, late-fusion, and
#' cooperative-fusion analyses.
#'
#' @section Main workflows:
#' \itemize{
#'   \item [stabl_fit()] fits the core single-matrix STABL selector.
#'   \item [stabl_multiomic_train_validate()] runs train/validation multi-omic
#'     workflows with optional early, late, and cooperative fusion.
#'   \item [stabl_multiomic_cv()] runs outer cross-validation for named
#'     multi-omic inputs, preserving optional group structure.
#' }
#'
#' @section Learners and outcomes:
#' The core selector supports lasso, elastic net, adaptive lasso, and optional
#' sparse group lasso backends.  Supported `glmnet` families include Gaussian,
#' binomial, multinomial, and Cox where the selected backend supports them.
#'
#' @section Outputs:
#' Fitted objects expose stable S3 accessors for support masks, selected feature
#' names, stability scores, importances, cooperative-fusion features, and
#' cooperative diagnostics.  Plotting and export helpers provide stability
#' paths, FDP+ curves, ROC/PRC plots, selected-feature visualizations, and CSV
#' output.
#'
#' @seealso [stabl_fit()], [stabl_multiomic_train_validate()],
#'   [get_support()], [get_importances()], [plot_stabl_path()]
"_PACKAGE"

# ".data" is the rlang data-mask pronoun used inside ggplot2::aes(); rlang
# stays in Suggests (it arrives with ggplot2), so the pronoun is registered
# here instead of via importFrom(rlang, .data), which would make rlang a
# hard namespace dependency (R CMD check ERROR).
utils::globalVariables(c(
  ".data", "artificial", "FDR", "count", "feature", "fdp_plus", "fpr",
  "frequency", "kind", "lambda_idx", "multiview.cox.path", "outcome",
  "panel_size", "precision", "proportion", "real", "realized", "recall",
  "rep", "reported", "score", "threshold", "tpr", "value"
))

# Private helper so tests can mock Suggests-package availability.
.has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# Unified package-guard helper: stop with an actionable install message if
# `pkg` is not installed. `context` is appended after "is required" when given.
.require_pkg <- function(pkg, context = NULL) {
  if (.has_pkg(pkg)) return(invisible(NULL))
  what <- paste0("Package '", pkg, "' is required",
                 if (!is.null(context)) paste0(" ", context) else "", ".")
  stop(what, "\nInstall with: install.packages(\"", pkg, "\")", call. = FALSE)
}
