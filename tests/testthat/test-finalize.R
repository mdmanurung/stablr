# Tests for stabl_finalize() / predict.stabl_predictor() /
# validate_newdata_schema() (1.0.0 spec A6; Phase 13 PR-F).
#
# Coverage: saveRDS/readRDS roundtrip, column-order invariance, schema
# violations (missing / renamed / duplicate / extra / non-numeric),
# training-data-free serialization, engine ladder (cvglmnet, one-feature
# glm, intercept-only and prevalence fallbacks), stored-preprocessing
# application without refit, development-summary provenance from
# stabl_analysis inputs, event mapping, and print smoke.

# ---- shared fixtures --------------------------------------------------------

.make_fin_data <- function(n = 100L, p = 12L, seed = 21L) {
  set.seed(seed)
  ids <- paste0("s", seq_len(n))
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(ids, paste0("f", seq_len(p))))
  yg <- setNames(2 * x[, 1L] - 1.5 * x[, 2L] + x[, 3L] + rnorm(n, sd = 1.5),
                 ids)
  yb <- setNames(factor(ifelse(yg > stats::median(yg), "case", "ctrl")), ids)
  list(x = x, yg = yg, yb = yb, ids = ids)
}

grid3 <- data.frame(lambda = c(0.3, 0.2, 0.1))

.make_fin <- function(seed = 21L) {
  d <- .make_fin_data(seed = seed)
  fit_g <- stabl_fit(d$x, d$yg, lambda_grid = grid3, family = "gaussian",
                     n_bootstraps = 30L, random_state = seed)
  fit_b <- stabl_fit(d$x, d$yb, lambda_grid = grid3, family = "binomial",
                     n_bootstraps = 30L, random_state = seed)
  list(d = d, fit_g = fit_g, fit_b = fit_b)
}

# ---- core finalize / predict ------------------------------------------------

test_that("gaussian finalize predicts and roundtrips through saveRDS", {
  f <- .make_fin()
  fin <- stabl_finalize(f$fit_g, f$d$x, f$d$yg,
                        preprocessing = c("center", "scale"),
                        random_state = 21L)
  expect_s3_class(fin, "stabl_predictor")
  expect_identical(fin$family, "gaussian")
  expect_true(length(fin$selected_features) >= 1L)
  p1 <- predict(fin, f$d$x[1:8, ])
  expect_type(p1, "double")
  expect_length(p1, 8L)

  tmp <- tempfile(fileext = ".rds")
  saveRDS(fin, tmp)
  fin2 <- readRDS(tmp)
  expect_identical(predict(fin, f$d$x[1:8, ]), predict(fin2, f$d$x[1:8, ]))
})

test_that("binomial finalize gives consistent class and probability", {
  f <- .make_fin()
  fin <- stabl_finalize(f$fit_b, f$d$x, f$d$yb, random_state = 21L)
  expect_identical(fin$family, "binomial")
  xb <- f$d$x[, fin$required_columns, drop = FALSE]
  cl <- predict(fin, xb[1:10, ])
  pr <- predict(fin, xb[1:10, ], type = "probability")
  expect_s3_class(cl, "factor")
  expect_identical(levels(cl), fin$levels)
  expect_true(all(pr >= 0 & pr <= 1))
  ev <- fin$event_class
  other <- setdiff(fin$levels, ev)
  expect_identical(as.character(cl), ifelse(pr >= 0.5, ev, other))

  tmp <- tempfile(fileext = ".rds")
  saveRDS(fin, tmp)
  fin2 <- readRDS(tmp)
  expect_identical(predict(fin, xb[1:10, ]), predict(fin2, xb[1:10, ]))
  expect_identical(predict(fin, xb[1:10, ], type = "probability"),
                   predict(fin2, xb[1:10, ], type = "probability"))
})

test_that("non-default event_class is respected in probabilities", {
  f <- .make_fin()
  ev1 <- f$fit_b$family  # sanity
  lv <- levels(factor(f$d$yb))
  fin_a <- stabl_finalize(f$fit_b, f$d$x, f$d$yb, random_state = 21L)
  fin_b <- stabl_finalize(f$fit_b, f$d$x, f$d$yb, event_class = lv[1L],
                          random_state = 21L)
  xb <- f$d$x[, fin_a$required_columns, drop = FALSE]
  pa <- predict(fin_a, xb[1:10, ], type = "probability")
  pb <- predict(fin_b, xb[1:10, ], type = "probability")
  expect_equal(pa + pb, rep(1, 10L), tolerance = 1e-12)
  # Hard classes are identical: the event mapping does not move the 0.5
  # threshold on P(second level).
  expect_identical(predict(fin_a, xb[1:10, ]), predict(fin_b, xb[1:10, ]))
})

test_that("column order is irrelevant and extra columns need explicit drop", {
  f <- .make_fin()
  fin <- stabl_finalize(f$fit_g, f$d$x, f$d$yg,
                        preprocessing = c("center", "scale"),
                        random_state = 21L)
  perm <- sample(ncol(f$d$x))
  expect_identical(predict(fin, f$d$x[1:10, ]),
                   predict(fin, f$d$x[1:10, perm]))
  x_extra <- cbind(f$d$x[1:5, ], brand_new = rnorm(5L))
  expect_error(predict(fin, x_extra), "unexpected column", fixed = TRUE)
  expect_identical(predict(fin, x_extra, extra_columns = "drop"),
                   predict(fin, f$d$x[1:5, ]))
})

test_that("schema violations are fail-closed with clear messages", {
  f <- .make_fin()
  fin <- stabl_finalize(f$fit_g, f$d$x, f$d$yg, random_state = 21L)
  req <- fin$required_columns

  x_miss <- f$d$x[1:3, setdiff(req, req[1L]), drop = FALSE]
  expect_error(validate_newdata_schema(fin, x_miss),
               "missing", ignore.case = TRUE)

  x_dup <- cbind(f$d$x[1:3, req], f$d$x[1:3, req[1L], drop = FALSE])
  expect_error(validate_newdata_schema(fin, x_dup), "unique", fixed = TRUE)

  x_ren <- f$d$x[1:3, req]
  colnames(x_ren)[1L] <- "renamed_feature"
  expect_error(validate_newdata_schema(fin, x_ren),
               "missing", ignore.case = TRUE)

  df_bad <- as.data.frame(f$d$x[1:3, req])
  df_bad[[1L]] <- letters[1:3]
  expect_error(validate_newdata_schema(fin, df_bad), "numeric", fixed = TRUE)

  x_noname <- f$d$x[1:3, req]
  colnames(x_noname) <- NULL
  expect_error(validate_newdata_schema(fin, x_noname), "column", fixed = TRUE)

  expect_error(validate_newdata_schema(list(), f$d$x[1:3, req]),
               "stabl_predictor", fixed = TRUE)

  # Valid input returns the required columns in training order.
  shuffled <- f$d$x[1:4, rev(req), drop = FALSE]
  out <- validate_newdata_schema(fin, shuffled)
  expect_identical(colnames(out), req)
  expect_equal(nrow(out), 4L)
})

test_that("no training observations are serialized", {
  f <- .make_fin()
  x2 <- f$d$x
  x2[1L, 1L] <- 987654.321  # distinctive sentinel
  fit2 <- stabl_fit(x2, f$d$yg, lambda_grid = grid3, family = "gaussian",
                    n_bootstraps = 15L, random_state = 1L)
  fin <- stabl_finalize(fit2, x2, f$d$yg,
                        preprocessing = c("center", "scale"),
                        random_state = 1L)
  tmp <- tempfile(fileext = ".rds")
  saveRDS(fin, tmp, compress = FALSE)
  raw <- readBin(tmp, "raw", file.size(tmp))
  needle <- charToRaw(format(987654.321, digits = 12))
  hits <- vapply(
    seq_len(length(raw) - length(needle)),
    function(i) identical(raw[i:(i + length(needle) - 1L)], needle),
    logical(1L)
  )
  expect_false(any(hits))
  # Object size is independent of n: a 20x larger training set must not
  # grow the serialized predictor materially (stored training data would
  # scale linearly with n).
  d_big <- .make_fin_data(n = 2000L, p = 12L, seed = 21L)
  fit_big <- stabl_fit(d_big$x, d_big$yg, lambda_grid = grid3,
                       family = "gaussian", n_bootstraps = 15L,
                       random_state = 1L)
  fin_big <- stabl_finalize(fit_big, d_big$x, d_big$yg,
                            preprocessing = c("center", "scale"),
                            random_state = 1L)
  tmp_big <- tempfile(fileext = ".rds")
  saveRDS(fin_big, tmp_big, compress = FALSE)
  expect_lt(file.size(tmp_big), file.size(tmp) * 1.5)
})

test_that("stored preprocessing is applied, never refit, at predict time", {
  f <- .make_fin()
  fin <- stabl_finalize(f$fit_g, f$d$x, f$d$yg,
                        preprocessing = c("center", "scale"),
                        random_state = 21L)
  pp <- fin$preprocessing
  # Stored parameters are the TRAINING moments.
  expect_equal(unname(pp$params$center$means),
               unname(colMeans(f$d$x)[pp$retained_features]))
  # Shifted new data is transformed with the stored moments: compare
  # against manual application of the stored transform + stored model.
  x_new <- f$d$x[1:6, ] + 100  # distribution shift
  xt <- stabl_preprocess_apply(pp, x_new[, fin$required_columns])
  manual <- as.numeric(stats::predict(
    fin$model, newx = xt[, fin$selected_features, drop = FALSE],
    s = "lambda.min"
  )[, 1L])
  expect_equal(predict(fin, x_new), manual)
})

test_that("one-feature panel uses the packed glm coefficients", {
  # Single-signal problem so the top importance is unique; threshold
  # between the top two importances keeps exactly one feature.
  set.seed(77L)
  n <- 100L; p <- 10L
  ids <- paste0("s", seq_len(n))
  x1 <- matrix(rnorm(n * p), n, p,
               dimnames = list(ids, paste0("g", seq_len(p))))
  y1 <- setNames(3 * x1[, 1L] + rnorm(n, sd = 1), ids)
  fit0 <- stabl_fit(x1, y1, lambda_grid = grid3, family = "gaussian",
                    n_bootstraps = 30L, random_state = 77L)
  imp <- sort(get_importances(fit0), decreasing = TRUE)
  expect_gt(imp[1L], imp[2L])  # unique top score on this DGP
  thr <- mean(imp[1:2])
  fit1 <- stabl_fit(x1, y1, lambda_grid = grid3, family = "gaussian",
                    n_bootstraps = 30L, hard_threshold = thr,
                    random_state = 77L)
  expect_equal(sum(get_support(fit1)), 1L)
  fin <- stabl_finalize(fit1, x1, y1, random_state = 77L)
  expect_identical(fin$engine, "glm_gaussian")
  expect_length(fin$model$coefficients, 2L)
  sel <- fin$selected_features
  manual <- fin$model$coefficients[1L] +
    fin$model$coefficients[2L] * x1[1:5, sel]
  expect_equal(unname(predict(fin, x1[1:5, sel, drop = FALSE])),
               unname(manual))
  # Full-matrix input is rejected unless extras are explicitly dropped.
  expect_error(predict(fin, x1[1:5, ]), "unexpected column", fixed = TRUE)
  expect_equal(unname(predict(fin, x1[1:5, ], extra_columns = "drop")),
               unname(manual))
})

test_that("empty panel falls back to the intercept-only model", {
  f <- .make_fin()
  fit_e <- stabl_fit(f$d$x, f$d$yg, lambda_grid = grid3, family = "gaussian",
                     n_bootstraps = 15L, hard_threshold = 1,
                     random_state = 1L)
  fin <- stabl_finalize(fit_e, f$d$x, f$d$yg, random_state = 1L)
  expect_identical(fin$engine, "training_mean")
  expect_length(fin$selected_features, 0L)
  expect_length(fin$required_columns, 0L)
  preds <- predict(fin, f$d$x[1:7, ])
  expect_length(preds, 7L)
  expect_equal(unique(preds), mean(f$d$yg))

  fit_eb <- stabl_fit(f$d$x, f$d$yb, lambda_grid = grid3, family = "binomial",
                      n_bootstraps = 15L, hard_threshold = 1,
                      random_state = 1L)
  fin_b <- stabl_finalize(fit_eb, f$d$x, f$d$yb, random_state = 1L)
  expect_identical(fin_b$engine, "training_prevalence")
  cl <- predict(fin_b, f$d$x[1:7, ])
  expect_length(unique(as.character(cl)), 1L)
  pr <- predict(fin_b, f$d$x[1:7, ], type = "probability")
  expect_equal(unique(pr), mean(f$d$yb == fin_b$event_class))
})

test_that("validated stabl_analysis carries its CV summary, labelled", {
  d <- .make_fin_data(n = 90L, p = 15L)
  val <- stabl_analyze(
    d$x, d$yg, family = "gaussian", mode = "validated",
    preprocessing = c("center", "scale"), lambda_grid = grid3,
    outer_v = 3L, inner_v = 3L, n_bootstraps = 20L,
    random_state = 31L, n_repeats = 3L, final_refit = TRUE
  )
  fin <- stabl_finalize(val, d$x, d$yg,
                        preprocessing = c("center", "scale"),
                        random_state = 31L)
  expect_identical(fin$development$source, "stabl_analysis_validated")
  expect_identical(fin$development$evidence_type, "full_data_refit")
  expect_identical(fin$development$validation_evidence_type,
                   "cross_validated")
  expect_true(!is.null(fin$development$performance))
  expect_true(all(c("mae", "rmse", "r2") %in%
                    names(fin$development$performance$pooled)))
  expect_match(fin$development$note, "PROCEDURE", fixed = TRUE)
  expect_length(predict(fin, d$x[1:4, ]), 4L)
})

test_that("validated analysis without final_refit cannot be finalized", {
  d <- .make_fin_data(n = 60L, p = 10L)
  val <- stabl_analyze(
    d$x, d$yg, family = "gaussian", mode = "validated",
    lambda_grid = grid3, outer_v = 3L, inner_v = 3L,
    n_bootstraps = 10L, random_state = 31L, n_repeats = 2L
  )
  expect_error(stabl_finalize(val, d$x, d$yg), "final_refit", fixed = TRUE)
})

test_that("exploratory stabl_analysis finalizes with apparent labelling", {
  d <- .make_fin_data(n = 60L, p = 10L)
  expl <- suppressWarnings(stabl_analyze(
    d$x, d$yb, family = "binomial", mode = "exploratory",
    lambda_grid = grid3, n_bootstraps = 20L, random_state = 41L,
    n_repeats = 2L
  ))
  fin <- stabl_finalize(expl, d$x, d$yb, random_state = 41L)
  expect_identical(fin$development$source, "stabl_analysis_exploratory")
  expect_null(fin$development$performance)
  xb <- d$x[, fin$required_columns, drop = FALSE]
  expect_s3_class(predict(fin, xb[1:4, ]), "factor")
})

test_that("finalize guards are fail-closed", {
  f <- .make_fin()
  expect_error(stabl_finalize(list(a = 1), f$d$x, f$d$yg),
               "stabl_fit", fixed = TRUE)
  expect_error(
    stabl_finalize(f$fit_g, f$d$x, f$d$yg, event_class = "case"),
    "`event_class` is only defined for family", fixed = TRUE
  )
  # Unnamed y of the wrong length, and wrongly named y, are rejected.
  expect_error(
    stabl_finalize(f$fit_g, f$d$x, unname(f$d$yg)[1:10]),
    "same length", ignore.case = TRUE
  )
  y_bad <- f$d$yg
  names(y_bad)[1L] <- "not_a_sample"
  expect_error(
    stabl_finalize(f$fit_g, f$d$x, y_bad),
    "cover exactly", ignore.case = TRUE
  )
  # Unsupported family
  fit_fake <- f$fit_g
  fit_fake$family <- "cox"
  expect_error(stabl_finalize(fit_fake, f$d$x, f$d$yg),
               "binomial and gaussian", fixed = TRUE)
  # Selected feature absent from x
  expect_error(
    stabl_finalize(f$fit_g, f$d$x[, 1:2], f$d$yg),
    "not present|missing", ignore.case = TRUE
  )
  # type mismatch
  fin <- stabl_finalize(f$fit_g, f$d$x, f$d$yg, random_state = 21L)
  expect_error(predict(fin, f$d$x[1:2, ], type = "class"), "should be")
  fin_b <- stabl_finalize(f$fit_b, f$d$x, f$d$yb, random_state = 21L)
  expect_error(
    predict(fin_b, f$d$x[1:2, fin_b$required_columns], type = "response"),
    "should be"
  )
})

test_that("recipes preprocessing has its training template stripped", {
  skip_if_not_installed("recipes")
  f <- .make_fin()
  rec <- recipes::recipe(~ ., data = as.data.frame(f$d$x))
  rec <- recipes::step_center(rec, recipes::all_predictors())
  fin <- stabl_finalize(f$fit_g, f$d$x, f$d$yg, preprocessing = rec,
                        random_state = 21L)
  expect_null(fin$preprocessing$prepped$template)
  expect_length(predict(fin, f$d$x[1:5, ]), 5L)
})

test_that("print.stabl_predictor is concise and honest", {
  f <- .make_fin()
  fin <- stabl_finalize(f$fit_g, f$d$x, f$d$yg, random_state = 21L)
  out <- capture.output(print(fin))
  expect_true(any(grepl("<stabl_predictor: gaussian >", out, fixed = TRUE)))
  expect_true(any(grepl("full-data refit; not cross-validated", out,
                        fixed = TRUE)))
  expect_identical(invisible(print(fin)), fin)
})
