# WS-8 follow-up: does per-cohort recalibration close the LOCO slope gap?
#
# The LOCO predictions are already out-of-cohort. Within each held-out
# cohort we cross-validate the recalibrator itself (5-fold x 5 repeats):
# fit the calibration map on 4 folds, apply to the 5th, and evaluate
# out-of-fold. Two maps are compared:
#   linear   : observed ~ predicted (least squares)
#   isotonic : monotone map via stats::isoreg + stepfun interpolation
# Metrics: calibration slope of lm(prediction ~ observed) and MAE, before
# (raw LOCO predictions) and after (out-of-fold recalibrated predictions).
#
# Run from the repository root:
#   R_LIBS=<rlib> Rscript benchmarks/ws8_recalibration.R

obj <- readRDS(file.path(
  Sys.getenv("WS8_OUT_DIR", "/workspace/ws8_out"), "ws8_objects.rds"
))

suppressPackageStartupMessages(library(ggplot2))

collapse_donors <- function(l) {
  data.frame(
    cohort    = as.character(l$cohort),
    observed  = as.numeric(unlist(l$obs)),
    predicted = as.numeric(unlist(l$pred)),
    stringsAsFactors = FALSE
  )
}

recal_cv <- function(d, reps = 5L, k = 5L) {
  n <- nrow(d)
  slope_pre <- unname(coef(lm(predicted ~ observed, data = d))[2])
  mae_pre   <- mean(abs(d$predicted - d$observed))
  per_rep <- lapply(seq_len(reps), function(r) {
    set.seed(20261013 + r)
    folds <- sample(rep(seq_len(k), length.out = n))
    pr_lin <- pr_iso <- rep(NA_real_, n)
    for (f in seq_len(k)) {
      tr <- folds != f
      te <- folds == f
      fit <- lm(observed ~ predicted, data = d[tr, ])
      pr_lin[te] <- predict(fit, newdata = d[te, ])
      iso <- stats::isoreg(d$predicted[tr], d$observed[tr])
      iso_fun <- as.stepfun(iso)
      pr_iso[te] <- iso_fun(d$predicted[te])
    }
    c(
      mae_lin   = mean(abs(pr_lin - d$observed)),
      mae_iso   = mean(abs(pr_iso - d$observed)),
      slope_lin = unname(coef(lm(pr_lin ~ d$observed))[2]),
      slope_iso = unname(coef(lm(pr_iso ~ d$observed))[2])
    )
  })
  m <- do.call(rbind, per_rep)
  r <- suppressWarnings(cor(d$observed, d$predicted))
  data.frame(
    mae_pre       = mae_pre,
    slope_pre     = slope_pre,
    pearson_r     = r,
    r_squared     = r^2,
    mae_lin       = mean(m[, "mae_lin"]),
    mae_iso       = mean(m[, "mae_iso"]),
    slope_lin     = mean(m[, "slope_lin"]),
    slope_lin_sd  = stats::sd(m[, "slope_lin"]),
    slope_iso     = mean(m[, "slope_iso"]),
    slope_iso_sd  = stats::sd(m[, "slope_iso"]),
    stringsAsFactors = FALSE
  )
  # Identity check: out-of-fold slope after linear recalibration equals r^2
  # in expectation (recalibration removes bias, not noise).
}

cal <- collapse_donors(obj$loco_sel)
res <- do.call(rbind, lapply(sort(unique(cal$cohort)), function(cc) {
  d <- cal[cal$cohort == cc, ]
  cbind(data.frame(cohort = cc, n = nrow(d), stringsAsFactors = FALSE),
        recal_cv(d))
}))
write.csv(res, "benchmarks/ws8_recalibration.csv", row.names = FALSE)

print(res, digits = 3)

# ---- Figure: calibration slope before vs after linear recalibration --------
fp <- rbind(
  data.frame(cohort = res$cohort, stage = "before recalibration",
             slope = res$slope_pre, stringsAsFactors = FALSE),
  data.frame(cohort = res$cohort, stage = "after linear recalibration",
             slope = res$slope_lin, stringsAsFactors = FALSE)
)
fp <- rbind(fp,
  data.frame(cohort = res$cohort, stage = "r-squared identity",
             slope = res$r_squared, stringsAsFactors = FALSE))
fp$cohort <- factor(fp$cohort, levels = res$cohort[order(res$slope_pre)])
fp$stage  <- factor(fp$stage, levels = c("before recalibration",
                                         "after linear recalibration",
                                         "r-squared identity"))

theme_set(theme_minimal(base_size = 11))
p <- ggplot(fp, aes(y = cohort)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
  geom_line(data = fp[fp$stage != "r-squared identity", ],
            aes(x = slope, group = cohort), color = "grey70",
            linewidth = 0.5) +
  geom_point(aes(x = slope, color = stage, shape = stage), size = 2.6) +
  scale_color_manual(values = c("before recalibration" = "#9aa5b1",
                                "after linear recalibration" = "#0279EE",
                                "r-squared identity" = "#D55E00"),
                     name = NULL) +
  scale_shape_manual(values = c("before recalibration" = 16,
                                "after linear recalibration" = 17,
                                "r-squared identity" = 5),
                     name = NULL) +
  scale_x_continuous(limits = c(0, 1.25)) +
  labs(x = "Calibration slope of LOCO predictions (ideal = 1)",
       y = "Held-out cohort") +
  theme(legend.position = "top",
        text = element_text(family = "Liberation Sans"),
        panel.grid.minor = element_blank())
dir.create("vignettes/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("vignettes/figures/ws8_recalibration.png", p,
       width = 6.6, height = 4.0, dpi = 300)

cat("recalibration table and figure written\n")
