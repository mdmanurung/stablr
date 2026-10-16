# WS-8 showcase figures and statistics (PNG) for the aging case-study vignette.
#
# Reads the saved result objects written by ws8_aging_showcase.R
# (WS8_OUT_DIR/ws8_objects.rds) and writes four PNG files plus two CSV
# statistics tables into the vignette figures folder. Run from the
# repository root:
#
#   R_LIBS=<rlib> Rscript benchmarks/ws8_figures.R
#
# The figures present an application case study. They are NOT calibration
# evidence.

obj <- readRDS(file.path(
  Sys.getenv("WS8_OUT_DIR", "/workspace/ws8_out"), "ws8_objects.rds"
))

fig_dir <- "vignettes/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
})

# Colorblind-safe settings
theme_set(theme_minimal(base_size = 11))
theme_update(
  text = element_text(family = "Liberation Sans"),
  panel.grid.minor = element_blank()
)
sel_col <- "#0279EE"
uns_col <- "#9aa5b1"
coh_cols <- c("#0279EE", "#FF9400", "#75A025", "#E9ED4C",
              "#FD9BED", "#00BFC4", "#D55E00", "#6A3D9A")

env_rf   <- obj$env_rf
sel_rf   <- obj$sel_rf
cor_mat  <- obj$cor_mat
loco_all <- obj$loco_all
loco_sel <- obj$loco_sel
n_boot   <- 200L

paper_consistent <- c("CT:MAIT", "CT:CD4+ RTE", "CT:CD8+ RTE",
                      "CT:CD8+ T Naive", "CT:CD8+ Temra", "CT:CD4+ Temra",
                      "CT:CD8+ Tem GZMK+", "CT:CD8+ Tem GZMB+")

# ---- Uncertainty on the min-agg selection frequency --------------------------
# The per-cohort stability score is a selection frequency over n_boot
# bootstraps. For each feature, the binding cohort is the one that attains
# the minimum score; its exact (Clopper-Pearson) 95% interval is a
# conservative uncertainty band for the min-agg score.
env_scores <- env_rf$env_scores
cp_ci <- function(k, n) {
  bt <- stats::binom.test(max(0L, k), n)
  as.numeric(bt$conf.int)
}
feats <- rownames(env_scores)
ci_list <- lapply(feats, function(f) {
  v <- env_scores[f, ]
  j <- which.min(v)
  cp_ci(round(v[j] * n_boot), n_boot)
})
ci_df <- do.call(rbind, ci_list)
colnames(ci_df) <- c("ci_low", "ci_high")
rank_df <- data.frame(
  feature = feats,
  min_score = apply(env_scores, 1, min),
  mean_score = apply(env_scores, 1, mean),
  ci_df,
  stringsAsFactors = FALSE
)
rank_df$selected <- rank_df$min_score >= 0.5
rank_df$paper <- rank_df$feature %in% paper_consistent

# ---- Sign-consistency test per feature ---------------------------------------
# k = number of cohorts sharing the dominant correlation sign; exact
# two-sided binomial test against p = 0.5.
sign_test_p <- function(v) {
  k <- max(sum(v < 0), sum(v > 0))
  min(1, 2 * sum(stats::dbinom(k:8, 8, 0.5)))
}
rank_df$sign_k <- vapply(feats, function(f) {
  v <- cor_mat[f, ]
  max(sum(v < 0), sum(v > 0))
}, integer(1))
rank_df$sign_p <- vapply(feats, function(f) sign_test_p(cor_mat[f, ]),
                         numeric(1))
rank_df$mean_rho <- vapply(feats, function(f) mean(cor_mat[f, ]),
                           numeric(1))

# ---- MAE bootstrap confidence intervals --------------------------------------
# Donors are resampled within each held-out cohort (2000 replicates,
# percentile intervals). The saved LOCO objects carry one row per donor.
collapse_donors <- function(l) {
  data.frame(
    cohort = as.character(l$cohort),
    observed = as.numeric(unlist(l$obs)),
    predicted = as.numeric(unlist(l$pred)),
    stringsAsFactors = FALSE
  )
}
mae_ci_table <- function(l, reps = 2000L) {
  d <- collapse_donors(l)
  d$abs_err <- abs(d$predicted - d$observed)
  do.call(rbind, lapply(sort(unique(d$cohort)), function(cc) {
    e <- d$abs_err[d$cohort == cc]
    n <- length(e)
    bs <- replicate(reps, mean(sample(e, n, replace = TRUE)))
    data.frame(cohort = cc, n = n, mae = mean(e),
               ci_low = quantile(bs, 0.025, names = FALSE),
               ci_high = quantile(bs, 0.975, names = FALSE),
               stringsAsFactors = FALSE)
  }))
}
mae_all <- mae_ci_table(loco_all)
mae_sel <- mae_ci_table(loco_sel)
names(mae_all)[4:5] <- c("ci_low_all", "ci_high_all")
mae_tab <- data.frame(
  cohort = mae_all$cohort, n = mae_all$n,
  mae_all = mae_all$mae, lo_all = mae_all$ci_low_all, hi_all = mae_all$ci_high_all,
  mae_sel = mae_sel$mae, lo_sel = mae_sel$ci_low, hi_sel = mae_sel$ci_high,
  stringsAsFactors = FALSE
)

# ---- Calibration statistics per cohort (selected model) ----------------------
cal <- collapse_donors(loco_sel)
cal_stats <- do.call(rbind, lapply(sort(unique(cal$cohort)), function(cc) {
  d <- cal[cal$cohort == cc, ]
  fit <- lm(predicted ~ observed, data = d)
  data.frame(cohort = cc, n = nrow(d),
             slope = unname(coef(fit)[2]),
             pearson_r = suppressWarnings(cor(d$observed, d$predicted)),
             stringsAsFactors = FALSE)
}))

# ---- Write statistics tables --------------------------------------------------
write.csv(rank_df, "benchmarks/ws8_feature_stats.csv", row.names = FALSE)
write.csv(mae_tab, "benchmarks/ws8_loco_mae_ci.csv", row.names = FALSE)
write.csv(cal_stats, "benchmarks/ws8_loco_calibration.csv", row.names = FALSE)

# ---- Figure 1: min-agg score ranking with 95% Clopper-Pearson bands ----------
rank_df$label <- sub("^CT:", "", rank_df$feature)
rank_df$sel_f <- factor(rank_df$selected, levels = c(TRUE, FALSE))
p1 <- ggplot(rank_df, aes(x = min_score, y = reorder(label, min_score))) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0,
                 color = "grey55", linewidth = 0.28) +
  geom_col(aes(fill = sel_f), width = 0.72) +
  geom_point(data = subset(rank_df, paper), aes(x = min_score, y = label),
             shape = 4, size = 2.0, stroke = 0.7, color = "black") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "#D55E00") +
  scale_fill_manual(values = c(`TRUE` = sel_col, `FALSE` = uns_col),
                    labels = c(`TRUE` = "selected (min-agg >= 0.5)",
                               `FALSE` = "not selected"), name = NULL) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Minimum cross-cohort selection frequency (95% exact CI)",
       y = NULL) +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 6.4),
        panel.grid.major.y = element_blank())
ggsave(file.path(fig_dir, "ws8_min_scores.png"), p1,
       width = 6.4, height = 8.6, dpi = 300)

# ---- Figure 2: consistency heatmap with annotated correlations ---------------
selected <- rank_df$feature[rank_df$selected]
hm <- cor_mat[selected, , drop = FALSE]
feat_order <- selected[order(rowMeans(hm), decreasing = TRUE)]
hm <- hm[feat_order, , drop = FALSE]
hm_df <- data.frame(
  feature = factor(sub("^CT:", "", rownames(hm)),
                   levels = sub("^CT:", "", feat_order)),
  cohort = factor(rep(colnames(hm), each = nrow(hm)), levels = colnames(hm)),
  rho = as.vector(hm),
  stringsAsFactors = FALSE
)
p2 <- ggplot(hm_df, aes(x = cohort, y = feature, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", rho),
                color = abs(rho) > 0.45), size = 2.1) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"),
                     guide = "none") +
  scale_fill_gradient2(low = "#D55E00", mid = "#F7F7F7", high = "#0279EE",
                       midpoint = 0, limits = c(-1, 1),
                       name = expression(rho)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6.8),
        panel.grid = element_blank())
ggsave(file.path(fig_dir, "ws8_consistency_heatmap.png"), p2,
       width = 6.6, height = 7.2, dpi = 300)

# ---- Figure 3: LOCO calibration scatter with per-cohort statistics ----------
coh_order <- names(sort(table(cal$cohort), decreasing = TRUE))
cal$cohort <- factor(cal$cohort, levels = coh_order)
cal_stats$cohort <- factor(cal_stats$cohort, levels = coh_order)
ann <- data.frame(
  cohort = cal_stats$cohort,
  label = sprintf("n = %d   slope = %.2f   r = %.2f",
                  cal_stats$n, cal_stats$slope, cal_stats$pearson_r),
  stringsAsFactors = FALSE
)
p3 <- ggplot(cal, aes(x = observed, y = predicted, color = cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40") +
  geom_point(alpha = 0.35, size = 1.1, shape = 16) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.4, size = 2.6, color = "grey20",
            inherit.aes = FALSE) +
  facet_wrap(~cohort, nrow = 2) +
  scale_color_manual(values = coh_cols) +
  coord_cartesian(xlim = c(15, 100), ylim = c(0, 105)) +
  labs(x = "Observed age (years)", y = "LOCO predicted age (years)") +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "ws8_loco_calibration.png"), p3,
       width = 8.6, height = 4.2, dpi = 300)

# ---- Figure 4: MAE forest plot with bootstrap CIs ----------------------------
fp <- rbind(
  data.frame(mae_tab[, c("cohort", "n")], model = "all 60 features",
             mae = mae_tab$mae_all, lo = mae_tab$lo_all,
             hi = mae_tab$hi_all, stringsAsFactors = FALSE),
  data.frame(mae_tab[, c("cohort", "n")], model = "24 selected features",
             mae = mae_tab$mae_sel, lo = mae_tab$lo_sel,
             hi = mae_tab$hi_sel, stringsAsFactors = FALSE)
)
fp$cohort <- factor(fp$cohort, levels = mae_tab$cohort)
fp$model <- factor(fp$model, levels = c("all 60 features",
                                        "24 selected features"))
p4 <- ggplot(fp, aes(x = mae, y = cohort, color = model)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey80") +
  geom_point(aes(shape = model), size = 2.4,
             position = position_dodge(width = 0.55)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, color = model), height = 0,
                 linewidth = 0.5,
                 position = position_dodge(width = 0.55)) +
  scale_color_manual(values = c("all 60 features" = "#9aa5b1",
                                "24 selected features" = sel_col),
                     name = NULL) +
  scale_shape_manual(values = c("all 60 features" = 16,
                                "24 selected features" = 17), guide = "none") +
  labs(x = "Leave-one-cohort-out age MAE (years, 95% bootstrap CI)",
       y = "Held-out cohort") +
  theme(legend.position = "top")
ggsave(file.path(fig_dir, "ws8_mae_forest.png"), p4,
       width = 7.0, height = 4.2, dpi = 300)

cat("figures and statistics written\n")
