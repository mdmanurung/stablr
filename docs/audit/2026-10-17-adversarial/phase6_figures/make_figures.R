# Phase 6 figures — four PNGs from the extracted curve records
suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(RColorBrewer)
})
themeset <- theme_minimal(base_size = 11, base_family = "Liberation Sans") +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom", legend.box = "vertical")
okabe <- c(RP = "#0072B2", KE = "#D55E00", null = "#CC79A7",
           signal = "#009E73", corr = "#E69F00", bound = "#000000")
dir.create("/mnt/results/phase6_figures", showWarnings = FALSE)

curves <- readRDS("/workspace/phase6_out/curves.rds")
res <- read.csv("/workspace/phase6_out/results.csv")
thr <- seq(0, 0.99, by = 0.01)

curve_df <- bind_rows(lapply(curves, function(r)
  data.frame(design = r$design, artificial = r$artificial, seed = r$seed,
             threshold = thr, fdp = pmin(r$FDRs, 2.5))))
curve_df$art <- ifelse(curve_df$artificial == "random_permutation", "RP", "KE")
curve_df$design_lab <- c(D1_sparse = "D1: sparse signal (indep.)",
                         D2_null = "D2: pure null",
                         D3_correlated = "D3: block-correlated")[curve_df$design]

# ---- Fig 1: FDP+ curves overlaid, null vs signal seeds ----------------------
# Per-threshold medians of jumpy integer-count ratios oscillate (comb artifact);
# summarize across seeds with a 10-90% quantile ribbon + median line instead.
qdf <- curve_df %>%
  group_by(design_lab, art, threshold) %>%
  summarise(lo = quantile(fdp, 0.10), med = median(fdp),
            hi = quantile(fdp, 0.90), .groups = "drop")
p1 <- ggplot(qdf, aes(threshold, med)) +
  geom_hline(yintercept = 1, linetype = 3, colour = "grey40") +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = design_lab), alpha = 0.28, colour = NA) +
  geom_line(aes(colour = design_lab), linewidth = 0.9) +
  facet_grid(art ~ ., labeller = labeller(art = c(RP = "random permutation",
                                                  KE = "knockoff (equi)"))) +
  scale_colour_manual(values = c("D2: pure null" = okabe[["null"]],
                                 "D1: sparse signal (indep.)" = okabe[["signal"]],
                                 "D3: block-correlated" = okabe[["corr"]]), name = NULL) +
  scale_fill_manual(values = c("D2: pure null" = okabe[["null"]],
                               "D1: sparse signal (indep.)" = okabe[["signal"]],
                               "D3: block-correlated" = okabe[["corr"]]), name = NULL) +
  coord_cartesian(ylim = c(0, 2.5)) +
  labs(x = "stability-score threshold", y = "FDP+ estimate (capped at 2.5)",
       title = "FDP+ curves: the argmin hunts the luckiest crossing",
       subtitle = "median across 20 seeds; ribbon = 10-90% across seeds; dotted = 1.0") +
  themeset
ggsave("/mnt/results/phase6_figures/fig1_fdp_curves_null_vs_signal.png", p1,
       width = 8, height = 6, dpi = 200)

# ---- Fig 2: calibration scatter ---------------------------------------------
res2 <- subset(res, is.na(error) | error == "NA")
res2$design_lab <- c(D1_sparse = "D1: sparse signal (indep.)",
                     D2_null = "D2: pure null",
                     D3_correlated = "D3: block-correlated")[res2$design]
res2$artificial_lab <- ifelse(res2$artificial == "random_permutation",
                              "random permutation", "knockoff (equi)")
p2 <- ggplot(res2, aes(reported_min_fdr, realized_fdp)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_abline(slope = 2, intercept = 0, linetype = 3, colour = "grey60") +
  geom_point(aes(colour = design_lab, shape = artificial_lab), alpha = 0.75, size = 2) +
  scale_colour_manual(values = c("D2: pure null" = okabe[["null"]],
                                 "D1: sparse signal (indep.)" = okabe[["signal"]],
                                 "D3: block-correlated" = okabe[["corr"]]), name = NULL) +
  scale_shape_manual(values = c("random permutation" = 16, "knockoff (equi)" = 17), name = NULL) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "reported min_fdr_", y = "realized FDP of selected support",
       title = "Calibration: realized FDP vs reported min_fdr_",
       subtitle = "dashed = perfect calibration; dotted = 2x anti-conservatism") +
  themeset
ggsave("/mnt/results/phase6_figures/fig2_calibration_scatter.png", p2,
       width = 8, height = 6, dpi = 200)

# ---- Fig 3: score distributions, null vs signal seed ------------------------
score_df <- bind_rows(lapply(curves, function(r)
  data.frame(design = r$design, artificial = r$artificial, seed = r$seed,
             kind = rep(c("real", "artificial"), c(length(r$scores_real), length(r$scores_art))),
             score = c(r$scores_real, r$scores_art))))
score_df <- subset(score_df, artificial == "knockoff_equi" & seed %in% 1:5)
score_df$panel <- ifelse(score_df$design == "D2_null", "D2: pure null (seed 1-5)",
                         ifelse(score_df$design == "D1_sparse", "D1: sparse signal (seed 1-5)", NA))
score_df <- subset(score_df, !is.na(panel))
p3 <- ggplot(score_df, aes(score, fill = kind)) +
  geom_density(alpha = 0.45, colour = NA) +
  facet_wrap(~ panel) +
  scale_fill_manual(values = c(real = okabe[["signal"]], artificial = okabe[["null"]]), name = NULL) +
  labs(x = "max-over-lambda stability score", y = "density",
       title = "Why nulls select features: real vs artificial score distributions",
       subtitle = "knockoff (equi); the argmin threshold lands where the real bulk outnumbers the artificial bulk by chance") +
  themeset
ggsave("/mnt/results/phase6_figures/fig3_score_distributions.png", p3,
       width = 8, height = 4, dpi = 200)

# ---- Fig 4: empirical tail probability vs Markov bound ----------------------
d2 <- subset(res2, design == "D2_null")
deltas <- seq(0.1, 5, by = 0.1)
tail_df <- data.frame(delta = deltas,
  empirical = sapply(deltas, function(d)
    mean(d2$realized_fdp >= (1 + d) * d2$reported_min_fdr)),
  bound = 1 / (1 + deltas))
p4 <- ggplot(tail_df, aes(delta, empirical)) +
  geom_line(aes(colour = "empirical (D2 null, 40 runs)"), linewidth = 1) +
  geom_line(aes(y = bound, colour = "Markov bound 1/(1+delta)"), linetype = 2, linewidth = 1) +
  scale_colour_manual(values = c("empirical (D2 null, 40 runs)" = okabe[["RP"]],
                                 "Markov bound 1/(1+delta)" = okabe[["bound"]]), name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = expression(delta), y = "P( FDP >= (1+delta) * q+ )",
       title = "The F-001 bound is nearly tight at small delta",
       subtitle = "the observed anti-conservatism is within theory at every delta") +
  themeset
ggsave("/mnt/results/phase6_figures/fig4_markov_bound.png", p4,
       width = 7, height = 5, dpi = 200)

cat("figures written\n")
