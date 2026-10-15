# WS-8 v2 — Figures from saved results only.
#
# Reads committed CSVs under results/ and writes PNGs to results/figures/.
# No model fitting here; every panel is traceable to a committed table.
#
# Usage: Rscript benchmarks/ws8_v2/figures.R   # run from the repo root

suppressPackageStartupMessages(library(ggplot2))

base <- "benchmarks/ws8_v2/results"
fig_dir <- file.path(base, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7",
               "#56B4E9", "#F0E442", "#000000")
method_cols <- setNames(okabe_ito[1:6], paste0("M", 0:5))
method_lbls <- c(M0 = "M0 meta", M1 = "M1 vote", M2 = "M2 lasso",
                 M3 = "M3 pooled", M4 = "M4 env-min", M5 = "M5 env-maj")

theme_set(theme_bw(base_size = 11) +
            theme(text = element_text(family = "Liberation Sans"),
                  panel.grid.minor = element_blank()))

save_png <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, name), p, width = w, height = h, dpi = 300,
         units = "in", bg = "white")
}

# ---- Fig 1: primary endpoint, per-cohort R@10 by method -----------------------
rep <- read.csv(file.path(base, "run", "heldout_replication.csv"))
rep10 <- subset(rep, K == 10)
rep10$method <- factor(rep10$method, levels = paste0("M", 0:5))
macro10 <- aggregate(R ~ method, rep10, mean)

p1 <- ggplot(rep10, aes(method, R, colour = method)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey60") +
  geom_point(position = position_jitter(width = 0.15, seed = 1),
             size = 1.8, alpha = 0.75) +
  geom_point(data = macro10, aes(method, R), shape = 18, size = 4,
             colour = "black") +
  scale_colour_manual(values = method_cols, labels = method_lbls, guide = "none") +
  scale_x_discrete(labels = method_lbls) +
  labs(x = NULL, y = "Held-out signed replication R (K = 10)") +
  annotate("text", x = 0.6, y = max(rep10$R), label = "diamond = macro mean",
           hjust = 0, size = 3, colour = "grey30")
save_png(p1, "fig1_replication_k10.png", 6.5, 4)

# ---- Fig 2: K sensitivity -----------------------------------------------------
rep$method <- factor(rep$method, levels = paste0("M", 0:5))
macroK <- aggregate(R ~ method + K, rep, mean)

p2 <- ggplot(macroK, aes(K, R, colour = method, group = method)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey60") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_vline(xintercept = 10, linetype = 2, colour = "grey75") +
  scale_colour_manual(values = method_cols, labels = method_lbls, name = NULL) +
  scale_x_continuous(breaks = c(5, 10, 15, 20)) +
  labs(x = "Panel size K", y = "Macro held-out signed replication R(K)") +
  theme(legend.position = "bottom")
save_png(p2, "fig2_k_sensitivity.png", 6, 4.2)

# ---- Fig 3: predictive transportability (M4, macro MAE by panel) --------------
pm <- read.csv(file.path(base, "run_m4b", "prediction_metrics.csv"))
panel_lbls <- c(topk = "M4 top-10", native = "M4 native",
                all59 = "All 59", median = "Median age")
pm$panel <- factor(pm$panel, levels = names(panel_lbls))
macro_mae <- aggregate(mae ~ panel, pm, mean)

p3 <- ggplot(pm, aes(panel, mae, fill = panel)) +
  geom_point(position = position_jitter(width = 0.15, seed = 1),
             size = 1.8, alpha = 0.7, colour = "grey25") +
  geom_point(data = macro_mae, aes(panel, mae), shape = 18, size = 4) +
  scale_x_discrete(labels = panel_lbls) +
  scale_fill_manual(values = okabe_ito[c(1, 2, 3, 7)], guide = "none") +
  labs(x = NULL, y = "Held-out age prediction MAE (years)")
save_png(p3, "fig3_prediction_mae.png", 5.5, 4)

# ---- Fig 4: S0 null probe — selection count under the null --------------------
null <- read.csv(file.path(base, "null", "null_results.csv"))
# Real-data comparator from the committed same-design artifact (review F2).
ref <- read.csv(file.path(base, "null", "real_reference_scores.csv"))
real_n <- sum(ref$native)

p4 <- ggplot(null, aes(perm, n_selected)) +
  geom_hline(yintercept = real_n, linetype = 2, colour = "#D55E00",
             linewidth = 0.9) +
  geom_point(size = 2.5, colour = "#0072B2") +
  annotate("text", x = 1, y = real_n - 2.2,
           label = sprintf("real data: %d features", real_n),
           hjust = 0, size = 3.2, colour = "#D55E00") +
  scale_x_continuous(breaks = seq(0, 50, 10)) +
  scale_y_continuous(limits = c(0, 62)) +
  labs(x = "Within-cohort age permutation",
       y = "Features selected (min-agg >= 0.5)")
save_png(p4, "fig4_null_probe.png", 6, 3.8)

# ---- Fig 5: reference populations, per-cohort Spearman rho --------------------
rho <- read.csv(file.path(base, "biology", "per_cohort_rho.csv"),
                row.names = 1, check.names = FALSE)
ref <- c("CT:CD8+ T Naive", "CT:CD8+ RTE", "CT:CD4+ RTE", "CT:CD4+ Temra",
         "CT:MAIT", "CT:CD8+ Temra", "CT:CD8+ Tem GZMK+", "CT:CD8+ Tem GZMB+")
ref <- intersect(ref, rownames(rho))
rho_df <- data.frame(
  feature = factor(rep(ref, ncol(rho)), levels = rev(ref)),
  cohort  = rep(colnames(rho), each = length(ref)),
  rho     = as.vector(unlist(rho[ref, ])))

p5 <- ggplot(rho_df, aes(cohort, feature, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.6) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
save_png(p5, "fig5_reference_populations.png", 7, 4.5)

cat("Figures written to", fig_dir, "\n")
