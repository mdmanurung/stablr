# WS-8 v3 — Figures for the candidate benchmark (PNG, Liberation Sans).
# Reads assembled tables (run assemble_v3.R first; env WS8_V3_TAB / WS8_V3_FIG).
suppressPackageStartupMessages({library(ggplot2)})

v3_dir <- Sys.getenv("WS8_V3_DIR", "benchmarks/ws8_v3")
tab <- Sys.getenv("WS8_V3_TAB", file.path(v3_dir, "results", "tables"))
fig <- Sys.getenv("WS8_V3_FIG", file.path(v3_dir, "results", "figures"))
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

library(scales)
theme_set(theme_bw(base_size = 11, base_family = "Liberation Sans") +
            theme(panel.grid.minor = element_blank(),
                  legend.position = "bottom"))
pal <- colorRampPalette(c("#0279EE", "#75A025", "#FF9400", "#FD9BED",
                          "#7c5cc4", "#8c8c8c"))(13)

METHODS <- c("M0","M1","M2","M3","M4","M5","M6","M7","M7c","M8","M9","M9s","M9b")
rep_df <- read.csv(file.path(tab, "heldout_replication_v3.csv"))
macro <- read.csv(file.path(tab, "macro_R_by_K_v3.csv"), row.names = 1)
rec_df <- read.csv(file.path(tab, "panel_recovery_per_fold_v3.csv"))

ord <- rownames(macro)[order(-macro[, "K10"])]
meth_f <- factor(rep_df$method, levels = ord)

# ---- Fig 1: per-fold R@10 by method -------------------------------------------
d10 <- rep_df[rep_df$K == 10, ]
d10$method <- factor(d10$method, levels = ord)
mu <- aggregate(R ~ method, d10, mean)
p1 <- ggplot(d10, aes(method, R)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey60") +
  geom_jitter(aes(colour = method), width = 0.15, height = 0, size = 1.8,
              alpha = 0.75, show.legend = FALSE) +
  geom_crossbar(data = mu, aes(method, R, ymin = R, ymax = R),
                width = 0.55, fatten = 1.2, colour = "black") +
  scale_colour_manual(values = setNames(pal, ord)) +
  labs(x = NULL, y = "Held-out signed-rho replication R@10") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig, "fig1_r10_by_method.png"), p1, width = 7.5, height = 4.2,
       dpi = 200)

# ---- Fig 2: macro R across K ----------------------------------------------------
mk <- data.frame(method = rep(rownames(macro), ncol(macro)),
                 K = rep(as.integer(sub("K", "", colnames(macro))),
                         each = nrow(macro)),
                 R = as.vector(as.matrix(macro)))
mk$method <- factor(mk$method, levels = ord)
p2 <- ggplot(mk, aes(K, R, colour = method, group = method)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  scale_colour_manual(values = setNames(pal, ord)) +
  scale_x_continuous(breaks = c(5, 10, 15, 20)) +
  labs(x = "Panel size K", y = "Macro signed-rho replication R(K)",
       colour = NULL) +
  guides(colour = guide_legend(nrow = 2))
ggsave(file.path(fig, "fig2_macro_R_by_K.png"), p2, width = 7, height = 4.5,
       dpi = 200)

# ---- Fig 3: P20 panel recovery (co-primary) -------------------------------------
r10 <- rec_df[rec_df$K == 10, ]
r10$method <- factor(r10$method, levels = ord)
rmu <- aggregate(recovery ~ method, r10, mean)
p3 <- ggplot(r10, aes(method, recovery)) +
  geom_jitter(aes(colour = method), width = 0.15, height = 0, size = 1.8,
              alpha = 0.75, show.legend = FALSE) +
  geom_crossbar(data = rmu, aes(method, recovery, ymin = recovery,
                                ymax = recovery),
                width = 0.55, fatten = 1.2, colour = "black") +
  scale_colour_manual(values = setNames(pal, ord)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Recovery of paper 20-population clock panel (top-10)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig, "fig3_p20_recovery.png"), p3, width = 7.5, height = 4.2,
       dpi = 200)

# ---- Fig 4: null probes for M7/M8 ------------------------------------------------
null_f <- file.path(tab, "null_summary_v3.csv")
null_scores <- file.path(v3_dir, "results", "null_v3", "null_scores_v3.csv")
if (file.exists(null_f) && file.exists(null_scores)) {
  ns <- read.csv(null_scores)
  ns <- ns[ns$perm > 0, ]
  real <- read.csv(null_f)
  p4 <- ggplot(ns, aes(score)) +
    geom_histogram(bins = 40, fill = "#8c8c8c", colour = "white") +
    geom_vline(data = real, aes(xintercept = real_max_score),
               colour = "#D55E00", linewidth = 1) +
    facet_wrap(~method, scales = "free") +
    labs(x = "Stability score", y = "Features x permutations",
         caption = "Grey: within-cohort age-permutation null. Orange: real-data max score.")
  ggsave(file.path(fig, "fig4_null_m7_m8.png"), p4, width = 8, height = 3.8,
         dpi = 200)
}

# ---- Fig 5: prediction MAE (top-K panel) -----------------------------------------
met <- read.csv(file.path(tab, "prediction_summary_v3.csv"))
mt <- met[met$panel == "topk", ]
mt$method <- factor(mt$method, levels = ord)
mb <- met[met$panel == "median", ]
p5 <- ggplot(mt, aes(method, mae, fill = method)) +
  geom_col(show.legend = FALSE) +
  geom_hline(yintercept = mean(mb$mae), linetype = 2, colour = "grey40") +
  scale_fill_manual(values = setNames(pal, ord)) +
  labs(x = NULL, y = "Held-out macro MAE, years (top-10 panel)",
       caption = "Dashed: training-median baseline.") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig, "fig5_prediction_mae.png"), p5, width = 7.5, height = 4.2,
       dpi = 200)

cat("Figures written to ", fig, "\n")
