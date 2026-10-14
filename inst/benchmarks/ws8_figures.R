# WS-8 showcase figures (PNG) for the aging case-study vignette.
#
# Reads the saved result objects written by ws8_aging_showcase.R
# (WS8_OUT_DIR/ws8_objects.rds) and writes three PNG files into the
# vignette figures folder. Run from the repository root:
#
#   R_LIBS=<rlib> Rscript inst/benchmarks/ws8_figures.R
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

# Phylo-friendly, colorblind-safe settings
theme_set(theme_minimal(base_size = 11))
theme_update(
  text = element_text(family = "Liberation Sans"),
  panel.grid.minor = element_blank()
)
sel_col <- "#0279EE"   # selected (min-agg >= 0.5)
uns_col <- "#9aa5b1"   # not selected
coh_cols <- c("#0279EE", "#FF9400", "#75A025", "#E9ED4C",
              "#FD9BED", "#00BFC4", "#D55E00", "#6A3D9A")

sel_rf <- obj$sel_rf
cor_mat <- obj$cor_mat
loco_sel <- obj$loco_sel

paper_consistent <- c("CT:MAIT", "CT:CD4+ RTE", "CT:CD8+ RTE",
                      "CT:CD8+ T Naive", "CT:CD8+ Temra", "CT:CD4+ Temra",
                      "CT:CD8+ Tem GZMK+", "CT:CD8+ Tem GZMB+")

# ---- Figure 1: min-agg score ranking ----------------------------------------
rank_df <- sel_rf[order(sel_rf$min_score, sel_rf$mean_score), ]
rank_df$label <- sub("^CT:", "", rank_df$feature)
rank_df$selected <- factor(rank_df$min_score >= 0.5, levels = c(TRUE, FALSE))
rank_df$paper <- rank_df$feature %in% paper_consistent

p1 <- ggplot(rank_df, aes(x = min_score, y = reorder(label, min_score))) +
  geom_col(aes(fill = selected), width = 0.72) +
  geom_point(data = subset(rank_df, paper), aes(x = min_score, y = label),
             shape = 4, size = 1.9, stroke = 0.7, color = "black") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "#D55E00") +
  scale_fill_manual(values = c(`TRUE` = sel_col, `FALSE` = uns_col),
                    labels = c(`TRUE` = "selected (min-agg >= 0.5)",
                               `FALSE` = "not selected"), name = NULL) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Minimum cross-cohort stability score",
       y = NULL,
       shape = "Paper-consistent") +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 6.4),
        panel.grid.major.y = element_blank())
ggsave(file.path(fig_dir, "ws8_min_scores.png"), p1,
       width = 6.2, height = 8.6, dpi = 200)

# ---- Figure 2: LOCO calibration scatter (selected-feature model) ------------
# unlist() keeps the AsIs class from the saved objects; strip it so ggplot2
# treats the columns as plain numerics.
cal <- data.frame(
  cohort = as.character(loco_sel$cohort),
  observed = as.numeric(unlist(loco_sel$obs)),
  predicted = as.numeric(unlist(loco_sel$pred)),
  stringsAsFactors = FALSE
)
coh_order <- names(sort(table(cal$cohort), decreasing = TRUE))
cal$cohort <- factor(cal$cohort, levels = coh_order)

p2 <- ggplot(cal, aes(x = observed, y = predicted, color = cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40") +
  geom_point(alpha = 0.35, size = 1.1, shape = 16) +
  facet_wrap(~cohort, nrow = 2) +
  scale_color_manual(values = coh_cols) +
  labs(x = "Observed age (years)", y = "LOCO predicted age (years)") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "ws8_loco_calibration.png"), p2,
       width = 8.4, height = 4.0, dpi = 200)

# ---- Figure 3: consistency heatmap (selected features x cohorts) ------------
selected <- sel_rf$feature[sel_rf$selected_minagg]
hm <- cor_mat[selected, , drop = FALSE]
feat_order <- selected[order(rowMeans(hm), decreasing = TRUE)]
hm <- hm[feat_order, , drop = FALSE]
hm_df <- data.frame(
  feature = factor(sub("^CT:", "", rownames(hm)), levels = sub("^CT:", "", feat_order)),
  cohort = factor(rep(colnames(hm), each = nrow(hm)), levels = colnames(hm)),
  rho = as.vector(hm),
  stringsAsFactors = FALSE
)

p3 <- ggplot(hm_df, aes(x = cohort, y = feature, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#D55E00", mid = "#F7F7F7", high = "#0279EE",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6.8),
        panel.grid = element_blank())
ggsave(file.path(fig_dir, "ws8_consistency_heatmap.png"), p3,
       width = 6.4, height = 7.2, dpi = 200)

cat("figures written to", fig_dir, "\n")
