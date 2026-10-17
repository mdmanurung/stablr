# Workstream C analysis: figures + tables for the all-learner comparison.
# Reads /mnt/results/learner_comparison/{lb_per_seed.csv,lb_summary.csv},
# writes figures + a guide-ready table into the same folder.

suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})
OUT <- "/mnt/results/learner_comparison"
ps <- read.csv(file.path(OUT, "lb_per_seed.csv"))
summ <- read.csv(file.path(OUT, "lb_summary.csv"))

themeset <- theme_minimal(base_size = 11, base_family = "Liberation Sans") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 11),
        strip.text = element_text(face = "bold"))

ps$learner <- factor(ps$learner, levels = c(
  "lasso", "elastic_net", "adaptive_lasso", "sparse_group_lasso",
  "splsda", "rf", "oblique_rf"))
ps$pairing <- ifelse(ps$art == "random_permutation", "RP", "KE")
ps$dgp_lab <- factor(ps$dgp, levels = c("independent", "block"),
                     labels = c("Independent features", "Block-correlated (AR 0.6, blocks of 10)"))

# ---- Fig 1: TPR vs FDP, per learner x pairing, faceted by DGP ---------------
means <- ps %>%
  group_by(dgp_lab, learner, pairing) %>%
  summarise(mean_tpr = mean(tpr), mean_fdp = mean(fdp), .groups = "drop")

p1 <- ggplot(ps, aes(fdp, tpr, color = pairing, shape = pairing)) +
  geom_point(alpha = 0.35, size = 1.6, position = position_jitter(width = 0.02, height = 0.02)) +
  geom_point(data = means, aes(x = mean_fdp, y = mean_tpr, fill = pairing),
             inherit.aes = FALSE, size = 3.4,
             shape = 21, color = "black", stroke = 0.6) +
  facet_grid(dgp_lab ~ learner) +
  scale_color_manual(values = c(RP = "#0072B2", KE = "#D55E00")) +
  scale_fill_manual(values = c(RP = "#0072B2", KE = "#D55E00")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Realized FDP", y = "TPR (10 planted signals)",
       title = "All-learner comparison: selection quality by learner and artificial pairing",
       subtitle = "Small points = seeds; large points = means. n=200, p=150, 10 signals, 100 bootstraps") +
  themeset +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7))
ggsave(file.path(OUT, "figA_tpr_fdp.png"), p1, width = 12, height = 5.2, dpi = 200)

# ---- Fig 2: calibration ratio ------------------------------------------------
p2 <- ggplot(ps, aes(learner, calib, fill = pairing)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  stat_summary(fun = median, geom = "col", position = position_dodge(0.7), width = 0.62,
               color = "black", linewidth = 0.3, na.rm = TRUE) +
  stat_summary(aes(group = interaction(pairing, learner)), fun = median,
               geom = "point", position = position_dodge(0.7), size = 1.4, na.rm = TRUE) +
  facet_grid(~ dgp_lab) +
  scale_fill_manual(values = c(RP = "#0072B2", KE = "#D55E00"),
                    labels = c(RP = "random permutation", KE = "knockoff_equi")) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = "Realized FDP / reported FDP+ (median, seeds)",
       title = "FDP+ calibration by learner and pairing",
       subtitle = "Dashed line = perfect calibration; bars = median over seeds") +
  themeset +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 8))
ggsave(file.path(OUT, "figB_calibration.png"), p2, width = 10, height = 4.6, dpi = 200)

# ---- Fig 3: stability (pairwise Jaccard) + runtime ---------------------------
summ$learner <- factor(summ$learner, levels = c(
  "lasso", "elastic_net", "adaptive_lasso", "sparse_group_lasso",
  "splsda", "rf", "oblique_rf"))
summ$pairing <- ifelse(summ$art == "random_permutation", "RP", "KE")
summ$dgp_lab <- factor(summ$dgp, levels = c("independent", "block"),
                       labels = c("Independent", "Block-correlated"))

p3 <- ggplot(summ, aes(learner, mean_jaccard, fill = pairing)) +
  geom_col(position = position_dodge(0.7), width = 0.62,
           color = "black", linewidth = 0.3) +
  facet_grid(~ dgp_lab) +
  scale_fill_manual(values = c(RP = "#0072B2", KE = "#D55E00"),
                    labels = c(RP = "random permutation", KE = "knockoff_equi")) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = "Mean pairwise Jaccard across seeds",
       title = "Selection stability (higher = more reproducible support)",
       subtitle = "Mean pairwise Jaccard of selected supports across seeds") +
  themeset +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 8))
ggsave(file.path(OUT, "figC_stability.png"), p3, width = 10, height = 4.2, dpi = 200)

p4 <- ggplot(summ, aes(learner, mean_secs, fill = pairing)) +
  geom_col(position = position_dodge(0.7), width = 0.62,
           color = "black", linewidth = 0.3) +
  facet_grid(~ dgp_lab, scales = "free_y") +
  scale_fill_manual(values = c(RP = "#0072B2", KE = "#D55E00"),
                    labels = c(RP = "random permutation", KE = "knockoff_equi")) +
  scale_y_continuous(trans = "log10") +
  labs(x = NULL, y = "Mean runtime per fit (seconds, log)",
       title = "Runtime per stabl_fit call (100 bootstraps, n=200, p=150)",
       subtitle = "Single core per fit; rf/oblique_rf at 5 seeds, others 10") +
  themeset +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 8))
ggsave(file.path(OUT, "figD_runtime.png"), p4, width = 10, height = 4.2, dpi = 200)

# ---- Guide-ready wide table ---------------------------------------------------
guide <- summ %>%
  select(dgp, learner, pairing, n_seeds, mean_tpr, mean_fdp, mean_n_sel,
         median_reported, median_calib, mean_jaccard, mean_secs) %>%
  arrange(dgp, learner, pairing)
write.csv(guide, file.path(OUT, "guide_table.csv"), row.names = FALSE)
cat("figures + guide_table.csv written\n")
print(guide, digits = 2, row.names = FALSE)
