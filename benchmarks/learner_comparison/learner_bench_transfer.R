# Transfer analysis: do the binomial learner rankings transfer to a
# Gaussian outcome? Reads both regime summaries, compares rankings and
# calibration patterns, writes figures + verdict table.

suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})
OUT <- "/mnt/results/learner_comparison"
bin_ <- read.csv(file.path(OUT, "lb_summary.csv"))
gau <- read.csv(file.path(OUT, "lb_gauss_summary.csv"))

bin_$regime <- "binomial"; gau$regime <- "gaussian"
both <- rbind(bin_, gau)
both$pairing <- ifelse(both$art == "random_permutation", "RP", "KE")

themeset <- theme_minimal(base_size = 11, base_family = "Liberation Sans") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 11),
        strip.text = element_text(face = "bold"))

# ---- 1. Ranking transfer (RP arms, 6 common learners) ------------------------
rp <- both %>%
  filter(pairing == "RP", learner != "splsda") %>%
  group_by(regime, dgp, learner) %>%
  summarise(mean_tpr = mean(mean_tpr), mean_fdp = mean(mean_fdp),
            median_calib = median(median_calib), .groups = "drop")

ranks <- rp %>%
  group_by(regime, dgp) %>%
  mutate(rank_tpr = rank(-mean_tpr), rank_fdp = rank(mean_fdp)) %>%
  ungroup()

transfer <- ranks %>%
  select(dgp, learner, regime, rank_tpr, rank_fdp) %>%
  pivot_wider(names_from = regime, values_from = c(rank_tpr, rank_fdp))
cor_tpr <- cor(transfer$rank_tpr_binomial, transfer$rank_tpr_gaussian,
               method = "spearman")
cor_fdp <- cor(transfer$rank_fdp_binomial, transfer$rank_fdp_gaussian,
               method = "spearman")
cat("Spearman rank correlation (RP arms, 6 learners, pooled over DGPs):\n")
cat("  TPR ranks:", round(cor_tpr, 2), " FDP ranks:", round(cor_fdp, 2), "\n")

# ---- 2. Calibration pattern transfer -----------------------------------------
calib_pat <- both %>%
  group_by(regime, learner, pairing) %>%
  summarise(median_calib = median(median_calib), .groups = "drop") %>%
  pivot_wider(names_from = regime, values_from = median_calib) %>%
  mutate(pattern =
           case_when(
             binomial <= 1.15 & gaussian <= 1.15 ~ "calibrated in both",
             binomial > 1.15 & gaussian > 1.15 ~ "anti-conservative in both",
             binomial < 0.85 & gaussian < 0.85 ~ "conservative in both",
             TRUE ~ "REGIME-DEPENDENT"))
write.csv(calib_pat, file.path(OUT, "lb_transfer_verdict.csv"), row.names = FALSE)
cat("\nCalibration pattern transfer (median calib, pooled over DGPs):\n")
print(as.data.frame(calib_pat), digits = 2, row.names = FALSE)

# ---- Fig E: calibration side-by-side ------------------------------------------
both$learner <- factor(both$learner, levels = c(
  "lasso", "elastic_net", "adaptive_lasso", "sparse_group_lasso",
  "splsda", "rf", "oblique_rf"))
pE <- ggplot(both, aes(learner, median_calib, fill = pairing)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_col(position = position_dodge(0.7), width = 0.62,
           color = "black", linewidth = 0.3) +
  facet_grid(~ regime) +
  scale_fill_manual(values = c(RP = "#0072B2", KE = "#D55E00"),
                    labels = c(RP = "random permutation", KE = "knockoff_equi")) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = "Realized FDP / reported FDP+ (median, seeds)",
       title = "FDP+ calibration: binomial vs Gaussian outcome",
       subtitle = "Dashed line = perfect calibration; sPLS-DA has no Gaussian variant") +
  themeset +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 8))
ggsave(file.path(OUT, "figE_calibration_transfer.png"), pE,
       width = 10, height = 4.6, dpi = 200)

# ---- Fig F: TPR slope chart (RP arms) -----------------------------------------
rp$learner <- factor(rp$learner, levels = c(
  "lasso", "elastic_net", "adaptive_lasso", "sparse_group_lasso",
  "rf", "oblique_rf"))
pF <- ggplot(rp, aes(regime, mean_tpr, group = learner, color = learner)) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 2.6) +
  facet_grid(~ dgp) +
  scale_color_manual(values = c(
    lasso = "#0072B2", elastic_net = "#E69F00", adaptive_lasso = "#56B4E9",
    sparse_group_lasso = "#009E73", rf = "#D55E00", oblique_rf = "#CC79A7")) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = "Mean TPR (RP pairing, 10 planted signals)",
       title = "Learner power ranking: binomial vs Gaussian outcome",
       subtitle = "Each line is one learner; parallel lines = ranking transfers") +
  themeset
ggsave(file.path(OUT, "figF_tpr_transfer.png"), pF, width = 9, height = 4.4, dpi = 200)

cat("\nfigures figE/figF written\n")
