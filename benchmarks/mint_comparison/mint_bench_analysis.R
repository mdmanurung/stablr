# Analysis + figures for the stablr vs MINT benchmark
suppressMessages({
  library(ggplot2); library(dplyr); library(matrixStats)
})
themeset <- theme_minimal(base_size = 11, base_family = "Liberation Sans") +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom", legend.box = "vertical")
okabe <- c(stablr_rp = "#0072B2", stablr_ke = "#D55E00", wrap_rp = "#009E73",
           wrap_ke = "#CC79A7", mint = "#E69F00", splsda = "#56B4E9")
meth_lab <- c(stablr_rp = "stablr lasso+RP", stablr_ke = "stablr lasso+KE",
              wrap_rp = "sPLS-DA stabl+RP", wrap_ke = "sPLS-DA stabl+KE",
              mint = "mint.splsda (tuned)", splsda = "sPLS-DA pooled (tuned)")
ck <- "/mnt/shared-workspace/shared/mint_bench"
out <- "/mnt/results/mint_benchmark"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

filesA <- sort(list.files(ck, pattern = "^seedA_", full.names = TRUE))
filesB <- sort(list.files(ck, pattern = "^seedB_", full.names = TRUE))
A <- lapply(filesA, readRDS); B <- lapply(filesB, readRDS)
SIG <- 1:10

# ---- Part A per-seed table --------------------------------------------------
rowsA <- bind_rows(lapply(A, function(r) {
  bind_rows(lapply(c("stablr_rp", "stablr_ke", "wrap_rp", "wrap_ke", "mint", "splsda"),
    function(m) {
      sel <- r[[m]]$sel
      time_field <- if (grepl("^stablr", m)) "stablr_time" else
        if (grepl("^wrap", m)) "wrap_time" else "mint_time"
      data.frame(seed = r$seed, method = m,
                 tpr = mean(SIG %in% sel),
                 fdp = if (length(sel)) mean(!sel %in% SIG) else NA_real_,
                 n_sel = length(sel),
                 min_fdr = if (is.null(r[[m]]$min_fdr)) NA_real_ else r[[m]]$min_fdr,
                 time = r[[time_field]])
    }))}))
rowsA$method <- factor(rowsA$method, levels = names(okabe), labels = meth_lab)

summA <- rowsA %>% group_by(method) %>%
  summarise(tpr = mean(tpr), fdp = mean(fdp, na.rm = TRUE),
            n_sel = mean(n_sel),
            min_fdr = mean(min_fdr, na.rm = TRUE),
            calib_ratio = mean(fdp, na.rm = TRUE) / mean(min_fdr, na.rm = TRUE),
            time_s = mean(time), .groups = "drop")
# stability: mean pairwise Jaccard of selected sets across seeds
stab <- sapply(names(okabe), function(m) {
  sels <- lapply(A, function(r) r[[m]]$sel)
  js <- c()
  for (i in seq_along(sels)) for (j in (i + 1):length(sels)) {
    if (j > length(sels)) break
    js <- c(js, if (length(sels[[i]]) == 0 && length(sels[[j]]) == 0) NA
              else length(intersect(sels[[i]], sels[[j]])) /
                   length(union(sels[[i]], sels[[j]]))) }
  mean(js, na.rm = TRUE)})
summA$jaccard <- stab[as.integer(summA$method)]
print(summA, digits = 2, width = 200)
write.csv(summA, file.path(out, "partA_summary.csv"), row.names = FALSE)
write.csv(rowsA, file.path(out, "partA_per_seed.csv"), row.names = FALSE)

# ---- Part B summary ---------------------------------------------------------
rowsB <- bind_rows(lapply(B, function(r) {
  bind_rows(lapply(1:3, function(k) {
    data.frame(seed = r$seed, cohort = paste0("C", k),
               stablr_rp = r[[paste0("stablr_rp_C", k)]],
               wrap_rp = r[[paste0("wrap_rp_C", k)]],
               mint = r[[paste0("mint_C", k)]],
               splsda = r[[paste0("splsda_C", k)]])}))}))
summB <- rowsB %>% summarise(across(-c(seed, cohort), mean, na.rm = TRUE)) %>%
  tidyr::pivot_longer(everything(), names_to = "method", values_to = "ber_mean")
summB$method <- factor(summB$method, levels = names(okabe), labels = meth_lab)
print(summB, digits = 3)
write.csv(rowsB, file.path(out, "partB_per_seed.csv"), row.names = FALSE)
write.csv(summB, file.path(out, "partB_summary.csv"), row.names = FALSE)

# ---- Figures ----------------------------------------------------------------
pA <- ggplot(rowsA, aes(tpr, fdp, colour = method)) +
  geom_point(aes(shape = method), alpha = 0.55, size = 2, position = position_jitter(width = 0.01, height = 0.01)) +
  geom_point(data = summA, aes(tpr, fdp), size = 4.5, shape = 21,
             fill = "white", stroke = 1.2, show.legend = FALSE) +
  scale_colour_manual(values = setNames(okabe, meth_lab), name = NULL) +
  scale_shape_manual(values = c(16, 17, 15, 18, 8, 4), name = NULL) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "TPR (10 true features)", y = "realized FDP of selected support",
       title = "Selection quality: stablr vs MINT sPLS-DA",
       subtitle = "multi-cohort DGP (3 cohorts, batch shifts + prevalence offsets); large dots = means; 20 seeds") +
  themeset
ggsave(file.path(out, "figA_selection_quality.png"), pA, width = 8, height = 6, dpi = 200)

ctrl <- subset(rowsA, method %in% c("stablr lasso+RP", "stablr lasso+KE",
                                    "sPLS-DA stabl+RP", "sPLS-DA stabl+KE"))
pB <- ggplot(ctrl, aes(min_fdr, fdp, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(shape = method), alpha = 0.6, size = 2) +
  scale_colour_manual(values = setNames(okabe[names(okabe) %in% c("stablr_rp", "stablr_ke", "wrap_rp", "wrap_ke")], meth_lab[1:4]), name = NULL) +
  scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "reported min_fdr_ (FDP+ estimate)", y = "realized FDP",
       title = "Calibration of the FDP+ estimate across learners",
       subtitle = "dashed = perfect calibration; the estimator is learner-agnostic") +
  themeset
ggsave(file.path(out, "figB_calibration.png"), pB, width = 8, height = 6, dpi = 200)

summA_long <- summA %>% select(method, jaccard, n_sel) %>%
  tidyr::pivot_longer(c(jaccard, n_sel), names_to = "metric", values_to = "value")
pC <- ggplot(summA_long, aes(method, value, fill = method)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = setNames(okabe, meth_lab)) +
  labs(x = NULL, y = NULL,
       title = "Selection stability (mean pairwise Jaccard) and set size",
       subtitle = "20 seeds; higher Jaccard = more reproducible support") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8)) +
  themeset
ggsave(file.path(out, "figC_stability.png"), pC, width = 9, height = 5, dpi = 200)

pD <- ggplot(rowsB %>% tidyr::pivot_longer(-c(seed, cohort), names_to = "method", values_to = "ber") %>%
               mutate(method = factor(method, levels = names(okabe), labels = meth_lab)),
             aes(method, ber, fill = method)) +
  geom_boxplot(width = 0.55, alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5, colour = "grey30") +
  geom_hline(yintercept = 0.5, linetype = 3, colour = "grey40") +
  scale_fill_manual(values = setNames(okabe, meth_lab)) +
  labs(x = NULL, y = "balanced error rate (held-out cohort)",
       title = "LOCO prediction: train on 2 cohorts, predict the third",
       subtitle = "10 seeds x 3 folds; dotted = chance") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  themeset
ggsave(file.path(out, "figD_loco_ber.png"), pD, width = 8, height = 5.5, dpi = 200)

cat("analysis complete\n")
