# WS-8 v3 — Per-feature heterogeneity map underlying M6's native panel.
# Full-atlas M6 fit (all 8 cohorts): per-feature RE-meta rho, BH-adjusted pooled
# p, Cochran's Q p, tau2, I2, homogeneity pass, and per-cohort Spearman rhos.
# Outputs: tables/heterogeneity_m6_fullatlas.csv + figures/fig_m6_heterogeneity_map.png
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table)
})
engine_dir <- Sys.getenv("WS8_ENGINE", v2_dir)
source(file.path(engine_dir, "metrics.R"))
source(file.path(engine_dir, "loco.R"))
source(file.path(engine_dir, "selectors.R"))
source(file.path(v3_dir, "selectors_v3.R"))

# Path config (repo-relative defaults; override via env vars).
v3_dir <- Sys.getenv("WS8_V3_DIR", "benchmarks/ws8_v3")
v2_dir <- Sys.getenv("WS8_V2_DIR", "benchmarks/ws8_v2")
in_dir <- Sys.getenv("WS8_V2_IN", file.path(v2_dir, "results", "m1"))
out <- Sys.getenv("WS8_V3_FIGOUT", file.path(v3_dir, "results"))
dir.create(file.path(out, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out, "figures"), showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(file.path(in_dir, "donor_tables.rds"))
prop <- (obj$ct_mat + 0.5) / (obj$total + 0.5)
X <- log(prop / (1 - prop)); colnames(X) <- paste0("CT:", colnames(X))
meta <- obj$meta[match(rownames(X), obj$meta$donor), ]
y <- setNames(as.numeric(meta$age), meta$donor)
cohort <- setNames(as.character(meta$dataset), meta$donor)
recipe <- fit_fold_transform(X, cohort)
Xf <- apply_fold_transform(recipe, X)

m6 <- sel_m6_meta_homog(Xf, y, cohort)
ex <- m6$extra
rho_c <- m6$detail                       # features x cohorts Spearman rho
stopifnot(identical(rownames(rho_c), names(m6$scores)))

tab <- data.frame(
  feature = names(m6$scores),
  rho_re = as.numeric(m6$scores),
  pooled_padj = as.numeric(ex$pooled_padj[names(m6$scores)]),
  q_p = as.numeric(ex$q_p[names(m6$scores)]),
  tau2 = as.numeric(ex$tau2[names(m6$scores)]),
  I2 = as.numeric(ex$I2[names(m6$scores)]),
  homog_pass = as.logical(ex$homog_pass[names(m6$scores)]),
  native = names(m6$scores) %in% m6$native,
  sign_discordant = sapply(names(m6$scores), function(f)
    any(rho_c[f, ] > 0.05) && any(rho_c[f, ] < -0.05)),
  stringsAsFactors = FALSE)
tab <- cbind(tab, as.data.frame(rho_c[tab$feature, , drop = FALSE]))
tab <- tab[order(-abs(tab$rho_re)), ]
write.csv(tab, file.path(out, "tables", "heterogeneity_m6_fullatlas.csv"),
          row.names = FALSE)

cat(sprintf("features: %d | native (BH<0.05 & Q p>=0.10): %d [%s]\n",
            nrow(tab), sum(tab$native),
            paste(tab$feature[tab$native], collapse = ", ")))
cat(sprintf("BH<0.05: %d | of those, heterogeneous (Q p<0.10): %d\n",
            sum(tab$pooled_padj < 0.05),
            sum(tab$pooled_padj < 0.05 & !tab$homog_pass)))
cat(sprintf("median I2: %.2f | sign-discordant across cohorts: %d features\n",
            median(tab$I2), sum(tab$sign_discordant)))

# ---- Figure -------------------------------------------------------------------
tab$feature <- factor(tab$feature, levels = rev(tab$feature))  # top = largest |rho|
tab$category <- ifelse(tab$native, "native (sig + homogeneous)",
                ifelse(tab$pooled_padj < 0.05, "significant, heterogeneous",
                ifelse(tab$homog_pass, "n.s., homogeneous",
                       "n.s., heterogeneous")))
tab$category <- factor(tab$category,
                       levels = c("native (sig + homogeneous)",
                                  "significant, heterogeneous",
                                  "n.s., homogeneous", "n.s., heterogeneous"))
cat_cols <- c("native (sig + homogeneous)" = "#75A025",
              "significant, heterogeneous" = "#FF9400",
              "n.s., homogeneous" = "grey55", "n.s., heterogeneous" = "grey85")

rho_long <- data.frame(feature = rep(tab$feature, ncol(rho_c)),
                       cohort = rep(colnames(rho_c), each = nrow(tab)),
                       rho = as.numeric(as.matrix(rho_c[levels(tab$feature), ])))
# fix ordering: as.matrix on reordered frame
rho_long <- do.call(rbind, lapply(colnames(rho_c), function(cn)
  data.frame(feature = tab$feature, cohort = cn,
             rho = rho_c[as.character(tab$feature), cn])))

th <- theme_bw(base_size = 10, base_family = "Liberation Sans") +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6.5))

p_heat <- ggplot(rho_long, aes(cohort, feature, fill = rho)) +
  geom_tile() +
  scale_fill_gradient2(low = "#0279EE", mid = "white", high = "#FF9400",
                       midpoint = 0, limits = c(-0.8, 0.8), oob = scales::squish,
                       name = "Spearman rho") +
  labs(x = NULL, y = NULL) +
  th + theme(axis.text.y = element_text(
    size = 6.5,
    face = ifelse(levels(tab$feature) %in% tab$feature[tab$native],
                  "bold", "plain")),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5))

p_cat <- ggplot(tab, aes("status", feature, fill = category)) +
  geom_tile() +
  scale_fill_manual(values = cat_cols, name = NULL, drop = FALSE) +
  labs(x = NULL, y = NULL) +
  th + theme(axis.text = element_blank(), axis.ticks = element_blank(),
             legend.position = "bottom",
             legend.text = element_text(size = 7))

p_i2 <- ggplot(tab, aes(I2, feature)) +
  geom_col(aes(fill = category), width = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = cat_cols) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "I2", y = NULL) +
  th + theme(axis.text.y = element_blank())

p_tau <- ggplot(tab, aes(tau2, feature)) +
  geom_col(aes(fill = category), width = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = cat_cols) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = "tau2 (Fisher-z)", y = NULL) +
  th + theme(axis.text.y = element_blank())

fig <- p_cat + p_heat + p_i2 + p_tau +
  plot_layout(widths = c(0.25, 2.6, 1.0, 1.0), guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(out, "figures", "fig_m6_heterogeneity_map.png"), fig,
       width = 10.5, height = 9.5, dpi = 200)
cat("figure written\n")
