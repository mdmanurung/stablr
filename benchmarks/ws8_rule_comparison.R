# WS-8 follow-up: how many min-agg selected features survive the source
# publication's ">= 5/8 cohorts" consistency rule?
#
# Two approximations of the paper's rule are computed for all 60 features:
#   (a) significance vote: per-cohort Spearman correlation with age,
#       BH-adjusted within cohort, counted as consistent when adj. p < 0.05
#       in >= 5 of 8 cohorts (raw p < 0.05 reported as sensitivity);
#   (b) frequency vote: per-cohort STABL selection frequency >= 0.5 in
#       >= 5 of 8 cohorts (the vote relaxation of the min-agg rule).
# The paper's exact per-cohort test is not public in our extract, so (a) is
# an approximation; both directions of disagreement are reported.
#
# Run from the repository root:
#   R_LIBS=<rlib> Rscript benchmarks/ws8_rule_comparison.R

obj <- readRDS(file.path(
  Sys.getenv("WS8_OUT_DIR", "/workspace/ws8_out"), "ws8_objects.rds"
))

X <- obj$X
y <- obj$y_age
env <- obj$env
cohorts <- sort(unique(env))

# Per-cohort Spearman p-values (asymptotic; ties present in age).
p_raw <- matrix(NA_real_, ncol(X), length(cohorts),
                dimnames = list(colnames(X), cohorts))
for (cv in cohorts) {
  idx <- env == cv
  p_raw[, cv] <- vapply(
    seq_len(ncol(X)),
    function(j) suppressWarnings(
      stats::cor.test(X[idx, j], y[idx], method = "spearman",
                      exact = FALSE)$p.value
    ),
    numeric(1)
  )
}
p_bh <- apply(p_raw, 2, p.adjust, method = "BH")

sig_count_bh  <- rowSums(p_bh < 0.05)
sig_count_raw <- rowSums(p_raw < 0.05)
freq_count    <- rowSums(obj$env_rf$env_scores >= 0.5)

out <- data.frame(
  feature          = colnames(X),
  min_score        = obj$sel_rf$min_score[match(colnames(X), obj$sel_rf$feature)],
  selected_minagg  = colnames(X) %in% obj$sel_rf$feature[obj$sel_rf$selected_minagg],
  freq_ge50_count  = as.integer(freq_count),
  sig_count_bh     = as.integer(sig_count_bh),
  sig_count_raw    = as.integer(sig_count_raw),
  pass_5of8_freq   = freq_count >= 5L,
  pass_5of8_sig_bh = sig_count_bh >= 5L,
  stringsAsFactors = FALSE
)
write.csv(out, "benchmarks/ws8_rule_comparison.csv", row.names = FALSE)

sel <- out[out$selected_minagg, ]
cat("== Rule comparison on the 24 min-agg selected features ==\n")
cat(sprintf("survive >=5/8 significance vote (BH):  %d/24\n",
            sum(sel$pass_5of8_sig_bh)))
cat(sprintf("survive >=5/8 frequency vote:          %d/24\n",
            sum(sel$pass_5of8_freq)))
cat("\nMin-agg features FAILING the >=5/8 significance vote:\n")
print(sel[!sel$pass_5of8_sig_bh,
          c("feature", "min_score", "sig_count_bh", "freq_ge50_count")],
      row.names = FALSE)

cat("\n== Set sizes under each rule ==\n")
cat(sprintf("min-agg (8/8 frequency >= 0.5):        %d features\n",
            sum(out$selected_minagg)))
cat(sprintf(">=5/8 frequency vote:                  %d features\n",
            sum(out$pass_5of8_freq)))
cat(sprintf(">=5/8 significance vote (BH):          %d features\n",
            sum(out$pass_5of8_sig_bh)))
cat(sprintf(">=5/8 significance vote (raw p):       %d features\n",
            sum(sig_count_raw >= 5L)))

cat("\n== Features passing >=5/8 significance but NOT min-agg ==\n")
extra <- out[out$pass_5of8_sig_bh & !out$selected_minagg, ]
cat(sprintf("count: %d\n", nrow(extra)))
print(extra[order(-extra$sig_count_bh),
            c("feature", "min_score", "sig_count_bh", "freq_ge50_count")],
      row.names = FALSE)
