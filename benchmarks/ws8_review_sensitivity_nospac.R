# WS-8 review follow-up 2: binding-cohort disclosure and min-agg sensitivity
# excluding SPAC.
#
# The min-agg score is the minimum per-cohort selection frequency.  This
# script (a) reports which cohort sets that minimum for each of the 24
# primary-selected features, and (b) recomputes the min-agg rule over the
# seven cohorts excluding SPAC (n = 33, the only cohort with fewer donors
# than features).  Reads the objects saved by ws8_aging_showcase.R; no
# refitting.
#
# Usage: run from the repo root, e.g.
#   Rscript benchmarks/ws8_review_sensitivity_nospac.R
#
# Application case study — NOT calibration evidence.

pkgload::load_all(".", compile = FALSE, quiet = TRUE)

obj_path <- Sys.getenv("WS8_OBJECTS", "/workspace/ws8_out/ws8_objects.rds")
out_dir  <- Sys.getenv("WS8_OUT_DIR", "benchmarks")
obj <- readRDS(obj_path)

es  <- obj$env_rf$env_scores   # features x cohorts, max-over-lambda scores
sel <- obj$sel_rf
stopifnot(identical(rownames(es), sel$feature))

# ---- (a) Binding cohort for each of the 24 primary-selected features --------
sel24 <- sel$feature[sel$selected_minagg]
M <- es[sel24, , drop = FALSE]
bind_first <- apply(M, 1, function(r) colnames(es)[which.min(r)])
bind_ties  <- apply(M, 1, function(r) paste(colnames(es)[r == min(r)], collapse = "+"))
binding <- data.frame(
  feature         = sel24,
  min_score       = M[cbind(sel24, bind_first)],
  binding_cohort  = bind_first,
  binding_ties    = bind_ties,
  stringsAsFactors = FALSE
)
write.csv(binding, file.path(out_dir, "ws8_binding_cohorts.csv"), row.names = FALSE)

cat("Binding cohort (first argmin) counts for the 24 primary features:\n")
print(table(binding$binding_cohort))
cat("Features where SPAC is among the minima (ties included):",
    sum(grepl("SPAC", bind_ties)), "of", length(sel24), "\n\n")

# ---- (b) Min-agg sensitivity: exclude SPAC ----------------------------------
es_ns   <- es[, setdiff(colnames(es), "SPAC"), drop = FALSE]
min_ns  <- setNames(apply(es_ns, 1, min), rownames(es))
sens <- data.frame(
  feature             = sel$feature,
  min_score_all8      = sel$min_score,
  selected_all8       = sel$selected_minagg,
  min_score_excl_SPAC = as.numeric(min_ns[sel$feature]),
  selected_excl_SPAC  = as.numeric(min_ns[sel$feature]) >= 0.5,
  stringsAsFactors    = FALSE
)
write.csv(sens, file.path(out_dir, "ws8_minagg_sensitivity_nospac.csv"),
          row.names = FALSE)

cat(sprintf("Selected: all-8 rule %d; excluding SPAC %d (of %d features)\n",
            sum(sens$selected_all8), sum(sens$selected_excl_SPAC), nrow(sens)))
cat("Leave the panel:", paste(setdiff(sel24, sens$feature[sens$selected_excl_SPAC]), collapse = ", "), "\n")
cat("Join the panel:",  paste(setdiff(sens$feature[sens$selected_excl_SPAC], sel24), collapse = ", "), "\n")

paper8 <- c("CT:CD8+ T Naive", "CT:CD8+ RTE", "CT:CD4+ RTE", "CT:CD4+ Temra",
            "CT:MAIT", "CT:CD8+ Temra", "CT:CD8+ Tem GZMK+", "CT:CD8+ Tem GZMB+")
cat("\nPredeclared-8 min scores, all-8 vs excluding SPAC:\n")
print(data.frame(all8 = sel$min_score[match(paper8, sel$feature)],
                 excl_SPAC = round(min_ns[paper8], 3), row.names = paper8))
