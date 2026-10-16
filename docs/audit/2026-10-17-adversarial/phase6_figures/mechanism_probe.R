# Probe: why does the RP null curve hug 1.0 while RP signal curves dip deepest?
# Quantifies (a) real-vs-artificial count curves per regime, (b) what the argmin
# threshold actually lands on (the "notch" lottery), (c) winner's-curse depth.
suppressMessages(library(dplyr))
curves <- readRDS("/workspace/phase6_out/curves.rds")
thr <- seq(0, 0.99, by = 0.01)

cnt <- bind_rows(lapply(curves, function(r) {
  nr <- colSums(outer(r$scores_real, thr, ">"))  # strict exceedance, matches findInterval usage
  na <- colSums(outer(r$scores_art,  thr, ">"))
  data.frame(design = r$design, art = r$artificial, seed = r$seed,
             threshold = thr, nr = nr, na = na, fdp = pmin(r$FDRs, 2.5),
             min_fdr = r$min_fdr, thr_star = r$thr_star)
}))

cat("=== median counts n_real / n_art and median FDP+ at fixed thresholds ===\n")
print(as.data.frame(
  cnt %>% filter(threshold %in% c(0.3, 0.5, 0.7, 0.9)) %>%
    group_by(design, art, threshold) %>%
    summarise(nr = median(nr), na = median(na), fdp_med = round(median(fdp), 2), .groups = "drop")
))

cat("\n=== where the argmin lands (per-seed notch, averaged over seeds) ===\n")
notch <- bind_rows(lapply(curves, function(r) {
  i <- which.min(r$FDRs)  # first-index-wins, matches the fit's argmin
  data.frame(design = r$design, art = r$artificial, seed = r$seed,
             thr = thr[i], fdp = r$FDRs[i],
             nr = sum(r$scores_real > thr[i]),
             na = sum(r$scores_art > thr[i]))
}))
print(as.data.frame(
  notch %>% group_by(design, art) %>%
    summarise(thr = round(mean(thr), 2), fdp = round(mean(fdp), 3),
              nr = round(mean(nr), 1), na = round(mean(na), 1), .groups = "drop")
))

cat("\n=== winner's curse: mean reported min_fdr_ vs mean median-of-curve FDP+ ===\n")
print(as.data.frame(
  cnt %>% group_by(design, art, seed) %>%
    summarise(min_fdr = first(min_fdr), med_fdp = median(fdp), .groups = "drop") %>%
    group_by(design, art) %>%
    summarise(reported = round(mean(min_fdr), 3),
              curve_median = round(mean(med_fdp), 3), .groups = "drop")
))

cat("\n=== mean max-over-lambda score, real vs artificial (exchangeability check) ===\n")
print(as.data.frame(
  bind_rows(lapply(curves, function(r)
    data.frame(design = r$design, art = r$artificial,
               msr = mean(r$scores_real), msa = mean(r$scores_art)))) %>%
    group_by(design, art) %>%
    summarise(mean_real = round(mean(msr), 4), mean_art = round(mean(msa), 4),
              diff = round(mean(msr - msa), 4), .groups = "drop")
))
