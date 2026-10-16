# stablr (commit 5cf603b) — Phase 1 Adversarial Audit: Specification-to-Implementation Traceability Matrix

Auditor: independent statistical methodologist (Phase 1 of adversarial audit).
Repo: `/workspace/integrate/stablr` (read-only; no repo file modified).
Scope: reconstruct the mathematical specification of each core algorithm from the primary
literature, then trace it into the implementation with file:line evidence, classify each
component, and rank the most scientifically consequential gaps.

---

## (i) Reference list (primary sources, with DOIs)

| # | Reference | DOI / identifier | Role in audit |
|---|-----------|------------------|---------------|
| R1 | Hédou J, Marić I, Bellan G, Einhaus J, Gaudillière DK, Ladant F-X, Verdonk F, Stelzer IA. **"Discovery of sparse, reliable omic biomarkers with Stabl."** *Nature Biotechnology* 42:208–219 (2024). | 10.1038/s41587-023-02033-x | Stabl algorithm, FDP+ estimator, reliability threshold θ, Theoretical guarantee (Lemma 1, Theorem 1, Assumption 1) |
| R2 | Meinshausen N, Bühlmann P. **"Stability selection."** *JRSS-B* 72(4):417–473 (2010). | 10.1111/j.1467-9868.2010.00740.x (arXiv:0809.2932) | Subsampling scheme, stability path, PFER bound E(V) ≤ q²/(p(2π_thr−1)) |
| R3 | Ding DY, Li S, Narasimhan B, Tibshirani RJ. **"Cooperative learning for multiview analysis."** *PNAS* 119(38):e2202113119 (2022). | 10.1073/pnas.2202113119 | Cooperative-learning objective, agreement penalty, block coordinate descent |
| R4 | Barber RF, Candès E. **"Controlling the false discovery rate via knockoffs."** *Annals of Statistics* 43(5):2055–2085 (2015). | 10.1214/15-AOS1337 | Fixed-X knockoffs, knockoff+ FDP estimator, stopping-time structure of Eq. (3)-style thresholds |
| R5 | Candès E, Fan Y, Janson L, Lv J. **"A note on the knockoff filter and the knockoff+ construction."** *Annals of Statistics* 46(3):940–968 (2018) ("Panning for gold: model-X knockoffs"). | 10.1214/17-AOS1588 | Model-X knockoff construction requirements: exchangeability (pairwise + joint with response), sufficiency, antiseptic/constructive |
| R6 | Spector P, Janson L. **"Powerful knockoffs via minimizing reconstructability."** *Annals of Statistics* 50(4):2022–2047 (2022); arXiv:2011.14625. | 10.48550/arXiv.2011.14625 | MVR S-matrix loss L(S) = tr((2Σ−S)⁻¹) + tr(S⁻¹) |
| R7 | Shah RD, Samworth RJ. **"Variable selection with error control: another look at stability selection."** *JRSS-B* 75(1):55–80 (2013). | 10.1111/rssb.12015 | Complementary-pairs variant; assumption-free bounds; critique of exchangeability |
| R8 | Peters J, Bühlmann P, Meinshausen N. **"Causal inference by using invariant prediction."** *JRSS-B* 78(5):947–1012 (2016). | 10.1111/rssb.12167 | ICP — conceptual motivation for environment stability (no causal claims in stablr) |
| R9 | Reference implementation: Bellan G, Hédou J et al., **Python `Stabl` library**, https://github.com/gregbellan/Stabl (`stabl/stabl.py`). | — | Cross-language parity target for stablr |
| R10 | Emami H, Emami M, et al. / `knockpy` package (Spector & Janson), https://github.com/amspector100/knockpy (`mrc.py`). | — | Source of the vendored MVR solver logic |

---

## (ii) Formal specifications (what the algorithms SHOULD compute)

### (a) Stabl bootstrap stability selection (R1, Methods; R2)

Data (Y, X) ∈ ℝⁿ × ℝⁿˣᵖ, p ≫ n regime. Let an SRM (sparse regularized model) be e.g. the
lasso β̂_Lasso(λ) = argmin_b {‖Y − Xb‖₂² + λ‖b‖₁}.

1. **Augmentation.** Generate p artificial features X̃ = (X̃₁,…,X̃_p) ∈ ℝⁿˣᵖ, uninformative by
   construction (MX knockoffs or random permutations), and set 𝕏 = [X | X̃] ∈ ℝⁿˣ²ᵖ.
   O = {1,…,p} (original), A = {p+1,…,2p} (artificial).
2. **Resampling.** B iterations k = 1,…,B; at each, a subsample (Y,X)ₖ of size ⌊n/2⌋
   (the SS scheme of R2, Definition 1: drawn **without replacement**; M&B choose ⌊n/2⌋
   because it "resembles most closely the bootstrap" while remaining a U-statistic scheme).
3. **Path traversal.** Fit the SRM on each subsample over a grid λ ∈ Λ ⊆ ℝ≥₀, giving
   β̂(k,λ) = (β̂₁(k,λ),…,β̂₂ₚ(k,λ))ᵀ. For models with two hyperparameters (EN, SGL) the grid
   is 2-D and the max below is taken over the full hyperparameter set.
4. **Selection frequency (stability path).**
   f_i(λ) = (1/B) Σ_{k=1}^B 1[β̂_i(k,λ) ≠ 0];  the "Stabl path"/stability graph plots f_i(λ)
   against 1/λ. The path-level score is the **max over the grid**:
   f_j = max_{λ∈Λ} f_j(λ).
5. **Selection at threshold t.** Ŝ(t) = {j ∈ O : f_j ≥ t} (paper uses ≥; see Gap G2).
6. **Final model.** β̂_Stabl = argmin_b ‖Y − Xb‖₂² s.t. b_i = 0 for i ∉ Ŝ(θ) — a refit on the
   selected support on the **original** (non-augmented) X.

Key assumptions: exchangeability (Assumption 1, below); the SRM's selection event is
well-defined per (k,λ); B large enough that Monte Carlo noise in f_j is small (paper:
>500 artificial features improves FDP+ estimation; Extended Data Fig. 3).

### (b) Artificial-feature proportional FDP control — the FDP+ estimator (R1)

With f_j as above and threshold t ∈ [0,1]:

**FDP+(t) = (1 + Σ_{j∈A} 1[f_j ≥ t]) / (Σ_{j∈O} 1[f_j ≥ t] ∨ 1)**   (R1, Methods, Eq. before (1))

**Reliability threshold:** θ ∈ argmin_{t∈[0,1]} FDP+(t)  (R1, Eq. 1), and the achieved
surrogate q+ = FDP+(θ) = min_{t∈[0,1]} FDP+(t)  (R1, Eq. 2).

**Alpha semantics — the crucial contrast (R1, "Theoretical guarantee" section).** Prior
knockoff-style procedures fix a target FDR q ∈ (0,1) and set
t̂ = min{t ∈ (0,1) : FDP+(t) ≤ q}  (R1, Eq. 3),
which is a **stopping time** w.r.t. a suitable filtration, yielding
FDR_t̂ = E{FDP+(t̂)} ≤ q  (R1, Eq. 4).
Stabl instead **minimizes** FDP+ over t; the paper explicitly states θ "is not a stopping
time, and, therefore, a similarly simple argument is not available" and that q+ "is itself
random." Consequently:

- There is **no user-facing α/q** in Stabl: the operating level is the random q+ = min_t FDP+(t).
- The guarantee is only **relative and stochastic** (R1, Lemma 1 / Theorem 1):
  under Assumption 1, for any Δ > 0,
  **P(FDP(θ) ≥ (1+Δ)·q+) ≤ 1/(1+Δ)**  (R1, Eq. 7);
  and if |S| ≤ p/2, with M = |N∩Ŝ(θ)| + |A∩Ŝ(θ)| false discoveries,
  **P(FDP(θ) ≥ (1+Δ)q+ and M ≥ m) ≤ 2e^{−c* m Δ²}**  (R1, Eq. 8),
  i.e. FDP(θ) ≤ (1+O_P(m^{−1/2}))·q+  (R1, Eq. 9).
- **Assumption 1 (exchangeability of the extended null set).** For any permutation π of
  N ∪ A: (Y, X_S, X_{N∪A}^π) ≐ (Y, X_S, X_{N∪A})  (R1, Eq. 6). The proof (R1 appendix)
  builds Z_k = n_k/(a_k+1) as a martingale w.r.t. F_k = σ({n_i,a_i : i≤k} ∪ {f_j : j∈S}).

So the correct reading: FDP+ is a *surrogate* whose **minimum** is taken ex ante; the
resulting q+ is an estimated upper-bound surrogate, not a controlled FDR level. Any
documentation or accessor that invites reading `min_fdr_` as "FDR ≤ min_fdr_" overstates
the theory.

### (c) Model-X knockoff construction requirements (R4, R5)

A valid model-X knockoff matrix X̃ ∈ ℝⁿˣᵖ must satisfy, **conditional on X** (and hence for
every realization of X):

1. **Pairwise exchangeability (sufficiency):** X̃_j ⫫ X_{−j}, X̃_{−j} | X  and
   (X, X̃)_{swap(j)} ≐ (X, X̃) for each j — i.e. (X_j, X̃_j) can be swapped conditional on
   all other columns. For Gaussian designs this reduces to the Gram conditions
   X̃ᵀX̃ = XᵀX and XᵀX̃ = XᵀX − diag(s), s ∈ [0, 2·eig(Σ)]ᵀ (R4 for fixed-X; R5 for model-X).
2. **Joint exchangeability with the response (model-X):** (X, X̃, Y) ≐ (X_{swap(j)}, X̃, Y)
   — knockoff validity is a property of the **conditional law of X given everything else**,
   not of the regression; Y may be adversarial.
3. **Antiseptic/constructive requirement:** knockoffs must be constructed **without looking
   at Y** (model-X) or, for fixed-X, must satisfy the Gram condition on the realized design;
   and each knockoff column must remain paired with its own original column (the "swap"
   is per-column). Sub-sampling rows of a valid construction preserves model-X validity
   (rows i.i.d. from the fitted model), but **fixed-X** validity (Gram condition on the
   realized matrix) does **not** survive arbitrary row deletion.
4. **Gaussian second-order knockoffs** (what stablr uses): approximate X ~ N(μ, Σ) with
   μ̂ = colMeans, Σ̂ = sample covariance (optionally shrunk), then draw X̃ | X from the
   conditional Gaussian with the equi (S = γΣ̂) or MVR (S = argmin tr((2Σ−S)⁻¹)+tr(S⁻¹))
   choice. This is exact model-X **only when X is truly Gaussian**; otherwise it is a
   second-moment approximation, and Assumption 1 of R1 inherits that approximation.

### (d) Cooperative learning objective (R3)

For M views X^(1),…,X^(M) (standardized columns, centered y), fixed ρ ≥ 0, common λ:

min over β^(1),…,β^(M):
  (1/2)‖y − Σ_m X^(m)β^(m)‖² + λ Σ_m [ α‖β^(m)‖₁ + (1−α)/2 ‖β^(m)‖₂² ]
  + ρ Σ_{m<m′} w_{mm′} ‖X^(m)β^(m) − X^(m′)β^(m′)‖₂²   (R3, Eqs. 3–5, 9, 22)

- ρ = 0 ⇒ early fusion (concatenated elastic net); ρ = 1 with ℓ₂-only ⇒ late fusion
  (average of marginal fits) (R3, "Relation to Early/Late Fusion").
- **Direct algorithm:** stack augmented rows — for each pair (m,m′) add rows
  √(ρw_{mm′})·[X^(m), −X^(m′)] — turning Eq. (22) into a single lasso/EN on 2n·C(M,2) + n
  pseudo-observations (R3, Eqs. 6–7, 23–24).
- **One-at-a-time algorithm (block coordinate descent):** update view m by fitting an
  elastic net on the penalty-adjusted partial residual
  (y − Σ_{m′≠m} X^(m′)β^(m′)) with the agreement term folded in (R3, Eqs. 10, 21, 25);
  for GLMs, an outer IRLS loop with an inner WLS-EN solve (R3, Eq. 13 discussion).
- Weights w_{mm′}: paper uses uniform in experiments; view-pair-specific weights are
  allowed (R3, after Eq. 21: "We can also have different weights on the agreement penalties
  for distinct pairs of data views").
- CV folds are formed on rows of X, Z (not of the augmented matrix) (R3, Remark after Alg. 1).

### (e) Stability selection theory (R2)

- **Scheme:** I ⊂ {1,…,n} random subsample of size ⌊n/2⌋, **without replacement** (R2, Def. 1).
- **Selection probability / stability path:** Π̂_k(λ) = P*(k ∈ Ŝ_λ(I)); stable set
  Ŝ^stable = {k : max_{λ∈Λ} Π̂_k(λ) ≥ π_thr}  (R2, Defs. 1–2, Eq. 7).
- **Theorem 1 (error control).** Assume (i) exchangeability of {1{k∈Ŝ_λ}, k∈N} for all λ,
  and (ii) "not worse than random guessing": E(|S∩Ŝ_λ|)/E(|N∩Ŝ_λ|) ≥ |S|/|N|  (R2, Eq. 8).
  Then for π_thr ∈ (1/2, 1):
  **E(V) ≤ (1/(2π_thr − 1)) · q_Λ²/p**,  q_Λ = E(|Ŝ_Λ(I)|)  (R2, Eq. 9).
- **Error-rate semantics:** the bound is on the **per-family error rate (PFER)** — the
  *expected number* of false selections — not on FDR/FDP; E(V)/p is the per-comparison
  error rate (R2, after Eq. 9). This is a *per-family* guarantee; it does not give a
  proportional (FDR) statement, which is precisely the gap Stabl's FDP+ machinery (and its
  own Lemma 1/Theorem 1) is designed to fill.
- R7 (CPSS) removes/weakens exchangeability via complementary pairs; neither Stabl nor
  stablr uses complementary pairs.

---

## (iii) Traceability matrix (specification → implementation)

Line numbers refer to the working tree at commit 5cf603b (`/workspace/integrate/stablr`).
Classifications: **FT** faithful translation · **NE** numerically equivalent alternative ·
**IE** intentional (documented) extension · **UE** unverified extension ·
**PD** potential implementation defect · **SCD** scientifically consequential divergence.

| # | Component | Formal estimand / objective (spec) | Key assumptions | Implementation location (file:line) | Classification | Independent validation procedure |
|---|-----------|-----------------------------------|-----------------|-------------------------------------|----------------|--------------------------------|
| 1 | Subsampling scheme | B subsamples of size ⌊n/2⌋ without replacement (R1 §Methods; R2 Def. 1) | i.i.d. rows; ⌊n/2⌋ for U-statistic interpretation | `R/stabl_fit.R:402` (`n_subsamples <- as.integer(floor(sample_fraction * n_samples))`), defaults `sample_fraction=0.5`, `replace=FALSE` (`R/stabl_fit.R:205-206`); index draws `R/bootstrap_helpers.R:63-135` (`sample_bootstrap_indices`) | **FT** (matches R1/R2 and Python defaults `sample_fraction=0.5, replace=False`, R9) | Unit test: empirical inclusion probabilities of each unit across B draws ≈ n_subsamples/n; verify `floor()` semantics and no-replacement; cross-run Python STABL with same seed-family and compare index distributions (KS test) |
| 2 | Selection event per (k,λ) | 1[β̂_j(k,λ) ≠ 0] (R1) | Exact zeros meaningful | `R/stabl_fit.R:679` (`.BOOTSTRAP_COEF_THRESHOLD <- 1e-5`); masks `abs(coef) > 1e-5` at `R/learner_adapters.R:515,532`, `R/learner_ordinal.R:253` | **NE** (numerical tolerance for "≠0"; identical to Python's l1 `bootstrap_threshold=1e-5`, R9) | Sensitivity analysis: recompute stability paths at thresholds {0, 1e-8, 1e-5, 1e-3}; confirm support and FDP+ invariant for well-conditioned problems |
| 3 | Stability path & score | f_j(λ) = B⁻¹Σ_k 1[·]; f_j = max_λ f_j(λ) (R1) | Grid Λ covers informative λ range | Accumulation `R/stabl_fit.R:560-608` (`.run_bootstrap_accumulation`: sum masks, divide by `n_bootstraps`); max-over-λ in `R/fdp_control.R:117-118` (`rowMaxs`) and `R/stabl_accessors.R:58` (`get_importances`) | **FT** (identical to Python `np.vstack(...).mean(axis=0)` + `np.max(axis=1)`, R9) | Equivalence test vs. literal per-(k,λ) loop on small simulated data; verify division by total B (see row 4) |
| 4 | Denominator under failed fits | f_j(λ) = (1/B)Σ over **all** B draws (R1 silent on failures) | All B fits succeed | Errors from glmnet are converted to `stablr_numerical_infeasibility` and **re-thrown** (`R/conditions.R:83-100`, `.with_glmnet_numerical_infeasibility`); no tryCatch in `.run_bootstrap_accumulation` (`R/stabl_fit.R:560-608`) → whole `stabl_fit()` aborts (fail-closed). Classification bootstraps: infeasible class-composition draws are **redrawn** (≤1000) then fall back to stratified sampling (`R/stabl_fit.R:558-660`, `.classification_bootstrap_sampler`) | **FT** for failures (fail-closed; denominators never silently altered — no divergence from paper because paper assumes all fits succeed); **IE** for the classification redraw (conditional sampling, documented) | (a) Inject a deterministic glmnet failure (e.g., constant y) and assert `stabl_fit` errors rather than returning scores with denominator < B. (b) For the classifier path, empirically compare selection frequencies under redraw-on-infeasible vs. rejection-free sampling to bound the conditioning effect |
| 5 | FDP+ estimator | FDP+(t) = (1 + Σ_{j∈A}1[f_j≥t]) / (Σ_{j∈O}1[f_j≥t] ∨ 1) (R1) | |A| = p (paper); exchangeability of N∪A | `R/fdp_control.R:119-121`: `FDRs <- (inv_prop * n_art + 1.0) / pmax(1.0, n_real)` with `inv_prop = 1/artificial_proportion`; strict-`>` counting `R/fdp_control.R:131-149` | **FT** at π=1; **NE/IE** generalization: (1/π)·n_art is the standard p/m scaling of the knockoff FDP estimator when m ≠ p (R4/R5); identical to Python `_compute_FDRc` (R9) which also carries the 1/π factor | Exact-match test vs. Python `_compute_FDRc` on identical score matrices; brute-force recount with `sum(max_scores > t)`; verify π=1 reduces to the paper's formula symbol-by-symbol |
| 6 | Threshold grid & argmin rule | θ ∈ argmin_{t∈[0,1]} FDP+(t); q+ = min_t FDP+(t) (R1 Eqs. 1–2) | Continuous t in theory; unique minimizer in practice | `R/fdp_control.R:122-126`: grid `seq(0, 0.99, by=0.01)`; `min_fdr <- min(FDRs)`; `fdr_min_threshold <- min(fdr_threshold_range[which.min(FDRs)], 1.0)` with first-index tie-breaking; cap at 1.0 when min_fdr > 1 (`R/fdp_control.R:126-129`) — all mirroring Python (R9) | **FT** (Python-parity); grid discretization of a continuous argmin is a benign, unanalyzed approximation (see Gap G3) | Dense-grid study: recompute θ on grid 0.001 and compare selected sets and q+; verify the paper's Lemma-1 bound P(FDP(θ) ≥ (1+Δ)q+) ≤ 1/(1+Δ) empirically in simulation with known S |
| 7 | Alpha / FDR-target semantics | No user α; q+ random; only relative stochastic bound (R1, Eq. 7) | Users must not read q+ as controlled FDR | No α/q parameter exists anywhere in `stabl_fit`/`compute_fdp_plus`; `min_fdr_` returned as diagnostic (`R/stabl_fit.R:497-500`); docstrings state "The fitting rule chooses the FDP+ minimizer; it is not a universal false-discovery guarantee" (`R/fdp_control.R:1-9`, `R/artificial_features.R` dispatcher docs) | **FT** of the paper's semantics (correctly *not* implementing Eq. (3)-style fixed-q thresholding); residual risk is interpretive (Gap G1) | Documentation audit + user-facing test: assert no code path treats `min_fdr_` as a nominal level; simulation: empirical FDR of Ŝ(θ) vs. distribution of q+ across replicates (should satisfy Eq. 7, not FDR ≤ q+) |
| 8 | Support rule at θ | Ŝ(θ) = {j ∈ O : f_j ≥ θ} (R1, "≥") | Ties measure-zero for continuous scores | `R/stabl_accessors.R:58-59`: `mask <- max_scores > threshold` (strict `>`), matching Python `_get_support_mask` (`max_scores > final_cutoff`, R9) | **NE** (strict `>` vs paper's `≥`); with B=1000 and a 0.01 grid, exact ties (e.g., f_j = 0.50, θ = 0.50) are possible and are **excluded** by both implementations, contrary to the paper's definition | Targeted test with B chosen so scores hit grid values exactly (e.g., B=100, θ=0.5, f_j=0.5): paper counts selected, stablr/Python do not; quantify impact on FDP+ numerator/denominator |
| 9 | Artificial features — random permutation | Column copies, independently shuffled (R1: "random permutation of the original features") | Exchangeability of permuted nulls with real nulls (only plausible under weak dependence) | `R/artificial_features.R:16-37` (`make_rp_features`: `sample.int` column choice, per-column `sample()` shuffle); generated **once** before the bootstrap loop (`R/stabl_fit.R:351-360`) | **FT** (Python-parity, R9) | Distributional test: marginal quantiles of permuted vs. source columns; joint-independence check (permuted columns must be asymptotically uncorrelated with Y and with each other's ordering) |
| 10 | Artificial features — model-X equi knockoffs | Second-order Gaussian MX knockoffs, S = γΣ̂ (R5; R1 uses `GaussianSampler(method='equicorrelated')`) | Gaussian approximation to X; Σ̂ ≈ Σ | `R/artificial_features.R:427-514` (`make_knockoff_equi_features`): μ = colMeans (l.438), Σ̂ = cov with corpcor `make.positive.definite` fallback (`R/artificial_features.R:395-408`), `knockoff::create.gaussian(x, mu, Sigma, method="equi")` (l.446); generated once, reused across all B bootstraps | **FT** of Python parity; **IE** as second-order approximation — exact model-X only under Gaussianity (Gap G4) | (a) Verify swap exchangeability empirically: marginal/conditional moments of (X_j, X̃_j); (b) simulation with Gaussian X: FDP+ calibration of the full Stabl pipeline; (c) non-Gaussian margins (NB, zero-inflated): quantify FDP+ miscalibration |
| 11 | Artificial features — MVR knockoffs | S = argmin tr((2Σ−S)⁻¹) + tr(S⁻¹) s.t. PSD constraints (R6) | Same Gaussian approximation | `R/knockoff_mvr.R:1-233` (`solve_mvr` coordinate descent, `.mvr_loss` l.60-72, quadratic solve l.88-130, PSD guards l.137-176); sampling via `knockoff::create.gaussian(..., diag_s = S)` (`R/artificial_features.R:598-612`); equi fallback on solver failure (l.558-596) | **IE** (novel vs. Python STABL and vs. R1 — documented as such); **NE** to knockpy's solver (port, R6/R10) | (a) Numerical: compare S and loss against `knockpy.mrc` on identical Σ; (b) statistical: FDP+ calibration of Stabl+MVR vs. Stabl+equi in simulations with correlated designs (MVR should have ≥ power at equal FDP) |
| 12 | Fixed-X knockoff mode | Gram condition X̃ᵀX̃ = XᵀX, XᵀX̃ = XᵀX − diag(s) on the realized design (R4) | Validity tied to the realized X matrix | `R/artificial_features.R` (`make_knockoff_features`, `knockoff::create.fixed(x, sigma=1)` with rescale-to-input-scale contract, documented at the function's Rd); explicit caveat that fixed-X exchangeability "does not survive the row subsampling of the bootstrap loop" | **IE** with documented invalidity under the bootstrap (the package itself recommends model-X modes); kept for backward compatibility | Verify the Gram identities numerically post-rescale; demonstrate FDP+ miscalibration under `sample_fraction<1` to confirm the caveat is not merely theoretical |
| 13 | Knockoff↔feature pairing & subsetting | Swap is per-column; knockoff j pairs with original j (R5) | Positional alignment of score rows | `R/artificial_features.R:494-506` (aligned, sorted `sel_idx`; comment explains positional accumulation) | **FT** | Test: `noise_col_indices` maps each artificial column to its source; permuted-label regression check that artificial score rows align with the intended sources |
| 14 | Dimensional guard & fallback policy | — (not in R1) | FDP+ calibration target must not silently change | `R/artificial_features.R:268-278` (`.assert_supported_knockoff_dimension`, 3000-feature cap, fail-closed); `knockoff_fallback = "stop"` default (`R/stabl_fit.R:224`); strict-fallback abort `R/artificial_features.R:381-391` | **IE** (defensive; prevents a scientifically consequential silent fallback to permutation nulls, which do not calibrate FDP+ under dependence) | Audit trail test: force a knockoff failure and assert error under `"stop"`, warning + provenance row under `"allow"` |
| 15 | Final refit on selected support | β̂_Stabl = argmin ‖Y−Xb‖² s.t. b_i = 0 off Ŝ(θ) (R1, final equation) | Refit on original X (no artificial columns) | Not in `stabl_fit` (returns scores/support only); explicit, labelled full-development refit offered by the analysis workflow (`R/analysis_workflow.R:18-19,66-86`, `final_refit = TRUE`) | **IE** (refit is user-invoked and labelled; consistent with the paper's separation of selection and refit) | End-to-end test: run workflow with `final_refit=TRUE` and confirm the design matrix excludes artificial columns and rows are the full development set |
| 16 | Cooperative learning engine | Objective of R3 (Eq. 22) with agreement penalty ρΣ_{m<m′} w_{mm′}‖X^(m)β^(m) − X^(m′)β^(m′)‖² | Standardized views; gaussian/binomial only in stablr | Vendored fork of CRAN multiview 1.0: `.mv_multiview` (`R/cooperative-multiview.R:208`), pair-row construction `rows <- lapply(pairs, make_row, ...)` (`R/cooperative-multiview.path.R:648`), objective evaluation with `view_components, rho` (l.595, 694), inner EN solve `elnet.fit` (l.677) and the glmnetpp WLS point solver `src/wls_exp.cpp:9-77` (`ElnetPoint`, gaussian WLS) / `spwls_exp.cpp` for sparse; divergence log `inst/PATCHES.cooperative.md` (D1 renames, D2 family restriction — both "No" statistical change under stablr usage); wrappers `R/cooperative-backend.R:26-63` (unit weights, no offsets, ρ passthrough) | **FT** of R3's block-coordinate/augmented-row scheme for gaussian/binomial (restricted family is an **IE**); ρ=0 reduces to early fusion as in R3 | (a) Objective-value test: at ρ>0, evaluate R3 Eq. (22) at the returned solution and compare against the augmented-matrix direct solve (R3 Eqs. 23–24) on small problems — must match to tolerance; (b) ρ=0 vs. glmnet on concatenated views: bit-level agreement; (c) verify CV fold construction on original rows (R3 remark) |
| 17 | Environment-aware stability (min/mean aggregation) | **No literature estimand.** Per-environment Stabl fits; per-feature score = min_e (max_λ f^{(e)}_j(λ)) (and mean) | Environments are labels; no causal claim; no FDP control at the aggregated level | `R/environment_stability.R:107-160` (`stabl_environment_stability`: per-level `stabl_fit` with seeds `random_state + i`, l.117-124); aggregation `min_scores <- apply(env_scores, 1, min)`, `mean_scores <- rowMeans(...)` (l.155-157); ICP motivation and explicit "no causal claims" disclaimer (l.1-13, Rd §Details); **no quantile aggregation exists in the package** (grep across `R/` returns nothing) | **UE** — novel, unvalidated extension. The min-over-environments statistic has no FDP/FDR theory; per-environment FDP+ thresholds are computed but the aggregated score is never thresholded through FDP+; min-aggregation of max-over-λ scores also compounds the leniency of max-over-λ (acknowledged in the Rd) | (a) Simulation with environment-specific nulls: empirical FDR of "min_score ≥ θ_env" rules vs. nominal; (b) permutation test: permute environment labels and confirm min-scores of true nulls are stochastically ≤ those of signals; (c) compare against a pooled-data Stabl fit with environment covariates to quantify power/conservativeness trade-off |
| 18 | Learner extensions (rf, oblique_rf, ordinal, SGL, AL) | R1 validates Stabl with Lasso/EN/AL/SGL; RF-family and ordinal are not in R1 | FDP+ calibration is learner-dependent | Adapters `R/learner_adapters.R` (glmnet/AL/SGL), `R/learner_rf.R`, `R/learner_oblique_rf.R`, `R/learner_ordinal.R`; guardrail warnings that RF-family FDP+ is "empirical, not guaranteed" and must pair with model-X knockoffs (`R/stabl_fit.R:362-377`, `.warn_rf_random_permutation_nulls` / `.warn_rf_fixedx_knockoff_nulls`); `stabl()` beginner default is `knockoff_equi` (`R/stabl.R:20-27`) | **IE/UE** — documented extensions; the RF-family pairing is explicitly flagged as empirically validated only | Simulation-based FDP+ calibration curves per learner × artificial-type cell (the package's own RF validation study protocol); require empirical FDR ≤ q+ + tolerance before promoting from UE |
| 19 | n_injected rounding & effective proportion | Paper: p artificial features (m = p) | π = 1 in paper | `R/stabl_fit.R:337-346`: `n_injected = min(p, max(1, round(p·π)))`, `effective_artificial_proportion = n_injected/p` used in FDP+ (l.486); Python uses `int(p·π)` (truncation) and the **nominal** π in `_compute_FDRc` (R9) | **NE** (rounding vs. truncation; effective vs. nominal π) — identical at π = 1, sub-1% count differences otherwise | Cross-language test at π ∈ {0.3, 0.5, 0.7}: compare n_injected and FDP+ curves |
| 20 | PFER vs FDR framing of the underlying SS theory | R2 Theorem 1 bounds **E(V)** (PFER), not FDR | π_thr ∈ (1/2,1); exchangeability; not-worse-than-random-guessing | stablr makes no PFER claim; FDP+ path is the only error-rate surface (`R/fdp_control.R`); docstrings consistently say "not a universal false-discovery guarantee" | **FT** (no overclaim found in code/docs) | Documentation review; confirm no accessor or vignette quotes the M&B bound as an FDR bound |

---

## (iv) Ranked list of the most scientifically consequential gaps / unverified extensions

**G1 (rank 1) — Interpretive risk around `min_fdr_` / q+ ("alpha semantics").**
Stabl has no user-facing FDR target; θ = argmin_t FDP+(t) and q+ = min_t FDP+(t) is a
**random** quantity whose relationship to the true FDP is only
P(FDP(θ) ≥ (1+Δ)q+) ≤ 1/(1+Δ) (R1, Eq. 7). stablr implements this faithfully and its
docstrings disclaim universal control, but the returned `min_fdr_` invites a
"controlled-FDR = min_fdr_" misreading in downstream tooling. *Consequence:* users may
report "FDR ≈ 0.1" when the theory supports only a stochastic relative bound.
*Validation:* simulation study reporting the joint distribution of (q+, realized FDP)
and empirical coverage of Eq. (7).

**G2 (rank 2) — Environment-aware min/mean aggregation is an unvalidated extension with
no error control.** `stabl_environment_stability` (R/environment_stability.R:107-160)
aggregates max-over-λ scores by min/mean across per-environment fits. Nothing in R1–R8
provides an estimand or error guarantee for this statistic; the per-environment FDP+
thresholds are not propagated to the aggregated score; and the min statistic is
anti-conservative under the documented leniency of max-over-λ scores. The task brief's
"min/quantile aggregation" is not implemented — only min/mean exist. *Consequence:*
"transportable biomarkers" claims rest on an uncalibrated statistic. *Validation:*
label-permutation and environment-shift simulations (row 17).

**G3 (rank 3) — Threshold-minimization optimism is covered by the paper's theory only
under Assumption 1, and the grid discretization is unanalyzed.** The implementation
reproduces the paper's argmin rule (R/fdp_control.R:122-129), and the paper itself
acknowledges θ is not a stopping time (unlike Eq. (3)-style thresholds, R4). The
supermartingale bounds (R1 Lemma 1/Theorem 1) require exchangeability of N∪A; the
0.01-grid argmin and the first-index tie-break are benign but unanalyzed departures from
the continuous argmin. *Validation:* dense-grid sensitivity + simulation coverage of
Eq. (7) under both exchangeable and correlated null designs.

**G4 (rank 4) — Second-order Gaussian knockoffs: model-X validity is an approximation.**
`knockoff_equi`/`knockoff_mvr` estimate (μ, Σ) once from the data
(R/artificial_features.R:395-408, 438-446) and draw i.i.d. knockoff rows; this is exact
model-X only when X is (approximately) multivariate Gaussian, which omic data typically
are not (R5's requirements are on the conditional law of X, not just moments). The
single-knockoff-realization reuse across all B bootstraps is Python-parity but means FDP+
carries unaveraged Monte Carlo noise from one draw. *Validation:* calibration under
non-Gaussian margins; repeated-knockoff-realization variance study.

**G5 (rank 5) — Random-permutation nulls do not satisfy Assumption 1 under dependence.**
Permuting columns destroys correlation structure, so permuted features are exchangeable
with real nulls only if nulls are (near-)independent; under strong omic correlation the
FDP+ numerator is biased. stablr warns only for the RF learners
(R/stabl_fit.R:369-377); for glmnet learners the default `artificial_type =
"random_permutation"` (R/stabl_fit.R:203) silently carries the risk, while the beginner
`stabl()` defaults to knockoffs (R/stabl.R:20-27). *Validation:* FDP+ calibration under
block-correlated null designs, permutation vs. knockoff.

**G6 (rank 6) — Strict `>` support rule vs. the paper's `≥` (and FDP+ counts).**
Both stablr (R/stabl_accessors.R:58-59) and Python use strict `>` where R1 defines Ŝ(t)
with `≥`. With discrete stability scores (multiples of 1/B) and a 0.01 threshold grid,
exact ties occur and features at exactly f_j = θ are excluded — a systematic
anti-conservative bias in |Ŝ(θ)| and in the FDP+ counts at grid points. Small effect,
but trivially fixable and a genuine spec deviation. *Validation:* targeted tie test
(row 8).

**G7 (rank 7) — RF / oblique-RF / ordinal learners: FDP+ calibration is empirical only.**
These learners are outside R1; the package itself warns that FDP+ is "empirical, not
guaranteed" for them and restricts them to model-X knockoffs
(R/stabl_fit.R:362-377, Rd of `stabl_fit`). Importance-threshold "paths" for RF replace
the λ path, so the max-over-grid score has a different (untheorized) null distribution.
*Validation:* the package's own RF validation protocol, extended to oblique RF and
ordinal; treat as research features until calibrated.

**G8 (rank 8) — Classification bootstrap repair changes the subsampling distribution.**
When a drawn subsample has <2 per class, stablr redraws (≤1000) and finally falls back to
stratified sampling (R/stabl_fit.R:558-660). This conditions the resampling distribution
on class feasibility — a deviation from the i.i.d. subsampling scheme underlying both R1
and R2 (benign in expectation for balanced problems, non-trivial for rare classes).
Failures of the fit itself, by contrast, fail closed (no denominator distortion — the
paper is silent and stablr's choice is the safe one). *Validation:* compare selection
frequencies with/without repair on imbalanced simulations.

**G9 (rank 9) — Cooperative engine: faithful for gaussian/binomial, but family-restricted
and ρ-semantics rely on the vendored fork.** The vendored multiview engine implements
R3's objective via augmented pair-rows and block coordinate descent with the glmnetpp WLS
point solver (src/wls_exp.cpp), with a machine-readable divergence log
(inst/PATCHES.cooperative.md) asserting no statistical change under stablr's usage (unit
weights, no offsets, gaussian/binomial only). The restriction and the untested
poisson/cox dead code are documented; the objective itself should be verified
numerically against the direct augmented-matrix solve (R3 Eqs. 23–24). *Validation:*
row 16's objective-value test.

**G10 (rank 10) — Minor cross-language numeric differences.** `round` vs. `int` for
n_injected; effective vs. nominal artificial_proportion in the FDP+ scaling
(R/stabl_fit.R:337-346 vs. Python `_compute_FDRc`); auto λ grids are glmnet-derived
(R/learner_adapters.R, `auto_lambda_grid`) rather than Python's `np.linspace(0.01, 1, 30)`
on the inverse-penalty C — same *kind* of object (a data-driven path), different
calibration; results are comparable only via the stability-path/FDP+ layer, not by
λ-value matching. *Validation:* cross-language equivalence harness on identical data.

---

### Summary judgment

The core Stabl machinery in stablr (subsampling, path accumulation, FDP+ formula, argmin
threshold rule, knockoff plumbing) is a **faithful, Python-parity translation** of R1/R9
with a small number of numerically equivalent or documented generalizations (1/π scaling,
rounding, strict-`>` support rule). The two scientifically consequential **unverified
extensions** are the environment-aware min/mean aggregation (no theory, no FDP control at
the aggregated level) and the RF-family/ordinal learners (empirical calibration only).
The most dangerous **interpretive** gap is the semantics of `min_fdr_` (q+ is random and
only relatively bounded, never a controlled FDR level). Failed bootstrap fits fail closed
and therefore cannot distort denominators; the classification-path feasibility repair is
the only place the resampling distribution itself is altered.
