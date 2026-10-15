# WS-4 — Fusion benchmark summary

- Predeclared: MATURITY_PREDESIGN.md; 5 seeds; 100 bootstraps; lasso; random_permutation decoys.
- Held-out metric: BER (binomial) / MSE (gaussian) on the 50-sample validation split.
- Evaluation rule: late/cooperative use native valid_predictions; early uses a downstream
   stats::glm refit on the concatenated selected features (stabl_fit stores no reusable model).

## Mean validation performance and selection stability

           arm   family        mode performance          perf_source
 complementary binomial       early   0.2342814 early_downstream_glm
 complementary binomial        late   0.2642015          late_native
 complementary binomial cooperative   0.2430181   cooperative_native
 complementary gaussian       early   0.2634792 early_downstream_glm
 complementary gaussian        late   1.0898308          late_native
 complementary gaussian cooperative   0.2628399   cooperative_native
     redundant binomial       early   0.3548425 early_downstream_glm
     redundant binomial        late   0.2997018          late_native
     redundant binomial cooperative   0.3086285   cooperative_native
     redundant gaussian       early   0.2740361 early_downstream_glm
     redundant gaussian        late   0.2912674          late_native
     redundant gaussian cooperative   0.2579450   cooperative_native
 mean_seed_jaccard view2_share note
         0.7500000  0.50000000     
         0.2104985  0.51724138     
         0.3225595  0.52941176     
         0.5966667  0.50000000     
         0.2296967  0.39024390     
         0.2337588  0.53571429     
         0.1532345  0.50000000     
         0.1203109  0.41176471     
         0.2330159  0.09090909     
         0.2230303  0.57142857     
         0.2862034  0.50000000     
         0.1883333  0.42857143     

## Per-seed rows

           arm   family        mode performance          perf_source
 complementary binomial       early   0.2750000 early_downstream_glm
 complementary binomial       early   0.2152778 early_downstream_glm
 complementary binomial       early   0.2756410 early_downstream_glm
 complementary binomial       early   0.1683087 early_downstream_glm
 complementary binomial       early   0.2371795 early_downstream_glm
 complementary binomial        late   0.3583333          late_native
 complementary binomial        late   0.2118056          late_native
 complementary binomial        late   0.3125000          late_native
 complementary binomial        late   0.2027915          late_native
 complementary binomial        late   0.2355769          late_native
 complementary binomial cooperative   0.2750000   cooperative_native
 complementary binomial cooperative   0.1840278   cooperative_native
 complementary binomial cooperative   0.2948718   cooperative_native
 complementary binomial cooperative   0.1855501   cooperative_native
 complementary binomial cooperative   0.2756410   cooperative_native
 complementary gaussian       early   0.2214186 early_downstream_glm
 complementary gaussian       early   0.3704356 early_downstream_glm
 complementary gaussian       early   0.2707881 early_downstream_glm
 complementary gaussian       early   0.1915814 early_downstream_glm
 complementary gaussian       early   0.2631724 early_downstream_glm
 complementary gaussian        late   0.8760204          late_native
 complementary gaussian        late   1.4877049          late_native
 complementary gaussian        late   0.9216950          late_native
 complementary gaussian        late   1.2158186          late_native
 complementary gaussian        late   0.9479152          late_native
 complementary gaussian cooperative   0.2253609   cooperative_native
 complementary gaussian cooperative   0.3317087   cooperative_native
 complementary gaussian cooperative   0.2744905   cooperative_native
 complementary gaussian cooperative   0.2272005   cooperative_native
 complementary gaussian cooperative   0.2554389   cooperative_native
     redundant binomial       early   0.4423077 early_downstream_glm
     redundant binomial       early   0.2200000 early_downstream_glm
     redundant binomial       early   0.3636364 early_downstream_glm
     redundant binomial       early   0.3149351 early_downstream_glm
     redundant binomial       early   0.4333333 early_downstream_glm
     redundant binomial        late   0.2964744          late_native
     redundant binomial        late   0.2200000          late_native
     redundant binomial        late   0.3457792          late_native
     redundant binomial        late   0.2029221          late_native
     redundant binomial        late   0.4333333          late_native
     redundant binomial cooperative   0.2996795   cooperative_native
     redundant binomial cooperative   0.2400000   cooperative_native
     redundant binomial cooperative   0.3814935   cooperative_native
     redundant binomial cooperative   0.2386364   cooperative_native
     redundant binomial cooperative   0.3833333   cooperative_native
     redundant gaussian       early   0.2160151 early_downstream_glm
     redundant gaussian       early   0.3768704 early_downstream_glm
     redundant gaussian       early   0.2777841 early_downstream_glm
     redundant gaussian       early   0.2124875 early_downstream_glm
     redundant gaussian       early   0.2870233 early_downstream_glm
     redundant gaussian        late   0.2351725          late_native
     redundant gaussian        late   0.4119059          late_native
     redundant gaussian        late   0.3095566          late_native
     redundant gaussian        late   0.2178414          late_native
     redundant gaussian        late   0.2818607          late_native
     redundant gaussian cooperative   0.2192647   cooperative_native
     redundant gaussian cooperative   0.3192937   cooperative_native
     redundant gaussian cooperative   0.2930188   cooperative_native
     redundant gaussian cooperative   0.2122491   cooperative_native
     redundant gaussian cooperative   0.2458986   cooperative_native
 mean_seed_jaccard view2_share note
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>
                NA          NA <NA>

## Caveats

- Simulation-only; two-view designs; view-2 noise in the redundant arm is the deliberate null probe.
 - Cooperative fusion: built-in vendored engine, gaussian/binomial only (Cox deferred — see roadmap WS-4).
 - Early-fusion held-out performance depends on the documented downstream glm refit; it is not a
   package-native prediction path (stabl_fit stores no model).
