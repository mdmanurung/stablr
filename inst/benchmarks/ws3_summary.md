# WS-3 — Oblique RF fix-or-fence summary

- Predeclared: MATURITY_PREDESIGN.md; nominal FDP+ 0.1; success criterion mean FDP <= 0.15 on BOTH B and B_het.

## Diagnosis: decoy selection-frequency noise proxy (design B, knockoff_equi, 10 seeds)

     learner decoy_mean_imp decoy_max_imp       fdp
1 oblique_rf      0.2621429         0.601 0.2421429
2         rf      0.2277143         0.486 0.0650000

## Remedies (oblique_rf + knockoff_equi)

  variant         remedy n_selected       fdp       tpr decoy_mean_imp
1       B       baseline        3.4 0.2421429 0.7666667      0.2621429
2   B_het       baseline        4.5 0.4504762 0.7333333      0.3138571
3       B R1_standardize        3.4 0.2421429 0.7666667      0.2621429
4   B_het R1_standardize        4.5 0.4504762 0.7333333      0.3138571
5       B     R2_boot200        3.1 0.1838095 0.7333333      0.2640714
6   B_het     R2_boot200        4.2 0.3921429 0.7333333      0.3183571
7       B        R3_both        3.1 0.1838095 0.7333333      0.2640714
8   B_het        R3_both        4.2 0.3921429 0.7333333      0.3183571
  decoy_max_imp
1        0.6010
2        0.6450
3        0.6010
4        0.6450
5        0.6010
6        0.6555
7        0.6010
8        0.6555

## Outcome

NO predeclared remedy met the criterion — adapter is FENCED: not-recommended scope goes into man pages/pkgdown; adverse result retained. 

## Caveats

- Remedies tested at benchmark level (input-level transforms), not runtime changes.
 - FDP+ is a diagnostic minimizer, not a universal FDR guarantee. Simulation-only.
