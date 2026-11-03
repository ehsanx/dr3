# dr3 0.0.0.9001

* **Correction.** Double cross-fitting (`DC-AIPW`, `DC-TMLE`) now uses one fixed rotation
  of roles, as in Zivich & Breskin (2021): observations in fold v are scored with the
  propensity model fit on fold v+1 alone and the outcome model fit on fold v+2 alone.
  Version 0.0.0.9000 fitted both role assignments and averaged their predictions, so the
  averaged propensity and outcome models each used both training folds. Recompute
  double-cross-fit results obtained with 0.0.0.9000; the other estimators are unchanged.
* `dr_nuisance(strategy = "double")` requires `Vd = 3`, which the rotation needs.
* `dr_tmle()` stops with an error when the outcome range is zero, since the outcome cannot
  then be rescaled to [0, 1].
* Tests check the double cross-fit's splitting rule without fitting a model, and the
  zero-range stop.
* The documentation explains seeding: full-sample fits use the current random-number
  state, and `dr_folds()` calls `set.seed()`. `?dr_rhc` describes how the Hmisc copy of
  the RHC data differs from the one analysed in the paper.
