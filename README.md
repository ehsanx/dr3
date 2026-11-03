# dr3

<!-- badges: start -->
<!-- badges: end -->

**Doubly-robust ATE estimators — three routes to one answer.**

`dr3` is a transparent, teaching-oriented implementation of the average
treatment effect (ATE) estimators reconciled in *Three Routes to One Answer:
Reconciling AIPW, TMLE, and Double Machine Learning for Applied Researchers*.
AIPW, TMLE, and double machine learning all target the **same** efficient
influence function

```
D(O) = [A/g(X) - (1-A)/(1-g(X))] (Y - Q_A(X)) + (Q_1(X) - Q_0(X)) - psi
```

and differ only in how they solve it (average the score vs. targeting) and how
the nuisances `g(X) = P(A=1|X)` and `Q_a(X) = E[Y|A=a,X]` are fit (full-sample,
cross-fit, or double cross-fit). `dr3` runs the whole family on **one shared
Super Learner library**, so differences between estimators come from the
combination rule and the fitting scheme (and, for split-based estimators, from
the random split), not from separate implementations.

Two choices generate the eight-estimator panel:

|                | full sample | single cross-fit | double cross-fit |
|----------------|-------------|------------------|------------------|
| **one-step**   | AIPW        | DML              | DC-AIPW          |
| **targeting**  | TMLE        | CV-TMLE          | DC-TMLE          |

plus the singly-robust baselines (g-computation, IPTW).

**Correction (0.0.0.9001).** Double cross-fitting now uses one fixed rotation of
roles, so the propensity and outcome models are trained on disjoint folds.
DC-AIPW and DC-TMLE results from 0.0.0.9000 should be recomputed; see
[NEWS](NEWS.md).

## Scope

`dr3` is for **teaching and reproducibility**. It is a readable reference engine,
whose tests check it against the established
[`tmle`](https://CRAN.R-project.org/package=tmle) and
[`AIPW`](https://CRAN.R-project.org/package=AIPW) packages; it is not a
replacement for them. It targets the ATE for a **binary treatment at a single
time point**; it is not a general causal-inference framework.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ehsanx/dr3", dependencies = TRUE)
```

## Quick start

```r
library(dr3)

d   <- dr_nhefs()                 # ships in causaldata; nothing to download
set.seed(1)                       # the full-sample fits use the current RNG state
fit <- dr_panel(d$Y, d$A, d$X)    # one call -> all eight estimators
fit
plot(fit)
```

`dr_panel()` is a thin wrapper: fit the shared nuisances yourself and hand them
to any estimator — pairing `dr_aipw()` with cross-fit nuisances *is* DML.

```r
nu <- dr_nuisance(d$Y, d$A, d$X, strategy = "crossfit")
dr_aipw(d$Y, d$A, nu)   # DML
dr_tmle(d$Y, d$A, nu)   # CV-TMLE
```

## Datasets for the tutorials

dr3 bundles no data: NHEFS comes from the `causaldata` package, and RHC is
downloaded from Vanderbilt's biostatistics repository.

| Loader | Data | Source |
|---|---|---|
| `dr_nhefs()` | Smoking cessation → weight change (continuous) | `causaldata` |
| `dr_rhc()`   | Right-heart catheterization → 30-day mortality (binary) | Vanderbilt's R-format export, downloaded by `Hmisc::getHdata("rhc")`; the paper analyses the CSV export, which differs slightly (see `?dr_rhc`) |

## Related

- Paper companion (regenerates the manuscript's reported results from archived outputs):
  <https://ehsanx.github.io/dml-alternatives/>
- Equivalence with the paper's engine: the companion repository's
  `extras/dr3_equivalence.R` runs `dr3` and the paper's hand-coded engine side by
  side on NHEFS and compares their estimates; see
  [its README](https://github.com/ehsanx/dml-alternatives#checking-dr3-against-the-engine).
- Other implementations: the tests check `dr3` against `tmle` and `AIPW`; the
  software-bridge article also runs `DoubleML`; the paper companion also compares
  `tmle3`. dr3's own nuisance fits all use `SuperLearner`; `DoubleML` fits ranger
  learners through mlr3, and `tmle3` an sl3 stack.

## Citation

Please cite the paper, [arXiv:2609.26142](https://arxiv.org/abs/2609.26142), and the
package, <https://github.com/ehsanx/dr3>; `citation("dr3")` gives both entries.

## Licence

MIT; see [LICENSE.md](LICENSE.md).
