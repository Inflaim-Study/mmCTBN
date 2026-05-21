# mmCTBN

> **Continuous-Time Bayesian Networks for Multimorbidity in Electronic Health Records**

An R package that unifies fast Maximum-A-Posteriori (MAP / Laplace) and
fully Bayesian (Stan / MCMC) estimation of Continuous-Time Bayesian
Networks (CTBNs) for the analysis of co-occurring chronic conditions in
longitudinal Electronic Health Record (EHR) data.

## Why `mmCTBN`?

Most existing CTBN tools handle only main effects and assume a
fixed number of nodes. `mmCTBN` is designed for **modern multimorbidity
research** and ships with **five interchangeable estimation backends**
under one common API:

  * **`map`** — fast regularised Poisson MAP (L-BFGS-B + Laplace SEs),
    with four prior families (`spike_slab`, `structured`, `lasso`,
    `horseshoe`) and up to 5-way interactions.
  * **`stan`** — full Bayesian MCMC via \pkg{rstan} (same priors).
  * **`classic`** — the original Nodelman et al. (2002) CTBN:
    closed-form MLE per (target, parent-configuration). No covariates.
  * **`fctbn`** — Faruqui et al. (2021) Functional CTBN. Adaptive
    group lasso (FISTA) with multiplicative parent effects and an
    optional Gaussian-mixture early-stop.
  * **`cph`** — Guillamet et al. (2025) CTBN-CPH. Classical baseline
    intensities modulated per edge by Cox proportional-hazards
    regressions, yielding individualised, covariate-dependent rates.

All five share:

  * **Any number of nodes** — `simulate_ctbn_data()` and the random
    network generator are size-agnostic.
  * **Patient-level k-fold CV** — `ctbn_cv()` plus
    `summarise_cv()` and `plot_tdauc()` for time-dependent AUC.
  * **Parallel everywhere** — `future.apply` for CV folds, parallel
    targets, and parallel patient simulation.
  * **Uniform downstream tooling** — `coef()`, `summary()`,
    `compute_F_m()`, network/heatmap/PIP/synergy plots.

## Installation

```r
# install.packages("remotes")
remotes::install_github("InflAim-Study/mmCTBN")
```

For the Stan backend you also need:

```r
install.packages("rstan")
```

## Quick start

```r
library(mmCTBN)

# 1. Simulate a 6-node network with up to 3-way interactions
nw <- make_random_network(n_nodes = 6, seed = 2026)
DT <- simulate_ctbn_data(network = nw,
                          n_patients         = 2000,
                          t_horizon          = 10,
                          interaction_order  = 3,
                          parallel = TRUE, n_cores = 4)

# 2. Fast MAP fit with spike-and-slab selection
fit <- ctbn_fit(DT,
                method            = "map",
                prior             = "spike_slab",
                max_order         = 2,
                fixed_covs        = c("sex_male"),
                time_varying_covs = c("age","smk_current", "smk_former"),
                parallel          = TRUE)
summary(fit)
plot(fit, type = "network")

# 3. Inference: cumulative incidence and synergistic effects
F5  <- compute_F_m(fit, target = "C1",
                    patient_dt = DT[eid == DT$eid[1]], tau = 5)
syn <- compute_interaction_effects(fit, DT, tau = 5)
plot_synergy_forest(syn, top_n = 15)

# 4. Cross-validation across priors
fit_fns <- list(
  ss  = function(dt, ...) ctbn_map(dt, prior = "spike_slab", max_order = 2,
                                    verbose = FALSE, ...),
  hs  = function(dt, ...) ctbn_map(dt, prior = "horseshoe",  max_order = 2,
                                    verbose = FALSE, ...))
future::plan(future::multisession, workers = 4)
cv <- ctbn_cv(DT, fit_fns = fit_fns, k_folds = 3,
              eval_times = c(2, 5, 8))
plot_tdauc(cv)
```

## Function reference

| Family               | Functions |
|----------------------|-----------|
| Fitting              | `ctbn_fit()`, `ctbn_map()`, `ctbn_stan()`, `ctbn_classic()`, `ctbn_fctbn()`, `ctbn_cph()` |
| Simulation           | `make_random_network()`, `make_default_network()`, `simulate_ctbn_data()`, `generate_true_interactions()` |
| Inference            | `compute_F_m()`, `compute_interaction_effects()`, `coef()`, `summary()` |
| Cross-validation     | `ctbn_cv()`, `make_patient_folds()`, `summarise_cv()`, `plot_tdauc()` |
| Recovery / selection | `compute_recovery_metrics()`, `compute_selection_metrics()`, `compute_oracle_metrics()`, `compute_pred_metrics()` |
| Plotting             | `plot_network()`, `plot_heatmap()`, `plot_pip_heatmap()`, `plot_synergy_forest()` |
| Utilities            | `prepare_wide()`, `build_design_matrix()`, `build_interaction_cols()` |

## Method comparison

| Method      | Covariates | Interactions | Selection            | Best when                                       |
|-------------|------------|--------------|----------------------|-------------------------------------------------|
| `map`       | yes        | up to 5-way  | PIP / kappa          | high-dim, want sparse synergy decomposition     |
| `stan`      | yes        | up to 5-way  | PIP / kappa          | need full posteriors / credible intervals       |
| `classic`   | no         | none         | none (uses topology) | small problems; closed-form, no tuning          |
| `fctbn`     | yes        | none         | adaptive group lasso | Faruqui-style structure learning on EHR data    |
| `cph`       | yes        | none         | none (uses topology) | per-patient individualised risk via Cox-PH      |

## Citation

If you use `mmCTBN` in published work, please cite the underlying
methodology papers and this package. See `citation("mmCTBN")`.

## License

GPL3.
