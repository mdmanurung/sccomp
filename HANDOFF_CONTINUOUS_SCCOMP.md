# Handoff: production continuous-outcome sccomp

Date: 2026-08-21
Repository: `/exports/para-lipg-hpc/mdmanurung/sccomp`
Branch at handoff: `master`
Runtime: use `/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript`

## Next-session objective

Turn the current uncommitted Student-t proof of concept into a production-quality Bayesian continuous-outcome model for sccomp. The target is a robust Bayesian analogue of limma/dream for transformed continuous assays, with exact limma-family residual-variance moderation, optional repeated-measure random effects, and no learned shrinkage of biological effects in v1.

The next session should implement the staged design below, not merely polish the current prototype. Preserve the released beta/beta-binomial behavior and prove compatibility before promotion.

## Promotion state

Current state: **prototype executed, not scientifically validated, not promotion-ready**.

- The worktree began clean; all modifications listed below arose during this task and remain uncommitted.
- The current Student-t Stan likelihood compiles and has run fixed-effect, donor-random-intercept, prediction, residual, testing, and outlier-refit smoke paths.
- A small existing-family beta regression smoke path also ran.
- Those HMC smoke chains were deliberately short and reached maximum tree depth; they establish routing/execution only.
- The complete testthat/package suite has not run because `instantiate >= 0.2.3` is unavailable in `R4_51`. The installed `cmdstanr` is also older than the package's documented `>= 0.9.0` requirement.

Do not represent compilation, a completed chain, or the current tests as statistical validation.

## Locked product decisions

These decisions were explicitly selected during planning:

| Decision | Locked v1 choice |
|---|---|
| Continuous families | Gaussian and Student-t |
| Default family | Student-t |
| Main information sharing | Residual-variance pooling only |
| Biological-effect pooling | Excluded from v1 |
| Residual-variance prior | Fully Bayesian scaled inverse-chi-squared family |
| Differential variability | Retained as covariate-dependent log variance ratios |
| Mean-variance trend | Optional; off by default |
| Repeated measures | Random intercepts and correlated random slopes |
| Student-t degrees of freedom | One shared learned value; allow numeric fixed value for sensitivity tests |
| Input scaling | One global robust affine rescaling, reversed in all user-facing outputs |
| Public entrypoint | New `sccomp_estimate_continuous()` wrapper |
| User terminology | Location and variance; no composition naming |
| Result shape | Long component table |
| Authoritative inference | HMC |
| Variational inference | Exploratory only and explicitly labeled |
| Promotion envelope | Up to 1,000 features; four-chain HMC within 8 hours and 64 GB |

## Scientific contract

### Supported data

- Long `data.frame` input with `sample`, `feature`, and continuous `value` columns.
- Exactly one finite value for every sample-feature pair; reject duplicate or missing pairs in v1.
- Allow negative transformed values.
- Covariates and grouping variables must be constant within a sample across feature rows.
- Require at least two features; warn below ten because variance-hyperparameter learning will be weak.
- Require a full-rank fixed-effect design and identifiable random slopes.

The caller remains responsible for an assay-appropriate transformation. Raw RNA-seq counts, censoring, hurdle outcomes, exact-boundary proportions, and untransformed compositions are out of scope. Cell-frequency applications must supply scientifically justified coordinates such as ILR coordinates rather than treating marginal logits as generically non-compositional.

### Global transformation

With `scaling = "global"`:

1. Set `global_center` to the median of all finite input values.
2. Compute a MAD within each non-constant feature using the usual consistency factor 1.4826.
3. Set `global_scale` to the median positive feature MAD.
4. Fail if no positive scale exists; do not silently substitute feature-wise z-scores.
5. Fit `y_star = (y - global_center) / global_scale` and store both constants.

This is one affine transformation for the whole matrix, so relative feature variances remain intact. Back-transform every location, residual SD/variance, prediction, and residual before returning it. Log variance ratios and Student-t degrees of freedom are invariant.

With `scaling = "none"`, require an explicitly supplied prior object; do not apply the standardized defaults silently.

### Generative model

For sample `i` and feature `g`:

```text
mu_ig = X_i beta_g + Z_i b_g
log(v_ig) = log(sigma_g^2) + W_i gamma_g
```

`W` contains no intercept, uses centered numeric predictors and sum contrasts for factors, and is column-centered. Therefore `sigma_g^2` is the design-geometric-mean residual variance, not a reference-level artifact. Random effects are permitted only in `formula_location` in v1; reject random terms in `formula_variance`.

Gaussian family:

```text
y*_ig ~ Normal(mu_ig, sqrt(v_ig))
```

Student-t family:

```text
y*_ig ~ StudentT(nu_obs, mu_ig, sqrt(v_ig * (nu_obs - 2) / nu_obs))
```

The correction makes `v_ig` the actual residual variance rather than the Stan Student-t scale squared. Enforce `nu_obs > 2`.

### Exact limma-family variance moderation

For each feature:

```text
sigma_g^2 ~ InvGamma(nu0 / 2, nu0 * s0_g^2 / 2)
```

This is the scaled inverse-chi-squared family used by limma. Estimate its hyperparameters jointly:

```text
log(s0) ~ Normal(0, 1)
nu0 = 2 + exp(eta0)
eta0 ~ Normal(log(8), 1)
```

When `variance_trend = "none"`, set `s0_g = s0`. When `variance_trend = "linear"`, compute a fixed, centered, unit-SD feature-grand-mean covariate `t_g` on the globally transformed input and use:

```text
log(s0_g) = log(s0) + trend_slope * t_g
trend_slope ~ Normal(0, 0.5)
```

Document the trend as limma-trend-style conditional moderation because `t_g` is derived from the observed outcomes. Do not reuse or map arbitrary location coefficients to variance coefficients.

Variance effects use fixed weak priors, not learned cross-feature effect shrinkage:

```text
gamma_kg ~ Normal(0, 1)
```

A requested variance contrast is reported as a log variance ratio and optionally exponentiated to a variance ratio.

### Location and random-effect priors

On the globally scaled outcome:

```text
beta_intercept,g ~ Normal(0, 5)
beta_nonintercept,kg ~ Normal(0, 2.5)
```

These hyperparameters are fixed in v1. Do not learn `mu_beta` or `tau_beta`; that is the deliberately deferred effect-pooling extension.

For each random-effect coefficient `q` and feature `g`:

```text
log(omega_gq) = mu_omega,q + tau_omega,q * z_gq
z_gq ~ Normal(0, 1)
mu_omega,q ~ Normal(log(0.5), 1)
tau_omega,q ~ HalfNormal(0, 0.5)
```

Use per-feature Cholesky correlation factors with `LKJ(2)` for correlated intercept/slope clauses. `tau_omega` must be constrained nonnegative. The current shared implementation leaves its analogous parameter unconstrained, creating a sign-symmetric posterior; repair this before mixed-effect benchmarking.

Support the existing maximum of four random-effect clauses. Error for fewer than two grouping levels or a random slope with no within-group variation; warn for fewer than five grouping levels.

### Shared Student-t degrees of freedom

The default is learned globally:

```text
nu_obs = 2 + exp(eta_obs)
eta_obs ~ Normal(log(8), 1)
```

Also accept a numeric `student_t_df > 2` to enable fixed-df sensitivity and Gaussian-limit tests. Return a posterior summary for learned `nu_obs`.

Features remain conditionally independent given shared hyperparameters. Do not claim that v1 models cross-feature residual covariance.

## Public API contract

Add this entrypoint without changing released count-model behavior:

```r
sccomp_estimate_continuous(
  .data,
  formula_location = ~ 1,
  formula_variance = ~ 1,
  sample,
  feature,
  value,
  family = c("student_t", "gaussian"),
  variance_trend = c("none", "linear"),
  scaling = c("global", "none"),
  student_t_df = "estimate",
  prior = sccomp_continuous_prior(),
  inference_method = c("hmc", "variational"),
  ...
)
```

Add `sccomp_continuous_prior()` with the defaults stated in the scientific contract, including separate log-prior-df and observation-df parameters. Validate every SD as finite and positive and every LKJ shape as at least one.

Return class `c("sccomp_continuous_tbl", "sccomp_tbl", ...)`, but dispatch continuous methods before legacy methods. Refactor S3 generics to accept `...` where necessary so the following semantic arguments do not leak composition terminology:

- `sccomp_test(fit, contrasts, effect_threshold, variance_ratio_threshold = NULL, posterior_fdr = 0.05, ...)`
- `sccomp_predict(fit, formula_location = NULL, new_data = NULL, ...)`
- `sccomp_calculate_residuals(fit)`
- `sccomp_replicate(fit, formula_location = NULL, formula_variance = NULL, ...)`

The estimate and test tables are long-format with these stable columns:

```text
feature, component, parameter, estimate, lower, upper,
rhat, ess_bulk, ess_tail
```

Allowed `component` values are `location`, `log_variance`, `residual_variance`, `random_sd`, `random_correlation`, and `student_t_df`. A global parameter uses `feature = NA`.

Testing adds `threshold`, `direction`, `p_null`, `posterior_fdr`, and `selected`. `effect_threshold` is required and expressed in original input units. Test variance effects only when `variance_ratio_threshold > 1` is supplied, using `log(variance_ratio_threshold)` internally. Compute posterior FDR separately across features for each contrast and component; never combine location and variance discoveries.

For VI fits, `sccomp_test()` must require `allow_exploratory = TRUE`, retain an `inference_status = "exploratory"` marker, and never silently present VI selection as confirmatory.

Predictions return posterior locations and generated continuous values in original units. Residual output includes observed-minus-location residuals and standardized residuals based on the posterior residual SD. Replicates return `generated_value`, `location`, `residual_variance`, and `.draw`.

`sccomp_remove_outliers()` must fail clearly for continuous fits in v1. Student-t robustness plus posterior-driven removal would otherwise double-use observations. `sccomp_proportional_fold_change()` remains unavailable for continuous fits.

## Implementation sequence

### 1. Recover a trustworthy baseline

- Restore `instantiate >= 0.2.3` and `cmdstanr >= 0.9.0` in `R4_51`; confirm CmdStan 2.39.0 or the package-qualified version.
- Run the complete existing test suite before refactoring and preserve its receipt.
- Remove the untracked generated `glm_multi_student_t*.hpp` files and ignore generated Stan C++ artifacts; source `.stan` remains authoritative.
- Preserve all released `multi_beta_binomial` APIs and tests. Do not wholesale-revert the dirty worktree; it contains the useful prototype routing and tests.
- Write the above contract into package developer documentation before changing the model.

Gate: no implementation work is promoted until the full pre-change suite passes or every pre-existing failure is recorded and shown unrelated.

### 2. Replace the prototype with a nested continuous Stan core

- Replace the current log-normal log-precision hierarchy with the exact inverse-gamma representation above.
- Use one `glm_multi_continuous.stan` model with a family flag and zero-length Student-t parameters under Gaussian, so adjacent family fits differ only in their likelihood.
- Use one matching generated-quantities model for `mu`, `residual_variance`, `y_rep`, and pointwise `log_lik`; keep parameter blocks exactly aligned.
- Move reusable scale correction, variance-prior, inclusion-mask, and random-effect functions into a continuous Stan include rather than duplicating family implementations.
- Retain `reduce_sum`; benchmark grainsize rather than hard-coding a value optimized for the tiny fixture.

Gate: both family branches pass `stanc`, full C++ compilation, prior-only sampling, generated-quantities parameter matching, and analytic density spot checks.

### 3. Build the continuous R boundary

- Add a dedicated long-data validator and model-input builder; do not route continuous values through count/proportion inference or `variability_to_composition_map`.
- Implement the wrapper, prior constructor, semantic fit class, model routing, initialization, back-transformation, cache fingerprinting, printing, prediction, replication, residuals, testing, and plotting.
- Preserve the original input, scaling constants, formulas, contrasts, priors, family, df mode, seed, inference method, Stan source hash, CmdStan version, and model diagnostics as fit metadata.
- Make cache identity depend on package version, Stan source/include hashes, family, and compile flags so a stale executable cannot be reused after model changes.
- Bump the package minor version from 2.5.0 to 2.6.0 only after compatibility and scientific gates pass.

Gate: legacy count fixtures are byte/schema-compatible where expected, and continuous outputs never expose `composition`, `cell_group`, `c_`, or logit-fold-change terminology.

### 4. Verify statistical behavior before optimizing speed

- Unit-test validation, global scaling/back-transformation, design coding, factor releveling, random-effect clauses, prior validation, df modes, missing/duplicate rejection, long output schema, FDR grouping, and legacy dispatch.
- Test limiting behavior: Gaussian estimates agree with ordinary feature-wise `lm` at large sample size; large fixed Student-t df approaches Gaussian; releveling preserves requested contrasts; disabling variance covariates recovers homoscedastic residuals.
- Run simulation-based calibration for beta, sigma-squared, variance effects, variance hyperparameters, random-effect scales, and shared Student-t df.
- Benchmark against `limma` under independent Gaussian designs and `variancePartition::dream` under repeated Gaussian designs; add contaminated/heavy-tail scenarios where Student-t has a predeclared opportunity to help.
- Include misspecification scenarios with correlated features, hypervariable features, unbalanced groups, small donor counts, and non-inverse-chi-squared variance distributions.

HMC promotion gates on every representative reference fit:

```text
0 divergences
0 maximum-treedepth transitions
R-hat < 1.01 for all decision and hyperparameters
bulk and tail ESS >= 400 for all decision and hyperparameters
E-BFMI > 0.3 in every chain
```

SBC/coverage gates:

- Use at least 200 fixed-effect simulations and 100 mixed-effect simulations.
- Require 90% and 95% interval coverage to lie inside the precomputed binomial Monte Carlo acceptance interval.
- Inspect SBC rank histograms and require no systematic location or dispersion pattern; do not rely on one omnibus p-value.
- Under null simulations, require realized false-discovery and false-sign behavior to be compatible with the nominal posterior decision target within Monte Carlo uncertainty.
- Under correctly specified Gaussian data, Bayesian RMSE may not exceed the comparator by more than 5%; Student-t interval width may not exceed Gaussian by more than 10% unless coverage improves.

### 5. Qualify performance, documentation, and release

- Benchmark four HMC chains on 1,000 features, 40 samples, 20 donors with two repeated observations, one condition effect, and a donor random intercept.
- Promotion ceiling: wall time no more than 8 hours and peak aggregate memory no more than 64 GB in the recorded CPU-HPC environment.
- Record scaling curves for 100, 500, and 1,000 features and for fixed effects, random intercepts, and random slopes separately.
- Document supported transformations, priors, original-unit output semantics, HMC versus VI status, compositional exclusions, and a complete paired-donor example.
- Regenerate Rd/NAMESPACE with roxygen, run `R CMD check`, the complete testthat suite, Stan semantic checks, downstream beta-family regression fixtures, and `git diff --check`.

Final promotion requires separate evidence for artifact presence, execution, technical correctness, scientific calibration, performance qualification, and human acceptance. Do not infer one state from another.

## Current worktree map

Modified prototype files include:

```text
R/model_fitting.R
R/print_method.R
R/sccomp_calculate_residuals.R
R/sccomp_estimate.R
R/sccomp_predict.R
R/sccomp_proportional_fold_change.R
R/sccomp_remove_outliers.R
R/sccomp_replicate.R
R/sccomp_test.R
R/utilities.R
README.md
man/sccomp_calculate_residuals.Rd
man/sccomp_estimate.Rd
man/sccomp_predict.Rd
man/sccomp_test.Rd
vignettes/introduction.Rmd
```

Untracked prototype/artifact files:

```text
inst/stan/glm_multi_student_t.stan
inst/stan/glm_multi_student_t_generate_data.stan
inst/stan/glm_multi_student_t.hpp                       # generated artifact; remove
inst/stan/glm_multi_student_t_generate_data.hpp         # generated artifact; remove
tests/testthat/test-student-t-noise-model.R
```

The current Student-t prototype uses a Normal hierarchy on log precision and fixed `student_t_df = 4`. It is therefore not the locked production model. Reuse its continuous input validation, family routing, prediction, residual, and compatibility tests where they remain valid; replace its variance and df hierarchy.

The inherited random-effect helper in `inst/stan/common_functions.stan` uses an unconstrained `random_effect_sigma_sigma` and divides the log scale by a hard-coded `3.0`. Replace this for the continuous model with the explicit positive hierarchy in the contract. Avoid changing the beta-binomial implementation until separate compatibility evidence justifies a shared repair.

## Known risks and fail-closed boundaries

- Full Bayesian HMC may not meet the 1,000-feature envelope; if it fails, report the measured boundary rather than silently promoting VI.
- Exact inverse-chi-squared moderation plus heteroscedastic variance effects is an extension of ordinary limma. Keep the baseline variance and log-ratio effects separately identifiable and test factor-coding invariance.
- Learned Student-t df, residual variance, and random-effect variance can confound in sparse repeated designs. Use prior predictive checks and the specified grouping warnings.
- Global rescaling does not make fundamentally incomparable assay units exchangeable. Documentation must require one coherent assay scale; assay-class-specific variance priors are deferred.
- Cross-feature residual correlation is absent. Correlated-feature simulations test robustness but do not turn v1 into a multivariate residual model.
- Student-t outlier refitting is out of scope. Do not reuse the prototype's continuous outlier-removal route in production v1.

## Recommended skills for the next session

- `bioinformatics-scientific-workflow` for the scientific contract and layered promotion gates.
- `scientific-critical-thinking` for simulation design, calibration claims, and adversarial review.
- `andrej-karpathy-skills:karpathy-guidelines` for bounded implementation changes and verification.
- `r-modern-guide` and `r-style-guide` for the new R API and S3 methods.

## Exact next action

Restore the `R4_51` package environment and run the untouched full test suite, saving the command, runtime versions, and failures before refactoring any Stan or R source.
