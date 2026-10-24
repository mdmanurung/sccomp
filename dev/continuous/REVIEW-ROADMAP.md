# Bayesian review incorporation: approved staged roadmap

Approved in the current session: qualify v1 first; then opt-in regularized
horseshoe biological-effect shrinkage. Qualify variance pooling before extending
it. Jointly model measured nuisance covariates using the existing formulas.
No latent factors, normalization-uncertainty model, feature strata or variance
mixture are included. This roadmap supersedes the unrelated feature proposals
in the older HANDOFF.md; it does not delete their historical record.

Source: Lewin, Bottolo and Richardson, *Bayesian methods for gene expression
analysis* (2018), sections 2.2, 2.3, 4.1 and 5.2:
https://drive.google.com/file/d/1467VwBxjs4Mg85mDvRELKzAs2dBtQgYJ/view.
Regularized horseshoe reference: Piironen and Vehtari (2017),
https://doi.org/10.1214/17-EJS1337SI.

## Stage 1: preserve v1 and explain its decisions

- [x] Freeze inherited source before edits: review-baseline-20261013-v1.
- [x] Explain minimum meaningful effects, posterior null/sign/directional error,
  tie handling and per-component/contrast discovery families in the vignette.
- [x] Test that a precise small effect remains unselected.
- [x] Finish current-source package/software qualification and classify all
  warnings and skips; existing receipts remain scoped to their frozen source.
  (2026-10-17: full check of frozen review-check-20261017-v1 = exit 0, 1 WARNING
  - environment-only TeX `inconsolata.sty` - and 0 NOTEs; test warnings/skips
  classified in STATUS.md.)

No public API, count model, continuous likelihood or prior changes in this stage.

## Stage 2: variance and measured-nuisance studies

`review-study.R` extends the independent fixed-truth generator, while the
existing source freeze, worker, bounds and Snakefile own reproducible execution.
All variants of one dataset use identical values, metadata and original-unit
thresholds. Corrected values are saved as an explicitly labeled diagnostic
analysis; the shared uncorrected input is never overwritten.

Variance panel: fractions 0%, 5%, 20% of features have tenfold residual variance;
Gaussian or standardized Student-t(df=5) observations; fixed, random-intercept
and identifiable correlated-slope designs. Fit both observation families with
and without the existing observed-mean trend. Hypervariable genes are fixed
across repeated datasets and recorded in the truth artifact.

Nuisance panel: Gaussian observations, fixed/intercept/slope designs, balanced
and partially confounded measured batch effects in both location and variance.
Compare joint batch adjustment, omitted adjustment, and batch-only least-squares
subtraction treated as fixed. The latter deliberately omits treatment when
estimating batch; it is a diagnostic of confounding and plug-in uncertainty,
not a recommendation or a claim that every two-stage method behaves this way.
Complete treatment-batch confounding is tested only for explicit rank rejection.
The airway example keeps its actual metadata; no batch is fabricated.

Comparators: plain, robust and trend/robust limma/dream minimum-effect tests;
`variancePartition::diffVar(method="AD", scale="leverage")` spread association
with BH adjustment. Spread results remain separate from minimum-variance-ratio
discoveries. Each comparator failure remains visible, including boundary warnings
not covered by the already approved narrow Hessian-warning handler.

The review pilot has 27 datasets and 114 fit/design checks; six are expected
rank rejections. The full supplemental grid has 2,700 datasets and 11,400
fit/design checks (100 datasets per setting). Seeds differ from prior campaigns.
Four chains, 1,500 warmup and sampling iterations, adapt_delta=.995, max_treedepth=12;
location threshold=1, residual variance ratio=1.5, posterior error cutoff=.05.
Freeze binomial bounds before fitting. Bootstrap independent datasets, never
correlated features; pair methods on shared datasets. Pilot results cannot qualify.

- [x] Implement the manifest, simulation, formula routing and comparison variants.
- [x] Record condition bias/RMSE, interval coverage, discovery errors and batch
  coefficient recovery; preserve expected rejections and attempted denominators.
- [ ] Validate all new scripts and frozen-worker routing.
- [ ] Run representative pilot targets and cost estimates before a full launch.
- [ ] Complete frozen v1 SBC, fixed-truth, real-data and performance gates.
- [ ] Complete supplemental qualification and inspect precision/calibration.

`assess-review.R` reports v1-supplemental evidence with scientific_qualified=FALSE.
Omitted/precorrected fits are misspecification controls, not promotion candidates.
If variance pooling fails qualification, stop promotion and document the cause;
do not silently adopt a broader variance family or weaken an acceptance gate.

## Stage 3: biological-effect shrinkage, gated on v1 qualification

Implementation must wait until all v1 software/scientific/application/performance
gates pass. Add estimator arguments `effect_shrinkage=c("none","regularized_horseshoe")`
and `shrinkage_coefficients=NULL`; default remains none. Enabled shrinkage requires
explicit fixed-design column names, excluding intercepts. Retain fixed priors
for unselected biological/nuisance coefficients, variance effects and random effects.

Extend the prior constructor with optional positive `effect_global_scale` and
`effect_slab_scale`, both explicitly required when enabling shrinkage and expressed
in the fitted model's outcome scale. For each selected coefficient k, share tau_k
across genes, with independent gene-specific lambda_kg:

```
beta_kg = tau_k * lambda_tilde_kg * z_kg
lambda_tilde_kg^2 = c^2 * lambda_kg^2 / (c^2 + tau_k^2 * lambda_kg^2)
z_kg ~ Normal(0,1)
lambda_kg ~ HalfCauchy(0,1)
tau_k ~ HalfCauchy(0,effect_global_scale)
c = effect_slab_scale  # fixed finite slab width
```

Use noncentered coordinates and matching generated quantities. Preserve beta
draw names, output schema, existing predictions/residuals/testing, legacy prior
objects and saved fits (interpreted as shrinkage-disabled). Record coefficient
selection and priors in metadata/cache identity; qualify new sampling parameters
and unit-restored scale summaries. Reject unknown/duplicate/absent coefficient
names. Record coding: recoding that changes the prior is not invariant by fiat.

## Stage 4: independent shrinkage qualification

Independent prior simulation/density/GQ tests; exact-model SBC with no empirical
scaling; practical fixed-truth panels with empirical scaling. Signal fractions
0%, 5%, 20%, 50%; include large and near-threshold effects, correlated genes,
and all qualified fixed/intercept/slope designs. Benchmark global scale=.25,
slab=2.5 on the globally scaled model; halve/double each separately for sensitivity.
These are frozen benchmark choices, not implicit production priors.

Reuse v1 replication minima, HMC and precision gates. Require paired sparse-panel
RMSE improvement with an uncertainty interval excluding equality, dense-panel
RMSE no more than 5% above the disabled model, and retained calibration/coverage/
large-effect recovery. Report failures, missing discoveries and complete runtime.
No release, version promotion, commit or push is authorized by this local work.
