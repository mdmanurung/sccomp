# Continuous inference qualification protocol

Supplemental review protocol: [REVIEW-ROADMAP.md](REVIEW-ROADMAP.md), `review-v1`.
It introduces independent variance/batch stress panels and comparison variants,
without changing v1 priors, likelihood, thresholds or acceptance criteria.

Protocol v1, frozen before continuous-model fitting on 2026-10-10. Targets derive
from `HANDOFF_CONTINUOUS_SCCOMP.md`. Amendments must be dated and distinguish
implementation corrections from changes to scientific choices or acceptance.

## Two separate questions

1. Does inference reproduce the posterior of the implemented generative model?
   Run SBC using independently written R simulation and `scaling="none"` with
   explicit matching priors. Trend is off for exact-model SBC. Known affine
   changes are tested with correctly transformed priors and results.
2. Does the practical workflow give useful uncertainty and discoveries for the
   intended assay/design? Run fixed-truth repeated-sampling experiments including
   the actual empirical global scaling and, separately, the observed-mean trend.
   These do not inherit a claim of exact prior-predictive SBC calibration.

No simulator calls a model-generated RNG or extracts its truth from fitted
outputs. Analytic density and simple-model reference checks are independent of
SBC to detect shared mistakes in simulator and likelihood.

## Fixed identities and replication

- Seed base: 10102026. Dataset identities encode family, design, scenario,
  replication, and protocol version. Pilot uses base+1,000,000 and qualification
  uses base+2,000,000, keeping both separate from development reference seeds.
  Distinct deterministic seeds are saved in
  a manifest before jobs run. Pilot identities are separate from qualification.
- SBC initial campaign: for each family, 200 fixed-effect datasets and 100
  random-intercept datasets. At least one correlated intercept/slope panel is
  required before claiming slope qualification; its replication is chosen from
  pilot precision/cost estimates before launching that panel.
- Default small calibration panel: 20 features, 40 samples, balanced condition;
  repeated panel has 20 donors with two observations each. Broader simulation
  panels cover 5/10/20 donors and 10/50 features. Feature-specific rank gates
  preselect features 1 and 10, avoiding treating dependent features as independent
  replicates. Report the remaining features as supplementary diagnostics.
- Retain all attempted fits. Do not replace a failed dataset with another seed.
  Retry only for a documented infrastructure failure using the identical input.
  Sampler failures remain in the denominator and require remediation.

## SBC and HMC

Store truth and thinned rank statistics for location coefficients, sigma_sq,
variance effects, log(s0), nu0, learned observation df, and random scales and
correlations when present. Assess both individual parameters and informative
derived quantities. Preserve unthinned draws for convergence checks.

Use four chains. Every representative reference fit must have zero divergences,
zero maximum-treedepth transitions, Rhat<1.01, bulk/tail ESS>=400 for all decision
and hyperparameters, and E-BFMI>0.3 in every chain. Sampling effort/adaptation may
increase to satisfy gates; record the change and rerun affected panels. Never
substitute VI or suppress difficult parameters. Tail probability decisions also
need adequate Monte Carlo precision near the selection boundary.

Inspect rank ECDFs/histograms against simultaneous simulation envelopes;
investigate systematic location/dispersion patterns. One uniformity p-value is
insufficient. Calculate 90% and 95% coverage separately for each predeclared
parameter/feature stratum across independent datasets. Precompute binomial
acceptance intervals using the nominal coverage and campaign size. Where many
strata are judged jointly, use a declared simultaneous diagnostic envelope.
Increase dataset replication if Monte Carlo uncertainty prevents interpretation;
do not treat a wide acceptance interval as strong precision evidence.

## Fixed-truth performance and robustness

Compare Gaussian fixed-effect estimates to limma and repeated-effect estimates
to variancePartition::dream, with identical continuous inputs, covariates,
contrasts, and independent units. Match the model's estimand; posterior FDR and
frequentist adjusted p-values are different decision rules. Report both effect
accuracy and their own discovery operating characteristics. Do not use different
preprocessing to give one method an advantage.

Required panels: independent Gaussian, Student-t heavy tails, sparse gross
contamination, correlated features, hypervariable feature subsets, unbalanced
conditions, few donors, and non-inverse-chi-squared residual variances. Include
Gaussian-limit Student-t fits, variance-effect null/non-null panels, no-trend and
linear-trend panels, and feature/factor-reordering invariance checks.

The unbalanced panel retains 40 samples/20 donors: ten donors have a matched
control/treated pair and ten have two controls (30 control/10 treated samples).
Within-donor slope variation remains in the ten paired donors. All default
fixed-truth panels except the five-donor scenario retain the 40-sample envelope.

Begin with 100 independent fixed-truth replicates per scenario after a pilot
cost estimate. Expand where the Monte Carlo confidence intervals cannot assess
the following targets. For clean Gaussian data, RMSE <=1.05 times the matched
comparator. Student-t interval widths <=1.10 times Gaussian unless coverage
improves. Report metrics in original input units, with uncertainty over datasets.

Predeclare original-unit location thresholds and variance-ratio thresholds in
the manifest before fitting. Thresholds describe the generating assay scale;
do not derive them from fitted effects. Discovery families are features within
component and contrast. Report posterior decision quantities and realized FDP
and false-sign proportion separately, setting FDP=0 for no discoveries.

For global-null panels, FDR equals the fraction of datasets with any discovery;
use precomputed binomial Monte Carlo acceptance intervals around the nominal
decision target. For mixed null/non-null panels, quantify mean FDP and false-sign
uncertainty by resampling independent datasets, not individual correlated
features. Retain counts of discoveries, true/false discoveries, no-discovery
datasets, failed fits, and near-threshold Monte Carlo instability. Report the
uncertainty interval and compatibility with the target; do not claim universal
frequentist FDR control from a posterior FDR value or one simulation panel.

## Application and compute

A real-data example must establish assay transformation, independent unit,
pairing, covariates, input identity/license/provenance, feature-selection rule,
contrasts, thresholds, priors, and missing-data/exclusion handling before fitting.
Evaluate prior and posterior predictive discrepancies for location, dispersion,
tails, condition/donor structure, and influential features/donors. Perform
sensitivity to reasonable prior changes, Gaussian/Student-t family, fixed/learned
df, and influential donors. Predictive checks alone do not establish discovery
calibration.

Four-chain CPU benchmarks: 100/500/1000 features with fixed, random-intercept, and
correlated random-slope designs. Required headline envelope: 1000 features,
40 samples, 20 donors with paired observations, condition plus donor intercept,
<=8h and <=64GB peak aggregate RAM. Measure the complete fit and result extraction
and report compilation separately. Record CPU/memory environment and chain/thread
allocation; pilots determine aggregate campaign cost before submission.

Snakemake owns dependency/staleness for campaigns. Versioned run roots freeze
the source and manifest; mutable current-checkout edits must not alter running
jobs. Persist logs, failures, draws or archival summaries with their source
identity, diagnostic receipts, runtime, metric tables, and gate evaluations.

## Identifiability amendment before qualification launch

2026-10-10: the matched dream pilot rejected the binary correlated-slope design
with two observations per donor because there are as many random coefficients
as observations. An independent covariance identity establishes structural
confounding: for each binary pair, `I = K_intercept - K_cross + 2*K_slope`.
See check-covariance-design.R and the preserved failed slope comparator receipt
under .local/continuous/campaigns/fixed-pilot-2026-10-10-v1/comparators.
Proper priors and passing HMC do not establish likelihood identifiability.

The locked contract already requires identifiable slopes. The continuous input
boundary now checks the rank of the combined random/residual covariance basis.
Future slope panels use four observations per donor, two per condition (80
samples for 20 donors, 20 for the five-donor stress panel). Fixed and intercept
panels retain two observations per donor and the 40-sample envelope. The slope
unbalanced panel uses 60 controls/20 treated observations. Earlier wording of
40 samples for all designs is superseded for slopes. Priors, thresholds,
replication minima and acceptance criteria are unchanged. Old two-observation
slope runs remain technical target/sampler evidence, not supported-design
qualification; their snapshots and failures are never overwritten.

The fixed-truth manifest predeclares 3,600 shared datasets and 8,100 fits:
12 scenarios x 3 designs x 100 replicates, both families; matched Gaussian-limit
Student-t fits for Gaussian data; trend-on fits alongside trend-off fits for the
trend scenario. The initial pilot has six fits on three Gaussian datasets.
Gaussian/Student-t share each input file and original-unit thresholds (location
1, residual variance ratio 1.5). The variance contrast is explicitly
`-2*condition1`, treated/control under centered sum coding. Matched limma/dream
receive the same values, metadata, minimum-effect threshold and independent units.

Coverage bounds are frozen for features 1/10 and both components, separately at
90/95%, across 324 strata using alpha/(4*324) tails. Null any-discovery bounds
are precomputed binomial upper counts with simultaneous tails across the 81
model-setting cells. Mean FDP, false-sign, interval-width and accuracy intervals
bootstrap independent datasets (2,000 resamples, seed 90102026); paired method
comparisons resample matched datasets jointly. Failed fits remain in attempted
denominators. Compatibility and adequate precision are reported separately;
uncertainty width >.05 for discovery error requires more replicates. Selection
changes under +/- two autocorrelation-aware probability MCSE are a separate
gate. Pilot size cannot qualify a panel.

2026-10-10: tied posterior directional-error probabilities share the cumulative
mean at the end of their tie group. This conservative correction prevents
feature order from determining selection. New campaigns use it; historical
pilot decisions predate it and are not discovery qualification evidence.

## References

### Implementation detail frozen before qualification launch

2026-10-10: numerical remediation preserves the physical model. Log-variance
deviations are scaled by sqrt(2/nu0), retaining the inverse-gamma Jacobian.
Random effects use an invertible Gaussian-based conditional coordinate map
with its exact Jacobian; prior-only draws retain the original whitened map.
Independent multivariate-normal and Student-t density checks cover the map.
HMC uses a narrower default unconstrained initialization range of .5, overridable
through CmdStan arguments. Earlier failed pilots/refits remain in their original
receipts; these changes do not amend the HMC or scientific acceptance thresholds.

2026-10-10: primary SBC strata are features 1/10 for each parameter plus shared
hyperparameters, with coefficient slots kept separate. This gives 91 strata
over the two families and three designs. Coverage bounds use binomial quantiles
at alpha/(4*91) and 1-alpha/(4*91), simultaneously over both 90/95% coverage
levels. Rank ECDFs use discrete uniform CDFs over 0..100 and the simultaneous
DKW bound sqrt(log(2*91/.05)/(2*n)). These conservative envelopes are diagnostics,
not strong precision evidence or replacements for visual inspection. Exact
95% binomial intervals for observed conditional coverage are reported alongside
attempted-fit denominators. Any failed fit or incomplete stratum closes the SBC
promotion gate. Pilot size never passes a qualification gate.

Coverage and envelope critical values must be archived with the qualification
manifest before its first fit. Each chain contributes 25 spaced draws to ranks;
spacing/random starts are saved. Inspect residual autocorrelation of these
thinned draws before interpreting ranks; increase sampling/thinning on the same
inputs if they are dependent. Discovery calculations use unthinned draws.

- Stan SBC: https://mc-stan.org/docs/stan-users-guide/simulation-based-calibration.html
- Stan diagnostics: https://mc-stan.org/learn-stan/diagnostics-warnings.html
- Stan predictive checks: https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html

## Comparator boundary-fit amendment

2026-10-10: the frozen fixed-pilot-2026-10-10-v2 slope comparator receipt is
`fit_failed`: `dream` dropped 2 of 20 features. dream sets `options(warn = 2)`, so
lmerTest's Satterthwaite warning "Model failed to converge with 1 negative
eigenvalue" on a singular (boundary) REML fit became an error, and the
completeness assertion rejected the dataset. `comparators.R` now muffles only that
warning, keeps the estimates, and flags every boundary fit in a `boundary` column;
all other dream checks and warnings are unchanged. The original receipt is
preserved in the denominator. Rerun on the same frozen inputs:
`runs/comparator-rerun-2026-10-10-boundary` (links.csv points to each original
receipt). Fixed and intercept results match the originals (all.equal, 1e-12).
The slope rerun has 20/20 valid rows, **11/20 boundary fits** (intercept design:
3/20). Satterthwaite df and intervals on boundary fits are less reliable, so
report the boundary rate with every dream comparison. Qualification freezes use
the amended comparator; this rerun is not itself a qualification receipt.

2026-10-10: `fit_continuous_comparator(trend =, robust =)` adds limma/dream variance-moderation
variants (`eBayes(trend = TRUE, robust = TRUE)`; the trend covariate is the feature mean and
`getTreat`/`treat` inherit the per-feature df). Estimates are unchanged; only intervals,
p-values and df move. Reruns on the frozen v2 inputs: `runs/comparator-rerun-2026-10-10-robust`
and `...-trend_robust`. On the clean Gaussian pilot they match the plain comparator (95% width
ratio 1.005 for fixed), so they matter only in the hypervariable, non-inverse-chi-squared and trend
scenarios of the full grid. The worker still runs the plain comparator; add the variants to the
qualification freeze before launching those scenarios.
