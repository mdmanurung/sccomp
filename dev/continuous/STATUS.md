# Continuous sccomp status

Goal unfinished; no publication-readiness claim.

## Commit preparation and current evidence (2026-10-24)

- Version remains 2.5.0; the continuous API remains under qualification.
- The completed 2026-10-18 pilot snapshot records 108/114 HMC-passed fits
  and six expected confounded-design rejections. Comparator bundles record
  16 passes, eight failures (the known slope diffVar assertion), and three
  expected design rejections. No pilot fits remain pending. Scientific
  calibration and performance qualification remain open.
- Current DESCRIPTION, NAMESPACE, R, Rd, tests, Stan sources and vignette
  sources (103 files) match the frozen review-check-20261017-v1 snapshot
  byte for byte. Its full check records 804 passes, zero failures, seven
  skips and 246 test warnings; the sole check warning is missing
  `inconsolata.sty`. That full check was not rerun today.
- Fresh focused continuous tests passed in
  `runs/commit-check-20261024-v1`; `git diff --check` passed. Run outputs
  and local campaign snapshots remain ignored and are not part of delivery.
- Eight simulator/assessment self-checks passed. `comparators.R --self-check`
  reproduced the known `all(features %in% rownames(table)) is not TRUE`
  failure after dream reported seven failed responses. The comparator
  remains an open development defect; it was not repaired for this commit.

Next: repair the slope diffVar comparator's dropped-feature handling, then
freeze a new campaign identity before rerunning affected comparisons.
Full calibration/performance campaigns still need separate authorization.

## Package-check remediation and fixed-source verification (2026-10-17)

This session reviewed the package and the uncommitted continuous extension
toward publication readiness: quick-check regressions were fixed, a full check
was completed against a frozen source, and bounded pilot/cost evidence was
prepared for the full-grid decision. Nothing is committed; version stays 2.5.0.

Fixes (working tree only):
- S3 consistency (`tools::checkS3methods`): legacy methods
  `sccomp_predict/replicate/test.sccomp_tbl` now accept `...`; the continuous
  methods keep the legacy argument positions and route semantic arguments with
  fail-loud guards (`formula_composition`→`formula_location`,
  `formula_variability`→`formula_variance`, `robust`, and
  `percent_false_positive`/`test_composition_above_logit_fold_change`→
  `effect_threshold`/`variance_ratio_threshold`).
- `gt(0)` typo fixed in an `sccomp_estimate` uniqueness guard
  (R/sccomp_estimate.R:1101 → `nrow() > 0`).
- Dead `get_model_from_data` (undeclared rstan `stan_model`) removed from
  R/model_fitting.R.
- Package-level imports (`_PACKAGE` block) completed; new R/globals.R;
  DESCRIPTION imports `methods` and drops unused `boot`/`cli`/`fansi`;
  `^LICENSE$` added to `.Rbuildignore`.
- `man/sccomp_test.sccomp_continuous_tbl.Rd` and NAMESPACE regenerated with
  roxygen2 8.1.0 (canonical style, strict superset: +10 importFrom, no
  losses); five new guard assertions in tests/testthat/test-continuous-methods.R.

Verification (receipts preserved; no attempt overwritten):
- runs/review-package-check-20261017-v1/-v2 (quick static checks): 3 WARNINGs
  (including code/documentation mismatch) reduced to 2 WARNINGs, 0 NOTEs; the
  remaining two were `--no-build-vignettes` artifacts and vanish in the full run.
- runs/review-package-check-20261017-v3 (full build + check of frozen
  review-check-20261017-v1; check phase 8,716 s): exit 0; Status: 1 WARNING,
  0 NOTEs. The only warning is environmental: the system TeX Live 2018 lacks
  `inconsolata.sty`, which R's own Rd.sty requests for the indexed PDF manual
  (the manual without index is OK). Tests: 804 pass / 0 fail / 0 error /
  7 skips (six empty tests, one fold-change guard) / 246 warnings dominated by
  posterior ESS-cap notices, plus instantiate/posterior deprecation notices and
  one treedepth event; tests took 5,313 s. Examples and built plus re-built
  vignette outputs OK.

Review pilot evidence (receipts remain scoped to their frozen source):
- review-pilot-20261013-v1: 4/114 fits hmc_passed (fixed design, 96.7-239.8 s).
- review-pilot-20261013-v2: 6/114 fits hmc_passed plus 6/114 expected confounded
  rejections (`design_rejected` in 0.43-1.18 s), so `expected_rejections_passed`
  is TRUE; joint/omit/precorrect shows the intended log-variance discovery-FDP
  separation.
- Snapshots runs/review-progress-20261017-v1/-v2 keep
  scientific_qualified=FALSE (n=1 per cell; gates need >=100).

Cost accounting before any full-grid launch (2,700 datasets / 11,400 checks):
- Matching-dimension proxies from fixed-pilot-2026-10-10-v2 (20x20; each fit
  runs four chains on four cores): fixed 86 s, intercept 792 s, slope 2,080 s
  per fit; review fixed-design fits ran about 1.6x slower. Slope timing on the
  current source is unmeasured, so a bounded fixed-source slope/intercept batch
  must run before a full-grid decision.
- Projection: >=12,500 core-hours for the full grid (~8 days at 64 cores,
  ~4 days at 128); comparator bundles and simulation are negligible.

Full pilot completion run on the frozen source (launched 2026-10-17):
- Phase 1 serial precompile (--cores 4): 4/4 fixed-design fits hmc_passed;
  only two compiled model keys are created (keys are family-scoped; trend
  none/linear share a key).
- Phase 2 (--cores 8, two concurrent 4-core fits; PID 2954278) runs the
  remaining 110 fits plus summarize. At 15:33 CEST: 11/114 fits. Measured
  walls: fixed 105-261 s (n=4), slope 1972-3489 s (n=7, roughly 25-45 min
  each under a loaded shared node); intercept not yet observed. Projected
  remaining wall ~12-21 h (central ~17 h) from per-design means; completion
  expected 2026-10-18 morning. A detached orchestrator (PID 3714086) runs
  review-progress.R on completion and writes a uniquely named snapshot.
- Partial snapshot dev/continuous/runs/review-progress-20261017-v3 validates
  the assessment pipeline on partial state: 10 hmc_passed / 104 pending,
  Partial TRUE, scientific_qualified FALSE (by design).
- Finding: 8 of 24 non-confounded comparator bundles failed - every
  slope-design dataset - with `all(features %in% rownames(table)) is not TRUE`
  at comparators.R:97. Cause is code-side: variancePartition diffVar's
  per-feature AD/leverage refits legitimately drop non-converged
  (boundary/singular random-slope) features, so topTable returns fewer rows
  (14/20 or 9/20) and the assertion aborts that adjustment-diffVar variant.
  All location/robust/trend_robust variants and their fit files are retained;
  only the spread tables are missing. The assertion predates this campaign and
  slope diffVar comparators had never actually run. A fix (reindex to
  features, NA/valid flags instead of abort) belongs in comparators.R and
  would require a new campaign identity for any rerun; neither was changed
  here. Decide before any full-grid launch whether to fix and re-freeze.

## Bayesian review incorporation (2026-10-13)

Approved roadmap: [REVIEW-ROADMAP.md](REVIEW-ROADMAP.md). Older HANDOFF.md feature
proposals are historical; this session selected staged v1 qualification, followed
by optional regularized horseshoe effects. V1 likelihood/priors and released count
code have not changed in this incorporation work. No horseshoe API is implemented:
its prerequisite qualification remains open.

- Baseline frozen before edits: `.local/continuous/campaigns/review-baseline-20261013-v1`.
- Decision/nuisance documentation and the precise-small-effect regression test
  are implemented. `runs/review-unit-20261013-v1`: 128 assertions, zero failures,
  errors or skips.
- Variance/batch manifests, independent simulation, known-covariate formula routing,
  robust/trend comparator variants and separate diffVar spread comparisons are
  implemented. Self-checks: `.local/continuous/review-self-check-v2.log`,
  `review-comparator-check-v1.log`, `review-assessment-check-v2.log`; the latter's
  fixture row-name warnings were corrected and its adjustment-specific pairing
  passes with warnings treated as errors (`review-assessment-strict.log`).
- New pilot frozen: `.local/continuous/campaigns/review-pilot-20261013-v1`, 27
  datasets/114 fit-design checks. The first bounded target batch is four fixed
  Gaussian/Student-t fits on 0%/20% hypervariable Gaussian datasets, plus their
  comparison bundles. Dry-run passed. Two comparison bundles passed (plain,
  robust, trend/robust and diffVar); HMC batch started, not yet qualified.
- Full package build attempt `runs/review-package-check-20261013-v1` failed:
  PATH found system Pandoc 2.0.6, below rmarkdown's >=2.8 requirement. A fresh
  attempt v2 uses existing R4_51 Pandoc 3.4 via PATH; no installation or global
  environment change. Both attempts are preserved.
- No full SBC/fixed-truth or supplemental qualification launch. The supplemental
  full grid contains 2,700 datasets/11,400 fit-design checks; pilot diagnostics
  and compute estimates must precede launch. Package/application/performance
  gates and all staged model-extension work remain open.

Live execution logs: `.local/continuous/review-pilot-v1-execution.log` and
`.local/continuous/review-package-check-v2.log`. Frozen workers use their source
snapshots, not the mutable checkout. Later changes to assessment self-checks and
the package-check runner are outside the pilot source snapshot and do not change
the frozen likelihood, inputs or workers.

| Gate | State | Evidence / next action |
| --- | --- | --- |
| Baseline / counts | Inherited suite classified; count fixtures exact | Attempt 5: 663 passing assertions, six attribute-only failures, zero errors, seven skips (six empty tests after execution). See BASELINE.md. Three inherited packaged-example count fixtures reproduced exact inputs, draw matrices, estimates and test statistics at runs/count-regression-2026-10-10-v2. Their source provenance is undocumented; this compares the inherited dirty baseline rather than a clean released checkout or external cohort. Full corrected suite against frozen v2 source finished (runs/full-suite-2026-10-10-v2): 995 assertions, 0 failed, 0 errors, 7 skipped, 213 warnings (e.g. ESS capped); warnings and skips not yet individually inspected. |
| Model | Density proof and current-coordinate prior/API/GQ checks passed | Density receipt runs/densities-2026-10-10-conditional-random-v3 verifies likelihood, inverse-gamma/df-coordinate Jacobians, random-coordinate prior/Jacobian, exact Gaussian marginal identity, Student-t joint density and prior-only identity. Earlier density receipts had a matrix-shape assertion defect and are not valid automated evidence. Current-coordinate core routing passed both families and fixed/slope designs. Integration at runs/integration-2026-10-10-conditional-random passed both HMC references and API/GQ comparisons; the receipt audit confirms finite errors and 12 compared features per family. |
| R boundary | 115 unit assertions passed; both current-coordinate HMC references passed | runs/unit-2026-10-10-attempt-12, zero failures/errors/skips; source hashes are recorded, with a later label-only clarification of hyperparameter units. Covers original-unit hyperparameter summaries and rejection of poorly mixed multi-coefficient contrasts despite passing individual coefficients. Integration passed Gaussian fixed and Student-t correlated slopes, zero divergences/treedepth and all Rhat/ESS gates; E-BFMI .849-.973. Maximum API/GQ location difference 1.78e-15 and variance difference 4.44e-16. Integration predates metadata-unit/contrast-guard additions, which have focused tests. |
| Calibration | v3 pilot: all 6 fits HMC-passed; fixed-truth pilot v2 running | v3 pilot (immutable, corrected coordinates): all six fits passed every HMC gate; elapsed 103-3193 s (slope fits longest). Assessment runs/pilot-v3-assessment-2026-10-10-v2*: n=1 per cell so coverage/ECDF bounds are degenerate and sbc_campaign_passed=FALSE as expected; rank ESS 33-150 of 100 draws and rank ACF up to .45, so decide rank thinning and cost before freezing 800-dataset SBC. Fixed-truth pilot v2 (campaign fixed-pilot-2026-10-10-v2) finished 14/14 steps: all six sccomp fits hmc_passed (75-2151 s); assessment runs/fixed-pilot-v2-assessment-2026-10-10.rds (all_hmc TRUE, qualified FALSE; all_comparators FALSE because the original slope comparator receipt is fit_failed, preserved). Amended comparator reruns on the frozen inputs: runs/comparator-rerun-2026-10-10-boundary (11/20 slope, 3/20 intercept boundary fits), -robust and -trend_robust (identical to plain on clean Gaussian data). Pilot n=1 per cell, location RMSE sccomp Gaussian/Student-t vs limma/dream: fixed .456/.450 vs .473; intercept .298/.292 vs .305; slope .227/.230 vs .225. Student-t/Gaussian 95% width ratio 1.00 (location), 1.03-1.06 (log variance); estimated nu_obs 29-37 (posterior 9-108). Pilot only: gates need >=100 datasets per cell and the contamination/student_t/hypervariable/non_ig/trend scenarios. No full SBC/fixed-truth campaign launched. |
| Application/performance | Input frozen; fitting/qualification pending | airway 1.30.0 restored locally. Protocol APPLICATION.md; frozen 8-sample/4-cell-line/100-feature input at runs/airway-2026-10-10-v1/input.rds with archive SHA256/metadata/selection receipts. No application fit or 100/500/1000-feature benchmark yet. |

Runtime: `/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript`.
Contract: `HANDOFF_CONTINUOUS_SCCOMP.md`; execution plan: `dev/continuous/PLAN.md`.

Dependencies restored in .local/continuous/library, including cmdstanr 0.9.0
and instantiate 0.2.4. Shared CmdStan 2.39.0 is read-only; use
MAKEFLAGS='PRECOMPILED_HEADERS=false' for compilation. Never restart the active
baseline until its exact session/log has been checked. No version bump or release.

Earlier core attempt 4 used stale includes and is invalid for current-source
validation. Source hashes and forced recompilation corrected the core runner.

Current numerical remediation: df-scaled log variance plus exact Gaussian-based
conditional random-effect coordinates, preserving priors/likelihood through
Jacobians. Details in inst/continuous-model.md. The initial df-scaled refit
(96651) finished with an analysis-script syntax error after sampling;
completed chains include a divergence and E-BFMI=.290, so it cannot qualify.
Preserve its draws/log and do not replace its input/seed.
The successful conditional refit started before a parenthesization correction that
ensures symmetric precision updates at large magnitudes; its source receipt
must be interpreted accordingly. Repeat current-source matched GQ/prior routing
and sampler gates before scientific qualification. The corrected prior/GQ routing at
runs/core-2026-10-10-conditional-random passed fixed and slope cases for both
families: all 16 chains had zero divergences/treedepth, E-BFMI .844-1.14. This is
technical prior routing, not posterior/SBC qualification. No full SBC launch
until current integration and pilot diagnostics/cost are assessed.
