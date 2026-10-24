# Continuous sccomp implementation and qualification

Started: 2026-10-10. Goal: publication-ready continuous Gaussian/Student-t
analysis, preserving released count models. Scientific contract:
`HANDOFF_CONTINUOUS_SCCOMP.md`, sections Locked product decisions through Public
API contract. This plan does not change those decisions.

Approved staged extension: [REVIEW-ROADMAP.md](REVIEW-ROADMAP.md). First strengthen
decisions and qualify v1 variance/nuisance handling; implement optional regularized
horseshoe effects only after v1 qualification. The review roadmap supersedes
the separate feature proposals recorded in the older HANDOFF.md.

## Scope and ownership

The starting dirty tree is the inherited Student-t prototype. Preserve it in an
immutable baseline snapshot before edits. Work is local; no commits, pushes,
release, shared-environment changes, or publication are authorized.

Changes must reach continuous Stan source and generated quantities, dedicated R
validation/model builder, prior constructor and estimator, S3 generics and
continuous methods, namespace/Rd, unit and integration tests, scientific
validation scripts/workflow, real-data provenance/example, and documentation.
Legacy beta/beta-binomial likelihoods and shared random-effect code are protected
until compatibility evidence supports any change.

## Ordered gates

1. **Baseline**: snapshot source and dirty diff; restore instantiate >=0.2.3 and
   cmdstanr >=0.9.0 in a local library; record runtime; execute complete existing
   suite and preserve results. Classify every failure/skip. Capture count-model
   schema/numerical fixtures. Baseline failures do not count as new successes.
2. **Model**: separate continuous core, exact inverse-gamma variance family,
   centered variance design, fixed effect priors, optional observed-mean trend,
   positive random-scale hierarchy and per-feature correlations, shared learned
   or fixed df, corrected Student-t scale, nested Gaussian branch. Pass stanc,
   compilation, analytic densities, prior sampling, and generated-quantity
   parameter matching. No residual cross-feature covariance claim.
3. **R boundary**: dedicated complete long-data validation; global affine
   scaling/back-transformation; continuous estimator/prior/class; source-aware
   cache; prediction/replication/residual/test/print/plot dispatch. Reject outlier
   refit and proportional fold changes. VI selections require exploratory opt-in.
   All continuous methods use location and variance semantics. Test duplicates,
   missingness, design rank, grouping and slopes, factor coding, units, schema,
   contrasts and posterior FDR families; preserve count regression fixtures.
4. **Scientific calibration**: freeze protocol before campaigns. Independent
   simulator; model-based SBC separate from fixed-truth coverage/discovery under
   misspecification. Explicitly handle observed-data scaling/trend conditioning.
   At least 200 fixed/100 mixed simulations, expanded when uncertainty requires.
   All attempted/failed fits retained. HMC: four chains, no divergence or
   treedepth saturation, Rhat<1.01, bulk/tail ESS>=400, E-BFMI>0.3. Coverage at
   90/95% within precomputed binomial Monte Carlo intervals. Inspect SBC ranks.
   Null FDR/false-sign behavior compatible with nominal targets within Monte
   Carlo uncertainty. Gaussian RMSE <=1.05 times matched limma/dream comparator;
   Student-t interval width <=1.10 times Gaussian unless coverage improves.
   Include contamination, correlated/hypervariable features, small donor counts,
   unbalanced groups, non-inverse-chi-squared variance distributions.
5. **Application and qualification**: provenance-backed real continuous assay,
   paired example, prior/posterior checks and sensitivity to donors/priors/df.
   Benchmark 100/500/1000 features and fixed/intercept/slope designs. Required
   envelope: 1000 features, 40 samples, 20 donors, four HMC chains, <=8h and
   <=64GB aggregate on recorded hardware. Snakemake owns campaign staleness;
   pilot measurements inform total compute before large submission. Complete
   test suite, R CMD check, Stan checks, regression fixtures, diff checks,
   documented scope/limitations and exact recovery commands. Version promotion
   to 2.6.0 only after all compatibility/scientific gates pass.

## Evidence rules

Artifacts, execution, technical correctness, scientific qualification,
performance, and acceptance are separate states. Failures remain failures; do
not weaken criteria or promote VI. A scientific-contract change requires a
concrete diagnosis and user review. No suitable real data leaves that gate open.
Archive commands, input/source hashes, runtime, priors/settings/seeds, diagnostics,
fit failures, metrics and wall time in immutable run directories under
`dev/continuous/runs/`. Large caches/libraries/build products stay local and out
of package source. Track current state in `dev/continuous/STATUS.md`.
