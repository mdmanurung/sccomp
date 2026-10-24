# Current continuation: commit delivery and pilot evidence (2026-10-24)

The user authorized committing and pushing the existing changes on 2026-10-24.
This supersedes the historical no-commit/no-push notes below. Version stays 2.5.0;
no release or new scientific campaign is authorized.

The pilot finished: the 2026-10-18 snapshot records 108 HMC-passed fits and
six expected confounded-design rejections. Eight slope comparator bundles
failed at the known diffVar assertion; 16 passed and three were rejected as
expected. No pilot fits remain pending. Current package sources (103 files)
match the frozen source used for the full 2026-10-17 check. Fresh focused
continuous tests passed on 2026-10-24. See dev/continuous/STATUS.md for evidence.

Next: repair dropped-feature handling in the slope diffVar comparator and
freeze a new campaign identity before rerunning affected comparisons.
Scientific calibration and performance qualification remain open.

---

# Historical continuation: review fixes, full check, pilot cost accounting (2026-10-17)

Scope: publication-readiness review of the package and its uncommitted
continuous-outcome extension (REVIEW-ROADMAP.md stages 1-2). Older feature
proposals below are historical and must not replace it on resume. Details and
evidence for everything below: dev/continuous/STATUS.md section
"Package-check remediation and fixed-source verification (2026-10-17)".

Done this session (nothing committed; version stays 2.5.0):
- Fixed the quick-check regressions: S3 method signatures (legacy methods accept
  `...`; continuous methods keep legacy positions with fail-loud semantic-argument
  guards), `gt(0)` guard typo (R/sccomp_estimate.R:1101), dead
  `get_model_from_data` removed, failed imports/globals/Rbuildignore/DESCRIPTION,
  one Rd plus NAMESPACE regenerated, five new guard assertions.
- Full R CMD check of frozen `review-check-20261017-v1`
  (runs/review-package-check-20261017-v3, check phase 8,716 s): exit 0,
  1 WARNING, 0 NOTEs. The warning is environmental only (system TeX Live lacks
  inconsolata.sty for the indexed PDF manual). Tests 804 pass / 0 fail /
  7 skips; vignettes and examples OK.
- Pilot assessments: v1 4/114 hmc fits; v2 6/114 hmc plus 6/114 expected
  confounded design_rejected (expected_rejections_passed TRUE). Snapshots
  runs/review-progress-20261017-v1/-v2 (scientific_qualified FALSE by design).
- Cost proxies (matching 20x20 single-dataset, four cores per fit): fixed 86 s,
  intercept 792 s, slope 2,080 s; full grid >=12,500 core-hours. Slope timing
  on the current source is unmeasured.

Next: the user chose "complete the frozen-source pilot first"; it is running.
- Phase 2 snakemake PID 2954278, log
  `.local/continuous/review-check-v1-pilot-execution-phase2.log`, 11/114 fits
  at 15:33 CEST, ETA ~12-21 h (central ~17 h). A detached orchestrator
  (PID 3714086, `/tmp/opencode/wait-assess-pilot.sh`, log
  `.local/continuous/review-check-v1-pilot-assessment.log`) runs
  review-progress.R into a fresh date-unique snapshot when report.rds appears.
- Open finding: all slope-design diffVar comparators aborted at
  comparators.R:97 (`all(features %in% rownames(table))`); diffVar drops
  non-converged features. Location/robust/trend_robust variants retained.
  Decide (with the user) whether to fix comparators.R and re-freeze before
  the full grid; never rerun over frozen receipts.
- Other open gates (priced): full grid ~12.2k core-h (~8 d at 64 cores);
  SBC v1 800 fits ~710 core-h at pilot effort (budget undocumented); perf and
  airway benchmarks small but extrapolated. Do not launch without the user's
  choice.

Runtime: `/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript`;
`R_LIBS` under `.local/continuous/library`; PATH prefixed with R4_51/bin
(Pandoc 3.4); `MAKEFLAGS='PRECOMPILED_HEADERS=false'`;
`SCCOMP_STAN_CACHE_DIR=.local/continuous/baseline-cache`; XDG_CACHE_HOME for
snakemake. No commits, pushes, release, global installs or memory writes.

---

# Handoff — continuous sccomp: qualification evidence and "beyond limma" feature plan

**Date:** 2026-10-11 · **Branch:** master @ f340569 · **Status:** pilot v2 done and assessed; plan revised after 3 Biomni reviews; no new feature code written; nothing committed

## Goal
Publication-ready continuous Gaussian/Student-t sccomp (contract: `HANDOFF_CONTINUOUS_SCCOMP.md`;
gates: `dev/continuous/PLAN.md`), plus features that justify it over limma/dream for cytometry
cluster frequencies and scRNA pseudobulk. Version 2.6 needs ALL scientific/software/performance gates.

## Next action
User picks: "start feature 1" (R-side `known_variance` validation, `sccomp_asinsqrt_variance`,
unit tests; no Stan yet) or "3a first" (cost scaling at 20/50/100 features + marginalisation memo).
Default if unanswered: Feature 1. Plan: `/home/mdmanurung/.claude/plans/plan-implementation-and-comparisons-warm-hopcroft.md`.

## State
- Uncommitted: whole continuous extension (45 paths in `git status`), including this session's
  edits to `dev/continuous/comparators.R`, `fixed-truth-worker.R`, `run-comparator-rerun.R` (new),
  `CALIBRATION.md`, `STATUS.md`.
- Nothing of ours is running. (A `flashier-utils` snakemake in `ps` is unrelated.)
- Self-checks pass: `comparators.R --self-check`, `discovery-metrics.R --self-check`.
- Fixed-truth pilot v2 (`.local/continuous/campaigns/fixed-pilot-2026-10-10-v2`): 14/14 steps, all 6 sccomp
  fits `hmc_passed`; assessment `dev/continuous/runs/fixed-pilot-v2-assessment-2026-10-10.rds`
  (`all_comparators` FALSE because the original slope comparator receipt is `fit_failed`, preserved).
- Pilot result (n=1 per cell): location RMSE sccomp Gaussian/Student-t vs limma/dream: fixed .456/.450 vs .473;
  intercept .298/.292 vs .305; slope .227/.230 vs .225; Student-t/Gaussian location width ratio ~1.0,
  nu_obs est. 29-37. Parity, not superiority. sccomp 75-2151 s per fit vs limma/dream 1-10 s.
- Comparator reruns on frozen v2 inputs: `dev/continuous/runs/comparator-rerun-2026-10-10-{boundary,robust,trend_robust}`.
  dream slope has 11/20 boundary fits (intercept 3/20); robust/trend variants identical to plain on clean Gaussian data.
- Biomni review project `prj_0104acnWZY5weunZIigtijHI` (R1 `tsk_0159aBUvdCGiWmlm1HmvBdPj`,
  R2 `tsk_016w0359rSq28xq70Q4PG6eR`, R3 `tsk_0150vv3xbWVVqZd3k7c9Xtun`), all completed. Their numbers are
  unverified illustrations.

## Locked decisions
- Comparator accepts dream boundary fits by muffling only lmerTest's "negative eigenvalue" warning; flags them in a
  `boundary` column — user, 2026-10-10. Validation was a comparator-only rerun, not a v3 pilot.
- Defaults of any new feature must reproduce the current model; no edits to frozen campaigns — qualification rule.
- Competitors for the new features are crumblr+dream (frequencies) and dreamlet (pseudobulk), not plain limma — user, 2026-10-10.
- Benefit scenarios use count-level data (multinomial/Dirichlet-multinomial), never data from the fitted model — Biomni R1.
- Closure via ilr/CLR option, not factors; factors only for pseudobulk, deferred — Biomni R2.
- Pseudobulk scale: marginalise feature effects (HMC on hyperparameters) primary, two-stage empirical Bayes fallback;
  Pathfinder/Laplace/ADVI diagnostics only — Biomni R3 (needs the 3a memo to confirm for Student-t).

## Dead ends — do not redo
- Editing dream's `vpcontrol` (`check.conv.hess = "ignore"`): the warning comes from lmerTest's Satterthwaite step and dream
  sets `options(warn = 2)`. Use `withCallingHandlers` as in `comparators.R`.
- Binary slope design with 2 obs/donor: unidentifiable (as many random coefficients as observations); v2 uses 4 obs/donor.
- Plain asin-sqrt with `1/k` weights as the crumblr analogue: crumblr uses CLR with Dirichlet-multinomial variances.

## Read first
- `dev/continuous/STATUS.md` — evidence table (updated this session).
- `dev/continuous/CALIBRATION.md` — comparator amendments at the end.
- `dev/continuous/comparators.R:1-60` — trend/robust/boundary comparator.
- `inst/continuous-model.md:45-70` — likelihood and priors; model changes must be documented here.
- `inst/stan/continuous/functions.stan:26-70` — random coordinates and likelihood; `continuous/transformed_parameters.stan:17-22`.

## Open items
- Feature 1 vs "3a first" — user choice.
- Rank-thinning for SBC: keep 25 draws/chain (stride 60); max-ACF and rank ESS match the iid null (median .36 / 110); re-check after ~100 datasets. Full 800-dataset SBC not launched.
- Inspect 213 warnings and 7 skips of the full suite (`dev/continuous/runs/full-suite-2026-10-10-v2/test-results.csv`).
- Add robust/trend+robust comparator variants to the worker before a qualification freeze (currently plain only).
- Application (airway, frozen input `runs/airway-2026-10-10-v1/input.rds`) and 100/500/1000-feature benchmarks not started.

## Constraints
No commits, pushes, version bump, release or memory writes. Rscript:
`/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript`; `R_LIBS=$PWD/.local/continuous/library`;
`MAKEFLAGS='PRECOMPILED_HEADERS=false'`; `XDG_CACHE_HOME=$PWD/.local/continuous/snakemake-cache`;
`pkgload::load_all(..., helpers=FALSE)` for focused work.
