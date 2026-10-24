# Frozen campaign execution

Review extensions are specified in `REVIEW-ROADMAP.md`. `review_pilot` and
`review_truth` use the same freeze/worker/Snakefile interface, adding variance
stress and measured-batch panels. Assess a completed review campaign with
`assess-review.R CAMPAIGN NEW.rds`; expected rank rejections are tracked separately
from sampler failures, and a pilot cannot establish qualification. The full
review grid has 11,400 fit/design checks: inspect pilot cost before launching.

For package builds, prepend the R4_51 environment's `bin` directory to PATH.
It contains Pandoc 3.4; the system Pandoc 2.0.6 cannot build these vignettes.
The package-check receipt records PATH and detected Pandoc. Preserve failed
builds and use a fresh run ID after correcting the runtime.

Use R4_51 and a writable local dependency library. `freeze-campaign.R` snapshots
tracked and nonignored source, installs that snapshot into an isolated library,
records hashes/dependencies and freezes the seed/settings manifest. Workers
verify the snapshot hashes and load the installed snapshot, never the mutable
checkout. A campaign directory is never overwritten.

```bash
Rscript dev/continuous/freeze-campaign.R pilot-YYYY-MM-DD-v1 pilot
XDG_CACHE_HOME="$PWD/.local/continuous/snakemake-cache" \
  snakemake -s .local/continuous/campaigns/pilot-YYYY-MM-DD-v1/source/sccomp/dev/continuous/Snakefile \
  --config campaign=.local/continuous/campaigns/pilot-YYYY-MM-DD-v1 --cores 4 -n
XDG_CACHE_HOME="$PWD/.local/continuous/snakemake-cache" \
  snakemake -s .local/continuous/campaigns/pilot-YYYY-MM-DD-v1/source/sccomp/dev/continuous/Snakefile \
  --config campaign=.local/continuous/campaigns/pilot-YYYY-MM-DD-v1 --cores 4 --resources mem_mb=16000
```

In this workspace use the absolute Rscript and Snakemake paths from STATUS.md
and the R4_51 environment. Shared CmdStan is read-only; rule fit supplies
`MAKEFLAGS=PRECOMPILED_HEADERS=false`. Source-cache writes require the explicit
local `XDG_CACHE_HOME`; never redirect HOME. Four cores intentionally serialize
fits and initial compilation, avoiding simultaneous writes to a new cache key.
Precompile all cache keys before allowing multiple simultaneous fits.

The pilot has six fits: two families and three designs, each with 20 features,
20 donors, variance covariates, four chains, 1500 warmup/1500 sampling,
adapt_delta=.995. Source/seed identities differ from qualification identities.
Future correlated-slope panels have four observations per donor; two-observation
binary slope designs failed identifiability and are retained only as historical
technical evidence. Fixed/intercept panels keep two observations per donor.
`sbc` mode creates the initial qualification manifest: 200 fixed, 100 intercept
and 100 slope datasets per family. Do not launch it until pilot correctness,
diagnostics and cost have been inspected. This work has not yet qualified that
campaign or its decisions.

`fixed_pilot` and `fixed_truth` modes reuse the source freeze and Snakefile.
Simulation is shared across Gaussian/Student-t fits; comparator runs once per
dataset. The six-fit pilot precedes the 8,100-fit full grid and its compute
estimate. Assess completed runs with `assess-fixed-truth.R CAMPAIGN NEW.rds`.
Use absolute paths for explicit Snakemake targets because rule outputs are
absolute. Never rerun over a frozen receipt or suppress a comparator failure.

`run-package-check.R SNAPSHOT NEW-RUN-ID` builds vignettes and runs R CMD check
with tests/examples/manual enabled. `SCCOMP_STAN_CACHE_DIR` selects a writable
default before package loading; HOME and shared CmdStan remain unchanged.

Simulation and fit are separate rules. A fit receipt stores source/input hashes,
settings, runtime, HMC diagnostics, full fit, coverage/rank metrics, elapsed time
and errors. Sampler or analysis failures produce a failure receipt and remain
in the denominator. Snakemake benchmarks aggregate the worker and descendants'
RSS; report PSS as well when available. Its report is explicitly unqualified
until scientific gates are separately evaluated.

No failed fit is reseeded or replaced. Same-input sampling-effort increases get
new run identities and link to the original failed receipt. Infrastructure
recovery must preserve the original log/input/settings and explain the reason.
Never edit a running snapshot, touch outputs to mask staleness, or force a rerun
over frozen receipts. Check the live session and receipts before resuming.

The exact-model SBC path disables empirical scaling/trend and supplies matching
priors. Fixed-truth coverage/discovery, comparators, misspecification, application
sensitivity and 100/500/1000-feature benchmarks require separate campaign panels
and gate reports. A completed Snakemake graph alone is not publication readiness.
