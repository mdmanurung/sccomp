# Inherited source baseline

The pre-implementation dirty checkout is archived at
`runs/baseline-2026-10-10/source.tar.gz`, with its SHA256 manifest, patch and git
state. Corrected runner attempt 5 executed the full inherited suite under R4_51:
663 passing assertions, six failed assertions, zero errors, seven skips.
Receipts: `runs/baseline-2026-10-10-attempt-5/test-results.{rds,csv}`.

All six failures concern assertions against attributes, not numeric differences:

- Two reproducibility assertions compared inherited fit tables with live CmdStan
  objects attached. Timestamps, temporary files and runtimes differ between fits.
  The separate exact draw-matrix assertions passed for pathfinder and HMC.
  Current assertions compare the named column lists, keeping exact value/type
  identity and the independent exact draw check.
- One prototype assertion compared a named matrix column with an unnamed vector.
  Current assertion checks sample row names explicitly and then the values.
- Three variability-design assertions expected vectors, whereas the inherited
  implementation returns one-dimensional integer arrays for Stan. Current
  expectations retain the array shape and exact mapping values.

Six skips are testthat's "empty test" notices after paths executed without an
expectation: outliers, multilevel Seurat, unwanted effects, SCE, metadata, and
variability plotting. One is the explicit identical-condition fold-change skip.
These do not establish assertion coverage. No CmdStan/dependency skip occurred.

Earlier attempts are retained as unsuccessful runner/dependency recovery evidence;
they do not supersede attempt 5. Frozen current-source full-suite execution and
exact count fixture regression are separate gates. No legacy likelihood or
count-model numerical implementation was changed to address these failures.

The three exact regression fixtures derive from packaged `seurat_obj` examples.
Its documentation does not establish source provenance. Regression receipts
compare the inherited dirty baseline with the modified API, not an external
cohort or a clean released checkout. Earlier receipt wording of "real-data"
was unsupported; immutable historical receipts are retained with this correction.
