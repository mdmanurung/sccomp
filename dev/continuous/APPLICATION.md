# Paired real-data example protocol

Frozen preprocessing/analysis choices before fitting, 2026-10-10.

Use `airway` 1.30.0 (Bioconductor 3.22, compatible with R 4.5), LGPL,
downloaded from
https://bioconductor.org/packages/3.22/data/experiment/src/contrib/airway_1.30.0.tar.gz.
Official source metadata:
https://bioconductor.org/packages/3.22/data/experiment/html/airway.html.
Experiment: Himes et al., PLoS ONE 2014, doi:10.1371/journal.pone.0099625,
GEO GSE52778. Eight libraries comprise four primary human airway smooth-muscle
cell lines measured with and without dexamethasone. The independent unit is
the cell line; preserve the within-line treatment pairing.

This example analyzes a continuous transformation of RNA-seq measurements.
It is an approximation for sufficiently expressed genes, not a raw-count
likelihood or a replacement for a validated RNA-seq differential-expression
workflow. Four independent units cannot establish broad mixed-model adequacy.
The example is explicitly separate from the 20-donor simulation/benchmark gates.

## Preprocessing and estimands

- Use the package's `airway` count assay and its `cell`/`dex` sample metadata.
  Archive package archive SHA256, assay/sample identifiers, package citation,
  metadata, library sizes, selected features and transformed matrix.
- TMM normalize library sizes with edgeR using all nonzero genes; calculate
  CPM for filtering and log2 CPM with prior.count=2 for continuous values.
- Require CPM >=1 and raw count >=10 in every one of the eight samples. Select
  the first 100 eligible Ensembl identifiers in lexical order, without looking
  at treatment effects or fitted results. Freeze the resulting identifiers.
- No samples excluded, no missing-value imputation. Assert exactly two samples
  per cell line and one treated/one untreated sample per line.
- Primary location contrast: treated minus untreated log2 CPM;
  `formula_location = ~ condition + (1 | donor)`.
  Primary variance model: `~1`. Its baseline variance is descriptive;
  no treatment-dependent variance discovery with four pairs.
- Student-t with learned shared df is primary; Gaussian is a sensitivity fit.
  Global scaling/default priors; original-unit location threshold=1 log2 CPM.
  Posterior expected error cutoff=5%, features within the treatment contrast.
  Selection remains blocked until that fit passes all HMC gates.

## Sensitivity and checks

Fit four chains with initially 1500 warmup/1500 sampling, adapt_delta=.995,
max_treedepth=12, fixed recorded seeds. Increase effort on the same dataset if
needed, preserving the failed receipt; never omit difficult features silently.
Compare Gaussian and Student-t effects/intervals, fixed Student-t df=8, wider
location-coefficient prior SD=5, narrower log-scale prior SD=.5, and each
leave-one-line-out fit (three independent units; explicitly weak evidence).
Report location, dispersion, tails, within-line treatment differences, and
largest standardized residual discrepancies in prior/posterior replication.
Plot/inspect influential genes and lines. Sensitivities do not establish FDR
calibration. Do not describe any selected gene as a newly established biological
discovery from this software example.
