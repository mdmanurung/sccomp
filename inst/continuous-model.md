# Continuous-outcome model contract

Status: development contract; not evidence of publication readiness. Qualification
is tracked in `dev/continuous/STATUS.md`. Released composition models retain
their own assumptions and interpretation.

Implementation uses a df-scaled log-variance coordinate to reduce a sampler
funnel: log(sigma_sq) = 2*log(s0) + sqrt(2/nu0)*raw, with the corresponding
log Jacobian in the prior. This is an exact reparameterization of the stated
inverse-gamma family. The optional trend adds 2*slope*t to log(sigma_sq).

Random effects use exact Gaussian conditional coordinates during posterior
sampling. For each feature/block/group, let S = diag(omega)*L,
P = I + S'Z'V^-1 Z S, h = S'Z'V^-1(y-X*beta), and U = chol(P) (upper).
The physical latent z is solve(P,h) + solve(U,raw), and b = S*z. Its prior is
standard normal on z plus the log Jacobian -sum(log(diag(U))). The likelihood
and physical random-effect priors are unchanged, including for Student-t.
For multiple clauses this is a Gaussian-based coordinate proposal per block;
it does not integrate away other clauses. Prior-only sampling uses z=raw and
zero Jacobian. Shared includes keep estimation and generated quantities aligned.

## Data and independent units

Input is one finite continuous value per sample-feature pair, with no duplicate
or missing pairs. Negative transformed values are supported. Covariates and
grouping variables must be constant within each sample. At least two features
are required; fewer than ten warrants a variance-learning warning. The fixed
design must have full column rank. The sample/donor design determines the
independent units; repeated observations from one donor are not independent.

Supply one coherent, assay-appropriate transformed scale. Raw counts, censored
or hurdle outcomes, boundary proportions, and untransformed compositions are
outside v1. Compositional applications require justified coordinates such as
ILR, not arbitrary independent marginal logits.

## Scaling

Global scaling subtracts the median of all input values and divides by the
median positive within-feature MAD (consistency factor 1.4826). Fail if no
positive MAD exists. Do not standardize individual features. Restore original
units in location, prediction, residual, random SD, and residual-variance
outputs; log variance ratios and degrees of freedom are invariant. Scaling
`none` requires explicitly supplied priors appropriate to the assay units.

## Likelihood and variance

For sample i and feature g, mu[i,g] = X[i,] beta[,g] + Z[i,] b[,g], and
log(v[i,g]) = log(sigma_sq[g]) + W[i,] gamma[,g]. W has no intercept, uses sum
contrasts for factors, and its columns are centered. The baseline sigma_sq is
the design-geometric-mean residual variance. Random terms are permitted only
in the location formula.

Gaussian observations have mean mu and SD sqrt(v). Student-t observations have
location mu, shared df nu_obs > 2, and scale sqrt(v * (nu_obs - 2) / nu_obs).
Thus v is the actual residual variance in both families. Features are
conditionally independent given shared hyperparameters; cross-feature residual
covariance is not modeled.

## Priors on the globally scaled outcome

- sigma_sq[g] ~ InvGamma(nu0/2, nu0*s0[g]^2/2).
- log(s0) ~ Normal(0,1); nu0 = 2 + exp(eta0), eta0 ~ Normal(log(8),1).
- No trend: s0[g] = s0. Linear trend: log(s0[g]) = log(s0) + slope*t[g],
  slope ~ Normal(0,0.5), where t is the centered, unit-SD observed feature mean.
  This is conditional, limma-trend-style moderation, not a fully generative
  model for the trend covariate. Sharing the variance-prior family does not
  make inference identical to empirical-Bayes limma.
- Location intercept ~ Normal(0,5); other coefficients ~ Normal(0,2.5).
  These prior widths are fixed; biological-effect pooling is excluded in v1.
- Variance coefficients gamma ~ Normal(0,1), with no learned effect pooling.
- Random SD: log(omega[g,q]) = mu_omega[q] + tau_omega[q]*z[g,q], with
  z ~ Normal(0,1), mu_omega ~ Normal(log(0.5),1), tau_omega ~ HalfNormal(0,0.5).
  Per-feature random-effect correlation matrices use LKJ(2) Cholesky factors.
- Learned observation df: nu_obs = 2 + exp(eta_obs),
  eta_obs ~ Normal(log(8),1). A fixed numeric df > 2 is supported.

Allow up to four random-effect clauses. Reject fewer than two grouping levels
and random slopes without within-group variation; warn below five levels.

## Inference, decisions, and output

HMC is authoritative. VI is exploratory, labeled accordingly, and selection
requires explicit opt-in. Report feature, component, parameter, estimate,
lower, upper, rhat, ess_bulk, ess_tail. Components are location, log_variance,
residual_variance, random_sd, random_correlation, student_t_df. Global
parameters have feature NA. Preserve prior/hyperparameter summaries in fit
metadata in addition to the public component table.
Metadata hyperparameter summaries restore log(s0) and random log-SD means to
original input units; df, log-SD dispersion and trend slopes are unit invariant.
Retained raw Stan draws and the prior object describe the scaled model, with
the global transformation stored explicitly.

Require a positive location effect threshold in original units. Test variance
effects only with a variance-ratio threshold > 1, converted to log units.
Posterior FDR families are features within each component and requested
contrast; they do not imply study-wide control across additional contrasts.
Testing adds threshold, direction, p_null, posterior_fdr, selected.
Multi-coefficient HMC contrasts must also satisfy Rhat<1.01 and bulk/tail
ESS>=400; passing individual coefficients does not guarantee a well-mixed
linear combination. Nonfinite contrast values are rejected.
Tied posterior directional-error probabilities use the cumulative mean through
the entire tied group. Feature ordering cannot choose a subset of equal-risk
features. Simulation assessments also report decision changes under a two-MCSE
perturbation of those probabilities, separately from realized discovery errors.

The input boundary checks the rank of the marginal covariance basis, including
residual variance directions and all random covariance components. Correlated
intercept/binary-slope clauses with one control and one treated observation per
group cannot separate residual variance from random covariance and are rejected.
Four observations per group (two per condition) pass this necessary check.
Covariance-basis rank alone does not establish adequate statistical information;
priors and sparse-group warnings still require sensitivity analysis.

Predictions and replicates report continuous generated values, locations, and
residual variances in original units. Residuals include observed-minus-location
and posterior-SD standardized residuals. Continuous outlier refitting and
proportional fold changes are unavailable in v1.

## Qualification and reproducibility

Record original input, scaling, designs/formulas/contrasts, priors, family/df,
seed/settings, Stan source/include hashes, compiler/runtime identity, and
diagnostics. Compiled-cache identity includes source/include hashes, package
version, family and compile flags.

Qualification requires separate software, model, calibration, performance, and
application evidence. Exact-model SBC must generate parameters and observations
from the same fixed prior/likelihood, using an independent simulator. Do not
claim exact SBC for an empirical preprocessing workflow merely because the
underlying fixed-scale Stan core passes SBC.
