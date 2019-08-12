//
// Hierarchical shrinkage prior on regression coeffs (logistic regression)
//
// This implements the regularized horseshoe prior according to the alternative
// parametrization presented in:
//
// Piironen, Vehtari (2017), Sparsity information and regularization in the
// horseshoe and other shrinkage priors, Electronic Journal of Statistics

functions {
#include /chunks/hs.fun
}

data {

  // number of columns in model matrix
  int P;

  // number of unpenalized columns in model matrix
  int U;

  // number of observations
  int N;

  // design matrix
  matrix[N, P] X;

  // binary response variable
  int<lower=0, upper=1> y[N];

  // prior standard deviation for the unpenalised variables
  real<lower=0> scale_u;

  // whether the regularized horseshoe should be used
  int<lower=0, upper=1> regularized;

  // degrees of freedom for the half-t priors on lambda
  real<lower=1> nu;

  // scale for the half-t prior on tau
  real<lower=0> global_scale;

  // degrees of freedom for the half-t prior on tau
  real<lower=1> global_df;

  // slab scale for the regularized horseshoe
  real<lower=0> slab_scale;

  // slab degrees of freedom for the regularized horseshoe
  real<lower=0> slab_df;
}

parameters {

  // unpenalized regression parameters
  vector[U] beta_u;

  // auxiliary variables
  vector[P-U] z;
  real<lower=0> r1_global;
  real<lower=0> r2_global;
  vector<lower=0>[P-U] r1_local;
  vector<lower=0>[P-U] r2_local;
  real<lower=0> c2;
}

transformed parameters {

  // penalized regression parameters
  vector[P-U] beta_p;

  if (regularized)
    beta_p = reg_hs(z, r1_local, r2_local, r1_global, r2_global,
                    global_scale, slab_scale * c2);
  else
    beta_p = hs(z, r1_local, r2_local, r1_global, r2_global);
}

model {

  // linear predictor
  vector[N] mu = X[, 1:U] * beta_u + X[, (U+1):P] * beta_p;

  // half t-priors for lambdas and tau
  z ~ std_normal();
  r1_local ~ std_normal();
  r2_local ~ inv_gamma(0.5 * nu, 0.5 * nu);
  r1_global ~ std_normal();
  r2_global ~ inv_gamma(0.5 * global_df, 0.5 * global_df);

  // inverse-gamma prior for c^2
  c2 ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);

  // unpenalized coefficients including intercept
  beta_u ~ normal(0, scale_u);

  // likelihood
  y ~ bernoulli_logit(mu);
}
