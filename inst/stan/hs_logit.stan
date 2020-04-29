//
// Hierarchical shrinkage prior on regression coeffs (logistic regression)
//
// This implements the regularized horseshoe prior according to the simple
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

  // global shrinkage parameter
  real<lower=0> tau;

  // local shrinkage parameter
  vector<lower=0>[P-U] lambda;

  // auxiliary variables
  vector[P-U] z;
  real<lower=0> c2;
}

transformed parameters {

  // penalized regression parameters
  vector[P-U] beta_p;

  if (regularized)
    beta_p = reg_hs(z, lambda, tau, slab_scale^2 * c2);
  else
    beta_p = hs(z, lambda, tau);
}

model {

  // regression coefficients
  vector[P] beta = append_row(beta_u, beta_p);

  // half t-priors for lambdas and tau
  z ~ std_normal();
  lambda ~ student_t(nu, 0, 1);
  tau ~ student_t(global_df, 0, global_scale);

  // inverse-gamma prior for c^2
  c2 ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);

  // unpenalized coefficients including intercept
  beta_u ~ normal(0, scale_u);

  // likelihood
  y ~ bernoulli_logit_glm(X, 0, beta);
}
