//
// Hierarchical shrinkage prior on regression coeffs (logistic regression)
//

data {

  // number of columns in model matrix
  int P;

  // number of unpenalized columns in model matrix
  int U;

  // number of training observations
  int N_train;

  // number of test observations
  int N_test;

  // prior standard deviation for the unpenalised variables
  int <lower=0> scale_u;

  // degrees of freedom of t distribution
  real <lower=1> nu;

  // X matrix for training data
  matrix[N_train, P] X_train;

  // X matrix for test data
  matrix[N_test, P] X_test;

  // binary response variable
  int<lower=0, upper=1> y_train[N_train];
}

parameters {

  // unpenalized regression parameters
  vector[U] beta_u;

  // auxiliary variable
  vector[P-U] z;

  // t priors are specified as mixtures: normal scaled by sqrt(~ inv_gamma)
  // gives better sampler than specifying t priors directly
  real<lower=0> r1_global;
  real<lower=0> r2_global;
  vector<lower=0>[P-U] r1_local;
  vector<lower=0>[P-U] r2_local;
}

transformed parameters {

  // penalized regression parameters
  vector[P-U] beta_p;

  // nested block to declare local variables
  {
    // global penalty parameter
    real tau;

    // local penalty parameters
    vector[P-U] lambda;

    tau = r1_global * sqrt(r2_global);
    lambda = r1_local .* sqrt(r2_local);
    beta_p = z .* lambda * tau;
  }
}

model {

  // linear predictor
  vector[N_train] mu;
  mu = X_train[, 1:U] * beta_u + X_train[, (U+1):P] * beta_p;

  // half t-priors for lambdas (nu = 1 corresponds to horseshoe)
  z ~ normal(0, 1);
  r1_local ~ normal(0.0, 1.0);
  r2_local ~ inv_gamma(0.5 * nu, 0.5 * nu);

  // half cauchy for tau
  r1_global ~ normal(0.0, 1.0);
  r2_global ~ inv_gamma(0.5, 0.5);

  // unpenalized coefficients including intercept
  beta_u ~ normal(0, scale_u);

  // likelihood
  y_train ~ bernoulli_logit(mu);
}
