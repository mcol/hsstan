//
// Hierarchical shrinkage prior on regression coeffs (linear regression)
//

data {

  // number of columns in model matrix
  int P;

  // number of unpenalized columns in model matrix
  int U;

  // degrees of freedom of t distribution
  real <lower=1> nu;

  // number of training observations
  int N_train;

  // number of test observations
  int N_test;

  // X matrix for training data
  matrix[N_train, P] X_train;

  // X matrix for test data
  matrix[N_test, P] X_test;

  // continuous response variable
  real y_train[N_train];
}

parameters {

  // unpenalized regression parameters
  vector[U] beta_u;

  // auxiliary variable
  vector[P-U] z;

  // residual standard deviation
  real <lower=0> sigma;

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

  // global penalty parameter
  real <lower=0> tau;

  // local penalty parameters
  vector <lower=0>[P-U] lambda;

  // linear predictor
  vector[N_train] theta_train;

  tau = r1_global * sqrt(r2_global);
  lambda = r1_local .* sqrt(r2_local);
  beta_p = z .* lambda * tau;
  theta_train = X_train[, 1:U] * beta_u + X_train[, (U+1):P] * beta_p;
}

model {

  // half t-priors for lambdas (nu = 1 corresponds to horseshoe)
  z ~ normal(0, 1);
  r1_local ~ normal(0.0, 1.0);
  r2_local ~ inv_gamma(0.5 * nu, 0.5 * nu);

  // half cauchy for tau
  r1_global ~ normal(0.0, 1.0);
  r2_global ~ inv_gamma(0.5, 0.5);

  // unpenalized coefficients including intercept
  beta_u ~ normal(0, 1000);

  // noninformative gamma priors on scale parameter are not advised
  sigma ~ inv_gamma(1, 1);

  y_train ~ normal(theta_train, sigma);
}

generated quantities {

  // predicted outcome on widthdrawn data
  vector[N_test] y_pred = X_test[, 1:U] * beta_u + X_test[, (U+1):P] * beta_p;
}
