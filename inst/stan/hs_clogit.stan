//
// Gaussian prior on regression coeffs (conditional logistic regression with 1 case per stratum)
//

functions {
#include /chunks/hs.fun
}

data {
  // prior standard deviation for the unpenalised variables
  real <lower=0> scale_u;

  // number of columns in model matrix
  int P;

  int N; // number of observations
  int U; // number of unpenalized columns in model matrix
  int S;  // number of strata

  // store ragged array as flat array
  int starts[S];
  int stops[S];

  matrix[N, P] X;   // design matrix

  // categorical response variable taking values in 1:lengths[s]
  int<lower=1> ycat[S];

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

transformed data {

  // categorical response variable taking values in 1:lengths[s]
  //  int<lower=1> ycat[S];

  //  for(s in 1:S) {
  //   ycat[s] = which(y[starts[s]:stops[s]]);
  // }

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

  // local variables
  vector[P] beta = append_row(beta_u, beta_p); // regression coefficients
  vector[N] X_beta = X * beta; // linear predictors

  // half t-priors for lambdas and tau
  z ~ std_normal();
  lambda ~ student_t(nu, 0, 1);
  tau ~ student_t(global_df, 0, global_scale);
  
  // inverse-gamma prior for c^2
  c2 ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  
  // unpenalized coefficients including intercept
  beta_u ~ normal(0, scale_u);

  // conditional likelihood is evaluated by looping over strata
  for (s in 1:S) {
    ycat[s] ~ categorical_logit(X_beta[starts[s]:stops[s]]);
  }

}
