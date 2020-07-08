//
// Gaussian prior on regression coeffs (conditional logistic regression with 1 case per stratum)
//

functions {

  int which(int[ ] x) { // return first subscript for which x[i] > 0)
    int ycat;
    int K; 
    int i;
    ycat = K + 1; // assign to an out of range value
    K = size(x);
    i = 1;
    while (ycat == K + 1 && i <= K) {
      if(x[i] > 0) {
	ycat = i;
      }
      i = i + 1;
    }
    return ycat;
  }

}

data {
  // prior standard deviation for the unpenalised variables
  real <lower=0> scale_u;

  int N; // number of observations
  int U; // number of unpenalized columns in model matrix
   int S;  // number of strata

  // store ragged array as flat array
  int starts[S];
  int stops[S];

  matrix[N, U] X;   // design matrix
  // int <lower=0, upper=1> y[N]; // binary outcome

  // categorical response variable taking values in 1:lengths[s]
  int<lower=1> ycat[S];

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
}

transformed parameters {
  vector[N] X_beta;
  
  X_beta = X * beta_u;
  
}

model {
  // unpenalized coefficients including intercept
  beta_u ~ normal(0, scale_u);
  // conditional likelihood is evaluated by looping over strata
  for (s in 1:S) {
    ycat[s] ~ categorical_logit(X_beta[starts[s]:stops[s]]);
  }
}
