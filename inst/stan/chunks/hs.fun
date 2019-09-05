vector hs(vector z, vector r1_local, vector r2_local,
          real r1_global, real r2_global) {

    // number of penalized parameters
    int K = rows(z);

    // global shrinkage parameter
    real tau = r1_global * sqrt(r2_global);

    // local shrinkage parameters
    vector[K] lambda = r1_local .* sqrt(r2_local);

    // penalized regression coefficients
    return z .* lambda * tau;
}

vector reg_hs(vector z, vector r1_local, vector r2_local,
              real r1_global, real r2_global, real global_scale, real c2) {

    // number of penalized parameters
    int K = rows(z);

    // global shrinkage parameter
    real tau = r1_global * sqrt(r2_global) * global_scale;

    // local shrinkage parameters
    vector[K] lambda2 = square(r1_local .* sqrt(r2_local));

    // truncated local shrinkage parameter
    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));

    // penalized regression coefficients
    return z .* lambda_tilde * tau;
}
