vector hs(vector z, vector lambda, real tau) {

    // penalized regression coefficients
    return z .* lambda * tau;
}

vector reg_hs(vector z, vector lambda, real tau, real c2) {

    // number of penalized parameters
    int K = rows(z);

    // local shrinkage parameters
    vector[K] lambda2 = square(lambda);

    // truncated local shrinkage parameter
    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));

    // penalized regression coefficients
    return z .* lambda_tilde * tau;
}
