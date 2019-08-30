# Current development version

### Major Changes

- Add the `kfold()` and  `posterior_summary()` functions.
- Remove the deprecated `sample.stan()` and `sample.stan.cv()`.
- Replace `get.cv.performance()` with `posterior_performance()`.
- Report the intercept-only results from `projsel()`.
- Add options to `plot.projsel()` for choosing the number of points to plot and
  whether to show a point for the null model.

### Smaller Changes and Bug Fixes

- Cap to 4 the number of cores used by default when loading the package.
- Don't change an already set `mc.cores` option when loading the package.
- Drop the internal horseshoe parameters from the stanfit object by default.
- Speed up the parallel loops in the projection methods.
- Evaluate the full model in `projsel()` only if selection stopped early.
- Rename the `max.num.pred` argument of `projsel()` to `max.iters`.
- Expand the documentation and add examples.

# hsstan 0.5 (11 August 2019)

### Major Changes

- Update the interface to `hsstan()`.
- Don't standardize the data inside `hsstan()`.
- Implement the thin QR decomposition and use it by default.
- Replace uses of `foreach()`/`%dopar%` with `parallel::mclapply()`.
- Add the `posterior_interval()`, `posterior_linpred()`, `posterior_predict()`
  `log_lik()`, `bayes_R2()`, `loo_R2()` and `waic()` functions.
- Change the folds format from a list of indices to a vector of fold numbers.

### Smaller Changes and Bug Fixes

- Add the `nsamples()` and `sampler.stats()` functions.
- Use `crossprod()`/`tcrossprod()` instead of matrix multiplications.
- Don't return the posterior mean of sigma in the hsstan object.
- Store covariates and biomarkers in the hsstan object.
- Remove option for using variational Bayes.
- Add option to control the number of Markov chains run.
- Fix computation of fitted values for logistic regression.
- Fix two errors in the computation of the elpd in `fit.submodel()`.
- Store the original data in the hsstan object.
- Use `log_lik()` instead of computing and storing the log-likelihood in Stan.
- Allow the use of regular expressions for `pars` in `summary.hsstan()`.

# hsstan 0.4 (24 July 2019)

### Major Changes

- Merge `sample.stan()` and `sample.stan.cv()` into `hsstan()`.
- Implement the regularized horseshoe prior.
- Add a `loo()` method for hsstan objects.
- Change the default `adapt.delta` argument for base models from 0.99 to 0.95.
- Decrease the default `scale.u` from 20 to 2.

### Smaller Changes and Bug Fixes

- Add option to set the seed of the random number generator.
- Add computation of log-likelihoods in the generated quantities.
- Use `scale()` to standardize the data in `sample.stan.cv()`.
- Remove the standardize option so that data is always standardized.
- Remove option to create a png file from `plot.projsel()`.
- Make `get.cv.performance()` work also on a non-cross-validated hsstan object.
- Add `print()` and `summary()` functions for hsstan objects.
- Add options for horizontal and vertical label adjustment in `plot.projsel()`.

# hsstan 0.3 (4 July 2019)

### Major Changes

- Add option to set the `adapt_delta` parameter and change the default for all
  models from 0.95 to 0.99.
- Allow to control the prior scale for the unpenalised variables.

### Smaller Changes and Bug Fixes

- Add option to control the number of iterations.
- Compute the elpd instead of the mlpd in the projection.
- Fix bug in the assignment of readable variable names.
- Don't compute the predicted outcome in the generated quantities block.

# hsstan 0.2 (13 November 2018)

### Major Changes

- Switch to `doParallel` since `doMC` is not packaged for Windows.

### Smaller Changes and Bug Fixes

- Enforce the direction when computing the AUC.
- Check that there are no missing values in the design matrix.
- Remove code to disable clipping of text labels from `plot.projsel()`.

### Notes

- This version was used in Colombo, M., Valo, E., McGurnaghan, S.J. et al.
  Diabetologia (2019) 62: 1616. https://doi.org/10.1007/s00125-019-4915-0.

# hsstan 0.1 (14 June 2018)

- First release.
