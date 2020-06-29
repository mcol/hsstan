# hsstan 0.8 (29 June 2020)

### Major Changes

- Add the `sub.idx` option to `posterior_performance()` to select the
  observations to be used in the computation of the performance measures.
- Add the `start.from` option to run `projsel()` to start the selection
  procedure from a submodel different from the set of unpenalized covariates.
- Allow interaction terms in the formula for unpenalized covariates.
- Speed up matrix multiplications in `posterior_linpred()` and `projsel()`:
  this also benefits all other functions that use `posterior_linpred()`, such
  as `log_lik()`, `posterior_predict()`, `posterior_performance()` and others.

### Smaller Changes and Bug Fixes

- Fix parallelized loop boundaries in `posterior_performance()` for Windows.
- Speed up `posterior_performance()` for gaussian models.
- Handle correctly the case in which a variable is mentioned both among the
  unpenalized covariates and the penalized predictors.
- Fix bug in handling of a factor variable with multiple levels in the set of
  penalized predictors.
- Use the correct sigma term in the computation of the elpd for gaussian models.
- Allow running `projsel()` on models with no penalized predictors.

# hsstan 0.7 (1 May 2020)

### Major Changes

- Speed up all models up to 4-5 times by using Stan's `normal_id_glm()` and
  `bernoulli_logit_glm()`.
- Use a simpler parametrization of the regularized horseshoe prior.

### Smaller Changes and Bug Fixes

- Allow using the `iter` and `warmup` options in `kfold()`.
- Switch to `rstantools` 2.0.0.
- Fix bug in the use of the `slab.scale` parameter of `hsstan()`, as it was not
  squared in the computation of the slab component of the regularized horseshoe
  prior. The default value of 2 in the current version corresponds to using the
  value 4 in versions 0.6 and earlier.

# hsstan 0.6 (14 September 2019)

### Major Changes

- First version to be available on CRAN.
- Add the `kfold()` and  `posterior_summary()` functions.
- Implement parallelization on Windows using `parallel::parLapply()`.
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
- Validate the options passed to `rstan::sampling()`.
- Expand the documentation and add examples.

### Notes

- This version was used in:
  - [M. Colombo][mcol], S.J. McGurnaghan, L.A.K. Blackbourn et al.,
    Comparison of serum and urinary biomarker panels with albumin creatinin
    ratio in the prediction of renal function decline in type 1 diabetes,
    _Diabetologia_ (2020): 63 (4) 788-798.
    https://doi.org/10.1007/s00125-019-05081-8

# hsstan 0.5 (11 August 2019)

### Major Changes

- Update the interface of `hsstan()`.
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
- Allow to control the prior scale for the unpenalized variables.

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

- This version was used in:
  - [M. Colombo][mcol], E. Valo, S.J. McGurnaghan et al.,
    Biomarkers associated with progression of renal disease in type 1 diabetes,
    _Diabetologia_ (2019) 62 (9): 1616-1627.
    https://doi.org/10.1007/s00125-019-4915-0
  - [A. Spiliopoulou][athina], [M. Colombo][mcol], D. Plant et al.,
    Association of response to TNF inhibitors in rheumatoid arthritis with
    quantitative trait loci for CD40 and CD39,
    _Annals of the Rheumatic Diseases_ (2019) 78: 1055-1061.
    https://doi.org/10.1136/annrheumdis-2018-214877

# hsstan 0.1 (14 June 2018)

- First release.

[mcol]:   https://pm2.phs.ed.ac.uk/~mcolombo/
[athina]: http://www.homepages.ed.ac.uk/aspiliop/
