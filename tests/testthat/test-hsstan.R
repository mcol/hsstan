test_that("hsstan",
{
    expect_s3_class(hs.gauss,
                    "hsstan")
    expect_s4_class(hs.gauss$stanfit,
                    "stanfit")
    expect_named(hs.gauss,
                 c("stanfit", "betas", "call", "data", "model.terms", "family",
                   "hsstan.settings"))
    expect_false("lambda[1]" %in% names(hs.gauss$stanfit))
    expect_equal(hs.gauss$family,
                 gaussian())
    expect_named(hs.gauss$betas$unpenalized,
                 c("(Intercept)", "X1b", "X1c", "X2", "X3"))
    expect_named(hs.gauss$betas$penalized,
                 hs.gauss$model.terms$penalized)
    expect_length(hs.gauss$betas, 2)
    expect_named(hs.gauss$hsstan.settings,
                 c("adapt.delta", "qr", "seed", "scale.u", "regularized", "nu",
                   "par.ratio", "global.scale", "global.df", "slab.scale",
                   "slab.df"))
    expect_equal(hs.gauss$hsstan.settings$global.scale,
                 0.007071068)
    expect_equal(hs.gauss$hsstan.settings$adapt.delta,
                 0.99)
})

test_that("hsstan with no penalized predictors",
{
    expect_null(hs.base$betas$penalized)
    expect_length(hs.base$penalized,
                  0)
    expect_named(hs.base$hsstan.settings,
                 c("adapt.delta", "qr", "seed", "scale.u"))
    expect_equal(hs.base$hsstan.settings$adapt.delta,
                 0.95)
})

test_that("hsstan handles categorical variables in the penalized predictors",
{
    SW({
        hs.1 <- hsstan(df, y.gauss ~ X2 + X3, "X1", iter=250,
                       keep.hs.pars=TRUE, refresh=0)
    })
    expect_named(hs.1$betas$unpenalized,
                 c("(Intercept)", "X2", "X3"))
    expect_named(hs.1$betas$penalized,
                 c("X1b", "X1c"))
})

test_that("hsstan handles penalized predictors appearing in the formula",
{
    SW({
        hs.1 <- hsstan(df, y.gauss ~ X1 + X2 + X3, "X2", iter=250,
                       keep.hs.pars=TRUE, refresh=0)
        hs.2 <- hsstan(df, y.gauss ~ X1 + X3, "X2", iter=250,
                       keep.hs.pars=TRUE, refresh=0)
    })
    for (val in c("beta", "data", "model.terms"))
        expect_equal(hs.1[[val]],
                     hs.2[[val]])
})

test_that("hsstan handles interaction terms correctly",
{
    SW({
        hs.int.2 <- hs(y.gauss ~ X1 + X3 + X2 + X1b_X3 + X1c_X3 + X3_X2, gaussian)
    })
    expect_equal(names(hs.inter$betas$unpenalized),
                 gsub("_", ":", names(hs.int.2$betas$unpenalized)))
    expect_equivalent(hs.inter$betas$unpenalized,
                      hs.int.2$betas$unpenalized)
    expect_equal(hs.inter$betas$penalized,
                 hs.int.2$betas$penalized)
})

test_that("hsstan handles interaction term without main effects",
{
    SW({
        hs.int.0 <- hsstan(df, y.gauss ~ X1:X3, iter=200, refresh=0)
    })
    expect_named(hs.int.0$betas$unpenalized,
                 c("(Intercept)", "X1a:X3", "X1b:X3", "X1c:X3"))
})

test_that("hsstan doesn't use the QR decomposition if P > N",
{
    SW({
        hs.noqr <- hsstan(df[1:5, ], mod.gauss, pen, iter=100, qr=TRUE,
                          keep.hs.pars=TRUE, refresh=0)
    })
    expect_false(hs.noqr$hsstan.settings$qr)
    expect_match(names(hs.noqr$stanfit),
                 "lambda", all=FALSE)
    expect_match(names(hs.noqr$stanfit),
                 "tau", all=FALSE)
})

test_that("kfold",
{
    expect_s3_class(cv.gauss,
                    c("kfold", "loo"))
    expect_output(print(cv.gauss),
                  "Based on 2-fold cross-validation")
    expect_named(cv.gauss,
                 c("estimates", "pointwise", "fits", "data"))
    expect_equal(rownames(cv.gauss$estimates),
                 c("elpd_kfold", "p_kfold", "kfoldic"))
    expect_equal(colnames(cv.gauss$estimates),
                 c("Estimate", "SE"))
    expect_equal(nrow(cv.gauss$pointwise),
                 N)
    expect_equal(colnames(cv.gauss$pointwise),
                 c("elpd_kfold", "p_kfold", "kfoldic"))
    expect_true(all(is.na(cv.gauss$pointwise[, "p_kfold"])))
    expect_length(cv.gauss$fits[[1]]$stanfit@stan_args,
                  2)

    expect_named(cv.binom,
                 c("estimates", "pointwise", "fits", "data"))
    for (i in 1:max(folds))
        expect_s3_class(cv.binom$fits[[i]],
                        "hsstan")
    expect_silent(validate.samples(cv.binom$fits[[1]]))
    expect_equal(nrow(cv.binom$fits),
                 2)
    expect_length(cv.binom$fits[[1]]$stanfit@stan_args,
                  1)
    expect_length(cv.binom$fits,
                  max(folds) * 2)
    expect_equal(cv.binom$fits[[max(folds) + 1]],
                 which(folds == 1))
})

test_that("kfold with store.fits=FALSE",
{
    expect_named(cv.nofit,
                 c("estimates", "pointwise"))
})

test_that("hsstan with invalid inputs",
{
    expect_error(hsstan(df, mod.gauss, adapt.delta=1),
                 "'adapt.delta' must be less than 1")

    expect_error(hsstan(df, mod.gauss, iter=0),
                 "'iter' must be a positive integer")
    expect_error(hsstan(df, mod.gauss, iter=-1),
                 "'iter' must be a positive integer")

    expect_error(hsstan(df, mod.gauss, warmup=0),
                 "'warmup' must be a positive integer")
    expect_error(hsstan(df, mod.gauss, warmup=-1),
                 "'warmup' must be a positive integer")

    expect_error(hsstan(df, mod.gauss, iter=1000, warmup=1000),
                 "'warmup' must be smaller than 'iter'")

    expect_error(hsstan(df, mod.gauss, chains=0),
                 "rstan::sampling failed")
})

test_that("kfold with invalid inputs",
{
    expect_error(kfold(hs.gauss, folds, chains=0),
                 "'chains' must be a positive integer")
    expect_error(kfold(hs.gauss, folds, chains=4.4),
                 "'chains' must be a positive integer")
})
