test_that("log_lik",
{
    out <- log_lik(hs.gauss)
    expect_equal(nrow(out), iters)
    expect_equal(ncol(out), N)

    expect_equal(log_lik(hs.binom),
                 log_lik(hs.binom, newdata=df))

    expect_equal(log_lik(cv.gauss$fits[[1]]),
                 log_lik(cv.gauss$fits[[1]], newdata=df[folds == 2, ]))
    expect_equal(log_lik(cv.binom$fits[[2]]),
                 log_lik(cv.binom$fits[[2]], newdata=df[folds == 1, ]))
})

test_that("posterior_interval",
{
    expect_error(posterior_interval(hs.gauss, pars=1:3),
                 "'pars' must be a character vector")
    expect_error(posterior_interval(hs.gauss, pars="zzz"),
                 "No pattern in 'pars' matches parameter names")
    expect_error(posterior_interval(hs.gauss, prob=0),
                 "'prob' must be a single value between 0 and 1")

    out <- posterior_interval(hs.gauss)
    expect_is(out, "matrix")
    expect_equal(nrow(out),
                 P + 1 + 1) # intercept and extra factor level for X1
    expect_equal(colnames(out), c("2.5%", "97.5%"))

    out <- posterior_interval(hs.gauss, pars="X1", prob=0.5)
    expect_equal(nrow(out),
                 3) # X1b, X1c, X10
    expect_equal(colnames(out), c("25%", "75%"))
})

test_that("posterior_linpred",
{
    expect_equal(posterior_linpred(hs.gauss),
                 posterior_linpred(hs.gauss, newdata=df))
    expect_equal(posterior_linpred(hs.binom, transform=TRUE),
                 posterior_linpred(hs.binom, transform=TRUE, newdata=df))
})

test_that("posterior_predict",
{
    expect_error(posterior_predict(hs.gauss, nsamples=0),
                 "'nsamples' must be a positive integer")
    expect_equal(posterior_predict(hs.gauss, seed=1),
                 posterior_predict(hs.gauss, seed=1, newdata=df))
    expect_equal(posterior_predict(hs.binom, seed=1),
                 posterior_predict(hs.binom, seed=1, newdata=df))
})

test_that("posterior_performance",
{
    expect_error(posterior_performance(x),
                 "Not an 'hsstan' or 'kfold' object")
    expect_error(posterior_performance(cv.nofit),
                 "No fitted models found, run 'kfold' with store.fits=TRUE")

    out <- posterior_performance(cv.gauss)
    expect_is(out,
              "matrix")
    expect_equal(rownames(out),
                 c("r2", "llk"))
    expect_equal(colnames(out),
                 c("mean", "sd", "2.5%", "97.5%"))
    expect_named(attributes(out),
                 c("dim", "dimnames", "type"))
    expect_equal(attributes(out)$type,
                 "cross-validated")
    expect_equivalent(out["r2", ],
                      c(0.00497231, 0.01238094, 0.00000000, 0.04022560),
                      tolerance=tol)
    expect_equivalent(out["llk", ],
                      c(-142.41211783, 23.39468621, -202.41609434, -119.25843542),
                      tolerance=tol)

    out <- posterior_performance(hs.binom, prob=0.89)
    expect_equal(rownames(out),
                 c("auc", "llk"))
    expect_equal(colnames(out),
                 c("mean", "sd", "5.5%", "94.5%"))
    expect_equal(attributes(out)$type,
                 "non cross-validated")
    expect_equivalent(out["auc", ],
                      c(0.66507840, 0.08794497, 0.54062400, 0.80817600),
                      tolerance=tol)
    expect_equivalent(out["llk", ],
                      c(-33.1964259, 3.22441630, -37.54544278, -27.68313103),
                      tolerance=tol)

    out <- posterior_performance(cv.binom, summary=FALSE)
    expect_equal(nrow(out),
                 nsamples(cv.binom$fits[[1]]))
    expect_equal(ncol(out), 2)
    expect_equal(attributes(out)$type,
                 "cross-validated")
    expect_equivalent(posterior_summary(out),
                      posterior_performance(cv.binom))
})

test_that("loo",
{
    out <- loo(hs.gauss)
    expect_s3_class(out, "loo")
    expect_equivalent(out$estimates[1:2, "Estimate"],
                      c(-114.087049, 7.010475),
                      tolerance=tol)

    out <- loo(hs.binom)
    expect_equivalent(out$estimates[1:2, "Estimate"],
                      c(-39.24420, 10.01289),
                      tolerance=tol)
})

test_that("waic",
{
    out <- waic(hs.gauss)
    expect_s3_class(out, c("waic", "loo"))
    expect_equivalent(out$estimates[1:2, "Estimate"],
                      c(-113.9345867, 6.8580124),
                      tolerance=tol)

    out <- waic(hs.binom)
    expect_equivalent(out$estimates[1:2, "Estimate"],
                      c(-38.7986589, 9.5673524),
                      tolerance=tol)
})

test_that("bayes_R2",
{
    expect_error(bayes_R2(hs.gauss, prob=1),
                 "'prob' must be a single value between 0 and 1")
    expect_error(bayes_R2(hs.gauss, prob=c(0.2, 0.5)),
                 "'prob' must be a single value between 0 and 1")

    out <- bayes_R2(hs.gauss)
    expect_is(out, "numeric")
    expect_named(out,
                 c("mean", "sd", "2.5%", "97.5%"))

    out <- bayes_R2(hs.binom, summary=FALSE)
    expect_is(out, "numeric")
    expect_length(out, iters * chains / 2)
})

test_that("loo_R2",
{
    out <- loo_R2(hs.gauss, summary=FALSE)
    expect_is(out, "numeric")
    expect_length(out, iters * chains / 2)

    out <- loo_R2(hs.binom)
    expect_is(out, "numeric")
    expect_named(out,
                 c("mean", "sd", "2.5%", "97.5%"))
})
