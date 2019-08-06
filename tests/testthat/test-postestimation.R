test_that("log_lik",
{
    out <- log_lik(hs.gauss)
    expect_equal(nrow(out), iters)
    expect_equal(ncol(out), N)

    expect_equal(log_lik(hs.binom),
                 log_lik(hs.binom, newdata=df))

    expect_equal(log_lik(cv.gauss[[1]]),
                 log_lik(cv.gauss[[1]], newdata=df[folds[[2]], ]))
})

test_that("posterior_interval",
{
    expect_error(posterior_interval(hs.gauss, pars=1:3),
                 "'pars' must be a character vector")
    expect_error(posterior_interval(hs.gauss, pars="zzz"),
                 "No pattern in 'pars' matches parameter names")
    expect_error(posterior_interval(hs.gauss, prob=0),
                 "'prob' should be a single number greater than 0 and less than 1")

    out <- posterior_interval(hs.gauss)
    expect_is(out, "matrix")
    expect_equal(nrow(out), P + 1)
    expect_equal(colnames(out), c("2.5%", "97.5%"))

    out <- posterior_interval(hs.gauss, pars="X1", prob=0.5)
    expect_equal(nrow(out), 2)
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

test_that("loo",
{
    out <- loo(hs.gauss)
    expect_s3_class(out, "loo")
})
