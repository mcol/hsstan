test_that("hsstan",
{
    expect_error(hsstan(df, mod.gauss, chains=0),
                 "'chains' must be a positive integer")
    expect_error(hsstan(df, mod.gauss, adapt.delta=1),
                 "'adapt.delta' must be less than 1")

    expect_s3_class(hs.gauss,
                    "hsstan")
    expect_s4_class(hs.gauss$stanfit,
                    "stanfit")
    expect_equal(hs.gauss$family,
                 gaussian())
    expect_equal(names(hs.gauss$betas$unpenalized),
                 c("(Intercept)", "X1b", "X1c", "X2", "X3"))
    expect_equal(names(hs.gauss$betas$penalized),
                 hs.gauss$model.terms$penalized)
    expect_length(hs.gauss$betas, 2)

    expect_is(cv.gauss, "list")
    expect_length(cv.gauss, 2)
    expect_s3_class(cv.gauss[[1]],
                    "hsstan")
    expect_silent(validate.samples(cv.gauss[[1]]))
})

test_that("sample.stan",
{
    SW({
        ss.gauss <- sample.stan(df, df$y.gauss, unp, pen,
                                iter=iters, chains=chains)
    })
    expect_equal(names(ss.gauss),
                 names(hs.gauss))
    expect_equal(dim(ss.gauss$data), # sample.stan adds the hsstan_y_ column
                 dim(hs.gauss$data) + c(0, 1))
    expect_equal(ss.gauss$model.terms[3:4],
                 hs.gauss$model.terms[3:4])
    skip.check <- c("stanfit", "data", "model.terms")
    for (field in setdiff(names(hs.gauss), skip.check))
        expect_equal(ss.gauss[[field]],
                     hs.gauss[[field]])
})

test_that("sample.stan.cv",
{
    SW({
        sv.binom <- sample.stan.cv(df, df$y.binom, unp, pen, logit=TRUE,
                                   iter=iters, chains=chains, folds=folds)
    })
    expect_equal(names(cv.binom[[1]]),
                 names(sv.binom[[1]]))
    expect_equal(dim(sv.binom[[1]]$data), # sample.stan adds the hsstan_y_ column
                 dim(cv.binom[[1]]$data) + c(0, 1))
    expect_equal(sv.binom[[1]]$model.terms[3:4],
                 cv.binom[[1]]$model.terms[3:4])
    skip.check <- c("stanfit", "data", "model.terms")
    for (field in setdiff(names(cv.binom[[1]]), skip.check))
        expect_equal(cv.binom[[1]][[field]],
                     sv.binom[[1]][[field]])
})
