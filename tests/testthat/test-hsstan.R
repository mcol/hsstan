test_that("hsstan",
{
    expect_error(hsstan(df, mod.gauss, adapt.delta=1),
                 "'adapt.delta' must be less than 1")
    expect_message(out <- hsstan(df, mod.gauss, chains=0),
                  "the number of chains is less than 1")
    expect_null(out)

    expect_s3_class(hs.gauss,
                    "hsstan")
    expect_s4_class(hs.gauss$stanfit,
                    "stanfit")
    expect_false("r1_local[1]" %in% names(hs.gauss$stanfit))
    expect_equal(hs.gauss$family,
                 gaussian())
    expect_equal(names(hs.gauss$betas$unpenalized),
                 c("(Intercept)", "X1b", "X1c", "X2", "X3"))
    expect_equal(names(hs.gauss$betas$penalized),
                 hs.gauss$model.terms$penalized)
    expect_length(hs.gauss$betas, 2)
})

test_that("hsstan doesn't use the QR decomposition if P > N",
{
    SW({
        hs.noqr <- hsstan(df[1:5, ], mod.gauss, pen, iter=100, qr=TRUE,
                          keep.hs.pars=TRUE)
    })
    expect_false(hs.noqr$qr)
    expect_match(names(hs.noqr$stanfit),
                 "r1_local", all=FALSE)
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

    expect_named(cv.binom,
                 c("estimates", "pointwise", "fits", "data"))
    for (i in 1:max(folds))
        expect_s3_class(cv.binom$fits[[i]],
                        "hsstan")
    expect_silent(validate.samples(cv.binom$fits[[1]]))
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
