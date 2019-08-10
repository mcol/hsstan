test_that("summary.hsstan",
{
    expect_error(summary(hs.gauss, prob=c(0.2, 0.8)),
                 "'prob' must be a single value between 0 and 1")
    expect_error(summary(hs.gauss, prob=1),
                 "'prob' must be a single value between 0 and 1")

    out <- summary(hs.gauss)
    expect_is(out, "matrix")
    expect_equal(colnames(out),
                 c("mean", "sd", "2.5%", "97.5%", "n_eff", "Rhat"))
    expect_equal(nrow(out),
                 P + 1 + 1) # intercept and extra factor level for X1

    out <- summary(hs.gauss, pars="X1")
    expect_equal(rownames(out),
                 c("X1b", "X1c", "X10"))

    expect_equal(summary(hs.gauss, max.rows=0),
                 summary(hs.gauss))
    expect_equal(nrow(summary(hs.gauss, max.rows=5)), 5)

    out <- summary(hs.gauss, sort="n_eff", decreasing=FALSE)
    expect_true(all(diff(out[, "n_eff"]) > 0))
})

test_that("print.hsstan",
{
    expect_output(print(hs.gauss))
})

test_that("get.cv.performance works for non-crossvalidated models",
{
    out <- get.cv.performance(hs.gauss)
    expect_is(out, "data.frame")
    expect_equal(colnames(out),
                 c("set", "test.llk", "r2"))
    expect_equal(out$set, "Non cross-validated")
    expect_equal(out$test.llk,
                 -110.3405, tolerance=tol)
    expect_equal(out$r2,
                 0.2214606, tolerance=tol)

    out <- get.cv.performance(hs.binom)
    expect_equal(colnames(out),
                 c("set", "test.llk", "auc"))
    expect_equal(out$test.llk,
                 -32.66869, tolerance=tol)
    expect_equal(out$auc,
                 0.800, tolerance=tol)
})

test_that("get.cv.performance works for cross-validated models",
{
    out <- get.cv.performance(cv.gauss, out.csv="out.csv")
    expect_equal(nrow(out), max(folds) + 1)
    expect_equal(out$set,
                 c(paste("Fold", 1:max(folds)), "Overall"))
    expect_equal(out$test.llk,
                 c(-69.29462, -70.00014, -139.29476), tolerance=tol)
    expect_equal(out$r2,
                 c(0.01770498, 0, 0), tolerance=tol)
    expect_true(file.exists("out.csv"))
    unlink("out.csv")

    out <- get.cv.performance(cv.binom)
    expect_equal(out$test.llk,
                 c(-23.51937, -27.46995, -50.98932), tolerance=tol)
    expect_equal(out$auc,
                 c(0.5909091, 0.5649351, 0.5152000), tolerance=tol)
})

test_that("get.cv.performance works for one fold of a cross-validated model",
{
    out <- get.cv.performance(cv.gauss[[1]])
    expect_equal(out$set,
                 "Fold 1")
    expect_equal(out$test.llk,
                 get.cv.performance(cv.gauss)$test.llk[1])
    expect_equal(out$r2,
                 get.cv.performance(cv.gauss)$r2[1])
})

test_that("sampler.stats",
{
    out <- sampler.stats(hs.gauss)
    expect_is(out, "matrix")
    expect_equal(colnames(out),
                 c("accept.stat", "stepsize", "divergences", "treedepth",
                   "gradients", "warmup", "sample"))
    expect_equal(rownames(out),
                 c(paste0("chain:", 1:chains), "all"))
})

test_that("nsamples",
{
    expect_equal(nsamples(hs.gauss), iters * chains / 2)
})
