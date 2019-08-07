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
                 c("(Intercept)", hs.gauss$model.terms$unpenalized))
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

test_that("get.cv.performance works for non-crossvalidated models",
{
    out <- get.cv.performance(hs.gauss)
    expect_is(out, "data.frame")
    expect_equal(colnames(out),
                 c("set", "test.llk", "r2"))
    expect_equal(out$set, "Non cross-validated")
    expect_equal(out$test.llk,
                 -108.7622, tolerance=tol)
    expect_equal(out$r2,
                 0.2737298, tolerance=tol)

    out <- get.cv.performance(hs.binom)
    expect_equal(colnames(out),
                 c("set", "test.llk", "auc"))
    expect_equal(out$test.llk,
                 -33.28791, tolerance=tol)
    expect_equal(out$auc,
                 0.736, tolerance=tol)
})

test_that("get.cv.performance works for cross-validated models",
{
    out <- get.cv.performance(cv.gauss, out.csv="out.csv")
    expect_equal(nrow(out), length(folds) + 1)
    expect_equal(out$set,
                 c(paste("Fold", 1:length(folds)), "Overall"))
    expect_equal(out$test.llk,
                 c(-67.06889, -59.52169, -126.59058), tolerance=tol)
    expect_equal(out$r2,
                 c(0.016581920, 0.182042541, 0.005568317), tolerance=tol)
    expect_true(file.exists("out.csv"))
    unlink("out.csv")

    out <- get.cv.performance(cv.binom)
    expect_equal(out$test.llk,
                 c(-28.85088, -27.24616, -56.09704), tolerance=tol)
    expect_equal(out$auc,
                 c(0.4935065, 0.4935065, 0.4640000), tolerance=tol)
})

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
    expect_equal(nrow(out), P + 1)

    expect_equal(summary(hs.gauss, max.rows=0), out)
    expect_equal(nrow(summary(hs.gauss, max.rows=5)), 5)

    out <- summary(hs.gauss, sort="n_eff", decreasing=FALSE)
    expect_true(all(diff(out[, "n_eff"]) > 0))
})

test_that("print.hsstan",
{
    expect_output(print(hs.gauss))
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
