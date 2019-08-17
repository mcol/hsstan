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

test_that("posterior_summary",
{
    out <- posterior_summary(1:100)
    expect_is(out,
              "matrix")
    expect_equal(nrow(out), 1)
    expect_equal(colnames(out),
                 c("mean", "sd", "2.5%", "97.5%"))
    expect_equal(as.numeric(out),
                 c(50.5, 29.01149198, 3.475, 97.525))

    out <- posterior_summary(hs.binom)
    expect_equal(rownames(out),
                 names(c(hs.binom$betas$unpenalized, hs.binom$betas$penalized)))
    expect_equal(out[, "mean"],
                 c(hs.binom$betas$unpenalized, hs.binom$betas$penalized))
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
