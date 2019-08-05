library(hsstan)

## silence output and warnings
SW <- function(expr) capture.output(suppressWarnings(expr))

library(doParallel)
registerDoParallel(cores=1)

set.seed(1)
N <- 50
P <- 10
U <- 3
x <- matrix(rnorm(N * P), nrow=N, ncol=P)
b <- runif(P) - 0.5
y.gauss <- rnorm(N, mean=x %*% b, sd=runif(1, 1, 2))
y.binom <- rbinom(N, 1, binomial()$linkinv(x %*% b))
df <- data.frame(x, y.gauss=y.gauss, y.binom=y.binom)
folds <- list(1:25, 25:N)
unp <- paste0("X", 1:U)
pen <- setdiff(paste0("X", 1:P), unp)
mod.gauss <- reformulate(unp, "y.gauss")
mod.binom <- reformulate(unp, "y.binom")

SW({
    hs.gauss <- hsstan(df, mod.gauss, pen, iter=500, chains=2, family=gaussian)
    hs.binom <- hsstan(df, mod.binom, pen, iter=500, chains=2, family=binomial)
    cv.gauss <- hsstan(df, mod.gauss, pen, iter=500, chains=2, family=gaussian,
                       folds=folds)
    cv.binom <- hsstan(df, mod.binom, pen, iter=500, chains=2, family=binomial,
                       folds=folds)
})

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
        ss.gauss <- sample.stan(df, df$y.gauss, unp, pen, iter=500, chains=2)
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
        sv.binom <- sample.stan.cv(df, df$y.binom, unp, pen, iter=500, chains=2,
                                   family=binomial, folds=folds)
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

test_that("get.cv.performance",
{
    tol <- 0.000001
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

    out <- get.cv.performance(cv.gauss, out.csv="out.csv")
    expect_equal(nrow(out), length(folds) + 1)
    expect_equal(out$set,
                 c(paste("Fold", 1:length(folds)), "Overall"))
    expect_equal(out$test.llk,
                 c(-67.06889, -63.05863, -130.12752), tolerance=tol)
    expect_equal(out$r2,
                 c(0.016581920, 0.175190218, 0.003522558), tolerance=tol)
    expect_true(file.exists("out.csv"))
    unlink("out.csv")

    out <- get.cv.performance(cv.binom)
    expect_equal(out$test.llk,
                 c(-28.85088, -27.02137, -55.87226), tolerance=tol)
    expect_equal(out$auc,
                 c(0.4935065, 0.4821429, 0.4692308), tolerance=tol)
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

test_that("nsamples",
{
    expect_equal(nsamples(hs.gauss), 500)
})