## silence output and warnings
SW <- function(expr) capture.output(suppressWarnings(expr))

set.seed(1)
N <- 50
P <- 10
U <- 3
x <- matrix(rnorm(N * P), nrow=N, ncol=P)
b <- runif(P) - 0.5
y.gauss <- rnorm(N, mean=x %*% b, sd=runif(1, 1, 2))
y.binom <- rbinom(N, 1, binomial()$linkinv(x %*% b))
df <- data.frame(x, y.gauss=y.gauss, y.binom=y.binom)
unp <- paste0("X", 1:U)
pen <- setdiff(paste0("X", 1:P), unp)

mod.gauss <- reformulate(unp, "y.gauss")
mod.binom <- reformulate(unp, "y.binom")

SW({
    hs.gauss <- hsstan(df, mod.gauss, pen, iter=500, chains=2, family=gaussian)
    hs.binom <- hsstan(df, mod.binom, pen, iter=500, chains=2, family=binomial)
})

test_that("posterior_interval",
{
    expect_error(posterior_interval(hs.gauss, pars=1:3),
                 "'pars' must be a character vector")
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
