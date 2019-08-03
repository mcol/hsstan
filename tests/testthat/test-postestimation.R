## silence output and warnings
SW <- function(expr) capture.output(suppressWarnings(expr))

set.seed(1)
N <- 50
P <- 10
x <- matrix(rnorm(N * P), nrow=N, ncol=P)
b <- runif(P) - 0.5
y.gauss <- rnorm(N, mean=x %*% b, sd=runif(1, 1, 2))
y.binom <- rbinom(N, 1, binomial()$linkinv(x %*% b))
df <- data.frame(x, y.gauss=y.gauss, y.binom=y.binom)

mod.gauss <- reformulate(paste0("X", 1:P), "y.gauss")
mod.binom <- reformulate(paste0("X", 1:P), "y.binom")

SW({
    hs.gauss <- hsstan(df, mod.gauss, family=gaussian, iter=500, chains=2)
    hs.binom <- hsstan(df, mod.binom, family=binomial, iter=500, chains=2)
})

test_that("posterior_linpred",
{
    expect_equal(posterior_linpred(hs.gauss),
                 posterior_linpred(hs.gauss, newdata=df))
})

test_that("posterior_predict",
{
    expect_equal(posterior_predict(hs.gauss, seed=1),
                 posterior_predict(hs.gauss, seed=1, newdata=df))
})
