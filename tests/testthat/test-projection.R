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
unp <- paste0("X", 1:U)
pen <- setdiff(paste0("X", 1:P), unp)

mod.gauss <- reformulate(unp, "y.gauss")
mod.binom <- reformulate(unp, "y.binom")
num.sel <- 5

test_that("projsel",
{
    SW({
        hs.nopen <- hsstan(df, mod.gauss, iter=50, chains=1, family=gaussian)
    })
    expect_error(projsel(hs.nopen),
                 "Model doesn't contain penalized predictors")
})

test_that("projsel for gaussian family",
{
    SW({
        hs.gauss <- hsstan(df, mod.gauss, pen, iter=300, chains=2, family=gaussian)
        sel.gauss <- projsel(hs.gauss, out.csv="out.csv")
    })

    tol <- 0.000001
    expect_s3_class(sel.gauss,
                    "projsel")
    expect_equal(colnames(sel.gauss),
                 c("var", "kl", "elpd", "delta.elpd"))
    expect_equal(nrow(sel.gauss),
                 length(pen) + 1)
    expect_equal(sel.gauss$var[1],
                 "Initial set of covariates")
    expect_equal(sel.gauss$var[-1],
                 paste0("X", c(9, 4, 8, 6, 10, 7, 5)))
    expect_equal(sel.gauss$kl[1],
                 0.0416673097, tolerance=tol)
    expect_equal(sel.gauss$kl[length(pen) + 1],
                 0)
    expect_equal(sel.gauss$elpd[1],
                 -111.93956, tolerance=tol)
    expect_equal(sel.gauss$elpd[length(pen) + 1],
                 -109.46688, tolerance=tol)
    expect_equal(sel.gauss$delta.elpd[1],
                 -2.47267928, tolerance=tol)
    expect_true(all(diff(sel.gauss$kl) < 0))
    expect_true(file.exists("out.csv"))
    unlink("out.csv")
})

test_that("projsel for binomial family",
{
    SW({
        hs.binom <- hsstan(df, mod.binom, pen, iter=300, chains=2, family=binomial)
        sel.binom <- projsel(hs.binom, max.num.pred=num.sel)
    })

    tol <- 0.000001
    expect_equal(nrow(sel.binom),
                 num.sel + 1)
    expect_equal(sel.binom$var[-1],
                 paste0("X", c(6, 9, 5, 10, 8)))
    expect_equal(sel.binom$kl[1],
                 0.024903468, tolerance=tol)
    expect_equal(sel.binom$kl[num.sel + 1],
                 0.0016025, tolerance=tol)
    expect_equal(sel.binom$elpd[1],
                 -34.07116, tolerance=tol)
    expect_equal(sel.binom$elpd[num.sel + 1],
                 -32.87313, tolerance=tol)
    expect_equal(sel.binom$delta.elpd[1],
                 -1.172055517, tolerance=tol)
    expect_true(all(diff(sel.binom$kl) < 0))
})

test_that("projsel for a cross-validated object",
{
    SW({
        folds <- list(1:25, 26:N)
        cv.gauss <- hsstan(df, mod.gauss, pen, iter=200, chains=2,
                           folds=folds, family=gaussian)
        sel.gauss <- projsel(cv.gauss[[1]])
    })

    expect_equal(sel.gauss$elpd[length(pen) + 1],
                 sum(cv.gauss[[1]]$loglik))
})
