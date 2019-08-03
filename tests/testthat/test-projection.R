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
num.sel <- 5

SW({
    hs.nopen <- hsstan(df, mod.gauss, iter=100, chains=2, family=gaussian)
    hs.gauss <- hsstan(df, mod.gauss, pen, iter=300, chains=2, family=gaussian)
    hs.binom <- hsstan(df, mod.binom, pen, iter=300, chains=2, family=binomial)

    sel.gauss <- projsel(hs.gauss, out.csv="out.csv")
    sel.binom <- projsel(hs.binom, max.num.pred=num.sel)
})

test_that("projsel",
{
    expect_error(projsel(hs.nopen),
                 "Model doesn't contain penalized predictors")

    expect_s3_class(sel.gauss,
                    "projsel")
    expect_equal(colnames(sel.gauss),
                 c("var", "kl", "elpd", "delta.elpd"))
    expect_equal(nrow(sel.gauss),
                 length(pen) + 1)
    expect_equal(as.character(sel.gauss$var[1]),
                 "Initial set of covariates")
    expect_true(all(diff(sel.gauss$kl) < 0))
    expect_true(file.exists("out.csv"))
    unlink("out.csv")

    expect_equal(nrow(sel.binom),
                 num.sel + 1)
    expect_true(all(diff(sel.binom$kl) < 0))
})
