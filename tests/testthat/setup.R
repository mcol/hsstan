library(doParallel)
registerDoParallel(cores=1)

## silence output and warnings
SW <- function(expr) capture.output(suppressWarnings(expr))

## dataset
set.seed(1)
N <- 50
P <- 10
U <- 3
x <- matrix(rnorm(N * P), nrow=N, ncol=P)
b <- runif(P) - 0.5
y.gauss <- rnorm(N, mean=x %*% b, sd=runif(1, 1, 2))
y.binom <- factor(rbinom(N, 1, binomial()$linkinv(x %*% b)))
df <- data.frame(x, y.gauss=y.gauss, y.binom=y.binom)

## model options
unp <- paste0("X", 1:U)
pen <- setdiff(paste0("X", 1:P), unp)
mod.gauss <- reformulate(unp, "y.gauss")
mod.binom <- reformulate(unp, "y.binom")
folds <- list(1:25, 26:N)
iters <- 500
chains <- 2

## numerical tolerance
tol <- 0.000001

## wrapper to set commonly used options
hs <- function(model, family, ...)
    hsstan(df, model, pen, iter=iters, chains=chains, family=family, ...)

message("Running hsstan models...")
SW({
    hs.gauss <- hs(mod.gauss, gaussian)
    hs.binom <- hs(mod.binom, binomial)
    cv.gauss <- hs(mod.gauss, gaussian, folds=folds)
    cv.binom <- hs(mod.binom, binomial, folds=folds)
})
