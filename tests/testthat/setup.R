## silence output and warnings
SW <- function(expr) suppressMessages(suppressWarnings(expr))

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
df[, 1] <- factor(letters[rbinom(N, 2, 0.5) + 1])
df$X1b_X3 <- df$X3 * (df$X1 == "b")
df$X1c_X3 <- df$X3 * (df$X1 == "c")
df$X3_X2 <- df$X3 * df$X2

## model options
unp <- paste0("X", 1:U)
pen <- setdiff(paste0("X", 1:P), unp)
mod.gauss <- reformulate(unp, "y.gauss")
mod.binom <- reformulate(unp, "y.binom")
folds <- c(rep(1, N / 2), rep(2, N / 2))
iters <- 500
chains <- 2

## numerical tolerance
tol <- 0.000001

## wrapper to set commonly used options
hs <- function(model, family, ...)
    hsstan(df, model, pen, iter=iters, chains=chains, family=family,
           refresh=0, ...)

message("Running hsstan models...")
SW({
    hs.gauss <- hs(mod.gauss, gaussian)
    hs.binom <- hs(mod.binom, binomial)
    hs.inter <- hs(y.gauss ~ X1 * X3 + X2 * X3, gaussian)
})

message("Running cross-validated hsstan models...")
SW({
    cv.gauss <- kfold(hs.gauss, folds=folds, chains=2)
    cv.binom <- kfold(hs.binom, folds=folds)
    cv.nofit <- cv.gauss
    cv.nofit$fits <- cv.nofit$data <- NULL
})
