N <- 10
P <- 6
x <- data.frame(matrix(rnorm(N * P), ncol=P))
x[1, P - 1] <- NA
colnames(x)[c(P - 1, P)] <- c("X.NA", "y")
y.gauss <- rnorm(N)
y.binom <- rbinom(N, 1, 0.5)
y.factr <- factor(y.binom, labels=c("No", "Yes"))


test_that("valid objects",
{
    expect_error(validate.hsstan(x),
                 "Not an object of class 'hsstan'")
    expect_error(validate.samples(list(stanfit=NA)),
                 "No valid posterior samples stored in the 'hsstan' object")
})

test_that("newdata",
{
    expect_error(validate.newdata(list(data=x), y.gauss),
                 "'newdata' must be a data frame or a matrix")
    expect_error(validate.newdata(list(data=matrix(nrow=0, ncol=P)), NULL),
                 "'newdata' contains no rows or no columns")
    expect_error(validate.newdata(list(data=matrix(nrow=P, ncol=0)), NULL),
                 "'newdata' contains no rows or no columns")

    mt <- list(outcome="y", unpenalized="X1", penalized="X.NA")
    expect_error(validate.newdata(list(model.terms=mt), x),
                 "'newdata' contains missing values")

    mt <- list(outcome="y", unpenalized="X1", penalized="X2")
    expect_error(validate.newdata(list(model.terms=mt), x[, 3:4]),
                 "object 'X1' not found")
    expect_equal(validate.newdata(list(model.terms=mt), x),
                 model.matrix(reformulate(c("X1", "X2")), x))
})

test_that("model formula",
{
    expect_error(validate.model(NULL),
                 "argument is not a valid model")
    expect_error(validate.model(c("y", "X1")),
                 "Model formula specified incorrectly")
    expect_error(validate.model(~ 1),
                 "No outcome variable specified in the model")
    expect_error(validate.model(y ~ 0),
                 "Models with no intercept are not supported")
    expect_error(validate.model(y ~ X1 * X2),
                 "Interaction terms are not supported")
    expect_error(validate.model(y ~ X1 + X2, 1:5),
                 "'penalized' must be a character vector")

    model <- validate.model(y ~ X1 + X2, c("X3", "X4"))
    expect_type(model, "list")
    expect_named(model,
                 c("outcome", "unpenalized", "penalized"))
    expect_equal(model$outcome,
                 "y")
    expect_equal(model$unpenalized,
                 c("X1", "X2"))
    expect_equal(model$penalized,
                 c("X3", "X4"))

    model <- validate.model(y ~ X1 + X2 + X3, c("X2", "X4"))
    expect_equal(model$unpenalized,
                 c("X1", "X3"))
    expect_equal(model$penalized,
                 c("X2", "X4"))

    model <- validate.model(y ~ X1 + X2 + X3, c("X2", "X2"))
    expect_equal(model$unpenalized,
                 c("X1", "X3"))
    expect_equal(model$penalized,
                 "X2")

    model <- validate.model(y ~ x1 + x2, c())
    expect_equal(validate.model(y ~ x1 + x2, NULL),
                 model)
    expect_length(model$penalized, 0)

    model <- validate.model(y ~ 1, c())
    expect_equal(model$outcome,
                 "y")
    expect_length(model$unpenalized, 0)
    expect_length(model$penalized, 0)
})

test_that("model data",
{
    expect_error(validate.data(NULL, y ~ 1),
                 "'x' must be a data frame or a matrix")
    expect_error(validate.data(1:N, y ~ 1),
                 "'x' must be a data frame or a matrix")
    expect_equal(validate.data(as.matrix(x), validate.model(y ~ X1, "X3")),
                 x)
    expect_error(validate.variables(x, c()),
                 "No predictors present in the model")
    expect_error(validate.variables(x, c("X1", "zzz")),
                 "'zzz' not present in 'x'")
    expect_error(validate.variables(x, "X.NA"),
                 "Model variables contain missing values")
})

test_that("outcome variable",
{
    expect_error(validate.outcome("zzz"),
                 "Outcome variable of invalid type")
    expect_error(validate.outcome(as.factor(y.gauss)),
                 "A factor outcome variable can only have two levels")
    expect_equal(validate.outcome(y.factr),
                 y.binom)
    expect_equal(validate.outcome(y.factr == "Yes"),
                 y.binom)
})

test_that("family is valid",
{
    expect_error(validate.family(),
                 "Argument of 'family' is missing")
    expect_error(validate.family("zzz"),
                 "'zzz' is not a valid family")
    expect_error(validate.family(NULL),
                 "Argument of 'family' is not a valid family")
    expect_error(validate.family(x),
                 "Argument of 'family' is not a valid family")
    expect_error(validate.family(poisson),
                 "Only 'gaussian' and 'binomial' are supported families")
    expect_error(validate.family(binomial()),
                 "is missing, with no default")
})

test_that("invalid family inputs",
{
    expect_error(validate.family(binomial(), y.gauss),
                 "must contain two classes")
    expect_error(validate.family(binomial(), y.binom + 1),
                 "must contain 0-1 values")
    expect_error(validate.family(binomial(), y.binom - 1),
                 "must contain 0-1 values")
})

test_that("valid family inputs",
{
    expect_is(validate.family(gaussian()),
              "family")
    expect_equal(validate.family("binomial", y.binom)$family,
                 "binomial")
    expect_equal(validate.family("binomial", y.factr)$family,
                 "binomial")
    expect_equal(validate.family("gaussian")$family,
                 "gaussian")
    expect_equal(validate.family(gaussian())$family,
                 "gaussian")
    expect_equal(validate.family(gaussian)$family,
                 "gaussian")
})

test_that("validate.indices",
{
    expect_error(validate.indices(c(1:N, NA), N, "test"),
                 "'test' contains missing values")
    expect_error(validate.indices(letters, N, "test"),
                 "'test' must be an integer vector")
    expect_error(validate.indices(list(1:N), N, "test"),
                 "'test' must be an integer vector")
    expect_error(validate.indices(matrix(1:9, 3, 3), N, "test"),
                 "'test' must be an integer vector")
    expect_error(validate.indices(c(1:9, 10.5), N, "test"),
                 "'test' must be an integer vector")
    expect_error(validate.indices(1, N, "test"),
                 "'test' must contain at least two elements")
    expect_error(validate.indices(1:N, 5, "test"),
                 "'test' contains out of bounds indices")
    expect_error(validate.indices(0:N, N, "test"),
                 "'test' contains out of bounds indices")
    expect_error(validate.indices(rep(2, N), N, "test"),
                 "'test' contains duplicate indices")

    expect_silent(validate.indices(rep(2, N), N, "test", throw.duplicates=FALSE))
})

test_that("validate.folds",
{
    expect_error(validate.folds(c(1:N, NA), N),
                 "'folds' contains missing values")
    expect_error(validate.folds(list(1:N), N),
                 "'folds' must be an integer vector")
    expect_error(validate.folds(c(1:9, 10.5), N),
                 "'folds' must be an integer vector")
    expect_error(validate.folds(1:N, N + 1),
                 "'folds' should have length")
    expect_error(validate.folds(rep(2, N), N),
                 "'folds' must contain all indices up to")
    expect_error(validate.folds(sample(2:4, N, replace=TRUE), N),
                 "'folds' must contain all indices up to")

    folds <- rep(1:3, length.out=N)
    expect_equal(validate.folds(NULL, N),
                 rep(1, N))
    expect_equal(validate.folds(rep(1, N), N),
                 rep(1, N))
    expect_equal(validate.folds(folds, N),
                 folds)
    expect_type(validate.folds(folds * 1.0, N),
                "integer")
})

test_that("validate.start.from",
{
    expect_error(validate.start.from(hs.gauss, c("X1", NA)),
                 "'start.from' contains missing values")
    expect_error(validate.start.from(hs.gauss, 1:3),
                 "contains '1', '2', '3', which cannot be matched")
    expect_error(validate.start.from(hs.gauss, "a"),
                 "contains 'a', which cannot be matched")
    expect_error(validate.start.from(hs.gauss, "Intercept"),
                 "contains 'Intercept', which cannot be matched")
    expect_error(validate.start.from(hs.gauss, "X1b"),
                 "contains 'X1b', which cannot be matched")
    expect_equal(validate.start.from(hs.gauss, character(0)),
                 1)

    idx <- validate.start.from(hs.gauss, c("X2", "X1"))
    expect_equal(idx,
                 validate.start.from(hs.gauss, c("X1", "X2")))
    expect_equal(idx,
                 validate.start.from(hs.gauss, c("X2", "X1", "X1", "X2")))

    idx <- validate.start.from(hs.gauss, NULL)
    expect_equal(idx,
                 seq_along(hs.gauss$betas$unpenalized))
    expect_type(idx,
                "integer")
})

test_that("validate.adapt.delta",
{
    expect_error(validate.adapt.delta("a"),
                 "must be a single numerical value")
    expect_error(validate.adapt.delta(1:3),
                 "must be a single numerical value")
    expect_error(validate.adapt.delta(0.5),
                 "must be at least 0.8")
    expect_error(validate.adapt.delta(1.5),
                 "must be less than 1")
})

test_that("validate.probability",
{
    expect_error(validate.probability(NULL),
                 "'prob' must be a single value between 0 and 1")
    expect_error(validate.probability("a"),
                 "'prob' must be a single value between 0 and 1")
    expect_error(validate.probability(0),
                 "'prob' must be a single value between 0 and 1")
    expect_error(validate.probability(1),
                 "'prob' must be a single value between 0 and 1")
    expect_error(validate.probability(c(0.2, 0.8)),
                 "'prob' must be a single value between 0 and 1")
})

test_that("validate.rstan.args",
{
    expect_error(validate.rstan.args(iters=200),
                 "Argument 'iters' not recognized")
})

test_that("get.pars",
{
    expect_error(get.pars(hs.gauss, NA),
                 "'pars' must be a character vector")
    expect_error(get.pars(hs.gauss, "zzz"),
                 "No pattern in 'pars' matches parameter names")
    expect_equal(get.pars(hs.gauss, NULL),
                 c("beta_u", "beta_p"))
    expect_equal(get.pars(hs.gauss, "X1"),
                 c("X1b", "X1c", "X10"))
    expect_equal(get.pars(hs.gauss, c("X[345]", "9$")),
                 c("X3", "X4", "X5", "X9"))
})

test_that("expand.terms",
{
    expect_error(expand.terms(df, ""),
                 "subscript out of bounds")
    expect_error(expand.terms(df, "zzz"),
                 "object 'zzz' not found")
    expect_equal(expand.terms(df, NULL),
                 character(0))
    expect_equal(expand.terms(df, c()),
                 character(0))
    expect_equal(expand.terms(df, character(0)),
                 character(0))
    expect_equal(expand.terms(df, "X1"),
                 c("(Intercept)", "X1b", "X1c"))
})
