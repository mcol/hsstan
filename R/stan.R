##=============================================================================
##
## Copyright (c) 2017 Marco Colombo and Paul McKeigue
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##=============================================================================


#' Returns posterior means and 95\% credible intervals
#'
#' @param samples An object of class \code{stanfit}.
#' @param pars Vector of variables names to be extracted.
#' @param varnames Vector of variable names.
#'
#' @importFrom rstan As.mcmc.list
summarize.params <- function(samples, pars, varnames) {
    post.params <- As.mcmc.list(object=samples, pars)
    params.summary <- summary(post.params)
    params.summary <- data.frame(mean=params.summary$statistics[, 1],
                                 params.summary$quantiles[, c(1, 5)])
    colnames(params.summary) <- c("mean", "centile2.5", "centile97.5")
    rownames(params.summary) <- varnames
    return(params.summary)
}

#' Returns the posterior means for the specified variables
#'
#' @param samples An object of class \code{stanfit}.
#' @param varnames Vector of variable names to be extracted.
#'
#' @importMethodsFrom rstan extract
posterior.means <- function(samples, varnames) {
    unlist(lapply(extract(samples, pars=varnames),
                  function(z) if (is.matrix(z)) colMeans(z) else mean(z)))
}

#' Returns the posterior means of the regression coefficients
#'
#' @param samples An object of class \code{stanfit}.
#' @param coeff.names Vector of names for the coefficients.
get.coefficients <- function(samples, coeff.names) {
    beta.u <- posterior.means(samples, "beta_u")
    beta.p <- posterior.means(samples, "beta_p")
    stopifnot(length(c(beta.u, beta.p)) == length(coeff.names))
    u.idx <- 1:length(beta.u)
    names(beta.u) <- coeff.names[u.idx]
    names(beta.p) <- coeff.names[-u.idx]

    return(list(unpenalized=beta.u, penalized=beta.p))
}

#' Returns the posterior means of the unpenalized regression coefficients
#'
#' @param samples An object of class \code{stanfit}.
#' @param coeff.names Vector of names for the coefficients.
get.unpenalized.coefficients <- function(samples, coeff.names) {
    beta.u <- posterior.means(samples, "beta_u")
    stopifnot(length(beta.u) == length(coeff.names))
    names(beta.u) <- coeff.names

    return(list(unpenalized=beta.u))
}

#' Runs a cross-validated Stan model
#'
#' Runs either with Hamiltonian Monte Carlo or variational Bayes over the
#' cross-validation folds.
#'
#' @param stan.file Path to a Stan model file.
#' @param x Data.frame of predictors.
#' @param y Vector of outcomes. For a logistic regression model, this is
#'        expected to contain only \code{0-1} entries.
#' @param covariates Names of the variables to be used as (unpenalized)
#'        predictors.
#' @param biomarkers Names of the variables to be used as penalized predictors.
#'        If it is specified as an empty vector, a model with only unpenalized
#'        covariates is used.
#' @param folds List of cross-validation folds, where each element contains
#'        the indices of the test observations.
#' @param standardize Whether the design matrix should be standardized
#'        (\code{TRUE} by default).
#' @param store.samples Whether the posterior samples should be saved
#'        (\code{TRUE} by default).
#' @param nu Number of degrees of freedom for the half-Student-t priors.
#' @param model.type Either \code{"mc"} for Hamiltonian Monte Carlo, or
#'        \code{"vb"} for variational Bayes.
#'
#' @importFrom rstan stan_model
#' @importMethodsFrom rstan sampling vb
#' @export
sample.stan.cv <- function(stan.file, x, y, covariates, biomarkers, folds,
                           standardize=TRUE, store.samples=TRUE,
                           nu=3, model.type=c("mc", "vb")) {

    stopifnot(nrow(x) == length(y))
    stopifnot(all(covariates %in% colnames(x)))
    stopifnot(all(biomarkers %in% colnames(x)))
    model.type <- match.arg(model.type)

    ## create the design matrix
    model <- paste("~", paste(c(covariates, biomarkers), collapse=" + "))
    X <- model.matrix(as.formula(model), data=x)
    P <- ncol(X)
    U <- P - length(biomarkers)
    which.unpenalized <- 1:U
    which.penalized <- setdiff(1:P, which.unpenalized)
    X <- X[, c(which.unpenalized, which.penalized)]
    N <- nrow(X)
    num.folds <- length(folds)

    ## names of variables to be standardized
    stand.cols <- colnames(x)[!sapply(x, class) %in% c("factor", "character")]

    ## compile the model
    model <- stan_model(file=stan.file)

    ## cross-validation
    cv <- foreach(fold=1:num.folds) %dopar% {
        test <- 1:N %in% folds[[fold]]
        train <- !test
        X_train <- X[train, ]
        y_train <- y[train]
        N_train <- nrow(X_train)
        X_test <- X[test, ]
        y_test <- y[test]
        N_test <- nrow(X_test)

        ## standardize continuous outcome
        if (length(unique(y_train)) > 2) {
            train.mean <- mean(y_train, na.rm=TRUE)
            train.sd <- sd(y_train, na.rm=TRUE)
            y_train <- (y_train - train.mean) / train.sd
            y_test <- (y_test - train.mean) / train.sd
        }

        ## standardize all columns corresponding to numerical variables: this
        ## excludes those generated from factor/character variables and the
        ## intercept column
        if (standardize) {
            train.mean <- colMeans(X_train, na.rm=TRUE)
            train.sd <- apply(X_train, 2, sd, na.rm=TRUE)
            train.mean[!colnames(X_train) %in% stand.cols] <- 0
            train.sd[!colnames(X_train) %in% stand.cols] <- 1
            X_train <- ((X_train - rep(train.mean, each=N_train))
                         %*% diag(1 / train.sd))
            X_test <- ((X_test - rep(train.mean, each=N_test))
                        %*% diag(1 / train.sd))
        }

        ## parameters not used by a model are ignored
        data.input <- list(N_train=N_train, N_test=N_test,
                           y_train=y_train, y_test=y_test,
                           X_train=X_train, X_test=X_test,
                           P=P, U=U, nu=nu)

        if (model.type == "mc") {
            samples <- sampling(model, data=data.input,
                                chains=4, iter=1000, warmup=500,
                                seed=123, control=list(adapt_delta=0.8))
        }
        else {
            samples <- vb(model, data=data.input,
                          iter=50000, output_samples=2000,
                          init="random", algorithm="meanfield")
        }

        if (P == U) {
            betas <- get.unpenalized.coefficients(samples, colnames(X))
        } else {
            betas <- get.coefficients(samples, colnames(X))
        }

        ## linear predictor of test data and residual standard deviation
        y_pred <- posterior.means(samples, "y_pred")
        sigma <- tryCatch(posterior.means(samples, "sigma"),
                          error=function(e) return(1))

        if (!store.samples) samples <- NA
        list(samples=samples, betas=betas,
             y_pred=y_pred, sigma=sigma,
             X_train=X_train, X_test=X_test,
             y_train=y_train, y_test=y_test, train=train, test=test)
    }

    return(cv)
}

#' Runs a Stan model on all available data.
#'
#' @param stan.file Path to a Stan model file.
#' @param x Data.frame of predictors.
#' @param y Vector of outcomes. For a logistic regression model, this is
#'        expected to contain only \code{0-1} entries.
#' @param covariates Names of the variables to be used as (unpenalized)
#'        predictors.
#' @param biomarkers Names of the variables to be used as penalized predictors.
#'        If it is specified as an empty vector, a model with only unpenalized
#'        covariates is used.
#' @param standardize Whether the design matrix should be standardized.
#' @param nu Number of degrees of freedom for the half-Student-t priors.
#' @param model.type Either \code{"mc"} for Hamiltonian Monte Carlo, or
#'        \code{"vb"} for variational Bayes.
#'
#' @importFrom rstan stan_model
#' @importMethodsFrom rstan sampling vb
#' @export
sample.stan <- function(stan.file, x, y, covariates, biomarkers,
                        standardize=TRUE,
                        nu=3, model.type=c("mc", "vb")) {

    stopifnot(nrow(x) == length(y))
    stopifnot(all(covariates %in% colnames(x)))
    stopifnot(all(biomarkers %in% colnames(x)))
    model.type <- match.arg(model.type)

    ## create the design matrix
    model <- paste("~", paste(c(covariates, biomarkers), collapse=" + "))
    X <- model.matrix(as.formula(model), data=x)
    N <- nrow(X)
    P <- ncol(X)
    U <- P - length(biomarkers)
    which.unpenalized <- 1:U
    which.penalized <- setdiff(1:P, which.unpenalized)
    X <- X[, c(which.unpenalized, which.penalized)]

    ## standardize continuous outcome
    if (length(unique(y)) > 2) {
        y <- as.numeric(scale(y))
    }

    ## standardize all columns corresponding to numerical variables: this
    ## excludes those generated from factor/character variables as well as
    ## the intercept column
    if (standardize) {
        stand.cols <- colnames(x)[!sapply(x, class) %in% c("factor", "character")]
        stand.idx <- which(colnames(X) %in% stand.cols)
        X[, stand.idx] <- scale(X[, stand.idx])
    }

    ## use all available data for both training and testing: this effectively
    ## computes the fit of the model (y_pred) for all observations
    train <- test <- rep(TRUE, N)

    ## parameters not used by a model are ignored
    data.input <- list(N_train=N, N_test=N,
                       y_train=y, y_test=y,
                       X_train=X, X_test=X,
                       P=P, U=U, nu=nu)

    ## compile the model
    cat("Compiling STAN model", stan.file, "...")
    model <- stan_model(file=stan.file)
    if (model.type == "mc") {
        samples <- sampling(model, data=data.input, iter=1000, warmup=500,
                            chains=4, seed=123, control=list(adapt_delta=0.8))
    }
    else {
        samples <- vb(model, data=data.input, iter=50000, output_samples=2000,
                      init="random", algorithm="meanfield")
    }

    if (P == U) {
        betas <- get.unpenalized.coefficients(samples, colnames(X))
    } else {
        betas <- get.coefficients(samples, colnames(X))
    }

    ## linear predictor of test data and residual standard deviation
    y_pred <- posterior.means(samples, "y_pred")
    sigma <- tryCatch(posterior.means(samples, "sigma"),
                      error=function(e) return(1))

    return(list(samples=samples, betas=betas, train=train, test=test,
                y_pred=y_pred, sigma=sigma,
                data=X, y=y))
}

#' Extract measures of performance from the cross-validation results
#'
#' @param hs.cv Cross-validated Stan model.
#' @param out.csv Optional name of a file where the output can be saved in
#'        CSV format.
#'
#' @importFrom pROC roc
#' @export
get.cv.performance <- function(hs.cv, out.csv=NULL) {

    gaussian.llk <- function(y.pred, y.obs, disp)
        -0.5 * sum(log(disp) + (y.pred - y.obs)^2 / disp)
    binomial.llk <- function(y.pred, y.obs)
        sum(log(y.pred[y.obs == 1])) + sum(log(1 - y.pred[y.obs == 0]))
    r2 <- function(y.pred, y.obs) {
        corr <- cor(y.pred, y.obs)
        if (corr < 0)
            return(0)
        return(corr^2)
    }
    to.prob <- function(lin.pred) 1 / (1 + exp(-lin.pred))
    auc <- function(y.pred, y.obs) as.numeric(roc(y.obs, y.pred)$auc)
    loglik.ratio <- function(y.pred, y.obs, prop.cases)
        (2 * y.obs - 1) * (log(y.pred) - log(1 - y.pred)) -
        (2 * y.obs - 1) * (log(prop.cases) - log(1 - prop.cases))

    y.obs.all <- y.pred.hs.all <- NULL
    sigma.hs.all <- NULL
    num.folds <- length(hs.cv)
    llk <- perf <- rep(NA, num.folds)
    llk.ratio <- llk.ratio.var <- rep(NA, num.folds)

    ## loop over the folds
    for (fold in 1:num.folds) {
        y.obs <- hs.cv[[fold]]$y_test

        y.pred.hs <- hs.cv[[fold]]$y_pred
        sigma.hs <- hs.cv[[fold]]$sigma
        is.logistic <- length(sigma.hs) == 1 && sigma.hs == 1

        ## logistic regression
        if (is.logistic) {

            ## proportion of cases in the training fold (prior probability)
            prop.cases <- sum(hs.cv[[fold]]$y_train) / sum(hs.cv[[fold]]$train)

            y.pred.hs <- to.prob(y.pred.hs)
            llk[fold] <- binomial.llk(y.pred.hs, y.obs)
            perf[fold] <- auc(y.pred.hs, y.obs)
            llkr <- loglik.ratio(y.pred.hs, y.obs, prop.cases)
            llk.ratio[fold] <- mean(llkr)
            llk.ratio.var[fold] <- var(llkr)
        }

        ## linear regression
        else {
            llk[fold] <- gaussian.llk(y.pred.hs, y.obs, sigma.hs)
            perf[fold] <- r2(y.pred.hs, y.obs)
        }

        sigma.hs.all <- c(sigma.hs.all, sigma.hs)
        y.obs.all <- c(y.obs.all, y.obs)
        y.pred.hs.all <- c(y.pred.hs.all, y.pred.hs)
    }
    set <- paste("Fold", 1:num.folds)

    ## compute log-likelihood and performance measure of the full model
    ## on the full vector of withdrawn observations
    set <- c(set, "Overall")
    llk <- c(llk, ifelse(is.logistic,
                         binomial.llk(y.pred.hs.all, y.obs.all),
                         gaussian.llk(y.pred.hs.all, y.obs.all,
                                      mean(sigma.hs.all))))
    perf <- c(perf, ifelse(is.logistic,
                           auc(y.pred.hs.all, y.obs.all),
                           r2(y.pred.hs.all, y.obs.all)))

    res <- data.frame(set=set, test.llk=llk, perf=perf)
    colnames(res)[3] <- gsub("perf", ifelse(is.logistic, "auc", "r2"),
                             colnames(res)[3])
    if (is.logistic) {
        prop.cases <- sum(y.obs.all == 1) / length(y.obs.all)
        llkr <- loglik.ratio(y.pred.hs.all, y.obs.all, prop.cases)
        res$llk.ratio <- c(llk.ratio, mean(llkr))
        res$llk.ratio.var <- c(llk.ratio.var, var(llkr))
    }

    if (!is.null(out.csv))
        write.csv(file=out.csv, res, row.names=FALSE)
    return(res)
}