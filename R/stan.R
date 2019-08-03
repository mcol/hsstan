##=============================================================================
##
## Copyright (c) 2017-2019 Marco Colombo and Paul McKeigue
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


#' Return the posterior means of the regression coefficients
#'
#' @param samples An object of class \code{stanfit}.
#' @param coeff.names Vector of names for the coefficients.
#'
#' @return
#' A list with two named elements: \var{unpenalized} for the posterior means
#' of the unpenalized covariates, and \var{penalized} for the posterior means
#' of the penalized predictors (which can be \code{NULL} for baseline models).
get.coefficients <- function(samples, coeff.names) {
    beta.u <- colMeans(as.matrix(samples, pars="beta_u"))
    beta.p <- tryCatch(colMeans(as.matrix(samples, pars="beta_p")),
                       error=function(e) return(NULL))
    stopifnot(length(c(beta.u, beta.p)) == length(coeff.names))
    u.idx <- 1:length(beta.u)
    names(beta.u) <- coeff.names[u.idx]
    if (!is.null(beta.p))
        names(beta.p) <- coeff.names[-u.idx]

    return(list(unpenalized=beta.u, penalized=beta.p))
}

#' Fit a plain or cross-validated Stan model
#'
#' Runs either with Hamiltonian Monte Carlo or variational Bayes over the
#' cross-validation folds (if specified) or over all the data.
#'
#' @param x Data frame containing outcome, covariates and penalized predictors.
#'        Continuous predictors and outcome variable should be standardized
#'        before fitting the models as priors assume them to have mean zero and
#'        have the same scale.
#' @param covs.model Formula containing the unpenalized covariates.
#' @param penalized Names of the variables to be used as penalized predictors.
#'        If \code{NULL} or an empty vector, a model with only unpenalized
#'        covariates is fitted.
#' @param family Type of model fitted: either \code{gaussian()} for linear
#'        regression (default) or \code{binomial()} for logistic regression.
#' @param folds List of cross-validation folds, where each element contains
#'        the indices of the test observations. If \code{NULL} (default), no
#'        cross-validation is performed.
#' @param store.samples Whether the posterior samples should be saved
#'        (by default, \code{FALSE} for cross-validation and \code{TRUE}
#'        otherwise).
#' @param chains Number of Markov chains to run (4 by default).
#' @param seed Integer defining the seed for the pseudo-random number generator.
#' @param adapt.delta Target average proposal acceptance probability for
#'        adaptation, a value between 0.8 and 1 (excluded). If unspecified,
#'        it's set to 0.99 for hierarchical shrinkage models and to 0.95 for
#'        base models.
#' @param iter Total number of iterations in each chain, including warmup
#'        (by default, 1000 iterations for cross-validation and 2000 otherwise).
#' @param warmup Number of warmup iterations per chain (by default, half the
#'        total number of iterations).
#' @param scale.u Prior scale (standard deviation) for the unpenalised
#'        covariates.
#' @param regularized If \code{TRUE} (default), the regularized horseshoe prior
#'        is used as opposed to the original horseshoe prior.
#' @param nu Number of degrees of freedom of the half-Student-t prior on the
#'        local shrinkage parameters (by default, 1 if \code{regularized=TRUE}
#'        and 3 otherwise).
#' @param par.ratio Expected ratio of non-zero to zero coefficients (ignored
#'        if \code{regularized=FALSE}). The scale of the global shrinkage
#'        parameter corresponds to \code{par.ratio} divided by the square root
#'        of the number of observations; for linear regression only, it's further
#'        multiplied by the residual standard deviation \code{sigma}.
#' @param global.df Number of degrees of freedom for the global shrinkage
#'        parameter (ignored if \code{regularized=FALSE}).
#' @param slab.scale Scale of the regularization parameter (ignored if
#'        \code{regularized=FALSE}).
#' @param slab.df Number of degrees of freedom of the regularization parameter
#'        (ignored if \code{regularized=FALSE}).
#'
#' @importFrom foreach %dopar%
#' @importFrom rstan stan_model
#' @importFrom stats gaussian model.matrix reformulate
#' @export
hsstan <- function(x, covs.model, penalized=NULL, family=gaussian, folds=NULL,
                   iter=ifelse(is.null(folds), 2000, 1000), warmup=iter / 2,
                   scale.u=2, regularized=TRUE, nu=ifelse(regularized, 1, 3),
                   par.ratio=0.05, global.df=1, slab.scale=2, slab.df=4,
                   store.samples=is.null(folds), chains=4, seed=123,
                   adapt.delta=NULL) {

    model.terms <- validate.model(covs.model, penalized)
    x <- validate.data(x, model.terms)
    y <- validate.outcome(x[[model.terms$outcome]])
    family <- validate.family(family, y)
    regularized <- as.integer(regularized)
    if (chains < 1)
        stop("'chains' must be a positive integer.")

    ## choose the model to be fitted
    model <- ifelse(length(penalized) == 0, "base", "hs")
    if (family$family == "binomial") model <- paste0(model, "_logit")

    ## set or check adapt.delta
    if (is.null(adapt.delta)) {
        adapt.delta <- ifelse(grepl("hs", model), 0.99, 0.95)
    } else {
        validate.adapt.delta(adapt.delta)
    }

    ## create the design matrix
    unpenalized <- model.terms$unpenalized
    X <- model.matrix(reformulate(c(unpenalized, penalized)), data=x)
    P <- ncol(X)
    U <- P - length(penalized)
    which.unpenalized <- 1:U
    which.penalized <- setdiff(1:P, which.unpenalized)
    X <- X[, c(which.unpenalized, which.penalized)]
    N <- nrow(X)
    num.folds <- max(length(folds), 1)

    ## whether it's a proper cross-validation
    is.cross.validation <- num.folds > 1

    ## loop over the cross-validation folds
    fold <- NULL   # silence a note raised by R CMD check
    cv <- foreach(fold=1:num.folds) %dopar% {

        ## create a proper training/test split
        if (is.cross.validation) {
            test <- 1:N %in% folds[[fold]]
            train <- !test
        }

        ## use all available data for both training and testing: this
        ## effectively computes the fit of the model (y_pred) for all
        ## observations
        else {
            train <- test <- rep(TRUE, N)
        }

        X_train <- X[train, ]
        y_train <- y[train]
        N_train <- nrow(X_train)
        X_test <- X[test, ]
        y_test <- y[test]
        N_test <- nrow(X_test)

        ## global scale for regularized horseshoe prior
        global.scale <- if (regularized) par.ratio / sqrt(N_train) else 1

        ## parameters not used by a model are ignored
        data.input <- list(N_train=N_train, N_test=N_test,
                           y_train=y_train, y_test=y_test,
                           X_train=X_train, X_test=X_test,
                           P=P, U=U, scale_u=scale.u,
                           regularized=regularized, nu=nu,
                           global_scale=global.scale, global_df=global.df,
                           slab_scale=slab.scale, slab_df=slab.df)

        ## run the stan model
        samples <- rstan::sampling(stanmodels[[model]], data=data.input,
                                   iter=iter, warmup=warmup,
                                   chains=chains, seed=seed,
                                   control=list(adapt_delta=adapt.delta))

        ## assign proper names
        par.idx <- grep("^beta_[up]", names(samples))
        stopifnot(length(par.idx) == ncol(X))
        names(samples)[par.idx] <- colnames(X)

        ## linear predictor of test data, regression coefficients and
        ## residual standard deviation
        post.matrix <- as.matrix(samples)
        y_pred <- colMeans(tcrossprod(post.matrix[, par.idx], X_test))
        fitted <- family$linkinv(y_pred)
        betas <- get.coefficients(samples, colnames(X))
        coefs <- c(betas$unpenalized, betas$penalized)

        ## test log-likelihood
        loglik <- colMeans(post.matrix[, grep("log_lik", colnames(post.matrix))])

        if (!store.samples) samples <- NA
        obj <- list(stanfit=samples, betas=betas, coefficients=coefs,
                    model.terms=model.terms,
                    linear.predictors=y_pred, fitted.values=fitted, family=family,
                    loglik=loglik, data=X_train, y=y_train)
        if (is.cross.validation)
            obj <- c(obj, list(withdrawn.data=X_test, y_test=y_test,
                               train=train, test=test))
        class(obj) <- "hsstan"
        return(obj)
    }

    return(if(is.cross.validation) cv else cv[[1]])
}

#' Deprecated functions to fit hierarchical shrinkage models
#'
#' @inheritParams hsstan
#' @param x Data frame of predictors.
#' @param y Vector of outcomes. For a logistic regression model, this is
#'        expected to contain only \code{0-1} entries.
#' @param covariates Names of the variables to be used as (unpenalized)
#'        predictors.
#' @param biomarkers Names of the variables to be used as penalized predictors.
#'        If it is specified as an empty vector, a model with only unpenalized
#'        covariates is used.
#' @param logit \code{FALSE} for linear regression (default), \code{TRUE} for
#'        logistic regression.
#' @param ... Further options passed to \code{hsstan}.
#'
#' @section Note:
#' \code{sample.stan} and \code{sample.stan.cv} are thin wrappers around
#' \code{hsstan} provided for backward compatibility. They are considered
#' deprecated.
#'
#' @export
sample.stan <- function(x, y, covariates, biomarkers=NULL,
                        logit=FALSE, ...) {
    default.args <- list(x=cbind(x, hsstan_y_=y), penalized=biomarkers,
                         covs.model=reformulate(covariates, "hsstan_y_"),
                         folds=NULL, family=ifelse(logit, binomial, gaussian),
                         iter=2000)
    input.args <- list(...)
    default.args[names(input.args)] <- input.args
    input.args[names(default.args)] <- NULL
    do.call(hsstan, c(default.args, input.args))
}

#' @rdname sample.stan
#' @export
sample.stan.cv <- function(x, y, covariates, biomarkers=NULL, folds,
                           logit=FALSE, iter=1000, ...) {
    sample.stan(x, y, covariates, biomarkers, logit, folds=folds, iter=iter, ...)
}

#' Extract measures of performance from the cross-validation results
#'
#' @param hs.cv Cross-validated Stan model.
#' @param out.csv Optional name of a file where the output can be saved in
#'        CSV format.
#'
#' @importFrom pROC roc
#' @importFrom stats cor var
#' @importFrom utils write.csv
#' @export
get.cv.performance <- function(hs.cv, out.csv=NULL) {

    r2 <- function(y.pred, y.obs) {
        corr <- cor(y.pred, y.obs)
        if (corr < 0)
            return(0)
        return(corr^2)
    }
    auc <- function(y.pred, y.obs)
        as.numeric(roc(y.obs, y.pred, direction="<", quiet=TRUE)$auc)
    loglik.ratio <- function(y.pred, y.obs, prop.cases)
        (2 * y.obs - 1) * (log(y.pred) - log(1 - y.pred)) -
        (2 * y.obs - 1) * (log(prop.cases) - log(1 - prop.cases))

    if (inherits(hs.cv, "hsstan")) {
        hs.cv <- list(hs.cv)
        num.folds <- 1
    } else {
        num.folds <- length(hs.cv)
    }
    y.obs.all <- y.pred.hs.all <- NULL
    llk <- perf <- rep(NA, num.folds)
    llk.ratio <- llk.ratio.var <- rep(NA, num.folds)
    is.logistic <- is.logistic(hs.cv[[1]])
    perf.fun <- ifelse(is.logistic, auc, r2)

    ## loop over the folds
    for (fold in 1:num.folds) {
        y.obs <- hs.cv[[fold]]$y_test
        if (is.null(y.obs))
            y.obs <- hs.cv[[fold]]$y

        y.pred.hs <- hs.cv[[fold]]$fitted.values
        perf[fold] <- perf.fun(y.pred.hs, y.obs)

        ## log-likelihood ratio
        if (is.logistic) {
            ## proportion of cases in the training fold (prior probability)
            prop.cases <- sum(hs.cv[[fold]]$y) / length(hs.cv[[fold]]$y)
            llkr <- loglik.ratio(y.pred.hs, y.obs, prop.cases)
            llk.ratio[fold] <- mean(llkr)
            llk.ratio.var[fold] <- var(llkr)
        }

        llk[fold] <- sum(hs.cv[[fold]]$loglik)
        y.obs.all <- c(y.obs.all, y.obs)
        y.pred.hs.all <- c(y.pred.hs.all, y.pred.hs)
    }
    set <- paste("Fold", 1:num.folds)

    ## compute log-likelihood and performance measure of the full model
    ## on the full vector of withdrawn observations
    set <- c(set, "Overall")
    llk <- c(llk, sum(sapply(hs.cv, function(z) sum(z$loglik))))
    perf <- c(perf, perf.fun(y.pred.hs.all, y.obs.all))

    res <- data.frame(set=set, test.llk=llk, perf=perf)
    colnames(res)[3] <- gsub("perf", ifelse(is.logistic, "auc", "r2"),
                             colnames(res)[3])
    if (is.logistic) {
        prop.cases <- sum(y.obs.all == 1) / length(y.obs.all)
        llkr <- loglik.ratio(y.pred.hs.all, y.obs.all, prop.cases)
        res$llk.ratio <- c(llk.ratio, mean(llkr))
        res$llk.ratio.var <- c(llk.ratio.var, var(llkr))
    }

    if (num.folds == 1) {
        res <- res[1, ]
        if (is.null(hs.cv[[fold]]$y_test))
            res$set <- "Non cross-validated"
    }

    if (!is.null(out.csv))
        write.csv(file=out.csv, res, row.names=FALSE)
    return(res)
}

#' Approximate leave-one-out cross-validation
#'
#' This computes an efficient approximate leave-one-out cross-validation
#' using Pareto smoothed importance sampling (PSIS-LOO).
#'
#' @param x An object of class \code{hsstan}.
#' @param save.psis Whether intermediate results should be saved in the
#'        returned object (\code{FALSE} by default).
#' @param cores Number of cores used for parallelisation (the value of
#'        the \code{mc.cores} option by default).
#'
#' @return
#' A \code{loo} object.
#'
#' @importMethodsFrom rstan loo
#' @method loo hsstan
#' @aliases loo
#' @export loo
#' @export
loo.hsstan <- function(x, save.psis=FALSE, cores=getOption("mc.cores")) {
    validate.samples(x)
    suppressWarnings(rstan::loo(x$stanfit, save_psis=save.psis, cores=cores))
}
