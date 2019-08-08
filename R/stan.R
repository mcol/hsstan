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
#' Run the No-U-Turn Sampler (NUTS) as implemented in Stan over all the data
#' or within the cross-validation folds.
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
#' @param chains Number of Markov chains to run (4 by default).
#' @param seed Optional integer defining the seed for the pseudo-random number
#'        generator.
#' @param qr Whether the QR decomposition should be used to decorrelate the
#'        predictors (\code{TRUE} by default). This is silently set to
#'        \code{FALSE} if there are more predictors than observations.
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
#' @importFrom rstan stan_model
#' @importFrom stats gaussian model.matrix reformulate
#' @export
hsstan <- function(x, covs.model, penalized=NULL, family=gaussian, folds=NULL,
                   iter=ifelse(is.null(folds), 2000, 1000), warmup=iter / 2,
                   scale.u=2, regularized=TRUE, nu=ifelse(regularized, 1, 3),
                   par.ratio=0.05, global.df=1, slab.scale=2, slab.df=4,
                   chains=4, seed=123, qr=TRUE, adapt.delta=NULL) {

    model.terms <- validate.model(covs.model, penalized)
    x <- validate.data(x, model.terms)
    y <- x[[model.terms$outcome]]
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
    N <- nrow(X)
    P <- ncol(X)
    U <- P - length(penalized)
    num.folds <- max(length(folds), 1)

    ## thin QR decomposition
    if (P > N) qr <- FALSE
    if (qr) {
        qr.dec <- qr(X)
        Q.qr <- qr.Q(qr.dec)
        R.inv <- qr.solve(qr.dec, Q.qr) * sqrt(N - 1)
        Q.qr <- Q.qr * sqrt(N - 1)
    }

    ## whether it's a proper cross-validation
    is.cross.validation <- num.folds > 1

    ## loop over the cross-validation folds
    cv <- parallel::mclapply(X=1:num.folds, mc.preschedule=FALSE,
                             FUN=function(fold) {

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

        X_train <- if (qr) Q.qr[train, ] else X[train, ]
        y_train <- y[train]
        N_train <- nrow(X_train)
        X_test <- if (qr) Q.qr[test, ] else X[test, ]
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

        if (qr) {
            pars <- grep("beta_", samples@sim$pars_oi, value=TRUE)
            stopifnot(pars[1] == "beta_u")
            beta.tilde <- rstan::extract(samples, pars=pars,
                                         inc_warmup=TRUE, permuted=FALSE)
            B <- apply(beta.tilde, 1:2, FUN=function(z) R.inv %*% z)
            for (chain in 1:chains) {
                for (p in 1:P)
                    samples@sim$samples[[chain]][[par.idx[p]]] <- B[p, , chain]
            }
        }

        betas <- get.coefficients(samples, colnames(X))
        obj <- list(stanfit=samples, betas=betas,
                    model.terms=model.terms, family=family, qr=qr,
                    data=x, in.train=train, in.test=test)
        class(obj) <- "hsstan"
        return(obj)
    })

    return(if (is.cross.validation) cv else cv[[1]])
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
