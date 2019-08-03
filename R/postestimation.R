##=============================================================================
##
## Copyright (c) 2019 Marco Colombo
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


#' Posterior uncertainty intervals
#'
#' Compute posterior uncertainty intervals for \code{hsstan} objects.
#'
#' @param object An object of class \code{hsstan}.
#' @param pars Names of parameters for which posterior intervals should be
#'        returned, which can be specified as regular expressions. If \code{NULL}
#'        (default) then this refers to the set of predictors used in the model.
#' @param prob A value between 0 and 1 indicating the desired probability
#'        to be covered by the uncertainty intervals (0.95, by default).
#' @param ... Currently ignored.
#'
#' @return
#' A matrix with lower and upper interval bounds as columns and as many rows
#' as selected parameters.
#'
#' @importFrom rstantools posterior_interval
#' @method posterior_interval hsstan
#' @aliases posterior_interval
#' @export posterior_interval
#' @export
posterior_interval.hsstan <- function(object, pars=NULL, prob=0.95, ...) {
    validate.samples(object)
    if (is.null(pars))
        pars <- grep("^beta_", object$stanfit@model_pars, value=TRUE)
    else {
        if (!is.character(pars))
            stop("'pars' must be a character vector.")
        get.pars <- function(x) grep(x, object$stanfit@sim$fnames_oi, value=TRUE)
        pars <- unlist(sapply(pars, get.pars))
    }
    post.matrix <- as.matrix(object$stanfit, pars=pars)
    rstantools::posterior_interval(post.matrix, prob=prob)
}

#' Posterior distribution of the linear predictor
#'
#' Extract the posterior draws of the linear predictor, possibly transformed
#' by the inverse-link function.
#'
#' @param object An object of class \code{hsstan}.
#' @param transform Whether the linear predictor should be transformed using
#'        the inverse-link function (\code{FALSE} by default).
#' @param newdata An optional data frame containing the variables used to
#'        predict. If \code{NULL} (default), the model matrix is used. If
#'        specified, its continuous variables should be standardized, since
#'        the model coefficients are learnt on standardized data.
#' @param ... Currently ignored.
#'
#' @return
#' A matrix of size \code{S} by \code{N}, where \code{S} is the number of draws
#' from the posterior distribution of the (transformed) linear predictor, and
#' \code{N} is the number of data points.
#'
#' @importFrom rstantools posterior_linpred
#' @method posterior_linpred hsstan
#' @aliases posterior_linpred
#' @export posterior_linpred
#' @export
posterior_linpred.hsstan <- function(object, transform=FALSE,
                                     newdata=NULL, ...) {
    validate.samples(object)
    newdata <- validate.newdata(object, newdata)
    pars <- grep("^beta_", object$stanfit@model_pars, value=TRUE)
    post.matrix <- as.matrix(object$stanfit, pars=pars)
    linear.predictor <- tcrossprod(post.matrix, newdata)
    if (!transform)
        return(linear.predictor)
    return(object$family$linkinv(linear.predictor))
}

#' Posterior predictive distribution
#'
#' Draw from the posterior predictive distribution of the outcome.
#'
#' @param object An object of class \code{hsstan}.
#' @param newdata An optional data frame containing the variables used to
#'        predict. If \code{NULL} (default), the model matrix is used. If
#'        specified, its continuous variables should be standardized, since
#'        the model coefficients are learnt on standardized data.
#' @param nsamples A positive integer indicating the number of posterior samples
#'        to use. If \code{NULL} (default) all samples are used.
#' @param seed Optional integer defining the seed for the pseudo-random number
#'        generator.
#' @param ... Currently ignored.
#'
#' @return
#' A matrix of size \code{S} by \code{N}, where \code{S} is the number of
#' simulations from the posterior predictive distribution, and \code{N} is the
#' number of data points.
#'
#' @importFrom rstantools posterior_predict
#' @importFrom stats rbinom rnorm
#' @method posterior_predict hsstan
#' @aliases posterior_predict
#' @export posterior_predict
#' @export
posterior_predict.hsstan <- function(object, newdata=NULL, nsamples=NULL,
                                     seed=NULL, ...) {
    validate.samples(object)

    ## extract a random subset of the posterior samples
    if (!is.null(seed))
        set.seed(seed)
    num.samples <- nsamples(object)
    samp <- sample(num.samples, min(num.samples, nsamples))
    nsamples <- length(samp)
    if (nsamples == 0)
        stop("'nsamples' must be a positive integer.")

    ## generate the posterior predictions
    mu <- posterior_linpred(object, newdata, transform=TRUE)[samp, , drop=FALSE]
    nobs <- ncol(mu)
    if (is.logistic(object))
        pp <- t(sapply(1:nsamples, function(z) rbinom(nobs, 1, mu[z, ])))
    else {
        sigma <- as.matrix(object$stanfit, pars="sigma")[samp, , drop=FALSE]
        pp <- t(sapply(1:nsamples, function(z) rnorm(nobs, mu[z, ], sigma[z])))
    }
    return(pp)
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
