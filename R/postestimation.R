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


#' Pointwise log-likelihood
#'
#' Compute the pointwise log-likelihood.
#'
#' @param object An object of class \code{hsstan}.
#' @param newdata Optional data frame to use to evaluate the log-likelihood.
#'        If \code{NULL} (default), the model matrix is used.
#' @param ... Currently ignored.
#'
#' @return
#' A matrix of size \code{S} by \code{N}, where \code{S} is number of draws
#' from the posterior distribution, and \code{N} is the number of data points.
#'
#' @importFrom rstantools log_lik
#' @importFrom stats rbinom rnorm
#' @method log_lik hsstan
#' @aliases log_lik
#' @export log_lik
#' @export
log_lik.hsstan <- function(object, newdata=NULL, ...) {
    if (is.null(newdata))
        newdata <- object$data[object$in.train, ]
    y <- validate.outcome(newdata[[object$model.terms$outcome]])
    mu <- posterior_linpred(object, newdata, transform=TRUE)
    if (!is.logistic(object))
        sigma <- as.matrix(object$stanfit, pars="sigma")
    llkfun <- ifelse(is.logistic(object),
                     function(x) dbinom(y[x], 1, mu[, x], log=TRUE),
                     function(x) dnorm(y[x], mu[, x], sigma, log=TRUE))
    llk <- cbind(sapply(1:ncol(mu), llkfun))
    return(llk)
}

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
    validate.probability(prob)
    pars <- get.pars(object, pars)
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
#' @param newdata Optional data frame containing the variables to use to
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
#' @param newdata Optional data frame containing the variables to use to
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
#' Compute an efficient approximate leave-one-out cross-validation
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

#' Bayesian and LOO-adjusted R-squared
#'
#' Compute the Bayesian and the LOO-adjusted R-squared from the posterior
#' samples. For Bayesian R-squared it uses the modelled residual variance
#' (rather than the variance of the posterior distribution of the residuals).
#' The LOO-adjusted R-squared uses Pareto smoothed importance sampling LOO
#' residuals and Bayesian bootstrap.
#'
#' @param object An object of class \code{hsstan}.
#' @param prob Width of the posterior interval (0.95, by default). It is
#'        ignored if \code{summary=FALSE}.
#' @param summary Whether a summary of the distribution of the R-squared
#'        should be returned rather than the pointwise values (\code{TRUE} by
#'        default).
#' @param ... Currently ignored.
#'
#' @return
#' The mean, standard deviation and posterior interval of R-squared if
#' \code{summary=TRUE}, or a vector of R-squared values with length equal to
#' the size of the posterior sample if \code{summary=FALSE}.
#'
#' @references
#' Andrew Gelman, Ben Goodrich, Jonah Gabry and Aki Vehtari (2019),
#' R-squared for Bayesian regression models, \emph{The American Statistician},
#' \url{https://doi.org/10.1080/00031305.2018.1549100}.
#'
#' Aki Vehtari, Andrew Gelman, Ben Goodrich and Jonah Gabry (2019),
#' Bayesian R2 and LOO-R2,
#' \url{https://avehtari.github.io/bayes_R2/bayes_R2.html}.
#'
#' @importFrom rstantools bayes_R2
#' @method bayes_R2 hsstan
#' @aliases bayes_R2
#' @export bayes_R2
#' @export
bayes_R2.hsstan <- function(object, prob=0.95, summary=TRUE, ...) {
    validate.samples(object)
    validate.probability(prob)
    mu <- posterior_linpred(object, transform=TRUE)
    var.mu <- apply(mu, 1, var)
    if (is.logistic(object))
        sigma2 <- rowMeans(mu * (1 - mu))
    else
        sigma2 <- drop(as.matrix(object$stanfit, pars="sigma"))^2
    R2 <- var.mu / (var.mu + sigma2)
    if (summary) {
        low <- (1 - prob) / 2
        upp <- 1 - low
        R2 <- c(mean=mean(R2), sd=stats::sd(R2),
                stats::quantile(R2, c(low, upp)))
    }
    return(R2)
}

#' @rdname bayes_R2.hsstan
#' @importFrom rstantools loo_R2
#' @method loo_R2 hsstan
#' @aliases loo_R2
#' @export loo_R2
#' @export
loo_R2.hsstan <- function(object, prob=0.95, summary=TRUE, ...) {
    validate.samples(object)
    validate.probability(prob)
    y <- object$data[object$in.train, object$model.terms$outcome]
    mu <- posterior_linpred(object, transform=TRUE)
    ll <- log_lik(object)
    S <- nrow(mu)
    N <- ncol(mu)
    chains <- object$stanfit@sim$chains
    r.eff <- loo::relative_eff(exp(ll), chain_id=rep(1:chains, each=S / chains))
    psis <- suppressWarnings(loo::psis(-ll, r_eff=r.eff))
    mu.loo <- loo::E_loo(mu, psis, log_ratios=-ll)$value
    err.loo <- mu.loo - y

    ## set the random seed as the seed used in the first chain and ensure
    ## the old RNG state is restored on exit
    rng.state.old <- .Random.seed
    on.exit(assign(".Random.seed", rng.state.old, envir=.GlobalEnv))
    set.seed(object$stanfit@stan_args[[1]]$seed)

    ## dirichlet weights for bayesian bootstrap
    exp.draws <- matrix(stats::rexp(S * N, rate=1), nrow=S, ncol=N)
    dw <- exp.draws / rowSums(exp.draws)

    var.y <- (rowSums(sweep(dw, 2, y^2, FUN = "*")) -
              rowSums(sweep(dw, 2, y, FUN = "*"))^2) * (N / (N - 1))
    var.err.loo <- (rowSums(sweep(dw, 2, err.loo^2, FUN = "*")) -
                    rowSums(sweep(dw, 2, err.loo, FUN = "*")^2)) * (N / (N - 1))

    R2 <- 1 - var.err.loo / var.y
    R2[R2 < -1] <- -1
    R2[R2 > 1] <- 1

    if (summary) {
        low <- (1 - prob) / 2
        upp <- 1 - low
        R2 <- c(mean=mean(R2), sd=stats::sd(R2),
                stats::quantile(R2, c(low, upp)))
    }
    return(R2)
}
