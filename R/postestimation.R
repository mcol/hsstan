##=============================================================================
##
## Copyright (c) 2019 Marco Colombo
##
## Parts of the code are based on https://avehtari.github.io/bayes_R2/bayes_R2.html
## Portions copyright (c) 2019 Aki Vehtari, Andrew Gelman, Ben Goodrich, Jonah Gabry
##
## Parts of the code are based on https://github.com/stan-dev/rstanarm
## Portions copyright (c) 2015, 2016, 2017 Trustees of Columbia University
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
#' @param object An object of class `hsstan`.
#' @param newdata Optional data frame to use to evaluate the log-likelihood.
#'        If `NULL` (default), the model matrix is used.
#' @param ... Currently ignored.
#'
#' @return
#' A matrix of size `S` by `N`, where `S` is number of draws from the posterior
#' distribution, and `N` is the number of data points.
#'
#' @examples
#' \dontshow{utils::example("hsstan", echo=FALSE)}
#' # continued from ?hsstan
#' log_lik(hs.biom)
#'
#' @importFrom rstantools log_lik
#' @importFrom stats rbinom rnorm
#' @method log_lik hsstan
#' @aliases log_lik
#' @export log_lik
#' @export
log_lik.hsstan <- function(object, newdata=NULL, ...) {
    if (is.null(newdata))
        newdata <- object$data
    y <- validate.outcome(newdata[[object$model.terms$outcome]])
    mu <- posterior_linpred(object, newdata, transform=TRUE)
    if (!is.logistic(object))
        sigma <- as.matrix(object$stanfit, pars="sigma")
    llkfun <- ifelse(is.logistic(object),
                     function(x) stats::dbinom(y[x], 1, mu[, x], log=TRUE),
                     function(x) stats::dnorm(y[x], mu[, x], sigma, log=TRUE))
    llk <- cbind(sapply(1:ncol(mu), llkfun))
    return(llk)
}

#' Posterior uncertainty intervals
#'
#' Compute posterior uncertainty intervals for `hsstan` objects.
#'
#' @param object An object of class `hsstan`.
#' @param pars Names of parameters for which posterior intervals should be
#'        returned, which can be specified as regular expressions. If `NULL`
#'        (default) then this refers to the set of predictors used in the model.
#' @param prob A value between 0 and 1 indicating the desired probability
#'        to be covered by the uncertainty intervals (0.95, by default).
#' @param ... Currently ignored.
#'
#' @return
#' A matrix with lower and upper interval bounds as columns and as many rows
#' as selected parameters.
#'
#' @examples
#' \dontshow{utils::example("hsstan", echo=FALSE)}
#' # continued from ?hsstan
#' posterior_interval(hs.biom)
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
#' @param object An object of class `hsstan`.
#' @param transform Whether the linear predictor should be transformed using
#'        the inverse-link function (`FALSE` by default).
#' @param newdata Optional data frame containing the variables to use to
#'        predict. If `NULL` (default), the model matrix is used. If
#'        specified, its continuous variables should be standardized, since
#'        the model coefficients are learnt on standardized data.
#' @param ... Currently ignored.
#'
#' @return
#' A matrix of size `S` by `N`, where `S` is the number of draws from the
#' posterior distribution of the (transformed) linear predictor, and `N` is
#' the number of data points.
#'
#' @examples
#' \dontshow{utils::example("hsstan", echo=FALSE)}
#' # continued from ?hsstan
#' posterior_linpred(hs.biom)
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
    linear.predictor <- multiplyABt(post.matrix, newdata)
    if (!transform)
        return(linear.predictor)
    return(object$family$linkinv(linear.predictor))
}

#' Posterior predictive distribution
#'
#' Draw from the posterior predictive distribution of the outcome.
#'
#' @param object An object of class `hsstan`.
#' @param newdata Optional data frame containing the variables to use to
#'        predict. If `NULL` (default), the model matrix is used. If
#'        specified, its continuous variables should be standardized, since
#'        the model coefficients are learnt on standardized data.
#' @param nsamples A positive integer indicating the number of posterior samples
#'        to use. If `NULL` (default) all samples are used.
#' @param seed Optional integer defining the seed for the pseudo-random number
#'        generator.
#' @param ... Currently ignored.
#'
#' @return
#' A matrix of size `S` by `N`, where `S` is the number of simulations from
#' the posterior predictive distribution, and `N` is the number of data points.
#'
#' @examples
#' \dontshow{utils::example("hsstan", echo=FALSE)}
#' # continued from ?hsstan
#' posterior_predict(hs.biom)
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

#' Posterior measures of performance
#'
#' Compute the log-likelihood and a relevant measure of performance (R-squared
#' or AUC) from the posterior samples.
#'
#' @param obj An object of class `hsstan` or `kfold`.
#' @param prob Width of the posterior interval (0.95, by default). It is
#'        ignored if `summary=FALSE`.
#' @param sub.idx Vector of indices of observations in the dataset to be used
#'        in computing the performance measures. If `NULL` (default), all
#'        observations in the dataset are used.
#' @param summary Whether a summary of the distribution of the performance
#'        measure should be returned rather than the pointwise values
#'        (`TRUE` by default).
#' @param cores Number of cores to use for parallelization (the value of
#'        `options("mc.cores")` by default).
#'
#' @return
#' The mean, standard deviation and posterior interval of the performance
#' measure (R-squared or AUC) if `summary=TRUE`, or a vector of values
#' of the performance measure with length equal to the size of the posterior
#' sample if `summary=FALSE`. Attribute `type` reports whether the performance
#' measures are cross-validated or not. If `sub.idx` is not `NULL`, attribute
#' `subset` reports the index of observations used in the computations.
#'
#' @examples
#' \donttest{
#' \dontshow{utils::example("hsstan", echo=FALSE)}
#' # continued from ?hsstan
#' posterior_performance(hs.biom, cores=1)
#' }
#'
#' @export
posterior_performance <- function(obj, prob=0.95, sub.idx=NULL, summary=TRUE,
                                  cores=getOption("mc.cores", 1)) {
    if (inherits(obj, "hsstan")) {
        obj <- list(fits=array(list(fit=obj, test.idx=1:nrow(obj$data)), c(1, 2)),
                    data=obj$data)
        colnames(obj$fits) <- c("fit", "test.idx")
    } else if (inherits(obj, c("kfold", "loo"))) {
        if (is.null(obj[["fits"]]))
            stop("No fitted models found, run 'kfold' with store.fits=TRUE.")
    } else
        stop("Not an 'hsstan' or 'kfold' object.")

    if (is.null(sub.idx)) {
        sub.idx <- 1:nrow(obj$data)
        used.subset <- FALSE
    } else {
        validate.indices(sub.idx, nrow(obj$data), "sub.idx")
        sub.idx <- sort(sub.idx)
        used.subset <- length(sub.idx) < nrow(obj$data)
    }

    validate.samples(obj$fits[[1]])
    validate.probability(prob)
    logistic <- is.logistic(obj$fits[[1]])
    num.folds <- nrow(obj$fits)

    ## loop over the folds
    y <- mu <- llk <- NULL
    for (fold in 1:num.folds) {
        hs <- obj$fits[[fold]]
        test.idx <- intersect(obj$fits[, "test.idx"][[fold]], sub.idx)
        if (length(test.idx) == 0)
            next
        newdata <- obj$data[test.idx, ]
        y <- c(y, obj$data[test.idx, hs$model.terms$outcome])
        mu <- cbind(mu, posterior_linpred(hs, newdata=newdata, transform=TRUE))
        llk <- cbind(llk, log_lik(hs, newdata=newdata))
    }

    if (used.subset && logistic && length(unique(y)) != 2)
        stop("'sub.idx' must contain both outcome classes.")

    if (logistic) {
        par.auc <- function(i)
            as.numeric(pROC::roc(y, mu[i, ], direction="<", quiet=TRUE)$auc)
        if (.Platform$OS.type != "windows") {
            out <- parallel::mclapply(X=1:nrow(mu), mc.cores=cores,
                                      mc.preschedule=TRUE, FUN=par.auc)
        } else { # windows
            cl <- parallel::makePSOCKcluster(cores)
            on.exit(parallel::stopCluster(cl))
            out <- parallel::parLapply(X=1:nrow(mu), cl=cl, fun=par.auc)
        }
    } else {
        out <- pmax(fastCor(y, mu), 0)^2
    }

    out <- cbind(perf=unlist(out), llk=rowSums(llk))
    colnames(out)[1] <- ifelse(logistic, "auc", "r2")
    if (summary)
        out <- posterior_summary(out, prob)
    attr(out, "type") <- paste0(if(num.folds == 1) "non ", "cross-validated")
    if (used.subset)
        attr(out, "subset") <- sub.idx
    return(out)
}

#' Predictive information criteria for Bayesian models
#'
#' Compute an efficient approximate leave-one-out cross-validation
#' using Pareto smoothed importance sampling (PSIS-LOO), or the widely
#' applicable information criterion (WAIC), also known as the Watanabe-Akaike
#' information criterion.
#'
#' @param x An object of class `hsstan`.
#' @param cores Number of cores used for parallelisation (the value of
#'        `options("mc.cores")` by default).
#' @param ... Currently ignored.
#'
#' @return
#' A `loo` object.
#'
#' @examples
#' \dontshow{utils::example("hsstan", echo=FALSE)}
#' \dontshow{oldopts <- options(mc.cores=2)}
#' # continued from ?hsstan
#' loo(hs.biom)
#' waic(hs.biom)
#' \dontshow{options(oldopts)}
#'
#' @importFrom loo loo
#' @method loo hsstan
#' @aliases loo
#' @export loo
#' @export
loo.hsstan <- function(x, cores=getOption("mc.cores"), ...) {
    validate.samples(x)
    llk <- log_lik(x)
    chain.id <- rep(1:ncol(x$stanfit), each=nrow(x$stanfit))
    r.eff <- loo::relative_eff(exp(llk), chain_id=chain.id, cores=cores)
    loo <- suppressWarnings(loo::loo(llk, r_eff=r.eff, cores=cores))
    return(loo)
}

#' @rdname loo.hsstan
#' @importFrom loo waic
#' @method waic hsstan
#' @aliases waic
#' @export waic
#' @export
waic.hsstan <- function(x, cores=getOption("mc.cores"), ...) {
    validate.samples(x)
    llk <- log_lik(x)
    waic <- suppressWarnings(loo::waic(llk, cores=cores))
    return(waic)
}

#' Bayesian and LOO-adjusted R-squared
#'
#' Compute the Bayesian and the LOO-adjusted R-squared from the posterior
#' samples. For Bayesian R-squared it uses the modelled residual variance
#' (rather than the variance of the posterior distribution of the residuals).
#' The LOO-adjusted R-squared uses Pareto smoothed importance sampling LOO
#' residuals and Bayesian bootstrap.
#'
#' @param object An object of class `hsstan`.
#' @param prob Width of the posterior interval (0.95, by default). It is
#'        ignored if `summary=FALSE`.
#' @param summary Whether a summary of the distribution of the R-squared
#'        should be returned rather than the pointwise values (`TRUE` by
#'        default).
#' @param ... Currently ignored.
#'
#' @return
#' The mean, standard deviation and posterior interval of R-squared if
#' `summary=TRUE`, or a vector of R-squared values with length equal to
#' the size of the posterior sample if `summary=FALSE`.
#'
#' @references
#' Andrew Gelman, Ben Goodrich, Jonah Gabry and Aki Vehtari (2019),
#' R-squared for Bayesian regression models,
#' _The American Statistician_, 73 (3), 307-309.
#' \url{https://doi.org/10.1080/00031305.2018.1549100}
#'
#' Aki Vehtari, Andrew Gelman, Ben Goodrich and Jonah Gabry (2019),
#' Bayesian R2 and LOO-R2.
#' \url{https://avehtari.github.io/bayes_R2/bayes_R2.html}
#'
#' @examples
#' \dontshow{utils::example("hsstan", echo=FALSE)}
#' \dontshow{oldopts <- options(mc.cores=2)}
#' # continued from ?hsstan
#' bayes_R2(hs.biom)
#' loo_R2(hs.biom)
#' \dontshow{options(oldopts)}
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
    var.mu <- apply(mu, 1, stats::var)
    if (is.logistic(object))
        sigma2 <- rowMeans(mu * (1 - mu))
    else
        sigma2 <- drop(as.matrix(object$stanfit, pars="sigma"))^2
    R2 <- var.mu / (var.mu + sigma2)
    if (summary)
        R2 <- vector.summary(R2, prob)
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
    y <- object$data[[object$model.terms$outcome]]
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

    if (summary)
        R2 <- vector.summary(R2, prob)
    return(R2)
}
