##=============================================================================
##
## Copyright (c) 2017-2019 Marco Colombo and Paul McKeigue
##
## kfold.hsstan() is based on code from https://github.com/stan-dev/rstanarm
## Portions copyright (C) 2015, 2016, 2017 Trustees of Columbia University
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


#' Hierarchical shrinkage models
#'
#' Run the No-U-Turn Sampler (NUTS) as implemented in Stan to fit a hierarchical
#' shrinkage model.
#'
#' @param x Data frame containing outcome, covariates and penalized predictors.
#'        Continuous predictors and outcome variable should be standardized
#'        before fitting the models as priors assume them to have mean zero and
#'        unit variance.
#' @param covs.model Formula containing the unpenalized covariates.
#' @param penalized Names of the variables to be used as penalized predictors.
#'        Any variable that is already part of the `covs.model` formula will be
#'        penalized. If `NULL` or an empty vector, a model with only unpenalized
#'        covariates is fitted.
#' @param family Type of model fitted: either `gaussian()` for linear regression
#'        (default) or `binomial()` for logistic regression.
#' @param seed Optional integer defining the seed for the pseudo-random number
#'        generator.
#' @param qr Whether the thin QR decomposition should be used to decorrelate the
#'        predictors (`TRUE` by default). This is silently set to `FALSE` if
#'        there are more predictors than observations.
#' @param adapt.delta Target average proposal acceptance probability for
#'        adaptation, a value between 0.8 and 1 (excluded). If unspecified,
#'        it's set to 0.99 for hierarchical shrinkage models and to 0.95 for
#'        base models.
#' @param iter Total number of iterations in each chain, including warmup
#'        (2000 by default).
#' @param warmup Number of warmup iterations per chain (by default, half the
#'        total number of iterations).
#' @param scale.u Prior scale (standard deviation) for the unpenalized
#'        covariates.
#' @param regularized If `TRUE` (default), the regularized horseshoe prior
#'        is used as opposed to the original horseshoe prior.
#' @param nu Number of degrees of freedom of the half-Student-t prior on the
#'        local shrinkage parameters (by default, 1 if `regularized=TRUE`
#'        and 3 otherwise).
#' @param par.ratio Expected ratio of non-zero to zero coefficients (ignored
#'        if `regularized=FALSE`). The scale of the global shrinkage parameter
#'        corresponds to `par.ratio` divided by the square root of the number of
#'        observations; for linear regression only, it's further multiplied by
#'        the residual standard deviation `sigma`.
#' @param global.df Number of degrees of freedom for the global shrinkage
#'        parameter (ignored if `regularized=FALSE`). Larger values induce more
#'        shrinkage.
#' @param slab.scale Scale of the regularization parameter (ignored if
#'        `regularized=FALSE`).
#' @param slab.df Number of degrees of freedom of the regularization parameter
#'        (ignored if `regularized=FALSE`).
#' @param keep.hs.pars Whether the parameters for the horseshoe prior should be
#'        kept in the `stanfit` object returned (`FALSE` by default).
#' @param ... Further arguments passed to [rstan::sampling()],
#'        such as `chains` (4 by default), `cores` (the value of
#'        `options("mc.cores")` by default), `refresh` (`iter / 10` by default).
#'
#' @return
#' An object of class `hsstan` containing the following fields:
#' \item{stanfit}{an object of class `stanfit` containing the output
#'       produced by Stan, including posterior samples and diagnostic summaries.
#'       It can be manipulated using methods from the **rstan** package.}
#' \item{betas}{posterior means of the unpenalized and penalized regression
#'       parameters.}
#' \item{call}{the matched call.}
#' \item{data}{the dataset used in fitting the model.}
#' \item{model.terms}{a list of names for the outcome variable, the unpenalized
#'       covariates and the penalized predictors.}
#' \item{family}{the `family` object used.}
#' \item{hsstan.settings}{the optional settings used in the model.}
#'
#' @seealso
#' [kfold()] for cross-validating a fitted object.
#'
#' @examples
#' \dontshow{oldopts <- options(mc.cores=2)}
#' data(diabetes)
#'
#' # non-default settings for speed of the example
#' df <- diabetes[1:50, ]
#' hs.biom <- hsstan(df, Y ~ age + sex, penalized=colnames(df)[5:10],
#'                   chains=2, iter=250)
#' \dontshow{options(oldopts)}
#'
#' @importFrom stats gaussian
#' @export
hsstan <- function(x, covs.model, penalized=NULL, family=gaussian,
                   iter=2000, warmup=floor(iter / 2),
                   scale.u=2, regularized=TRUE, nu=ifelse(regularized, 1, 3),
                   par.ratio=0.05, global.df=1, slab.scale=2, slab.df=4,
                   qr=TRUE, seed=123, adapt.delta=NULL,
                   keep.hs.pars=FALSE, ...) {

    model.terms <- validate.model(covs.model, penalized)
    x <- validate.data(x, model.terms)
    y <- x[[model.terms$outcome]]
    family <- validate.family(family, y)
    regularized <- as.integer(regularized)
    validate.positive.scalar(iter, "iter", int=TRUE)
    validate.positive.scalar(warmup, "warmup", int=TRUE)
    if (iter <= warmup)
        stop("'warmup' must be smaller than 'iter'.", call.=FALSE)

    ## stop if options to be passed to rstan::sampling are not valid, as
    ## to work around rstan issue #681
    validate.rstan.args(...)

    ## parameter names not to include by default in the stanfit object
    hs.pars <- c("lambda", "tau", "z", "c2")
    if (keep.hs.pars)
        hs.pars <- NA

    ## retrieve the call and its actual argument values
    call <- match.call(expand.dots=TRUE)
    args <- c(as.list(environment()), list(...))
    for (nm in names(call)[-c(1:2)]) # exclude "" and "x"
        call[[nm]] <- args[[nm]]

    ## choose the model to be fitted
    model <- ifelse(length(penalized) == 0, "base0", "hs")
    if (family$family == "binomial") model <- paste0(model, "_logit")

    ## set or check adapt.delta
    if (is.null(adapt.delta)) {
        adapt.delta <- ifelse(grepl("hs", model), 0.99, 0.95)
    } else {
        validate.adapt.delta(adapt.delta)
    }

    ## create the design matrix
    X <- ordered.model.matrix(x, model.terms$unpenalized, model.terms$penalized)
    N <- nrow(X)
    P <- ncol(X)

    ## number of penalized and unpenalized columns in the design matrix
    K <- length(expand.terms(x, model.terms$penalized)[-1])
    U <- P - K

    ## thin QR decomposition
    if (P > N) qr <- FALSE
    if (qr) {
        qr.dec <- qr(X)
        Q.qr <- qr.Q(qr.dec)
        R.inv <- qr.solve(qr.dec, Q.qr) * sqrt(N - 1)
        Q.qr <- Q.qr * sqrt(N - 1)
    }

    ## global scale for regularized horseshoe prior
    global.scale <- if (regularized) par.ratio / sqrt(N) else 1

    ## core block
    {
        ## parameters not used by a model are ignored
        data.input <- list(X=if (qr) Q.qr else X, y=y, N=N,
                           P=P, U=U, scale_u=scale.u,
                           regularized=regularized, nu=nu,
                           global_scale=global.scale, global_df=global.df,
                           slab_scale=slab.scale, slab_df=slab.df)

        ## run the stan model
        samples <- rstan::sampling(stanmodels[[model]], data=data.input,
                                   iter=iter, warmup=warmup, seed=seed, ...,
                                   pars=hs.pars, include=keep.hs.pars,
                                   control=list(adapt_delta=adapt.delta))
        if (is.na(nrow(samples)))
            stop("rstan::sampling failed, see error message above.", call.=FALSE)

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
            chains <- ncol(beta.tilde)
            for (chain in 1:chains) {
                for (p in 1:P)
                    samples@sim$samples[[chain]][[par.idx[p]]] <- B[p, , chain]
            }
        }

        ## store the hierarchical shrinkage settings
        opts <- list(adapt.delta=adapt.delta, qr=qr, seed=seed, scale.u=scale.u)
        if (K > 0)
            opts <- c(opts, regularized=regularized, nu=nu, par.ratio=par.ratio,
                      global.scale=global.scale, global.df=global.df,
                      slab.scale=slab.scale, slab.df=slab.df)

        ## compute the posterior means of the regression coefficients
        betas <- list(unpenalized=colMeans(as.matrix(samples, pars="beta_u")),
                      penalized=tryCatch(colMeans(as.matrix(samples,
                                                            pars="beta_p")),
                                         error=function(e) NULL))
        obj <- list(stanfit=samples, betas=betas, call=call, data=x,
                    model.terms=model.terms, family=family, hsstan.settings=opts)
        class(obj) <- "hsstan"
    }

    return(obj)
}

#' K-fold cross-validation
#'
#' Perform K-fold cross-validation using the same settings used when fitting
#' the model on the whole data.
#'
#' @param x An object of class `hsstan`.
#' @param folds Integer vector with one element per observation indicating the
#'        cross-validation fold in which the observation should be withdrawn.
#' @param chains Number of Markov chains to run. By default this is set to 1,
#'        independently of the number of chains used for `x`.
#' @param store.fits Whether the fitted models for each fold should be stored
#'        in the returned object (`TRUE` by default).
#' @param cores Number of cores to use for parallelization (the value of
#'        `options("mc.cores")` by default). The cross-validation folds will
#'        be distributed to the available cores, and the Markov chains for each
#'        model will be run sequentially.
#' @param ... Further arguments passed to [rstan::sampling()].
#'
#' @return
#' An object with classes `kfold` and `loo` that has a similar structure as the
#' objects returned by [loo()] and [waic()] and is compatible with the
#' \code{\link[loo]{loo_compare}} function for
#' comparing models. The object contains the following fields:
#' \item{estimates}{a matrix containing point estimates and standard errors of
#'       the expected log pointwise predictive density ("elpd_kfold"),
#'       the effective number of parameters ("p_kfold", always `NA`) and the
#'       K-fold information criterion "kfoldic" (which is `-2 * elpd_kfold`,
#'       i.e., converted to the deviance scale).}
#' \item{pointwise}{a matrix containing the pointwise contributions of
#'       "elpd_kfold", "p_kfold" and "kfoldic".}
#' \item{fits}{a matrix with two columns and number of rows equal to the number
#'       of cross-validation folds. Column `fit` contains the fitted
#'       `hsstan` objects for each fold, and column `test.idx` contains
#'       the indices of the withdrawn observations for each fold. This is not
#'       present if `store.fits=FALSE`.}
#' \item{data}{the dataset used in fitting the model (before withdrawing
#'       observations). This is not present if `store.fits=FALSE`.}
#'
#' @examples
#' \donttest{
#' \dontshow{utils::example("hsstan", echo=FALSE)}
#' # continued from ?hsstan
#' # only 2 folds for speed of example
#' folds <- rep(1:2, length.out=length(df$Y))
#' cv.biom <- kfold(hs.biom, folds=folds, cores=2)
#' }
#'
#' @importFrom loo kfold
#' @method kfold hsstan
#' @aliases kfold
#' @export kfold
#' @export
kfold.hsstan <- function(x, folds, chains=1, store.fits=TRUE,
                         cores=getOption("mc.cores", 1), ...) {
    data <- x$data
    N <- nrow(data)
    folds <- validate.folds(folds, N)
    num.folds <- max(folds)
    validate.positive.scalar(chains, "chains", int=TRUE)
    validate.rstan.args(...)

    ## collect the list of calls to be evaluated in parallel
    calls <- list()
    for (fold in 1:num.folds) {
        test.idx <- which(folds == fold)
        fit.call <- stats::update(object=x, x=data[-test.idx, , drop=FALSE],
                                  chains=chains, cores=1, refresh=0,
                                  open_progress=FALSE, evaluate=FALSE, ...)
        fit.call$x <- eval(fit.call$x)
        calls[[fold]] <- fit.call
    }

    ## evaluate the models
    message("Fitting ", num.folds, " models using ",
            min(cores, num.folds), " cores")
    par.fun <- function(fold) {
        fit <- eval(calls[[fold]])

        ## log pointwise predictive densities (pointwise test log-likelihood)
        lppd <- log_lik(fit, newdata=data[which(folds == fold), , drop=FALSE])
        return(list(lppd=lppd, fit=if (store.fits) fit else NULL))
    }
    if (.Platform$OS.type != "windows") {
        cv <- parallel::mclapply(X=1:num.folds, mc.cores=cores,
                                 mc.preschedule=FALSE, FUN=par.fun)
    } else { # windows
        cl <- parallel::makePSOCKcluster(cores)
        on.exit(parallel::stopCluster(cl))
        cv <- parallel::parLapply(X=1:num.folds, cl=cl, fun=par.fun)
    }

    ## expected log predictive densities
    elpds.unord <- unlist(lapply(cv, function(z) apply(z$lppd, 2, logMeanExp)))
    obs.idx <- unlist(lapply(1:num.folds, function(z) which(folds == z)))
    elpds <- elpds.unord[obs.idx]

    pointwise <- cbind(elpd_kfold=elpds, p_kfold=NA, kfoldic=-2 * elpds)
    estimates <- colSums(pointwise)
    se.est <- sqrt(N * apply(pointwise, 2, stats::var))
    out <- list(estimates=cbind(Estimate=estimates, SE=se.est),
                pointwise=pointwise)
    rownames(out$estimates) <- colnames(pointwise)
    if (store.fits) {
        fits <- array(list(), c(num.folds, 2), list(NULL, c("fit", "test.idx")))
        for (fold in 1:num.folds)
            fits[fold, ] <- list(fit=cv[[fold]][["fit"]],
                                 test.idx=which(folds == fold))
        out$fits <- fits
        out$data <- data
    }
    attr(out, "K") <- num.folds
    class(out) <- c("kfold", "loo")
    return(out)
}
