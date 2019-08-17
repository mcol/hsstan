##=============================================================================
##
## Copyright (c) 2019 Marco Colombo and Paul McKeigue
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


#' Summary for the fitted model
#'
#' @param object An object of class \code{hsstan}.
#' @param pars Vector of parameter names to be extracted. If \code{NULL}
#'        (default) then this refers to the set of predictors used in the
#'        model.
#' @param prob Width of the posterior intervals (0.95, by default).
#' @param digits Number of decimal places to be reported (2 by default).
#' @param sort Column name used to sort the results according to the absolute
#'        value of the column. If \code{NULL} (default) or the column name cannot
#'        be found, no sorting occurs.
#' @param decreasing Whether the results should be sorted in decreasing order
#'        when a valid column name is specified in \code{sort} (\code{TRUE} by
#'        default).
#' @param max.rows Maximum number of rows to be returned. If \code{NULL}
#'        (default) or 0, all results are returned.
#' @param ... Currently ignored.
#'
#' @return
#' A matrix with summaries from the posterior distribution of the parameters
#' of interest.
#'
#' @importMethodsFrom rstan summary
#' @method summary hsstan
#' @export
summary.hsstan <- function(object, pars=NULL, prob=0.95, digits=2,
                          sort=NULL, decreasing=TRUE, max.rows=NULL, ...) {
    validate.samples(object)
    validate.probability(prob)
    pars <- get.pars(object, pars)
    low <- (1 - prob) / 2
    upp <- 1 - low
    summary <- summary(object$stanfit, pars=pars, probs=c(low, upp))$summary
    summary[, "n_eff"] <- round(summary[, "n_eff"], 0)
    if (!is.null(sort) && sort %in% colnames(summary)) {
        ord <- order(abs(summary[, sort]), decreasing=decreasing)
        summary <- summary[ord, ]
    }
    if (!is.null(max.rows)) {
        if (max.rows > 0 && max.rows < nrow(summary))
            summary <- summary[1:max.rows, , drop=FALSE]
    }
    round(summary[, -match("se_mean", colnames(summary)), drop=FALSE],
                digits)
}

#' Print a summary for the fitted model
#'
#' @param x An object of class \code{hsstan}.
#' @param ... Further arguments to \code{\link{summary.hsstan}}.
#'
#' @export
print.hsstan <- function(x, ...) {
    print(summary(x, ...))
}

#' Posterior summary
#'
#' Produce a summary of the posterior samples for the quantities of interest.
#'
#' @param x An object containing or representing posterior samples. If \code{x}
#'        is a matrix, it should have size \code{S} by \code{Q}, where \code{S}
#'        is the number of posterior samples, and \code{Q} is the number of
#'        quantities of interest.
#' @param prob Width of the posterior intervals (0.95, by default).
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' A matrix with columns containing mean, standard deviation and posterior
#' intervals for the given quantities.
#'
#' @seealso
#' \code{\link{summary.hsstan}} to produce summaries of \code{hsstan} objects
#' that include the number of effective samples and the split-Rhat diagnostic.
#'
#' @export
posterior_summary <- function(x, prob, ...) {
    UseMethod("posterior_summary")
}

#' @rdname posterior_summary
#' @export
posterior_summary.default <- function(x, prob=0.95, ...) {
    validate.probability(prob)
    if (NCOL(x) == 1)
        x <- as.matrix(x)
    t(apply(x, 2, vector.summary, prob))
}

#' @rdname posterior_summary
#' @param pars Vector of parameter names to be extracted. If \code{NULL}
#'        (default) then this refers to the set of predictors used in the
#'        model.
#'
#' @export
posterior_summary.hsstan <- function(x, prob=0.95, pars=NULL, ...) {
    validate.samples(x)
    pars <- get.pars(x, pars)
    posterior_summary(as.matrix(x$stanfit, pars=pars), prob, ...)
}

#' Performance measures
#'
#' Extract measures of performance from plain or cross-validated results.
#'
#' @param obj An object or a list of objects of class \code{hsstan}.
#' @param out.csv Optional name of a file where the output can be saved in
#'        CSV format.
#'
#' @importFrom pROC roc
#' @importFrom stats cor var
#' @importFrom utils write.csv
#' @export
get.cv.performance <- function(obj, out.csv=NULL) {

    r2 <- function(y.pred, y.obs) {
        corr <- cor(y.pred, y.obs)
        if (corr < 0)
            return(0)
        return(corr^2)
    }
    auc <- function(y.pred, y.obs)
        as.numeric(roc(y.obs, y.pred, direction="<", quiet=TRUE)$auc)

    if (inherits(obj, "hsstan")) {
        obj <- list(fits=array(list(fit=obj, test.idx=1:nrow(obj$data)), c(1, 2)),
                    data=obj$data)
        colnames(obj$fits) <- c("fit", "test.idx")
    } else if (inherits(obj, c("kfold", "loo"))) {
        if (is.null(obj[["fits"]]))
            stop("No fitted models found, run 'kfold' with store.fits=TRUE.")
    } else {
        stop("Not an 'hsstan' or 'kfold' object.")
    }

    num.folds <- nrow(obj$fits)
    y.obs.all <- y.pred.all <- NULL
    llk <- perf <- rep(NA, num.folds)
    llk.all <- 0
    is.logistic <- is.logistic(obj$fits[[1]])
    perf.fun <- ifelse(is.logistic, auc, r2)

    ## loop over the folds
    for (fold in 1:num.folds) {
        hs.fold <- obj$fits[[fold]]
        test.idx <- obj$fits[, "test.idx"][[fold]]
        y.obs <- obj$data[test.idx, hs.fold$model.terms$outcome]

        ## retrieve the fitted values on test data
        newdata <- obj$data[test.idx, ]
        y.pred <- colMeans(posterior_linpred(hs.fold, newdata=newdata,
                                             transform=TRUE))
        llk[fold] <- sum(colMeans(log_lik(hs.fold, newdata=newdata)))
        perf[fold] <- perf.fun(y.pred, y.obs)
        y.obs.all <- c(y.obs.all, y.obs)
        y.pred.all <- c(y.pred.all, y.pred)
        llk.all <- llk.all + llk[fold]
    }
    set <- paste("Fold", 1:num.folds)

    ## compute log-likelihood and performance measure of the full model
    ## on the full vector of withdrawn observations
    set <- c(set, "Overall")
    llk <- c(llk, llk.all)
    perf <- c(perf, perf.fun(y.pred.all, y.obs.all))

    res <- data.frame(set=set, test.llk=llk, perf=perf,
                      stringsAsFactors=FALSE)
    colnames(res)[3] <- gsub("perf", ifelse(is.logistic, "auc", "r2"),
                             colnames(res)[3])

    if (num.folds == 1) {
        res <- res[1, ]
        if (!inherits(obj, c("kfold", "loo")))
            res$set <- "Non cross-validated"
    }

    if (!is.null(out.csv))
        write.csv(file=out.csv, res, row.names=FALSE)
    return(res)
}

#' Sampler statistics
#'
#' Report statistics on the parameters used in the sampler, the sampler
#' behaviour and the sampling time.
#'
#' @param object An object of class \code{hsstan}.
#'
#' @return
#' A matrix with \code{C + 1} rows, where \code{C} is the number of Markov
#' chains, reporting average acceptance probability, average stepsize, number
#' of divergent transitions, maximum tree depth, total number of gradient
#' evaluations, warmup and sample times in seconds.
#'
#' @export
sampler.stats <- function(object) {
    validate.hsstan(object)
    validate.samples(object)
    sp <- rstan::get_sampler_params(object$stanfit, inc_warmup=FALSE)
    accept.stat <- sapply(sp, function(x) mean(x[, "accept_stat__"]))
    stepsize <- sapply(sp, function(x) mean(x[, "stepsize__"]))
    divergences <- sapply(sp, function(x) sum(x[, "divergent__"]))
    treedepth <- sapply(sp, function(x) max(x[, "treedepth__"]))
    gradients <- sapply(sp, function(z) sum(z[, "n_leapfrog__"]))
    et <- round(rstan::get_elapsed_time(object$stanfit), 2)
    res <- cbind(accept.stat, stepsize, divergences, treedepth, gradients, et)
    avg <- colMeans(res)
    tot <- colSums(res)
    round(rbind(res, all=c(avg[1:2], tot[3], max(res[, 4]), tot[5:7])), 4)
}

#' Number of posterior samples
#'
#' Extracts the number of posterior samples stored in a fitted model.
#'
#' @param object An object of class \code{hsstan}.
#' @param ... Currently ignored.
#'
#' @return
#' The total number of posterior samples across the chains after discarding
#' the warmup iterations.
#'
#' @importFrom rstantools nsamples
#' @method nsamples hsstan
#' @aliases nsamples
#' @export nsamples
#' @export
nsamples.hsstan <- function(object, ...) {
    validate.samples(object)
    with(object$stanfit@sim, sum(n_save - warmup2))
}
