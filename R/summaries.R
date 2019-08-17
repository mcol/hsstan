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
