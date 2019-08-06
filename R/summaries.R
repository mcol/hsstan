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
#'        (default) then this refers to the set of covariates and biomarkers
#'        (if available).
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
    if (is.null(pars))
        pars <- grep("^beta_", object$stanfit@model_pars, value=TRUE)
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

#' Sampler parameters
#'
#' Report the parameters used for the sampler.
#'
#' @param object An object of class \code{hsstan}.
#'
#' @return
#' A matrix with as many rows as number of Markov chains reporting the average
#' acceptance probability, the total number of divergent transitions, the maximum
#' tree depth and the total number of gradient evaluations.
#'
#' @export
sampler.params <- function(object) {
    validate.hsstan(object)
    validate.samples(object)
    sp <- rstan::get_sampler_params(object$stanfit, inc_warmup=FALSE)
    accept.stat <- sapply(sp, function(x) mean(x[, "accept_stat__"]))
    divergences <- sapply(sp, function(x) sum(x[, "divergent__"]))
    treedepth <- sapply(sp, function(x) max(x[, "treedepth__"]))
    n.gradients <- sapply(sp, function(z) sum(z[, "n_leapfrog__"]))
    round(cbind(accept.stat, divergences, treedepth, n.gradients), 4)
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
