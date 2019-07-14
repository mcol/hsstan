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


#' Prints a summary for the fitted model
#'
#' @param x An object of class \code{hsstan}.
#' @param pars Vector of parameter names to be extracted. If \code{NULL}
#'        (default) then this refers to the set of covariates and biomarkers
#'        (if available).
#' @param probs Quantiles to be printed.
#' @param digits Number of significant digits to be printed.
#' @param ... Currently ignored.
#'
#' @importMethodsFrom rstan summary
print.hsstan <- function(x, pars=NULL, probs=c(0.025, 0.5, 0.975),
                         digits=2, ...) {
    if (!inherits(x$stanfit, "stanfit")) {
        cat("No posterior samples found, run 'hsstan' with store.samples=TRUE.\n")
        return(invisible(NULL))
    }
    if (is.null(pars))
        pars <- grep("^beta_", x$stanfit@model_pars, value=TRUE)
    summary <- summary(x$stanfit, pars=pars, probs=probs)$summary
    summary[, "n_eff"] <- round(summary[, "n_eff"], 0)
    print(round(summary[, -match(c("se_mean", "sd"), colnames(summary))], digits))
    return(invisible(summary))
}
