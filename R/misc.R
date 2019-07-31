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


#' Validate hsstan object
#'
#' Checks that the object has been created by \code{\link{hsstan}}.
#'
#' @param obj An object to be checked.
#'
#' @return
#' Throws an error if the object is not an \code{hsstan} object.
#'
#' @noRd
validate.hsstan <- function(obj) {
    if (!inherits(obj, "hsstan")) {
        stop("Not an object of class 'hsstan'.")
    }
}

#' Validate the posterior samples
#'
#' Checks that the object contains valid posterior samples in the
#' \code{stanfit} field.
#'
#' @param obj An object of class \code{hsstan}.
#'
#' @return
#' Throws an error if the object does not contain posterior samples.
#'
#' @noRd
validate.samples <- function(obj) {
    if (!inherits(obj$stanfit, "stanfit")) {
        stop("No posterior samples found, run 'hsstan' with store.samples=TRUE.")
    }
}

#' Validate new data
#'
#' Checks that the new data contains all variables used in the model with no
#' missing values, and generates the corresponding model matrix.
#'
#' @param obj Object of class \code{hsstan}.
#' @param newdata Optional data frame containing the variables used in the
#'        model. If \code{NULL}, the model matrix used when fitting the model
#'        is returned.
#'
#' @return
#' A model matrix corresponding to the variables used in the model.
#'
#' @importFrom stats model.matrix reformulate
#' @noRd
validate.newdata <- function(obj, newdata) {

    if (is.null(newdata))
        return(obj$data)

    ## only check for NAs in the variables used in the model
    vars <- c(obj$covariates, obj$biomarkers)
    newdata <- newdata[, colnames(newdata) %in% vars, drop=FALSE]
    if (any(is.na(newdata)))
        stop("NAs are not allowed in 'newdata'.")

    ## this adds the intercept column back
    newdata <- model.matrix(reformulate(vars), newdata)

    return(newdata)
}

#' Validate the family argument
#'
#' Ensure that the family argument has been specified correctly.
#' This is inspired by code in \code{\link{glm}}.
#'
#' @param family Family argument to test.
#' @param y Outcome variable.
#'
#' @return
#' A valid family. The function throws an error if the family argument cannot
#' be used.
#'
#' @importFrom methods is
#' @noRd
validate.family <- function(family, y) {
    if (missing(family))
        stop("Argument of 'family' is missing.")
    if (is.character(family))
        tryCatch(
            family <- get(family, mode="function", envir=parent.frame(2)),
                          error=function(e)
                              stop("'", family, "' is not a valid family.")
        )
    if (is.function(family))
        family <- family()
    if (!is(family, "family"))
        stop("Argument of 'family' is not a valid family.")
    if (!family$family %in% c("gaussian", "binomial"))
        stop("Only 'gaussian' and 'binomial' are supported families.")

    if (family$family == "binomial") {
        if (length(table(y)) != 2)
            stop("y must contain two classes with family=binomial().")
      if (!is.factor(y) && any(y < 0 | y > 1))
          stop("y must contain 0-1 values with family=binomial().")
    }

    return(family)
}

#' Validate adapt.delta
#'
#' Checks that the \code{adapt.delta} value is valid.
#'
#' @param adapt.delta Value to be checked.
#'
#' @return
#' A valid acceptance probability for adaptation.
#' Throws an error if the value passed it not a valid acceptance probability
#' for adaptation.
#'
#' @noRd
validate.adapt.delta <- function(adapt.delta) {
    if (!is.numeric(adapt.delta) || length(adapt.delta) != 1) {
        stop("'adapt.delta' must be a single numeric value.")
    }
    if (adapt.delta < 0.8) {
        stop("'adapt.delta' must be at least 0.8.")
    }
    if (adapt.delta >= 1) {
        stop("'adapt.delta' must be below 1.")
    }
}
