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


#' Validate an hsstan object
#'
#' Check that the object has been created by [hsstan()].
#'
#' @param obj An object to be checked.
#'
#' @return
#' Throws an error if the object is not an `hsstan` object.
#'
#' @noRd
validate.hsstan <- function(obj) {
    if (!inherits(obj, "hsstan")) {
        stop("Not an object of class 'hsstan'.")
    }
}

#' Validate the posterior samples
#'
#' Check that the object contains valid posterior samples in the
#' `stanfit` field.
#'
#' @param obj An object of class `hsstan`.
#'
#' @return
#' Throws an error if the object does not contain posterior samples.
#'
#' @noRd
validate.samples <- function(obj) {
    if (!inherits(obj$stanfit, "stanfit")) {
        stop("No valid posterior samples stored in the 'hsstan' object.")
    }
}

#' Validate new data
#'
#' Check that the new data contains all variables used in the model and no
#' missing values, and generate the corresponding model matrix.
#'
#' @param obj Object of class `hsstan`.
#' @param newdata Optional data frame containing the variables used in the
#'        model. If `NULL`, the model matrix used when fitting the model
#'        is returned.
#'
#' @return
#' A design matrix corresponding to the variables used in the model.
#'
#' @noRd
validate.newdata <- function(obj, newdata) {
    if (is.null(newdata))
        newdata <- obj$data
    else if (!inherits(newdata, c("data.frame", "matrix")))
        stop("'newdata' must be a data frame or a matrix.")
    if (nrow(newdata) == 0 || ncol(newdata) == 0)
        stop("'newdata' contains no rows or no columns.")

    ## only check for NAs in the variables used in the model
    vars <- with(obj$model.terms, c(outcome, unpenalized, penalized))
    newdata <- newdata[, colnames(newdata) %in% vars, drop=FALSE]
    if (any(is.na(newdata)))
        stop("'newdata' contains missing values.")

    ## this adds the intercept column back
    ordered.model.matrix(as.data.frame(newdata),
                         obj$model.terms$unpenalized,
                         obj$model.terms$penalized)
}

#' Validate a model formula
#'
#' Check that the formula that specifies a model contains all required elements.
#'
#' @param model Formula to be checked.
#' @param penalized Vector of names for the penalized predictors.
#'
#' @return
#' A list containing the formula representing the covariates model, the name of
#' the outcome variable, the names of the upenalized and penalized predictors.
#'
#' @importFrom stats as.formula terms
#' @noRd
validate.model <- function(model, penalized) {
    if (is.character(model) && length(model) > 1)
        stop("Model formula specified incorrectly.")
    model <- as.formula(model)
    tt <- terms(model)
    if (attr(tt, "response") == 0)
        stop("No outcome variable specified in the model.")
    if (attr(tt, "intercept") == 0)
        stop("Models with no intercept are not supported.")
    if (length(penalized) > 0 && !is.character(penalized))
        stop("'penalized' must be a character vector.")
    if (any(grepl("[:*]", penalized)))
        stop("Interaction terms in penalized predictors are not supported.")
    penalized <- setdiff(unique(trimws(penalized)), "")
    return(list(outcome=as.character(model)[2],
                unpenalized=setdiff(attr(tt, "term.labels"), penalized),
                penalized=penalized))
}

#' Validate the model data
#'
#' Check if the model data can be used with the given model formula and
#' penalized predictors.
#'
#' @param x An object to be checked.
#' @param model Validated model formula.
#'
#' @return
#' A data frame containing the model data. A factor or logical outcome variable
#' is replaced by its numeric equivalent.
#'
#' @noRd
validate.data <- function(x, model) {
    if (!inherits(x, c("data.frame", "matrix")))
        stop("'x' must be a data frame or a matrix.")
    x <- as.data.frame(x)
    validate.variables(x, model$outcome)
    validate.variables(x, c(model$unpenalized, model$penalized))
    x[[model$outcome]] <- validate.outcome(x[[model$outcome]])
    return(x)
}

#' Validate variables
#'
#' Check that the required variables are in the dataset.
#'
#' @param x Data frame containing the variables of interest.
#' @param variables Vector of variable names.
#'
#' @return
#' Throws if variables are not present in the dataset or contain missing values.
#'
#' @noRd
validate.variables <- function(x, variables) {
    ## unpack interaction terms
    variables <- unique(unlist(strsplit(as.character(variables), ":")))
    if (length(variables) == 0)
        stop("No predictors present in the model.")
    var.match <- match(variables, colnames(x))
    if (anyNA(var.match))
        stop(collapse(variables[is.na(var.match)]), " not present in 'x'.")
    if (anyNA(x[, variables]))
        stop("Model variables contain missing values.")
}

#' Validate the outcome variable
#'
#' Check that the outcome variable can be converted to a valid numerical
#' vector.
#'
#' @param y Outcome vector to be checked.
#'
#' @return
#' A numeric vector.
#'
#' @noRd
validate.outcome <- function(y) {
    if (is.factor(y)) {
        if (nlevels(y) != 2)
            stop("A factor outcome variable can only have two levels.")
        y <- as.integer(y) - 1
    }
    if (!(is.numeric(y) || is.logical(y)))
        stop("Outcome variable of invalid type.")
    return(as.numeric(y))
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
            stop("Outcome variable must contain two classes with family=binomial.")
        if (!is.factor(y) && any(y < 0 | y > 1))
            stop("Outcome variable must contain 0-1 values with family=binomial.")
    }

    return(family)
}

#' Validate a vector of indices
#'
#' @param x Vector to be checked.
#' @param N Maximum valid index.
#' @param name Name of the vector to report in error messages.
#' @param throw.duplicates Whether the function should throw if the vector
#'        contains duplicate elements (`TRUE` by default).
#'
#' @return
#' Throws an error if the given vector is not an integer vector or contains
#' missing, out of bounds or duplicate indices (if `throw.duplicates` is `TRUE`).
#'
#' @noRd
validate.indices <- function(x, N, name, throw.duplicates=TRUE) {
    if (anyNA(x))
        stop("'", name, "' contains missing values.")
    if (!is.numeric(x) || NCOL(x) > 1 || any(x != as.integer(x)))
        stop("'", name, "' must be an integer vector.")
    if (length(x) < 2)
        stop("'", name, "' must contain at least two elements.")
    if (any(x < 1 | x > N))
        stop("'", name, "' contains out of bounds indices.")
    if (throw.duplicates && any(duplicated(x)))
        stop("'", name, "' contains duplicate indices.")
}

#' Validate the cross-validation folds
#'
#' @param folds Folds to be checked or `NULL`.
#' @param N Number of observations.
#'
#' @return
#' An integer vector with one element per observation indicating the
#' cross-validation fold in which the observation should be withdrawn.
#'
#' @noRd
validate.folds <- function(folds, N) {
    if (is.null(folds))
        return(rep(1, N))
    validate.indices(folds, N, "folds", throw.duplicates=FALSE)
    if (length(folds) != N)
        stop("'folds' should have length ", N, ".")
    K <- length(unique(folds))
    if (!all(1:K %in% folds))
        stop("'folds' must contain all indices up to ", K, ".")
    folds <- as.integer(folds)
}

#' Validate start.from
#'
#' Check that the predictor names provided is a valid subset of the variables
#' used in the model.
#'
#' @param obj An object of class `hsstan`.
#' @param start.from Vector to be checked.
#'
#' @return
#' A list of two elements: the names of the model terms matching `start.from`
#' and a vector of indices corresponding to the names listed in `start.from`.
#' Throws an error if any of the names mentioned does not match those available
#' in the model terms.
#'
#' @noRd
validate.start.from <- function(obj, start.from) {
    unp.terms <- obj$model.terms$unpenalized
    unp.betas <- names(obj$betas$unpenalized)
    mod.terms <- c(unp.terms, obj$model.terms$penalized)
    mod.betas <- c(unp.betas, names(obj$betas$penalized))
    start.from <- setdiff(start.from, "")
    if (is.null(start.from)) {
        if (length(obj$model.terms$penalized) > 0)
            return(list(start.from=unp.terms, idx=seq_along(unp.betas)))
        else
            return(list(start.from=character(0), idx=1))
    }
    if (length(start.from) == 0)
        return(list(start.from=character(0), idx=1))
    if (anyNA(start.from))
        stop("'start.from' contains missing values.")
    var.match <- match(start.from, mod.terms)
    if (anyNA(var.match))
        stop("'start.from' contains ", collapse(start.from[is.na(var.match)]),
             ", which cannot be matched.")

    ## unpack interaction terms so that also main effects are matched
    start.from <- mod.terms[mod.terms %in%
                            c(start.from, unlist(strsplit(start.from, ":")))]
    chosen <- expand.terms(obj$data, start.from)

    ## also consider interaction terms in reverse order
    chosen <- c(chosen, sapply(strsplit(chosen[grep(":", chosen)], ":"),
                               function(z) c(z, paste(rev(z), collapse=":"))))
    return(list(start.from=start.from, idx=which(mod.betas %in% chosen)))
}

#' Validate adapt.delta
#'
#' Check that an adaptation acceptance probability is valid.
#'
#' @param adapt.delta Value to be checked.
#'
#' @return
#' Throws an error if the given value is not a valid acceptance probability
#' for adaptation.
#'
#' @noRd
validate.adapt.delta <- function(adapt.delta) {
    if (!is.numeric(adapt.delta) || length(adapt.delta) != 1) {
        stop("'adapt.delta' must be a single numerical value.")
    }
    if (adapt.delta < 0.8) {
        stop("'adapt.delta' must be at least 0.8.")
    }
    if (adapt.delta >= 1) {
        stop("'adapt.delta' must be less than 1.")
    }
}

#' Validate a probability
#'
#' Check that a probability value is valid.
#'
#' @param prob Value to be checked.
#'
#' @return
#' Throws an error if the given value is not a valid probability.
#'
#' @noRd
validate.probability <- function(prob) {
    if (length(prob) != 1 || prob <= 0 || prob >= 1)
        stop("'prob' must be a single value between 0 and 1.\n")
}

#' Validate arguments passed to rstan
#'
#' Ensure that the options to be passed to \code{\link[rstan]{sampling}} are
#' valid, as to work around rstan issue #681.
#'
#' @param ... List of arguments to be checked.
#'
#' @return
#' Throws an error if any argument is not valid for \code{\link[rstan]{sampling}}.
#'
#' @noRd
validate.rstan.args <- function(...) {
    valid.args <- c("chains", "cores", "pars", "thin", "init", "check_data",
                    "sample_file", "diagnostic_file", "verbose", "algorithm",
                    "control", "open_progress", "show_messages", "chain_id",
                    "init_r", "test_grad", "append_samples", "refresh",
                    "save_warmup", "enable_random_init", "iter", "warmup")
    dots <- list(...)
    for (arg in names(dots))
        if (!arg %in% valid.args)
           stop("Argument '", arg, "' not recognized.")
}

#' Parameter names
#'
#' Get the parameter names corresponding to the regression coefficients or
#' matching a regular expression.
#'
#' @param obj An object of class `hsstan`.
#' @param pars Regular expression to match a parameter name, or `NULL`
#'        to retrieve the names of all regression coefficients.
#'
#' @return
#' A character vector.
#'
#' @noRd
get.pars <- function(object, pars) {
    if (is.null(pars))
        pars <- grep("^beta_", object$stanfit@model_pars, value=TRUE)
    else {
        if (!is.character(pars))
            stop("'pars' must be a character vector.")
        get.pars <- function(x) grep(x, object$stanfit@sim$fnames_oi, value=TRUE)
        pars <- unlist(lapply(pars, get.pars))
        if (length(pars) == 0)
            stop("No pattern in 'pars' matches parameter names.")
    }
    return(pars)
}

#' Create a design matrix with all unpenalized predictors first
#'
#' This is required as `model.matrix` puts the interaction terms after the
#' penalized predictors, but the Stan models expects all unpenalized terms to
#' appear before the penalized ones.
#'
#' @param x Data frame containing the variables of interest.
#' @param unpenalized Vector of variable names for the unpenalized covariates.
#' @param penalized Vector of variable names for the penalized predictors.
#'
#' @return
#' A design matrix with all unpenalized covariates (including interaction terms)
#' before the penalized predictors.
#'
#' @importFrom stats model.matrix reformulate
#' @noRd
ordered.model.matrix <- function(x, unpenalized, penalized) {
    X <- model.matrix(reformulate(c(unpenalized, penalized)), data=x)
    if (any(grepl("[:*]", unpenalized)))
        X <- X[, c(expand.terms(x, unpenalized), expand.terms(x, penalized)[-1])]
    return(X)
}

#' Expand variable names into formula terms
#'
#' @param x Data frame containing the variables of interest.
#' @param variables Vector of variable names.
#'
#' @return
#' A vector of variable names expanded by factor levels and interaction terms.
#'
#' @importFrom stats model.matrix reformulate
#' @noRd
expand.terms <- function(x, variables) {
    if (length(variables) == 0)
        return(character(0))
    colnames(model.matrix(reformulate(variables), x[1, ]))
}

#' Summarize a vector
#'
#' @param x A numerical vector.
#' @param prob Width of the interval between quantiles.
#'
#' @return
#' The mean, standard deviation and quantiles for the input vector.
#'
#' @noRd
vector.summary <- function(x, prob) {
    lower <- (1 - prob) / 2
    upper <- 1 - lower
    c(mean=mean(x), sd=stats::sd(x), stats::quantile(x, c(lower, upper)))
}

#' Check whether the model fitted is a logistic regression model.
#'
#' @param obj An object of class `hsstan`.
#'
#' @return
#' `TRUE` for logistic regression models, `FALSE` otherwise.
#'
#' @noRd
is.logistic <- function(obj) {
    obj$family$family == "binomial"
}

#' Comma-separated string concatenation
#'
#' Collapse the elements of a character vector into a comma-separated string.
#'
#' @param x Character vector.
#'
#' @return
#' A comma-separated string where each element of the original vector is
#' surrounded by single quotes.
#'
#' @noRd
collapse <- function(x) {
    paste0("'", x, "'", collapse=", ")
}

#' Fast computation of correlations
#'
#' This provides a loopless version of the computation of the correlation
#' coefficient between observed and predicted outcomes.
#'
#' @param y Vector of observed outcome.
#' @param x Matrix with as many columns as the number of elements in `y`,
#'          where each row corresponds to a predicted outcome.
#'
#' @return
#' A vector of correlations with as many elements as the number of rows in `x`.
#'
#' @noRd
fastCor <- function(y, x) {
    yx <- rbind(y, x)
    if (.Machine$sizeof.pointer == 8) {
        yx <- yx - rowMeans(yx)
        yx <- yx / sqrt(rowSums(yx^2))
        corr <- tcrossprod(yx, yx)
    } else {
        corr <- stats::cor(yx)
    }
    return(corr[-1, 1])
}

#' Log of sum of exponentials
#'
#' @noRd
logSumExp <- function(x) {
    xmax <- max(x)
    xmax + log(sum(exp(x - xmax)))
}

#' Log of average of exponentials
#'
#' @noRd
logMeanExp <- function(x) {
    logSumExp(x) - log(length(x))
}
