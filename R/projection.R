##=============================================================================
##
## Copyright (c) 2017-2019 Marco Colombo and Paul McKeigue
##
## Parts of the code are based on https://github.com/jpiironen/rstan-varsel
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


#' Compute projections of full predictors on to subspace of predictors
#'
#' @param x Design matrix.
#' @param fit Matrix of fitted values for the full model.
#' @param sigma2 Residual variance (1 for logistic regression).
#' @param indproj Vector of indices of the columns of \var{x} that form the
#'        projection subspace.
#' @param is.logistic Set to \code{TRUE} for a logistic regression model.
#'
#' @importFrom stats binomial glm.fit quasibinomial
#' @keywords internal
lm_proj <- function(x, fit, sigma2, indproj, is.logistic) {

    P <- ncol(x)
    S <- ncol(fit)

    ## pick the Q columns of x that form the projection subspace
    xp <- x[, indproj] # matrix of dimension N x Q

    ## logistic regression model
    if (is.logistic) {

        ## compute the projection for each sample
        wp <- parallel::mclapply(X=1:S, mc.preschedule=TRUE,
                                 FUN=function(z) {
                                     glm.fit(xp, fit[, z],
                                             family=quasibinomial())$coefficients
                                 })
        wp <- matrix(unlist(wp, use.names=FALSE), ncol=S)

        ## estimate the KL divergence between full and projected model
        fitp <- binomial()$linkinv(xp %*% wp)
        kl <- 0.5 * mean(colMeans(binomial()$dev.resids(fit, fitp, 1)))
        sigma2p <- 1
    }

    ## linear regression model
    else {
        ## solve the projection equations
        wp <- solve(crossprod(xp, xp), crossprod(xp, fit)) # Q x S matrix

        ## fit of the projected model
        fitp <- xp %*% wp

        ## estimate the KL divergence between full and projected model
        sigma2p <- sigma2 + colMeans((fit - fitp)^2)
        kl <- mean(0.5 * log(sigma2p / sigma2))
    }

    ## reshape wp so that it has same dimension (P x S) as t(x), and zeros for
    ## those variables that are not included in the projected model
    wp.all <- matrix(0, P, S)
    wp.all[indproj, ] <- wp
    return(list(w=wp.all, sigma2=sigma2p, fit=fitp, kl=kl))
}

#' Next variable to enter the current submodel
#'
#' Return the index of the variable that should be added to the current model
#' according to the smallest KL-divergence (linear regression) or the largest
#' score test (logistic regression).
#'
#' @param x Design matrix.
#' @param sigma2 Residual variance (1 for logistic regression).
#' @param fit Matrix of fitted values for the full model.
#' @param fitp Matrix of fitted values for the projected model.
#' @param chosen Vector of indices of the columns of \var{x} in the current
#'        submodel.
#' @param is.logistic Set to \code{TRUE} for a logistic regression model.
#'
#' @keywords internal
choose.next <- function(x, sigma2, fit, fitp, chosen, is.logistic) {
    notchosen <- setdiff(1:ncol(x), chosen)
    if (!is.logistic) {
        kl <- parallel::mclapply(X=1:length(notchosen), mc.preschedule=TRUE,
                                 FUN=function(z) {
                                     ind <- sort(c(chosen, notchosen[z]))
                                     lm_proj(x, fit, sigma2, ind, FALSE)$kl
                                 })
        idx.selected <- which.min(unlist(kl))
        return(notchosen[idx.selected])
    }

    ## score test
    dinv.link <- t(fitp * (1 - fitp))
    yminusexp <- t(fit - fitp)
    U <- colMeans(yminusexp %*% x[, notchosen] / sigma2)
    V <- colMeans(dinv.link %*% x[, notchosen]^2 / sigma2)
    idx.selected <- which.max(U^2 / V)
    return(notchosen[idx.selected])
}

#' Forward selection minimizing KL-divergence in projection
#'
#' @param samples Object produced by \code{\link{hsstan}}.
#' @param max.num.pred Maximum number of predictors after which the selection
#'        procedure should stop.
#' @param out.csv If not \code{NULL}, the name of a CSV file to save the
#'        output to.
#'
#' @importFrom stats dbinom dnorm
#' @importFrom utils write.csv
#' @export
projsel <- function(samples, max.num.pred=30, out.csv=NULL) {

    fit.submodel <- function(x, sigma2, fit, chosen, xt, yt, is.logistic) {

        ## projected parameters
        submodel <- lm_proj(x, fit, sigma2, chosen, is.logistic)
        eta <- xt %*% submodel$w

        ## expected log predictive density on test set
        if (is.logistic) {
            pd <- dbinom(yt, 1, binomial()$linkinv(eta), log=TRUE)
        }
        else {
            pd <- t(cbind(sapply(1:nrow(eta),
                                 function(z) dnorm(yt[z], eta[z, ],
                                                   sqrt(sigma2), log=TRUE))))
        }
        elpd <- sum(rowMeans(pd))

        return(list(fit=submodel$fit, kl=submodel$kl, elpd=elpd))
    }

    validate.hsstan(samples)
    validate.samples(samples)

    ## check that the model contains penalized predictors
    if (length(samples$model.terms$penalized) == 0) {
        stop("Model doesn't contain penalized predictors.")
    }

    x <- xt <- validate.newdata(samples, samples$data)
    yt <- samples$data[[samples$model.terms$outcome]]
    stanfit <- samples$stanfit

    is.logistic <- is.logistic(samples)
    sigma2 <- if (is.logistic) 1 else as.matrix(stanfit, pars="sigma")^2

    ## fit of the full model (matrix of dimension N x S)
    fit <- t(posterior_linpred(samples, transform=TRUE))

    ## U is number of unpenalized variables (always chosen) including intercept
    P <- length(c(samples$betas$unpenalized, samples$betas$penalized))
    U <- length(samples$betas$unpenalized)
    kl.elpd <- NULL
    cat(sprintf("%58s  %8s %11s\n", "Model", "KL", "ELPD"))
    report.iter <- function(msg, kl, elpd)
        cat(sprintf("%58s  %8.5f  %8.5f\n", substr(msg, 1, 55), kl, elpd))

    ## start from the intercept only model
    sub <- fit.submodel(x, sigma2, fit, 1, xt, yt, is.logistic)
    kl.elpd <- rbind(kl.elpd, c(sub$kl, sub$elpd))
    report.iter("Intercept only", sub$kl, sub$elpd)

    ## start from the model having only unpenalized variables
    chosen <- 1:U
    notchosen <- setdiff(1:P, chosen)
    sub <- fit.submodel(x, sigma2, fit, chosen, xt, yt, is.logistic)
    fitp <- sub$fit
    kl.elpd <- rbind(kl.elpd, c(sub$kl, sub$elpd))
    report.iter("Unpenalized covariates", sub$kl, sub$elpd)

    ## add variables one at a time
    for (k in 1:(P - U)) {
        sel.idx <- choose.next(x, sigma2, fit, fitp, chosen, is.logistic)
        chosen <- c(chosen, sel.idx)

        ## evaluate current submodel according to projected parameters
        sub <- fit.submodel(x, sigma2, fit, chosen, xt, yt, is.logistic)
        fitp <- sub$fit
        kl.elpd <- rbind(kl.elpd, c(sub$kl, sub$elpd))
        report.iter(colnames(x)[sel.idx], sub$kl, sub$elpd)

        if (length(chosen) - U == max.num.pred)
            break
    }

    ## evaluate the full model
    full <- fit.submodel(x, sigma2, fit, 1:P, xt, yt, is.logistic)

    res <- data.frame(var=c("Intercept only",
                            "Unpenalized covariates",
                            colnames(x)[setdiff(chosen, 1:U)]),
                      kl=kl.elpd[, 1], elpd=kl.elpd[, 2],
                      delta.elpd=kl.elpd[, 2] - full$elpd,
                      stringsAsFactors=FALSE)
    if (!is.null(out.csv))
        write.csv(file=out.csv, res, row.names=FALSE)

    class(res) <- c("projsel", "data.frame")
    return(res)
}

#' Plot of relative explanatory power of predictors
#'
#' The relative explanatory power of predictors is computed according to the KL
#' divergence from the full model to each submodel, scaled in such a way that
#' the baseline set of covariates are at 0, while the full model is at 1.
#'
#' A function of name \code{getfullname} to match variable names to full
#' names is searched on the current workspace, and if found it is used to
#' transform the labels before creating the plot.
#'
#' @param x A data frame created by \code{\link{projsel}}.
#' @param title Title of the plot. If \code{NULL}, no title is displayed.
#' @param max.labels Maximum number of points to be labelled. If \code{NULL},
#'        all those present in the \code{x} file are displayed.
#' @param font.size Size of the textual elements (labels and axes).
#' @param hadj,vadj Horizontal and vertical adjustment for the labels.
#' @param ... Currently ignored.
#'
#' @return
#' A \pkg{ggplot2} object showing the relative incremental contribution of each
#' predictor starting from the initial set of unpenalized covariates.
#'
#' @import ggplot2
#' @method plot projsel
#' @export
plot.projsel <- function(x, title=NULL, max.labels=NULL, font.size=12,
                         hadj=0.05, vadj=0, ...) {

    ## get full variable names if possible
    sel <- x
    labs <- tryCatch(get("getfullname")(sel$var),
                     error=function(e) sel$var)
    labs <- gsub(" \\(.*\\)$", "", labs)
    if (!is.null(max.labels)) {
        labs[-c(1:(max.labels + 1))] <- ""
    }

    ## convert from points to millimetres
    geom.text.size <- font.size * 25.4 / 72

    ## relative explanatory power
    sel$rel <- 1 - sel$kl / sel$kl[1]

    x <- seq(nrow(sel)) - 1
    text_idx <- x < mean(x) | x - floor(x / 2) * 2 == 1
    p <- ggplot(data=sel, aes(x=x, y=rel, label=labs)) +
      coord_cartesian(ylim=range(c(0, 1))) +
      geom_line() + geom_point(size=geom.text.size / 3) +
      geom_text(aes(x=x + ifelse(text_idx, hadj, -hadj),
                    y=rel + ifelse(text_idx, -vadj, vadj)),
                size=geom.text.size,
                hjust=ifelse(text_idx, "left", "right")) +
      xlab("Number of biomarkers") +
      ylab("Relative explanatory power") +
      theme(text=element_text(size=font.size))

    if (!is.null(title))
        p <- p + ggtitle(title)

    return(p)
}
