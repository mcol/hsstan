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

    Q <- length(indproj)
    N <- nrow(x)
    P <- ncol(x)
    S <- ncol(fit)

    ## pick the Q columns of x that form the projection subspace
    xp <- x[, indproj] # matrix of dimension N x Q

    ## logistic regression model
    if (is.logistic) {

        ## compute the projection for each sample
        s <- NULL   # silence a note raised by R CMD check
        wp <- foreach(s=1:S, .combine=cbind) %dopar% {
            glm.fit(xp, fit[, s], family=quasibinomial())$coefficients
        }

        ## estimate the KL divergence between full and projected model
        fitp <- 1 / (1 + exp(-xp %*% wp))
        kl <- 0.5 * mean(colMeans(binomial()$dev.resids(fit, fitp, 1)))
        sigma2p <- 1
    }

    ## linear regression model
    else {
        ## solve the projection equations
        wp <- solve(t(xp) %*% xp, t(xp) %*% fit) # matrix of dimension Q x S

        ## fit of the projected model
        fitp <- xp %*% wp

        ## estimate the KL divergence between full and projected model
        sigma2p <- sigma2 + colMeans((fit - fitp)^2)
        kl <- mean(0.5 * log(sigma2p / sigma2))
    }

    ## reshape wp so that it has same dimension (P x N) as t(x), and zeros for
    ## those variables that are not included in the projected model
    wptemp <- matrix(0, P, S)
    wptemp[indproj, ] <- wp
    wp <- wptemp
    return(list(w=wp, sigma2=sigma2p, fit=fitp, kl=kl))
}

#' Return the next variable to be added to the current submodel
#'
#' Returns the index of the variable that should be added to the current model
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
#' @importFrom foreach foreach %dopar%
#' @keywords internal
choose.next <- function(x, sigma2, fit, fitp, chosen, is.logistic) {
    notchosen <- setdiff(1:ncol(x), chosen)
    if (is.logistic) {
        dinv.link <- t(fitp * (1 - fitp))
    }
    else {
        i <- NULL   # silence a note raised by R CMD check
        kl <- foreach(i=1:length(notchosen), .combine=c) %dopar% {
            ind <- sort(c(chosen, notchosen[i]))
            lm_proj(x, fit, sigma2, ind, FALSE)$kl
        }
        idx.selected <- which.min(kl)
        return(notchosen[idx.selected])
    }
    yminusexp <- t(fit - fitp)

    ## score test
    U <- colMeans(yminusexp %*% x[, notchosen] / sigma2)
    V <- colMeans(dinv.link %*% x[, notchosen]^2 / sigma2)
    idx.selected <- which.max(U^2 / V)
    return(notchosen[idx.selected])
}

#' Forward selection minimizing KL-divergence in projection
#'
#' @param samples Object produced by \code{\link{sample.stan}}.
#' @param max.num.pred Maximum number of predictors after which the selection
#'        procedure should stop.
#' @param out.csv If not \code{NULL}, the name of a CSV file to save the
#'        output to.
#'
#' @importFrom stats dbinom dnorm
#' @importFrom utils write.csv
#' @importMethodsFrom rstan extract
#' @export
projsel <- function(samples, max.num.pred=30, out.csv=NULL) {

    fit.submodel <- function(x, w, sigma2, fit, chosen, xt, yt, is.logistic) {

        ## projected parameters
        submodel <- lm_proj(x, fit, sigma2, chosen, is.logistic)
        wp <- submodel$w
        sigma2p <- submodel$sigma2

        ## mean log predictive density on test set (mean test log-likelihood
        ## per observation)
        if (is.logistic) {
            pd <- dbinom(yt, 1, 1 / (1 + exp(-xt %*% wp)))
        }
        else {
            pd <- dnorm(yt, xt %*% wp, sqrt(sigma2p))
        }
        elpd <- sum(log(rowMeans(pd)))

        return(list(fit=submodel$fit, kl=submodel$kl, elpd=elpd))
    }

    validate.hsstan(samples)
    if (!inherits(samples$stanfit, "stanfit"))
        stop("No posterior samples found: run 'hsstan' with store.samples=TRUE.")

    ## check that the model contains penalized predictors
    if (is.null(samples$betas$penalized)) {
        stop("Model doesn't contain penalized predictors.")
    }

    x <- samples$data
    if (is.null(samples$withdrawn.data)) {
        xt <- samples$data
        yt <- samples$y
    }
    else {
        xt <- samples$withdrawn.data
        yt <- samples$y_test
    }
    stanfit <- samples$stanfit

    beta.samples <- extract(stanfit, pars=c("beta_u", "beta_p"))
    w <- t(cbind(beta.samples$beta_u, beta.samples$beta_p)) # P x S
    sigma2 <- tryCatch(unlist(extract(stanfit, pars=c("sigma")))^2,
                       error=function(e) return(1))
    is.logistic <- length(sigma2) == 1 && sigma2 == 1

    ## fit of the full model (matrix of dimension N x S)
    fit <- x %*% w
    if (is.logistic)
        fit <- 1 / (1 + exp(-fit))

    ## U is number of unpenalized variables (always chosen) including intercept
    P <- ncol(x)
    U <- ncol(beta.samples$beta_u)
    kl <- elpd <- rep(0, P - U + 1)
    cat(sprintf("%58s  %8s %9s\n", "Model", "KL", "ELPD"))

    ## start from the model having only unpenalized variables
    chosen <- 1:U
    notchosen <- setdiff(1:P, chosen)
    sub <- fit.submodel(x, w, sigma2, fit, chosen, xt, yt, is.logistic)
    fitp <- sub$fit
    kl[1] <- sub$kl
    elpd[1] <- sub$elpd
    cat(sprintf("%58s  %8.5f  %8.5f\n",
                "Initial set of covariates", kl[1], elpd[1]))

    ## add variables one at a time
    for (k in 2:(P - U + 1)) {
        sel.idx <- choose.next(x, sigma2, fit, fitp, chosen, is.logistic)
        chosen <- c(chosen, sel.idx)

        ## evaluate current submodel according to projected parameters
        sub <- fit.submodel(x, w, sigma2, fit, chosen, xt, yt, is.logistic)
        fitp <- sub$fit
        kl[k] <- sub$kl
        elpd[k] <- sub$elpd
        cat(sprintf(" + %55s  %8.5f  %8.5f\n",
                    colnames(x)[chosen[U + k - 1]], kl[k], elpd[k]))

        if (length(chosen) - U == max.num.pred)
            break
    }

    ## evaluate the full model
    full <- fit.submodel(x, w, sigma2, fit, 1:P, xt, yt, is.logistic)

    ## remove trailing zeros
    len <- length(chosen) - U + 1
    kl <- kl[1:len]
    elpd <- elpd[1:len]
    delta.elpd <- elpd - full$elpd

    res <- data.frame(var=c("Initial set of covariates",
                            colnames(x)[setdiff(chosen, 1:U)]),
                      kl=kl, elpd=elpd, delta.elpd=delta.elpd)
    if (!is.null(out.csv))
        write.csv(file=out.csv, res, row.names=FALSE)

    class(res) <- c("projsel", "data.frame")
    return(res)
}

#' Plot of the incremental contribution of each predictor
#'
#' A function of name \code{getfullname} to match variable names to full
#' names is searched on the current workspace, and if found it is used to
#' transform the labels before creating the plot.
#'
#' @param x A data frame created by \code{\link{projsel}}.
#' @param title Title of the plot. If \code{NULL}, no title is displayed.
#' @param max.labels Maximum number of points to be labelled. If \code{NULL},
#'        all those present in the \var{sel} file are displayed.
#' @param font.size Font size used to scale all text in the plot.
#' @param hadj,vadj Horizontal and vertical adjustment for the labels.
#' @param ... Other options to plot() (currently ignored).
#'
#' @return
#' A \pkg{ggplot2} object showing the relative incremental contribution of each
#' predictor starting from the initial set of unpenalized covariates.
#'
#' @import ggplot2
#' @method plot projsel
#' @export
plot.projsel <- function(x, title=NULL, max.labels=NULL, font.size=12,
                         hadj=0.05, vadj=0.02, ...) {

    ## get full variable names if possible
    sel <- x
    labs <- tryCatch(get("getfullname")(as.character(sel$var)),
                     error=function(e) as.character(sel$var))
    labs <- gsub(" \\(.*\\)$", "", labs)
    if (!is.null(max.labels)) {
        labs[-c(1:(max.labels + 1))] <- ""
    }

    geom.text.size <- font.size * 5 / 14

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
