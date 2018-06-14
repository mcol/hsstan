##=============================================================================
##
## Copyright (c) 2017-2018 Marco Colombo and Paul McKeigue
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


#' Computes projections of full predictors on to subspace of predictors
#'
#' @param x Design matrix.
#' @param fit Matrix of fitted values for the full model.
#' @param sigma2 Residual variance (1 for logistic regression).
#' @param indproj Vector of indices of the columns of \var{x} that form the
#'        projection subspace.
#' @param is.logistic Set to \code{TRUE} for a logistic regression model.
#'
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
        wp <- foreach(s=1:S, .combine=cbind) %dopar% {
            ## XXX: what about the intercept?
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

#' Returns the next variable to be added to the current submodel
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
        kl <- foreach(i=1:length(notchosen), .combine=c) %dopar% {
            ind <- sort(c(chosen, notchosen[i]))
            lm_proj(x, fit, sigma2, ind, FALSE)$kl
        }
        idx.selected <- which.min(kl)
        return(notchosen[idx.selected])
        sigma2 <- sigma2 + colMeans((fit - fitp)^2)
        dinv.link <- matrix(1, ncol(fit), nrow(fit))
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
#' @importMethodsFrom rstan extract
#' @export
lm_fprojsel <- function(samples, max.num.pred=30, out.csv=NULL) {

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
        mlpd <- mean(log(rowMeans(pd)))

        return(list(fit=submodel$fit, kl=submodel$kl, mlpd=mlpd))
    }

    ## check that the model contains penalized predictors
    if (is.null(samples$betas$penalized)) {
        stop("Model doesn't contain penalized predictors.")
    }

    x <- samples$data[samples$train, ]
    xt <- samples$data[samples$test, ]
    yt <- samples$y[samples$test]
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
    kl <- mlpd <- rep(0, P - U + 1)
    cat(sprintf("%58s  %8s %9s\n", "Model", "KL", "MLPD"))

    ## start from the model having only unpenalized variables
    chosen <- 1:U
    notchosen <- setdiff(1:P, chosen)
    sub <- fit.submodel(x, w, sigma2, fit, chosen, xt, yt, is.logistic)
    fitp <- sub$fit
    kl[1] <- sub$kl
    mlpd[1] <- sub$mlpd
    cat(sprintf("%58s  %8.5f  %8.5f\n",
                "Initial set of covariates", kl[1], mlpd[1]))

    ## add variables one at a time
    for (k in 2:(P - U + 1)) {
        sel.idx <- choose.next(x, sigma2, fit, fitp, chosen, is.logistic)
        chosen <- c(chosen, sel.idx)

        ## evaluate current submodel according to projected parameters
        sub <- fit.submodel(x, w, sigma2, fit, chosen, xt, yt, is.logistic)
        fitp <- sub$fit
        kl[k] <- sub$kl
        mlpd[k] <- sub$mlpd
        cat(sprintf(" + %55s  %8.5f  %8.5f\n",
                    colnames(x)[chosen[U + k - 1]], kl[k], mlpd[k]))

        if (length(chosen) - U == max.num.pred)
            break
    }

    ## evaluate the full model
    full <- fit.submodel(x, w, sigma2, fit, 1:P, xt, yt, is.logistic)

    ## remove trailing zeros
    len <- length(chosen) - U + 1
    kl <- kl[1:len]
    mlpd <- mlpd[1:len]
    delta.mlpd <- mlpd - full$mlpd

    res <- data.frame(var=c("Initial set of covariates",
                            colnames(x)[setdiff(chosen, 1:U)]),
                      kl=kl, mlpd=mlpd, delta.mlpd=delta.mlpd)
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
#' @param sel A data.frame created by \code{\link{lm_fprojsel}}.
#' @param title Title of the plot. If \code{NULL}, no title is displayed.
#' @param filename Optional name of a file where the plot is saved in png
#'        format. If \code{NULL}, no png file is produced.
#' @param max.labels Maximum number of points to be labelled. If \code{NULL},
#'        all those present in the \var{sel} file are displayed.
#' @param font.size Font size used to scale all text in the plot.
#' @param width Width of the png file in pixels.
#' @param height Height of the png file in pixels.
#' @param ... Other options to plot() (currently ignored).
#'
#' @import ggplot2
#' @method plot projsel
#' @export
plot.projsel <- function(sel, title=NULL, filename=NULL, max.labels=NULL,
                         font.size=12, width=800, height=800, ...) {

    ## get full variable names if possible
    labs <- tryCatch(get("getfullname")(as.character(sel$var)),
                     error=function(e) as.character(sel$var))
    labs <- gsub(" \\(.*\\)$", "", labs)
    if (!is.null(max.labels)) {
        labs[-c(1:(max.labels + 1))] <- ""
    }

    geom.text.size <- font.size * 5 / 14

    sel$rel <- 1 - sel$kl / sel$kl[1]
    x <- seq(nrow(sel)) - 1
    p <- ggplot(data=sel, aes(x=x, y=rel, label=labs)) +
      coord_cartesian(ylim=range(c(0, 1))) +
      geom_line() + geom_point(size=geom.text.size / 3) +
      geom_text(aes(x=x + ifelse(x < mean(x), 0.3, -0.3)),
                size=geom.text.size,
                hjust=ifelse(x < mean(x), "left", "right")) +
      xlab("Number of biomarkers") +
      ylab("Relative explanatory power") +
      theme(text=element_text(size=font.size))

    if (!is.null(title))
        p <- p + ggtitle(title)

    if (!is.null(filename)) {
        png(filename, width=width, height=height)
        print(p)
        dev.off()
    }
    else {
        ## disable clipping of the text labels
        p <- ggplot_gtable(ggplot_build(p))
        p$layout$clip[p$layout$name == "panel"] <- "off"
        p
    }
}
