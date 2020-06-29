##=============================================================================
##
## Copyright (c) 2017-2019 Marco Colombo and Paul McKeigue
##
## Parts of the code are based on https://github.com/jpiironen/rstan-varsel
## Portions copyright (c) 2015-2017 Juho Piironen
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
#' @param indproj Vector of indices of the columns of `x` that form the
#'        projection subspace.
#' @param is.logistic Set to `TRUE` for a logistic regression model.
#'
#' @importFrom stats binomial glm.fit quasibinomial
#' @keywords internal
lm_proj <- function(x, fit, sigma2, indproj, is.logistic) {

    P <- ncol(x)
    S <- ncol(fit)

    ## pick the Q columns of x that form the projection subspace
    xp <- x[, indproj, drop=FALSE] # matrix of dimension N x Q

    ## logistic regression model
    if (is.logistic) {

        ## compute the projection for each sample
        par.fun <- function(s) glm.fit(xp, fit[, s],
                                       family=quasibinomial())$coefficients
        if (.Platform$OS.type != "windows") {
            wp <- parallel::mclapply(X=1:S, mc.preschedule=TRUE, FUN=par.fun)
        } else { # windows
            cl <- parallel::makePSOCKcluster(getOption("mc.cores", 1))
            on.exit(parallel::stopCluster(cl))
            wp <- parallel::parLapply(X=1:S, cl=cl, fun=par.fun)
        }
        wp <- matrix(unlist(wp, use.names=FALSE), ncol=S)

        ## estimate the KL divergence between full and projected model
        fitp <- binomial()$linkinv(multiplyAB(xp, wp))
        kl <- 0.5 * mean(colMeans(binomial()$dev.resids(fit, fitp, 1)))
        sigma2p <- 1
    }

    ## linear regression model
    else {
        ## solve the projection equations
        wp <- solve(crossprod(xp, xp), multiplyAtB(xp, fit)) # Q x S matrix

        ## fit of the projected model
        fitp <- multiplyAB(xp, wp)

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
#' @param chosen Vector of indices of the columns of `x` in the current submodel.
#' @param is.logistic Set to `TRUE` for a logistic regression model.
#'
#' @keywords internal
choose.next <- function(x, sigma2, fit, fitp, chosen, is.logistic) {
    notchosen <- setdiff(1:ncol(x), chosen)
    par.fun <- function(idx) {
        lm_proj(x, fit, sigma2, sort(c(chosen, notchosen[idx])), FALSE)$kl
    }
    if (!is.logistic) {
        if (.Platform$OS.type != "windows") {
            kl <- parallel::mclapply(X=1:length(notchosen), mc.preschedule=TRUE,
                                     FUN=par.fun)
        } else { # windows
            cl <- parallel::makePSOCKcluster(getOption("mc.cores", 1))
            on.exit(parallel::stopCluster(cl))
            kl <- parallel::parLapply(X=1:length(notchosen), cl=cl, fun=par.fun)
        }
        idx.selected <- which.min(unlist(kl))
        return(notchosen[idx.selected])
    }

    ## score test
    yminusexp <- fit - fitp
    dinv.link <- fitp * (1 - fitp)
    U <- colMeans(multiplyAtB(yminusexp, x[, notchosen]))
    V <- colMeans(multiplyAtB(dinv.link, x[, notchosen]^2))
    idx.selected <- which.max(U^2 / V)
    return(notchosen[idx.selected])
}

#' Compute the fit of a submodel
#'
#' @param x Design matrix.
#' @param sigma2 Residual variance (1 for logistic regression).
#' @param mu Matrix of fitted values for the full model.
#' @param chosen Vector of indices of the columns of `x` in the current submodel.
#' @param xt Design matrix for test data.
#' @param yt Outcome variable for test data.
#' @param logistic Set to `TRUE` for a logistic regression model.
#'
#' @keywords internal
fit.submodel <- function(x, sigma2, mu, chosen, xt, yt, logistic) {
    ## projected parameters
    submodel <- lm_proj(x, mu, sigma2, chosen, logistic)
    eta <- multiplyAB(xt, submodel$w)
    return(c(submodel, elpd=elpd(yt, eta, submodel$sigma2, logistic)))
}

#' Expected log predictive density
#'
#' @param yt Outcome variable.
#' @param eta Matrix of linear predictors, with as many rows as the number of
#'        observations.
#' @param sigma2 Residual variance (1 for logistic regression).
#' @param logistic Set to `TRUE` for a logistic regression model.
#'
#' @noRd
elpd <- function(yt, eta, sigma2, logistic) {
    if (logistic)
        pd <- stats::dbinom(yt, 1, binomial()$linkinv(eta), log=TRUE)
    else {
        sigma <- sqrt(sigma2)
        pd <- t(cbind(sapply(1:nrow(eta),
                             function(z) stats::dnorm(yt[z], eta[z, ],
                                                      sigma, log=TRUE))))
    }
    sum(rowMeans(pd))
}

#' Forward selection minimizing KL-divergence in projection
#'
#' @param obj Object of class `hsstan`.
#' @param max.iters Maximum number of iterations (number of predictors selected)
#'        after which the selection procedure should stop.
#' @param start.from Vector of variable names to be used in the starting
#'        submodel. If `NULL` (default), selection starts from the set of
#'        unpenalized covariates if the model contains penalized predictors,
#'        otherwise selection starts from the intercept-only model.
#' @param out.csv If not `NULL`, the name of a CSV file to save the
#'        output to.
#'
#' @return
#' A data frame of class `projsel` where each row corresponds to a
#' forward-selected submodel that contains all variables listed up to that row.
#' Attribute `start.from` reports the predictors in the initial model.
#' The data frame contains the following columns:
#' \item{var}{names of the variables selected.}
#' \item{kl}{KL-divergence from the full model to the submodel.}
#' \item{rel.kl.null}{relative explanatory power of predictors starting from the
#'       intercept-only model.}
#' \item{rel.kl}{relative explanatory power of predictors starting from the
#'       initial submodel.}
#' \item{elpd}{the expected log predictive density of the submodels.}
#' \item{delta.elpd}{the difference in elpd from the full model.}
#'
#' @examples
#' \donttest{
#' \dontshow{utils::example("hsstan", echo=FALSE)}
#' \dontshow{oldopts <- options(mc.cores=2)}
#' # continued from ?hsstan
#' sel <- projsel(hs.biom, max.iters=3)
#' plot(sel)
#' \dontshow{options(oldopts)}
#' }
#'
#' @importFrom utils write.csv
#' @export
projsel <- function(obj, max.iters=30, start.from=NULL,
                    out.csv=NULL) {
    validate.hsstan(obj)
    validate.samples(obj)

    x <- xt <- validate.newdata(obj, obj$data)
    yt <- obj$data[[obj$model.terms$outcome]]
    is.logistic <- is.logistic(obj)
    sigma2 <- if (is.logistic) 1 else as.matrix(obj$stanfit, pars="sigma")^2

    ## set of variables in the initial submodel
    vsf <- validate.start.from(obj, start.from)
    start.from <- vsf$start.from
    chosen <- start.idx <- vsf$idx

    ## number of model parameters (including intercept)
    P <- sum(sapply(obj$betas, length))

    ## number of iterations to run
    K <- min(P - length(start.idx), max.iters)

    message(sprintf("%58s  %8s %11s", "Model", "KL", "ELPD"))
    report.iter <- function(msg, kl, elpd)
        message(sprintf("%58s  %8.5f  %8.5f", substr(msg, 1, 55), kl, elpd))

    ## fitted values for the full model (N x S)
    fit <- t(posterior_linpred(obj, transform=TRUE))

    ## intercept only model
    sub <- fit.submodel(x, sigma2, fit, 1, xt, yt, is.logistic)
    kl.elpd <- cbind(sub$kl, sub$elpd)
    report.iter("Intercept only", sub$kl, sub$elpd)

    ## initial submodel with the set of chosen covariates
    if (length(start.idx) > 1) {
        sub <- fit.submodel(x, sigma2, fit, chosen, xt, yt, is.logistic)
        kl.elpd <- rbind(kl.elpd, c(sub$kl, sub$elpd))
        report.iter("Initial submodel", sub$kl, sub$elpd)
    }

    ## add variables one at a time
    for (iter in seq_len(K)) {
        sel.idx <- choose.next(x, sigma2, fit, sub$fit, chosen, is.logistic)
        chosen <- c(chosen, sel.idx)

        ## evaluate current submodel according to projected parameters
        sub <- fit.submodel(x, sigma2, fit, chosen, xt, yt, is.logistic)
        kl.elpd <- rbind(kl.elpd, c(sub$kl, sub$elpd))
        report.iter(colnames(x)[sel.idx], sub$kl, sub$elpd)
    }

    ## evaluate the full model if selection stopped before reaching it
    full.elpd <- if (length(chosen) == P)
        sub$elpd else elpd(yt, t(posterior_linpred(obj)), sigma2, is.logistic)

    kl <- kl.elpd[, 1]
    rel.kl <- c(NA, if (length(kl) > 1) 1 - kl[-1] / kl[2])
    if (length(rel.kl) == 2 && is.nan(rel.kl[2]))
        rel.kl[2] <- 1

    res <- data.frame(var=c("Intercept only",
                            if (length(start.idx) > 1) "Initial submodel",
                            colnames(x)[setdiff(chosen, start.idx)]),
                      kl=kl,
                      rel.kl.null=1 - kl / kl[1],
                      rel.kl=rel.kl,
                      elpd=kl.elpd[, 2],
                      delta.elpd=kl.elpd[, 2] - full.elpd,
                      stringsAsFactors=FALSE, row.names=NULL)
    attr(res, "start.from") <- start.from
    if (!is.null(out.csv))
        write.csv(file=out.csv, res, row.names=FALSE)

    class(res) <- c("projsel", "data.frame")
    return(res)
}

#' Plot of relative explanatory power of predictors
#'
#' The function plots the relative explanatory power of each predictor in order
#' of selection. The relative explanatory power of predictors is computed
#' according to the KL divergence from the full model to each submodel, scaled
#' in such a way that the baseline model (either the null model or the model
#' containing only unpenalized covariates) is at 0, while the full model is at 1.
#'
#' @param x A data frame created by [projsel()].
#' @param title Title of the plot. If `NULL`, no title is displayed.
#' @param max.points Maximum number of predictors to be plotted. If `NULL`
#'        (default) or 0, all points are plotted.
#' @param max.labels Maximum number of predictors to be labelled. If `NULL`
#'        (default), all predictor labels present in `x` are displayed, which
#'        may result in overprinting.
#' @param from.covariates Whether the plotting should start from the unpenalized
#'        covariates (`TRUE` by default). If set to `FALSE`, the plot includes a
#'        point for the null (intercept-only) model.
#' @param font.size Size of the textual elements (labels and axes).
#' @param hadj,vadj Horizontal and vertical adjustment for the labels.
#' @param ... Currently ignored.
#'
#' @return
#' A **ggplot2** object showing the relative incremental contribution of each
#' predictor starting from the initial set of unpenalized covariates.
#'
#' @import ggplot2
#' @method plot projsel
#' @export
plot.projsel <- function(x, title=NULL, max.points=NULL, max.labels=NULL,
                         from.covariates=TRUE,
                         font.size=12, hadj=0.05, vadj=0, ...) {

    ## prepare the set of points to plot
    sel <- x
    if (!is.null(max.points) && max.points > 0)
        sel <- utils::head(sel, n=max.points + 2)
    if (!is.null(max.labels))
        sel$var[-c(1:(max.labels + 2))] <- ""
    if (from.covariates)
        sel <- sel[-1, ]
    labs <- sel$var

    ## convert from points to millimetres
    geom.text.size <- font.size * 25.4 / 72

    x <- seq(nrow(sel)) - 1
    text_idx <- x < mean(x) | x - floor(x / 2) * 2 == 1
    rel <- if (from.covariates) sel$rel.kl else sel$rel.kl.null
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
