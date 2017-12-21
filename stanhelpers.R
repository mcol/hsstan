library(doMC)
require(rstan)
require(coda)
library(ggplot2)
library(pROC)

summarize.params <- function(samples, pars, varnames) { # returns posterior means and 95% CIs
    post.params <- As.mcmc.list(object=samples, pars)
    params.summary <- summary(post.params)
    params.summary <- data.frame(mean=params.summary$statistics[, 1],
                                 params.summary$quantiles[, c(1, 5)])
    colnames(params.summary) <- c("mean", "centile2.5", "centile97.5")
    rownames(params.summary) <- varnames
    return(params.summary)
}

posterior.means <- function(samples, varnames) { # returns posterior means
    unlist(lapply(extract(samples, pars=varnames),
                  function(z) if (is.matrix(z)) colMeans(z) else mean(z)))
}

## returns the posterior means of the coefficients
get.coefficients <- function(samples, coeff.names) {
    beta.u <- posterior.means(samples, "beta_u")
    beta.p <- posterior.means(samples, "beta_p")
    stopifnot(length(c(beta.u, beta.p)) == length(coeff.names))
    u.idx <- 1:length(beta.u)
    names(beta.u) <- coeff.names[u.idx]
    names(beta.p) <- coeff.names[-u.idx]

    return(list(unpenalized=beta.u, penalized=beta.p))
}

## returns the posterior means of the unpenalized coefficients
get.unpenalized.coefficients <- function(samples, coeff.names) {
    beta.u <- posterior.means(samples, "beta_u")
    stopifnot(length(beta.u) == length(coeff.names))
    names(beta.u) <- coeff.names

    return(list(unpenalized=beta.u))
}

## runs a stan model (either with hamiltonian montecarlo or variational bayes)
## over the cross-validation folds
sample.stan.cv <- function(stan.file, x, y, covariates, biomarkers, folds,
                           standardize=TRUE,
                           nu=3, model.type=c("mc", "vb")) {

    stopifnot(nrow(x) == length(y))
    stopifnot(all(covariates %in% colnames(x)))
    stopifnot(all(biomarkers %in% colnames(x)))
    model.type <- match.arg(model.type)

    ## create the design matrix
    model <- paste("~", paste(c(covariates, biomarkers), collapse=" + "))
    X <- model.matrix(as.formula(model), data=x)
    P <- ncol(X)
    U <- P - length(biomarkers)
    which.unpenalized <- 1:U
    which.penalized <- setdiff(1:P, which.unpenalized)
    X <- X[, c(which.unpenalized, which.penalized)]
    N <- nrow(X)
    num.folds <- length(folds)

    ## names of variables to be standardized
    stand.cols <- colnames(x)[!sapply(x, class) %in% c("factor", "character")]

    ## compile the model
    model <- stan_model(file=stan.file)

    ## cross-validation
    cv <- foreach(fold=1:num.folds) %dopar% {
        test <- 1:N %in% folds[[fold]]
        train <- !test
        X_train <- X[train, ]
        y_train <- y[train]
        N_train <- nrow(X_train)
        X_test <- X[test, ]
        y_test <- y[test]
        N_test <- nrow(X_test)

        ## standardize continuous outcome
        if (length(unique(y_train)) > 2) {
            train.mean <- mean(y_train, na.rm=TRUE)
            train.sd <- sd(y_train, na.rm=TRUE)
            y_train <- (y_train - train.mean) / train.sd
            y_test <- (y_test - train.mean) / train.sd
        }

        ## standardize all columns corresponding to numerical variables: this
        ## excludes those generated from factor/character variables and the
        ## intercept column
        if (standardize) {
            train.mean <- colMeans(X_train, na.rm=TRUE)
            train.sd <- apply(X_train, 2, sd, na.rm=TRUE)
            train.mean[!colnames(X_train) %in% stand.cols] <- 0
            train.sd[!colnames(X_train) %in% stand.cols] <- 1
            X_train <- ((X_train - rep(train.mean, each=N_train))
                         %*% diag(1 / train.sd))
            X_test <- ((X_test - rep(train.mean, each=N_test))
                        %*% diag(1 / train.sd))
        }

        ## parameters not used by a model are ignored
        data.input <- list(N_train=N_train, N_test=N_test,
                           y_train=y_train, y_test=y_test,
                           X_train=X_train, X_test=X_test,
                           P=P, U=U, nu=nu)

        if (model.type == "mc") {
            samples <- sampling(model, data=data.input,
                                chains=4, iter=1000, warmup=500,
                                seed=123, control=list(adapt_delta=0.8))
        }
        else {
            samples <- vb(model, data=data.input,
                          iter=50000, output_samples=2000,
                          init="random", algorithm="meanfield")
        }

        if (P == U) {
            betas <- get.unpenalized.coefficients(samples, colnames(X))
        } else {
            betas <- get.coefficients(samples, colnames(X))
        }

        list(samples=samples, betas=betas,
             X_train=X_train, X_test=X_test,
             y_train=y_train, y_test=y_test, train=train, test=test)
    }

    return(cv)
}

## runs a stan model (either with hamiltonian montecarlo or variational bayes)
sample.stan <- function(stan.file, x, y, covariates, biomarkers,
                        standardize=TRUE,
                        nu=3, model.type=c("mc", "vb")) {

    stopifnot(nrow(x) == length(y))
    stopifnot(all(covariates %in% colnames(x)))
    stopifnot(all(biomarkers %in% colnames(x)))
    model.type <- match.arg(model.type)

    ## create the design matrix
    model <- paste("~", paste(c(covariates, biomarkers), collapse=" + "))
    X <- model.matrix(as.formula(model), data=x)
    N <- nrow(X)
    P <- ncol(X)
    U <- P - length(biomarkers)
    which.unpenalized <- 1:U
    which.penalized <- setdiff(1:P, which.unpenalized)
    X <- X[, c(which.unpenalized, which.penalized)]

    ## standardize continuous outcome
    if (length(unique(y)) > 2) {
        y <- as.numeric(scale(y))
    }

    ## standardize all columns corresponding to numerical variables: this
    ## excludes those generated from factor/character variables as well as
    ## the intercept column
    if (standardize) {
        stand.cols <- colnames(x)[!sapply(x, class) %in% c("factor", "character")]
        stand.idx <- which(colnames(X) %in% stand.cols)
        X[, stand.idx] <- scale(X[, stand.idx])
    }

    ## use all available data for both training and testing: this effectively
    ## computes the fit of the model (y_pred) for all observations
    train <- test <- rep(TRUE, N)

    ## parameters not used by a model are ignored
    data.input <- list(N_train=N, N_test=N,
                       y_train=y, y_test=y,
                       X_train=X, X_test=X,
                       P=P, U=U, nu=nu)

    ## compile the model
    cat("Compiling STAN model", stan.file, "...")
    model <- stan_model(file=stan.file)
    if (model.type == "mc") {
        samples <- sampling(model, data=data.input, iter=1000, warmup=500,
                            chains=4, seed=123, control=list(adapt_delta=0.8))
    }
    else {
        samples <- vb(model, data=data.input, iter=50000, output_samples=2000,
                      init="random", algorithm="meanfield")
    }

    if (P == U) {
        betas <- get.unpenalized.coefficients(samples, colnames(X))
    } else {
        betas <- get.coefficients(samples, colnames(X))
    }

    return(list(samples=samples, betas=betas, train=train, test=test,
                data=X, y=y))
}

## returns projections of full predictors on to subspace of predictors,
## KL divergence between full predictor and projection
lm_proj <- function(x, fit, sigma2, indproj, is.logistic) {
    ## x is model matrix of dimension N x P, including column of ones,
    ## fit is the N x S matrix of fitted values for full the model,
    ## sigma2 is residual variance (1 for logistic regression)
    ## indproj is vector of columns of x that form the projection subspace
    ##
    ## returns the projected samples and estimated kl-divergence

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

## return the index of the variable that should be added to the current model
## according to the smallest kl-divergence (linear regression) or the largest
## score test (logistic regression)
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

## forward selection minimizing KL-divergence in projection
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

    x <- samples$data[samples$train, ]
    xt <- samples$data[samples$test, ]
    yt <- samples$y[samples$test]
    samples <- samples$samples

    beta.samples <- extract(samples, pars=c("beta_u", "beta_p"))
    w <- t(cbind(beta.samples$beta_u, beta.samples$beta_p)) # P x S
    sigma2 <- tryCatch(unlist(extract(samples, pars=c("sigma")))^2,
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
    return(res)
}

## plot of the incremental contribution of each biomarker
plot.fprojsel <- function(sel, title, filename=NULL, max.labels=NULL,
                          font.size=12, width=800, height=800) {

    ## get full variable names if possible
    labs <- tryCatch(getfullname(as.character(sel$var)),
                     error=function(e) as.character(sel$var))
    labs <- gsub(" \\(.*\\)$", "", labs)
    if (!is.null(max.labels)) {
        labs[-c(1:max.labels)] <- ""
    }

    geom.text.size <- font.size * 5 / 14

    x <- seq(nrow(sel)) - 1
    p <- ggplot(data=sel, aes(x=x, y=kl, label=labs)) +
      coord_cartesian(ylim=range(c(0, sel$kl))) +
      scale_y_continuous(sec.axis=sec_axis(~ 1 - . / sel$kl[1],
                                           name="Relative explanatory power")) +
      geom_hline(yintercept=0, linetype=2) +
      geom_line() + geom_point(size=geom.text.size / 3) +
      geom_text(aes(x=x + ifelse(x < mean(x), 0.3, -0.3)),
                size=geom.text.size,
                hjust=ifelse(x < mean(x), "left", "right")) +
      ggtitle(title) +
      xlab("Number of biomarkers") +
      ylab("KL divergence (nats)") +
      theme(text=element_text(size=font.size))

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

## extract measures of performance from the cross-validation results
get.cv.performance <- function(hs.cv, base.cv, out.csv=NULL) {

    gaussian.llk <- function(y.pred, y.obs, disp)
        -0.5 * sum(log(disp) + (y.pred - y.obs)^2 / disp)
    binomial.llk <- function(y.pred, y.obs)
        sum(log(y.pred[y.obs == 1])) + sum(log(1 - y.pred[y.obs == 0]))
    r2 <- function(y.pred, y.obs) {
        corr <- cor(y.pred, y.obs)
        if (corr < 0)
            return(0)
        return(corr^2)
    }
    to.prob <- function(lin.pred) 1 / (1 + exp(-lin.pred))
    auc <- function(y.pred, y.obs) as.numeric(roc(y.obs, y.pred)$auc)
    loglik.ratio <- function(y.pred, y.obs, prop.cases)
        (2 * y.obs - 1) * (log(y.pred) - log(1 - y.pred)) -
        (2 * y.obs - 1) * (log(prop.cases) - log(1 - prop.cases))

    y.obs.all <- y.pred.hs.all <- y.pred.base.all <- NULL
    sigma.hs.all <- sigma.base.all <- NULL
    num.folds <- length(hs.cv)
    llk <- perf <- rep(NA, num.folds)
    llk.ratio <- llk.ratio.var <- rep(NA, num.folds)
    llk.base <- perf.base <- rep(NA, num.folds)
    llk.ratio.base <- llk.ratio.var.base <- rep(NA, num.folds)

    ## loop over the folds
    for (fold in 1:num.folds) {
        y.obs <- hs.cv[[fold]]$y_test

        y.pred.hs <- posterior.means(hs.cv[[fold]]$samples, "y_pred")
        sigma.hs <- tryCatch(summary(As.mcmc.list(hs.cv[[fold]]$samples,
                                                  pars="sigma"))$statistics[1],
                             error=function(e) return(NA))

        y.pred.base <- posterior.means(base.cv[[fold]]$samples, "y_pred")
        sigma.base <- tryCatch(summary(As.mcmc.list(base.cv[[fold]]$samples,
                                                    pars="sigma"))$statistics[1],
                               error=function(e) return(NA))

        ## logistic regression
        if (is.na(sigma.hs)) {
            y.pred.hs <- to.prob(y.pred.hs)
            llk[fold] <- binomial.llk(y.pred.hs, y.obs)
            perf[fold] <- auc(y.pred.hs, y.obs)
            prop.cases <- sum(y.obs) / sum(hs.cv[[fold]]$train)
            llkr <- loglik.ratio(y.pred.hs, y.obs, prop.cases)
            llk.ratio[fold] <- mean(llkr)
            llk.ratio.var[fold] <- var(llkr)

            y.pred.base <- to.prob(y.pred.base)
            llk.base[fold] <- binomial.llk(y.pred.base, y.obs)
            perf.base[fold] <- auc(y.pred.base, y.obs)
            llkr.base <- loglik.ratio(y.pred.base, y.obs, prop.cases)
            llk.ratio.base[fold] <- mean(llkr.base)
            llk.ratio.var.base[fold] <- var(llkr.base)
        }

        ## linear regression
        else {
            llk[fold] <- gaussian.llk(y.pred.hs, y.obs, sigma.hs)
            perf[fold] <- r2(y.pred.hs, y.obs)

            llk.base[fold] <- gaussian.llk(y.pred.base, y.obs, sigma.base)
            perf.base[fold] <- r2(y.pred.base, y.obs)
        }

        sigma.hs.all <- c(sigma.hs.all, sigma.hs)
        y.obs.all <- c(y.obs.all, y.obs)
        y.pred.hs.all <- c(y.pred.hs.all, y.pred.hs)
        sigma.base.all <- c(sigma.base.all, sigma.base)
        y.pred.base.all <- c(y.pred.base.all, y.pred.base)
        cat(".")
    }
    cat("\n")
    set <- paste("Fold", 1:num.folds)

    ## compute log-likelihood and performance measure of the full model
    ## on the full vector of withdrawn observations
    set <- c(set, "Overall")
    llk <- c(llk, ifelse(is.na(sigma.hs),
                         binomial.llk(y.pred.hs.all, y.obs.all),
                         gaussian.llk(y.pred.hs.all, y.obs.all,
                                      mean(sigma.hs.all))))
    perf <- c(perf, ifelse(is.na(sigma.hs),
                           auc(y.pred.hs.all, y.obs.all),
                           r2(y.pred.hs.all, y.obs.all)))

    ## compute log-likelihood and performance measure of the baseline model
    ## on the full vector of withdrawn observations
    llk.base <- c(llk.base, ifelse(is.na(sigma.base),
                                   binomial.llk(y.pred.base.all, y.obs.all),
                                   gaussian.llk(y.pred.base.all, y.obs.all,
                                                mean(sigma.base.all))))
    perf.base <- c(perf.base, ifelse(is.na(sigma.base),
                                     auc(y.pred.base.all, y.obs.all),
                                     r2(y.pred.base.all, y.obs.all)))

    res <- data.frame(set=set, test.llk=llk, perf=perf,
                      test.llk.base=llk.base, perf.base=perf.base)
    colnames(res)[c(3, 5)] <- gsub("perf", ifelse(is.na(sigma.hs), "auc", "r2"),
                                   colnames(res)[c(3, 5)])
    if (is.na(sigma.hs)) {
        prop.cases <- sum(y.obs.all == 1) / length(y.obs.all)
        llkr <- loglik.ratio(y.pred.hs.all, y.obs.all, prop.cases)
        llkr.base <- loglik.ratio(y.pred.base.all, y.obs.all, prop.cases)
        res$llk.ratio <- c(llk.ratio, mean(llkr))
        res$llk.ratio.var <- c(llk.ratio.var, var(llkr))
        res$llk.ratio.base <- c(llk.ratio.base, mean(llkr.base))
        res$llk.ratio.var.base <- c(llk.ratio.var.base, var(llkr.base))
    }

    if (!is.null(out.csv))
        write.csv(file=out.csv, res, row.names=FALSE)
    return(res)
}
