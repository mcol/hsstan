test_that("projsel for a model with no penalized predictors",
{
    SW({
        sel.base.00 <- projsel(hs.base, max.iters=0)
        sel.base.02 <- projsel(hs.base, max.iters=0, start.from="X2")
        sel.base.1 <- projsel(hs.base)
        sel.base.2 <- projsel(hs.base, start.from="X2")
    })

    expect_equal(nrow(sel.base.00),
                 1)
    expect_true(is.na(sel.base.00$rel.kl))
    expect_equal(sel.base.00$rel.kl.null,
                 0)

    expect_equal(nrow(sel.base.02),
                 2)
    expect_equal(sel.base.02$rel.kl[2],
                 0)
    expect_equal(sel.base.02$rel.kl.null[1],
                 0)
    expect_equal(sel.base.02$rel.kl.null[2],
                 0.38193515, tolerance=tol)

    expect_equal(nrow(sel.base.1),
                 length(hs.base$betas$unpenalized))
    expect_equal(attr(sel.base.1, "start.from"),
                 character(0))

    expect_equal(nrow(sel.base.2),
                 length(hs.base$betas$unpenalized))
    expect_equal(sel.base.2$var[2],
                 "Initial submodel")
    expect_equal(attr(sel.base.2, "start.from"),
                 "X2")
})

test_that("projsel for gaussian family",
{
    SW({
        sel.gauss <- projsel(hs.gauss, out.csv="out.csv")
    })

    expect_s3_class(sel.gauss,
                    "projsel")
    expect_equal(colnames(sel.gauss),
                 c("var", "kl", "rel.kl.null", "rel.kl", "elpd", "delta.elpd"))
    expect_equal(nrow(sel.gauss),
                 length(pen) + 2)
    expect_equal(attr(sel.gauss, "start.from"),
                 hs.gauss$model.terms$unpenalized)
    expect_equal(sel.gauss$var[1:2],
                 c("Intercept only", "Initial submodel"))
    expect_equal(sel.gauss$var[-c(1:2)],
                 paste0("X", c(9, 4, 8, 7, 5, 6, 10)))
    expect_equal(sel.gauss$kl[2],
                 0.02948656, tolerance=tol)
    expect_equal(tail(sel.gauss$kl, n=1),
                 0)
    expect_equal(sel.gauss$rel.kl.null[1],
                 0)
    expect_equal(tail(sel.gauss$rel.kl.null, n=1),
                 1)
    expect_true(is.na(sel.gauss$rel.kl[1]))
    expect_equal(sel.gauss$rel.kl[2],
                 0)
    expect_equal(tail(sel.gauss$rel.kl, n=1),
                 1)
    expect_equal(sel.gauss$elpd[2],
                 -111.370949, tolerance=tol)
    expect_equal(tail(sel.gauss$elpd, n=1),
                 -109.695740, tolerance=tol)
    expect_equal(sel.gauss$delta.elpd[2],
                 -1.67520810, tolerance=tol)
    expect_true(all(diff(sel.gauss$kl) < 0))
    expect_true(file.exists("out.csv"))
    unlink("out.csv")

    expect_s3_class(plot(sel.gauss, from.covariates=FALSE, max.points=3),
                    "ggplot")
    expect_silent(print(plot(sel.gauss, title="Test", max.labels=3)))
    unlink("Rplots.pdf")
})

test_that("projsel for binomial family",
{
    num.sel <- 4
    SW({
        sel.binom <- projsel(hs.binom, max.iters=num.sel)
    })

    expect_equal(nrow(sel.binom),
                 num.sel + 2)
    expect_equal(sel.binom$var[-c(1:2)],
                 paste0("X", c(6, 9, 5, 8)))
    expect_equal(sel.binom$kl[2],
                 0.03245156, tolerance=tol)
    expect_equal(tail(sel.binom$kl, n=1),
                 0.00478297, tolerance=tol)
    expect_equal(sel.binom$rel.kl.null[1],
                 0)
    expect_true(is.na(sel.binom$rel.kl[1]))
    expect_lt(tail(sel.binom$rel.kl, n=1),
              1)
    expect_lt(tail(sel.binom$rel.kl.null, n=1),
              1)
    expect_equal(sel.binom$elpd[2],
                 -35.3793754, tolerance=tol)
    expect_equal(tail(sel.binom$elpd, n=1),
                 -34.0411588, tolerance=tol)
    expect_equal(sel.binom$delta.elpd[2],
                 -1.39197661, tolerance=tol)
    expect_true(all(diff(sel.binom$kl) < 0))
})

test_that("projsel for a model with interaction terms",
{
    SW({
        sel.inter <- projsel(hs.inter, start.from="X1:X3")
    })
    expect_equal(attr(sel.inter, "start.from"),
                 c("X1", "X3", "X1:X3"))
})

test_that("projsel from the intercept-only model",
{
    SW({
        sel.gauss.1 <- projsel(hs.gauss, start.from=character(0))
        sel.gauss.2 <- projsel(hs.gauss, start.from="X2")
    })

    expect_equal(sel.gauss.1$var[1:2],
                 c("Intercept only", "X2"))
    expect_equal(sel.gauss.2$var[1:2],
                 c("Intercept only", "Initial submodel"))
    expect_equivalent(sel.gauss.1[-2, ],
                      sel.gauss.2[-2, ])
    expect_equal(nrow(sel.gauss.1),
                 length(c(hs.gauss$betas$unpenalized, hs.gauss$betas$penalized)))
    expect_equal(attr(sel.gauss.1, "start.from"),
                 character(0))
    expect_equal(attr(sel.gauss.2, "start.from"),
                 "X2")
})

test_that("projsel from a non-default starting submodel",
{
    SW({
        sel.gauss <- projsel(hs.gauss, start.from="X1")
    })

    expect_equal(nrow(sel.gauss),
                 length(c(grep("X1", names(hs.gauss$beta$unpenalized), inv=TRUE),
                          hs.gauss$betas$penalized)) + 1)
})

test_that("projsel from a submodel that includes all variables",
{
    SW({
        sel.gauss <- projsel(hs.gauss, start.from=paste0("X", 1:P))
    })

    expect_equal(sel.gauss$var,
                 c("Intercept only", "Initial submodel"))
    expect_equal(sel.gauss$rel.kl.null,
                 c(0, 1))
    expect_equal(sel.gauss$rel.kl,
                 c(NA, 1))
    expect_equal(sel.gauss$delta.elpd[2],
                 0)
    expect_equal(attr(sel.gauss, "row.names"),
                 1:2)
    expect_equal(attr(sel.gauss, "start.from"),
                 paste0("X", 1:P))
})

test_that("projsel for a model with few observations",
{
    SW({
        hs.small <- hsstan(df[1:5, ], y.gauss ~ X2 + X3, pen,
                           iter=iters,  chains=chains, refresh=0)
        expect_message(sel <- projsel(hs.small),
                      "Fully saturated model reached")
    })
    expect_equal(tail(sel$kl, 1),
                 0)
})

test_that("projsel for a cross-validated object",
{
    SW({
        sel.gauss <- projsel(cv.gauss$fits[[1]])
    })

    expect_equal(tail(sel.gauss$elpd, n=1),
                 sum(colMeans(log_lik(cv.gauss$fits[[1]],
                                      newdata=df[folds == 2, ]))))
})

test_that("projsel with invalid inputs",
{
    expect_error(projsel(hs.gauss, max.iters=-1),
                 "must be a non-negative integer")
    expect_error(projsel(hs.gauss, max.iters=2.3),
                 "must be a non-negative integer")
    expect_error(projsel(hs.gauss, max.iters=iris),
                 "must be a non-negative integer")
})
