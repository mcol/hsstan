test_that("projsel",
{
    SW({
        hs.nopen <- hsstan(df, mod.gauss, iter=50, chains=1, family=gaussian)
    })

    expect_error(projsel(hs.nopen),
                 "Model doesn't contain penalized predictors")
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
    expect_equal(sel.gauss$var[1:2],
                 c("Intercept only", "Unpenalized covariates"))
    expect_equal(sel.gauss$var[-c(1:2)],
                 paste0("X", c(9, 4, 8, 5, 7, 6, 10)))
    expect_equal(sel.gauss$kl[2],
                 0.024897581, tolerance=tol)
    expect_equal(tail(sel.gauss$kl, n=1),
                 0)
    expect_equal(sel.gauss$rel.kl.null[1],
                 0)
    expect_equal(tail(sel.gauss$rel.kl.null, n=1),
                 1)
    expect_true(is.na(sel.gauss$rel.kl[1]))
    expect_equal(tail(sel.gauss$rel.kl, n=1),
                 1)
    expect_equal(sel.gauss$elpd[2],
                 -111.717561, tolerance=tol)
    expect_equal(tail(sel.gauss$elpd, n=1),
                 -110.340511, tolerance=tol)
    expect_equal(sel.gauss$delta.elpd[2],
                 -1.37705019, tolerance=tol)
    expect_true(all(diff(sel.gauss$kl) < 0))
    expect_true(file.exists("out.csv"))
    unlink("out.csv")

    expect_s3_class(plot(sel.gauss),
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
                 0.054780181, tolerance=tol)
    expect_equal(tail(sel.binom$kl, n=1),
                 0.008075134, tolerance=tol)
    expect_equal(sel.binom$rel.kl.null[1],
                 0)
    expect_true(tail(sel.binom$rel.kl.null, n=1) < 1)
    expect_true(is.na(sel.binom$rel.kl[1]))
    expect_true(tail(sel.binom$rel.kl, n=1) < 1)
    expect_equal(sel.binom$elpd[2],
                 -35.2036857, tolerance=tol)
    expect_equal(tail(sel.binom$elpd, n=1),
                 -32.8004470, tolerance=tol)
    expect_equal(sel.binom$delta.elpd[2],
                 -2.53499502, tolerance=tol)
    expect_true(all(diff(sel.binom$kl) < 0))
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
