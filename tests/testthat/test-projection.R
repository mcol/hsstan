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
                 c("var", "kl", "elpd", "delta.elpd"))
    expect_equal(nrow(sel.gauss),
                 length(pen) + 1)
    expect_equal(sel.gauss$var[1],
                 "Unpenalized covariates")
    expect_equal(sel.gauss$var[-1],
                 paste0("X", c(9, 4, 8, 5, 7, 6, 10)))
    expect_equal(sel.gauss$kl[1],
                 0.024897581, tolerance=tol)
    expect_equal(sel.gauss$kl[length(pen) + 1],
                 0)
    expect_equal(sel.gauss$elpd[1],
                 -111.717561, tolerance=tol)
    expect_equal(sel.gauss$elpd[length(pen) + 1],
                 -110.340511, tolerance=tol)
    expect_equal(sel.gauss$delta.elpd[1],
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
        sel.binom <- projsel(hs.binom, max.num.pred=num.sel)
    })

    expect_equal(nrow(sel.binom),
                 num.sel + 1)
    expect_equal(sel.binom$var[-1],
                 paste0("X", c(6, 9, 5, 8)))
    expect_equal(sel.binom$kl[1],
                 0.054780181, tolerance=tol)
    expect_equal(sel.binom$kl[num.sel + 1],
                 0.008075134, tolerance=tol)
    expect_equal(sel.binom$elpd[1],
                 -35.2036857, tolerance=tol)
    expect_equal(sel.binom$elpd[num.sel + 1],
                 -32.8004470, tolerance=tol)
    expect_equal(sel.binom$delta.elpd[1],
                 -2.53499502, tolerance=tol)
    expect_true(all(diff(sel.binom$kl) < 0))
})

test_that("projsel for a cross-validated object",
{
    SW({
        sel.gauss <- projsel(cv.gauss$fits[[1]])
    })

    expect_equal(sel.gauss$elpd[length(pen) + 1],
                 sum(colMeans(log_lik(cv.gauss$fits[[1]],
                                      newdata=df[folds == 2, ]))))
})
