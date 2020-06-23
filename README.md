# Hierarchical Shrinkage Stan Models for Biomarker Selection

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/hsstan)](https://cran.r-project.org/package=hsstan)
[![CRAN\_Downloads\_Badge](https://cranlogs.r-pkg.org/badges/hsstan)](https://cran.r-project.org/package=hsstan)

The **hsstan** package provides linear and logistic regression models penalized
with hierarchical shrinkage priors for selection of biomarkers. Models are
fitted with [Stan](https://mc-stan.org), which allows to perform full Bayesian
inference ([Carpenter et al. (2017)](https://doi.org/10.18637/jss.v076.i01)).

It implements the horseshoe and regularized horseshoe priors [(Piironen and
Vehtari (2017)](https://doi.org/10.1214/17-EJS1337SI)), and the projection
predictive selection approach to recover a sparse set of predictive biomarkers
([Piironen, Paasiniemi and Vehtari (2020)](https://doi.org/10.1214/20-EJS1711)).

The approach is particularly suited to selection from high-dimensional panels
of biomarkers, such as those that can be measured by MSMS or similar technologies.

### Example

```r
library(hsstan)
data(diabetes)

## allow using as many cores as cross-validation folds
options(mc.cores=10)

## baseline model with only clinical covariates
hs.base <- hsstan(diabetes, Y ~ age + sex)

## model with additional predictors
hs.biom <- hsstan(diabetes, Y ~ age + sex, penalized=colnames(diabetes)[3:10])
print(hs.biom)
#              mean   sd  2.5% 97.5% n_eff Rhat
# (Intercept)  0.00 0.03 -0.07  0.06  5535    1
# age          0.00 0.04 -0.07  0.08  5618    1
# sex         -0.15 0.04 -0.22 -0.08  5161    1
# bmi          0.33 0.04  0.25  0.41  4230    1
# map          0.20 0.04  0.12  0.27  3725    1
# tc          -0.45 0.26 -0.95  0.05  3283    1
# ldl          0.28 0.21 -0.13  0.69  3283    1
# hdl          0.01 0.13 -0.23  0.27  3532    1
# tch          0.07 0.08 -0.06  0.25  4062    1
# ltg          0.43 0.11  0.22  0.65  3377    1
# glu          0.02 0.03 -0.03  0.10  3313    1

## run cross-validation
set.seed(1)
folds <- caret::createFolds(diabetes$Y, k=10, list=FALSE)
cv.base <- kfold(hs.base, folds=folds)
cv.biom <- kfold(hs.biom, folds=folds)

## cross-validated performance
round(posterior_performance(cv.base), 2)
#        mean   sd    2.5%   97.5%
# r2     0.02 0.00    0.01    0.03
# llk -623.14 1.67 -626.61 -620.13
# attr(,"type")
# [1] "cross-validated"

round(posterior_performance(cv.biom), 2)
#        mean   sd    2.5%   97.5%
# r2     0.48 0.01    0.47    0.50
# llk -483.18 3.75 -490.51 -476.24
# attr(,"type")
# [1] "cross-validated"

## projection predictive selection
sel.biom <- projsel(hs.biom)
print(sel.biom, digits=4)
#                       var       kl rel.kl.null rel.kl   elpd delta.elpd
# 1          Intercept only 0.352649     0.00000     NA -697.4 -225.90205
# 2  Unpenalized covariates 0.333616     0.05397 0.0000 -682.6 -211.09068
# 3                     bmi 0.139287     0.60503 0.5825 -542.5  -70.98338
# 4                     ltg 0.058546     0.83398 0.8245 -493.5  -21.99851
# 5                     map 0.035878     0.89826 0.8925 -482.9  -11.42176
# 6                     hdl 0.010386     0.97055 0.9689 -473.8   -2.37619
# 7                      tc 0.005371     0.98477 0.9839 -472.2   -0.69300
# 8                     ldl 0.002491     0.99294 0.9925 -471.9   -0.38562
# 9                     tch 0.001153     0.99673 0.9965 -471.5   -0.07225
# 10                    glu 0.000000     1.00000 1.0000 -471.5    0.00000
```

### References

* [M. Colombo][mcol], S.J. McGurnaghan, L.A.K. Blackbourn et al.,
  Comparison of serum and urinary biomarker panels with albumin creatinin
  ratio in the prediction of renal function decline in type 1 diabetes,
  _Diabetologia_ (2020) 63 (4): 788-798.
  https://doi.org/10.1007/s00125-019-05081-8

* [M. Colombo][mcol], E. Valo, S.J. McGurnaghan et al.,
  Biomarkers associated with progression of renal disease in type 1 diabetes,
  _Diabetologia_ (2019) 62 (9): 1616-1627.
  https://doi.org/10.1007/s00125-019-4915-0

[mcol]: https://pm2.phs.ed.ac.uk/~mcolombo/
