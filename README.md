# Hierarchical Shrinkage Stan Models for Biomarker Selection

The **hsstan** package provides linear and logistic regression models penalized
with hierarchical shrinkage priors for selection of biomarkers. Models are
fitted with [Stan](https://mc-stan.org), which allows to perform full Bayesian
inference ([Carpenter et al. (2017)](https://doi.org/10.18637/jss.v076.i01)).

It implements the horseshoe and regularized horseshoe priors [(Piironen and
Vehtari (2017)](https://doi.org/10.1214/17-EJS1337SI)), and the projection
predictive selection approach to recover a sparse set of predictive biomarkers
([Piironen, Paasiniemi and Vehtari (2018)](https://arxiv.org/abs/1810.02406)).

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
print(sel.biom)
#                       var          kl      elpd    delta.elpd
# 1          Intercept only 0.352649015 -697.3736 -225.90204686
# 2  Unpenalized covariates 0.333616277 -682.5622 -211.09067799
# 3                     bmi 0.139286593 -542.4549  -70.98338127
# 4                     ltg 0.058546210 -493.4700  -21.99851021
# 5                     map 0.035877800 -482.8933  -11.42175893
# 6                     hdl 0.010385908 -473.8477   -2.37619296
# 7                      tc 0.005370697 -472.1645   -0.69299636
# 8                     ldl 0.002490666 -471.8571   -0.38562094
# 9                     tch 0.001153068 -471.5438   -0.07224846
# 10                    glu 0.000000000 -471.4715    0.00000000
```

### References

* [M. Colombo][mcol], E. Valo, S.J. McGurnaghan et al.,
  Biomarkers associated with progression of renal disease in type 1 diabetes,
  _Diabetologia_ (2019) 62 (9): 1616-1627.
  https://doi.org/10.1007/s00125-019-4915-0

[mcol]: https://pm2.phs.ed.ac.uk/~mcolombo/
