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

## if possible, allow using as many cores as cross-validation folds
options(mc.cores=10)

## baseline model with only clinical covariates
hs.base <- hsstan(diabetes, Y ~ age + sex)

## model with additional predictors
hs.biom <- hsstan(diabetes, Y ~ age + sex, penalized=colnames(diabetes)[3:10])
print(hs.biom)
#              mean   sd  2.5% 97.5% n_eff Rhat
# (Intercept)  0.00 0.03 -0.07  0.07  4483    1
# age          0.00 0.04 -0.07  0.08  4706    1
# sex         -0.15 0.04 -0.22 -0.08  5148    1
# bmi          0.33 0.04  0.25  0.41  4228    1
# map          0.20 0.04  0.12  0.28  3571    1
# tc          -0.45 0.25 -0.94  0.04  3713    1
# ldl          0.28 0.20 -0.12  0.68  3674    1
# hdl          0.01 0.12 -0.23  0.25  3761    1
# tch          0.07 0.08 -0.06  0.25  4358    1
# ltg          0.43 0.11  0.22  0.64  3690    1
# glu          0.02 0.03 -0.03  0.10  3034    1

## behaviour of the sampler
sampler.stats(hs.base)
#         accept.stat stepsize divergences treedepth gradients warmup sample
# chain:1      0.9497   0.5723           0         3      6320   0.09   0.08
# chain:2      0.9357   0.6480           0         3      5938   0.09   0.08
# chain:3      0.9455   0.6014           0         3      6112   0.09   0.08
# chain:4      0.9488   0.5932           0         3      6238   0.09   0.08
# all          0.9449   0.6037           0         3     24608   0.36   0.32

sampler.stats(hs.biom)
#         accept.stat stepsize divergences treedepth gradients warmup sample
# chain:1      0.9821   0.0191           0         8    233656   5.04   4.28
# chain:2      0.9891   0.0158           1         8    255994   5.88   4.72
# chain:3      0.9908   0.0143           0         9    274328   5.77   5.14
# chain:4      0.9933   0.0121           0         9    344984   5.98   6.70
# all          0.9888   0.0153           1         9   1108962  22.67  20.84

## approximate leave-one-out cross-validation with Pareto smoothed
## importance sampling
loo(hs.base)
# Computed from 4000 by 442 log-likelihood matrix
#          Estimate   SE
# elpd_loo   -622.4 11.4
# p_loo         3.4  0.2
# looic      1244.9 22.7
# ------
# Monte Carlo SE of elpd_loo is 0.0.
#
# All Pareto k estimates are good (k < 0.5).

loo(hs.biom)
# Computed from 4000 by 442 log-likelihood matrix
#          Estimate   SE
# elpd_loo   -476.5 13.7
# p_loo         9.8  0.7
# looic       953.0 27.5
# ------
# Monte Carlo SE of elpd_loo is 0.1.
#
# All Pareto k estimates are good (k < 0.5).

## run 10-folds cross-validation
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
# llk -482.86 3.76 -490.45 -476.56
# attr(,"type")
# [1] "cross-validated"

## projection predictive selection
sel.biom <- projsel(hs.biom)
print(sel.biom, digits=4)
#                 var       kl rel.kl.null rel.kl   elpd delta.elpd
# 1    Intercept only 0.352283     0.00000     NA -627.3 -155.84260
# 2  Initial submodel 0.333156     0.05429 0.0000 -619.8 -148.39729
# 3               bmi 0.138629     0.60648 0.5839 -533.1  -61.69199
# 4               ltg 0.058441     0.83411 0.8246 -492.5  -21.09681
# 5               map 0.035970     0.89789 0.8920 -482.7  -11.25515
# 6               hdl 0.010304     0.97075 0.9691 -473.9   -2.41192
# 7                tc 0.005292     0.98498 0.9841 -472.2   -0.72490
# 8               ldl 0.002444     0.99306 0.9927 -471.8   -0.38292
# 9               tch 0.001105     0.99686 0.9967 -471.5   -0.07819
# 10              glu 0.000000     1.00000 1.0000 -471.4    0.00000
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
