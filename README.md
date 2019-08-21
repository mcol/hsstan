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

### References:

* [M. Colombo][mcol], E. Valo, S.J. McGurnaghan et al.,
  Biomarkers associated with progression of renal disease in type 1 diabetes,
  _Diabetologia_ (2019) 62 (9): 1616-1627.
  https://doi.org/10.1007/s00125-019-4915-0

[mcol]: https://pm2.phs.ed.ac.uk/~mcolombo/
