##=============================================================================
##
## Copyright (c) 2018-2019 Marco Colombo and Paul McKeigue
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

#' Hierarchical shrinkage Stan models for biomarker selection
#'
#' The **hsstan** package provides linear and logistic regression models
#' penalized with hierarchical shrinkage priors for selection of biomarkers.
#' Models are fitted with Stan (Carpenter et al. (2017)), which allows to
#' perform full Bayesian inference.
#'
#' The package implements the horseshoe and regularized horseshoe priors
#' (Piironen and Vehtari (2017)), and the projection predictive selection
#' approach to recover a sparse set of predictive biomarkers (Piironen,
#' Paasiniemi and Vehtari (2020)).
#'
#' The approach is particularly suited to selection from high-dimensional
#' panels of biomarkers, such as those that can be measured by MSMS or similar
#' technologies (Colombo, Valo, McGurnaghan et al. (2019), Colombo, McGurnaghan,
#' Blackbourn et al. (2020)).
#'
#' @docType package
#' @name hsstan-package
#' @import Rcpp
#' @import methods
#' @useDynLib hsstan, .registration = TRUE
#'
#' @references
#' B. Carpenter et al. (2017),
#' Stan: a probabilistic programming language,
#' _Journal of Statistical Software_, 76 (1).
#' \url{https://doi.org/10.18637/jss.v076.i01}
#'
#' J. Piironen and A. Vehtari (2017),
#' Sparsity information and regularization in the horseshoe and other shrinkage
#' priors, _Electronic Journal of Statistics_, 11 (2), 5018-5051.
#' \url{https://doi.org/10.1214/17-EJS1337SI}
#'
#' J. Piironen, M. Paasiniemi and A. Vehtari (2020),
#' Projective inference in high-dimensional problems: prediction and feature
#' selection, _Electronic Journal of Statistics_, 14 (1), 2155-2197.
#' \url{https://doi.org/10.1214/20-EJS1711}
#'
#' M. Colombo, E. Valo, S.J. McGurnaghan et al. (2019),
#' Biomarkers associated with progression of renal disease in type 1 diabetes,
#' _Diabetologia_, 62 (9), 1616-1627.
#' \url{https://doi.org/10.1007/s00125-019-4915-0}
#'
#' M. Colombo, S.J. McGurnaghan, L.A.K. Blackbourn et al. (2020),
#' Comparison of serum and urinary biomarker panels with albumin creatinin
#' ratio in the prediction of renal function decline in type 1 diabetes,
#' _Diabetologia_, 63 (4), 788-798.
#' \url{https://doi.org/10.1007/s00125-019-05081-8}
NULL

.onLoad <- function(libname, pkgname) { # nocov start
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
} # nocov end

.onAttach <- function(libname, pkgname) {

    ## number of cores used by default for sampling from the chains
    if (is.null(options()$mc.cores))
        options(mc.cores=min(ceiling(parallel::detectCores() / 2), 4)) # nocov

    packageStartupMessage("hsstan ", utils::packageVersion("hsstan"),
                          ": using ", options("mc.cores"),
                          " cores, set 'options(mc.cores)' to change.")
}
