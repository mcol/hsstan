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

#' The 'hsstan' package.
#'
#' Regression models with hierarchical shrinkage priors
#'
#' @docType package
#' @name hsstan-package
#' @import Rcpp
#' @import methods
#' @useDynLib hsstan, .registration = TRUE
#'
#' @references
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.17.3. http://mc-stan.org
NULL

.onLoad <- function(libname, pkgname) { # nocov start
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
} # nocov end

#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname) {

    ## number of cores used by default for sampling from the chains
    if (.Platform$OS.type != "windows")
        options(mc.cores=min(ceiling(parallel::detectCores() / 2), 4))
    else
        options(mc.cores=1) # nocov

    packageStartupMessage("hsstan ", packageVersion("hsstan"), ":")
    if (.Platform$OS.type != "windows")
        packageStartupMessage("  - Defaulting to ", options("mc.cores"),
                              " cores: use 'options(mc.cores=<cores>)' to change")
}
