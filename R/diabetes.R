##=============================================================================
##
## Copyright (c) 2019 Marco Colombo
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


#' Diabetes data with interaction terms
#'
#' The dataset consists of observations on 442 individuals for which a
#' quantitative measure of diabetes progression is recorded in variable `Y`.
#' Predictors include 10 baseline measurements, 45 interactions and 9 quadratic
#' terms, for a total of 64 variables for each individual. Each variable has
#' been standardized by subtracting the mean and then dividing it by its
#' standard deviation.
#'
#' @format
#' A data frame with 442 rows and 65 columns (centred and scaled).
#'
#' @source
#' B. Efron, T. Hastie, I. Johnstone and R. Tibshirani (2004),
#' Least angle regression, _The Annals of Statistics_, 32 (2), 407-499.
#' \url{https://doi.org/10.1214/009053604000000067}
#'
#' The original dataset is available from
#' \url{https://www.stanford.edu/~hastie/Papers/LARS/data64.txt}
#'
#' @examples
#' data(diabetes, package="hsstan")
#'
#' @name diabetes
#' @docType data
#' @keywords datasets
NULL
