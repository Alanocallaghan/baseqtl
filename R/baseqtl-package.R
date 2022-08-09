#' The 'baseqtl' package.
#'
#' @description Bayesian allele specific expression for eqtl analysis
#'
#' @docType package
#' @name baseqtl-package
#' @aliases baseqtl
#' @useDynLib baseqtl, .registration = TRUE
#' @import data.table
#' @import GUESSFM
#' @import addstrings
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom stats dt end fisher.test line qlogis reshape sd setNames start var
#' @importFrom utils capture.output head write.table
#'
#' @references
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.18.2. https://mc-stan.org
#'
NULL
