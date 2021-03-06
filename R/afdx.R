#' afdx: Diagnosis performance indicators from attributable fraction estimates.
#' 
#' The afdx package provides functions to estimate the attributable fraction 
#' using logit exponential model or  bayesian latent class model.
#' 
#' @section The logit exponential model:
#' 
#' The logitexp function estimated the logit exponential function fitting a
#' maximum likelihood model. The senspec() function estimate the sensitivity,
#' specificity, positive predicted value and negative predicted values for
#' the specified cut-off points.
#' 
#' @section The bayesian latent class model:
#' 
#' The get_latent_model() provides an rjags model template to estimate
#' the attributable fraction and the sensitivity, specificity, positive
#' predicted value and negative predicted value of the latent class model.
#'  
#'  @docType package
#'  @name afdx
"_PACKAGE"
