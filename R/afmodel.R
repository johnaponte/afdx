## Helper function to estimate maximum likelihood of model type 3
## Setup optimization by Orvalho Augusto
## 20200502 by JJAV
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


#' Exponential logit model for two variables
#' 
#' Fit a logit model of v.density on v.fever v.density with
#' a exponential coefficient for the v.density
#' 
#' logit(v.fever) ~ beta * (v. density ^ tau)
#' 
#' This corresponds to the  model 3 describe by 
#' Smith, T., Schellenberg, J.A., Hayes, R., 1994. 
#' Attributable fraction estimates and case definitions for malaria 
#' in endemic areas. Stat Med 13, 2345–2358.
#'
#' @param v.fever numeric vector of 0/1 indicating fever or equivalent
#' @param v.density numeric vector of values >= 0 indicating the density
#' 
#' @return S3 object of class afmodel with 4 components: data, model, coefficients 
#'        and the estimated attributable fraction.
#' 
#' @import maxLik
#' @importFrom stats AIC
#' @importFrom stats BIC
#' @importFrom stats binomial
#' @importFrom stats coef
#' @importFrom stats glm
#' @importFrom stats logLik
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats vcov
#' @importFrom magrittr `%>%`
#' @importFrom dplyr filter
#' @importFrom dplyr summarise
#' @export
logitexp <- function(v.fever, v.density) {
  
  stopifnot(v.fever %in% c(0,1))
  stopifnot(v.density >= 0)
  stopifnot(length(v.fever) == length(v.density))
  
  df <- data.frame(v.fever = v.fever, v.density = v.density)
  
 
  llk.af <- function(param, x = NA, y = NA) {
    alpha <- param[1]
    beta <- param[2]
    tau <- param[3]
    xb <- alpha + beta * ((x)^(tau))
    -sum(y * log(1 + exp(-xb)) + (1 - y) * log(1+exp(xb))) 
  }
  

  ls.result <-suppressWarnings(
    glm(v.fever ~ v.density, data = df, family = binomial(link = "logit")))
 
  stval <- c(coef(ls.result), 0.001)   
  names(stval) <- c("alpha", "beta", "tau")
  

  logit.result.opt <-
    maxLik::maxLik(
      start = stval,
      logLik = llk.af,
      method = "BFGS",
      x = df$v.density,
      y = df$v.fever
    )
  
  coefs <- coef(logit.result.opt)
  tau <-  coefs["tau"]
  names(tau) <- "tau"
  

  df$or = exp(coefs["beta"]*df$v.density ^tau)
  df$rr = (df$or-1)/df$or
  
  af <- 
    df %>% 
    dplyr::filter(v.fever == 1) %>%
    dplyr::summarise(af = mean(rr)) %>% 
    .$af %>%
    `names<-`("AF")

  structure(
    list(
      model = logit.result.opt,
      data = df,
      coef = coefs,
      tau = tau,
      af = af
    ),
    class = c("afmodel")
  )
}

#' @export
format.afmodel <- function(x, ...) {
  coef = coef(x)
  se = sqrt(diag(vcov(x)))
  cbind(
    coef = coef,
    se = se,
    lb = coef - qnorm(0.975) * se,
    ub = coef + qnorm(0.975) * se,
    z = (z <- coef/se),
    p.val = 2*pnorm(abs(z), lower.tail = FALSE)
  )    
}

#' @export
print.afmodel <- function(x,...) {
  # Results of the model
  print(format(x))
  cat(rep("=",10),"\n")
  #cat("Tau: ", x$tau, "\n")
  cat("AF: ", x$af, "\n")
  invisible(x)
}

#' @export
summary.afmodel <- function(object,...) {
  summary(object$model)
}

#' @export
logLik.afmodel <- function(object, ...){
  logLik(object$model)
}

#' @importFrom stats AIC
#' @export
AIC.afmodel <- function(object, ...){
  AIC(object$model)
}

#' @export
BIC.afmodel <- function(object, ...){
  BIC(object$model)
}

#' @export
coef.afmodel <- function(object,...) {
  object$coef
}

#' @export
vcov.afmodel <- function(object,...){
  vcov(object$model)
}
