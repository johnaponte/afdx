# Template for the bayesian latent class model to be use in rJags
# Code is based on the work from Vountasou and Smith, adapted
# by O. Augusto and optimized by J. Aponte.
# 20200502 by JJAV


#' Template for the bayesian latent class model
#' 
#' This function returns a template that can be use as model in an rJags model
#' it requires two vectors with the number of subjects in the symptoms, 
#' like fever in the case of malaria (n) and the number of non-symptomatic (m)
#' in each of the categories of results of the diagnostic test. The first 
#' category is reserved for the negatives by the diagnostic test (in the malaria
#' case those with asexual density 0) and the rest categories each one with 
#' higher values than the previous category.
#' 
#' See:
#' Smith T, Vounatsou P. Logistic regression and latent class models for 
#' estimating positivities in diagnostic assays with poor resolution. 
#' Communications in Statistics - Theory and Methods. 1997 Jan;26(7):1677–700. 
#'
#' Vounatsou P, Smith T, Smith AFM. Bayesian analysis of two-component mixture 
#' distributions applied to estimating malaria attributable fractions. Journal 
#' of the Royal Statistical Society: Series C (Applied Statistics). 
#' 1998;47(4):575–87. 
#'
#' Müller I, Genton B, Rare L, Kiniboro B, Kastens W, Zimmerman P, et al. 
#' Three different Plasmodium species show similar patterns of clinical 
#' tolerance of malaria infection. Malar J. 2009;8(1):158. 
#'  
#' Plucinski MM, Rogier E, Dimbu PR, Fortes F, Halsey ES, Aidoo M, et al.
#' Performance of Antigen Concentration Thresholds for Attributing Fever to 
#' Malaria among Outpatients in Angola. J Clin Microbiol. 2019;57(3). 
#' 
#' @return a string value
#' @export
get_latent_model <- function(){
  "
data {
  # Number of categories
  K = length(n)
  # Total events by group
  Sn <- sum(n[])
  Sm <- sum(m[]) 
}

model {
  for (i in 1:2){   
    z0[i] <- (i-1)*0.0001
    phi0[i] <- theta[i]*z0[i] 
  }
  theta[1] <- 1-St
  eltheta[1] ~ dgamma(1.0,1.0) 
  theta[2] <- eltheta[2]/(1+Sr) 
  eltheta[2] ~ dgamma(1.0E-3,1.0E-3) 
  for (i in 3:K){ 
    phi0[i] <- theta[i]*z0[i] 
   eltheta[i] ~ dgamma(1.0E-3,1.0E-3) 
    theta[i] <- eltheta[i]/(1+Sr) 
    z0[i] <- z0[i-1]/q[i]
  }
  
  Sr <- sum(eltheta[2:K]) 
  St <- sum(theta[2:K]) 
  Sp <-sum(p0[]) 
  Sphi0 <- sum(phi0[]) 
  for (i in 1:K){
    phi[i] <- phi0[i]/Sphi0 
    z[i] <- z0[i]/Sphi0 
    q[i] ~ dunif(0.001,0.999)
    p0[i] <- theta[i]*(1-lambda)+lambda*phi[i] 
    p[i] <- p0[i]/Sp 
    lami[i] <- lambda*phi[i]/p0[i]
    
    # Sens and specificity
    sens[i] <- sum(phi[i:K])
    P[i] <- (1-lami[i])*p0[i]
    spec[i] <- sum(P[1:i])/sum(P)
    Q[i] <- lami[i]*p0[i]
    ppv[i] <- sum(Q[i:K])/sum(p0[i:K])
    npv[i] = spec[i] * (1 - lambda) /(spec[i] * (1 - lambda) + ((1 - sens[i]) * lambda))
  }

  # m is afebrile (training sample)
  m[1:K] ~ dmulti(theta[1:K], Sm) 
  
  # n is febrile (mixtured sample)
  n[1:K] ~ dmulti(p[1:K], Sn) 
  
  # lambda is the fraction from the mixture that belong to the g2
  lambda~ dunif(0.00001,0.99999)
}
"
}
