# Template for the bayesian latent class model

This function returns a template that can be use as model in an `rjags`
model it requires two vectors with the number of subjects in the
symptoms, like fever in the case of malaria (n) and the number of
non-symptomatic (m) in each of the categories of results of the
diagnostic test. The first category is reserved for the negatives by the
diagnostic test (in the malaria case those with asexual density 0) and
the rest categories each one with higher values than the previous
category.

## Usage

``` r
get_latent_model()
```

## Value

a string value

## Details

See: Smith T, Vounatsou P. Logistic regression and latent class models
for estimating positivities in diagnostic assays with poor resolution.
Communications in Statistics - Theory and Methods. 1997
Jan;26(7):1677–700.

Vounatsou P, Smith T, Smith AFM. Bayesian analysis of two-component
mixture distributions applied to estimating malaria attributable
fractions. Journal of the Royal Statistical Society: Series C (Applied
Statistics). 1998;47(4):575–87.

Müller I, Genton B, Rare L, Kiniboro B, Kastens W, Zimmerman P, et al.
Three different Plasmodium species show similar patterns of clinical
tolerance of malaria infection. Malar J. 2009;8(1):158.

Plucinski MM, Rogier E, Dimbu PR, Fortes F, Halsey ES, Aidoo M, et al.
Performance of Antigen Concentration Thresholds for Attributing Fever to
Malaria among Outpatients in Angola. J Clin Microbiol. 2019;57(3).

## Examples

``` r
{
 get_latent_model()
}
#> [1] "\ndata {\n  # Number of categories\n  K = length(n)\n  # Total events by group\n  Sn <- sum(n[])\n  Sm <- sum(m[]) \n}\n\nmodel {\n  for (i in 1:2){   \n    z0[i] <- (i-1)*0.0001\n    phi0[i] <- theta[i]*z0[i] \n  }\n  theta[1] <- 1-St\n  eltheta[1] ~ dgamma(1.0,1.0) \n  theta[2] <- eltheta[2]/(1+Sr) \n  eltheta[2] ~ dgamma(1.0E-3,1.0E-3) \n  for (i in 3:K){ \n    phi0[i] <- theta[i]*z0[i] \n   eltheta[i] ~ dgamma(1.0E-3,1.0E-3) \n    theta[i] <- eltheta[i]/(1+Sr) \n    z0[i] <- z0[i-1]/q[i]\n  }\n  \n  Sr <- sum(eltheta[2:K]) \n  St <- sum(theta[2:K]) \n  Sp <-sum(p0[]) \n  Sphi0 <- sum(phi0[]) \n  for (i in 1:K){\n    phi[i] <- phi0[i]/Sphi0 \n    z[i] <- z0[i]/Sphi0 \n    q[i] ~ dunif(0.001,0.999)\n    p0[i] <- theta[i]*(1-lambda)+lambda*phi[i] \n    p[i] <- p0[i]/Sp \n    lami[i] <- lambda*phi[i]/p0[i]\n    \n    # Sens and specificity\n    sens[i] <- sum(phi[i:K])\n    P[i] <- (1-lami[i])*p0[i]\n    spec[i] <- sum(P[1:i])/sum(P)\n    Q[i] <- lami[i]*p0[i]\n    ppv[i] <- sum(Q[i:K])/sum(p0[i:K])\n    npv[i] = spec[i] * (1 - lambda) /(spec[i] * (1 - lambda) + ((1 - sens[i]) * lambda))\n  }\n\n  # m is afebrile (training sample)\n  m[1:K] ~ dmulti(theta[1:K], Sm) \n  \n  # n is febrile (mixtured sample)\n  n[1:K] ~ dmulti(p[1:K], Sn) \n  \n  # lambda is the fraction from the mixture that belong to the g2\n  lambda~ dunif(0.00001,0.99999)\n}\n"
```
