# Attributable fraction using a latent class model

By John Aponte and Orvalho Augusto.

## Introduction

In malaria endemic areas, asymptomatic carriage of malaria parasites
occurs frequently and the detection of malaria parasites in blood films
from a febrile individual does not necessarily indicate clinical
malaria.

A case definition for symptomatic malaria that is used widely in endemic
areas requires the presence of fever or history of fever together with a
parasite density above a specific cutoff. If the parasite density is
equal or higher than the cutoff point, the fever is considered due to
malaria.

How to estimate what is the sensitivity and specificity of the cutoff
point in the classification of the fever, without knowing what is the
true value of the fevers due to malaria?

Using the attributable fraction, one can estimate the expected number of
true cases due to malaria and with the positive predictive value
associated with a given cutoff point, we can estimate the expected
number of true cases among the fever cases that have a parasite density
higher or equal than the selected cutoff point.

In order to estimate the attributable fraction and the positive
predictive values, we follow the method proposed by Vounatsou et al (1)
fitting a bayesian latent class model.

The latent class model have several advantages over the logistic
exponential model: do not make any assumption on the shape of the risk
of the fever, only that at higher categories, higher the risk of fever
so there is little risk of bias, and provide direct estimation of the
attributable fraction and the other quantities like specificity and
sensitivity, allowing to estimate directly confidence intervals, however
the attributable fraction can be underestimated if the number of
categories is too high.

Here it is presented how to do it with the `afdx` package for the
R-software.

## Example using synthetic data

The data used in this example (malaria_df1) is a simulated data set as
seen frequently in malaria cross-sectionals where two main outcomes are
measured, the presence of fever or history of fever (fever column) and
the measured parasite density in parasites per $\mu l$ (density column).

``` r
library(afdx)
```

    #> 
    #> Attaching package: 'dplyr'
    #> The following objects are masked from 'package:stats':
    #> 
    #>     filter, lag
    #> The following objects are masked from 'package:base':
    #> 
    #>     intersect, setdiff, setequal, union
    #> 
    #> Attaching package: 'magrittr'
    #> The following object is masked from 'package:tidyr':
    #> 
    #>     extract
    #> 
    #> Attaching package: 'kableExtra'
    #> The following object is masked from 'package:dplyr':
    #> 
    #>     group_rows

| fever | density |
|------:|--------:|
|     1 |  475896 |
|     1 |   12008 |
|     0 |    1392 |
|     0 |    1664 |
|     0 |       0 |
|     1 |       0 |

head(malaria_df1) first 6 observations

In this simulation, there are 2000 observations, from which 785 have
fever or history of fever and 744 have a density of malaria greater than
0. A total of 437 have both fever and a malaria density higher than 0.

| k (category lower limit) | m (no fever) | n (fever) |
|:-------------------------|-------------:|----------:|
| 0                        |          908 |       348 |
| 1                        |           12 |         6 |
| 100                      |            8 |         4 |
| 200                      |           19 |         7 |
| 400                      |           37 |         8 |
| 800                      |           47 |        11 |
| 1600                     |           43 |        28 |
| 3200                     |           41 |        33 |
| 6400                     |           42 |        51 |
| 12800                    |           23 |        59 |
| 25600                    |           24 |        60 |
| 51200                    |           10 |        54 |
| 102400                   |            0 |        50 |
| 204800                   |            1 |        66 |

Distribution of fevers by density categories

## The latent class model

    * TODO: Include more Details on the model

    * TODO: Include Detail on the calculation of sensitivity, specificity and predictive values

## Estimating the bayesian latent class model

The `adfx` provide functions that facilitate the fitting of the bayesian
latent class model using the `rjags` package, but the user is
responsible to setup the appropriate bayesian workflow and confirm the
convergence of the model. Here we present one way to do it but there are
many other possibilities. We use a burn-in of 1000 iterations and
monitor samples from 4 chains 10000 iterations.

The function
[`get_latent_model()`](https://johnaponte.github.io/afdx/reference/get_latent_model.md)
provides a string with a model that can be use by `rjags`,

``` r
model <- get_latent_model()
cat(model)
#> 
#> data {
#>   # Number of categories
#>   K = length(n)
#>   # Total events by group
#>   Sn <- sum(n[])
#>   Sm <- sum(m[]) 
#> }
#> 
#> model {
#>   for (i in 1:2){   
#>     z0[i] <- (i-1)*0.0001
#>     phi0[i] <- theta[i]*z0[i] 
#>   }
#>   theta[1] <- 1-St
#>   eltheta[1] ~ dgamma(1.0,1.0) 
#>   theta[2] <- eltheta[2]/(1+Sr) 
#>   eltheta[2] ~ dgamma(1.0E-3,1.0E-3) 
#>   for (i in 3:K){ 
#>     phi0[i] <- theta[i]*z0[i] 
#>    eltheta[i] ~ dgamma(1.0E-3,1.0E-3) 
#>     theta[i] <- eltheta[i]/(1+Sr) 
#>     z0[i] <- z0[i-1]/q[i]
#>   }
#>   
#>   Sr <- sum(eltheta[2:K]) 
#>   St <- sum(theta[2:K]) 
#>   Sp <-sum(p0[]) 
#>   Sphi0 <- sum(phi0[]) 
#>   for (i in 1:K){
#>     phi[i] <- phi0[i]/Sphi0 
#>     z[i] <- z0[i]/Sphi0 
#>     q[i] ~ dunif(0.001,0.999)
#>     p0[i] <- theta[i]*(1-lambda)+lambda*phi[i] 
#>     p[i] <- p0[i]/Sp 
#>     lami[i] <- lambda*phi[i]/p0[i]
#>     
#>     # Sens and specificity
#>     sens[i] <- sum(phi[i:K])
#>     P[i] <- (1-lami[i])*p0[i]
#>     spec[i] <- sum(P[1:i])/sum(P)
#>     Q[i] <- lami[i]*p0[i]
#>     ppv[i] <- sum(Q[i:K])/sum(p0[i:K])
#>     npv[i] = spec[i] * (1 - lambda) /(spec[i] * (1 - lambda) + ((1 - sens[i]) * lambda))
#>   }
#> 
#>   # m is afebrile (training sample)
#>   m[1:K] ~ dmulti(theta[1:K], Sm) 
#>   
#>   # n is febrile (mixtured sample)
#>   n[1:K] ~ dmulti(p[1:K], Sn) 
#>   
#>   # lambda is the fraction from the mixture that belong to the g2
#>   lambda~ dunif(0.00001,0.99999)
#> }
```

Data in the model must be provided as a list with two vectors:

- `n` with the number of subjects with symptoms in the category

- `m` with the number of subjects without symptoms

The model calculate as data the number of categories (K) and the total
number of subjects in each group (Sn, Sm)

``` r
library(rjags)
library(coda)

# compile the model
af_latent <-
  jags.model(
    textConnection(get_latent_model()),
    data = list(n = data$`n (fever)`,
                m = data$`m (no fever)`),
    n.chains = 4,
    n.adapt = 1000,
    inits = list(
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 1111),
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 2222),
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 3333),
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 4444)
    )
  )

# Simulate the posterior
latent_sim <-
  coda.samples(
    model = af_latent,
    variable.names = c('lambda','sens','spec','ppv','npv'),
    n.thinning = 5,
    n.iter =  10000 )

# Extract and Analyze the posterior
latent_sum <-  summary(latent_sim)
latent_eff <-  effectiveSize(latent_sim)
```

``` r
# reformat to present the results
summary_table <-
  data.frame(latent_sum[[1]]) %>%
  bind_cols(data.frame(latent_sum[[2]])) %>%
  mutate(varname = row.names(latent_sum[[1]])) %>%
  mutate(cutoff  = c(NA, rep(cutoffs,4))) %>%
  select(varname, cutoff,Mean, X2.5., X50., X97.5.,Naive.SE ) %>%
  mutate(eff_size = floor(latent_eff)) %>%
  filter(is.na(cutoff) | cutoff != 0) 


mean_table <- summary_table %>%
  rename(point = Mean) %>%
  rename(lci = `X2.5.`) %>%
  rename(uci = `X97.5.`) %>%
  mutate(varname = gsub("\\[.*\\]","",varname)) %>%
  filter(varname != "lambda") %>%
  select(cutoff,varname, lci, uci,point) %>%
  pivot_longer(-c("cutoff","varname"),names_to = "xxv", values_to = "value") %>%
  unite("varx",varname,xxv ) %>%
  pivot_wider(names_from = "varx", values_from = "value") %>%
  select(cutoff,
         sens_point, 
         sens_lci, 
         sens_uci, 
         spec_point,
         spec_lci, 
         spec_uci,
         ppv_point, 
         ppv_lci, 
         ppv_uci,
         npv_point,
         npv_lci,
         npv_uci) %>%
  rename(sensitivity = sens_point) %>%
  rename(specificity = spec_point) %>%
  rename(ppv = ppv_point) %>%
  rename(npv = npv_point) %>%
  mutate_if(is.numeric, round,3)

# Lambda corresponds to the attributable fraction
afrow <- 
  summary_table %>%
  filter(varname == "lambda") %>% 
  mutate_if(is.numeric, round,3)
```

| cutoff | sensitivity | sens_lci | sens_uci | specificity | spec_lci | spec_uci |   ppv | ppv_lci | ppv_uci |   npv | npv_lci | npv_uci |
|-------:|------------:|---------:|---------:|------------:|---------:|---------:|------:|--------:|--------:|------:|--------:|--------:|
|      1 |       1.000 |    1.000 |    1.000 |       0.762 |    0.739 |    0.784 | 0.743 |   0.698 |   0.784 | 1.000 |   1.000 |   1.000 |
|    100 |       1.000 |    0.999 |    1.000 |       0.769 |    0.747 |    0.791 | 0.751 |   0.706 |   0.792 | 1.000 |   0.999 |   1.000 |
|    200 |       1.000 |    0.999 |    1.000 |       0.785 |    0.762 |    0.806 | 0.757 |   0.712 |   0.797 | 1.000 |   0.999 |   1.000 |
|    400 |       0.999 |    0.996 |    1.000 |       0.811 |    0.790 |    0.832 | 0.769 |   0.725 |   0.809 | 0.999 |   0.997 |   1.000 |
|    800 |       0.998 |    0.990 |    1.000 |       0.845 |    0.825 |    0.865 | 0.792 |   0.749 |   0.829 | 0.998 |   0.991 |   1.000 |
|   1600 |       0.993 |    0.976 |    1.000 |       0.883 |    0.864 |    0.900 | 0.821 |   0.782 |   0.857 | 0.994 |   0.980 |   1.000 |
|   3200 |       0.968 |    0.934 |    0.996 |       0.917 |    0.901 |    0.931 | 0.855 |   0.820 |   0.885 | 0.975 |   0.947 |   0.997 |
|   6400 |       0.916 |    0.867 |    0.964 |       0.951 |    0.938 |    0.963 | 0.887 |   0.857 |   0.913 | 0.940 |   0.901 |   0.976 |
|  12800 |       0.808 |    0.747 |    0.869 |       0.972 |    0.963 |    0.980 | 0.922 |   0.897 |   0.944 | 0.875 |   0.828 |   0.920 |
|  25600 |       0.666 |    0.604 |    0.729 |       0.990 |    0.984 |    0.995 | 0.945 |   0.924 |   0.962 | 0.804 |   0.755 |   0.852 |
|  51200 |       0.502 |    0.440 |    0.566 |       0.998 |    0.996 |    1.000 | 0.973 |   0.957 |   0.986 | 0.735 |   0.686 |   0.784 |
| 102400 |       0.348 |    0.293 |    0.406 |       0.999 |    0.998 |    1.000 | 0.994 |   0.982 |   0.999 | 0.680 |   0.631 |   0.729 |
| 204800 |       0.201 |    0.158 |    0.248 |       1.000 |    1.000 |    1.000 | 0.995 |   0.986 |   1.000 | 0.634 |   0.588 |   0.681 |

Summary of diagnostic characteristics at selected cutoff points

Attributable fraction: 0.419 95%CI(0.373, 0.463 )

## Bibliography

1.  Vounatsou P, Smith T, Smith AFM. Bayesian analysis of two-component
    mixture distributions applied to estimating malaria attributable
    fractions. Journal of the Royal Statistical Society: Series C
    (Applied Statistics). 1998;47(4):575–87. DOI:
    10.1111/1467-9876.00129

2.  Müller I, Genton B, Rare L, Kiniboro B, Kastens W, Zimmerman P, et
    al. Three different Plasmodium species show similar patterns of
    clinical tolerance of malaria infection. Malar J. 2009;8(1):158.
    DOI: 10.1186/1475-2875-8-158

3.  Plucinski MM, Rogier E, Dimbu PR, Fortes F, Halsey ES, Aidoo M, et
    al. Performance of Antigen Concentration Thresholds for Attributing
    Fever to Malaria among Outpatients in Angola. J Clin Microbiol.
    2019;57(3):e01901-18. DOI: 10.1128/JCM.01901-18
