# Creates synthetic data to test the afdx package
# inspired on real crossectional data
# 20210201 By John Aponte

library(tidyverse)


# Dataset 1
prev = 0.35 # Proportion of densities > 0 in the sample
meanlog10den = 4 # mean of the log10 densities
sdlog10den = 1 #  sd of the log10 denstities
ss = 2000 # samplesize 
# Parameters for the exponential logistic regression
alphacoef = -1
betacoef = 0.01
tau = 0.5

set.seed(12345)
malaria_df1 <-
  tibble(
    positive =
      rbinom(ss,1,prev)
  )  %>%
  mutate(density = round(10^(rnorm(ss,meanlog10den,sdlog10den)))*positive) %>%
  mutate(odds_fever = exp(alphacoef + (betacoef*(density^tau)))) %>%
  mutate(pfever = odds_fever/(1+odds_fever)) %>%
  mutate(fever = rbinom(ss,1,pfever)) %>%
  select(fever, density)


# Dataset 2
prev = 0.75
alphacoef = -2
betacoef = 0.001
tau = 0.8
meanlog10den = 3
sdlog10den = 0.5

set.seed(657891)
malaria_df2 <-
  tibble(
    positive =
      rbinom(ss,1,prev)
  )  %>%
  mutate(density = round(10^(rnorm(ss,meanlog10den,sdlog10den)))*positive) %>%
  mutate(odds_fever = exp(alphacoef + (betacoef*(density^tau)))) %>%
  mutate(pfever = odds_fever/(1+odds_fever)) %>%
  mutate(fever = rbinom(ss,1,pfever)) %>%
  select(fever, density)


usethis::use_data(malaria_df1, malaria_df2, overwrite = T)
