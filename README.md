# afdx: Diagnosis performance using attributable fraction

## Introduction

This R-package help on the estimation of diagnosis performance 
(Sensitivity, Specificity, Positive predictive value, Negative predicted value) 
of a diagnostic test where the golden standard can't be measured but can be
estimated using the attributable fraction

Two methods are presented with examples for Malaria diagnosis, using a maximum
likelihood estimated logistic exponential model  and using a bayesian latent 
class model.

To install the package from github use:

`devtools::install_github("johnaponte/afdx", build_manual = T, build_vignettes = T)`


