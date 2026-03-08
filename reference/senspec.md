# S3 methods to estimate diagnosis performance of an afmodel

Estimate sensitivity, specificity, positive predicted value and negative
predicted value negative predictive value from an afmodel. The estimated
"true" negative and "true" positive are estimated using the estimated
overall attributable fraction and the predictive positive value
associated with each cut-off point as described by Smith, T.,
Schellenberg, J.A., Hayes, R., 1994. Attributable fraction estimates and
case definitions for malaria in endemic areas. Stat Med 13, 2345–2358.

## Usage

``` r
senspec(object, ...)

# Default S3 method
senspec(object, ...)

# S3 method for class 'afmodel'
senspec(object, cutoff, ...)
```

## Arguments

- object:

  with the data to calculate the sensitivity and specificity

- ...:

  other parameters for the implementing functions

- cutoff:

  vector of cut-off points to make the estimations

## Value

a matrix with the columns sensitivity and specificity, ppv (positive
predicted value) and npv (negative predicted value)

No return value. Raise an error.

a matrix with the columns sensitivity and specificity, ppv (positive
predicted value) and npv (negative predicted value)

## See also

[`logitexp`](https://johnaponte.github.io/afdx/reference/logitexp.md)

## Examples

``` r
{
# Get the sample data
head(malaria_df1)
fit <- logitexp(malaria_df1$fever, malaria_df1$density)
fit
senspec(fit,  c(1,100,500,1000,2000,4000,8000,16000, 32000,54000,100000))
}
#>       cutoff sensitivity specificity       ppv       npv
#>  [1,]      1   1.0000000   0.7727793 0.7658522 1.0000000
#>  [2,]    100   0.9986640   0.7851102 0.7754762 0.9987369
#>  [3,]    500   0.9937948   0.8059183 0.7919063 0.9943103
#>  [4,]   1000   0.9861357   0.8224325 0.8049691 0.9876265
#>  [5,]   2000   0.9720090   0.8430224 0.8214885 0.9759178
#>  [6,]   4000   0.9250541   0.8858478 0.8576031 0.9408427
#>  [7,]   8000   0.8707224   0.9165290 0.8857480 0.9051178
#>  [8,]  16000   0.7372786   0.9594746 0.9311339 0.8309098
#>  [9,]  32000   0.5936829   0.9837722 0.9645255 0.7651379
#> [10,]  54000   0.4863330   0.9928155 0.9805100 0.7222734
#> [11,] 100000   0.3470470   0.9981097 0.9927245 0.6728613
```
