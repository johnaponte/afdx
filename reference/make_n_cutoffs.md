# Make a defined number of categories having similar number of positives in each category

Generate the categories in a way that each category have at least the
`mintot` number of observation. It generate all possible categories were
there is change and then collapse to have minimum number of observations
in each category

## Usage

``` r
make_n_cutoffs(v.fever, v.density, mintot, add1 = TRUE)
```

## Arguments

- v.fever:

  numeric vector of 0/1 indicating fever or equivalent

- v.density:

  numeric vector of values \>= 0 indicating the density

- mintot:

  minimum number of observations per category

- add1:

  a logical value to indicate the category started with 1 is included

## Value

a vector with the cutoff points

## Examples

``` r
{
make_n_cutoffs(malaria_df1$fever, malaria_df1$density, mintot=50)
}
#>  [1]      0      1    297    728   1360   2135   3446   5721   8477  12314
#> [11]  17572  27360  43152  68500  95757 269349
```
