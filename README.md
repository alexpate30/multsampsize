
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multsampsize

<!-- badges: start -->
<!-- badges: end -->

The goal of multsampsize is to …

## Installation

You can install the development version of multsampsize from
[GitHub](https://github.com/) with:

``` r
# pak::pak("alexpate30/multsampsize")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(multsampsize)
## basic example code
```

This can be used in the same way as `pmsampsize::pmsampsize` by
specifying `type = "m"`.

The number of outcome categories `K` must be specified.

The number of expected events in each outcome category must be
specified. Note, this is used to estimate the prevalence of each outcome
category relative to eachother. Therefore, the absolute values are not
important, just their relative values. These values can therefore be
estimated from a previously published study, even if this study is of a
different sample size to the study you are planning to implement.

One, and only one of nagrsquared, csrsquared and cstatistic must be
specified. These must be vectors of length equal to the number of pairs
of outcome categories. The order must be sequential, first in category
k, then category r. For example, when `K = 3`, the elements must be as
follows:

- Element 1: Pair(1,2); k = 1, r = 2.
- Element 2: Pair(1,3); k = 1, r = 3.
- Element 3: Pair(2,3); k = 2, r = 3.

When `K = 4`, the elements must be as follows:

- Element 1: Pair(1,2); k = 1, r = 2.
- Element 2: Pair(1,3); k = 1, r = 3.
- Element 3: Pair(1,4); k = 1, r = 4.
- Element 4: Pair(2,3); k = 2, r = 3.
- Element 5: Pair(2,4); k = 2, r = 4.
- Element 6: Pair(3,4); k = 3, r = 4.

One, and only one of `mult_nagrsquared_overall` and
`mult_rsquared_overall` must be specified. The rsquared of a multinomial
logistic regression model, `mult_rsquared_overall`, is rarely reported
and may be difficult to pre-specify. It will therefore be common to
specify `mult_nagrsquared_overall`

This first examples mimics the worked example from manuscript Pate et
al.

``` r
my_sampsize_calc <- pmsampsize_mult_general(type = "m",
                                            parameters = 17,
                                            shrinkage = 0.9,
                                            cstatistic = c(0.85, 0.92, 0.99, 0.95, 0.75, 0.95, 0.87, 0.87, 0.71, 0.82),
                                            ### New parameters
                                            K = 5,
                                            mult_n_events = c(2557, 186, 176, 467, 120),
                                            mult_nagrsquared_overall = 0.15)
```

The output contains the expected pairwise performance metrics
(pertaining to criterion i), as well as the sample size for each
criterion.

``` r
my_sampsize_calc
#> $results_criteria1
#>          CS_Rsq Max_Rsq Nag_Rsq targ_shrinkage
#> Pair 1,2  0.116   0.391   0.297          0.988
#> Pair 1,3  0.179   0.380   0.471          0.993
#> Pair 1,4  0.497   0.577   0.861          0.998
#> Pair 1,5  0.168   0.306   0.549          0.992
#> Pair 2,3  0.185   0.750   0.247          0.945
#> Pair 2,4  0.499   0.697   0.716          0.991
#> Pair 2,5  0.372   0.738   0.504          0.972
#> Pair 3,4  0.329   0.691   0.477          0.985
#> Pair 3,5  0.128   0.741   0.173          0.900
#> Pair 4,5  0.210   0.637   0.330          0.971
#> 
#> $results_all
#> Criteria 1 Criteria 2 Criteria 3      Final 
#>      13136       1477        524      13136 
#> 
#> $sample_size
#> [1] 13136
#> 
#> $parameters
#> [1] 17
#> 
#> $events_expected
#> [1] 9580.3628  696.8899  659.4227 1749.7182  449.6064
#> 
#> $type
#> [1] "multinomial"
```

Another random example:

``` r
my_sampsize_calc <- pmsampsize_mult_general(type = "m",
                                            nagrsquared = c(0.15,0.15,0.15),
                                            parameters = 17,
                                            shrinkage = 0.9,
                                            ### New parameters
                                            K = 3,
                                            mult_n_events = c(50,100,150),
                                            mult_nagrsquared_overall = 0.15)
my_sampsize_calc
#> $results_criteria1
#>          CS_Rsq Max_Rsq Nag_Rsq targ_shrinkage
#> Pair 1,2  0.108   0.720    0.15          0.900
#> Pair 1,3  0.101   0.675    0.15          0.920
#> Pair 2,3  0.111   0.740    0.15          0.942
#> 
#> $results_all
#> Criteria 1 Criteria 2 Criteria 3      Final 
#>       2660        714        664       2660 
#> 
#> $sample_size
#> [1] 2660
#> 
#> $parameters
#> [1] 17
#> 
#> $events_expected
#> [1]  443.3333  886.6667 1330.0000
#> 
#> $type
#> [1] "multinomial"
```
