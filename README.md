
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bfboin

<!-- badges: start -->
<!-- badges: end -->

The goal of bfboin is to provide an independent implementation of the
[BF-BOIN design](https://pubmed.ncbi.nlm.nih.gov/38048044/) aiming to
reproduce the functionality of the [shiny app by MD
Anderson](https://biostatistics.mdanderson.org/shinyapps/BF-BOIN/).

## Installation

You can install the development version of bfboin like so:

``` r
# install.packages("devtools")
devtools::install_gitlab("boin/bfboin")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(bfboin)
## basic example code
get.oc.bf(ntrial = 100,
          seed = 9,
          target = 0.25,
          p.true = c(0.1, 0.5),
          ncohort = 10,
          cohortsize = 3,
          n.earlystop = 9,
          startdose = 1,
          titration = FALSE,
          cutoff.eli = 0.95,
          extrasafe = TRUE,
          offset = 0.1,
          boundMTD=FALSE,
          n.cap = 12,
          end.backfill = TRUE,
          n.per.month = 1,
          dlt.window = 1,
          p.response.true = c(0.001, 0.001))
#> $selpercent
#> [1] 81 15
#> 
#> $npatients
#> [1] 10.58  7.08
#> 
#> $percentpatients
#> [1] 59.9094 40.0906
#> 
#> $ntox
#> [1] 1.18 3.57
#> 
#> $totaltox
#> [1] 4.75
#> 
#> $totaln
#> [1] 17.66
#> 
#> $percentstop
#> [1] 4
#> 
#> $duration
#> [1] 22.13623
```
