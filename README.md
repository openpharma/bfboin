
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bfboin

<!-- badges: start -->
<!-- badges: end -->

The goal of bfboin is to provide an independent implementation of the
[BF-BOIN design](https://pubmed.ncbi.nlm.nih.gov/38048044/) of Zhou et
al. (2024), aiming to reproduce the functionality of the [shiny app by
MD Anderson](https://biostatistics.mdanderson.org/shinyapps/BF-BOIN/).

Zhao, Y., Yuan, Y., Korn, E.L. and Freidlin, B., 2024. Backfilling
patients in phase I dose-escalation trials using Bayesian optimal
interval design (BOIN). Clinical Cancer Research, 30(4), pp.673-679.

## Installation

You can install the development version of bfboin like so:

``` r
# install.packages("devtools")
devtools::install_github("openpharma/bfboin")
```

## Example

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
          p.response.true = c(0.001, 0.001),
          accrual = "uniform")
#> $selpercent
#> [1] 75 18
#> 
#> $npatients
#> [1] 10.54  7.41
#> 
#> $percentpatients
#> [1] 58.71866 41.28134
#> 
#> $ntox
#> [1] 1.01 3.71
#> 
#> $totaltox
#> [1] 4.72
#> 
#> $totaln
#> [1] 17.95
#> 
#> $percentstop
#> [1] 7
#> 
#> $duration
#> [1] 20.52315
```
