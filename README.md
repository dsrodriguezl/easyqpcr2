
<!-- README.md is generated from README.Rmd. Please edit that file -->

# easyqpcr2

<!-- badges: start -->
<!-- badges: end -->

Rebuild of the abandoned EasyqpcR package, from Bioconductor (V 2.12).  
As EasyqpcR, this package is intended for the analysis of low-throughput
real-time quantitative PCR data, based on the qBase algorithms published
by Hellemans et al. in 2007.

Several functions were modyfied to be able to work with the newest
version of R, and/or to add extra functionalities to them. New versions
of the original functions state that they are modified of an original
EasyqpcR function, and of which, in their descriptions.

## Installation

You can install the development version of easyqpcr2 from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dsrodriguezl/easyqpcr2"
  , dependencies = TRUE
  , build_vignettes = TRUE)
```
