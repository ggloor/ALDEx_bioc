
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/michellepistner/ALDEx_bioc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/michellepistner/ALDEx_bioc/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Introduction

Welcome to a scale simulation within ALDEx2!

The `ALDEx2` package is a Bioconductor package for differential
abundance analysis across two or more conditions. It is useful for
analyzing data from standard RNA-seq or meta-RNA-seq assays as well as
selected and unselected values from in-vitro sequence selections. Unlike
other packages, `ALDEx2` uses a Dirichlet-multinomial model to infer
abundance from counts, optimized for three or more experimental
replicates. The method infers biological and sampling variation to
calculate the expected false discovery rate, given the variation, based
on a Wilcox rank test or Welch t-test (via aldex.ttest), or a glm and
Kruskal-Wallis test (via aldex.glm). The `ALDEx2` package reports
p-values and Benjamini-Hochberg corrected p-values. Effect sizes \> 1
are generally preferred metrics. This repository also allows for scale
simulation to be incorporated within ALDEx2.

## Quick start

You can install the developmental branch of `ALDEx2` plus scale
simulation from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("ggloor/ALDEx_bioc")
```

Getting started with `ALDEx2` is easy. All you need is a matrix (with
rows as variables and columns as samples) and a character vector of
group labels. Finally, use the `denom` argument to choose a set of
variables to use as the reference for the analysis. You can provide a
user-defined reference set (e.g., known house-keeping genes), or choose
a method that finds references from the data (`denom = "iqlr"` usually
performs well!).

``` r
library(ALDEx2)
#> Loading required package: zCompositions
#> Loading required package: MASS
#> Loading required package: NADA
#> Loading required package: survival
#> 
#> Attaching package: 'NADA'
#> The following object is masked from 'package:stats':
#> 
#>     cor
#> Loading required package: truncnorm
#> Loading required package: lattice
#> Loading required package: latticeExtra
data(selex)
group <- c(rep("A", 7), rep("B", 7))
res <- aldex(selex, group, denom = "iqlr")
#> aldex.clr: generating Monte-Carlo instances and clr values
#> operating in serial mode
#> computing iqlr centering
#> aldex.ttest: doing t-test
#> aldex.effect: calculating effect sizes
```

See the scaleSim or ALDEx2 vignettes for more details.
