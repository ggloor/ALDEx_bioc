
<!-- README.md is generated from README.Rmd. Please edit that file -->
Introduction
------------

Welcome to the **official active development branch** of the `ALDEx2` package!

The `ALDEx2` package is a Bioconductor package for differential abundance analysis across two or more conditions. It is useful for analyzing data from standard RNA-seq or meta-RNA-seq assays as well as selected and unselected values from in-vitro sequence selections. Unlike other packages, `ALDEx2` uses a Dirichlet-multinomial model to infer abundance from counts, optimized for three or more experimental replicates. The method infers biological and sampling variation to calculate the expected false discovery rate, given the variation, based on a Wilcox rank test or Welch t-test (via aldex.ttest), or a glm and Kruskal-Wallis test (via aldex.glm). The `ALDEx2` package reports p-values and Benjamini-Hochberg corrected p-values.

Quick start
-----------

You can install the developmental branch of `ALDEx2` from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("ggloor/ALDEx_bioc")
```

Getting started with `ALDEx2` is easy. All you need is a matrix (with rows as variables and columns as samples) and a character vector of group labels. Finally, use the `denom` argument to choose a set of variables to use as the reference for the analysis. You can provide a user-defined reference set (e.g., known house-keeping genes), or choose a method that finds references from the data (`denom = "iqlr"` usually performs well!).

``` r
library(ALDEx2)
data(selex)
group <- c(rep("A", 7), rep("B", 7))
res <- aldex(selex, group, denom = "iqlr")
#> [1] "aldex.clr: generating Monte-Carlo instances and clr values"
#> [1] "operating in serial mode"
#> [1] "computing iqlr centering"
#> [1] "aldex.ttest: doing t-test"
#> [1] "running tests for each MC instance:"
#> |------------(25%)----------(50%)----------(75%)----------|
#> [1] "aldex.effect: calculating effect sizes"
#> [1] "operating in serial mode"
```

See the vignette for more details.
