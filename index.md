# Adaptive gPCA

## Description

Adaptive gPCA, described in
[Fukuyama 2017](https://arxiv.org/abs/1702.00501), is a flexible
method for incorporating side information into principal components
analysis. It was developed for using information about the
phylogenetic structure of bacteria in microbiome data analysis, but it
is applicable to more general kinds of structure. 

## Installation

Adaptive gPCA is implemented in the R package adaptiveGPCA, which can
be installed either from CRAN or github. To install from CRAN, use
```r
install.packages("adaptiveGPCA")
```

To install from GitHub, first install devtools, and then use
```r
devtools::install_github("jfukuyama/adaptiveGPCA")
```
