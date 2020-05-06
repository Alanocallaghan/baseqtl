
<!-- README.md is generated from README.Rmd. Please edit that file -->

# baseqtl

BaseQTL is a suite of models to discover eQTLs combining between
individual variation (modelled by a negative binomial distribution) with
allele specific expression (ASE modelled by a beta-binomial
distribution). We descibe models for 4 scenarios:

1.  eQTL analysis when DNA-seq and RNA-seq data are available (with
    genotypes).
2.  eQTL analysis when only RNA-seq is available (unknown genotypes).
3.  Condition specific eQTL analysis in paired samples (two treatments
    on same samples) with genotypes.
4.  Condition specific eQTL analysis in different samples (diseases vs
    control) and unkown gentopyes.

## System requirements

    R versions >= 3.4.0.
    GNU make
    bcftools

## Instalations

bcftools can be download it from
[bcftools](http://www.htslib.org/download/)

Move the downloaded file to your bin directory to build

``` bash
# uncompress
tar xvfj bcftools-x.tar.bz2
cd bcftools-x
./configure --prefix=/where/to/intall
make
make install
# The executable programs will be installed to a bin subdirectory
# Add bcftools path to $PATH, otherwise R won't find bcftools
export PATH=/where/to/install/bin:$PATH
```

You can install the released version of baseqtl from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("baseqtl")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
