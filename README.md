
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BaseQTL

**B**ayesian **a**llele **s**pecific **e**xpression of **QTL**, is a
suite of models to discover molecular QTLs combining between individual
variation (modelled by a negative binomial distribution, NB) with allele
specific expression (ASE modelled by a beta-binomial distribution). We
descibe models for 4 scenarios:

1.  eQTL analysis when DNA-seq and RNA-seq data are available (with
    genotypes).
2.  eQTL analysis when only RNA-seq is available (unknown genotypes).
3.  Condition specific eQTL analysis in paired samples (two treatments
    on same samples) with genotypes.
4.  Condition specific eQTL analysis in different samples (diseases vs
    control) and unkown gentopyes.

## System requirements

    R versions >= 3.4.0.
    GNU make >= 3.82
    bcftools >= 1.3.1

## Instalations

### bcftools can be downloaded from [htslib](http://www.htslib.org/download/)

Move the downloaded file to your bin directory to build

``` bash
# uncompress
tar xvfj bcftools-x.tar.bz2
# build
cd bcftools-x
./configure --prefix=/where/to/intall
make
make install
# The executable programs will be installed to a bin subdirectory
# Add bcftools path to $PATH in your .bashrc or equivalent file, otherwise R won't find bcftools
export PATH=/where/to/install/bin:$PATH
```

### Install baseqtl from [GitLab](https://gitlab.com):

Installation has been tested on R 3.5.1. Installation time is estimated
as 2 minutes.

``` r
## Within R:
install.packages("devtools") # if you don't already have the package
library(devtools)
devtools::install_git(url = "https://gitlab.com/evigorito/addstrings.git") 
devtools::install_git(url = "https://gitlab.com/evigorito/baseqtl.git") 
```

In its current form baseqtl depends on GUESSFM, to install GUESSFM do:

``` r
library(devtools)
devtools::install_github("chr1swallace/GUESSFM", ref="groups")
## you may need to install additional R packages that GUESSFM requires
```

## Running baseqtl

Preparation of the required files to run these functions can be made as
described in
[baseqtl\_pipeline](https://gitlab.com/evigorito/baseqtl_pipeline).

### eQTL analysis when DNA-seq and RNA-seq data are available

The function to call is **baseqtl.gt**. It will test for eQTL effects in
one specified gene. To show an example we provide data files for gene
ENSG00000159958 on chromosome 22. Genome coordinates are in built 37.

**Arguments**

*gene*: ensembl gene id

*chr*: number of chromosome

*snps*: if numeric, it will test for eQTL effects within a cis-window
from start/end of the gene expressed in base pairs, defaults to
\(5*10^5\). For testing specific SNPs input a character vector with
pos:reference:alternative allele (example snps=c(“13444352:A:G”,
“13444567:T:C”)).

*counts.f*: full name of a txt file with total gene counts, first column
is gene\_id, followed by samples, details in
[snakefile](https://gitlab.com/evigorito/baseqtl_pipeline/-/blob/master/input/Snakefile)
output from rule total\_gene\_counts.

*covariates*: full name to rds file with a matrix of covariates, details
in
[snakefile](https://gitlab.com/evigorito/baseqtl_pipeline/-/blob/master/input/Snakefile)
output from rule total\_gene\_counts. For running the analysis without
covariates set covariates=1. You can add extra columns to the matrix for
additional covariates

*additional\_cov*: full name to file with covariates that are gene
independent, especially useful when using argument *covariates* with
library size adjusted by GC-content. Format is first column sample names
and additional columns with covariate information. Defaults to NULL

*e.snps*: full name of txt file with a list of exonic SNPS across genes,
details in
[snakefile](https://gitlab.com/evigorito/baseqtl_pipeline/-/blob/master/input/Snakefile)
output from rule fSNP\_gene (fsnps output).

*u.esnps*: optional argument, when strand information is not available
for RNA-seq some exonic SNPs could be shared between genes making
difficult to assess allele specific expression. In this case it is
recommended to provide a list of exonic SNPs uniquely mapping genes. In
this mode, e.snps will be used to imporve phasing accuracy and u.esnps
to compute allele specific expression, further deatails in
[snakefile](https://gitlab.com/evigorito/baseqtl_pipeline/-/blob/master/input/Snakefile)
output from rule fSNP\_gene (ufsnps output).

*gene.coord*: full name to file with gene id, gene sart and gene end,
details in
[snakefile](https://gitlab.com/evigorito/baseqtl_pipeline/-/blob/master/input/Snakefile)
output from rule exon\_by\_gene.

*vcf*: full name to vcf file with genotypes and allele specific
expression, details in
[snakefile](https://gitlab.com/evigorito/baseqtl_pipeline/-/blob/master/input/Snakefile)
output from rule merge\_vcf (source=“DNA”).

*le.file*: full name to legend file with external reference panel SNP
description. We use the [1000 Genomes
Phase3](https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html).

*h.file*: full name to hap file with haplotypes for external reference
panel.

*population*: one of AFR AMR EAS EUR SAS ALL, defaults to EUR. Used to
calculate EAF of eQTL SNP to include in output.

*nhets*: minimun number of hetrozygous individuals for the eQTL SNP to
run model, defaults to 5.

*min.ase*: minimun number of ASE counts for an individual to be included
in the model, deaults to 5.

*min.ase.het*: minumun number of heterozygous individuals for the eQTL
SNP with sufficient ASE counts in order to run the ASE component of the
model, defaults to 5.

*min.ase.n*: minimun number of individual (any genotype) with sufficient
ASE counts to run ASE model, defaults to 5.

*tag.threshold*:
![r^2](https://render.githubusercontent.com/render/math?math=r%5E2)
threshold for grouping SNPs to reduce the number of running test,
recommended when running with cis-windows, defaults to 0.9. To disable,
when running selected SNPS set to “no”.

*out*: path to output directory to write files, deafults to current
directory.

*prefix*: prefix to add to output files, defaults to gene id.

*model*: whether to run NB-ASE only (full model negative binomial and
allele specific counts),NB (negat ive binomial only) or both (NB-ASE
when sufficient information and NB for those associations with no ASE
information), defaults to “both”.

*stan.model*: optional parameter, defaults is to run built-in NB-ASE
model. When AI\_estimate is provided (below) the model corrects for
reference panel bias, otherwise id doesn’t. Alternatively, users can
provide their alternative compiled stanmodel object, whcih should use
the same input data.

*stan.negonly*: compiled stanmodel object with neg only side, deafults
to built-in model. model.

*prob*: number p∈(0,1) indicating the desired posterior probability mass
to report, defaults to 0.99 -0.95 quantiles.

*prior*: prior for eQTL effect. Defaults to a mixture of 2 Gaussians
learnt from external eQTL data.

*ex.fsnps*: for sensitivity analysis it is possible to exclude selected
exonic (feature) SNPs by proving a character vector with
psoition:reference:alternative allele, defaults to NULL (no exlcusion).

*AI\_estimate*: full name to txt file with pre-computed allelic
imbalance estimates for reference panel bias for exonic SNPs, details in
[snakefile](https://gitlab.com/evigorito/baseqtl_pipeline/-/blob/master/input/Snakefile)
output from rule get\_AI with source=“DNA”. Defaults to NULL for no
correction.

*pretotalReads*: cut-off for total initial reads to consider AI
estimates, defaults to 100, details in **input4baseqtl**. When reference
panel correction is not enabled this argument won’t be used.

``` r
## Retrive input files for running baseqtl.gt
## Most of the files contain information for the selected gene only to minimise their size

counts.f <- system.file("extdata/input", "counts.txt", package = "baseqtl", mustWork = TRUE)
covariates <- system.file("extdata/input", "lbsize_gc.rds", package = "baseqtl", mustWork = TRUE)
e.snps = system.file("extdata/input", "chr22.fSNPS.ENSG00000159958.txt", package = "baseqtl", mustWork = TRUE)
u.snps = system.file("extdata/input", "chr22.unique.fSNPS.ENSG00000159958.txt", package = "baseqtl", mustWork = TRUE)
gene.coord = system.file("extdata/input", "ENSG00000159958_data.txt", package = "baseqtl", mustWork = TRUE)
vcf = system.file("extdata/input", "chr22GT.86GEU.vcf.gz", package = "baseqtl", mustWork = TRUE)
le.file = system.file("extdata/input", "1000GP_Phase3_subset_chr22.legend.gz", package = "baseqtl", mustWork = TRUE)
h.file = system.file("extdata/input", "1000GP_Phase3_subset_chr22.hap.gz", package = "baseqtl", mustWork = TRUE)
AI_estimate = system.file("extdata/input", "AI_estimate.GT.txt", package = "baseqtl", mustWork = TRUE)
```

``` r
## Choose your output directory
out <- "path/to/output_dir"

## To minimise file sizes and computational load SNPs are within a 10^4 cis-window
## baseqtl is computational intense, it is recommended to be run with several cores
## cores are automatically detected by R

## Run baseqtl.gt:
baseqtl.gt(gene = "ENSG00000159958",
       chr = 22,
       snps = 10^4,
       counts.f = counts.f,
       covariates = covariates,
       additional_cov = NULL,
       e.snps = e.snps,
       u.esnps = u.snps,
       gene.coord = gene.coord,
       vcf = vcf,
       le.file = le.file,
       h.file = h.file,
       out = out,
       AI_estimate = AI_estimate)
```

  - The following messages will appear on screen:

  - max r2 is now0.894478527607361 \#\# relates to the tagging

  - \[1\] “Effective number of fSNPs: 1” (number of exonic SNPs used in
    ASE model)

  - \[1\] “Preparing stan inputs”

  - \[1\] “Running NB\_ASE model” (stan is running NB\_ASE model)

  - \[1\] “Running NB model” (stan is running negative binomial model)

**Output files**

With default ‘prefix’ argument you will find the following files in your
output directory:

1.  ENSG00000159958.GT.eqtl.tags.lookup.txt
2.  ENSG00000159958.GT.excluded.rsnps.txt
3.  ENSG00000159958.GT.stan.summary.txt

<!-- end list -->

``` r
## Inspecting output files
tags <- data.table::fread(system.file("extdata/output", "ENSG00000159958.GT.eqtl.tags.lookup.txt", package = "baseqtl", mustWork = TRUE))
excluded.rsnps <- data.table::fread(system.file("extdata/output", "ENSG00000159958.GT.excluded.rsnps.txt", package = "baseqtl", mustWork = TRUE))
results <- data.table::fread(system.file("extdata/output", "ENSG00000159958.GT.stan.summary.txt", package = "baseqtl", mustWork = TRUE))
```

**ENSG00000159958.GT.eqtl.tags.lookup.txt**

When tagging: column tag is the tag SNP, column SNP are the tagged SNPs.
SNPs are identified by position:reference allele:alternative allele

| Gene\_id        | tag          | SNP          |
| :-------------- | :----------- | :----------- |
| ENSG00000159958 | 42312807:G:A | 42311058:T:G |
| ENSG00000159958 | 42312807:G:A | 42312807:G:A |
| ENSG00000159958 | 42312807:G:A | 42312938:A:C |
| ENSG00000159958 | 42312807:G:A | 42313324:C:A |
| ENSG00000159958 | 42312807:G:A | 42313995:C:T |
| ENSG00000159958 | 42312807:G:A | 42314124:C:T |

**ENSG00000159958.GT.excluded.rsnps.txt**

The table details the SNPs that were excluded fron running NB and/or
NB\_ASE models and the reason for exclusion.

  - Reasons for total exclusion are:
    1.  Missing gentoypes in all samples or homozygous in all samples
    2.  Snp with zero variance
    3.  rsnp with less than ‘nhet’ het ind. ‘nhet’ is the nhet argument
        in baseqtl.gt
  - Reasons to exclude from NB\_ASE, if model = “both” those SNPs will
    be run with NB model
    1.  No entry for gene “gene” in e.snps
    2.  No rsnp in reference panel
    3.  Not enough individuals with ASE counts
    4.  Not unique fsnps in gene
    5.  No fSNPs with GT or AI estimates (when applying reference panel
        bias correction)

| id           | reason                         |
| :----------- | :----------------------------- |
| 42311062:A:C | Missing or homo GT all samples |
| 42311071:C:G | Missing or homo GT all samples |
| 42311196:C:G | Missing or homo GT all samples |
| 42311285:G:C | Missing or homo GT all samples |
| 42311309:A:G | Missing or homo GT all samples |
| 42311332:C:T | Missing or homo GT all samples |

**ENSG00000159958.GT.stan.summary.txt**

The summary file has the following information:

``` r
names(results)
#>  [1] "Gene_id"          "tag"              "log2_aFC_mean"    "log2_aFC_se_mean"
#>  [5] "log2_aFC_sd"      "log2_aFC_0.5%"    "log2_aFC_2.5%"    "log2_aFC_25%"    
#>  [9] "log2_aFC_50%"     "log2_aFC_75%"     "log2_aFC_97.5%"   "log2_aFC_99.5%"  
#> [13] "log2_aFC_d"       "null.99"          "Signif"           "n_eff"           
#> [17] "Rhat"             "model"            "nhets"            "ASE.hets"        
#> [21] "tag.EAF"          "n.fsnps"          "PEP"              "min_AI"
```

  - Description:
      - Gene\_id: ensembl gene id
      - tag or SNP id: position:reference:alternative allele, column
        name tag when tagging, SNP otherwise. Corresponds to the
        regulatory SNP (cis-SNP).
      - log2\_aFC\_mean: eQTL effect expressed as log2 allelic fold
        change (alternative/reference allele): corresponds to the
        posterior mean
      - log2\_aFC\_se\_mean: Monte Carlo standard error to assess model
        performance (see stan manual for more details)
      - log2\_aFC\_sd: eQTL posterior standard deviation
      - log2\_aFC\_0.5% - log2\_aFC\_99.5%: quantiles for the posterior
        distribution
      - null.99: “yes” when 0 (null) is within the 99% credible
        interval, “no” otherwise
      - Signif: “no” if 0 falls within 99% of the eQTL effect posterior
        distribution, “yes” otherwise
      - log2\_aFC\_d: when eQTL effect is significant is the distance
        from the closest quantile (0.5% or 99.5%) to the null. When eQTl
        effect is not significant is the width of the credible interval
      - n\_eff: effective sample size (see stan manual for more details)
      - Rhat: R-hat statisitc (it is recommended only results if
        *\(Rhat < 1.1\)*, see stan manual for more details)
      - model: whether NB or NB-ASE was run
      - nhets: numbr of hets for rSNP
      - ASE.hets: when running NB-ASE number of hets for rSNP with
        sufficient ASE counts, otherwise NA
      - tag.EAF/ SNP.EAF: EAF for rSNP based on the external reference
        panel for the population requested, NA if the rSNP is not in the
        reference panel
      - n.fsnps: number of exonic SNPs used in model, NA for NB model
        \*PEP: posterior exclusion probablity based on 4000 posterior
        draws. Gives the proportion of posterior draws of opposite sign
        to the posterior mean
      - min\_AI: when using reference panel bias correcion gives the
        most extreme allelic imbalance estimate for the exonic SNPs used
        in model (no imbalance corresponds to 0.5), NA
otherwise

| Gene\_id        | tag          | log2\_aFC\_mean | log2\_aFC\_se\_mean | log2\_aFC\_sd | log2\_aFC\_0.5% | log2\_aFC\_2.5% |
| :-------------- | :----------- | --------------: | ------------------: | ------------: | --------------: | --------------: |
| ENSG00000159958 | 42311097:G:A |     \-0.0029311 |           0.0005074 |     0.0403060 |     \-0.1129825 |     \-0.0842575 |
| ENSG00000159958 | 42311204:G:T |     \-0.0021680 |           0.0016226 |     0.0576564 |     \-0.1930193 |     \-0.0947769 |

| log2\_aFC\_25% | log2\_aFC\_50% | log2\_aFC\_75% | log2\_aFC\_97.5% | log2\_aFC\_99.5% | log2\_aFC\_d | null.99 | Signif |
| -------------: | -------------: | -------------: | ---------------: | ---------------: | -----------: | :------ | :----- |
|    \-0.0293565 |    \-0.0029378 |      0.0230779 |        0.0782694 |        0.1073190 |    0.2203016 | yes     | no     |
|    \-0.0320624 |    \-0.0020506 |      0.0299898 |        0.0902738 |        0.1317523 |    0.3247716 | yes     | no     |

|   n\_eff |      Rhat | model  | nhets | ASE.hets |   tag.EAF | n.fsnps |     PEP |   min\_AI |
| -------: | --------: | :----- | ----: | -------: | --------: | :------ | ------: | --------: |
| 6309.026 | 0.9998788 | NB-ASE |    27 |       14 | 0.2465209 | 1       | 0.46775 | 0.4878554 |
| 1262.605 | 1.0004557 | NB     |     8 |       NA | 0.0328032 | NA      | 0.48350 |        NA |

## eQTL analysis when only RNA-seq is available

The function to call is **baseqtl.nogt**. To show an example we provide
data files for gene ENSG00000159958 on chromosome 22. I will describe
below the arguments that differ from baseqtl.gt.

**Arguments**

*vcf*: full name for vcf file with genotypes and ASE counts for exonic
SNPS. This file is produced from RNA-seq data only as explained in
[snakefile](https://gitlab.com/evigorito/baseqtl_pipeline/-/blob/master/input/Snakefile)
output from rule merge\_vcf, source=“RNA”\*\*

*info*: cut-off based on the imputation quality of cis-SNP. By default
the model is only run for those SNPs with info\(>= 0.3\)

*ex.fsnp*: if numeric p-value cut-off for fSNPs to exclude based on
fisher test for suspected genotype error based on comparing the
proportion of hets in the sample and reference panel, defaults to 0.01.
If character: vector with pos:ref:alt for fsnps to exclude (same as
baseqtl.gt),if NULL no fSNP will be excluded.

*sample.file*: sample file for the reference panel (sample description),
to be used if ex.fsnp test is numeric and population is not the whole
reference panel

*save\_input*: whether to save input to stan model for QC purposes,
defaults to FALSE to save disk space. Object ending with
“noGT.stan.input.rds” is a named list with each element the inputs for
a cis-SNP. For each cis-SNP there is a list of 2 elements: “NB” and
“ase”. “NB” is a list with elements “counts” and “p.g”. “Counts” is
a data.table with columns sample names and one row corresponding to the
gene, values total read counts. “p.g” is a named list with each element
a sample. For each sample there is an array with names genotypes (0,1,2)
and values the genotype probabilities. For the “ase” list they are for
elements: “m” numeric vector with total ASE counts per sample. “g” list
with each element a sample and for each sample the genoptype of the cis
SNP coded as 0,1,2 and -1, with -1 indicating that the alternative
allele is in haplotype 1. “p” has the same structure as “g” and
indicates the probability for each genotype. “n” is similar to “g” and
“p” but contains the mapped reads to haplotype 2. The file ending with
“noGT.fsnps.counts.rds is a matrix with rows samples and columns fSNPS.
When a fSNPs ends with”.n" correspond to the counts matching the
alternative allele and “.m” indicates the total counts matching the SNP.

``` r
## Retrive input files for running baseqtl.nogt

counts.f <- system.file("extdata/input", "counts.txt", package = "baseqtl", mustWork = TRUE)
covariates <- system.file("extdata/input", "lbsize_gc.rds", package = "baseqtl", mustWork = TRUE)
e.snps <- system.file("extdata/input", "chr22.fSNPS.ENSG00000159958.txt", package = "baseqtl", mustWork = TRUE)
u.snps <- system.file("extdata/input", "chr22.unique.fSNPS.ENSG00000159958.txt", package = "baseqtl", mustWork = TRUE)
gene.coord <- system.file("extdata/input", "ENSG00000159958_data.txt", package = "baseqtl", mustWork = TRUE)
vcf = system.file("extdata/input", "chr22noGT.86GEU.vcf.gz", package = "baseqtl", mustWork = TRUE)
sample.f <- system.file("extdata/input", "1000GP_Phase3.sample", package = "baseqtl", mustWork = TRUE)
le.file <- system.file("extdata/input", "1000GP_Phase3_subset_chr22.legend.gz", package = "baseqtl", mustWork = TRUE)
h.file <- system.file("extdata/input", "1000GP_Phase3_subset_chr22.hap.gz", package = "baseqtl", mustWork = TRUE)
AI_estimate <- system.file("extdata/input", "AI_estimate.noGT.txt", package = "baseqtl", mustWork = TRUE)
```

``` r
## Choose your output directory
out <- "path/to/output_dir"

## To minimise file sizes and computational load SNPs are within a 10^4 cis-window
## baseqtl is computational intense, it is recommended to be run with several cores
## cores are automatically detected by R

## Run baseqtl.nogt:
baseqtl.nogt(gene = "ENSG00000159958",
       chr = 22,
       snps = 10^4,
       counts.f = counts.f,
       covariates = covariates,
       additional_cov = NULL,
       e.snps = e.snps,
       u.esnps = u.snps,
       gene.coord = gene.coord,
       vcf = vcf,
       sample.f=sample.f,
       le.file = le.file,
       h.file = h.file,
       out = out,
       AI_estimate = AI_estimate)
```

**Output files**

With default ‘prefix’ argument you will find the following files in your
output directory:

1.  ENSG00000159958.noGT.eqtl.tags.lookup.txt, same format as
    ENSG00000159958.GT.eqtl.tags.lookup.txt explained above.
2.  ENSG00000159958.noGT.excluded.rsnps.txt, same format as
    ENSG00000159958.GT.excluded.rsnps.txt explained above.
3.  ENSG00000159958.fsnps.het.fisher.test.txt
4.  ENSG00000159958.noGT.stan.summary.txt

<!-- end list -->

``` r
## Inspecting output files
fsnps <- data.table::fread(system.file("extdata/output", "ENSG00000159958.fsnps.het.fisher.test.txt", package = "baseqtl", mustWork = TRUE))
results <- data.table::fread(system.file("extdata/output", "ENSG00000159958.noGT.stan.summary.txt", package = "baseqtl", mustWork = TRUE))
```

**ENSG00000159958.fsnps.het.fisher.test.txt**

Table with feature SNP (exonic SNP) id (position:reference:altrnative
allele), odds ratio (OR) and pvalue testing the frequencey of
heterozygocity between sample and reference panel and ensembl gene
id.

``` r
fsnps <- data.table::fread(system.file("extdata/output", "ENSG00000159958.fsnps.het.fisher.test.txt", package = "baseqtl", mustWork = TRUE))
```

| fsnp         |       OR |    pvalue | gene\_id        |
| :----------- | -------: | --------: | :-------------- |
| 42321251:A:G | 1.291604 | 0.3259933 | ENSG00000159958 |

**ENSG00000159958.noGT.stan.summary.txt**

The summary file has the following information:

``` r
names(results)
#>  [1] "Gene_id"          "tag"              "log2_aFC_mean"    "log2_aFC_se_mean"
#>  [5] "log2_aFC_sd"      "log2_aFC_0.5%"    "log2_aFC_2.5%"    "log2_aFC_25%"    
#>  [9] "log2_aFC_50%"     "log2_aFC_75%"     "log2_aFC_97.5%"   "log2_aFC_99.5%"  
#> [13] "log2_aFC_d"       "null.99"          "Signif"           "n_eff"           
#> [17] "Rhat"             "model"            "nhets"            "ASE.hets"        
#> [21] "tag.EAF"          "info"             "n.fsnps"          "min.p.fsnp"      
#> [25] "PEP"              "min_AI"
```

  - Description: same as section **ENSG00000159958.GT.stan.summary.txt**
    except:
      - info: quality of imputation for rSNP.
      - min.p.fsnp: minumun pvalue for Fisher test of heterozygocity
        across all fSNPs (exonic
SNPS)

## eQTL analysis with paired samples (two treatments on same samples) with genotypes.

The function to call is **baseqtl.gt.paired**. It will test for eQTL
interaction effect between 2 conditions in a pair design. We use again
gene ENSG00000159958 on chromosome 22.

**Arguments**

Same as baseqtl.gt except:

*count.f*: vector with full name for total gene counts for each
treatment. First column gene id followed by samples. Samples in same
order and same name in both files, details in **input4baseqtl**.

*vcf*: vector with full name to vcf files with GT and ASE counts.
Genotype field should be the same in both files. Order of treatments
should be the same as in count.f, details in **input4baseqtl**.

``` r
## Retrive input files for running baseqtl.gt
## Most of the files contain information for the selected gene only to minimise their size
## For simplicity I will use the same counts.f and vcf files twice, but in reality each file
## will match a different treatment or condition

counts.f <- rep(system.file("extdata/input", "counts.txt", package = "baseqtl", mustWork = TRUE),2)
covariates <- system.file("extdata/input", "lbsize_gc.rds", package = "baseqtl", mustWork = TRUE)
e.snps = system.file("extdata/input", "chr22.fSNPS.ENSG00000159958.txt", package = "baseqtl", mustWork = TRUE)
u.snps = system.file("extdata/input", "chr22.unique.fSNPS.ENSG00000159958.txt", package = "baseqtl", mustWork = TRUE)
gene.coord = system.file("extdata/input", "ENSG00000159958_data.txt", package = "baseqtl", mustWork = TRUE)
vcf = rep(system.file("extdata/input", "chr22GT.86GEU.vcf.gz", package = "baseqtl", mustWork = TRUE), 2)
le.file = system.file("extdata/input", "1000GP_Phase3_subset_chr22.legend.gz", package = "baseqtl", mustWork = TRUE)
h.file = system.file("extdata/input", "1000GP_Phase3_subset_chr22.hap.gz", package = "baseqtl", mustWork = TRUE)
AI_estimate = system.file("extdata/input", "AI_estimate.GT.txt", package = "baseqtl", mustWork = TRUE)
```

``` r
## Choose your output directory
out <- "path/to/output_dir"

## To minimise file sizes and computational load SNPs are within a 10^4 cis-window
## baseqtl is computational intense, it is recommended to be run with several cores
## cores are automatically detected by R

## Run baseqtl.gt.paired:
baseqtl.gt.paired(gene = "ENSG00000159958",
       chr = 22,
       snps = 10^4,
       counts.f = counts.f,
       covariates = covariates,
       additional_cov = NULL,
       e.snps = e.snps,
       u.esnps = u.snps,
       gene.coord = gene.coord,
       vcf = vcf,
       le.file = le.file,
       h.file = h.file,
       out = out,
       AI_estimate = AI_estimate)
```

**Output files**

With default ‘prefix’ argument you will find the following files in your
output directory:

1.  ENSG00000159958.paired.GT.eqtl.tags.lookup.txt, same format as
    ENSG00000159958.GT.eqtl.tags.lookup.txt explained above.
2.  ENSG00000159958.paired.GT.excluded.rsnps.txt, same format as
    ENSG00000159958.GT.excluded.rsnps.txt explained above.
3.  ENSG00000159958.paired.GT.stan.summary.txt

<!-- end list -->

``` r
## Inspecting output files
results <- data.table::fread(system.file("extdata/output", "ENSG00000159958.paired.GT.stan.summary.txt", package = "baseqtl", mustWork = TRUE))
```

**ENSG00000159958.paired.GT.stan.summary.txt**

The summary file has the following information:

``` r
names(results)
#>  [1] "Gene_id"              "tag"                  "log2_aFC_mean.ba"    
#>  [4] "log2_aFC_se_mean.ba"  "log2_aFC_sd.ba"       "log2_aFC_0.5%.ba"    
#>  [7] "log2_aFC_2.5%.ba"     "log2_aFC_25%.ba"      "log2_aFC_50%.ba"     
#> [10] "log2_aFC_75%.ba"      "log2_aFC_97.5%.ba"    "log2_aFC_99.5%.ba"   
#> [13] "n_eff.ba"             "Rhat.ba"              "Signif.ba"           
#> [16] "PEP.ba"               "log2_aFC_mean.bd"     "log2_aFC_se_mean.bd" 
#> [19] "log2_aFC_sd.bd"       "log2_aFC_0.5%.bd"     "log2_aFC_2.5%.bd"    
#> [22] "log2_aFC_25%.bd"      "log2_aFC_50%.bd"      "log2_aFC_75%.bd"     
#> [25] "log2_aFC_97.5%.bd"    "log2_aFC_99.5%.bd"    "n_eff.bd"            
#> [28] "Rhat.bd"              "Signif.bd"            "PEP.bd"              
#> [31] "log2_aFC_mean.bt1"    "log2_aFC_se_mean.bt1" "log2_aFC_sd.bt1"     
#> [34] "log2_aFC_0.5%.bt1"    "log2_aFC_2.5%.bt1"    "log2_aFC_25%.bt1"    
#> [37] "log2_aFC_50%.bt1"     "log2_aFC_75%.bt1"     "log2_aFC_97.5%.bt1"  
#> [40] "log2_aFC_99.5%.bt1"   "n_eff.bt1"            "Rhat.bt1"            
#> [43] "Signif.bt1"           "PEP.bt1"              "log2_aFC_mean.bt2"   
#> [46] "log2_aFC_se_mean.bt2" "log2_aFC_sd.bt2"      "log2_aFC_0.5%.bt2"   
#> [49] "log2_aFC_2.5%.bt2"    "log2_aFC_25%.bt2"     "log2_aFC_50%.bt2"    
#> [52] "log2_aFC_75%.bt2"     "log2_aFC_97.5%.bt2"   "log2_aFC_99.5%.bt2"  
#> [55] "n_eff.bt2"            "Rhat.bt2"             "Signif.bt2"          
#> [58] "PEP.bt2"              "tag.EAF"              "n.fsnps"             
#> [61] "model"                "nhets"                "ASE.hets"
```

  - **Description**
    
      - Similar output as baseqtl.gt, except that now we look at 4
        coefficients: ba, bd, bt1 and bt2.
      - bt1 and bt2 corespond to the two treatments/condition
        respectively, as ordered in counts.f and vcf inputs.
      - ba corresponds to the ‘addition’ coefficient, \(ba = bt1 + bt2\)
      - bd corresponds to the ‘difference’ coefficient,
        \(bd = bt1 - bt2\)
      - We evaluate ‘ba’ and ‘bd’ to look for a condition specific
        effect
      - When ‘ba’ is significant implies a significant eQTL efect in one
        or both treatments. If ‘bd’ is significant, there is evidence
        for condition specific effect. The most common scenario is a
        significant eQTL effect in only one condition, look at ‘bt1’ and
        ‘bt2’.
      - If ‘ba’ is not significant but ‘bd’ is significant it could
        indicate eQTL effects in opposite directions, look at ‘bt1’ and
        ‘bt2’.
      - ASE.hets gives the number of hets individuals for the rSNP with
        sufficient ASE counts for each treatment as
‘14,10’

## eQTL analysis with samples from two conditions (diseases vs control) and unkown gentopyes.

In this example we are going to compare psoriasis vs normal skin,
RNA-seqdata from [psoriasis](ftp://ftp.sra.ebi.ac.uk/vol1/fastq) We call
**baseqtl2T.nogt**. I will describe below the arguments that differ from
baseqtl.nogt.

**Arguments**

*counts.f*: character vector with names of files with total gene counts
for each treatment

*covariates*:character vector with names of files with covariates for
each treatment

*additional\_cov*: character vector with names of files with covariates
for each treatment. See description above, defaults to NULL

*vcf*: character vector with names of fvcf iles with GT and ASE for
fSNPs for each treatment

*treatments* character vector with the 2 treatments (tissues, diseases)
to study

**File order in counts.f, covariates, vcf and treatment names must be
the same**

*fishjoin* whether to run Fisher test for heterozygocity frequency
between samples and reference panel fSNPs (ex.fsnp argument) jointly in
all samples or by treatment, defaults to jointly but for QC purposes, to
use the same fSNPs in ech treatment as the ones used when the treatments
were run using independent models select NULL.

``` r
## Retrive input files for running baseqtl2T.nogt

counts.f <- c(system.file("extdata/input", "counts_Psoriasis_skin.txt", package = "baseqtl", mustWork = TRUE),
     system.file("extdata/input", "counts_normal_skin.txt", package = "baseqtl", mustWork = TRUE))
     
covariates <- c(system.file("extdata/input", "Psoriasis_skin_gc_lib_size.rds", package = "baseqtl", mustWork = TRUE),
       system.file("extdata/input", "normal_skin_gc_lib_size.rds", package = "baseqtl", mustWork = TRUE))

vcf = c(system.file("extdata/input", "chr10.ASE.Psoriasis_skin.vcf.gz", package = "baseqtl", mustWork = TRUE),
    system.file("extdata/input", "chr10.ASE.normal_skin.vcf.gz", package = "baseqtl", mustWork = TRUE))
    
e.snps <- system.file("extdata/input", "chr10.fSNPS.ENSG00000178372.txt", package = "baseqtl", mustWork = TRUE)
u.snps <- system.file("extdata/input", "chr10.unique.fSNPS.ENSG00000178372.txt", package = "baseqtl", mustWork = TRUE)
gene.coord <- system.file("extdata/input", "ENSG00000178372_data.txt", package = "baseqtl", mustWork = TRUE)
sample.f <- system.file("extdata/input", "1000GP_Phase3.sample", package = "baseqtl", mustWork = TRUE)
le.file <- system.file("extdata/input", "1000GP_Phase3_subset_chr10.legend.gz", package = "baseqtl", mustWork = TRUE)
h.file <- system.file("extdata/input", "1000GP_Phase3_subset_chr10.hap.gz", package = "baseqtl", mustWork = TRUE)
AI_estimate <- system.file("extdata/input", "AI_estimate.psoriasis.txt", package = "baseqtl", mustWork = TRUE)
```

``` r
## Choose your output directory
out <- "path/to/output_dir"

## To minimise file sizes and computational load SNPs are within a 10^4 cis-window
## baseqtl is computational intense, it is recommended to be run with several cores
## cores are automatically detected by R

## Run baseqtl.n2Togt:
baseqtl2T.nogt(gene = "ENSG00000178372",
       chr = 10,
       snps = 10^4,
       counts.f = counts.f,
       covariates = covariates,
       additional_cov = NULL,
       e.snps = e.snps,
       u.esnps = u.snps,
       gene.coord = gene.coord,
       vcf = vcf,
       sample.f=sample.f,
       le.file = le.file,
       h.file = h.file,
       out = out,
       treatment=c("Psoriasis_skin","normal_skin"),
       AI_estimate = AI_estimate)
```

  - The following messages will show on the screen:
  - max r2 is now0.884460446665381
  - \[1\] “Preparing stan inputs”
  - \[1\] “Effective number of exonic SNPs: Psoriasis\_skin 3”
  - \[2\] “Effective number of exonic SNPs: normal\_skin 3”  
  - Running stan for 9 rSNPS

**Output files**

With default ‘prefix’ argument you will find the following files in your
output directory:

1.  ENSG00000178372.2T.noGT.eqtl.tags.lookup.txt, same format as
    ENSG00000159958.GT.eqtl.tags.lookup.txt explained above.
2.  ENSG00000178372.2T.noGT.excluded.snps.txt, same format as
    ENSG00000159958.GT.excluded.rsnps.txt explained above.
3.  ENSG00000178372.2T.fsnps.het.fisher.test.txt same format as
    ENSG00000159958.fsnps.het.fisher.test.txt
4.  ENSG00000178372.2T.noGT.summary.txt same format as
    ENSG00000159958.paired.GT.stan.summary.txt

All ouput files are available at “extdata/output”.
