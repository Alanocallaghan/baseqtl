% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_baseqtl_GT.R
\name{aux.snps}
\alias{aux.snps}
\title{get gt.ase info given gene and snps}
\usage{
aux.snps(gene, chr, snps, gcoord, gene.coord, vcf, rsnps.ex)
}
\arguments{
\item{gene}{gene id for the gene to run}

\item{chr}{chromosome where the gene is, example chr=22}

\item{snps}{either cis-window or character vector with pos:ref:alt allele for each snp, defaults to cis-window}

\item{gcoord}{data table from reading gene.coord}

\item{gene.coord}{path to file listing gene coordinates and exons, prepared in input.R}

\item{vcf}{path to vcf file with ASE and GT for the chromosome where the gene is}

\item{rsnps.ex}{data table with id and reason column to exclude rsnps from analysis}
}
\value{
list with gt.ase info for whole cis-window to contain fSNPs, rs is the genotype plus rsnp info, excluded rsnps and cis-window
}
\description{
get gt.ase info given gene and snps
}
