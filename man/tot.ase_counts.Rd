% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/real.data.R
\name{tot.ase_counts}
\alias{tot.ase_counts}
\title{Total gene and ASE counts, per fsnp, per individual}
\usage{
tot.ase_counts(x, y = NULL, z = NULL)
}
\arguments{
\item{x}{DT with ASE and GT created from reading vcf}

\item{y}{data table with total counts for samples}

\item{z}{data table with each row the genotype for 1 rsnp coded as 0,1,-1 or 2, output from a rec_mytrecase_rSNPs}
}
\value{
matrix if z=NULL or list of  data tables, each data table corresponds to each rsnp, cols are total counts (y), GT 0,1,-1,2 for the rSNP, ase counts per fsnps across samples
tot.ase_counts
}
\description{
Get total and AS counts per snp per individual
}
\keyword{ASE}
\keyword{counts}
\keyword{gene}
