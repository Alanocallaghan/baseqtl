% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stan.eff.R
\name{stan.paired.for}
\alias{stan.paired.for}
\title{Format stan matrix output for paired design}
\usage{
stan.paired.for(s, snp)
}
\arguments{
\item{s}{stan matrix}

\item{snp}{id for rsnp}
}
\value{
data table with corrected PEP based on eQTL effect sign
stan.paired.for()
}
\description{
This function allows further formatting of stan matrix output when dialing with multiple parameters
}
