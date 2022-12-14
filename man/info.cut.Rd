% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/real.data.R
\name{info.cut}
\alias{info.cut}
\title{wrap function to calculate r2=var(E(G))/var(G) to use for info input for baseqtl noGT}
\usage{
info.cut(stan.noGT, rp.r, info)
}
\arguments{
\item{stan.noGT}{input list for one gene}

\item{rp.r}{matrix with reference panel haplotype info, cols individuals, rows snps}

\item{info, }{cut-off to remove snps, when r2<info, remove}
}
\value{
named vector with r2 and names snp_id, r2 above threshold
info.cut()
}
\description{
This function allows you to to calculate select snps with r2=var(E(G))/var(G) above specific threshold
}
