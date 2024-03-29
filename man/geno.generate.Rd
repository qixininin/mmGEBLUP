% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geno.generate.R
\name{geno.generate}
\alias{geno.generate}
\title{geno.generate function
A simple function to generate homozygotic genotype data (coded with 1 or -1)
with respect to the number of individual, the number of markers, MAF range and chromosome.
Markers are independent.}
\usage{
geno.generate(indNum, snpNum, maf.min, maf.max, chr.snpNum)
}
\arguments{
\item{indNum}{The number of individual}

\item{snpNum}{The number of markers}

\item{maf.min}{The lower bound of MAF, so that MAF is sampled from U(maf.min, maf.max)}

\item{maf.max}{The upper bound of MAF, so that MAF is sampled from U(maf.min, maf.max)}

\item{chr.snpNum}{The number of snps on each chromosome}
}
\value{
geno_data a data frame, where the first three columns are $CHR, $SNP, and $BP
}
\description{
geno.generate function
A simple function to generate homozygotic genotype data (coded with 1 or -1)
with respect to the number of individual, the number of markers, MAF range and chromosome.
Markers are independent.
}
\examples{
geno.generate(100, 200, 0.05, 0.5, c(50,50,40,60))
}
