% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snp.effect.R
\name{snp.effect}
\alias{snp.effect}
\title{snp.effect function}
\usage{
snp.effect(
  snpNum,
  envNum,
  major_a_idx,
  major_ae_idx,
  variance_a_major,
  variance_ae_major,
  variance_a_minor,
  variance_ae_minor
)
}
\arguments{
\item{snpNum}{The number of SNPs}

\item{envNum}{The number of environments}

\item{major_a_idx}{The index for major additive SNP}

\item{major_ae_idx}{The index for major additive-by-environment SNP}

\item{variance_a_major}{The variance for major additive SNP effect}

\item{variance_ae_major}{The variance for major additive-by-environment SNP effect}

\item{variance_a_minor}{The variance for minor additive SNP effect}

\item{variance_ae_minor}{The variance for minor additive-by-environment SNP effect}
}
\value{
a list

$effects an envNum * snpNum matrix for total SNP effect
$main_effects a 1 * snpNum matrix for main SNP effect
$interaction_effects a envNum * snpNum matrix for interaction SNP effect
}
\description{
snp.effect function
}
\examples{
snp.effect(snpNum = 2000, envNum = 3,
                     major_a_idx = c(500, 750, 1000, 1250, 1500),
                     major_ae_idx = c(250, 500, 1000, 1500, 1750),
                     variance_a_major = 0.02,  variance_ae_major = 0.01,
                     variance_a_minor = 0.002,  variance_ae_minor = 0.001)
}
