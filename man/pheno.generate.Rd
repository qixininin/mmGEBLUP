% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pheno.generate.R
\name{pheno.generate}
\alias{pheno.generate}
\title{pheno.generate function}
\usage{
pheno.generate(genotypes, effects, envNum, indNum, sigma.error)
}
\arguments{
\item{genotypes}{an N*M genotype matrix, with N rows of individual and M columns of markers, individual ID and marker ID are given as rownames and colnames.}

\item{effects}{an H*M effect matrix, with H rows of environment and M columns of markers.}

\item{envNum}{the number of environments}

\item{indNum}{the number of individuals}

\item{sigma.error}{the variance for residual effects.}
}
\value{
pheno_data a data frame, with three columns $ENV, $GID, and $SimTrait
}
\description{
pheno.generate function
}
\examples{
\dontrun{pheno.generate(genotypes = t(geno_data[-c(1:3)]), effects = b,
                   envNum = envNum, indNum = indNum, sigma.error = sigma_error)}
}
