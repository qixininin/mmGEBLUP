% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmgblup.R
\name{mmgblup}
\alias{mmgblup}
\title{mmgblup function}
\usage{
mmgblup(data, Ka, AE)
}
\arguments{
\item{data}{a data frame}

\item{Ka}{additive genetic relationship matrix}

\item{AE}{additive-by-environment relationship matrix1}
}
\value{
list(mod, BV)
}
\description{
mmgblup function
}
\examples{
\dontrun{rst = mmgblup(data = cbind(dt, mmdata$Xa), Ka = mmdata$Ka, AE = mmdata$AE)}
}
