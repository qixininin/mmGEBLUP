#' calculateKaXa function
#'
#' @param geno_data genotype data
#' @param pheno_data phenotype data
#' @param A additive relationship matrix
#' @param site_qtl a vector, names of QTL sites
#' @param m_qtl the number of QTLs
#'
#' @return list(Ka = Ka,Xa = Xa)
#' @export
#' @import dplyr
#'
#' @examples \dontrun{KaXa = calculateKaXa(mmgeno_data, mmpheno_data, A, site_qtl, m_qtl)}
calculateKaXa <- function(geno_data, pheno_data, A, site_qtl, m_qtl)
{
  # extract qtl information from what we got from GWAS
  if(m_qtl>0)
  {
    Ka = A.mat(as.matrix(geno_data[,!colnames(geno_data) %in% site_qtl]))
    Xa = geno_data[as.character(pheno_data$GID), site_qtl]
  } else {
    Ka = A
    Xa = c()
  }

  colnames(Ka) = rownames(Ka) = rownames(geno_data)

  return(list(Ka = Ka,Xa = Xa))
}
