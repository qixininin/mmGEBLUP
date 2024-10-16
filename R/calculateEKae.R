#' calculateEKae function
#'
#' @param geno_data genotype data
#' @param E the environment covariance matrix
#' @param site_env_qtl_all a dataframe with all AE QTL
#'
#' @return a list of KaeE1 and KaeE2
#' @export
#'
#' @examples \dontrun{EKae = calculateEKae(geno_data = mmgeno_data, E, site_env_qtl_all)}
calculateEKae <- function(geno_data, E, site_env_qtl_all)
{

  site_qe = unique(site_env_qtl_all$QTL)

  ## Variance-covariance matrix in terms of small AE effects
  Kae_s = A.mat(as.matrix(geno_data[,!colnames(geno_data) %in% site_qe]))
  colnames(Kae_s) = rownames(Kae_s) = rownames(geno_data)
  KaeE_s = kronecker(Kae_s, E, make.dimnames = T)

  ## Variance-covariance matrix in terms of large AE effects (a list)
  # genotypes for large AE loci
  Xae = as.matrix(geno_data[, colnames(geno_data) %in% site_qe])
  # imputation with mean (- 2024/10/17)
  Xae = apply(Xae, 2, function(x) {
    if(any(is.na(x))){
      x[which(is.na(x))] = mean(x, na.rm = T)
    }
    return(x)
    })
  KaeE_l = list()
  for(t in 1:length(site_qe))
  {
    Kae_l = tcrossprod(Xae[,t])
    colnames(Kae_l) = rownames(Kae_l) = names(Xae[,t])
    KaeE_l[[t]] = kronecker(Kae_l, E, make.dimnames = T)
  }
  KaeE_l = as.list(KaeE_l)

  return(list(KaeE1 = KaeE_l, KaeE2 = KaeE_s))

}
