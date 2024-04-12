#' calculateEKae function
#'
#' @param geno_data genotype data
#' @param A additive relationship matrix
#' @param E environment variance-covariance matrix
#' @param EA additive-by-environment relationship matrix
#' @param site_env_qtl_all a dataframe with all AE QTL
#' @param m_env_qtl the number of AE QTL
#'
#' @return a list of EKae1 and EKae2
#' @export
#'
#' @examples \dontrun{EKae = calculateEKae(geno_data = mmgeno_data, A, E, EA, site_env_qtl_all, m_env_qtl)}
calculateEKae <- function(geno_data, A, E, EA, site_env_qtl_all, m_env_qtl)
{
  EKae_l = EA
  EKae_s = EA
  lineNum = nrow(A)
  trialName = rownames(E)
  trialNum = length(trialName)

  if(m_env_qtl>0)
  {
    for(e in 1:trialNum)
    {
      env = trialName[e]
      left = (e-1)*lineNum+1
      right = e*lineNum
      site_env_qtl = site_env_qtl_all %>% dplyr::filter(ENV == env)

      if(nrow(site_env_qtl)>0)
      {
        # adjusted additive-by-env relationship matrix
        EKae_l[left:right,left:right] = A.mat(as.matrix(geno_data[, colnames(geno_data) %in% site_env_qtl$QTL]))
        EKae_s[left:right,left:right] = A.mat(as.matrix(geno_data[,!colnames(geno_data) %in% site_env_qtl$QTL]))
      } else {
        EKae_l[left:right,left:right] = diag(lineNum)
        EKae_s[left:right,left:right] = A
      }
    }

    ## TEST
    site_qe = unique(site_env_qtl_all$QTL)
    Xae = as.matrix(geno_data[, colnames(geno_data) %in% site_qe])
    EKae_l = list()
    for(t in 1:length(site_qe))
    {
      Ka = tcrossprod(Xae[,t])
      colnames(Ka) = rownames(Ka) = names(Xae[,t])
      EKae_l[[t]] = kronecker(Ka, E, make.dimnames = T)
    }
    EKae_l = as.list(EKae_l)
    ## TEST END

    return(list(EKae1 = EKae_l,
                EKae2 = EKae_s))
  } else {
    return(list(EKae1 = EA,
                EKae2 = NULL))
  }
}
