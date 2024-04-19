#' mmdata function
#' Format data for mmGEBLUP analysis
#'
#' @param geno_data input genotype data
#' @param pheno_data input phenotype data
#' @param qtl_data input qtl data
#' @param qtl_env_data input environmental qtl data
#'
#' @return a list with QCed genotype and phenotype data, and other data required for GS analysis
#'         $ mmgeno_data an n*m matrix
#'         $ mmpheno_data an n_obs*3 data frame
#'         $ Ka additive relationship matrix for mmGEBLUP
#'         $ Xa additive fixed effect coefficient for mmGEBLUP
#'         $ EKae1 additive-by-environment relationship matrix 1 for mmGEBLUP
#'         $ EKae2 additive-by-environment relationship matrix 2 for mmGEBLUP
#'         $ A additive relationship matrix for GEBLUP
#'         $ AE additive-by-environment relationship matrix for GEBLUP
#'         $ mmsummary data summary
#' @export
#' @import dplyr
#'
#' @examples \dontrun{qcdata = dataqc(geno_data, pheno_data, qtl_data, qtl_env_data)}
mmdata <- function(geno_data, pheno_data, qtl_data, qtl_env_data)
{
  # Transform genotype
  rownames(geno_data) = geno_data[,2]
  mmgeno_data = t(geno_data[,-c(1:3)])

  # Summary
  traitName = colnames(pheno_data)[-c(1:2)]
  traitNum = length(traitName)
  trialName = unique(pheno_data$ENV)
  trialNum = length(trialName)
  lineName = intersect(unique(pheno_data$GID),rownames(mmgeno_data))
  lineNum = length(lineName)

  # Prepare genotype
  mmgeno_data = as.matrix(mmgeno_data[lineName, ])

  # Prepare phenotype
  mmpheno_data = pheno_data %>%
    dplyr::filter(GID %in% lineName) %>%
    droplevels()
  colnames(mmpheno_data)[colnames(mmpheno_data) == traitName] <- deparse(substitute(trait))

  # Prepare qtl
  if(!is.null(qtl_data)){
    site_qtl = qtl_data$QTL
    m_qtl = length(site_qtl)
  } else {
    m_qtl = 0
  }

  # Prepare qtl_env
  if(!is.null(qtl_data)){
    site_env_qtl_all = qtl_env_data %>% dplyr::filter(ENV %in% trialName)
    m_env_qtl = nrow(site_env_qtl_all)
  } else {
    m_env_qtl = 0
  }

  # Summary
  full_summary <- paste0("  -------------------------------------------------------------", "\n",
                         "                           Data Summary                        ", "\n",
                         "  -------------------------------------------------------------", "\n",
                         "  **The number of trials/environments: ", trialNum, "\n",
                         "  **The number of lines:               ", lineNum,  "\n",
                         "  **The number of markers:             ", ncol(mmgeno_data),  "\n",
                         "  **The number of observation:         ", nrow(mmpheno_data),  "\n","\n",

                         "  **Trait Name: ", traitName,  "\n",
                         "  **Trial Name: ", paste(trialName, collapse = ", "),  "\n",
                         "  **Line  Name: ", paste(lineName, collapse = ", "),  "\n"
  )



  # Additive relationship matrix
  A = A.mat(mmgeno_data)
  colnames(A) = rownames(A) = rownames(mmgeno_data)
  E = diag(trialNum)
  rownames(E) = colnames(E) = trialName
  AE = kronecker(A, E, make.dimnames = TRUE)

  # Calculate reduced additive relationship matrix
  KaXa = calculateKaXa(mmgeno_data, mmpheno_data, A, site_qtl, m_qtl)

  # Calculate two additive-by-environment matrices
  if(m_env_qtl>0){
    KaeE = calculateEKae(mmgeno_data, E, site_env_qtl_all)
  } else {
    KaeE = list(KaeE1 = AE, KaeE2 = NULL)
  }

  return(list(mmgeno_data = mmgeno_data,
              mmpheno_data = mmpheno_data,
              Ka = KaXa$Ka,
              Xa = KaXa$Xa,
              KaeE1 = KaeE$KaeE1,
              KaeE2 = KaeE$KaeE2,
              A = A,
              AE = AE,
              mmsummary = list(traitName = traitName,
                               trialName = trialName,
                               lineName  = lineName)))
}
