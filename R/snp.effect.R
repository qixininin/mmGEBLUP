#' snp.effect function
#'
#' @param snpNum The number of SNPs
#' @param envNum The number of environments
#' @param major_a_idx The index for major additive SNP
#' @param major_ae_idx The index for major additive-by-environment SNP
#' @param variance_a_major The variance for major additive SNP effect
#' @param variance_ae_major The variance for major additive-by-environment SNP effect
#' @param variance_a_minor The variance for minor additive SNP effect
#' @param variance_ae_minor The variance for minor additive-by-environment SNP effect
#'
#' @return a list
#'
#' $effects an envNum * snpNum matrix for total SNP effect
#' $main_effects a 1 * snpNum matrix for main SNP effect
#' $interaction_effects a envNum * snpNum matrix for interaction SNP effect
#'
#' @export
#' @importFrom stats rnorm
#'
#' @examples snp.effect(snpNum = 2000, envNum = 3,
#'                      major_a_idx = c(500, 750, 1000, 1250, 1500),
#'                      major_ae_idx = c(250, 500, 1000, 1500, 1750),
#'                      variance_a_major = 0.02,  variance_ae_major = 0.01,
#'                      variance_a_minor = 0.002,  variance_ae_minor = 0.001)
snp.effect <- function(snpNum, envNum, major_a_idx, major_ae_idx,
                       variance_a_major, variance_ae_major,
                       variance_a_minor, variance_ae_minor) {

  if(max(major_a_idx)>snpNum | max(major_ae_idx)>snpNum) {
    stop("Error: snp.effect(). The input major_a_idx or major_ae_idx is out of the snpNum range.")
  }

  # determine minor index
  minor_a_idx = setdiff(1:snpNum, major_a_idx)
  minor_ae_idx = setdiff(1:snpNum, major_ae_idx)

  # store main and interaction effect for all snps
  main_effects     <- matrix(0, nrow = 1     , ncol = snpNum)
  interact_effects <- matrix(0, nrow = envNum, ncol = snpNum)

  # Simulate main gene effect

  main_effects[1,major_a_idx] <- rnorm(length(major_a_idx), mean = 0, sd = sqrt(variance_a_major))
  main_effects[1,minor_a_idx] <- rnorm(length(minor_a_idx), mean = 0, sd = sqrt(variance_a_minor))


  for (i in 1:envNum) {
    interact_effects[i,major_ae_idx] <- rnorm(length(major_ae_idx), mean = 0, sd = sqrt(variance_ae_major))
    interact_effects[i,minor_ae_idx] <- rnorm(length(minor_ae_idx), mean = 0, sd = sqrt(variance_ae_minor))
  }

  effects = main_effects[rep(1,envNum),] + interact_effects

  # Return effects matrix
  return(list(effects = effects,
              main_effects = main_effects,
              interact_effects = interact_effects))
}
