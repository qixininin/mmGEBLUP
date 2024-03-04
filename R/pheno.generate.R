#' pheno.generate function
#'
#' @param genotypes an N*M genotype matrix, with N rows of individual and M columns of markers, individual ID and marker ID are given as rownames and colnames.
#' @param effects an H*M effect matrix, with H rows of environment and M columns of markers.
#' @param envNum the number of environments
#' @param indNum the number of individuals
#' @param sigma.error the variance for residual effects.
#'
#' @return pheno_data a data frame, with three columns $ENV, $GID, and $SimTrait
#' @export
#'
#' @examples \dontrun{pheno.generate(genotypes = t(geno_data[-c(1:3)]), effects = b,
#'                    envNum = envNum, indNum = indNum, sigma.error = sigma_error)}
pheno.generate <- function(genotypes, effects, envNum, indNum, sigma.error){

  if(nrow(genotypes)!=indNum)
  {
    stop("Error: pheno.generate() function, the number of rows of genotypes doesn't match with indNum. Make sure it is an N*M matrix.")
  }

  if(nrow(effects)!=envNum)
  {
    stop("Error: pheno.generate() function, the number of rows of effects doesn't match with envNum. Make sure it is an H*M matrix.")
  }

  pheno_data = data.frame(ENV = as.factor(rep(paste0("env",1:envNum), each=indNum)),
                          GID = as.factor(rep(rownames(genotypes), envNum)))

  g = tcrossprod(genotypes, b)
  g = as.vector(g)
  error = rnorm(indNum*envNum, mean = 0, sd = sqrt(sigma.error))
  pheno_data = cbind(pheno_data, g + error)

  colnames(pheno_data) = c("ENV","GID","SimTrait")

  pheno_data$ENV = as.factor(pheno_data$ENV)
  pheno_data$GID = as.factor(pheno_data$GID)

  return(pheno_data)
}
