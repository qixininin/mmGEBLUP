#' mmgeblup function
#'
#' @param data a data frame
#' @param Ka additive genetic relationship matrix
#' @param EKae1 additive-by-environment relationship matrix1
#' @param EKae2 additive-by-environment relationship matrix2
#'
#' @return list(mod, BV)
#' @export
#' @import sommer
#' @import dplyr
#'
#' @examples \dontrun{rst = mmgeblup(data = cbind(dt, mmdata$Xa), Ka = mmdata$Ka, EKae1 = mmdata$EKae1, EKae2 = mmdata$EKae2)}
mmgeblup <- function(data, Ka, EKae1, EKae2)
{
  ## Evaluate input
  if(missing(Ka) | missing(EKae1)){
    stop("Error: no input for additive kinship matrix or additive-by-environment kinship matrix.")
  }
  majorGE = TRUE
  if(missing(EKae2)){
    majorGE = FALSE
  }

  ## Receive fixed effect column names
  ## If no fixed effect columns in data, then only intercept is fixed effect
  fixColName = names(data)[!(names(data) %in% c("ENV","GID","trait"))]
  if(length(fixColName)==0){
    fixColName = "1"
  }

  ## Data prepare
  datafake = data %>% dplyr::mutate(GID1 = GID)

  ## TEST --------

  site_qe = unique(qtl_env_data$QTL)
  Xae = mmgeno_data[, colnames(mmgeno_data) %in% site_qe]
  EKae = apply(Xae, 2, simplify = F,  function(x)  {
    A = tcrossprod(x)
    colnames(A) = rownames(A) = names(x)
    kronecker(A, E, make.dimnames = T) } )

  datafake = data %>%
    dplyr::mutate(GID1 = GID) %>%
    dplyr::mutate(GID2 = GID) %>%
    dplyr::mutate(GID3 = GID) %>%
    dplyr::mutate(GID4 = GID)

  mod = mmer(reformulate(fixColName, "trait"),
             random = ~vsr(GID, Gu=Ka) + vsr(ENV) +
               vsr(GID1:ENV, Gu=EKae[[1]]) + vsr(GID1:ENV, Gu=EKae[[2]]) + vsr(GID1:ENV, Gu=EKae[[3]]) + vsr(GID1:ENV, Gu=EKae[[4]]) +
               vsr(ENV:GID, Gu=EKae2),
             rcov = ~units,
             data = datafake,
             verbose = FALSE, date.warning = FALSE)
  ## END TEST ---------

  ## Perform sommer models
  if(majorGE){ # If two AE kinship matrices are inputted
    mod = mmer(reformulate(fixColName, "trait"),
               random = ~vsr(GID, Gu=Ka) + vsr(ENV) + vsr(ENV:GID, Gu=EKae1) + vsr(ENV:GID1, Gu=EKae2),
               rcov = ~units,
               data = datafake,
               verbose = FALSE, date.warning = FALSE)
    ## Predict
    # mu = as.matrix(cbind(rep(1, nrow(data)), data[,as.vector(mod$Beta$Effect)[-1]])) %*% mod$Beta$Estimate
    BV = data.frame(data[,c("ENV","GID")])
    BV$mu = mod$Beta$Estimate[1]                                                                     # mu
    BV$A_l = as.matrix(data[,as.vector(mod$Beta$Effect)[-1]]) %*% as.vector(mod$Beta$Estimate)[-1]   # G-major
    BV$A_s = mod$U$`u:GID`$trait[BV$GID]                                                             # G-minor
    BV$AE_l = mod$U$`u:ENV:GID`$trait[paste0(BV$ENV,":",BV$GID)]                                     # GE-major
    BV$AE_s = mod$U$`u:ENV:GID1`$trait[paste0(BV$ENV,":",BV$GID)]                                    # GE-minor
    BV$E = mod$U$`u:ENV`$trait[BV$ENV]                                                               # E
    BV$pre = BV$mu + BV$A_l + BV$A_s + BV$AE_l + BV$AE_s + BV$E

  } else { # If only one AE kinship matrices is inputted
    mod = mmer(reformulate(fixColName, "trait"),
               random = ~vsr(GID, Gu=Ka) + vsr(ENV) + vsr(ENV:GID, Gu=EKae1),
               rcov = ~units,
               data = data,
               verbose = FALSE, date.warning = FALSE)
    ## Predict
    # mu = as.matrix(cbind(rep(1, nrow(data)), data[,as.vector(mod$Beta$Effect)[-1]])) %*% mod$Beta$Estimate
    BV = data.frame(data[,c("ENV","GID")])
    BV$mu = mod$Beta$Estimate[1]                                                                     # mu
    BV$A_l = as.matrix(data[,as.vector(mod$Beta$Effect)[-1]]) %*% as.vector(mod$Beta$Estimate)[-1]   # G-major
    BV$A_s = mod$U$`u:GID`$trait[BV$GID]                                                             # G-minor
    BV$AE = mod$U$`u:ENV:GID`$trait[paste0(BV$ENV,":",BV$GID)]                                       # GE
    BV$E = mod$U$`u:ENV`$trait[BV$ENV]                                                               # E
    BV$pre = BV$mu + BV$A_l + BV$A_s + BV$AE + BV$E
  }

  return(list(mod, BV))

}
