#' mmgeblup function
#'
#' @param data a data frame
#' @param Ka additive genetic relationship matrix
#' @param KaeE1 additive-by-environment relationship matrix1
#' @param KaeE2 additive-by-environment relationship matrix2
#'
#' @return list(mod, BV)
#' @export
#' @import sommer
#' @import dplyr
#'
#' @examples \dontrun{rst = mmgeblup(data = cbind(dt, mmdata$Xa), Ka = mmdata$Ka, KaeE1 = mmdata$KaeE1, KaeE2 = mmdata$KaeE2)}
mmgeblup <- function(data, Ka, KaeE1, KaeE2)
{
  ## Evaluate input
  if(missing(Ka) | missing(KaeE1)){
    stop("Error: no input for additive kinship matrix or additive-by-environment kinship matrix.")
  }
  majorGE = TRUE
  if(missing(KaeE2)){
    majorGE = FALSE
  }

  ## Receive fixed effect column names
  ## If no fixed effect columns in data, then only intercept is fixed effect
  fixColName = names(data)[!(names(data) %in% c("ENV","GID","trait"))]
  if(length(fixColName)==0){
    fixColName = "1"
  }

  ## Perform sommer models
  if(majorGE){ # If two AE kinship matrices are inputted

    m_env_qtl = length(KaeE1)
    datafake = data
    for(t in 1:m_env_qtl)
    {
      datafake = datafake %>% dplyr::mutate(!! paste0("GID",t) :=GID)
    }

    mod = mmer(reformulate(fixColName, "trait"),
               random = as.formula(paste0("~vsr(GID, Gu=Ka) + vsr(ENV) + ",
                                          paste0("vsr(GID", 1:m_env_qtl,":ENV, Gu=KaeE1[[", 1:m_env_qtl,"]])", collapse = " + "),
                                          " + vsr(GID:ENV, Gu=KaeE2)")),
               rcov = ~units,
               data = datafake,
               verbose = FALSE, date.warning = FALSE)

    BV = data.frame(data[,c("ENV","GID")])
    BV$mu = mod$Beta$Estimate[1]                                                                     # mu
    BV$A_l = as.matrix(data[,as.vector(mod$Beta$Effect)[-1]]) %*% as.vector(mod$Beta$Estimate)[-1]   # G-major
    BV$A_s = mod$U$`u:GID`$trait[BV$GID]                                                             # G-minor
    BV$AE_l = 0
    for(t in 1:m_env_qtl)                                                                            # GE-major
    {
      BV$AE_l = BV$AE_l + mod$U[[2+t]]$trait[paste0(BV$GID,":",BV$ENV)]
    }
    BV$AE_s = mod$U$`u:GID:ENV`$trait[paste0(BV$GID,":",BV$ENV)]                                     # GE-minor
    BV$E = mod$U$`u:ENV`$trait[BV$ENV]                                                               # E
    BV$pre = BV$mu + BV$A_l + BV$A_s + BV$AE_l + BV$AE_s + BV$E

  } else { # If only one AE kinship matrices is inputted
    mod = mmer(reformulate(fixColName, "trait"),
               random = ~vsr(GID, Gu=Ka) + vsr(ENV) + vsr(GID:ENV, Gu=KaeE1),
               rcov = ~units,
               data = data,
               verbose = FALSE, date.warning = FALSE)
    ## Predict
    # mu = as.matrix(cbind(rep(1, nrow(data)), data[,as.vector(mod$Beta$Effect)[-1]])) %*% mod$Beta$Estimate
    BV = data.frame(data[,c("ENV","GID")])
    BV$mu = mod$Beta$Estimate[1]                                                                     # mu
    BV$A_l = as.matrix(data[,as.vector(mod$Beta$Effect)[-1]]) %*% as.vector(mod$Beta$Estimate)[-1]   # G-major
    BV$A_s = mod$U$`u:GID`$trait[BV$GID]                                                             # G-minor
    BV$AE = mod$U$`u:GID:ENV`$trait[paste0(BV$GID,":",BV$ENV)]                                       # GE
    BV$E = mod$U$`u:ENV`$trait[BV$ENV]                                                               # E
    BV$pre = BV$mu + BV$A_l + BV$A_s + BV$AE + BV$E
  }

  return(list(mod, BV))

}
