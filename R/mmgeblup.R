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

  ## Perform sommer models
  if(majorGE){ # If two AE kinship matrices are inputted
    mod = mmer(reformulate(fixColName, "trait"),
               random = ~vsr(GID, Gu=Ka) + vsr(ENV) + vsr(ENV:GID, Gu=EKae1) + vsr(ENV:GID1, Gu=EKae2),
               rcov = ~units,
               data = datafake,
               verbose = FALSE, date.warning = FALSE)
    ## Predict
    mu = as.matrix(cbind(rep(1, nrow(data)), data[,as.vector(mod$Beta$Effect)[-1]])) %*% mod$Beta$Estimate
    BV = data.frame(data[,c("ENV","GID")])
    BV$pre = mu +                                          # mu
      mod$U$`u:GID`$trait[BV$GID]+                         # G
      mod$U$`u:ENV`$trait[BV$ENV]+                         # E
      mod$U$`u:ENV:GID`$trait[paste0(BV$ENV,":",BV$GID)]+  # GE-major
      mod$U$`u:ENV:GID1`$trait[paste0(BV$ENV,":",BV$GID)]  # GE-minor

  } else { # If only one AE kinship matrices is inputted
    mod = mmer(reformulate(fixColName, "trait"),
               random = ~vsr(GID, Gu=Ka) + vsr(ENV) + vsr(ENV:GID, Gu=EKae1),
               rcov = ~units,
               data = data,
               verbose = FALSE, date.warning = FALSE)
    ## Predict
    mu = as.matrix(cbind(rep(1, nrow(data)), data[,as.vector(mod$Beta$Effect)[-1]])) %*% mod$Beta$Estimate
    BV = data.frame(data[,c("ENV","GID")])
    BV$pre = mu +                                          # mu
      mod$U$`u:GID`$trait[BV$GID]+                         # G
      mod$U$`u:ENV`$trait[BV$ENV]+                         # E
      mod$U$`u:ENV:GID`$trait[paste0(BV$ENV,":",BV$GID)]   # GE
  }

  return(list(mod, BV))

}
