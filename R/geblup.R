#' geblup function
#'
#' @param data a data frame
#' @param A additive genetic relationship matrix
#' @param AE additive-by-environment relationship matrix1
#'
#' @return list(mod, BV)
#' @export
#' @import sommer
#' @import dplyr
#'
#' @examples \dontrun{rst = geblup(data = dt, A = mmdata$A, AE = mmdata$AE)}
geblup <- function(data, A, AE)
{
  ## Evaluate input
  if(missing(A) | missing(AE)){
    stop("Error: no input for additive kinship matrix or additive-by-environment covariance matrix")
  }

  mod = mmer(reformulate("1", "trait"),
             random = ~vsr(GID, Gu=A) + vsr(ENV) + vsr(GID:ENV, Gu=AE),
             rcov = ~units,
             data = data,
             verbose = FALSE, date.warning = FALSE)
  ## Predict
  BV = data.frame(data[,c("ENV","GID")])
  BV$mu = mod$Beta$Estimate[1]                                                                     # mu
  BV$A  = mod$U$`u:GID`$trait[BV$GID]                                                              # G
  BV$AE = mod$U$`u:GID:ENV`$trait[paste0(BV$GID,":",BV$ENV)]                                       # GE
  BV$E  = mod$U$`u:ENV`$trait[BV$ENV]                                                              # E
  BV$pre = BV$mu + BV$A + BV$AE + BV$E


  return(list(mod, BV))

}
