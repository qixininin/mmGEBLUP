#' gblup function
#'
#' @param data a data frame
#' @param A additive genetic relationship matrix
#'
#' @return list(mod, BV)
#' @export
#' @import sommer
#' @import dplyr
#'
#' @examples \dontrun{rst = gblup(data = dt, A = mmdata$A)}
gblup <- function(data, A)
{
  ## Evaluate input
  if(missing(A)){
    stop("Error: no input for additive kinship matrix")
  }

  mod = mmer(reformulate("1", "trait"),
             random = ~vsr(GID, Gu=A) + vsr(ENV),
             rcov = ~units,
             data = data,
             verbose = FALSE, date.warning = FALSE)
  ## Predict
  BV = data.frame(data[,c("ENV","GID")])
  BV$mu = mod$Beta$Estimate[1]                                                                     # mu
  BV$A  = mod$U$`u:GID`$trait[BV$GID]                                                              # G
  BV$E  = mod$U$`u:ENV`$trait[BV$ENV]                                                              # E
  BV$pre = BV$mu + BV$A + BV$E

  return(list(mod, BV))
}
