#' mmgblup function
#'
#' @param data a data frame
#' @param Ka additive genetic relationship matrix
#'
#' @return list(mod, BV)
#' @export
#' @import sommer
#' @import dplyr
#'
#' @examples \dontrun{rst = mmgblup(data = cbind(dt, mmdata$Xa), Ka = mmdata$Ka)}
mmgblup <- function(data, Ka)
{
  ## Evaluate input
  if(missing(Ka)){
    stop("Error: no input for additive kinship matrix")
  }

  ## Receive fixed effect column names
  ## If no fixed effect columns in data, then only intercept is fixed effect
  fixColName = names(data)[!(names(data) %in% c("ENV","GID","trait"))]
  if(length(fixColName)==0){
    fixColName = "1"
  }

  mod = mmer(reformulate(fixColName, "trait"),
             random = ~vsr(GID, Gu=Ka) + vsr(ENV) ,
             rcov = ~units,
             data = data,
             verbose = FALSE, date.warning = FALSE)

  BV = data.frame(data[,c("ENV","GID")])
  BV$mu  = mod$Beta$Estimate[1]                                                                     # mu
  BV$A_l = as.matrix(data[,as.vector(mod$Beta$Effect)[-1]]) %*% as.vector(mod$Beta$Estimate)[-1]    # G-major
  BV$A_s = mod$U$`u:GID`$trait[BV$GID]                                                              # G-minor
  BV$E   = mod$U$`u:ENV`$trait[BV$ENV]                                                              # E
  BV$pre = BV$mu + BV$A_l + BV$A_s +  BV$E

  return(list(mod, BV))

}
