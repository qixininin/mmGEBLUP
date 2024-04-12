#' snp.effect.plot function
#'
#' @param effects an envNum * snpNum matrix for total SNP effect
#'
#' @return a ggplot object
#' @export
#' @import ggplot2
#'
#' @examples \dontrun{p = snp.effect.plot(effects = b)}
snp.effect.plot <- function(effects)
{
  snpNum = ncol(effects)
  envNum = nrow(effects)
  dt = data.frame(x=rep(1:snpNum,envNum),
                  y=as.vector(t(effects)),
                  env=rep(paste0("Env",1:envNum),each=snpNum))
  lmt = max(max(dt$y), -min(dt$y)) * 1.1
  p = ggplot(dt,aes(x=x,y=y))+
    geom_point(cex = 0.8) +
    scale_y_continuous(limits = c(-lmt,lmt))+
    facet_grid(~env) +
    theme_bw() +
    labs(x="SNP",y="Overall effect")
  return(p)
}
