#' qtxnetwork.output.trans
#'
#' @param pheno_data a data frame of phenotypic data, one of the outputs from simulation
#' @param pre_file path and name for .pre file
#'
#' @return a list
#'         $ qtl_data a data frame with additive qtl information
#'         $ qtl_env_data a data frame with additive-by-environment qtl information
#' @export
#'
#' @examples \dontrun{qtl = qtxnetwork.output.trans(pheno_data, pre_file)
#'                    qtl$qtl_data
#'                    qtl$qtl_env_data}
qtxnetwork.output.trans <- function(pheno_data, pre_file)
{
  traitName = colnames(pheno_data)[-c(1:2)]
  traitNum = length(traitName)
  envName = as.character(unique(pheno_data$ENV))
  envNum = length(envName)

  start_title = "_1D_effect"
  end_title = "_1D_heritability"
  df.qtl = data.frame()
  for(c in 1:traitNum)
  {
    trait = traitName[c]
    file_lines = readLines(pre_file)
    # locate target lines
    start_line <- 0
    end_line <- 0
    for (i in 1:length(file_lines)) {
      line <- file_lines[i]
      if (line == start_title) {
        start_line <- i
        next
      }
      if (line == end_title) {
        end_line <- i
        break
      }
    }
    # store in data.frame
    if(start_line) {
      tmp <- as.data.frame(do.call(rbind, strsplit(file_lines[(start_line + 2):(end_line - 2)], "\\s+")))
      df.qtl = rbind(df.qtl, data.frame(trait, tmp))
    }
  }

  if(nrow(df.qtl)==0) {return(list(qtl_data = data.frame(TRAIT=character(), QTL=character()),
                                   qtl_env_data = data.frame(TRAIT=character(), ENV=character(), QTL=character())))}

  colnames(df.qtl) = c("TRAIT", "QTL", "SNPID", "A", "SE", "P-Value",
                       as.vector(rbind(paste0(envName),paste0("SE", 1:envNum),paste0("Pvalue", 1:envNum))))

  dt = df.qtl %>% dplyr::select(c("TRAIT","SNPID","A", envName))
  dt_reshape=reshape(dt,
                     idvar=c("TRAIT", "SNPID"),
                     varying=c("A",envName),
                     v.names="Effect",
                     timevar="ENV",
                     times=c("A",envName),
                     direction="long")
  # Additive qtl
  qtl_data = dt_reshape %>%
    dplyr::filter(ENV == "A") %>%
    dplyr::select(c("TRAIT", "SNPID")) %>%
    dplyr::rename(QTL = SNPID)
  # AE qtl
  qtl_env_data = dt_reshape %>%
    dplyr::filter(ENV != "A") %>%
    dplyr::filter(Effect != "---") %>%
    dplyr::select(c("TRAIT", "ENV", "SNPID")) %>%
    dplyr::rename(QTL = SNPID) %>%
    dplyr::arrange(factor(TRAIT, levels = traitName))

  if(nrow(qtl_data)>0) rownames(qtl_data) = 1:nrow(qtl_data)
  if(nrow(qtl_env_data)>0) rownames(qtl_env_data) = 1:nrow(qtl_env_data)

  return(list(qtl_data = qtl_data, qtl_env_data = qtl_env_data))
}
