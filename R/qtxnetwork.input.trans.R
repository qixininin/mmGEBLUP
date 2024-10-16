#' qtxnetwork.input.trans function
#'
#' @param geno_data genotype data, output from geno.generate()
#' @param pheno_data phentype data, output from pheno.generate()
#' @param geno_output_prefix output prefix for .gen file
#' @param pheno_output_prefix output prefix for .phe file
#'
#' @return NULL
#' @export
#' @importFrom utils write.table
#'
#' @examples \dontrun{qtxnetwork.input.trans(geno_data, pheno_data,
#'                                           geno_output_prefix = geno_output_prefix,
#'                                           pheno_output_prefix = pheno_output_prefix)}
qtxnetwork.input.trans <- function(geno_data, pheno_data, geno_output_prefix, pheno_output_prefix)
{

  # chromosome or linkage group summary
  chr = c(table(geno_data[,1]))
  chrName = names(chr)
  # rownames of geno_data
  rownames(geno_data) = geno_data[,2]
  geno_data = t(geno_data[,-c(1:3)])
  gid = rownames(geno_data)
  geno_data = cbind(gid, geno_data)

  ## phenotype
  # pheno_data = read.table(pheno_file, header = T, na.strings = ".")
  traitName = colnames(pheno_data)[-c(1:2)]
  traitNum = length(traitName)

  ## output ----------------------------------------------------------------------
  colnames(geno_data)[1] = "#Ind"
  colnames(pheno_data)[1:2] = c("Env#", "Geno#")

  # .gen
  print(paste0("Writing QTXNetwork2.0 \'.gen\' files at ", pheno_output_prefix, ".gen"))
  outputFile = paste0(geno_output_prefix,".gen")
  printer = file(outputFile, "w")
  write.table(geno_data, file = printer, append = T, na = ".",
              row.names = F, col.names = T, quote = F,sep = "\t", eol = "\n")
  close(printer)

  # .phe
  print(paste0("Writing ", traitNum, " traits into QTXNetwork2.0 \'.phe\' files at ", pheno_output_prefix, "_TRAIT.phe"))
  for(c in 1:traitNum)
  {
    trait = traitName[c]
    outputFile = paste0(pheno_output_prefix,"_", trait, ".phe")
    printer = file(outputFile, "w")
    write("_Population\tFFFFFFF", printer)
    write(paste0("_Genotypes\t", nrow(geno_data)), printer, append = T)
    write(paste0("_Observations\t", nrow(pheno_data)), printer, append = T)
    write("_Environments\tyes", printer, append = T)
    write("_Replications\tno", printer, append = T)
    write(paste0("_TraitNumber\t", 1), printer, append = T)
    write(paste0("_Chromosomes\t", length(chrName), "\t", paste(chrName, collapse = " ")), printer, append = T)
    write(paste0("_TotalMarker\t", sum(chr), "\t", paste(chr, collapse = " ")), printer, append = T)
    write("_MarkerCode\tP1=1\tP2=-1\tF1=0\n", printer, append = T)
    write("*TraitBegin*", printer, append = T)
    write.table(pheno_data[,c(1,2,c+2)], file = printer, append = T, na = ".",
                row.names = F, col.names = T, quote = F,sep = "\t", eol = ";\n")
    write("*TraitEnd*", printer, append = T)
    close(printer)
  }

}
