#' geno.generate function
#' A simple function to generate homozygotic genotype data (coded with 1 or -1)
#' with respect to the number of individual, the number of markers, MAF range and chromosome.
#' Markers are independent.
#'
#' @param indNum The number of individual
#' @param snpNum The number of markers
#' @param maf.min The lower bound of MAF, so that MAF is sampled from U(maf.min, maf.max)
#' @param maf.max The upper bound of MAF, so that MAF is sampled from U(maf.min, maf.max)
#' @param chr.snpNum The number of snps on each chromosome
#'
#' @return geno_data a data frame, where the first three columns are $CHR, $SNP, and $BP
#' @export
#' @importFrom stats rbinom runif
#'
#' @examples geno.generate(100, 200, 0.05, 0.5, c(50,50,40,60))
geno.generate <- function(indNum, snpNum, maf.min, maf.max, chr.snpNum)
{

  if(sum(chr.snpNum)!=snpNum){
    stop("Error: geno.generate() function: the summation of chr.snpNum is not equal to snpNum.")
  }

  Ga = matrix(NA, indNum, snpNum)
  freq = runif(snpNum, maf.min, maf.max)
  for(j in 1:snpNum)
  {
    # Ga[,j] = rbinom(indNum, 2, freq[j])-1              # additive(-1,0,1)
    Ga[,j] = ifelse(rbinom(indNum, 1, freq[j])==0,-1,1)  # additive(-1,1)
  }
  rownames(Ga) = paste0("gid", 1:indNum)
  colnames(Ga) = paste0("snp", 1:snpNum)

  chrName = paste0("chr", 1:length(chr.snpNum))


  geno_data = data.frame(CHR = as.vector(unlist(sapply(chrName, function(x) rep(x, chr.snpNum[which(chrName == x)])))),
                         SNP = colnames(Ga),
                         BP = as.vector(unlist(sapply(chr.snpNum, function(x) seq(1:x)))),
                         t(Ga))

  return(geno_data)
}

