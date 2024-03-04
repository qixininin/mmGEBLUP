dataqc <- function(geno_data, pheno_data)
{
  # Summary
  traitName = colnames(pheno_data)[-c(1:2)]
  traitNum = length(traitName)
  trialName = unique(pheno_data$ENV)
  trialNum = length(trialName)
  lineName = intersect(unique(pheno_data$GID),rownames(geno_data))
  lineNum = length(lineName)

  # Prepare genotype
  rownames(geno_data) = geno_data[,2]
  geno_data = t(geno_data[,-c(1:3)])
  geno_data = as.matrix(geno_data[lineName, ])

  # Prepare phenotype
  pheno_data = pheno_data[which(pheno_data$GID %in% lineName), ]

  return(list(geno_data_qc = geno_data,
              pheno_data_qc = pheno_data,
              summary = list()))


}
