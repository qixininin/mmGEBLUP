install.packages("/Users/zqx/ZJU-PhD/2-Work/Rpackage/mmGEBLUP_0.1.0.tar.gz", repos = NULL, type = "source")

library(mmGEBLUP)

# Parameter setting ------------------------------
# repNum = 1
envNum = 3
pheNum = 5
indNum = 1000
snpNum = 2000
rho_a = 0.5
rho_ae = 0.5
h2_a = 0.3
h2_ae = 0.3
sigma_a = h2_a
sigma_ae = h2_ae
sigma_error = 1-h2_a-h2_ae

major_a_idx = c(500, 750, 1000, 1250, 1500)
snpNum_a_major = length(major_a_idx)
snpNum_a_minor = snpNum-snpNum_a_major

major_ae_idx = c(250, 500, 1000, 1500, 1750)
snpNum_ae_major = length(major_ae_idx)
snpNum_ae_minor = snpNum-snpNum_ae_major

sigma_a_major = rho_a * sigma_a/ snpNum_a_major
sigma_a_minor = (1-rho_a) * sigma_a / snpNum_a_minor
sigma_ae_major = rho_ae * sigma_ae/ snpNum_ae_major
sigma_ae_minor = (1-rho_ae) * sigma_ae / snpNum_ae_minor

# Effect generation --------------------------------------------------------
eff_list = snp.effect(snpNum = snpNum, envNum = envNum,
                      major_a_idx = major_a_idx, major_ae_idx = major_ae_idx,
                      variance_a_major = sigma_a_major , variance_a_minor = sigma_a_minor,
                      variance_ae_major = sigma_ae_major, variance_ae_minor = sigma_ae_minor)
b0 = eff_list$main_effects
bh = eff_list$interact_effects
b = eff_list$effects

p = snp.effect.plot(effects = b)
p

# Genotype generation --------------------------------------------------------
Ga = matrix(NA, indNum, snpNum)
freq = runif(snpNum, 0.05, 0.5)
for(j in 1:snpNum)
{
  # Ga[,j] = rbinom(indNum, 2, freq[j])-1              # additive(-1,0,1)
  Ga[,j] = ifelse(rbinom(indNum, 1, freq[j])==0,-1,1)  # additive(-1,1)
}
rownames(Ga) = paste0("gid", 1:indNum)
colnames(Ga) = paste0("snp", 1:snpNum)

geno_data = data.frame(CHR = rep(c("chr1","chr2","chr3", "chr4"), each = snpNum/4),
                       SNP = colnames(Ga),
                       BP = rep(seq(1:(snpNum/4)),4),
                       t(Ga))
# write.csv(geno_data, file = "Simulation_geno.csv", quote = F, col.names = T, row.names = F)

# Phenotype generation ---------------------------
qtl = data.frame()
pheno_data = data.frame(ENV = as.factor(rep(paste0("env",1:envNum), each=indNum)),
                        GID = as.factor(rep(rownames(Ga), envNum)))
for(p in 1:pheNum)
{
  g = tcrossprod(Ga, b)
  g = as.vector(g)
  error = rnorm(indNum*envNum, mean = 0, sd = sqrt(sigma_error))
  pheno_data = cbind(pheno_data, g + error)
  qtl = rbind(qtl, data.frame(TRAIT = paste0("TRAIT",p),
                              QTL = paste0("snp:",major_idx)))
  print(var(g)/(var(error)+var(g)))
  # Examine rhoA
  # var_l = diag(var(tcrossprod(Ga[,major_idx], t(as.matrix(b0[major_idx])))))
  # var_s = diag(var(tcrossprod(Ga[,-major_idx], t(as.matrix(b0[-major_idx])))))
  # print(var_l/(var_l+var_s))
}


colnames(pheno_data) = c("ENV","GID",paste0("TRAIT",1:pheNum))

write.table(pheno_data, file = paste0(filePrefix, ".txt"), quote = F, col.names = T, row.names = F)


# QTSformat-input.R ------------------------------------------------------------
# This R script will help you transform into QTS .Gen and .Phe files
# Please be careful about the build-in header in .Phe files
setwd("/public3/zqx/GS-GE/2-QTS/")

# NOTES:
# 1. pheno_file should have at least three columns, with a header. The first two columns are "ENV" and "GID"
# 2. geno_file should be a csv file, with a header. The first three columns are "CHR", "SNP", and "BP.
geno_file = "/public3/zqx/GS-GE/1-Data/Simulation_geno.csv"
pheno_file = paste0("/public3/zqx/GS-GE/1-Data/",filePrefix,".txt")
geno_output_prefix = "Simulation"
pheno_output_prefix = filePrefix

## genotype --------------------------------------------------------------------
geno_data = read.csv(geno_file)
# chromosome or linkage group summary
chr = c(table(geno_data[,1]))
chrName = names(chr)
# rownames of geno_data
rownames(geno_data) = geno_data[,2]
geno_data = t(geno_data[,-c(1:3)])
gid = rownames(geno_data)
geno_data = cbind(gid, geno_data)

## phenotype -------------------------------------------------------------------
pheno_data = read.table(pheno_file, header = T, na.strings = ".")
traitName = colnames(pheno_data)[-c(1:2)]
traitNum = length(traitName)

## output ----------------------------------------------------------------------
colnames(geno_data)[1] = "#Ind"
colnames(pheno_data)[1:2] = c("Env#", "Geno#")

# .phe
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
  write.table(pheno_data[,c(1,2,c+2)], file = printer, append = T,
              row.names = F, col.names = T, quote = F,sep = "\t", eol = ";\n")
  write("*TraitEnd*", printer, append = T)
  close(printer)
}

# .gen
# outputFile = paste0(geno_output_prefix,".gen")
# printer = file(outputFile, "w")
# write.table(geno_data, file = printer, append = T,
#             row.names = F, col.names = T, quote = F,sep = "\t", eol = "\n")
# close(printer)

# Perform QTSNetwork -----------------------------------------------------------
data_dir <- "/public3/zqx/GS-GE/2-QTS"
phe_files <- list.files(data_dir, pattern = filePrefix, full.names = TRUE)

commands <- lapply(phe_files, function(phe_file) {
  paste0("/public3/zqx/QTXNetwork_4.0/build/QTXNetwork",
         " --map /public3/zqx/GS-GE/2-QTS/Simulation.gen",
         " --txt ", phe_file,
         " --out ", sub("\\.phe$", ".pre", phe_file),
         " --QTX 1 --only-1D 1")
})
results <- sapply(commands, system, wait = TRUE)
if (all(results == 0)) {
  cat("All commands executed successfully.\n")
} else {
  failed_commands <- commands[results != 0]
  cat("Some commands failed to execute:\n")
  cat(paste(failed_commands, collapse = "\n"))
}

# QTSformat-output.R -----------------------------------------------------------
geno_file = "/public3/zqx/GS-GE/1-Data/Simulation_geno.csv"
pheno_file = paste0("/public3/zqx/GS-GE/1-Data/",filePrefix, ".txt")
geno_data = read.csv(geno_file)
pheno_data = read.table(pheno_file, header = T, na.strings = ".")
traitName = colnames(pheno_data)[-c(1:2)]
traitNum = length(traitName)
envName = unique(pheno_data$ENV)
envNum = length(envName)

start_title = "_1D_effect"
end_title = "_1D_heritability"
df.qtl = data.frame()
for(c in 1:traitNum)
{
  trait = traitName[c]
  pre_file = list.files(".", pattern = paste0(filePrefix,"_",trait, ".*\\.pre$"), full.names = FALSE)
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

# GS prepare -------------------------------------------------------------------
# Include marker_data and geno_data
# prepare genotype matrix(rows: individual id/columns: marker id)
marker_data = geno_data[,1:3]
colnames(marker_data) = c("CHR","SNP","Distance")
rownames(geno_data) = geno_data[,2]
geno_data = t(geno_data[,-c(1:3)])

# Phenotype
pheno_data$ENV = as.factor(pheno_data$ENV)
pheno_data$GID = as.factor(pheno_data$GID)

save(geno_data, pheno_data, marker_data, qtl_data, qtl_env_data, file = paste0("/public3/zqx/GS-GE/3-GS/",filePrefix, "_forGS.Rdata"))

# Perform GS -------------------------------------------------------------------
setwd("/public3/zqx/GS-GE/3-GS/")

system(paste0("nohup /usr/local/R/bin/Rscript GS-GE-beta-2.1.R ",filePrefix," 1 > ",filePrefix,"-rep1.log 2>&1 &"))
system(paste0("nohup /usr/local/R/bin/Rscript GS-GE-beta-2.1.R ",filePrefix," 2 > ",filePrefix,"-rep2.log 2>&1 &"))
system(paste0("nohup /usr/local/R/bin/Rscript GS-GE-beta-2.1.R ",filePrefix," 3 > ",filePrefix,"-rep3.log 2>&1 &"))
system(paste0("nohup /usr/local/R/bin/Rscript GS-GE-beta-2.1.R ",filePrefix," 4 > ",filePrefix,"-rep4.log 2>&1 &"))
system(paste0("nohup /usr/local/R/bin/Rscript GS-GE-beta-2.1.R ",filePrefix," 5 > ",filePrefix,"-rep5.log 2>&1 &"))
system(paste0("nohup /usr/local/R/bin/Rscript GS-GE-beta-2.1.R ",filePrefix," 6 > ",filePrefix,"-rep6.log 2>&1 &"))
system(paste0("nohup /usr/local/R/bin/Rscript GS-GE-beta-2.1.R ",filePrefix," 7 > ",filePrefix,"-rep7.log 2>&1 &"))
system(paste0("nohup /usr/local/R/bin/Rscript GS-GE-beta-2.1.R ",filePrefix," 8 > ",filePrefix,"-rep8.log 2>&1 &"))
system(paste0("nohup /usr/local/R/bin/Rscript GS-GE-beta-2.1.R ",filePrefix," 9 > ",filePrefix,"-rep9.log 2>&1 &"))
system(paste0("nohup /usr/local/R/bin/Rscript GS-GE-beta-2.1.R ",filePrefix," 10 > ",filePrefix,"-rep10.log 2>&1 &"))
