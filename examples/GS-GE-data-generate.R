#### Install package -----------------------------------------------------------
## Install from github
library(devtools)
install_github("qixininin/deepKin")
## Install from source files
install.packages("/Users/zqx/ZJU-PhD/2-Work/Rpackage/mmGEBLUP_0.1.0.tar.gz", repos = NULL, type = "source")

#### Load package --------------------------------------------------------------
library(mmGEBLUP)

#### Data simulation -----------------------------------------------------------
# Parameter setting
# repNum = 1
envNum = 3
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

# Effect generation
eff_list = snp.effect(snpNum = snpNum, envNum = envNum,
                      major_a_idx = major_a_idx, major_ae_idx = major_ae_idx,
                      variance_a_major = sigma_a_major , variance_a_minor = sigma_a_minor,
                      variance_ae_major = sigma_ae_major, variance_ae_minor = sigma_ae_minor)
b = eff_list$effects
# Plot phenotypes
p = snp.effect.plot(effects = b)

# Genotype generation
geno_data = geno.generate(indNum = indNum, snpNum = snpNum, maf.min = 0.05, maf.max = 0.5,
                          chr.snpNum = c(500, 500, 500, 500))
# Phenotype generation
pheno_data = pheno.generate(genotypes = t(geno_data[-c(1:3)]), effects = b,
                            envNum = envNum, indNum = indNum, sigma.error = sigma_error)


save(geno_data, pheno_data, file = "./inst/Simulation.Rdata")

## QTXnetwork ------------------------------------------------------------------
load("./inst/Simulation.Rdata")
pheno_output_prefix = "./inst/Simulation"
geno_output_prefix  = "./inst/Simulation"

qtxnetwork.input.trans(geno_data, pheno_data,
                       geno_output_prefix = geno_output_prefix,
                       pheno_output_prefix = pheno_output_prefix)


# Perform QTSNetwork -----------------------------------------------------------
qtxnetwork_path = "/public3/zqx/QTXNetwork_4.0/build/QTXNetwork"

commands <- lapply(phe_files, function(phe_file) {
  paste0(qtxnetwork_path,
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
