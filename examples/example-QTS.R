## This example file is to perform QTS analysis based on genotype and phenotype data
## 1. If you have already obtained QTS results, you can skip this example.
## 2. If you have generated the simulation data by following the 'example-data-generate.R',
##    you can load 'Simulation.Rdata' directly.
## 3. If you want to input your own data for 'QTXnetwork' from R,
##    you can also perpare R-formated 'geno_data' and 'pheno_data' into 'Prefix.Rdata'

# Load packages ----------------------------------------------------------------
library(mmGEBLUP)
library(dplyr)

# Prepare QTXNetwork ------------------------------------------------------------------
load("./inst/Simulation-genphe.Rdata")
pheno_output_prefix = "./inst/Simulation"
geno_output_prefix  = "./inst/Simulation"

## The default output names are
## .gen = Simulation.gen
## .phe = Simulation_SimTrait.phe
qtxnetwork.input.trans(geno_data, pheno_data,
                       geno_output_prefix = geno_output_prefix,
                       pheno_output_prefix = pheno_output_prefix)


# Perform QTXNetwork -----------------------------------------------------------
qtxnetwork_path = "/public3/zqx/QTXNetwork_4.0/build/QTXNetwork"
gen_file = "./inst/Simulation.gen"
phe_file = "./inst/Simulation_SimTrait.phe"

system(paste0(qtxnetwork_path,
              " --map ", gen_file,
              " --txt ", phe_file,
              " --out ", sub("\\.phe$", ".pre", phe_file),
              " --QTX 1 --only-1D 1"))

# Extract QTXNetwork -----------------------------------------------------------
load("./inst/Simulation-genphe.Rdata")
pre_file = "./inst/Simulation_SimTrait.pre"
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

save(qtl_data, qtl_env_data, file = "./inst/Simulation-qtl.Rdata")
