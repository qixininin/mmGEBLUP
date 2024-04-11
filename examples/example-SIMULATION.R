#### Install package -----------------------------------------------------------
## Install from github
library(devtools)
install_github("qixininin/mmGEBLUP")
## Install from source files
install.packages("/public3/zqx/mmGEBLUP_0.1.0.tar.gz", repos = NULL, type = "source")

#### Load package --------------------------------------------------------------
library(mmGEBLUP)

#### Data simulation -----------------------------------------------------------
# Parameter setting
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


save(geno_data, pheno_data, file = "./inst/Simulation-genphe.Rdata")
