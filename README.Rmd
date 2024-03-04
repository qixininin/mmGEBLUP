# mmGEBLUP

mmGEBLUP: A systematic simulation program and a novel genomic prediction strategy based on integrated evaluation of major and minor gene and gene-by-environment effects    

In our package, we offer **Three** key functions in the field of quantatative genetics including **Genotype-by-Environment (GE)** interaction.   
- Trait simulation    
- GWAS analysis    
- Genomic prediction (a novel model named mmGEBLUP)     

## Package installation
The current GitHub version of **mmGEBLUP** can be installed via:
```
library(devtools)
install_github("qixininin/mmGEBLUP")
```

## Load the library
```
library(mmGEBLUP)
```

Correspondingly, we offer **Three** examples as an illustration.    

## Example 1 Trait simulation 

```
# Load packages 
library(mmGEBLUP)


# Set parameters
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
```


## Example 2 GWAS analysis
```
# Load packages 
library(mmGEBLUP)
library(dplyr)

# Prepare QTXNetwork
load("./inst/Simulation-genphe.Rdata")
pheno_output_prefix = "./inst/Simulation"
geno_output_prefix  = "./inst/Simulation"

## The default output names are
## .gen = Simulation.gen
## .phe = Simulation_SimTrait.phe
qtxnetwork.input.trans(geno_data, pheno_data,
                       geno_output_prefix = geno_output_prefix,
                       pheno_output_prefix = pheno_output_prefix)


# Perform QTXNetwork
qtxnetwork_path = "/public3/zqx/QTXNetwork_4.0/build/QTXNetwork"
gen_file = "./inst/Simulation.gen"
phe_file = "./inst/Simulation_SimTrait.phe"
pre_file = "./inst/Simulation_SimTrait.pre"

qtxnetwork.perform(qtxnetwork_path, gen_file, phe_file, pre_file)

# Extract QTXNetwork
qtl = qtxnetwork.output.trans(pheno_data, pre_file)
qtl_data = qtl$qtl_data
qtl_env_data = qtl$qtl_env_data

save(qtl_data, qtl_env_data, file = "./inst/Simulation-qtl.Rdata")
```


## Exmaple 3 Genomic prediction (a novel model named mmGEBLUP)

```
# Load packages 
library(mmGEBLUP)
library(dplyr)

# Load data 
load("./inst/Simulation-genphe.Rdata")
load("./inst/Simulation-qtl.Rdata")

# geno_data[1:5,1:10]
# pheno_data[1:5,]
# qtl_data[1:5,]
# qtl_env_data[1:5,]

# Data QC and Summary ----------------------------------------------------------
mmdata = mmdata(geno_data, pheno_data, qtl_data, qtl_env_data)
mmgeno_data = mmdata$mmgeno_data
mmpheno_data = mmdata$mmpheno_data

# mmgeno_data[1:5,1:10]
# mmpheno_data[1:5,]


# Set parameters ---------------------------------------------------------------
trainEnv = c("env1","env2")
validEnv = c("env3")
cvNum = 4
cvSet = matrix(sample(mmdata$mmsummary$lineName), nrow = cvNum)

# mmGEBLUP ------------------------------------------------------
a <- 0
mmgeblup_list <- list()
mmgeblup_bv_list <- list()
for(i in 1:cvNum) # loop for cross validation fold
{
  a <- a + 1

  # set phenotype in validation set and in validation env to be NA
  cv = unique(cvSet[i,])
  dt = mmpheno_data %>% dplyr::mutate(trait = ifelse(GID %in% cv & ENV %in% validEnv, NA, trait))
  dt = as.data.frame(dt)

  # mmGEBLUP model
  rst = mmgeblup(data = cbind(dt, mmdata$Xa), Ka = mmdata$Ka, EKae1 = mmdata$EKae1, EKae2 = mmdata$EKae2)
  BV = rst[[2]]

  # Calculate correlation
  cor <- BV %>% dplyr::mutate(obs = mmpheno_data$trait) %>%
    dplyr::filter(GID %in% cv) %>%
    dplyr::filter(ENV %in% validEnv) %>%
    dplyr::summarise(cor(obs,pre,use="pairwise.complete.obs")) %>%
    as.numeric()

  mmgeblup_list[[a]] = data.frame(TRAIT = mmdata$mmsummary$traitName, CV = i, COR = cor, R2 = cor^2)
  mmgeblup_bv_list[[a]] = BV %>% dplyr::filter(GID %in% cv) %>% dplyr::filter(ENV %in% validEnv)

  print(a)
}


rstmmGEBLUP <- dplyr::bind_rows(mmgeblup_list)

print("NOTE: mmGEBLUP DONE")

save(rstmmGEBLUP, file = "./inst/Simulation-GSresult.Rdata")
```
