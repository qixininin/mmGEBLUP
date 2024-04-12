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
load("./inst/Pivot-genphe.Rdata")
pheno_output_prefix = "./inst/Pivot"
geno_output_prefix  = "./inst/Pivot"

## The default output names are
## .gen = Simulation.gen
## .phe = Simulation_SimTrait.phe
qtxnetwork.input.trans(geno_data, pheno_data,
                       geno_output_prefix = geno_output_prefix,
                       pheno_output_prefix = pheno_output_prefix)


# Perform QTXNetwork -----------------------------------------------------------
qtxnetwork_path = "/public3/zqx/QTXNetwork_4.0/build/QTXNetwork"
gen_file = "./inst/Pivot.gen"
phe_file = "./inst/Pivot_SimTrait.phe"
pre_file = "./inst/Pivot_SimTrait.pre"

qtxnetwork.perform(qtxnetwork_path, gen_file, phe_file, pre_file)

# Extract QTXNetwork -----------------------------------------------------------
qtl = qtxnetwork.output.trans(pheno_data, pre_file)
qtl_data = qtl$qtl_data
qtl_env_data = qtl$qtl_env_data

save(qtl_data, qtl_env_data, file = "./inst/Pivot-qtl.Rdata")
