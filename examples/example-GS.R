# Load packages ----------------------------------------------------------------
library(mmGEBLUP)
library(sommer)
library(dplyr)

# Load data --------------------------------------------------------------------
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

  mmgeblup_list[[a]] = data.frame(TRAIT = trait, CV = i, COR = cor, R2 = cor^2)
  mmgeblup_bv_list[[a]] = BV %>% dplyr::filter(GID %in% cv) %>% dplyr::filter(ENV %in% validEnv)

  print(a)
}


rstmmGEBLUP <- dplyr::bind_rows(mmgeblup_list)

print("NOTE: mmGEBLUP DONE")

save(rstmmGEBLUP, file = paste0(prefix,"-rep", rep, "-rst.Rdata"))
