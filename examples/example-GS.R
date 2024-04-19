# Load packages ----------------------------------------------------------------
library(mmGEBLUP)
library(dplyr)
set.seed(215)

# Load data --------------------------------------------------------------------
load("./inst/Pivot-genphe.Rdata")
load("./inst/Pivot-qtl.Rdata")

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
naNum = ifelse(length(mmdata$mmsummary$lineName)%%cvNum==0, 0, cvNum-length(mmdata$mmsummary$lineName)%%cvNum)
cvSet = matrix(c(sample(mmdata$mmsummary$lineName), rep(NA, naNum)), nrow = cvNum)

# mmGEBLUP ------------------------------------------------------
a <- 0
mmgeblup_list <- list()
mmgeblup_bv_list <- list()
for(i in 1:cvNum) # loop for cross validation fold
{
  a <- a + 1

  # set phenotype in validation set and in validation env to be NA
  cv = as.vector(na.omit(unique(cvSet[i,])))
  dt = mmpheno_data %>% dplyr::mutate(trait = ifelse(GID %in% cv & ENV %in% validEnv, NA, trait))
  dt = as.data.frame(dt)

  # mmGEBLUP model
  rst = mmgeblup(data = cbind(dt, mmdata$Xa), Ka = mmdata$Ka, KaeE1 = mmdata$KaeE1, KaeE2 = mmdata$KaeE2)
  BV = rst[[2]]

  # Calculate correlation
  cor <- BV %>% dplyr::mutate(obs = mmpheno_data$trait) %>%
    dplyr::filter(GID %in% cv) %>%
    dplyr::filter(ENV %in% validEnv) %>%
    dplyr::summarise(cor(obs,pre,use="pairwise.complete.obs")) %>%
    as.numeric()

  mmgeblup_list[[a]] = data.frame(TRAIT = mmdata$mmsummary$traitName, CV = i, COR = cor, R2 = cor^2)
  mmgeblup_bv_list[[a]] = BV %>%
    dplyr::mutate(obs = mmpheno_data$trait) %>%
    dplyr::filter(GID %in% cv) %>%
    dplyr::filter(ENV %in% validEnv)

  print(a)
}

rstmmGEBLUP <- dplyr::bind_rows(mmgeblup_list)
rstBV <- dplyr::bind_rows(mmgeblup_bv_list)

save(mmdata, rstmmGEBLUP, rstBV, cvSet, file = "./inst/Simulation-GSresult.Rdata")

# GEBLUP ------------------------------------------------------
a <- 0
geblup_list <- list()
geblup_bv_list <- list()
for(i in 1:cvNum) # loop for cross validation fold
{
  a <- a + 1

  # set phenotype in validation set and in validation env to be NA
  cv = as.vector(na.omit(unique(cvSet[i,])))
  dt = mmpheno_data %>% dplyr::mutate(trait = ifelse(GID %in% cv & ENV %in% validEnv, NA, trait))
  dt = as.data.frame(dt)

  # GEBLUP model
  rst = geblup(data = dt, A = mmdata$A, AE = mmdata$AE)
  BV = rst[[2]]

  # Calculate correlation
  cor <- BV %>% dplyr::mutate(obs = mmpheno_data$trait) %>%
    dplyr::filter(GID %in% cv) %>%
    dplyr::filter(ENV %in% validEnv) %>%
    dplyr::summarise(cor(obs,pre,use="pairwise.complete.obs")) %>%
    as.numeric()

  geblup_list[[a]] = data.frame(TRAIT = mmdata$mmsummary$traitName, CV = i, COR = cor, R2 = cor^2)

  print(a)
}

rstGEBLUP <- dplyr::bind_rows(geblup_list)
