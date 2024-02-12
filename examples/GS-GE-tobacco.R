setwd("/public3/zqx/YK_RILs/")
source("/public3/zqx/GS-GE/GS-GE-function.R")

# Input arguments --------------------------------------------------------------
## Examples nohup Rscript GS-GE-tobacco.R 1 1 YK_RILs_smok > exp1-rep1.log 2>&1 &
arg = commandArgs(T)
expNo = as.numeric(arg[1])
rep = as.numeric(arg[2])
prefix = arg[3]
# expNo = 1
# rep = 1
# prefix = "YK_RILs_agro"

# Load packages ----------------------------------------------------------------
library(sommer)
library(dplyr)

# Load data --------------------------------------------------------------------
load(paste0(prefix, ".Rdata"))

# Set parameters ---------------------------------------------------------------
cvNum = 4
bin_size = 0
set.seed(12345+rep-1)

# Trait category ---------------------------------------------------------------
if(prefix == "YK_RILs_agro"){
  ## agronomic traits
  expList = list(list(trainEnv = c("2018_SL","2019_SL","2020_SL","2021_SL","2018_YH","2019_YH","2020_YH"), validEnv = c("2022_SL")), # Exp1
                 list(trainEnv = c("2018_SL","2019_SL","2020_SL","2021_SL")                              , validEnv = c("2022_SL")), # Exp2
                 list(trainEnv = c("2018_SL","2019_SL","2020_SL","2018_YH","2019_YH","2020_YH")          , validEnv = c("2021_SL")), # Exp3
                 list(trainEnv = c("2018_SL","2019_SL","2020_SL")                                        , validEnv = c("2021_SL")), # Exp4
                 list(trainEnv = c("2018_SL","2019_SL","2020_SL","2018_YH","2019_YH","2020_YH")          , validEnv = c("2021_SL","2022_SL")), # Exp5
                 list(trainEnv = c("2018_SL","2019_SL","2020_SL")                                        , validEnv = c("2021_SL","2022_SL"))) # Exp6
} else if (prefix == "YK_RILs_smok"){
  ## smoke traits
  expList = list(list(trainEnv = c("2018_SL","2019_SL"), validEnv = c("2020_SL")), # Exp1
                 list(trainEnv = c("2018_SL")          , validEnv = c("2020_SL")), # Exp2
                 list(trainEnv = c("2019_SL")          , validEnv = c("2020_SL")), # Exp3
                 list(trainEnv = c("2018_SL")          , validEnv = c("2019_SL")), # Exp4
                 list(trainEnv = c("2018_SL")          , validEnv = c("2019_SL","2020_SL"))) # Exp5
} else {
  ## chem traits
  expList = list(list(trainEnv = head(levels(pheno_data$ENV),-1), validEnv = levels(pheno_data$ENV)[length(levels(pheno_data$ENV))])) # Exp1
}


# Summary ----------------------------------------------------------------------
modelName = c("GBLUP","mmGBLUP","GEBLUP","mmGEBLUP")
modelNum = length(modelName)
traitName = colnames(pheno_data)[-c(1:2)]
traitNum = length(traitName)
trainEnv = expList[[expNo]]$trainEnv
validEnv = expList[[expNo]]$validEnv
trialName = c(trainEnv,validEnv)
trialNum = length(trialName)
lineName = intersect(unique(pheno_data$GID),rownames(geno_data))
lineNum = length(lineName)

# Prepare genotype
Ga = geno_data ; rm(geno_data)
Ga = as.matrix(Ga[lineName, ])

# Additive relationship matrix
A = A.mat(Ga)
colnames(A) = rownames(A) = rownames(Ga)
E = diag(trialNum)
rownames(E) = colnames(E) = trialName
EA = kronecker(E, A, make.dimnames = TRUE)

# Cross validation set
cvSet = matrix(sample(lineName), nrow = cvNum)

write.table(cvSet, file = paste0("cvFold-",prefix, "-rep",rep,".txt"), quote = F, row.names = F, col.names = F)

# GBLUP ------------------------------------------------------------------------
a <- 0
gblup_list <- list()
gblup_bv_list <- list()
for(c in 1:traitNum) # loop for traits
{
  trait = traitName[c]
  y = pheno_data %>% dplyr::filter(ENV %in% trialName) %>%
    dplyr::filter(GID %in% lineName) %>%
    dplyr::select(c("ENV","GID",trait)) %>%
    droplevels()
  # replace specific trait name with just word "trait"
  colnames(y)[colnames(y) == trait] <- deparse(substitute(trait))
  
  for(i in 1:cvNum)
  {
    a <- a + 1
    cv = unique(cvSet[i,])
    # set phenotype in validation set and in validation env to be NA
    dt = y %>% dplyr::mutate(trait = ifelse(GID %in% cv & ENV %in% validEnv, NA, trait))
    dt = as.data.frame(dt)
    
    # GBLUP model
    rst = GBLUP(data = dt, Ka = A)
    BV = rst[[2]]
    
    # Calculate correlation
    cor <- BV %>% dplyr::mutate(obs = y$trait) %>%
      dplyr::filter(GID %in% cv) %>%
      dplyr::filter(ENV %in% validEnv) %>%
      dplyr::summarise(cor(obs,pre,use="pairwise.complete.obs")) %>%
      as.numeric()
    
    gblup_list[[a]] = data.frame(TRAIT = trait, CV = i, COR = cor, R2 = cor^2)
    gblup_bv_list[[a]] = BV %>% dplyr::filter(GID %in% cv) %>% dplyr::filter(ENV %in% validEnv)
    
    print(a)
  }
}

rstGBLUP <- dplyr::bind_rows(gblup_list)
save(rstGBLUP, gblup_bv_list, file = paste0(prefix,"-exp", expNo, "-GBLUP-rep",rep,"-rst.Rdata"))

print("NOTE: GBLUP DONE")

# mmGBLUP ----------------------------------------------------------------------
a <- 0
mmgblup_list <- list()
mmgblup_bv_list <- list()
for(c in 1:traitNum) # loop for traits
{
  trait = traitName[c]
  y = pheno_data %>% dplyr::filter(ENV %in% trialName) %>%
    dplyr::filter(GID %in% lineName) %>%
    dplyr::select(c("ENV","GID",trait)) %>%
    droplevels()
  # replace specific trait name with just word "trait"
  colnames(y)[colnames(y) == trait] <- deparse(substitute(trait))
  
  for(i in 1:cvNum) # loop for cross validation fold
  {
    a <- a + 1
    # extract qtl information from what we got from GWAS
    if(nrow(qtl_data %>% dplyr::filter(TRAIT == trait))>0){
      site_qtl = qtl_data %>% dplyr::filter(TRAIT == trait)
      m_qtl = nrow(site_qtl)
    } else {
      m_qtl = 0
    }
    
    # set phenotype in validation set and in validation env to be NA
    cv = unique(cvSet[i,])
    dt = y %>% dplyr::mutate(trait = ifelse(GID %in% cv & ENV %in% validEnv, NA, trait))
    dt = as.data.frame(dt)
    
    if(m_qtl>0){
      # adjusted additive relationship matrix
      Ka = calculateKa(Ga = Ga, map_data = marker_data, markers = site_qtl$QTL, bin = bin_size)
    }
    
    # mmGBLUP model
    if(m_qtl>0){
      dtnew = cbind(dt, Ga[y$GID, site_qtl$QTL])
      rst = mmGBLUP(data = dtnew, Ka = Ka)
      BV = rst[[2]]
    } else {
      mmgblup_list[[a]] = rstGBLUP[a,]
      mmgblup_bv_list[[a]] = gblup_bv_list[[a]]
      
      print(a)
      next
    }
    
    # Calculate correlation
    cor <- BV %>% dplyr::mutate(obs = y$trait) %>%
      dplyr::filter(GID %in% cv) %>% 
      dplyr::filter(ENV %in% validEnv) %>% 
      dplyr::summarise(cor(obs,pre,use="pairwise.complete.obs")) %>% 
      as.numeric()
    
    mmgblup_list[[a]] = data.frame(TRAIT = trait, CV = i, COR = cor, R2 = cor^2)
    mmgblup_bv_list[[a]] = BV %>% dplyr::filter(GID %in% cv) %>% dplyr::filter(ENV %in% validEnv)
    
    print(a)
  }
}

rstmmGBLUP <- dplyr::bind_rows(mmgblup_list)
save(rstmmGBLUP, mmgblup_bv_list, file = paste0(prefix,"-exp", expNo, "-mmGBLUP-rep",rep,"-rst.Rdata"))

print("NOTE: mmGBLUP DONE")

# GEBLUP -----------------------------------------------------------------------
a <- 0
geblup_list <- list()
geblup_bv_list <- list()
for(c in 1:traitNum) # loop for traits
{
  trait = traitName[c]
  y = pheno_data %>% dplyr::filter(ENV %in% trialName) %>%
    dplyr::filter(GID %in% lineName) %>%
    dplyr::select(c("ENV","GID",trait)) %>%
    droplevels()
  # replace specific trait name with just word "trait"
  colnames(y)[colnames(y) == trait] <- deparse(substitute(trait))
  
  for(i in 1:cvNum)
  {
    a <- a + 1
    cv = unique(cvSet[i,])
    # set phenotype in validation set and in validation env to be NA
    dt = y %>% dplyr::mutate(trait = ifelse(GID %in% cv & ENV %in% validEnv, NA, trait))
    dt = as.data.frame(dt)
    
    # GEBLUP model
    rst = GEBLUP(data = dt, Ka = A, EKae = EA)
    BV = rst[[2]]
    
    # Calculate correlation
    cor <- BV %>% dplyr::mutate(obs = y$trait) %>%
      dplyr::filter(GID %in% cv) %>%
      dplyr::filter(ENV %in% validEnv) %>%
      dplyr::summarise(cor(obs,pre,use="pairwise.complete.obs")) %>%
      as.numeric()
    
    geblup_list[[a]] = data.frame(TRAIT = trait, CV = i, COR = cor, R2 = cor^2)
    geblup_bv_list[[a]] = BV %>% dplyr::filter(GID %in% cv) %>% dplyr::filter(ENV %in% validEnv)
    
    print(a)
  }
}

rstGEBLUP <- dplyr::bind_rows(geblup_list)
save(rstGEBLUP, geblup_bv_list, file = paste0(prefix,"-exp", expNo, "-GEBLUP-rep",rep,"-rst.Rdata"))

print("NOTE: GEBLUP DONE")

# mmGEBLUP ------------------------------------------------------
a <- 0
mmgeblup_list <- list()
mmgeblup_bv_list <- list()
for(c in 1:traitNum) # loop for traits
{
  trait = traitName[c]
  y = pheno_data %>% dplyr::filter(ENV %in% trialName) %>%
    dplyr::filter(GID %in% lineName) %>%
    dplyr::select(c("ENV","GID",trait)) %>%
    droplevels()
  # replace specific trait name with just word "trait"
  colnames(y)[colnames(y) == trait] <- deparse(substitute(trait))
  
  for(i in 1:cvNum) # loop for cross validation fold
  {
    a <- a + 1
    # extract qtl information from what we got from GWAS
    if(nrow(qtl_data %>% dplyr::filter(TRAIT == trait))>0){
      site_qtl = qtl_data %>% dplyr::filter(TRAIT == trait)
      m_qtl = nrow(site_qtl)
    } else {
      m_qtl = 0
    }
    if(nrow(qtl_env_data %>% dplyr::filter(TRAIT == trait))>0){
      site_env_qtl_all = qtl_env_data %>% dplyr::filter(TRAIT == trait) %>% dplyr::filter(ENV %in% trialName)
      m_env_qtl = nrow(site_env_qtl_all)
    } else {
      m_env_qtl = 0
    }
    
    # set phenotype in validation set and in validation env to be NA
    cv = unique(cvSet[i,])
    dt = y %>% dplyr::mutate(trait = ifelse(GID %in% cv & ENV %in% validEnv, NA, trait))
    dt = as.data.frame(dt)
    
    if(m_qtl>0){
      # adjusted additive relationship matrix
      Ka = calculateKa(Ga = Ga, map_data = marker_data, markers = site_qtl$QTL, bin = bin_size)
    }
    
    if(m_env_qtl>0){
      EKae = calculateEKae(Ga = Ga, EA = EA, A = A, site_env_qtl_all = site_env_qtl_all, trialName = trialName)
      EKae_l = EKae[[1]]
      EKae_s = EKae[[2]]
    }
    
    # mmGEBLUP model
    if(m_qtl>0 & m_env_qtl>0){
      dtnew = cbind(dt, Ga[y$GID, site_qtl$QTL])
      rst = mmGEBLUP(data = dtnew, Ka = Ka, EKae1 = EKae_l, EKae2 = EKae_s)
      BV = rst[[2]]
      
    } else if (m_qtl>0 & m_env_qtl==0) {
      dtnew = cbind(dt, Ga[y$GID, site_qtl$QTL])
      rst = mmGEBLUP(data = dtnew, Ka = Ka, EKae1 = EA)
      BV = rst[[2]]
      
    } else if (m_qtl==0 & m_env_qtl>0) {
      dtnew = dt
      rst = mmGEBLUP(data = dtnew, Ka = A, EKae1 = EKae_l, EKae2 = EKae_s)
      BV = rst[[2]]
      
    } else {
      mmgeblup_list[[a]] = rstGEBLUP[a,]
      mmgeblup_bv_list[[a]] = geblup_bv_list[[a]]
      
      print(a)
      next
    }
    
    # Calculate correlation
    cor <- BV %>% dplyr::mutate(obs = y$trait) %>%
      dplyr::filter(GID %in% cv) %>% 
      dplyr::filter(ENV %in% validEnv) %>% 
      dplyr::summarise(cor(obs,pre,use="pairwise.complete.obs")) %>% 
      as.numeric()
    
    mmgeblup_list[[a]] = data.frame(TRAIT = trait, CV = i, COR = cor, R2 = cor^2)
    mmgeblup_bv_list[[a]] = BV %>% dplyr::filter(GID %in% cv) %>% dplyr::filter(ENV %in% validEnv)
    
    print(a)
  }
}

rstmmGEBLUP <- dplyr::bind_rows(mmgeblup_list)
save(rstmmGEBLUP, mmgeblup_bv_list, file = paste0(prefix,"-exp", expNo, "-mmGEBLUP-rep",rep,"-rst.Rdata"))

print("NOTE: mmGEBLUP DONE")

df = data.frame(rbind(rstGBLUP,rstmmGBLUP,rstGEBLUP,rstmmGEBLUP))
df$METHOD = rep(modelName, each=traitNum*cvNum)
df$TRAIT = factor(df$TRAIT, levels = unique(df$TRAIT))
save(df, file = paste0(prefix,"-exp", expNo, "-rep",rep,"-rst.Rdata"))