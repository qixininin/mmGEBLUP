setwd("/public3/zqx/IRRI-A")
source("/public3/zqx/GS-GE/GS-GE-function.R")

# Input arguments --------------------------------------------
## Examples nohup Rscript GS-GE-IRRI-A.R 1 > irri-exp1.log 2>&1 &
arg = commandArgs(T)
expNo = as.numeric(arg[1])

# Load packages ----------------------------------------------
library(dplyr)
library(sommer)

# Load data --------------------------------------------------
prefix = "IRRI_RYT_8756mk"
load(paste0(prefix, "_forGS.Rdata")) ## Load geno_data, pheno_data
load("irri-a-cvFoldDesign.Rdata") ## Load cv fold information (cvSet)

# Set parameters ---------------------------------------------
cvNum = 5
alpha = 0.05
bin_size = 0
set.seed(12345)

expList = list(list(trainEnv = c("2009-1","2009-2","2010-1","2010-2","2011-1","2011-2")         , validEnv = c("2012-1")), # Exp1
               list(trainEnv = c("2009-1","2009-2","2010-1","2010-2","2011-1","2011-2","2012-1"), validEnv = c("2012-2")), # Exp2
               list(trainEnv = c("2009-1","2009-2","2010-1","2010-2","2011-1","2011-2")         , validEnv = c("2012-2")), # Exp3
               list(trainEnv = c("2009-1","2010-1","2011-1")                                    , validEnv = c("2012-1")), # Exp4
               list(trainEnv = c("2009-2","2010-2","2011-2")                                    , validEnv = c("2012-2")), # Exp5
               list(trainEnv = c("2011-1","2011-2")                                             , validEnv = c("2012-1")), # Exp6
               list(trainEnv = c("2011-1","2011-2")                                             , validEnv = c("2012-2")), # Exp7
               list(trainEnv = c("2011-1","2011-2","2012-1")                                    , validEnv = c("2012-2")), # Exp8
               list(trainEnv = c("2010-1","2010-2","2011-1","2011-2")                           , validEnv = c("2012-1")), # Exp9
               list(trainEnv = c("2010-1","2010-2","2011-1","2011-2")                           , validEnv = c("2012-2")), # Exp10
               list(trainEnv = c("2010-1","2010-2","2011-1","2011-2","2012-1")                  , validEnv = c("2012-2")), # Exp11
               list(trainEnv = c("2010-2","2011-1","2011-2")                                    , validEnv = c("2012-1")))  # Exp12

# summary
modelName = c("GBLUP","mmGBLUP","GEBLUP","mmGEBLUP")
modelNum = length(modelName)
trainEnv = expList[[expNo]]$trainEnv
validEnv = expList[[expNo]]$validEnv
trialName = c(trainEnv,validEnv)
trialNum = length(trialName)
if("2009-2" %in% trainEnv){ # 2009-2 miss trait - PH
  traitName = c("YLD","FL")
}else{ 
  traitName = c("YLD","PH","FL")
} 
traitNum = length(traitName)

# Unifying line id
lineName = intersect(unique(pheno_data$GID),rownames(geno_data))
lineNum = length(lineName)

# Factorize phenotype data
pheno_data$ENV = as.factor(pheno_data$ENV)
pheno_data$GID = as.factor(pheno_data$GID)

# Prepare genotype (-1,0,1)
Ga = geno_data ; rm(geno_data)
Ga = as.matrix(Ga[lineName, ])

# Additive relationship matrix
A = A.mat(Ga)
colnames(A) = rownames(A) = rownames(Ga)
E = diag(trialNum)
rownames(E) = colnames(E) = trialName
EA = kronecker(E, A, make.dimnames = TRUE)

# # GWAS ---------------------------------------------------------------
# a <- 0
# b <- 0
# qtl_list <- list()
# env_qtl_list <- list()
# pvalue_list <- list()
# env_pvalue_list <- list()
# for(c in 1:traitNum) # loop for traits
# {
#   trait = traitName[c]
#   y = pheno_data %>% dplyr::filter(ENV %in% trialName) %>%
#     dplyr::filter(GID %in% lineName) %>%
#     dplyr::select(c("ENV","GID",trait)) %>%
#     droplevels()
# 
#   # replace specific trait name with just word "trait"
#   colnames(y)[colnames(y) == trait] <- deparse(substitute(trait))
# 
#   for(i in 1:cvNum) # loop for cross validation fold
#   {
#     a <- a + 1
#     cv = unique(cvSet[i,])
#     # set phenotype in validation set and in validation env to be NA
#     dt = y %>% dplyr::mutate(trait = ifelse(GID %in% cv & ENV %in% validEnv, NA, trait))
#     dt = as.data.frame(dt)
#     # model
#     mod0 = GWAS(trait~ENV,
#                 random = ~vsr(GID, Gu=A) + vsr(ENV:GID, Gu=EA),
#                 rcov = ~units,
#                 data = dt,
#                 M = Ga,
#                 gTerm = "u:GID",
#                 method = "NR", min.MAF = 0.05, naMethodY="exclude", # default
#                 verbose = FALSE, date.warning = FALSE)
# 
#     # extract p-value
#     pvalue = data.frame(snp = rownames(mod0$scores),
#                         oriP = as.vector(ifelse(10^-mod0$scores==0,NA,10^-mod0$scores)))
#     # reorder p-value
#     pvalue = pvalue[order(pvalue$oriP),]
#     # FDR
#     pvalue$adjP = p.adjust(pvalue$oriP, method = "BH")
#     site_qtl = pvalue[which(pvalue$adjP<alpha),"snp"]
#     if(length(site_qtl)>0) qtl_list[[a]] = data.frame(TRAIT = trait,
#                                                       CV = i,
#                                                       QTL = site_qtl)
#     pvalue_list[[a]] = pvalue
# 
#     # additive-by-env gwas
#     for(e in 1:trialNum)
#     {
#       b <- b + 1
#       env = trialName[e]
#       dte = y %>% dplyr::filter(ENV==env) %>% dplyr::mutate(trait = ifelse(GID %in% cv & ENV %in% validEnv, NA, trait))
#       # model
#       mod0 = GWAS(trait~1,
#                   random = ~ vsr(GID, Gu=A),
#                   rcov = ~ units,
#                   data = dte,
#                   M = Ga,
#                   gTerm = "u:GID",
#                   verbose = FALSE, date.warning = FALSE)
# 
#       # extract p-value
#       pvalue = data.frame(snp = rownames(mod0$scores),
#                           oriP = as.vector(ifelse(10^-mod0$scores==0,NA,10^-mod0$scores)))
#       # reorder p-value
#       pvalue = pvalue[order(pvalue$oriP),]
#       # FDR
#       pvalue$adjP = p.adjust(pvalue$oriP, method = "BH")
#       site_env_qtl = pvalue[which(pvalue$adjP<alpha),"snp"]
#       if(length(site_env_qtl)>0) env_qtl_list[[b]] = data.frame(TRAIT = trait,
#                                                                 ENV = env,
#                                                                 CV = i,
#                                                                 QTL = site_env_qtl)
#       env_pvalue_list[[b]] = pvalue
#     }
#   }
# }
# 
# gwas_qtl <- dplyr::bind_rows(qtl_list)
# gwas_env_qtl <- dplyr::bind_rows(env_qtl_list)
# print("NOTE: GWAS DONE")
# 
# save(gwas_qtl,gwas_env_qtl,pvalue_list, env_pvalue_list, file = paste0(prefix,"-exp", expNo, "-qtl.Rdata"))

# # QTS selection ------------------------------------------------------
# load(paste0(prefix,"-exp", expNo, "-qtl.Rdata"))
# # Extract chromosome and position
# gwas_qtl$CHR=as.numeric(unlist(strsplit(gsub("S","",gwas_qtl$QTL),"_"))[seq(1,2*nrow(gwas_qtl),2)])
# gwas_qtl$BP=as.numeric(unlist(strsplit(gsub("S","",gwas_qtl$QTL),"_"))[seq(1,2*nrow(gwas_qtl),2)+1])
# gwas_env_qtl$CHR=as.numeric(unlist(strsplit(gsub("S","",gwas_env_qtl$QTL),"_"))[seq(1,2*nrow(gwas_env_qtl),2)])
# gwas_env_qtl$BP=as.numeric(unlist(strsplit(gsub("S","",gwas_env_qtl$QTL),"_"))[seq(1,2*nrow(gwas_env_qtl),2)+1])
# 
# # QTL that is treated as large-effect QTL
# le_qtl <-  gwas_qtl %>%
#   mutate(bin_id = floor(BP / bin_size)) %>%
#   group_by(TRAIT, CV, CHR, bin_id) %>%
#   filter(row_number() == 1) %>%
#   ungroup()
# 
# # QTL that is treated as large-effect envQTL
# le_env_qtl <-  gwas_env_qtl %>%
#   mutate(bin_id = floor(BP / bin_size)) %>%
#   group_by(TRAIT, CV, ENV, CHR, bin_id) %>%
#   filter(row_number() == 1) %>%
#   ungroup()


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
save(rstGBLUP, gblup_bv_list, file = paste0(prefix,"-exp", expNo, "-GBLUP-rst.Rdata"))

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
save(rstmmGBLUP, mmgblup_bv_list, file = paste0(prefix,"-exp", expNo, "-mmGBLUP-rst.Rdata"))

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
save(rstGEBLUP, geblup_bv_list, file = paste0(prefix,"-exp", expNo, "-GEBLUP-rst.Rdata"))

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
save(rstmmGEBLUP, mmgeblup_bv_list, file = paste0(prefix,"-exp", expNo, "-mmGEBLUP-rst.Rdata"))

df = data.frame(rbind(rstGBLUP,rstmmGBLUP,rstGEBLUP,rstmmGEBLUP))
df$METHOD = rep(modelName, each=traitNum*cvNum)
df$TRAIT = factor(df$TRAIT, levels = unique(df$TRAIT))

save(df, file = paste0(prefix,"-exp", expNo, "-rst.Rdata"))