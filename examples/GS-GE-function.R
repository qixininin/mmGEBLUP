# Simulation -------------------------------------------------------------------
# simulate_effects_with_GxE
simulate_effects_with_GxE <- function(snpNum, envNum, prop, mode, major_idx, 
                                      variance_a_major, variance_ae_major,
                                      variance_a_minor, variance_ae_minor) {
  
  # The mode indicate whether the major QTL has either of A or AE or both A+AE effects
  if (!mode %in% c("A","AE","A+AE") ) {
    stop("Invalid mode: choose 'A', 'AE', or 'A+AE'")
  }
  
  # Select major effects
  # major_idx <- sample(1:snpNum, size = round(snpNum * prop), replace = FALSE)
  
  
  main_effects     <- matrix(0, nrow = 1     , ncol = snpNum)
  interact_effects <- matrix(0, nrow = envNum, ncol = snpNum)
  
  # Simulate minor gene effect
  main_effects[1,-major_idx] <- rnorm(snpNum-length(major_idx), mean = 0, sd = sqrt(variance_a_minor))
  for (i in 1:envNum) {
    interact_effects[i,-major_idx] <- rnorm(snpNum-length(major_idx), mean = 0, sd = sqrt(variance_ae_minor))
  }
  
  # Simulate major gene effect
  if(mode=="A"){
    main_effects[1, major_idx] <- rnorm(length(major_idx), mean = 0, sd = sqrt(variance_a_major))
    for (i in 1:envNum) {
      interact_effects[i, major_idx] <- rnorm(length(major_idx), mean = 0, sd = sqrt(variance_ae_minor))
    }    
  } else if (mode=="AE") {
    main_effects[1, major_idx] <- rnorm(length(major_idx), mean = 0, sd = sqrt(variance_a_minor))
    for (i in 1:envNum) {
      interact_effects[i, major_idx] <- rnorm(length(major_idx), mean = 0, sd = sqrt(variance_ae_major))
    }
  } else {
    main_effects[1, major_idx] <- rnorm(length(major_idx), mean = 0, sd = sqrt(variance_a_major))
    for (i in 1:envNum) {
      interact_effects[i, major_idx] <- rnorm(length(major_idx), mean = 0, sd = sqrt(variance_ae_major))
    }
  }
  
  # Return effects matrix and major markers based on mode
  return(list(main_effects = main_effects, 
              interact_effects = interact_effects))
}

plot_effects_with_GxE <- function(dt, lmt)
{
  p1 = ggplot(dt,aes(x=x,y=y,color=grp))+
    geom_point(cex = 0.8) +
    scale_y_continuous(limits = c(-lmt,lmt))+
    scale_color_manual(values = c("orange","lightblue"), name = "Simulated effect size") + 
    facet_grid(~env) +
    theme_bw() + 
    labs(x="SNP",y="Main effect")
  dt = data.frame(x=rep(1:snpNum,envNum),y=as.vector(t(bh)),grp=rep(grp,envNum),env=rep(paste0("Env",1:envNum),each=snpNum))
  p2 = ggplot(dt,aes(x=x,y=y,color=grp))+
    geom_point(cex = 0.8) +
    scale_y_continuous(limits = c(-lmt,lmt))+
    scale_color_manual(values = c("orange","lightblue"), name = "Simulated effect size") + 
    facet_grid(~env) +
    theme_bw() +
    labs(x="SNP",y="Intersect effect")
  dt = data.frame(x=rep(1:snpNum,envNum),y=as.vector(t(b)),grp=rep(grp,envNum),env=rep(paste0("Env",1:envNum),each=snpNum))
  p3 = ggplot(dt,aes(x=x,y=y,color=grp))+
    geom_point(cex = 0.8) +
    scale_y_continuous(limits = c(-lmt,lmt))+
    scale_color_manual(values = c("orange","lightblue"), name = "Simulated effect size") + 
    facet_grid(~env) +
    theme_bw() + 
    labs(x="SNP",y="Overall effect")
  p = ggpubr::ggarrange(p1,p2,p3,nrow=3,common.legend = T, legend = "right")
  return(p)
}


# Model function ---------------------------------------------------------------
GBLUP <- function(data, Ka)
{
  ## Evaluate input
  if(missing(Ka)){
    stop("Error: no input for additive kinship matrix.")
  }
  
  mod = mmer(trait~1,
             random = ~vsr(GID, Gu=Ka) + vsr(ENV),
             rcov = ~units,
             data = data,
             verbose = FALSE, date.warning = FALSE)
  ## Predict
  mu = mod$Beta$Estimate 
  BV = data.frame(data[,c("ENV","GID")])
  BV$pre = mu +                                          # mu
    mod$U$`u:GID`$trait[BV$GID]+                         # G
    mod$U$`u:ENV`$trait[BV$ENV]                          # E
  
  return(list(mod, BV))
}

mmGBLUP <- function(data, Ka)
{
  ## Evaluate input
  if(missing(Ka)){
    stop("Error: no input for additive kinship matrix.")
  }
  
  ## Receive fixed effect column names
  ## If no fixed effect columns in data, then only intercept is fixed effect
  fixColName = names(data)[!(names(data) %in% c("ENV","GID","trait"))]
  if(length(fixColName)==0){
    fixColName = "1"
  }
  
  ## Perform sommer models
  mod = mmer(reformulate(fixColName, "trait"),
             random = ~vsr(GID, Gu=Ka) + vsr(ENV),
             rcov = ~units,
             data = data,
             verbose = FALSE, date.warning = FALSE)
  ## Predict
  mu = as.matrix(cbind(rep(1, nrow(data)), data[,as.vector(mod$Beta$Effect)[-1]])) %*% mod$Beta$Estimate 
  BV = data.frame(data[,c("ENV","GID")])
  BV$pre = mu +                                          # mu
    mod$U$`u:GID`$trait[BV$GID]+                         # G
    mod$U$`u:ENV`$trait[BV$ENV]                          # E
  
  return(list(mod, BV))
}

GEBLUP <- function(data, Ka, EKae)
{
  ## Evaluate input
  if(missing(Ka) | missing(EKae)){
    stop("Error: no input for additive kinship matrix or additive-by-environment kinship matrix.")
  }
  
  mod = mmer(trait~1,
             random = ~vsr(GID, Gu=Ka) + vsr(ENV) + vsr(ENV:GID, Gu=EKae),
             rcov = ~units,
             data = data,
             verbose = FALSE, date.warning = FALSE)
  ## Predict
  mu = mod$Beta$Estimate 
  BV = data.frame(data[,c("ENV","GID")])
  BV$pre = mu +                                          # mu
    mod$U$`u:GID`$trait[BV$GID]+                         # G
    mod$U$`u:ENV`$trait[BV$ENV]+                         # E
    mod$U$`u:ENV:GID`$trait[paste0(BV$ENV,":",BV$GID)]   # GE
  
  return(list(mod, BV))
}

mmGEBLUP <- function(data, Ka, EKae1, EKae2)
{
  ## Evaluate input
  if(missing(Ka) | missing(EKae1)){
    stop("Error: no input for additive kinship matrix or additive-by-environment kinship matrix.")
  }
  majorGE = TRUE
  if(missing(EKae2)){
    majorGE = FALSE
  }
  
  ## Receive fixed effect column names
  ## If no fixed effect columns in data, then only intercept is fixed effect
  fixColName = names(data)[!(names(data) %in% c("ENV","GID","trait"))]
  if(length(fixColName)==0){
    fixColName = "1"
  }
  
  ## Data prepare
  datafake = data %>% dplyr::mutate(GID1 = GID)
  
  ## Perform sommer models
  if(majorGE){ # If two AE kinship matrices are inputted
    mod = mmer(reformulate(fixColName, "trait"),
               random = ~vsr(GID, Gu=Ka) + vsr(ENV) + vsr(ENV:GID, Gu=EKae1) + vsr(ENV:GID1, Gu=EKae2),
               rcov = ~units,
               data = datafake,
               verbose = FALSE, date.warning = FALSE)
    ## Predict
    mu = as.matrix(cbind(rep(1, nrow(data)), data[,as.vector(mod$Beta$Effect)[-1]])) %*% mod$Beta$Estimate 
    BV = data.frame(data[,c("ENV","GID")])
    BV$pre = mu +                                          # mu
      mod$U$`u:GID`$trait[BV$GID]+                         # G
      mod$U$`u:ENV`$trait[BV$ENV]+                         # E
      mod$U$`u:ENV:GID`$trait[paste0(BV$ENV,":",BV$GID)]+  # GE-major
      mod$U$`u:ENV:GID1`$trait[paste0(BV$ENV,":",BV$GID)]  # GE-minor
    
  } else { # If only one AE kinship matrices is inputted
    mod = mmer(reformulate(fixColName, "trait"),
               random = ~vsr(GID, Gu=Ka) + vsr(ENV) + vsr(ENV:GID, Gu=EKae1),
               rcov = ~units,
               data = data,
               verbose = FALSE, date.warning = FALSE)
    ## Predict
    mu = as.matrix(cbind(rep(1, nrow(data)), data[,as.vector(mod$Beta$Effect)[-1]])) %*% mod$Beta$Estimate 
    BV = data.frame(data[,c("ENV","GID")])
    BV$pre = mu +                                          # mu
      mod$U$`u:GID`$trait[BV$GID]+                         # G
      mod$U$`u:ENV`$trait[BV$ENV]+                         # E
      mod$U$`u:ENV:GID`$trait[paste0(BV$ENV,":",BV$GID)]   # GE
  }
  
  return(list(mod, BV))
  
}



find_LD_markers <- function(map_data, markers, bin) {
  result <- data.frame() 
  
  for (i in 1:length(markers)) {
    
    index <- which(map_data$SNP == markers[i])
    
    if (length(index) > 0) {
      
      chromosome <- map_data$CHR[index]
      position <- map_data$Distance[index]
      
      nearby_markers <- map_data[map_data$CHR == chromosome & abs(map_data$Distance - position) <= bin, ]
      
      result <- rbind(result, nearby_markers)
    }
  }
  
  return(unique(result$SNP))
}

calculateKa <- function(Ga, map_data, markers, bin)
{
  # NEW KINSHIP MATRIX
  # !! REMOVE ASOCIATED MARKERS
  if(bin<0) {
    stop("Error: bin size should not be minus.")
  }
  
  if(bin>0)
  {
    # find all markers within bin windows
    remove_marker = find_LD_markers(map_data = map_data, markers = markers, bin = bin)
    Ka = A.mat(as.matrix(Ga[,!colnames(Ga) %in% remove_marker]))
  } else {   # if bin=0, then only QTS is removed in calculate Ka
    Ka = A.mat(as.matrix(Ga[,!colnames(Ga) %in% markers]))
  }
  colnames(Ka) = rownames(Ka) = rownames(Ga)
  
  return(Ka)
}

calculateEKae <- function(Ga, EA, A, site_env_qtl_all, trialName)
{
  EKae_l = EA
  EKae_s = EA
  trialNum = length(trialName)
  lineNum = nrow(A)
  
  for(e in 1:trialNum)
  {
    env = trialName[e]
    left = (e-1)*lineNum+1
    right = e*lineNum
    site_env_qtl = site_env_qtl_all %>% dplyr::filter(ENV == env)
    
    if(nrow(site_env_qtl)>0)
    {
      # adjusted additive-by-env relationship matrix
      EKae_l[left:right,left:right] = A.mat(as.matrix(Ga[, colnames(Ga) %in% site_env_qtl$QTL]))
      EKae_s[left:right,left:right] = A.mat(as.matrix(Ga[,!colnames(Ga) %in% site_env_qtl$QTL]))
    } else {
      EKae_l[left:right,left:right] = diag(lineNum)
      EKae_s[left:right,left:right] = A
    }
  }
  
  return(list(EKae_l,EKae_s))
}

