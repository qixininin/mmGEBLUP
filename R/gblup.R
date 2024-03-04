gblup <- function(geno_data, pheno_data, Ka)
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
