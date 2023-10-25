library("optparse")
option_list <- list(
  make_option(c("-r", "--row"), action="store_true", default=TRUE,
              help="r0w",type='numeric')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

row = opt$row

library(susieR)
library(MASS)
require(bigsnpr)
require(MOSTWAS)

###set parameter we want to vary
p.causal = c(0.01,0.05,.1)
eqtl.h2.loc = c(0.05,.1,.25)
eqtl_ratio = c(.1,.5,1,1.5,2)
params = expand.grid(
  p.causal,
  eqtl.h2.loc,
  eqtl_ratio)
colnames(params) = c('p.causal','eqtl.h2','ratio')


p.causal = params$p.causal[row]
eqtl.h2 = params$eqtl.h2[row]
ratio =params$ratio[row]


for (chr in 1:22) {
  

  


  #read and clean the UK biobank and 1000G data
  snps = snp_attach(snp_readBed2(paste0('/u/scratch/y/yhc1998/ukb_1k_filtersnp/chr_',chr,'.bed'),
                                 backingfile = tempfile()))
  
  
  snps_afr = snp_attach(snp_readBed2(paste0('/u/home/y/yhc1998/project-pasaniuc/yhc/AFR/chr',chr,'_new.bed'),
                                     backingfile = tempfile()))
  
  #let snps lists are same in each population
  eur_geno<-snps$map$physical.pos
  afr_geno<-snps_afr$map$physical.pos
  rm_list<-afr_geno[(afr_geno %in% eur_geno)]
  
  snps_afr = snp_attach(subset(snps_afr,
                               ind.col = which(snps_afr$map$physical.pos %in% rm_list)))
  


adjr_europe_adjR<-data.frame(matrix(ncol = 1000, nrow = 1))
adjr_afr_adjR<-data.frame(matrix(ncol = 1000, nrow = 1))

for (i in 1:1000){
  
  ### Simulate cis expression
  Z_qtl= as.data.frame(snps$genotypes[])
  
  # getting mean of each column using apply() 
  all_column_mean<- apply(Z_qtl, 2, mean, na.rm=TRUE)
  
  # imputing mean value with NA 
  for(k in colnames(Z_qtl))
    Z_qtl[,k][is.na(Z_qtl[,k])] <- all_column_mean[k]
  
  
  
  
  ## spliting training and testing set for ukbiobank
  train_ind <- sample(seq_len(nrow(Z_qtl)), size = 500)
  
  Z_qtl_train <-  as.matrix(Z_qtl[train_ind, ])
  Z_qtl_test <-  as.matrix(Z_qtl[-train_ind, ])
  
  #sample individual
  n = nrow(Z_qtl_train)
  p = ncol(Z_qtl_train)
  
  # GRM AND LD
  b = simBeta(p.causal,
              eqtl.h2,
              n_snps = ncol(Z_qtl_train))
  
  # sim training gene expression dataset
  gexpr_train = simTrait(Z_qtl_train%*% b, eqtl.h2)
  

  # sim testing gene expression dataset
  gexpr_test= simTrait(Z_qtl_test%*% b, eqtl.h2)
  
  
  ###African Testing Set (n=500)###
  Z_qtl_afr= as.data.frame(snps_afr$genotypes[])
  
  # getting mean of each column using apply() 
  all_column_mean_afr<- apply(Z_qtl_afr, 2, mean, na.rm=TRUE)
  
  # imputing mean value with NA 
  for(j in colnames(Z_qtl_afr))
    Z_qtl_afr[,j][is.na(Z_qtl_afr[,j])] <- all_column_mean_afr[j]
  
  
  ## spliting training and testing set for ukbiobank
  smp_size <- 500
  test_ind <- sample(seq_len(nrow(Z_qtl_afr)), size = smp_size)
  Z_qtl_afr_test <-  as.matrix(Z_qtl_afr[test_ind, ])
  
  # sim testing african gene expression dataset
  gexpr_afr_test = simTrait(Z_qtl_afr_test%*% b, eqtl.h2 * ratio)
  
  
  # run SuSiE model
  model_susie_European = susie(Matrix::Matrix(Z_qtl_train,sparse=T),
                               gexpr_train,
                               estimate_residual_variance =
                                 TRUE,
                               estimate_prior_variance =
                                 FALSE,
                               scaled_prior_variance = 0.1,
                               verbose = F)

  
  ### Input European Testing Set###
  
  df1<-predict(model_susie_European, newx =Matrix::Matrix(Z_qtl_test,sparse=T))
  df_eur_test<-data.frame(trainGE = gexpr_test,est_GE = df1[,1])
  df_eur_lm_test<-lm(df_eur_test$trainGE~df_eur_test$est_GE , data=df_eur_test)
  adjr_europe_test=summary(df_eur_lm_test)$adj.r.squared
  adjr_europe_adjR[1,i]=adjr_europe_test
  
  
  
  ### Input African Testing Set###
  df2<-predict(model_susie_European, newx =Matrix::Matrix(Z_qtl_afr_test,sparse=T))
  df_afr_test<-data.frame(trainGE = gexpr_afr_test,est_GE = df2[,1])
  df_afr_lm<-lm(trainGE ~ est_GE, data=df_afr_test)
  adjr_afr_test=summary(df_afr_lm)$adj.r.squared
  adjr_afr_adjR[1,i]=adjr_afr_test
  
  
}
diff_eur_afr = adjr_europe_adjR - adjr_afr_adjR
out.file1 <- paste("/u/scratch/y/yhc1998/h2_ratio/r2_chr",chr,"_", row, ".txt", sep="")

write.table(diff_eur_afr, file=out.file1)
}
