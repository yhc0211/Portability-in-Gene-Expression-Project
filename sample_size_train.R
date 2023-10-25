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
library(data.table)
#n.eqtl = c(500,1000)
p.causal = c(0.01)
eqtl.h2.loc = c(0.05)
train_sample_size = c(250, 500, 1000, 10000)
params = expand.grid(
  p.causal,
  eqtl.h2.loc,
  train_sample_size)
colnames(params) = c('p.causal','eqtl.h2','train_sample_size')



p.causal = params$p.causal[row]
eqtl.h2 = params$eqtl.h2[row]
sample.size = params$train_sample_size[row]

adjr_europe_adjR<-data.frame(matrix(ncol = 1000, nrow = 1))
adjr_afr_adjR<-data.frame(matrix(ncol = 1000, nrow = 1))



for (chr in 1:22){
  
  
  
  snps = snp_attach(snp_readBed2(paste0('/u/scratch/y/yhc1998/sample_size/ukb_10k_filtersnp/chr_',chr,'.bed'),
                                 backingfile = tempfile()))
  
  
  snps_afr = snp_attach(snp_readBed2(paste0('/u/home/y/yhc1998/project-pasaniuc/yhc/AFR/chr',chr,'_new.bed'),
                                     backingfile = tempfile()))
  
  
  eur_geno<-snps$map$physical.pos
  afr_geno<-snps_afr$map$physical.pos
  rm_list<-afr_geno[(afr_geno %in% eur_geno)]
  
  snps_afr = snp_attach(subset(snps_afr,
                               ind.col = which(snps_afr$map$physical.pos %in% rm_list)))
  
  
  
  
  
  for (i in 1:1000){
    
    ### Simulate cis expression
    Z_qtl= as.data.frame(snps$genotypes[])
    
    # getting mean of each column using apply() 
    all_column_mean<- apply(Z_qtl, 2, mean, na.rm=TRUE)
    
    # imputing mean value with NA 
    for(k in colnames(Z_qtl))
      Z_qtl[,k][is.na(Z_qtl[,k])] <- all_column_mean[k]
    
    
    
    
    ## spliting training and testing set for ukbiobank
    train_ind <- sample(seq_len(nrow(Z_qtl)), size = sample.size)
    
    Z_qtl_train <-  as.matrix(Z_qtl[train_ind, ])
    
    Z_qtl_test <-  as.matrix(Z_qtl[-train_ind, ])
    
    test_ind <- sample(seq_len(nrow(Z_qtl_test)), size = 500)
    Z_qtl_test <-  as.matrix(Z_qtl_test[test_ind, ])
    #sample individual
    n = nrow(Z_qtl_train)
    p = ncol(Z_qtl_train)
    
    # GRM AND LD
    b = simBeta(p.causal,
                eqtl.h2,
                n_snps = ncol(Z_qtl_train))
    
    # sim gene expression
    gexpr_train = simTrait(Z_qtl_train%*% b, eqtl.h2)
    
    
    
    ###European Testing Set (n=500)###
    
    # sim gene expression
    gexpr_test= simTrait(Z_qtl_test%*% b, eqtl.h2)
    
    
    ###African Testing Set (n=500)###
    
    Z_qtl_afr= as.data.frame(snps_afr$genotypes[])
    
    # getting mean of each column using apply() 
    all_column_mean_afr<- apply(Z_qtl_afr, 2, mean, na.rm=TRUE)
    
    # imputing mean value with NA 
    for(j in colnames(Z_qtl_afr))
      Z_qtl_afr[,j][is.na(Z_qtl_afr[,j])] <- all_column_mean_afr[j]
    
    
    
    test_ind_afr <- sample(seq_len(nrow(Z_qtl_afr)), size = 500)
    
    Z_qtl_afr_test <-  as.matrix(Z_qtl_afr[test_ind_afr, ])
    
    # sim gene expression
    gexpr_afr_test = simTrait(Z_qtl_afr_test%*% b, eqtl.h2)
    
    
    
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
  
  out.file1 <- paste("/u/scratch/y/yhc1998/sample_size/training_r2_result/r2_EUR_sample_chr",chr,"_", row, ".txt", sep="")
  out.file2 <- paste("/u/scratch/y/yhc1998/sample_size/training_r2_result/r2_AFR_sample_chr",chr,"_", row, ".txt", sep="")
  write.table(adjr_europe_adjR, file=out.file1)
  write.table(adjr_afr_adjR, file=out.file2)
}
