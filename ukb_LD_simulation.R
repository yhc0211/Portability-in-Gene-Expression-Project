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
library(mcreplicate)


infor<-fread("/u/home/y/yhc1998/project-pasaniuc/yhc/ukb_simulation/ukb_LD/genes_infor.txt")
chromosome.name<-infor$chromosome_name[row]


###European

snps = snp_attach(snp_readBed2(paste0('/u/scratch/y/yhc1998/ukb_LD_500/chr_list',row,'.bed'),
                               backingfile = tempfile()))
snps_afr = snp_attach(snp_readBed2(paste0('/u/project/pasaniuc/pasaniucdata/DATA/plink_1KG_bypop/AFR/chr',chromosome.name,".bed"),
                                   backingfile = tempfile()))

eur_geno<-snps$map$physical.pos


snps_afr = snp_attach(subset(snps_afr,
                             ind.col = which(snps_afr$map$physical.pos %in% eur_geno)))


p.causal = c(0.01,0.05,.1)
eqtl.h2.loc = c(0.05,.1,.25)
params = expand.grid(
  p.causal,
  eqtl.h2.loc)
colnames(params) = c('p.causal','eqtl.h2')



for (parameter in 1:9) {
  
  
  p.causal = params$p.causal[parameter]
  eqtl.h2 = params$eqtl.h2[parameter]
  
  
  imputeTrainSets<-function(snps,p.causal,eqtl.h2,snps_afr){
    
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
    
    # sim gene expression
    gexpr_train = simTrait(Z_qtl_train%*% b, eqtl.h2)
    
    

    
    # sim gene expression
    gexpr_test= simTrait(Z_qtl_test%*% b, eqtl.h2)
    
    
    
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
    
    
    
    
    
    df1<-predict(model_susie_European, newx =Matrix::Matrix(Z_qtl_test,sparse=T))
    df_eur_test<-data.frame(trainGE = gexpr_test,est_GE = df1[,1])
    df_eur_lm_test<-lm(df_eur_test$trainGE~df_eur_test$est_GE , data=df_eur_test)
    adjr_europe_test=summary(df_eur_lm_test)$adj.r.squared

    
    
    
    ### Input African Testing Set###
    df2<-predict(model_susie_European, newx =Matrix::Matrix(Z_qtl_afr_test,sparse=T))
    df_afr_test<-data.frame(trainGE = gexpr_afr_test,est_GE = df2[,1])
    df_afr_lm<-lm(trainGE ~ est_GE, data=df_afr_test)
    adjr_afr_test=summary(df_afr_lm)$adj.r.squared

    
    result<-c(europe=adjr_europe_test,african=adjr_afr_test)
    return(result)
  }
  result<-as.data.frame(mc_replicate(1000,imputeTrainSets(snps,p.causal,eqtl.h2,snps_afr)),row.names = NULL)
  result<-t(result)
  
  
  out.file1 <- paste("/u/scratch/y/yhc1998/ukb_500_simulation/r2_", row,'_',parameter, ".txt", sep="")
  write.table(result, file=out.file1)
}  
