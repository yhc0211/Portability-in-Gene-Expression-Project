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
runLD = function(row){  
  LD_norm_result<-data.frame(matrix(ncol = 6, nrow = 1))
  LD_norm_result<-setNames(LD_norm_result,c("chromosome_name","start_position","end_position","norm 1","norm 2","norm I"))
  
  
  results<-fread("/u/home/y/yhc1998/project-pasaniuc/yhc/ukb_simulation/ukb_LD/genes_infor.txt")
  chromosome.name<-results$chromosome_name[row]
  start.position<-results$start_position[row]
  end.position<-results$end_position[row]


  snps = snp_attach(snp_readBed2(paste0('/u/scratch/y/yhc1998/ukb_LD_500/chr_list',row,'.bed'),
                                 backingfile = tempfile()))
  snps_afr = snp_attach(snp_readBed2(paste0('/u/project/pasaniuc/pasaniucdata/DATA/plink_1KG_bypop/AFR/chr',chromosome.name,".bed"),
                                     backingfile = tempfile()))
  
  eur_geno<-snps$map$physical.pos
  
  
  snps_afr = snp_attach(subset(snps_afr,
                               ind.col = which(snps_afr$map$physical.pos %in% eur_geno)))
  
  LD_eur<-as.matrix(snp_cor(snps$genotypes))
  LD_afr<-as.matrix(snp_cor(snps_afr$genotypes))
  
  LD<- LD_eur - LD_afr
  
  norm1<-norm(LD,type = "O")
  norm2<-norm(LD,type = "2")
  normI<-norm(LD,type = "I")
  
  
  LD_norm_result[1,1]<-chromosome.name
  LD_norm_result[1,2]<-start.position
  LD_norm_result[1,3]<-end.position
  LD_norm_result[1,4]<-norm1
  LD_norm_result[1,5]<-norm2
  LD_norm_result[1,6]<-normI
  
  out.file1 <- paste("/u/scratch/y/yhc1998/ukb_LD_500_norm/LD_norm_", row, ".txt", sep="")
  write.table(LD_norm_result, file=out.file1,row.names=F)
}

require(pbmcapply)
tc_assoc = function(s){
  tryCatch(runLD(s),
           error = function(e) {
             print(paste0('error with ',s))
           }
  )
}


lapply(1:500,tc_assoc)