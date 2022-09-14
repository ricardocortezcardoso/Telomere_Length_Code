##This code was developed to run hyprcoloc

#installing hyprcoloc and plotting functions
#devtools::install_github("cnfoley/hyprcoloc", build_opts = c("--resave-data", "--no-manual"), build_vignettes = TRUE)

#imput: snps_list is a string vector with the RSIDs of interest to be queried
#imput: file_to_GWAS is a string vector with the path to the GWAS summary statitics. Order of traits is the same as the order of file_to_GWAS
#GWAS summary statistics should have the follwing header: rsid (SNP RSID),chr (chromossome),pos (genome position),beta (beta estimate of the original GWAS or log OR if binary trait),beta_se (standard error of beta estimate),a1 (effect allele in the original GWAS)
#more details about output check https://github.com/jrs95/hyprcoloc

#load R packages
library(tidyverse)
library(hyprcoloc)

#######Run hyprcoloc function

hyprcoloc_function<-function(snp_list,file_to_GWAS,output_folder,output_name){
  
  #import sumstats
  for (i in length(file_to_GWAS)){
      eval(parse(text = paste0("trait",i,"= bigreadr::fread2(","file_to_GWAS[",i,"])")))}
  
  n_traits=length(file_to_GWAS)
  
  #formatting the files for coloc:
  #get the SNPs +/- 75kb centered on the queried SNPs	
  res1=data.frame()
  for(i in snp_list){
    trait1 %>% filter(rsid==i) %>% mutate(Start=pos-75000, End=pos+75000, Start=ifelse(Start<0,0,Start)) %>%
      dplyr::select(rsid,chr,pos,Start,End) -> keep_regions #create file with regions to be excluded from sumstats
    keep_regions=as.data.frame(keep_regions)
    exc.rows <- NULL
    id <- which(trait1$chr == keep_regions$chr[1] & trait1$pos >= keep_regions$Start[1] & trait1$pos <= keep_regions$End[1])
    exc.rows <- c(exc.rows, id)
    filtered_sumstat <- trait1[exc.rows, ]
    filtered_sumstat %>% dplyr::select(rsid,chr,pos,beta,beta_se,a1) -> temp0
    colnames(temp0)=c('rsid','chr','pos','beta_1','beta_se_1','a1_1')
    
    #filter sumstats and get the a1 and effect related to trait1
    for(j in seq(2,n_traits)){
      temp1=get(paste0('trait',j))
      temp1 %>% dplyr::select(rsid,beta,beta_se,a1) -> temp1
      temp2=left_join(temp0,temp1,by='rsid')
      temp2 %>% mutate(beta=ifelse(a1==a1_1,beta,beta*(-1))) %>% dplyr::select(rsid,beta,beta_se)-> temp2
      colnames(temp2)=c('rsid',paste0('beta_',j),paste0('beta_se_',j)) 
      temp0=left_join(temp0,temp2,by='rsid')
    }
    
    
    temp0=na.omit(temp0) #remove rows with any NA from data.frame
    
    #matrix1: betas
    temp0 %>% dplyr::select(paste0('beta_',seq(1,n_traits))) -> betas
    colnames(betas)=c(paste0('T',seq(1,n_traits)))
    betas=as.matrix(betas)
    row.names(betas)=temp0$rsid
    
    #matrix2: ses (standard error)
    temp0 %>% dplyr::select(paste0('beta_se_',seq(1,n_traits))) -> ses
    colnames(ses)=c(paste0('T',seq(1,n_traits)))
    ses=as.matrix(ses)
    row.names(ses)=temp0$rsid
    
    #hyprcoloc
    traits <- c(paste0('T',seq(1,n_traits)));
    rsid <- rownames(betas);
    skip_to_next <- FALSE
    tryCatch(res <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid,snpscores = F),error = function(e) { skip_to_next <<- TRUE})
    if(skip_to_next) { next }
    res=as.data.frame(res$results)
    res$query_snp=i
    res1=rbind(res1,res)
  }
  write_csv(res1,paste0(output_folder,output_name,'.csv'))
}






