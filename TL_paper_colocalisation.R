####This code is developed to calculate the posterior probability of two traits within a 150kb locus###################
##########################################load the R packages #########################################################
library(tidyverse)
library(coloc)#version 5.1.0.1

##################Colocalisation function####################################### 
#####imputs: 1)Vector with genetic variant RSIDs (i.e, 144 genetic variant RSIDs used as LTL intruments in the MR analyses)
#####imputs: 2) Output folder as string
#####imputs: 3) Output name as string
#####prior chosen were: p1 = 1e-03,p2 = 1e-04,p12 = 1e-5
#####You must have a token to use LDlinkR::LDmatrix to generate the LD matrix
#####outputs:4) Output all 5 posterior probabilities, number of credible ste of variants; query_snp, and applied colocalisation method
#####more details about output check: https://chr1swallace.github.io/coloc/

coloc_func<-function(snp_list,output_folder,output_name){
  #import sumstats
  coloc_results=data.frame()
  trait1 <-bigreadr::fread2('PATH_TO_TRAIT1_FILE')
  trait2 <-bigreadr::fread2('PATH_TO_TRAIT2_FILE')
  
  #formatting the files for coloc:
  #get the all SNPs +/-75kb centered on the queried SNPs	
  for(i in snp_list){
    print(i)
    trait1 %>% filter(rsid==i) %>% mutate(Start=pos-75000, End=pos+75000, Start=ifelse(Start<0,0,Start)) %>%
      dplyr::select(rsid,chr,pos,Start,End) -> keep_regions #create file with regions to be excluded from sumstats
    keep_regions=as.data.frame(keep_regions)
    exc.rows <- NULL
    id <- which(trait1$chr == keep_regions$chr[1] & trait1$pos >= keep_regions$Start[1] & trait1$pos <= keep_regions$End[1])
    exc.rows <- c(exc.rows, id)
    filtered_sumstat <- trait1[exc.rows, ]
    filtered_sumstat %>% dplyr::select(rsid,beta,beta_se,a1,p,a1_freq) -> temp0
    colnames(temp0)=c('rsid','beta_1','beta_se_1','a1_1','p_1','a1_freq')
    
    #filter the sumstats and get the a1 and effect related to trait1
    temp1=trait2
    temp1 %>% dplyr::select(rsid,beta,beta_se,a1,p) -> temp1
    temp2=left_join(temp0,temp1,by='rsid')
    temp2 %>% mutate(beta=ifelse(a1==a1_1,beta,beta*(-1))) %>% dplyr::select(rsid,beta,beta_se,p)-> temp2
    colnames(temp2)=c('rsid','beta_2','beta_se_2','p_2') 
    temp0=left_join(temp0,temp2,by='rsid')
    
    temp0=na.omit(temp0) #remove rows with any NA from data.frame
    temp0 %>% mutate(z_1=beta_1/beta_se_1,z_2=beta_1/beta_se_2)-> temp0
    
    #ld matrix used by conditioning and mask methods
    skip_to_next <- FALSE
    ldmatrix<-LDlinkR::LDmatrix(temp0$rsid, pop = "EUR", r2d = "r2", token = 'Type of token', file = FALSE)
    tryCatch(ldmatrix %>% dplyr::select(i)-> a, error = function(e) { skip_to_next <<- TRUE})
    if(skip_to_next) { next }
    
    #break in the code to calculate the list of snps to filter the ldmatrix
    ldmatrix %>% mutate(na = rowSums(is.na(.))) %>% filter(na==nrow(ldmatrix)) %>% dplyr::select(RS_number,na)-> to_exclude
    snp<-ldmatrix$RS_number
    ldmatrix=ldmatrix[,2:ncol(ldmatrix)]
    row.names(ldmatrix)=snp
    ldmatrix=ldmatrix[! rownames(ldmatrix) %in% to_exclude$RS_number,]
    ldmatrix=ldmatrix[,rownames(ldmatrix)]
    temp0 %>% filter(rsid %in% colnames(ldmatrix)) -> temp0
    corr=ldmatrix[temp0$rsid,]
    corr=corr[,temp0$rsid]
    corr %>% mutate(across(everything(), sqrt)) -> corr
    
    ## Create the lists required for coloc
    D1=list(beta=temp0$beta_1,
            varbeta=temp0$beta_se_1^2,
            snp=temp0$rsid,
            pvalues=temp0$p_1,
            LD = as.matrix(corr),
            MAF=temp0$a1_freq,
            sdY=1,#
            N=464716, #(i.e, number of participant of LTL GWAS) 
            type="quant")
    
    D2=list(beta=temp0$beta_2,
            varbeta=temp0$beta_se_2^2,
            snp=temp0$rsid,
            pvalues=temp0$p_2,
            LD = as.matrix(corr),
            MAF=temp0$a1_freq,
            s=0.1968, # This is the case proportion. Make sure to specify if using a case-control GWAS (get this information from the paper
            N=61073, #number of participants
            type="cc")
    
    ########running different colocalisation methods
    
    #mask iterative
    res=coloc.signals(D1,D2,method = c("mask"),LD=NULL,mode = c("iterative"),p1 = 1e-03,p2 = 1e-04,p12 = 1e-5,
                      maxhits = 3,r2thr = 0.01,pthr = 1e-06)
    res=as.data.frame(res$summary)
    res <- mutate(res,query_snp=i,method='mask_iterative',avg_z1=temp0$avg_z1,avg_z2=temp0$avg_z2)
    coloc_results=rbind(coloc_results,res)
    
    #mask allbutone
    res=coloc.signals(D1,D2,method = c("mask"),LD=NULL,mode = c("allbutone"),p1 = 1e-03,p2 = 1e-04,p12 = 1e-5,
                      maxhits = 3,r2thr = 0.01,pthr = 1e-06)
    res=as.data.frame(res$summary)
    res <- mutate(res,query_snp=i,method='mask_allbutone',avg_z1=temp0$avg_z1,avg_z2=temp0$avg_z2)
    coloc_results=rbind(coloc_results,res)
    
    
    #conditional iterative
    res=coloc.signals(D1,D2,method = c('cond'),mode = c("iterative"),p1 = 1e-03,p2 = 1e-04,p12 = 1e-5,
                      maxhits = 3,r2thr = 0.01,pthr = 1e-06)
    
    res=as.data.frame(res$summary)
    res <- mutate(res,query_snp=i,method='cond_iterative',avg_z1=temp0$avg_z1,avg_z2=temp0$avg_z2)
    coloc_results=rbind(coloc_results,res)
    
    #default: single (original coloc)
    res=coloc.signals(D1,D2,method = c('single'),LD=NULL,mode = NULL,p1 = 1e-03,p2 = 1e-04,p12 = 1e-5,
                      maxhits = 3,r2thr = 0.01,pthr = 1e-06)
    res=as.data.frame(res$summary)
    res <- mutate(res,query_snp=i,method='single',avg_z1=temp0$avg_z1,avg_z2=temp0$avg_z2)
    coloc_results=rbind(coloc_results,res)
  }
  write_csv(coloc_results,paste0(output_folder,output_name,'.csv'))
}




