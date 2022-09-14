###########This is code is developed to perform Principal Component Analysis (PCA) based on RNA-sequencing data from TCGA
# Load R packages
library(tidyverse)
library(TCGAbiolinks)
library(knitr)
library(factoextra)
library(fpc)
library(smacof)
library(gridExtra)
library(MASS)
library(BBmisc)
library(stats)

#Query platform Illumina HiSeq with a list of barcodes
#For TL manuscript, we included primary lung adenocarcinoma tumors with both RNA-seq and germline data avaibale ('BARCODE_LIST')

#Query the samples of interest for a given platform in TCGA dataset using TCGAbiolinks 

query <- GDCquery(project = 'TCGA_PROJECT', #you can query any TCGA cohort of interest (i.e, TCGA-LUAD for the current study)
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  barcode= 'BARCODE_LIST',
                  legacy = TRUE)

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
TCGAbiolinks::GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
rnaseq_assay <- GDCprepare(query)

# normalization of genes
dataNorm <- TCGAanalyze_Normalization(tabDF = rnaseq_assay, geneInfo =  geneInfo)

# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

#formating expression data: Gene symbols as cols and barcode as rows
#log+1 transform RNA-seq data
dataFilt[dataFilt==0]<-1
dataFilt=log(dataFilt)
dataFilt_T<-t(dataFilt) #transpose matrix
a<-colnames(dataFilt)
b<-row.names(dataFilt)
dataFilt_T=as.data.frame(dataFilt_T)#convert matrix to dataframe format
colnames(dataFilt_T)=b

#detection of outliers in our analyses using Eucledian distance in two dimensions (mds1 and mds2)
distance <- dist(as.matrix(dataFilt_T[,1:nrow(dataFilt_T)])) 
mds<-cmdscale(distance, k=2)

#create a dataframe to check for outliers
temp_data=data.frame(barcode=rownames(dataFilt_T),mds1=mds=mds[,1],mds2=mds[,2])

#flagging samples detected as outliers (outside 1.5 times the interquartile range above the upper quartile and bellow the lower quartile)
#in mds1 or mds2 if any detected. In our analyses, only one sample was dropped.
temp_data<- mutate(temp_data,
  outlier=ifelse(mds1 %in% boxplot.stats(mds[,1])$out | mds2 %in% boxplot.stats(mds[,2])$out,'Yes','No'))

#Running PCA
data_matrix_1=as.matrix(dataFilt_T[,1:nrow(dataFilt_T)])#convert dataframe to matrix
data_matrix_1_PCA = prcomp(data_matrix_1,scale=T) #run PCA, scale gene expression to the mean before run
PCA_values<- as.data.frame(data_matrix_1_PCA$x) #extract eignvalues

#plot variance explained by the first 5 principal components (PC)
fviz_eig(data_matrix_1_PCA,barfill = 'steelblue3',ggtheme= theme_classic(base_size = 22) ,addlabels = T,
         hjust = 0.5,main = '',ncp = 5,xlab = 'Principal Components')




