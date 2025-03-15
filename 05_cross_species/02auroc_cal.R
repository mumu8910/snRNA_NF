#SCE -> auroc

# Input:  
# 1) --out: Output path; the main directory for cross_species will be created at this path.  
# 2) --A_dataset: Dataset name for species A.  
# 3) --A_col: Column name in the metadata for species A that contains cell type information.  
# 4) --B_dataset: Dataset name for species B.  
# 5) --B_col: Column name in the metadata for species B that contains cell type information.

#environment
rm(list = ls())
options(stringsAsFactors = F)

#package
library(MetaNeighbor)
library(SingleCellExperiment)
library(magrittr)
library(Matrix.utils)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
#
source("2017-08-28-runMN-US.R")
source("2017-08-28-runMN-US.pearson.R")

#####################
# Receive parameters
#####################
# 0: Accept parameters
parser = argparse::ArgumentParser(description="script to caculate AUROC")

# 1: Retrieve parameters  
#' 1) --out: Output path; a new folder named after the column will be created in this path,  
#'    and "pb.rds" will be saved in it, storing the SCE object.  
parser$add_argument('--out',help='output path')
# 2) --A_dataset: Dataset name for species A.  
parser$add_argument('--A_dataset',help='dataset in species A')
# 3) --A_col: Column name in the metadata for species A that contains cell type information.  
parser$add_argument('--A_col',help='the column name which describe celltype in species A')
# 4) --B_dataset: Dataset name for species B.  
parser$add_argument('--B_dataset',help='dataset in species B')
# 5) --B_col: Column name in the metadata for species B that contains cell type information.
parser$add_argument('--B_col',help='the column name which describe celltype in species B')


# 2: Read input parameters  
args = parser$parse_args()

# 3: Split parameters
path = args$out
dataset.A <- args$A_dataset
col.A <- args$A_col
dataset.B <- args$B_dataset
col.B <- args$B_col


####
#OUT
####
#workdir
oup_name <- paste0(dataset.A,'-',col.A,'_vs_',dataset.B,'-',col.B)
oup_path = paste0(path,'/metaNeighbor/',oup_name)
if (!dir.exists(oup_path)) {
  dir.create(oup_path)
}
setwd(oup_path)


#read SCE
#A
path.A <- paste0(path,'/data/',dataset.A,'/',col.A)

gene_file.A <- paste0(path.A,'/marker.tsv')
gene.A <- read.delim(gene_file.A)

inp_file.A <- paste0(path.A,'/pb.RDS')
sce.A <- readRDS(inp_file.A)
colnames(sce.A) <- paste0(dataset.A,'.',colnames(sce.A))

#B
path.B <- paste0(path,'/data/',dataset.B,'/',col.B)

gene_file.B <- paste0(path.B,'/marker.tsv')
gene.B <- read.delim(gene_file.B)

inp_file.B <- paste0(path.B,'/pb.RDS')
sce.B <- readRDS(inp_file.B)
colnames(sce.B) <- paste0(dataset.B,'.',colnames(sce.B))

#merge
data <- list(A=sce.A,B=sce.B)
data.merge <- mergeSCE(data)
data.merge$study_id <- str_split_fixed(colnames(data.merge),'\\.',2)[,1]
data.merge$cell_type <- paste0(data.merge$study_id,'.',data.merge$cell_type)
#check data
#dim(data.merge)
#head(data.merge)
#table(data.merge$cell_type, data.merge$study_id)

#data
# data: Genes as rows, cells as columns  
data <- counts(data.merge)

# data.scale: Genes as rows, cells as columns
data.scale <- apply(X=data,MARGIN=2,FUN = function(x){x*1000/sum(x)})

# data.zscore: Cells as rows, genes as columns
# The z-score calculation is done after merging, so gene normalization is performed across cells from both species.
data.zscore <- apply(X=data.scale,MARGIN=1,FUN = function(x){if (sd(x)!=0) {(x-mean(x))/sd(x)} else{x-mean(x)}})
# Therefore, the data needs to be transposed.
data.zscore <- t(data.zscore)

#pheno
pheno <- colData(data.merge) %>%
  as.data.frame() %>%
  rownames_to_column(var='Sample_ID') %>%
  select(Sample_ID,study_id,cell_type) %>%
  rename(Study_ID=study_id,Celltype=cell_type)

#celltypes
celltypes <- unique(pheno$Celltype)

#gene
A.gene <- gene.A %>% 
  filter(p_val_adj < 0.05,avg_log2FC>=1.25) %>% 
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) %>% as.data.frame() %>%
  select(gene)

B.gene <- gene.B %>% 
  filter(p_val_adj < 0.05,avg_log2FC>=1.25) %>% 
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) %>% as.data.frame() %>%
  select(gene)
marker.genes<-union(A.gene$gene,B.gene$gene)
#length(marker.genes)

#AUROC:zscore+maker
AUROC.matrix.s=run_MetaNeighbor_s(marker.genes, data.zscore, celltypes, pheno)
AUROC.matrix.p=run_MetaNeighbor_p(marker.genes, data.zscore, celltypes, pheno)
AUROC.data.s<-reshape2::melt(AUROC.matrix.s,value.name = "AUROCs")
AUROC.data.p<-reshape2::melt(AUROC.matrix.p,value.name = "AUROCp")
res<-merge(AUROC.data.p,AUROC.data.s,by=c("Var1","Var2"))
#save
write.table(res,"auroc.tsv",sep="\t",quote=F,row.names=F)


