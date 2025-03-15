# UMI + metadata -> SCE  
# Marker calculation using z-score normalized matrix  

# Input:  
# 1) --out: Output path; the script will generate "pb.rds" in this path, which stores the SCE object  

# Gene name conversion between species  
# 2) --B2A: (Optional) Gene ID mapping file  
# 3) --A: (Optional) Gene ID to gene symbol mapping for species A (Apis mellifera)  
# 4) --B: (Optional) Gene ID to gene symbol mapping for species B (another species)  

# 5) --col: Column name in metadata that contains cell type information  


#environment
rm(list = ls())
options(stringsAsFactors = F)

#package
library(Seurat)
library(SingleCellExperiment)
library(magrittr)
library(Matrix.utils)
library(dplyr)
library(stringr)
library(tibble)

##############
#接收parameter
##############
# 0: Accept parameters  
parser = argparse::ArgumentParser(description="script to tidy SCE")

# 1: Retrieve parameters  
#' 1) --out: Output path; a new folder named after the column will be created in this path,  
#'    and "pb.rds" will be saved in it, storing the SCE object.  
parser$add_argument('--out',help='output path')
#' 2) --B2A: (Optional) Gene ID mapping file.
parser$add_argument('--B2A',help='one-to-one orthologues gene_id between different species')
#' 3) --A: (Optional) Gene ID to gene symbol mapping for species A (Apis mellifera). 
parser$add_argument('--A',help='gene_id and gene_symbol in species A, especially Apis melifera')
#' 4) --B: (Optional) Gene ID to gene symbol mapping for species B (another species).
parser$add_argument('--B',help='gene_id and gene_symbol in species A, especially other species')
#' 5) --col: Column name in the metadata that contains cell type information.
parser$add_argument('--col',help='the column name which describe celltype')


# 2: Read input parameters  
args = parser$parse_args()

####
#OUT
####
# Required metadata columns (cell types for subsequent correlation analysis)
column.use = args$col

#workdir
path = args$out
oup_path = paste0(path,'/',column.use)
if (!dir.exists(oup_path)) {
  dir.create(oup_path)
}
setwd(oup_path)


######
#INPUT
######
#umi matrix and metadata
rds = paste0(path,'/seurat.RDS')
umi = paste0(path,'/matrix.RDS')
mdata = paste0(path,'/metadata.txt')

if (file.exists(rds)) {
  EC <- readRDS(rds)
  # Extract counts and metadata
  EC.counts <- EC@assays$RNA@counts
  EC.metadata <- EC@meta.data
  rm(EC)
  #gene
  EC.gene <- rownames(EC.counts)
}else if (file.exists(umi)&file.exists(mdata)) {
  EC.counts <- readRDS(umi)
  EC.metadata <- read.delim(mdata)
  EC.gene <- rownames(EC.counts)
  if ('CellID' %in% colnames(EC.metadata)) {
    rownames(EC.metadata) <- EC.metadata$CellID
  }else if ('barcode' %in% colnames(EC.metadata)) {
    rownames(EC.metadata) <- EC.metadata$barcode
  }else {
    return(message('Sorry~ please try another data format'))
  }
}else{
  return(message('Sorry~ please try another data format'))
}


#Homologous gene mapping  (A <- B)
B2A.file = args$B2A
A.file = args$A
B.file = args$B
B2A.check <- (if (is.null(B2A.file)) 'false' else 'true')


###################
#PROCESS
###################
#1 pseudo_group
###################
#add pseudo_group in metadta 
EC.metadata <- EC.metadata[column.use]
colnames(EC.metadata) <- c('cell_type')
EC.metadata <- mutate_all(EC.metadata,as.character)
cell_type <- unique(EC.metadata$cell_type)
EC.metadata.list <-
  lapply(X=cell_type, 
         FUN=function(x){
           data.use <- EC.metadata %>% filter(cell_type == x)
           cnt <- nrow(data.use)
           if (cnt <= 10) {
             data.use$pseudo_group <- rep(1,cnt)
           }else {
             data.use$pseudo_group <-sample(c(rep(seq(1,floor(cnt/10)),10),rep(0,cnt%%10)))
           }
           return(data.use)
         })
EC.metadata <- do.call(what = 'rbind', args = EC.metadata.list)
cell.keep <-rownames(EC.metadata)

###########
#2 ortholog gene
###########
# If the species is Apis mellifera, mapping is not needed.  
# For other species, first identify homologous genes and then find the intersection with existing genes.  
if (B2A.check == 'true') {
  B2A.df <- read.delim(B2A.file,header=FALSE) %>% select(V1,V2)
  colnames(B2A.df) <- c('gene_id_others','gene_id_apis')
  A.df <- read.delim(A.file) %>% select('gene_name','gene_id')
  B.df <- read.delim(B.file) %>% select('gene_name','gene_id')
  symbols.B <- B.df$gene_name[match(B2A.df$gene_id_others,B.df$gene_id)]
  symbols.A <- A.df$gene_name[match(B2A.df$gene_id_apis,A.df$gene_id)]
  gene.keep <- intersect(symbols.B,EC.gene)
  gene.rename <- symbols.A[match(gene.keep,symbols.B)]
}else{
  gene.keep <- EC.gene
}

###############
#3 Retain homologous genes
###############
# Adjust order! Ensure that counts and metadata correspond to the same cells.  
EC.counts <- EC.counts[,rownames(EC.metadata)]
#sum(rownames(EC.metadata) == colnames(EC.counts))
sce <- SingleCellExperiment(assays=list(counts=EC.counts),
                            colData=EC.metadata)
#colData
colData(sce) %>% 
  as.data.frame %>%
  mutate_all(as.factor) %>%
  set_rownames(colnames(sce)) %>%
  DataFrame -> colData(sce)
#head(colData(sce))

#filter
sce <- sce[gene.keep, cell.keep]
if (B2A.check == 'true') {
  gene.keep <- EC.gene
  rownames(sce) <- gene.rename
}

#########################
#4 pseudo_cell matrix
#########################
# 4.1 Pseudo-cell UMI matrix  
# Sum UMI counts for the same gene within the same pseudo-group  
# 'pb' represents pseudo-cells (rows), and genes are represented as columns 
groups <- colData(sce)[,c('cell_type','pseudo_group')]
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum")

# 4.2 CP10K: Normalization within a pseudo-cell  
# Perform library size normalization within each pseudo-cell (across rows)  
# 'pb.CP10K' represents genes (rows) and pseudo-cells (columns) 
pb.CP10K <- apply(X=pb,MARGIN=1,FUN = function(x){x*1000/sum(x)})

# 4.3 Z-score: Normalization within a gene  
# Perform z-score normalization for each gene (x - mean / standard deviation)  
# 'pb.norm' represents pseudo-cells (rows) and genes (columns)  
pb.norm <- apply(X=pb.CP10K,MARGIN=1,FUN = function(x){if (sd(x)!=0) {(x-mean(x))/sd(x)} else{x-mean(x)}})

#######################
#5 pseudo_cell metadata
#######################
pb.metadata <- colData(sce) %>%
  as.data.frame %>%
  mutate_all(as.character) %>%
  unique() %>%
  `row.names<-`(NULL)
pb.metadata <- pb.metadata %>%
  mutate(pseudo_cell=str_c(cell_type,pseudo_group,sep='_')) %>%
  column_to_rownames(var='pseudo_cell') %>%
  select(cell_type)

#######################
#6 save pseudo_cell sce
#######################
# Make sure to adjust the order!  
# Only after ensuring that the cell names in counts correspond to those in pb.metadata can the SCE object be created!
pb <- pb[rownames(pb.metadata),]
pb.sce <- SingleCellExperiment(assays=list(counts=t(pb)),
                               colData=pb.metadata)
#save
oup_file <- 'pb.RDS'
saveRDS(pb.sce, oup_file)

###################
#7 marker calculate
###################
EC <- CreateSeuratObject(counts=t(pb.norm),assay = "RNA")
EC@meta.data$orig.ident <- pb.metadata$cell_type[match(rownames(EC@meta.data),rownames(pb.metadata))]
ident <- EC@meta.data[,'orig.ident']
names(ident) <- row.names(EC@meta.data)
EC@active.ident <- factor(ident)

EC.markers <- FindAllMarkers(EC,test.use = "wilcox",only.pos = TRUE,min.pct = 0.2, logfc.threshold = 1.25)
markers <- select(EC.markers,c('cluster','gene','p_val_adj','p_val','avg_log2FC')) %>%
  filter(p_val_adj < 0.05)

#save
oup_file <- 'marker.tsv'
write.table(markers,oup_file,quote=FALSE,sep='\t',row.names=FALSE)




