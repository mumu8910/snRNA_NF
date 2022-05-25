# data handling
library(scater)
library(Seurat)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(data.table)
# visualzation
library(ComplexHeatmap)
library(ggplot2)
library(pheatmap)
library(scales)
# multiprocess
library(future)

#environment
rm(list = ls())
options(stringsAsFactors = F)

#parameter
##Get the parameters
parser = argparse::ArgumentParser(description="script to find DEG")
parser$add_argument('-c','--column',help='column to find DEG')
parser$add_argument('-i','--input',help='input rds')
parser$add_argument('-o','--path',help='output path')
parser$add_argument('-t','--thread', help='thread used')
parser$add_argument('-g','--ram', help='ram used(Gb)')
args = parser$parse_args()

# multiprocess
#thread <- 1; ram <- 0.5
thread <- as.numeric(if(!is.null(args$thread)) args$thread else 5)
ram <- as.numeric(if(!is.null(args$ram)) args$ram else 2)
plan("multiprocess",workers = thread)
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = ram*1024^3)

#result directory
path <- args$path
path.result <- paste0(path,'/01prepare_sce')
if (!dir.exists(path.result)) {
  dir.create(path.result)
}
setwd(path.result)

#input RDS
input <- args$input
EC <- readRDS(input)

#chose column
column <- args$column
check <- is.na(as.numeric(column))
EC.assay <- names(EC@assays)
if ('integrated' %in% EC.assay) {
  assay.m <- 'integrated'
}else if ('SCT' %in% EC.assay) {
  assay.m <- 'SCT'
}else{
  assay.m <- 'RNA'
}
column <- (if(check) column else paste0(assay.m,'_snn_res.',column))

EC.split <- EC
#counts and metadata
EC.split.counts <- EC.split@assays$RNA@counts
EC.split.metadata <- EC.split@meta.data
#select metadata columns
column.use <- c('orig.ident','split',column)
EC.split.metadata <- EC.split.metadata[column.use]
colnames(EC.split.metadata) <- c('orig.ident','split','cluster')
EC.split.metadata <- mutate_all(EC.split.metadata,as.character)

#creat SingCellExperiment object
sce <- SingleCellExperiment(assays=list(counts=EC.split.counts),
                            colData=EC.split.metadata)
#tidy metadata
colData(sce) %>% 
  as.data.frame %>% 
  transmute(
    group_id = orig.ident, 
    sample_id = split,
    cluster_id = cluster) %>%
  mutate_all(as.factor) %>% 
  set_rownames(colnames(sce)) %>% 
  DataFrame -> colData(sce)
head(colData(sce))

#cell type
nk <- length(kids <- set_names(levels(sce$cluster_id)))
#sample
ns <- length(sids <- set_names(levels(sce$sample_id)))

m <- match(sids, sce$sample_id)
n_cells <- as.numeric(table(sce$sample_id))
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"cluster_id")

# remove undetected genes
sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

sce <- addPerCellQC(sce)
sce <- addPerFeatureQC(sce)

# get cells w/ few/many detected genes
sce$is_outlier <- isOutlier(
  metric = sce$detected,
  nmads = 2, type = "both", log = TRUE)

# remove outlier cells
sce <- sce[, !sce$is_outlier]
dim(sce)

# compute normalized expression values for visualization
#sizeFactor, saved in matadata
sizeFactors(sce) <- librarySizeFactors(sce)
#logcounts, saved in assays
sce <- logNormCounts(sce)
#assayNames(sce)
#"counts"    "logcounts"
#range(logcounts(sce))

save(sce,
     nk,kids,
     ns,sids,
     ei,
     file = "01prepare_sce.Rdata")
