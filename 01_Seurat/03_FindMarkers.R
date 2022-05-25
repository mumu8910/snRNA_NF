#########
#prepare
#########
#package
library(Seurat);cat("Seurat:",as.character(packageVersion("Seurat")),"\n")
library(data.table);cat("data.table:",as.character(packageVersion("data.table")),"\n")
library(tidyr);cat("tidyr:",as.character(packageVersion("tidyr")),"\n")
library(dplyr);cat("dplyr:",as.character(packageVersion("dplyr")),"\n")
library(patchwork);cat("patchwork:",as.character(packageVersion("patchwork")),"\n")
library(ggpubr);cat("ggpubr:",as.character(packageVersion("ggpubr")),"\n")
library(future);cat("future:",as.character(packageVersion("future")),"\n")

#environment
rm(list = ls())
options(stringsAsFactors = F)

#parameter
##Get the parameters
parser = argparse::ArgumentParser(description="script to find DEG")
parser$add_argument('-r','--res_usage_i',help='resolution usage')
parser$add_argument('-o','--out',help='out directory')
parser$add_argument('-t','--thread', help='thread used')
parser$add_argument('-g','--ram', help='ram used(Gb)')
args = parser$parse_args()
##resolution
res.usage.i <- as.numeric(if(!is.null(args$res_usage_i)) args$res_usage_i else 0.8)
##workdir
path <- args$out
path.res <- paste0(path,'/','res_',res.usage.i)
if (!dir.exists(path.res)) {
  dir.create(path.res)
}
setwd(path.res)
##multiprocess
thread <- as.numeric(if(!is.null(args$thread)) args$thread else 5)
ram <- as.numeric(if(!is.null(args$ram)) args$ram else 2)
plan("multiprocess",workers = thread)
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = ram*1024^3)
##input RDS
EC <- readRDS(paste0(path,'/',"01_cluster.RDS"))
DefaultAssay(EC) <- 'RNA'

#print parameter
cat('resolution:',mode(res.usage.i),';',res.usage.i,"\n")
cat('out directory:',mode(path.res),';',path.res,"\n")
cat('thread:',mode(thread),';',thread,"\n")
cat('ram:',mode(ram),';',ram,"\n")

##########################
#diff analysis:find marker
##########################
EC.assay <- names(EC@assays)
if ('integrated' %in% EC.assay) {
  assay.m <- 'integrated'
}else if ('SCT' %in% EC.assay) {
  assay.m <- 'SCT'
}else{
  assay.m <- 'RNA'
}
res <- paste0(assay.m,'_snn_res.',res.usage.i)
ident <- EC@meta.data[,res]
names(ident) <- row.names(EC@meta.data)
EC@active.ident <- ident
#Find Markers
EC.markers <- FindAllMarkers(EC, test.use = "wilcox",only.pos = TRUE,min.pct = 0.25, logfc.threshold = 0.1)
#output
markers <- select(EC.markers,c('cluster','gene','p_val_adj','p_val','avg_log2FC'))
write.csv(markers,'./02_diff_marker.csv',row.names = F,quote = F)



