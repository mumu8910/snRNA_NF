#########
#prepare
#########
#package
library(Seurat);cat("Seurat:",as.character(packageVersion("Seurat")),"\n")
library(clustree);cat("clustree:",as.character(packageVersion("clustree")),"\n")
library(data.table);cat("data.table:",as.character(packageVersion("data.table")),"\n")
library(future);cat("future:",as.character(packageVersion("future")),"\n")

#environment
rm(list = ls())
options(stringsAsFactors = F)

##Get the parameters
parser = argparse::ArgumentParser(description="script to Cluster scRNA data")
parser$add_argument('-d','--dim',help='dim usage')
parser$add_argument('-k','--knn',help='defines k for the k-nearest neighbor algorithm')
parser$add_argument('-r','--res_file',help='resolution usage')
parser$add_argument('-i','--input', help='input merge seurat object rds')
parser$add_argument('-o','--out',help='out directory')
parser$add_argument('-t','--thread', help='thread used')
parser$add_argument('-g','--ram', help='ram used(Gb)')
args = parser$parse_args()
##dim and cluster
dim.usage <- as.numeric(if(!is.null(args$dim)) args$dim else 30)
k.usage <- as.numeric(if(!is.null(args$knn)) args$knn else 20)
res.usage.seq <- as.vector(as.matrix(read.table(args$res_file,header=F,stringsAsFactors=F)))
res.usage <- as.numeric(if(!is.null(res.usage.seq)) res.usage.seq else 0.8)
##data
input <- args$input
#workdir
path <- args$out
if (!dir.exists(path)) {
  dir.create(path)
}
setwd(path)
##multiprocess
thread <- as.numeric(if(!is.null(args$thread)) args$thread else 5)
ram <- as.numeric(if(!is.null(args$ram)) args$ram else 2)
plan("multiprocess",workers = thread)
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = ram*1024^3)

#print parameter
cat('dim usage:',mode(dim.usage),';',dim.usage,"\n")
cat('k for the k-nearest neighbor algorithm:',mode(k.usage),';',k.usage,"\n")
cat('resolution:',mode(res.usage),';',res.usage,"\n")
cat('input seurat rds:',mode(input),';',input,"\n")
cat('out directory:',mode(path),';',path,"\n")
cat('thread:',mode(thread),';',thread,"\n")
cat('ram:',mode(ram),';',ram,"\n")

##########
#read data
##########
EC <- readRDS(input)
EC.assay <- names(EC@assays)
if ('integrated' %in% EC.assay) {
  assay.m <- 'integrated'
}else if ('SCT' %in% EC.assay) {
  assay.m <- 'SCT'
}else{
  assay.m <- 'RNA'
}
DefaultAssay(EC) <- assay.m

##################
#standard analysis
##################
EC <- RunTSNE(EC, dims = 1:dim.usage)
#cluster
EC <- FindNeighbors(EC, reduction = "pca", dims = 1:dim.usage, k.param = k.usage, annoy.metric = 'euclidean')
EC <- FindClusters(EC, resolution = res.usage, algorithm = 1, random.seed = 0)

#RNA assay all genes scale
if ('SCT' %in% EC.assay) {
  assay <- 'RNA'
}else{
  assay <- 'RNA'
}
DefaultAssay(EC) <- assay
all.genes <- rownames(EC)
EC <- ScaleData(EC, features = all.genes)
saveRDS(EC,paste(path,'/',"01_cluster.RDS",sep=""))

##############
#clustree_plot
##############
p <- clustree(EC@meta.data,prefix = paste0(assay.m,'_snn_res.'),node_size = 10)
pdf(paste(path,'/',"01_clustree.pdf",sep=""),height = length(res.usage)*1.15,width = length(unique(EC@active.ident))*0.5)
print(p)
dev.off()

#####
#save
#####
saveRDS(EC,paste(path,'/',"01_cluster.RDS",sep=""))



