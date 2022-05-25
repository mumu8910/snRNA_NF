##################
#rpca+LogNormalize
##################
### package
library(Seurat);cat("Seurat:",as.character(packageVersion("Seurat")),"\n")
library(tidyr);cat("tidyr:",as.character(packageVersion("tidyr")),"\n")
library(future);cat("future:",as.character(packageVersion("future")),"\n")

### Get the parameters
parser = argparse::ArgumentParser(description='rpca+LogNormalize')
parser$add_argument('-f','--nfeatures', help='the number of variable features')
parser$add_argument('-o','--out', help='output directory')
parser$add_argument('-t','--thread', help='thread used')
parser$add_argument('-g','--ram', help='ram used(Gb)')
args = parser$parse_args()

nfeatures <- as.numeric(if(!is.null(args$nfeatures)) args$nfeatures else 2000)
out <- args$out
thread <- as.numeric(if(!is.null(args$thread)) args$thread else 10)
ram <- as.numeric(if(!is.null(args$ram)) args$ram else 50)
dim.usage <- 30
pc.usage <- 50
seed.usage <- 0
maxdim.usage <- "2L"

plan("multiprocess",workers = thread)
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = ram*1024^3)

###FindIntegrationAnchors
reduction.method <- 'rpca' #c("cca", "rpca", "rlsi")
norm.method <- 'LogNormalize' #c("LogNormalize", "SCT")
### input RDS files
input <- paste0(out,'/rds.list')
files <- as.vector(as.matrix(read.table(input,header=F,stringsAsFactors=F)))

### read data
objectlist <- list()
for(i in 1:length(files)){
  objectlist[[i]] <- readRDS(files[i])
}
cat("The number of RDS inputed:",length(files),"\n")

print('Objectlist creat begin!')
objectlist <- lapply(X = objectlist, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures)
  x <- ScaleData(x)
  x <- RunPCA(x)
  x <- RunUMAP(x,dims = 1:dim.usage)
})
saveRDS(objectlist,paste0(out,"/objectlist.RDS"))
#objectlist <- readRDS(paste0(out,"/objectlist.RDS"))
print('Objectlist done!')
cat('\n')

### Seurat Integrate
print('Integrate begin!')
object.anchors <- FindIntegrationAnchors(object.list = objectlist,
                                         anchor.features = nfeatures,
                                         dims = 1:dim.usage,
                                         reduction=reduction.method,
                                         normalization.method=norm.method)
object.combined <- IntegrateData(anchorset = object.anchors, dims = 1:dim.usage)
saveRDS(object.combined,paste0(out,"/merge.RDS"))
print('Integrate done!')
cat('\n')

print('Standard analysis begin!')
DefaultAssay(object.combined) <- "integrated"
object.combined <- ScaleData(object.combined, verbose = FALSE)
object.combined <- RunPCA(object.combined, npcs = pc.usage, seed.use = seed.usage,verbose = FALSE)
object.combined <- RunUMAP(object.combined, reduction = "pca", dims = 1:dim.usage, seed.use = seed.usage, max.dim = maxdim.usage)
#DefaultAssay(object.combined) <- "RNA"
saveRDS(object.combined,paste0(out,"/merge_final.RDS"))
print('Standard analysis done!')