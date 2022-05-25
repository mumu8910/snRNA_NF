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
library(Matrix.utils)
library(edgeR)
library(limma)
# visualzation
library(ComplexHeatmap)
library(ggplot2)
library(pheatmap)
library(scales)
library(UpSetR)
# multiprocess
library(future)

#environment
rm(list = ls())
options(stringsAsFactors = F)

#parameter
##Get the parameters
parser = argparse::ArgumentParser(description="script to find DEG")
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
#path <-'/Users/muxiaohuan/Project/Single_cell/5analysis/06dsda/test'
path <-args$path
path.result <- paste0(path,'/03pb_da')
if (!dir.exists(path.result)) {
  dir.create(path.result)
}
setwd(path.result)

#01prepare_sce.Rdata
rdata <- '../01prepare_sce/01prepare_sce.Rdata'
load(rdata)

###############
#1Visualization
###############
# calculate cluster-sample cell counts
n_cells <- table(sce$cluster_id, sce$sample_id)

# calculate cluster proportions across samples
freqs <- prop.table(n_cells, margin = 2)

# prep. data.frame for plotting
df <- data.frame(
  count = as.numeric(n_cells),
  frequency = as.numeric(freqs), 
  cluster_id = rep(kids, ns),
  sample_id = rep(sids, each = nk))
m <- match(df$sample_id, ei$sample_id)
df$group_id <- ei$group_id[m]

# barplot of relative cluster-abundances
p1 <- ggplot(df, aes(x = sample_id, y = frequency, fill = cluster_id)) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~ group_id, scales = "free_x") +
  theme_classic()  + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# boxplot of relative cluster-abundances
p2 <- ggplot(df, aes(x = group_id, y = frequency, color = group_id)) +
  geom_boxplot(outlier.colour = NA) +  geom_jitter() +
  facet_wrap(~ cluster_id, scales = "free_y", ncol = 4) +
  theme_classic()


###############
#Testing for DA
###############
# construct design & contrast matrix
design <- model.matrix(~ 0 + ei$group_id) %>% 
  set_rownames(ei$sample_id) %>% 
  set_colnames(levels(ei$group_id))

contrast.str <- paste0(levels(ei$group_id)[1],'-',levels(ei$group_id)[2])
contrast <- makeContrasts(contrast.str, levels = design)

# for ea. cluster, run edgeR w/ default parameters
y <- DGEList(counts = n_cells)
y <- estimateDisp(y, design, trend.method = "none")
fit <- glmFit(y, design)
fit <- glmLRT(fit, contrast = contrast)
#topTags(fit, n = Inf)$table %>% round(2)
res_abund <- topTags(fit, n = Inf)$table %>%
  dplyr::mutate(cluster_id = rownames(.))


#######
#output
#######
#abundance
write.table(df,'./03_1_cell_cnt_freq.txt',quote = F,sep = '\t',row.names = F)
#abundance plot
pdf('./03_1_cell_cnt_freq.pdf')
print(p1)
print(p2)
dev.off()
#different abundance
write.table(res_abund,'./03_2_diff_abund.txt',quote = F,sep = '\t',row.names = F)
#rdata
save(n_cells,
     freqs,
     df,
     design,
     contrast,
     file = "./03pb_da.Rdata")

