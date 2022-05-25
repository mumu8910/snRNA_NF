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
library(stringr)
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
path <-args$path
path.result <- paste0(path,'/02pb_ds')
if (!dir.exists(path.result)) {
  dir.create(path.result)
}
setwd(path.result)
if (!dir.exists('./02_2_diff_gene')) {
  dir.create('./02_2_diff_gene')
}
if (!dir.exists('./02_3_diff_gene_cnt')) {
  dir.create('./02_3_diff_gene_cnt')
}
#01prepare_sce.Rdata input
rdata <- '../01prepare_sce/01prepare_sce.Rdata'
load(rdata)

###########
#pseudobulk
###########
# aggregate by cluster-sample
groups <- colData(sce)[, c("cluster_id", "sample_id")]
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

# split by cluster, transform & rename columns
pb.split.cluster <- str_split_fixed(rownames(pb),'_',2)[,1]
pb.split.sample <- str_split_fixed(rownames(pb),'_',2)[,2]
pb <- split.data.frame(pb,pb.split.cluster) %>% 
  lapply(function(u){
    u.split = str_split_fixed(rownames(u),'_',2)[,2]
    set_colnames(t(u), u.split)
  })

##########################
#pseudobulk-level MDS plot
##########################
# compute MDS coordinates
mds <- pb %>% 
  lapply(as.data.frame.matrix) %>% 
  bind_cols %>% 
  DGEList(remove.zeros = TRUE) %>% 
  calcNormFactors %>% 
  plotMDS.DGEList(plot = FALSE)

# prep. data.frame for plotting
gg_df <- data.frame(mds[c("x", "y")],
                    cluster_id = pb.split.cluster,
                    sample_id = pb.split.sample,
                    group_id = ei$group_id[match(pb.split.sample, ei$sample_id)])

p <- ggplot(gg_df, aes(x, y, col = cluster_id, shape = group_id)) + 
  geom_point(size = 3, alpha = 0.8) +
  labs(x = "MDS dim. 1", y = "MDS dim. 2") + 
  theme(panel.grid.minor = element_blank()) +
  coord_fixed() + theme_classic()

#output
write.table(gg_df,'./02_1_pb_mds.txt',quote = F,sep = '\t',row.names = F)
pdf('./02_1_pb_mds.pdf',width = 5,height = 5)
print(p)
dev.off()

#####################################
#pseudobulk-level different expr gene
#####################################
# construct design & contrast matrix
design <- model.matrix(~ 0 + ei$group_id) %>% 
  set_rownames(ei$sample_id) %>% 
  set_colnames(levels(ei$group_id))

contrast.str <- paste0(levels(ei$group_id)[1],'-',levels(ei$group_id)[2])
contrast <- makeContrasts(contrast.str, levels = design)

# for ea. cluster, run edgeR w/ default parameters
res <- lapply(kids, function(k) {
  y <- pb[[k]]
  y <- DGEList(y, remove.zeros = TRUE)
  y <- calcNormFactors(y)
  y.design <- design[colnames(y),]
  y <- estimateDisp(y, y.design)
  fit <- glmQLFit(y, y.design)
  fit <- glmQLFTest(fit, contrast = contrast)
  topTags(fit, n = Inf, sort.by = "none")$table %>% 
    dplyr::mutate(gene = rownames(.), cluster_id = k) %>% 
    dplyr::rename(p_val = PValue, p_adj = FDR)
})

# filter FDR < 0.05, |logFC| > 1 & sort by FDR
res_fil <- lapply(res, 
                  function(u)  u %>% 
                    dplyr::filter(p_adj < 0.05, abs(logFC) > 1) %>% 
                    dplyr::arrange(p_adj))

###################################
#different expr gene stat & overlap
###################################
# nb. & % of DE genes per cluster
n_de <- vapply(res_fil, nrow, numeric(1))
n_de_stat.df <- data.frame(cbind(n_de, p_gs = n_de / nrow(sce) * 100))
n_de_stat.df$cluster_id <- rownames(n_de_stat.df)

# diff gene overlap
n_de.df <- fromList(map(res_fil, "gene"))

# Jaccard coefficient
n_de_jac <- mapply(
  FUN = function(ident1,ident2) {
    v1 <- n_de.df[ident1]=='1';v1.sum <- sum(v1)
    v2 <- n_de.df[ident2]=='1';v2.sum <- sum(v2)
    v <- v1&v2;v.sum <- sum(v)
    jac <- v.sum/(v1.sum+v2.sum-v.sum)
    return(list(cluster_id1 = ident1,cluster_id2 = ident2,jaccard = jac))
  },
  rep(kids,nk),
  rep(kids,times=rep(nk,nk)),
  SIMPLIFY = F,
  USE.NAMES = F
)
n_de_jac <- lapply(
  X = 1:length(n_de_jac),
  FUN = function(x) {
    data.use <- as.data.frame(x = n_de_jac[[x]])
    return(data.use)
  }
)
n_de_jac <- do.call(what = 'rbind', args = n_de_jac)

#######
#output
#######
#DEG list
lapply(res_fil, 
       function(u){
         write.table(u,paste0('./02_2_diff_gene/',unique(u$cluster_id),'.txt'),quote = F,sep = '\t',row.names = F)
       } )
#count for DEG
write.table(n_de_stat.df,'./02_3_diff_gene_cnt/diff_gene_cnt_stat.txt',quote = F,sep = '\t',row.names = F)
#DEG overlap:Jaccard coefficient
write.table(n_de_jac,'./02_3_diff_gene_cnt/diff_gene_jaccard.txt',quote = F,sep = '\t',row.names = F)
#rdata
save(pb,
     mds,
     groups,
     design,
     contrast,
     res,res_fil,
     file = "./02pb_ds.Rdata")
#DEG overlap
write.table(n_de.df,'./02_3_diff_gene_cnt/diff_gene_cnt.txt',quote = F,sep = '\t',row.names = F)
pdf('./02_3_diff_gene_cnt/diff_gene_cnt.pdf')
upset(n_de.df)
dev.off()

