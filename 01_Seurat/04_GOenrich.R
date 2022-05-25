########
#prepare
########
#environment
rm(list = ls())
options(stringsAsFactors = F)

#package
library(dplyr);cat("dplyr:",as.character(packageVersion("dplyr")),"\n")
library(stringr);cat("stringr:",as.character(packageVersion("stringr")),"\n")
library(clusterProfiler);cat("clusterProfiler:",as.character(packageVersion("clusterProfiler")),"\n")
library(ggplot2);cat("ggplot2:",as.character(packageVersion("ggplot2")),"\n")
library(ggpubr);cat("ggpubr:",as.character(packageVersion("ggpubr")),"\n")
library(patchwork);cat("patchwork:",as.character(packageVersion("patchwork")),"\n")
library(grid);cat("grid:",as.character(packageVersion("grid")),"\n")
library(future);cat("future:",as.character(packageVersion("future")),"\n")

##########
#parameter
##########
##Get the parameters
parser = argparse::ArgumentParser(description="script to analysis DEG go function")
parser$add_argument('-r','--res_usage_i',help='resolution usage')
parser$add_argument('-o','--out',help='out directory')
parser$add_argument('-a','--anno_config', help='annotation config directory')
parser$add_argument('-t','--thread', help='thread used')
parser$add_argument('-g','--ram', help='ram used(Gb)')
args = parser$parse_args()

##resolution
column <- args$res_usage_i
check <- is.na(as.numeric(column))
##workdir
path <- args$out
path.col <- (if(check) paste0(path,'/',column) else paste0(path,'/','res_',column))
if (!dir.exists(path.col)) {
  dir.create(path.col)
}
setwd(path.col)

##file
anno_file <- paste0(args$anno_config,'/apis_go.tsv')
gene_id_file <- paste0(args$anno_config,'/gene_gtf.tsv')
gene_file <- './02_diff_marker.csv'
##multiprocess
thread <- as.numeric(if(!is.null(args$thread)) args$thread else 5)
ram <- as.numeric(if(!is.null(args$ram)) args$ram else 2)
plan("multiprocess",workers = thread)
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = ram*1024^3)

###########
#annotation
###########
#annotation file
anno_file <- anno_file
anno_df <- read.delim(anno_file,header = T,sep = '\t')
#gene_id convert file
gene_id_file <- gene_id_file
gene_id_df <- read.delim(gene_id_file,header = T,sep = '\t')
#go_id|gene_name
term2gene.BP <- dplyr::filter(anno_df,go_ontology == 'BP') %>% select(c('go_id','gene_id')) 
term2gene.CC <- dplyr::filter(anno_df,go_ontology == 'CC') %>% select(c('go_id','gene_id')) 
term2gene.MF <- dplyr::filter(anno_df,go_ontology == 'MF') %>% select(c('go_id','gene_id')) 
#gene_id convert
term2gene.BP$gene_id <- gene_id_df$gene_name[match(term2gene.BP$gene_id,gene_id_df$gene_id)]
term2gene.CC$gene_id <- gene_id_df$gene_name[match(term2gene.CC$gene_id,gene_id_df$gene_id)]
term2gene.MF$gene_id <- gene_id_df$gene_name[match(term2gene.MF$gene_id,gene_id_df$gene_id)]
#go_id|go_term
term2name.BP <- dplyr::filter(anno_df,go_ontology == 'BP') %>% select(c('go_id','go_term'))
term2name.CC <- dplyr::filter(anno_df,go_ontology == 'CC') %>% select(c('go_id','go_term'))
term2name.MF <- dplyr::filter(anno_df,go_ontology == 'MF') %>% select(c('go_id','go_term'))

#####
#gene
#####
#gene file
gene_file <- gene_file
gene_df <- read.csv(gene_file,header = T,sep = ',')

###############################################
#select intrested genes and enrichment analysis
###############################################
#gene extract
oup_df.BP <- data.frame()
oup_df.CC <- data.frame()
oup_df.MF <- data.frame()
cluster_i <- unique(gene_df$cluster)
#i <- 6
for (i in seq(1,length(cluster_i))) {
  #sig:p_val<0.05
  gene_i <- dplyr::filter(gene_df,cluster==cluster_i[i],p_val<0.05)$gene
  #enrich
  #BP
  enrich_df.BP <- enricher(gene_i,TERM2GENE = term2gene.BP,TERM2NAME = term2name.BP,
                           pvalueCutoff = 1,pAdjustMethod = 'BH',qvalueCutoff = 0.5)
  enrich_df.result.BP <- enrich_df.BP@result
  enrich_df.result.BP.oup <- enrich_df.result.BP
  enrich_df.result.BP$cluster_id <- cluster_i[i]
  oup_df.BP <- rbind(oup_df.BP,enrich_df.result.BP)
  #CC
  enrich_df.CC <- enricher(gene_i,TERM2GENE = term2gene.CC,TERM2NAME = term2name.CC,
                           pvalueCutoff = 1,pAdjustMethod = 'BH',qvalueCutoff = 0.5)
  enrich_df.result.CC <- enrich_df.CC@result
  enrich_df.result.CC <- enrich_df.CC@result
  enrich_df.result.CC.oup <- enrich_df.result.CC
  enrich_df.result.CC$cluster_id <- cluster_i[i]
  oup_df.CC <- rbind(oup_df.CC,enrich_df.result.CC)
  #MF
  enrich_df.MF <- enricher(gene_i,TERM2GENE = term2gene.MF,TERM2NAME = term2name.MF,
                           pvalueCutoff = 1,pAdjustMethod = 'BH',qvalueCutoff = 0.5)
  enrich_df.result.MF <- enrich_df.MF@result
  enrich_df.result.MF <- enrich_df.MF@result
  enrich_df.result.MF.oup <- enrich_df.result.MF
  enrich_df.result.MF$cluster_id <- cluster_i[i]
  oup_df.MF <- rbind(oup_df.MF,enrich_df.result.MF)
}
write.table(oup_df.BP,'./cluster_go_bp.tsv',row.names = F,quote = F,sep = '\t')
write.table(oup_df.CC,'./cluster_go_cc.tsv',row.names = F,quote = F,sep = '\t')
write.table(oup_df.MF,'./cluster_go_mf.tsv',row.names = F,quote = F,sep = '\t')


