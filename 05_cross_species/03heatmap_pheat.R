require(dplyr)
require(tidyr)
require(stringr)
require(tibble)
require(pheatmap)
require(RColorBrewer)


#####################
# Receive parameters
#####################
# 0: Accept parameters
parser = argparse::ArgumentParser(description="script to caculate AUROC")

# 1: Retrieve parameters
#'1)--out: output path
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

#workdir
oup_name <- paste0(dataset.A,'-',col.A,'_vs_',dataset.B,'-',col.B)
oup_path = paste0(path,'/metaNeighbor/',oup_name)
if (!dir.exists(oup_path)) {
  dir.create(oup_path)
}
setwd(oup_path)

res.file <- 'auroc.tsv'
res <-  read.delim(res.file)

anno.A.file <- paste0(path,'/data/',dataset.A,'/annotation/',col.A,'.txt')
anno.A <- read.delim(anno.A.file,header=F)
colnames(anno.A) <- c('cluster','type')

anno.B.file <- paste0(path,'/data/',dataset.B,'/annotation/',col.B,'.txt')
anno.B <- read.delim(anno.B.file,header = F)
colnames(anno.B) <- c('cluster','type')

#####
#main
#####
df<-subset(res, grepl(dataset.B,Var1) & grepl(dataset.A,Var2)) %>%
  `row.names<-`(NULL)%>%
  separate(Var1,c('dataset.B','cluster.B'),'[.]',extra='merge') %>%
  separate(Var2,c('dataset.A','cluster.A'),'[.]',extra='merge') %>%
  mutate(type.B = anno.B$type[match(cluster.B,anno.B$cluster)]) %>%
  mutate(type.A = anno.A$type[match(cluster.A,anno.A$cluster)]) %>%
  mutate(Var1 = paste0(type.B,'(',cluster.B,')')) %>%
  mutate(Var2 = paste0(type.A,'(',cluster.A,')'))
write.table(df,'auroc.anno.tsv',sep='\t',quote=F,row.names=F)

#df -> matrix
auroc <- df %>%
  select(Var1,Var2,AUROCs) %>%
  spread(key=Var1,value=AUROCs) %>%
  column_to_rownames(var = 'Var2') %>%
  as.matrix()

#order
hc.row <- hclust(dist(auroc))
auroc <- auroc[hc.row$order,]
hc.col <- hclust(dist(t(auroc)))
auroc <- auroc[,hc.col$order]

#text
text <- matrix(0,nrow(auroc),ncol(auroc))
cutoff <- 0.5
for (i in seq(1,nrow(auroc))) {
  for (j in seq(1,ncol(auroc))) {
    if (auroc[i,j] < cutoff) {
      text[i,j] <- ''
    }else {
      text[i,j] <- as.character(signif(auroc[i,j],2))
    }
  }
}

#color.use
color.up.x <- brewer.pal(6,'Reds')
color.down.x <- brewer.pal(6,'Blues')
color.use <- colorRampPalette(c(rev(color.down.x),'white',color.up.x))(50)
#scales::show_col(color.use)

#plot
oup <- 'auroc.pdf'
main_title <- paste0(dataset.B,' (row) vs ',dataset.A,' (col)')
pheatmap(auroc,
         cluster_cols=T,
         cluster_rows=T,
         #scale = 'row',
         color = color.use, 
         na_col = 'white',
         cellwidth = 17,
         cellheight = 17,
         fontsize=8,
         fontsize_row=6,
         fontsize_col=6,
         display_numbers = text,
         number_color = 'black',
         show_rownames = T,
         show_colnames = T,
         angle_col = '45',
         main = main_title,
         filename = oup)

