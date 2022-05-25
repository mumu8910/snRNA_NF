#===================
#  01_PrepareSce.R
#===================
This script will convert the Seurat object to SingCellExperiment object

#===================
#  02_PseudoBulkDEG.R
#===================
The number of UMI for each gene from all cells in the cluster from each library were added together and normalized by Trimmed Mean of M-values (TMM) method. 
Differential gene expression analysis based on cluster-level pseudo bulk mainly using glmQLFit and glmQLFTest function in R package EdgeR.

