#==================
#  01_PrepareSce.R
#==================
This script will convert the Seurat object to SingCellExperiment object

#==================
#  02_CellRatioShifts.R
#==================
In this script, we do the following things:
(1) Calculate frequency of clusters across samples; 
(2) Cell ratio shifts in different groups were conducted by edgeR, mainly using function glmFit and glmLRT.
