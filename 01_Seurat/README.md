#================
#  01_Integration.R
#================
Using FindIntegrationAnchors and IntegrateData functions, all libraries were integrated into a single matrix with the first 30 principle components with the following parameter values: reduction = "rpca", normalization.method = "LogNormalize"

#================
#  02_FindCluster.R
#================
This script find clusters using FindNeighbors and FindClusters in Seurat with different resolution parameter (0.2, 0.5, 0.8, 1.2, 1.6), followed by clustree analysis.

#================
#  03_FindMarkers.R
#================
Markers was defined utilizing FindAllMarkers function in this script. Significantly overexpressed genes were defined based on the Wilcoxon rank-sum test with a P value < 0.01 and log2 FC > 0.1

#================
#  04_GOenrich.R
#================
Gene Ontology enrichment analyses were performed using the R package clusterProfiler in three levels (BP, Biological Process; CC, Cell Components; MF, Molecular Function). This script can describe the enrichment function of marker genes.
