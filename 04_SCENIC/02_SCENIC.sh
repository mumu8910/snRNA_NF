#Activate conda environment
source activate /dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/yangcan/KY_Project/Bee-brain/conda_environment/new_scenic/

###########
#pyscenic grn
###########
#Candidate regulatory modules are inferred from coexpression patterns between genes 

#Parameter input
tf_list=/dellfsqd3/C_OCEAN/USER/c-muxiaohuan/SingleCell_bee/motif_v2/bee_tf.txt
rna=/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/yangcan/KY_Project/Bee-brain/select_sample/loom/left/NF_filter.loom
out=/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/yangcan/KY_Project/Bee-brain/1_GRNBoost2_loose/NF/NF_adj.tsv

pyscenic grn \
--num_workers 20 \
--output ${out} \
--method grnboost2 \
${rna} \
${tf_list}

###########
#pyscenic ctx
###########
#Coexpression modules are refined by the elimination of indirect targets using TF motif information

#Parameter input
grn_out=/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/yangcan/KY_Project/Bee-brain/1_GRNBoost2_loose/NF/NF_adj.tsv
tf_rank=/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/yangcan/KY_Project/Bee-brain/Creat_database/loose_database/partial/Bee.regions_vs_motifs.rankings.rename.feather
tf_anno=/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/yangcan/KY_Project/Bee-brain/2_cisTarget_loose/motif/motifs-v9-nr.bee-m0.001-o0.1.tbl
rna=/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/yangcan/KY_Project/Bee-brain/select_sample/loom/NF_filter.loom
cis_out=/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/yangcan/KY_Project/Bee-brain/2_cisTarget_loose/NF_reg.csv

pyscenic ctx \
${grn_out} \
${tf_rank} \
--annotations_fname ${tf_anno} \
--expression_mtx_fname ${rna} \
--mode "dask_multiprocessing" \
--output ${cis_out} \
--num_workers 20 \
--mask_dropouts


#############
#pyscenic aucell
#############
#The activity of these discovered regulons is measured in each individual cell and used for clustering

#Parameter input
main=/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/yangcan/KY_Project/Bee-brain/
rna=${main}/select_sample/loom/NF.loom
cis_out=${main}/2_cisTarget_loose/NF_reg.csv
auc_out=${main}/3_AUCell_loose/NF_SCENIC.loom

pyscenic aucell \
${rna} \
${cis_out} \
--output ${auc_out} \
--num_workers 10