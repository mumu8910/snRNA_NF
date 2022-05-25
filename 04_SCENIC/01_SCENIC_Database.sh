#We build SCENIC database referring to the URL below:
#https://github.com/aertslab/create_cisTarget_databases

#Activate conda environment
source activate /dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/yangcan/KY_Project/Bee-brain/conda_environment/create_cistarget_databases

#Parameter input
main=/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/yangcan/KY_Project/Bee-brain/Creat_database/create_cisTarget_databases
fasta_filename=/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/yangcan/KY_Project/Bee-brain/Creat_database/database/fasta/gene_region.fasta
motifs_dir=/dellfsqd3/C_OCEAN/USER/c-muxiaohuan/SingleCell_bee/motif_v2/motif_cb
motifs_list_filename=/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/yangcan/KY_Project/Bee-brain/Creat_database/loose_database/motifs/motifs-v9-nr.bee-m0.001-o0.1.tbl
db_prefix=Bee
nbr_threads=20
nbr_parts=10

#Creat directory
mkdir /dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/yangcan/KY_Project/Bee-brain/Creat_database/database/partial

#create partial databases
for current_part in $( seq 1 10);do
    $main/create_cistarget_motif_databases.py \
    -f "${fasta_filename}" \
    -M "${motifs_dir}" \
    -m "${motifs_list_filename}" \
    -p "${current_part}" "10" \
    -o "partial/${db_prefix}" \
    -t "${nbr_threads}"
done

#Combine
$main/combine_partial_regions_or_genes_vs_motifs_or_tracks_scores_cistarget_dbs.py \
-i ./partial \
-o ./

#Create rankings from a complete cisTarget regions/genes vs motifs scores database
$main/convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py \
-i partial/Bee.motifs_vs_regions.scores.feather 