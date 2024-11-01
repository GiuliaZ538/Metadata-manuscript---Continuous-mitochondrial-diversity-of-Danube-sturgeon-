#MetaDMG analysis
##Creating config file
metaDMG config /path/to/files/*.bam --names path/to/ncbi_taxonomy/names.dmp --nodes path/to/ncbi_taxonomy/nodes.dmp --acc2tax path/to/ncbi_taxonomy/accession2taxid.gz

##Config file settings
samples: file_name_no_duplicates: /path/to/file_no_duplicates.bam
metaDMG_cpp: /path/to/metaDMG-cpp/metaDMG-cpp
names: path/to/ncbi_taxonomy/names.dmp
nodes: path/to/ncbi_taxonomy/nodes.dmp
acc2tax: path/to/ncbi_taxonomy/accession2taxid.gz #contains all accessions taxids used
min_similarity_score: 0.95
max_similarity_score: 1.0
min_edit_dist: 0
max_edit_dist: 10
min_mapping_quality: 0
lca_rank: ''
max_position: 15
weight_type: 1
custom_database: false
forward_only: true
output_dir: /path/to/results_folder
parallel_samples: 1
cores_per_sample: 1
bayesian: false
config_file: config.yaml
damage_mode: global
version: 0.24.7

##Generate final output
metaDMG convert --output /path/to/final_output.csv --results /path/to/results_folder/results
