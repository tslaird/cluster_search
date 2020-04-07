The cluster search pipeline is used for finding iac gene clusters within RefSeq genomes

The pipeline is implemented in Snakemake and involves the following steps:

Downloading bacterial genomes as ".gbff" files from the RefSeq database and checking their integrity based on md5 sums
Parsing the ".gbff" files to create ".fasta" protein files with informative headers
Combining the ".fasta" files and subsequently creating blast databases from those files
Running the blastp algorithm using iac genes as a query against the previously made databases
Parsing the resulting blastp output to identify iac gene clusters within genomes
Fetching metadata for genomes containing iac clusters
Downloading and cross-checking identified iac carrying genomes against the genome taxonomy database (GTDB)


The output of the results file includes:
accession
assembly
title
feature_count_nhbr
cluster_len_nhbr
synteny_nhbr
synteny_dir_dist_nhbr
synteny_dir_nhbr
cluster_number
adj_coord_list
tared_adj_coord_list
itol_cluster_string
nhbrhood_hit_list
nhbrhood_locus_tags
nhbrhood_old_locus_tags
nhbrhood_prot_ids
nhbrhood_prot_name
nhbrhood_prot_seq
clusterGC
genomeGC
diffGC
cluster_seq
genome_acc
filename
biosample
hits
cluster_length
synteny	synteny_dir_dist
synteny_dir_pident
synteny_dir
name
hit_list
old_locus_hit_list
protein_name_list
protein_id_list	pseudogene_list
query_list
coord_list
contig
complete_genome
plasmid
has_overlap
duplicated
isolation_src
env_biome
env_feature
host
host_sciname
strain
all_biosample_metadata
gtdb_tax
ncbi_tax
same_taxonomy
domain_gtdb
phylum_gtdb
class_gtdb
order_gtdb
family_gtdb
genus_gtdb
species_gtdb

