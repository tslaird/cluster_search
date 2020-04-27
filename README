The cluster search pipeline is used for finding iac gene clusters within RefSeq genomes

The pipeline is implemented in python using several python scripts ".py"

To download/unzip bacterial genomes as ".gbff" files from the RefSeq database and check their integrity based on md5 sums run:
```
python download_genomes.py
```

To parse the ".gbff" files and create ".fasta" protein files with informative headers run:
```
python parse_gbff_files.py
```
Combining the ".fasta" files and subsequently creating blast databases from those files:
```
python make_blastdb.py
```
Running the blastp algorithm using iac genes as a query against the previously made databases:
```
python run_blastp.py
```
Parsing the blastp output, identify iac gene clusters, and fetch appropriate metadata run:
```
python parse_blastp.py
```
The output of the results file is given in several file formats including a tsv file, xlsx file, and ".pickle" file (for use in python)

The table headers include:
accession:			the accession number of the contig containing the cluster
assembly:			the name of the RefSeq genome assembly the cluster is found in (e.g. "GCF_XXXXXXX.X)
title:				the name of the contig containing the cluster
feature_count_nhbr:		a count of genomic features in the cluster neighborhood
cluster_len_nhbr:		the length in nucleotides of the cluster neighborhood
synteny_nhbr:			
synteny_dir_dist_nhbr
synteny_dir_nhbr
cluster_number
adj_coord_list
tared_adj_coord_list
itol_cluster_string
nhbrhood_hit_list:		a list of 
nhbrhood_locus_tags:		a list of locus tags" from proteins in the cluster neighborhood
nhbrhood_old_locus_tags:	a list of "old locus tags" from proteins in the cluster neighborhood
nhbrhood_prot_ids		a list of protein ids from proteins in the cluster neighborhood
nhbrhood_prot_name		a list of protein name from proteins in the cluster neighborhood
nhbrhood_prot_seq:		a list of protein sequences from proteins in the cluster neighborhood
clusterGC:			the GC content of the cluster
genomeGC:			the GC content of the entire genome
diffGC:				the difference in GC content between the cluster and entire genome
cluster_seq:			the DNA sequence of the cluster
genome_acc
filename:			the local filename of the genome ".gbff" file
biosample:			the biosample id of the genome
hits
cluster_length
synteny	synteny_dir_dist
synteny_dir_pident
synteny_dir
name
hit_list			a list of the "locus tags" for a cluster parsed from the gbff file
old_locus_hit_list:		a list of the "old locus tags" for a cluster parsed from the gbff file
protein_name_list
protein_id_list	
#pseudogene_list
query_list:			the list of the 
coord_list:			the coordinate list of the cluster genes
contig:			
complete_genome:		whetehr or not the cluster comes from a complete genome based on pattern matching with the name
plasmid:			whether or not the cluster is on a plasmid based on pattern matching with the name
has_overlap:			whether or not genes in the cluster are annotated as overlapping
duplicated:			whether or not the cluster appears multiple times in a genome
isolation_src:			the isolation source fetched from the biosample report
env_biome:			the environmental biome fetched from the biosample report
env_feature:			the environmental feature fetched from the biosample report
host:				the host name fetched from the biosample report
host_sciname:			the host scientific name fetched from the biosample report
strain:				strain information fetched from the biosample report
all_biosample_metadata:		all the metadata fetched from the biosampple report
gtdb_tax:			all gtdb taxonomy category classifications according to the gtdb database file
ncbi_tax:			all ncbi taxonomy category classifications according to the gtdb database file
same_taxonomy:			whether or not the gtdb and ncbi taxonomies are the same
domain_gtdb:			the gtdb assigned domain
phylum_gtdb:			the gtdb assigned phylum
class_gtdb:			the gtdb assigned class
order_gtdb:			the gtdb assigned order
family_gtdb:			the gtdb assigned family
genus_gtdb:			the gtdb assigned genus
species_gtdb:			the gtdb assigned species

