
cluster_df<-read.csv('/home/tslaird/Downloads/FARM_out/iac_positive_all_data.tsv', sep='\t', stringsAsFactors = F)
cluster_df<-cluster_df[cluster_df$domain_gtdb!='',]
metadata_mapping<- read.csv('/home/tslaird/leveau_lab/cluster_search/Figures/metadata_mapping_iac_clusters_6.9.20 .tsv', sep='\t')
cluster_df_meta<-merge(cluster_df,metadata_mapping[,c('assembly_base','Main_category','Subcategory_1','Subcategory_2')], by= 'assembly_base',all.x=TRUE)
cluster_df_meta$synteny_alphabet_nhbr

synteny_metadata<-data.frame(table(cluster_df_meta$synteny_alphabet_nhbr, cluster_df_meta$Subcategory_1))
synteny_metadata=synteny_metadata[synteny_metadata$Freq!=0,]
metadata_mapping


