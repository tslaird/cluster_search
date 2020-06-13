
cluster_df<-read.csv('/home/tslaird/Downloads/FARM_out/iac_positive_all_data.tsv', sep='\t', stringsAsFactors = F)
cluster_df<-cluster_df[cluster_df$domain_gtdb!='',]
metadata_mapping<- read.csv('/home/tslaird/leveau_lab/cluster_search/Figures/metadata_mapping_iac_clusters_6.9.20 .tsv', sep='\t')
cluster_df_meta<-merge(cluster_df,metadata_mapping[,c('assembly_base','Main_category','Subcategory_1','Subcategory_2')], by= 'assembly_base',all.x=TRUE)
cluster_df_meta$synteny_alphabet_nhbr

synteny_metadata<-data.frame(table(cluster_df_meta$synteny_alphabet_nhbr, cluster_df_meta$Subcategory_1))
synteny_metadata=synteny_metadata[synteny_metadata$Freq!=0,]
metadata_mapping

#gene tally total
length(grep('A',cluster_df$synteny_nhbr))/length(cluster_df$synteny_nhbr)
length(grep('B',cluster_df$synteny_nhbr))/length(cluster_df$synteny_nhbr)
length(grep('C',cluster_df$synteny_nhbr))/length(cluster_df$synteny_nhbr)
length(grep('D',cluster_df$synteny_nhbr))/length(cluster_df$synteny_nhbr)
length(grep('E',cluster_df$synteny_nhbr))/length(cluster_df$synteny_nhbr)
length(grep('H',cluster_df$synteny_nhbr))/length(cluster_df$synteny_nhbr)
length(grep('I',cluster_df$synteny_nhbr))/length(cluster_df$synteny_nhbr)



#most represented cluster
data.frame(table(cluster_df$synteny_alphabet_nhbr))

#Acinetobacter species with HABICDE synteny
data.frame(table(cluster_df[cluster_df$genus_gtdb=='g__Acinetobacter',]$synteny_alphabet_nhbr))

#genera with HABICDE synteny
data.frame(table(cluster_df[cluster_df$synteny_alphabet_nhbr=='HABICDE',]$genus_gtdb))

#maximum cluster size
max(cluster_df$feature_count_nhbr)

#number of unique syntenies
length(unique(cluster_df[!grepl('nan',cluster_df$nhbrhood_prot_ids),]$synteny_alphabet_nhbr))
       
       