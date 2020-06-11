

cluster_df<-read.csv('/home/tslaird/Downloads/FARM_out/iac_positive_all_data.tsv', sep='\t', stringsAsFactors = F)
cluster_df<-cluster_df[cluster_df$domain_gtdb!='',]

mean(cluster_df[cluster_df$order_gtdb=='o__Enterobacterales',]$four_mer_distance)

mean(cluster_df[cluster_df$order_gtdb='o__Enterobacterales',]$four_mer_distance)

t.test(cluster_df[cluster_df$genus_gtdb=='g__Escherichia',]$four_mer_distance, 
       cluster_df[cluster_df$genus_gtdb!='g__Escherichia',]$four_mer_distance)


mean(cluster_df$diffGC)

library(plyr)
GC_tab<-ddply(cluster_df, .(genus_gtdb), summarise, mean(diffGC))

x<-cluster_df[cluster_df$genus_gtdb=='g__Escherichia',]
x<-cluster_df[cluster_df$clusterGC>0.5,]

#mean GC content of iac in genus Escherichia
mean(cluster_df[cluster_df$genus_gtdb=='g__Escherichia',]$diffGC)
mean(cluster_df[cluster_df$genus_gtdb=='g__Klebsiella',]$diffGC)

summary(lm(cluster_df$diffGC~cluster_df$four_mer_distance))


cluster_df[cluster_df$diffGC>0.07,]$order_gtdb
an<-aov(cluster_df$diffGC ~ cluster_df$order_gtdb)
TukeyHSD(an)

# identify iac clusters with mobilization genes nearby
cluster_df_10<-read.csv('/home/tslaird/leveau_lab/cluster_search/FARM_output/iac_positive_all_data_10_10.tsv',sep='\t')
mobilization<-cluster_df_10[grep('integrase|transposase|phage',cluster_df_10$nhbrhood_prot_name),]$accession
mobile_df<-cluster_df[cluster_df$accession %in% mobilization,]
length(unique(mobile_df$assembly))
mobile_df[mobile_df$diffGC>0.06,]$species_gtdb
table(mobile_df$family_gtdb)

#tabulate percent per taxonomy
GTDB_meta <- read.csv('/home/tslaird/leveau_lab/cluster_search/gtdb/bac120_metadata_r89.tsv', sep='\t')
GTDB_meta_RS<-GTDB_meta[grep('RS_',GTDB_meta$accession),]
library(dplyr)
library(tidyr)
GTDB_meta_RS<- GTDB_meta_RS %>% separate(gtdb_taxonomy,  c('domain_gtdb','phylum_gtdb','class_gtdb','order_gtdb','family_gtdb','genus_gtdb','species_gtdb'), ";")
GTDB_meta_RS$accession_base<-sapply(strsplit(as.character(GTDB_meta_RS$accession),'\\.|RS_'), function(x) x[2])

genus_tally<-rbind.data.frame(t(sapply(unique(cluster_df$genus_gtdb), function(x){
  check=x
  pos<-nrow(GTDB_meta_RS[((GTDB_meta_RS$accession_base %in% cluster_df$assembly_base) & (GTDB_meta_RS$genus_gtdb==check)),])
  neg<-nrow(GTDB_meta_RS[((!GTDB_meta_RS$accession_base %in% cluster_df$assembly_base) & (GTDB_meta_RS$genus_gtdb==check)),])
  c(check,pos,neg,(pos/(neg+pos))) })))
