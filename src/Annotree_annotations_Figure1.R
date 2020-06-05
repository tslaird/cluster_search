library("treeio")
library("ggtree")
library("phytools")
#must install magick
#sudo apt-get install -y libmagick++-dev
library(ggplot2)
library(phangorn)
library(phylobase)

#functions
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
get_cladelabel_position_3 <- function(data, node, angle = "auto", extend = 0) {
  if (length(extend) == 1) {
    extend = rep(extend, 2)
  }
  ## sp <- get.offspring.df(data, node)
  ## sp2 <- c(sp, node)
  ## sp.df <- data[match(sp2, data$node),]
  sp.df <- offspring.tbl_tree(data, node, self_include = TRUE)
  y <- sp.df$y
  y <- y[!is.na(y)]
  mx <- max(sp.df$x, na.rm=TRUE)
  d <- data.frame(x=mx, y=min(y) - extend[2], yend=max(y) + extend[1])
  if (missing(angle))
    return(d)
  if (angle == "auto") {
    d$angle <- mean(range(sp.df$angle))
  } else {
    d$angle <- angle
  }
  return(d)
}

# tree must be downloaded from AnnoTree website
tree <- read.tree("/home/tslaird/leveau_lab/cluster_search/gtdb/tree_of_life.newick")
#cluster_search output
cluster_df<-read.csv('/home/tslaird/leveau_lab/cluster_search/FARM_output/iac_positive_all_data.tsv', sep='\t', stringsAsFactors = F)
all_nodes<-c(unique(cluster_df$phylum_gtdb), unique(cluster_df$class_gtdb),unique(cluster_df$order_gtdb),unique(cluster_df$family_gtdb),unique(cluster_df$genus_gtdb))
all_nodes<-gsub("^\\w{3}",'',all_nodes)

#all_nodes_GTDB<-read.csv('/Users/tslaird/Box/leveau_lab/AnnoTree_trees/all_nodes_GTDB.txt',header=F)
all_indices<-unlist(sapply(all_nodes, function(x)  which(tree$node.label== x)))
all_tips<-unlist(sapply(all_nodes, function(x) which(tree$tip.label== x)))
all_ancestors=c(unlist(Ancestors(tree,all_indices,'all' )) ,all_indices )

tree_phylo4 = as(tree, 'phylo4')
#all_indices<-unlist(sapply(all_nodes, function(x) which(tree_phylo4@label== x)))
#all_ancestors=as.numeric(c(unlist(ancestors(tree_phylo4,all_indices,'all' )) ,all_indices ))

nodes=1:(nrow(tree$edge)+1)

#download from GTDB database
GTDB_meta <- read.csv('/home/tslaird/leveau_lab/cluster_search/gtdb/bac120_metadata_r89.tsv', sep='\t')

p<-ggtree(tree_phylo4,layout = 'circular') %<+% d + aes(color=I(color), size=I(size))
all_gtdb_classes=data.frame( 'class' = unlist(unique(stringr::str_extract_all(as.character(GTDB_meta$gtdb_taxonomy),'(?<=p__).+?(?=;)'))))
all_gtdb_classes$node = sapply(all_gtdb_classes$class, function(x) which(tree_phylo4@label==x)[1])
all_gtdb_classes$angle<- sapply(all_gtdb_classes$class, function(x) -p$data[p$data$label==x,]$angle[1] )
all_gtdb_classes_filtered <- all_gtdb_classes[(!is.na(all_gtdb_classes$node)) & (!grepl('^UBA|^UBP|-|\\d|^GCA|^GCF',all_gtdb_classes$class)),]

d = data.frame(node=1:(length(tree_phylo4@label)), color=c(ifelse((1:length(tree_phylo4@label)) %in% all_ancestors,'blue','darkgrey')), size=c(ifelse((1:length(tree_phylo4@label)) %in% all_ancestors,0.3,0.05)))

for(n in c(1:(nrow(all_gtdb_classes_filtered)))){
  clade_lab<-geom_cladelabel(node=all_gtdb_classes_filtered[n,2], label=all_gtdb_classes_filtered[n,1],angle='auto', hjust = 'right',geom='text',fontsize=3.2)
  calc<-get_cladelabel_position_3(data=p$data,node=all_gtdb_classes_filtered[n,2],angle = 'auto')
  calc_angle<-calc$angle
  calc_size<- calc$yend-calc$y
  #filter clades by size
  if(calc_size>30){
  clade_lab$angle<- ifelse((calc_angle <=90 | calc_angle>=270), calc_angle ,calc_angle +180)
  clade_lab$hjust<- ifelse((calc_angle <=90 | calc_angle>=270), 'left','right')
  p<-p+clade_lab
  }
}

p<-p+geom_hilight(node=9146, fill="#BB0044", alpha=.4)+
  geom_hilight(node=10056, fill="#FF8F3F", alpha=.4)+
  geom_hilight(node=6994, fill="#4BD0E4", alpha=.4)
p


ggsave(plot=p,'test.svg', device = 'svg', height = 10,width=10 )
