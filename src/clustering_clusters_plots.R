library(dendextend)

distmat_cluster<-read.csv('/Users/tslaird/Box/leveau_lab/iac_database/distmat_cluster.csv')
distmat_cluster2<-read.csv('/Users/tslaird/Box/leveau_lab/iac_database/distmat_cluster_twos.csv')
distmat_cluster3<-read.csv('/Users/tslaird/Box/leveau_lab/iac_database/distmat_cluster_threes.csv')
distmat_cluster4<-read.csv('/Users/tslaird/Box/leveau_lab/iac_database/distmat_cluster_fours.csv')
mat2<-distmat_cluster2[,2:ncol(distmat_cluster2)]
mat3<-distmat_cluster3[,2:ncol(distmat_cluster3)]
mat4<-distmat_cluster4[,2:ncol(distmat_cluster4)]

consensus_mat<-(mat2+mat3+mat4)/3
mat<-mat2


mat<-distmat_cluster2[,2:ncol(distmat_cluster2)]
mat_names<- gsub('〉','r',  distmat_cluster2[,1]) 
mat_names<- gsub('〈','l',  mat_names)
mat_names<- gsub(' ','',  mat_names)
mat_names<- gsub('|','I',  mat_names)
rownames(mat)<- mat_names
colnames(mat)<-mat_names

methods<-c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")
correlation<-c()
for(i in seq(1:length(methods))){
  hc<- hclust(as.dist(mat),method=methods[i])
  hc_coph<-cophenetic(hc)
  correlation[i]<- cor(as.dist(mat),hc_coph)
}
results<-cbind.data.frame(methods,correlation)
print(results[order(-results$correlation),])


hclust<-hclust(as.dist(mat), method='average')
groups<-cutree(hclust, k = 8, h = NULL)
names(groups)<-distmat_cluster[,1]
groups_df<-as.data.frame(groups)
par(mar = c(2, 2, 2, 25))
dend<-as.dendrogram(hclust)
hclust$labels <-rep(expression(symbol("\341")),118)
nodePar <- list(lab.cex = 0.5, pch = c(NA, 19), 
                cex = 0.7)
plot(dend,horiz=TRUE,nodePar=nodePar )

clus = cutree(hclust, h=0.3)
colors<-randomColor(max(clus))
plot(as.phylo(hclust), tip.color = colors[clus],
     label.offset = 0.01, cex = 0.5)

plot_dendro(dend,height = 600)
ggdendrogram(hclust)
plot(as.phylo(hclust), type = "unrooted", cex = 0.6,
     no.margin = TRUE)

hc <- hclust(dist(USArrests))
dend1 <- as.dendrogram(hc)
plot_dendro(dend ) %>% 
  hide_legend() %>% 
  highlight(persistent = TRUE, dynamic = TRUE)

hc <- hclust(dist(mtcars))
den<- as.phylo(hclust)
clus <- cutree(hc, 4)
g <- split(names(clus), clus)

p <- ggtree(den, linetype='dashed') + geom_tiplab()

geom_tiplab(aes(label=den$tip.label), size=3, hjust=.5, color='black')
p
clades <- sapply(g, function(n) MRCA(p, n))

p <- groupClade(p, clades, group_name='subtree') + aes(color=subtree)

d <- data.frame(label = names(clus), 
                cyl = mtcars[names(clus), "cyl"])

library(heatmaply)
x<-1-distmat_cluster2[,2:ncol(distmat_cluster2)]
rownames(x)<-NULL
colnames(x)<-NULL
label_names = gsub('〉|\\||〈| |〈','' , rownames(NMDS2D_co))
heatmaply(x,labRow =gsub('〉|\\||〈| |〈','' , rownames(NMDS2D_co)),
          labCol =gsub('〉|\\||〈| |〈','' , rownames(NMDS2D_co)),
          margins = c(60,100,40,20))

length(label_names)
NMDS_2D_co<-metaMDS(distmat_cluster[,2:ncol(distmat_cluster)],k =2,try=20,trymax=20)
NMDS2D_co <- data.frame(NMDS_1 = NMDS_2D_co$points[,1], NMDS_2 = NMDS_2D_co$points[,2]  )
rownames(NMDS2D_co)<-distmat_cluster[,1]
NMDS2D_co_plot <- plotly::plot_ly()
NMDS2D_co_plot<-add_markers(NMDS2D_co_plot,data=NMDS2D_co, x = NMDS2D_co$NMDS_1, y = NMDS2D_co$NMDS_2,color = NMDS2D_co$class_gtdb,
                              colors=pal,
                              marker=list(line=list(color='black', width=0.5),size=8,opacity = 1 ),
                              hoverinfo = 'text',  text = rownames(NMDS2D_co))
NMDS2D_co_plot


pcoa3D<- cmdscale(as.dist(mat),k=3)
pcoa3D<-as.data.frame(pcoa3D)
pcoa3D<- cbind.data.frame(pcoa3D)
pcoa3D_plot <- plotly::plot_ly(pcoa3D, x = ~V1, y = ~V2,z= ~V3,
                               colors=pal,
                               marker=list(line=list(color="black", width=1.0),size=5,opacity = 1 ),
                               hoverinfo = 'text',  text = rownames(NMDS2D_co))
pcoa3D_plot

plot_pcoa3d(distmat_cluster)

#----
library(dendextend)
distmat_cluster2<-read.csv('/Users/tslaird/Box/leveau_lab/iac_database/distmat_cluster_twos.csv')
mat<-distmat_cluster2[,2:ncol(distmat_cluster2)]
mat_names<- gsub('〉','r',  distmat_cluster2[,1]) 
mat_names<- gsub('〈','l',  mat_names)
mat_names<- gsub(' ','',  mat_names)
mat_names<- gsub('\\|','J',  mat_names)
rownames(mat)<- mat_names
colnames(mat)<-mat_names


df<- read.csv('/Users/tslaird/Box/leveau_lab/FARM_output/iac_search_09122019_no_pfam_filter/results_09122019_nopfam_filter/iac_positive_all_data.tsv', sep = '\t')

gtdb<- cluster_df[cluster_df$gtdb_tax !='' & !sapply(cluster_df$nhbrhood_prot_ids,function(x) grepl(x,'nan')),]
gtdb$synteny_dir_nhbr<- sapply(gtdb$synteny_dir_nhbr,function(x) gsub('〉','r',x))
gtdb$synteny_dir_nhbr<- sapply(gtdb$synteny_dir_nhbr,function(x) gsub('〈','l',x))
gtdb$synteny_dir_nhbr<- sapply(gtdb$synteny_dir_nhbr,function(x) gsub(' ','',x))
gtdb$synteny_dir_nhbr<- sapply(gtdb$synteny_dir_nhbr,function(x) gsub('\\|','J',x))


synteny_dist<-adist(unique(gtdb$synteny_alphabet_nhb))
colnames(synteny_dist)<-unique(gtdb$synteny_alphabet_nhb)
rownames(synteny_dist)<-unique(gtdb$synteny_alphabet_nhb)
plot(hclust(as.dist(synteny_dist)))


gtdb<-gtdb[gtdb$synteny_dir_nhbr %in% colnames(mat),]
crosstab<-as.data.frame.matrix(table(gtdb$synteny_alphabet_nhbr,gtdb$family_gtdb))[-1]
colnames(crosstab)<-gsub('f__','',colnames(crosstab))
#crosstab_rel<-t(apply(crosstab,1,function(x) x/sum(x)))
#hclust<-hclust(as.dist(mat), method='average')
hclust<-hclust(as.dist(synteny_dist))
#crosstab_rel<-crosstab_rel[hclust$labels[hclust$order],]
crosstab_rel_binary<-ifelse(crosstab>0,1,0)
x  <- as.matrix(crosstab_rel_binary[hclust$labels,])
x<-x[,unique(GTDB_taxonomy$Family)]
row_dend  <-hclust %>% as.dendrogram %>%
  set("branches_k_color", k = 10) %>% set("branches_lwd", c(0.5,0.5)) 

heatmaply(x, Rowv = row_dend, Colv = F,fontsize_row = 6,
          grid_color = 'white', grid_size = 0.01 ,
          colors = c('grey','darkgreen'))

library(pheatmap)

svglite::svglite('clusters_pheatmap_only.svg',width = 10, height=11)
pdf('clusters_pheatmap_bw.pdf',width = 10, height=11)
ph<-pheatmap(x, cluster_rows = hclust,cluster_cols = FALSE, cellwidth = 10, cellheight = 7,
         border_color = 'grey', labels_row = rep('',nrow(x)),
         color = c('white','black'))
ph$gtable$grobs[[1]]$gp <- gpar(lex=3.5)
ph
dev.off()
library(stringr)

color_map<-data.frame('gene'=c('A','B','C','D','E','H','I','x'),
                      'color'=c('red','orange','yellow','green','dodgerblue1','purple','violet','grey'))
library(genoPlotR)
clusters<-list()
for( i in hclust$labels[hclust$order]){
  names=unlist(str_extract_all(i,'[ABCDEHIx]'))
  names_dir=unlist(str_extract_all(i,'...'))
  starts<-c()
  ends<-c()
  strands<-c()
  cols<-c()
  for(n in c(1:(length(names)))){
    starts<-c(starts,n*-5)
    ends<-c(ends,(n+1)*-5)
    strands<-c(strands,ifelse(grepl('r',names_dir[n]),1,-1 ) )
    cols<-c(cols, toString(color_map[color_map$gene== names[n],]$color) )
  }
  clusters[[i]]<-dna_seg(data.frame(name=names, start=starts, end=ends,
                                          strand=strands, fill=cols, col='black')) 
}
cluster_tree<- hclust2phylog(hclust)
plot_gene_map(dna_segs=clusters, tree=cluster_tree,tree_branch_labels_cex=NULL,
              dna_seg_labels=NULL,
              seg_plot_height=70)

mid_pos <- middle(clusters[[1]])

annots<-list()
for(i in c(1:length(clusters))){
  mid_pos <- middle(clusters[[i]])
  annots[[i]]<- annotation(x1=mid_pos,text=paste(ifelse(clusters[[i]]$name !='x',clusters[[i]]$name,'')  ))
}


svglite::svglite('genoplotr_cluster_arrows.svg', width=10, height =22)
plot_gene_map(dna_segs = clusters, annotations = annots,tree=cluster_tree,annotation_height = 0.2,annotation_cex = 0.7)
dev.off()

pdf('test_genoplotr.pdf', width=10,height = 20)
plot_gene_map(dna_segs=clusters,tree=cluster_tree,seg_plot_height=5,annotation_height=5,fixed_gap_length=0.1)
dev.off()


cluster_tree$leaves
names(clusters)[1] == names(cluster_tree$leaves)[1]
plot(cluster_tree)

names(clusters)<-c(1:length(labels(clusters)))
cluster_tree$leaves<-as.character(c(1:length(labels(clusters))))




cluster_tree$tre

plot(cluster_tree)

sort(colSums(x))

unique(gtdb$family_gtdb)
library(pheatmap)
library(heatmaply)
crosstab_filtered<- crosstab[rowSums(crosstab)>=10,]
heatmaply(log10(crosstab_filtered+1),fontsize_row = 5)
mat<-mat2
mat_names<- distmat_cluster2[,1] 
rownames(mat)<- mat_names
colnames(mat)<- mat_names
hclust<-hclust(as.dist(mat), method='average')
dend<-as.dendrogram(hclust)
crosstab_rel<-crosstab_rel[hclust$labels[hclust$order],]
crosstab_rel_binary<-ifelse(crosstab_rel>0,1,0)
heatmaply(crosstab_rel,fontsize_row = 6, dendrogram = 'none' ,grid_color = 'white', grid_size = 0.01 )
heatmaply(crosstab_rel_binary,fontsize_row = 6,
          grid_color = 'white', grid_size = 0.01 ,colors = c('grey','darkgreen'),
          Rowv = dend, Colv= F)

x<-data.frame(table(gtdb$synteny_dir_nhbr))

hclust$labels

#compare cluster distance to IacD distance


#importing IacD distance matrix----
distmat <-read.table('/Users/tslaird/Box/leveau_lab/Analysis_no_pfam_filter/IacD_analysis/IacD_sequences_GTDB_09132019_NO_Sphingomonas_IacD1_.fa.trim.msa.uncorrected.distmat', fill = TRUE, skip = 7, sep ='\t')
distmat<-distmat[-c(1,ncol(distmat)-1)]
distmat[distmat==""]<- NA
rownames(distmat)<- distmat[1:nrow(distmat),ncol(distmat)]
colnames(distmat)<- distmat[1:nrow(distmat),ncol(distmat)]
distmat<-distmat[-c(ncol(distmat))]
distmat_tri<-as.dist(t(distmat))

df<- read.csv('/Users/tslaird/Box/leveau_lab/FARM_output/iac_search_09122019_no_pfam_filter/results_09122019_nopfam_filter/iac_positive_all_data.tsv', sep = '\t')
synteny_dist<-distmat_cluster2[,2:ncol(distmat_cluster2)]
synteny_clust<-hclust(as.dist(synteny_dist), method='average')
rownames(synteny_dist)<-distmat_cluster2[,1]
colnames(synteny_dist)<-distmat_cluster2[,1]

avg_dists<-c()
for(i in synteny_clust$labels[synteny_clust$order]){
  assemb<-(gsub('\\.','_',df[df$synteny_dir_nhbr == i,]$assembly_x))
  index<-c()
  for(a in assemb){
    index<-c(index,which(grepl(a,rownames(distmat))))
  }
  mat_subset<-distmat[index,index]
  if(!is.null(dim(mat_subset))){
  diag(mat_subset)<-NA
  }
  print(i)
  avg_dists<-c(avg_dists,c(na.omit(unlist(mat_subset))))
}

boxplot(avg_dists)\


names(avg_dists)<-synteny_clust$labels[synteny_clust$order]
#regression of cluster distance

synteny_subset_df<-df[df$synteny_dir_nhbr %in% rownames(synteny_dist),]
sampled_assemb<-synteny_subset_df[sample(nrow(synteny_subset_df), 30), ]$assembly_x
index<-c()
for(a in gsub('\\.','_',sampled_assemb)){
  index<-c(index,which(grepl(a,rownames(distmat))))
}
index<-sort(index)
length(index)
mat_subset<-distmat[index,index]
z<-dist2list(as.dist(t(mat_subset)))
z[,c(1,2)]<-apply(z[,c(1,2)],1, function(x) sort(x[c(1,2)]))
z<-z[z$col != z$row,]
z<-z[!duplicated(z[,c(1,2)] ),]
z$g1<- gsub('(?<=\\d)_','\\.',sapply(z$col, function(x) str_extract(x,'GCF_.+?_\\d')), perl = T)
z$g2<- gsub('(?<=\\d)_','\\.',sapply(z$row, function(x) str_extract(x,'GCF_.+?_\\d')), perl = T)

z$g1_syn<- sapply(z$g1, function(x) toString(synteny_subset_df[synteny_subset_df['assembly_x'] == x,]['synteny_dir_nhbr'][1,1] ) )
z$g2_syn<- sapply(z$g2, function(x) synteny_subset_df[synteny_subset_df['assembly_x'] == x,]['synteny_dir_nhbr'][1,1]  )
z$syn_dist<- apply(z,1, function(x) synteny_dist[x['g1_syn'],x['g2_syn']])

plot(z$value,z$syn_dist)
res=lm(z$syn_dist ~ z$value, data=z)
abline(res)

# test----
library(dendextend)

family_tree<-read.tree('/Users/tslaird/Box/leveau_lab/GTDB/bac120_msa_r89_selected_families.faa.tree')
plot(family_tree)
y<-lapply(unique(as.character(GTDB_meta[GTDB_meta$gtdb_taxonomy %in% gtdb$gtdb_tax,]$gtdb_taxonomy )), function(x) str_split(x,';'))

str_split(as.character(GTDB_meta$gtdb_taxonomy),';')[1:5]


paste0('f__',unique(GTDB_taxonomy$Family))
