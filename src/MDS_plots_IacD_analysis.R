library(plotly)
library(vegan)
library(ggtree)


#4mer distance matrix and tree
cluster_df<-read.csv('/home/tslaird/Downloads/FARM_out_061020/iac_positive_all_data.tsv', sep='\t', stringsAsFactors = F)
cluster_df<-cluster_df[cluster_df$domain_gtdb!='',]
four_mer_freq<-gsub('\\]|\\[|\\s,',"",cluster_df$four_mer_freq_cluster)
four_mer_freq_genome<-gsub('\\]|\\[|\\s,',"",cluster_df$four_mer_freq_genome )
four_mer_numeric<-sapply(four_mer_freq, function(x) as.numeric(strsplit(x,',')[[1]]) )
four_mer_numeric_genome<-sapply(four_mer_freq_genome, function(x) as.numeric(strsplit(x,',')[[1]]) )
colnames(four_mer_numeric)<-cluster_df$assembly
colnames(four_mer_numeric_genome)<-paste0(cluster_df$assembly,'_genome')
four_mer_all<-cbind(four_mer_numeric,four_mer_numeric_genome)

four_mer_dist<-dist(t(four_mer_numeric))
four_mer_mat<-as.matrix(four_mer_dist)

four_mer_dist_all<-dist(t(four_mer_all))
four_mer_mat<-as.matrix(four_mer_dist_all)
diag(four_mer_mat)<-1

four_mer_mat_masked<-four_mer_mat[colnames(four_mer_numeric), colnames(four_mer_numeric_genome)]

#closest_match<-cbind.data.frame(rownames(four_mer_mat_masked),(apply(four_mer_mat_masked,1, function(x) colnames(four_mer_mat_masked)[which.min(x)])))
closest_match<-cbind.data.frame(rownames(four_mer_mat_masked), t(apply(four_mer_mat_masked,1, function(x) names(sort(x)[1:5]))  ))

colnames(closest_match)<-c('cluster','genome1','genome2','genome3','genome4','genome5')
closest_match$cluster_tax<-sapply(closest_match$cluster, function(x) cluster_df[cluster_df$assembly==x,]$family_gtdb[1])
closest_match$genome1_tax<-sapply( gsub('_genome','',closest_match$genome1), function(x) cluster_df[cluster_df$assembly==x,]$genus_gtdb[1] )
closest_match$genome2_tax<-sapply( gsub('_genome','',closest_match$genome2), function(x) cluster_df[cluster_df$assembly==x,]$genus_gtdb[1] )
closest_match$genome3_tax<-sapply( gsub('_genome','',closest_match$genome3), function(x) cluster_df[cluster_df$assembly==x,]$genus_gtdbb[1] )
closest_match$genome4_tax<-sapply( gsub('_genome','',closest_match$genome4), function(x) cluster_df[cluster_df$assembly==x,]$genus_gtdb[1] )
closest_match$genome5_tax<-sapply( gsub('_genome','',closest_match$genome5), function(x) cluster_df[cluster_df$assembly==x,]$genus_gtdb[1] )

closest_match$cluster_tax<-sapply(closest_match$cluster, function(x) cluster_df[cluster_df$assembly==x,]$family_gtdb[1])
closest_match$synteny1<-sapply( gsub('_genome','',closest_match$genome1), function(x) cluster_df[cluster_df$assembly==x,]$synteny_alphabet_nhbr[1] )
closest_match$synteny2<-sapply( gsub('_genome','',closest_match$genome2), function(x) cluster_df[cluster_df$assembly==x,]$synteny_alphabet_nhbr[1] )
closest_match$synteny3<-sapply( gsub('_genome','',closest_match$genome3), function(x) cluster_df[cluster_df$assembly==x,]$synteny_alphabet_nhbr[1] )
closest_match$synteny4<-sapply( gsub('_genome','',closest_match$genome4), function(x) cluster_df[cluster_df$assembly==x,]$synteny_alphabet_nhbr[1] )
closest_match$synteny5<-sapply( gsub('_genome','',closest_match$genome5), function(x) cluster_df[cluster_df$assembly==x,]$synteny_alphabet_nhbr[1] )

closest_match$cluster_synteny<-sapply(closest_match$cluster, function(x) cluster_df[cluster_df$assembly==x,]$synteny_alphabet_nhbr[1])
closest_match$same_synteny<-apply( closest_match,1, function(x) sum( x == x['cluster_synteny'])-1)
closest_match$same_tax<-apply( closest_match,1, function(x) sum( x == x['cluster_tax'])-1)
closest_match$cluster_gtdb<-sapply(closest_match$cluster, function(x) cluster_df[cluster_df$assembly==x,]$species_gtdb[1])


closest_match$genome_gtdb<-sapply( gsub('_genome','',closest_match$genome), function(x) cluster_df[cluster_df$assembly==x,]$species_gtdb[1] )
closest_match_filtered<-closest_match[(closest_match$same_synteny>1 & closest_match$same_tax<3),]


cluster_df_10<-read.csv('/home/tslaird/leveau_lab/cluster_search/FARM_output/iac_positive_all_data_10_10.tsv',sep='\t')
mobilization<-cluster_df_10[grep('integrase|transposase|phage',cluster_df_10$nhbrhood_prot_name),]$assembly
length(mobilization)


closest_match_hgt<-closest_match[((closest_match$cluster_tax != closest_match$genome_tax) & (closest_match$cluster_synteny == closest_match$genome_synteny) ),]
closest_match_hgt<-closest_match[((closest_match$cluster_tax != closest_match$genome_tax)),]

closest_match_hgt_mobile<-closest_match_hgt[closest_match_hgt$cluster %in% mobilization,]


#comparing four_mer of clusters only
four_mer_dist<-dist(t(four_mer_numeric))
four_mer_mat<-as.matrix(four_mer_dist)

diag(four_mer_mat)<-1
closest_match_<-cbind.data.frame( rownames(four_mer_mat), t(apply(four_mer_mat,1, function(x) names(sort(x)[1:5]))  ))
colnames(closest_match)<-c('cluster','genome1','genome2','genome3','genome4','genome5')
closest_match_$genome1_tax<-sapply( gsub('_genome','',closest_match_$genome1), function(x) cluster_df[cluster_df$assembly==x,]$genus_gtdb[1] )
closest_match_$genome2_tax<-sapply( gsub('_genome','',closest_match_$genome2), function(x) cluster_df[cluster_df$assembly==x,]$genus_gtdb[1] )
closest_match_$genome3_tax<-sapply( gsub('_genome','',closest_match_$genome3), function(x) cluster_df[cluster_df$assembly==x,]$genus_gtdbb[1] )
closest_match_$genome4_tax<-sapply( gsub('_genome','',closest_match_$genome4), function(x) cluster_df[cluster_df$assembly==x,]$genus_gtdb[1] )
closest_match_$genome5_tax<-sapply( gsub('_genome','',closest_match_$genome5), function(x) cluster_df[cluster_df$assembly==x,]$genus_gtdb[1] )
closest_match_$cluster_tax<-sapply(closest_match$cluster, function(x) cluster_df[cluster_df$assembly==x,]$family_gtdb[1])




random_indices<-sample(1:ncol(four_mer_mat), 100, replace=FALSE)
four_mer_dist_subset<-four_mer_mat[random_indices,random_indices]

complete_genomes<-which(cluster_df$complete_genome=="True")
four_mer_dist_subset<-four_mer_mat[complete_genomes,complete_genomes]


h_clust<-hclust(as.dist(four_mer_dist_subset))
tree<-ggtree(h_clust, layout = 'rectangular')
tree_data<-data.frame('label' = h_clust$labels)
tree_data<-cbind.data.frame(tree_data,cluster_df[random_indices,])
tree_data<-cbind.data.frame(tree_data,cluster_df[complete_genomes,])
tree_data<-merge(tree_data,cluster_df_meta, by= )
ggtree(h_clust, layout ='circular',size=0.5) %<+% tree_data + 
  geom_tippoint(aes(color=factor(order_gtdb)),size=1) +
  geom_tiplab(aes(label=paste(synteny_alphabet_nhbr, species_gtdb)),size=2)

taxonomy<-cluster_df[,colnames(cluster_df) %in% c('domain_gtdb','phylum_gtdb','class_gtdb','order_gtdb','family_gtdb','genus_gtdb','species_gtdb')]

NMDS_2D<-metaMDS(four_mer_mat,k =2,try=1,trymax=1)
NMDS2D = data.frame(NMDS_1 = NMDS_2D$points[,1], NMDS_2 = NMDS_2D$points[,2]  )
NMDS2D<-cbind.data.frame(NMDS2D,taxonomy)
pal<-c("#4BD0E4","#FF8F3F","#BB0044")
p <- plotly::plot_ly()
NMDS2D_plot <- plotly::plot_ly()
NMDS2D_plot<-add_markers(NMDS2D_plot,data=NMDS2D, x = NMDS2D$NMDS_1, y = NMDS2D$NMDS_2,color = NMDS2D$class_gtdb,
                         colors=pal,
                         marker=list(line=list(color='black', width=0.5),size=8,opacity = 1 ),
                         hoverinfo = 'text',  text = paste(NMDS2D$species_gtdb, cluster_df$synteny_alphabet_nhbr))
NMDS2D_plot<-add_markers(NMDS2D_plot, x = NMDS2D$NMDS_1[grep('Pseudomonas',NMDS2D$genus_gtdb)], y = NMDS2D$NMDS_2[grep('Pseudomonas',NMDS2D$genus_gtdb)], data = NMDS2D,inherit = TRUE,
                         color = NMDS2D$class_gtdb[grep('Pseudomonas',NMDS2D$genus_gtdb)],colors=pal, text=NMDS2D[grep('Pseudomonas',NMDS2D$genus_gtdb),]$species_gtdb,
                         marker=list(line=list(color='yellow', width=0.8),size=8,opacity = 1 ))
NMDS2D_plot<-add_markers(NMDS2D_plot, x = NMDS2D$NMDS_1[grep('Halomonas',NMDS2D$genus_gtdb)], y = NMDS2D$NMDS_2[grep('Halomonas',NMDS2D$genus_gtdb)], data = NMDS2D,inherit = TRUE,
                         color = NMDS2D$class_gtdb[grep('Halomonas',NMDS2D$genus_gtdb)],colors=pal, text=rownames(NMDS2D[grep('Halomonas',NMDS2D$genus_gtdb),]),
                         marker=list(line=list(color='pink', width=0.8),size=8,opacity = 1 ))

NMDS2D_plot

#3D
NMDS_3D<-metaMDS(four_mer_mat,k =3,try=1,trymax=1)
NMDS3D = data.frame(NMDS_1 = NMDS_3D$points[,1], NMDS_2 = NMDS_3D$points[,2], NMDS_3 = NMDS_3D$points[,3]  )
NMDS3D<-cbind.data.frame(NMDS3D,taxonomy)
pal<-c("#4BD0E4","#FF8F3F","#BB0044")
p <- plotly::plot_ly()
NMDS3D_plot <- plotly::plot_ly()
NMDS3D_plot<-add_markers(NMDS3D_plot,data=NMDS3D, x = NMDS3D$NMDS_1, y = NMDS3D$NMDS_2, z=NMDS3D$NMDS_3, color = NMDS3D$class_gtdb,
                         colors=pal,
                         marker=list(line=list(color='black', width=0.5),size=3,opacity = 1 ),
                         hoverinfo = 'text',  text = paste(NMDS3D$species_gtdb, cluster_df$synteny_alphabet_nhbr, cluster_df$assembly))
NMDS3D_plot<-add_markers(NMDS3D_plot, x = NMDS3D$NMDS_1[grep('Pseudomonas',NMDS3D$genus_gtdb)], y = NMDS3D$NMDS_2[grep('Pseudomonas',NMDS3D$genus_gtdb)], z= NMDS3D$NMDS_3[grep('Pseudomonas',NMDS3D$genus_gtdb)],data = NMDS3D,inherit = TRUE,
                         color = NMDS3D$class_gtdb[grep('Pseudomonas',NMDS3D$genus_gtdb)],colors=pal, text=rownames(NMDS3D[grep('Pseudomonas',NMDS3D$genus_gtdb),]),
                         marker=list(line=list(color='lime', width=0.8),size=3,opacity = 1 ))
NMDS3D_plot

#implicit phylogeny
comb<-t(combn(colnames(four_mer_mat)[1:200],2))
z<-t(apply(comb,1, function(x) {
c=four_mer_mat[x[1],x[2]]
d=four_mer_mat[paste0(x[1],'_genome'),paste0(x[2],'_genome')]
return(cbind(c,d)) }))
imphy_df<-cbind.data.frame(comb,z)
colnames(imphy_df)<-c('g1','g2','dc','dg')

imphy_df$g1_tax<-sapply(imphy_df$g1, function(x) cluster_df[cluster_df$assembly==as.character(x),]$gtdb_tax)
imphy_df$g2_tax<-sapply(imphy_df$g1, function(x) cluster_df[cluster_df$assembly==as.character(x),]$gtdb_tax)

imphy_plot<-plotly::plot_ly()
imphy_plot<- add_markers(imphy_plot,data=imphy_df, x=imphy_df$dc, y=imphy_df$dg, text=paste(imphy_df$g1,imphy_df$g2 ))
imphy_plot

#four mer distance of ABICDE gene cluster
#get assemblys with ABICDE synteny
ABICDE_assemblies<-cluster_df[cluster_df$synteny_alphabet_nhbr=='ABICDE',]$assembly
indices<-colnames(four_mer_mat) %in% ABICDE_assemblies
ABICDE_mat<-four_mer_mat[indices,indices]
h_clust<-hclust(as.dist(four_mer_mat[indices,indices]))

tree_data<-cbind.data.frame(tree_data,cluster_df[indices,])
tree_data<-data.frame('label' = h_clust$labels)
ggtree(h_clust, layout ='rectangular',size=1) %<+% tree_data + 
  geom_tippoint(aes(color=factor(class_gtdb)),size=1,shape=3)+
  geom_tiplab(aes(label=paste(synteny_alphabet_nhbr, species_gtdb)),size=2)+
  geom_treescale()



#IacD----
#importing distance matrix
distmat <-read.table('/home/tslaird/Downloads/FARM_out/IacD_analysis/IacD_seqs.fa.trim.msa.uncorrected.distmat', fill = TRUE, skip = 7, sep ='\t')
distmat<-distmat[-c(1,ncol(distmat)-1)]
distmat[distmat==""]<- NA
rownames(distmat)<- distmat[1:nrow(distmat),ncol(distmat)]
colnames(distmat)<- distmat[1:nrow(distmat),ncol(distmat)]
distmat<-distmat[-c(ncol(distmat))]
distmat_tri<-as.dist(t(distmat))

#assign group colors IacD----
library(stringr)
library(randomcoloR)
#gtdb_tax<-read.table('/Users/tslaird/Box/leveau_lab/iac_database/iac_positive_gtdb_taxonomy.tsv',sep='\t', header = T, fill = T)
accessions<-unlist(lapply(rownames(distmat), function(x) str_extract(x,'GCF_\\d+_\\d')))
#taxonomy<-lapply(accessions, function(x) cluster_df[cluster_df$assembly== gsub( "(?<=\\d)_(?=\\d)","\\.",x, perl = T),]  )
#taxonomy<-do.call('rbind', taxonomy)

taxonomy<-lapply(accessions, function(x) cluster_df[cluster_df$assembly_base == gsub( "(?<=\\d)_\\d","",x, perl = T),][c('assembly','domain_gtdb','phylum_gtdb','class_gtdb','order_gtdb','family_gtdb','genus_gtdb','species_gtdb')]  )
taxonomy<-do.call('rbind', taxonomy)
taxonomy<-taxonomy[!duplicated(taxonomy$assembly),]


#NMDS_3D IacD
NMDS_3D<-metaMDS(distmat_tri,k =3,try=5,trymax=5)
NMDS3D = data.frame(NMDS_1 = NMDS_3D$points[,1], NMDS_2 = NMDS_3D$points[,2], NMDS_3 = NMDS_3D$points[,3]  )
NMDS3D<-cbind.data.frame(NMDS3D,taxonomy)
pal<-c("#4BD0E4","#FF8F3F","#BB0044")
p <- plotly::plot_ly()
NMDS3D_plot <- plotly::plot_ly()
NMDS3D_plot<-add_markers(NMDS3D_plot,data=NMDS3D, x = NMDS3D$NMDS_1, y = NMDS3D$NMDS_2, z=NMDS3D$NMDS_3, color = NMDS3D$class_gtdb,
                         colors=pal,
                         marker=list(line=list(color='black', width=0.5),size=8,opacity = 1 ),
                         hoverinfo = 'text',  text = rownames(NMDS3D))
NMDS3D_plot<-add_markers(NMDS3D_plot, x = NMDS3D$NMDS_1[grep('Pseudomonas',NMDS3D$genus_gtdb)], y = NMDS3D$NMDS_2[grep('Pseudomonas',NMDS3D$genus_gtdb)], z= NMDS3D$NMDS_3[grep('Pseudomonas',NMDS3D$genus_gtdb)],data = NMDS3D,inherit = TRUE,
                         color = NMDS3D$class_gtdb[grep('Pseudomonas',NMDS3D$genus_gtdb)],colors=pal, text=rownames(NMDS3D[grep('Pseudomonas',NMDS3D$genus_gtdb),]),
                         marker=list(line=list(color='lime', width=0.8),size=8,opacity = 1 ))
NMDS3D_plot

#2D NMDS IacD----
NMDS_2D<-metaMDS(distmat_tri,k =2,try=1,trymax=1)
NMDS2D = data.frame(NMDS_1 = NMDS_2D$points[,1], NMDS_2 = NMDS_2D$points[,2], assembly=  gsub( "(?<=\\d)_(?=\\d)","\\.",accessions, perl = T) )
NMDS2D<-merge(NMDS2D,taxonomy, by='assembly')
NMDS2D<-merge(NMDS2D, cluster_df, by='assembly')

pal<-c("#4BD0E4","#FF8F3F","#BB0044")
plot(NMDS2D$NMDS_1,NMDS2D$NMDS_2)
p <- plotly::plot_ly()
NMDS2D_plot <- plotly::plot_ly()
NMDS2D_plot<-add_markers(NMDS2D_plot,data=NMDS2D, x = NMDS2D$NMDS_1, y = NMDS2D$NMDS_2,color = NMDS2D$class_gtdb.x,
                         colors=pal,
                         marker=list(line=list(color='black', width=0.5),size=8,opacity = 1 ),
                         hoverinfo = 'text',  text = paste(NMDS2D$species_gtdb, NMDS2D$synteny_alphabet_nhbr))
NMDS2D_plot<-add_markers(NMDS2D_plot, x = NMDS2D$NMDS_1[grep('Pseudomonas',NMDS2D$genus_gtdb)], y = NMDS2D$NMDS_2[grep('Pseudomonas',NMDS2D$genus_gtdb)], data = NMDS2D,inherit = TRUE,
                         color = NMDS2D$class_gtdb[grep('Pseudomonas',NMDS2D$genus_gtdb)],colors=pal, text=NMDS2D[grep('Pseudomonas',NMDS2D$genus_gtdb),]$species_gtdb,
                         marker=list(line=list(color='yellow', width=0.8),size=8,opacity = 1 ))
NMDS2D_plot<-add_markers(NMDS2D_plot, x = NMDS2D$NMDS_1[grep('Halomonas',NMDS2D$genus_gtdb)], y = NMDS2D$NMDS_2[grep('Halomonas',NMDS2D$genus_gtdb)], data = NMDS2D,inherit = TRUE,
                         color = NMDS2D$class_gtdb[grep('Halomonas',NMDS2D$genus_gtdb)],colors=pal, text=rownames(NMDS2D[grep('Halomonas',NMDS2D$genus_gtdb),]),
                         marker=list(line=list(color='pink', width=0.8),size=8,opacity = 1 ))

NMDS2D_plot

#ggplot NMDS2D color by class
NMDS2Dgg<-ggplot(NMDS2D, aes(x=NMDS_1, y=NMDS_2)) + 
  geom_point(aes(fill=class_gtdb, color=class_gtdb,shape=class_gtdb), size=2)+
  theme_bw()+
  scale_color_manual(name='GTDB classification',labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria','test'),values=c("black","black","black",'pink'))+
  scale_fill_manual(name='GTDB classification',labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria','test'),values=c("#4BD0E4","#FF8F3F","#BB0044",'#BB0044'))+
  scale_shape_manual(name='GTDB classification',labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria','test'),values=c(21, 22, 24,24))+
  xlab("NMDS1") + ylab("NMDS2")+
  geom_point(data= NMDS2D[grep('Pseudomonas',NMDS2D$genus_gtdb),],
             aes(x=NMDS2D$NMDS_1[grep('Pseudomonas',NMDS2D$genus_gtdb)], y=NMDS2D$NMDS_2[grep('Pseudomonas',NMDS2D$genus_gtdb)]),
                 fill= '#BB0044', color ='pink',pch=24,size=2)+
  geom_point(data= NMDS2D[grep('Halomonas',NMDS2D$genus_gtdb),],
             aes(x=NMDS2D$NMDS_1[grep('Halomonas',NMDS2D$genus_gtdb)], y=NMDS2D$NMDS_2[grep('Halomonas',NMDS2D$genus_gtdb)],
                 fill= NMDS2D$class_gtdb[grep('Halomonas',NMDS2D$genus_gtdb)]),color ='yellow',pch=24,size=2,
             show.legend = FALSE)
  
NMDS2Dgg
ggsave(file="/home/tslaird/leveau_lab/cluster_search/Figures/IacD_NMDS.svg", plot=NMDS2Dgg, width=5, height=3, dpi=600 )



guides(color = guide_legend( override.aes = list(name='GTDB classification',
                            labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria','test'),
                            fill=c("#4BD0E4","#FF8F3F","#BB0044",'#BB0044'),
                            color='black','black','black,','pink') ) )

NMDS2Dgg +
geom_point(data= NMDS2D[grep('Halomonas',NMDS2D$genus_gtdb),],
           aes(x=NMDS2D$NMDS_1[grep('Halomonas',NMDS2D$genus_gtdb)], y=NMDS2D$NMDS_2[grep('Halomonas',NMDS2D$genus_gtdb)],
               fill= NMDS2D$class_gtdb[grep('Halomonas',NMDS2D$genus_gtdb)]),color ='yellow',pch=24,size=2)

  
geom_point(data= NMDS2D[grep('Pseudomonas',NMDS2D$genus_gtdb),],
             aes(x=NMDS2D$NMDS_1[grep('Pseudomonas',NMDS2D$genus_gtdb)], y=NMDS2D$NMDS_2[grep('Pseudomonas',NMDS2D$genus_gtdb)],
                 fill= NMDS2D$class_gtdb[grep('Pseudomonas',NMDS2D$genus_gtdb)]),color ='pink',pch=24,size=2,
             show.legend = FALSE)+
geom_point(data= NMDS2D[grep('Halomonas',NMDS2D$genus_gtdb),],
             aes(x=NMDS2D$NMDS_1[grep('Halomonas',NMDS2D$genus_gtdb)], y=NMDS2D$NMDS_2[grep('Halomonas',NMDS2D$genus_gtdb)],
                 fill= NMDS2D$class_gtdb[grep('Halomonas',NMDS2D$genus_gtdb)]),color ='yellow',pch=24,size=2,
             show.legend = FALSE)+
  #stat_ellipse(aes(color=NMDS2D$class_gtdb, values=c("#427df6","#00a94b","#c20049")),type = "norm")+
  scale_fill_manual(name='GTDB classification', labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria'),values=c("#4BD0E4","#FF8F3F","#BB0044"))+
  scale_shape_manual(name='GTDB classification',values=c(21, 22, 24))+
  xlab("NMDS1") + ylab("NMDS2")
#stress is 0.05906332
NMDS2Dgg
ggsave(file="/Users/tslaird/Box/leveau_lab/iac_database/IacD_analysis_101219/IacD_NMDS2D_ggplot.svg", plot=NMDS2Dgg, width=8, height=6)


#ggplot NMDS2D color by class highlight several

Pseudomonas <- NMDS2D[grep('Pseudomonas',NMDS2D$genus_gtdb),]
Halomonadacea<- NMDS2D[grep('Halomonadacea',NMDS2D$family_gtdb),]

NMDS2Dgg<- ggplot(NMDS2D, aes(x = NMDS_1, y =  NMDS_2, shape = class_gtdb, 
                fill = class_gtdb, color=class_gtdb))+  
  geom_point(size = 2,stroke = 0.3) +
  geom_point(data = Pseudomonas, aes(fill = "Pseudomonas", shape = "Pseudomonas",color="Pseudomonas"),size=2,stroke = 0.3) +
  geom_point(data = Halomonadacea, aes(fill = "Halomonadaceae", shape = "Halomonadaceae",color="Halomonadaceae"),size=2,stroke = 0.3) +
  scale_shape_manual(name='GTDB classification', labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria',"Halomonadaceae","Pseudomonas"),
                     values = c(21, 21, 21,21,21) ) +
  scale_fill_manual(name='GTDB classification', labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria',"Halomonadaceae","Pseudomonas"),
                    values =c("#4BD0E4","#FF8F3F","#BB0044",'#BB0044','#BB0044') )+
  scale_color_manual(name='GTDB classification', labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria',"Halomonadaceae","Pseudomonas"),
                     values =c('black',"black","black","green",'pink') )+
  theme_bw()
  
ggsave(file="/Users/tslaird/Box/leveau_lab/iac_database/IacD_analysis_101219/IacD_NMDS2D_ggplot_103019.svg", plot=NMDS2Dgg, width=8, height=6)
ggsave(file="/Users/tslaird/Box/leveau_lab/iac_database/IacD_analysis_101219/IacD_NMDS2D_ggplot_103019.pdf", plot=NMDS2Dgg, width=5, height=3, dpi=600 )

####

NMDS2D$ABICDE<- ifelse(NMDS2D$assembly_x %in% gtdb[gtdb$synteny_dir_nhbr =='JArJBrJIrJCrJDrJEr','assembly_x'],'Yes','No')

NMDS2Dgg_synteny<- ggplot(NMDS2D, aes(x = NMDS_1, y =  NMDS_2, shape = ABICDE, 
                              fill = class_gtdb, color=ABICDE, size=ABICDE ))+  
  geom_point(stroke=0.4) +
  scale_shape_manual(name='ABICDE synteny', labels=c('No','Yes'),
                     values = c(21, 24)) +
  scale_fill_manual(name='GTDB classification', labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria'),
                    values =c("#4BD0E4","#FF8F3F","#BB0044") )+
  scale_color_manual(name='ABICDE synteny', labels=c('No','Yes'),
                     values =c("#3d3d3d", 'black'))+
  scale_size_manual(name='ABICDE synteny', labels=c('No','Yes'),
                     values =c(1.8,2))+
  theme_bw()

NMDS2Dgg_synteny

ggsave(file="/Users/tslaird/Box/leveau_lab/iac_database/IacD_analysis_101219/IacD_NMDS2D_ggplot_synteny_updated102519.svg", plot=NMDS2Dgg_synteny, width=8, height=6)
ggsave(file="/Users/tslaird/Box/leveau_lab/iac_database/IacD_analysis_101219/IacD_NMDS2D_ggplot_synteny_updated102519.pdf", plot=NMDS2Dgg_synteny, width=5, height=3, dpi=600 )

###


#ggplot color by order
NMDS2Dgg<-ggplot(NMDS2D, aes(x=NMDS2D$NMDS_1, y=NMDS2D$NMDS_2)) + 
  geom_point(aes(fill=NMDS2D$order_gtdb),color="black",pch=21,size=2 )+
  theme_bw()+
  geom_point(data= NMDS2D[grep('Pseudomonas',NMDS2D$genus_gtdb),],
             aes(x=NMDS2D$NMDS_1[grep('Pseudomonas',NMDS2D$genus_gtdb)], y=NMDS2D$NMDS_2[grep('Pseudomonas',NMDS2D$genus_gtdb)],
                 fill= NMDS2D$order_gtdb[grep('Pseudomonas',NMDS2D$genus_gtdb)]),color ='#e38c00',pch=21,size=2,
             show.legend = FALSE)+
  #stat_ellipse(aes(color=NMDS2D$class_gtdb, values=c("#427df6","#00a94b","#c20049")),type = "norm")+
  #scale_fill_manual(name = 'GTDB_class',values=c("#427df6","#00a94b","#c20049"))+
  xlab("NMDS1") + ylab("NMDS2")
#stress is 0.05906332
NMDS2Dgg
ggsave(file="/Users/tslaird/Box/leveau_lab/iac_database/IacD_analysis/IacD_NMDS2D_ggplot_colorbyorder.svg", plot=NMDS2Dgg, width=8, height=6)



#pcoa 2D----
pcoa2D<- cmdscale(four_mer_mat ,k=2, eig = TRUE,list. = TRUE)
pcoa2D<-as.data.frame(pcoa2D$points)
pcoa2D<- cbind.data.frame(pcoa2D,taxonomy)
#pal<- distinctColorPalette(length(unique(pcoa2D$class_gtdb)),runTsne = T)
pal<-c("#427df6",
        "#00a94b",
        "#c20049")
pcoa2D_plot <- plotly::plot_ly()
pcoa2D_plot<-add_markers(pcoa2D_plot,data=pcoa2D, x = pcoa2D$V1, y = pcoa2D$V2,color = pcoa2D$class_gtdb,
                     colors=pal,
                     marker=list(line=list(color='black', width=0.5),size=8,opacity = 1 ),
                     hoverinfo = 'text',  text = pcoa2D$species_gtdb)
pcoa2D_plot<-add_markers(pcoa2D_plot, x = pcoa2D$V1[grep('Pseudomonas',pcoa2D$genus_gtdb)], y = pcoa2D$V2[grep('Pseudomonas',pcoa2D$genus_gtdb)], data = pcoa2D,inherit = TRUE,
                         color = pcoa2D$class_gtdb[grep('Pseudomonas',pcoa2D$genus_gtdb)],colors=pal, text=rownames(pcoa2D[grep('Pseudomonas',pcoa2D$genus_gtdb),]),
                         marker=list(line=list(color='orange', width=0.8),size=8,opacity = 1 ))
pcoa2D_plot

#3D separating Pseudomonas
pcoa3D_plot <- plotly::plot_ly()
pcoa3D_plot<-add_markers(pcoa3D_plot,data=pcoa3D, x = pcoa3D$V1, y = pcoa3D$V2,color = pcoa3D$class_gtdb,z=pcoa3D$V3,
                         colors=pal,
                         marker=list(line=list(color='black', width=0.5),size=7,opacity = 1 ),
                         hoverinfo = 'text',  text = rownames(distmat))
pcoa3D_plot<-add_markers(pcoa3D_plot, x = pcoa3D$V1[grep('Pseudomonas',pcoa3D$genus_gtdb)], y = pcoa3D$V2[grep('Pseudomonas',pcoa3D$genus_gtdb)], z= pcoa3D$V3[grep('Pseudomonas',pcoa3D$genus_gtdb)],
                         data = pcoa3D,inherit = TRUE,
                         color = pcoa3D$class_gtdb[grep('Pseudomonas',pcoa3D$genus_gtdb)],colors=pal, text=rownames(pcoa3D[grep('Pseudomonas',pcoa3D$genus_gtdb),]),
                         marker=list(line=list(color='orange', width=0.8),size=7,opacity = 1 ))
pcoa3D_plot

#pcoa 3D----
pcoa3D<- cmdscale(distmat_tri,k=3)
pcoa3D<-as.data.frame(pcoa3D)
pcoa3D<- cbind.data.frame(pcoa3D,taxonomy)
pal<- distinctColorPalette(length(unique(pcoa3D$order_gtdb)),runTsne = T)
pal<-c("#8961b3",
       "#969e4a",
       "#b8545c")
pcoa3D_plot <- plotly::plot_ly(pcoa3D, x = ~V1, y = ~V2,z= ~V3, color = pcoa3D$order_gtdb,
                     colors=pal,
                     marker=list(line=list(color="black", width=1.0),size=5,opacity = 1 ),
                     hoverinfo = 'text',  text = rownames(distmat))
pcoa3D_plot

#----
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
MDS3 = NMDS$points[,3]
MDS_plot = data.frame(MDS1 = MDS1, MDS2 = MDS2,MDS3 = MDS3 )
MDS_plot<-cbind.data.frame(MDS_plot,taxonomy)

set.seed(100)
pal<- distinctColorPalette(length(unique(MDS_plot$order_gtdb)))
p <- plotly::plot_ly(MDS_plot, x = ~MDS1, y = ~MDS2, z = ~MDS3,size = 2, color = ~order_gtdb,
                     colors=pal,
                     marker=list(line=list(color="black", width=0.5),opacity = 1 ),
                     hoverinfo = 'text',
                     text = rownames(distmat))
p

#PCoA
library(ape)
pcoa<- pcoa(D=distmat_tri, correction="lingoes")
biplot(x, Y=NULL, plot.axes = c(1,2), dir.axis1=1,
       dir.axis2=1, rn=NULL, main=NULL, ...)





#recA----
library(stringr)
library(vegan)
distmat_RecA <-read.table('/home/tslaird/Downloads/FARM_out/IacD_analysis/RecA_seqs.fa.trim.msa.uncorrected.distmat', fill = TRUE, skip = 7, sep ='\t')
distmat_RecA<-distmat_RecA[-c(1,ncol(distmat_RecA)-1)]
distmat_RecA[distmat_RecA==""]<- NA
rownames(distmat_RecA)<- distmat_RecA[1:nrow(distmat_RecA),ncol(distmat_RecA)]
colnames(distmat_RecA)<- distmat_RecA[1:nrow(distmat_RecA),ncol(distmat_RecA)]
distmat_RecA<-distmat_RecA[-c(ncol(distmat_RecA))]
accessions_RecA<-unlist(lapply(rownames(distmat_RecA), function(x) str_extract(x,'GCF_\\d+_\\d')))
keep<-which(accessions_RecA %in%  gsub( "(?<=\\d)\\.(?=\\d)","_",cluster_df$assembly, perl = T))
distmat_RecA<- distmat_RecA[c(keep),c(keep)]
distmat_RecA_tri<-as.dist(t(distmat_RecA))
NMDS_2D_RecA<-metaMDS(distmat_RecA_tri,k =2,try=1,trymax=1)

library(stringr)
library(randomcoloR)
library(ggplot2)
#gtdb_tax<-read.table('/Users/tslaird/Box/leveau_lab/iac_database/iac_positive_gtdb_taxonomy.tsv',sep='\t', header = T, fill = T)
accessions_RecA<-do.call('rbind', lapply(rownames(distmat_RecA), function(x) str_extract(x,'GCF_\\d+_\\d')))
taxonomy_RecA<-lapply(accessions_RecA, function(x) cluster_df[cluster_df$assembly_base == gsub( "(?<=\\d)_\\d","",x, perl = T),][c('assembly','domain_gtdb','phylum_gtdb','class_gtdb','order_gtdb','family_gtdb','genus_gtdb','species_gtdb')]  )
taxonomy_RecA<-do.call('rbind', taxonomy_RecA)
taxonomy_RecA<-taxonomy_RecA[!duplicated(taxonomy_RecA$assembly),]


NMDS2D_RecA <- data.frame(NMDS_1 = NMDS_2D_RecA$points[,1], NMDS_2 = NMDS_2D_RecA$points[,2]  )
NMDS2D_RecA<-cbind.data.frame(NMDS2D_RecA,taxonomy_RecA)

Pseudomonas <- NMDS2D_RecA[grep('Pseudomonas',NMDS2D_RecA$genus_gtdb),]
Halomonadacea<- NMDS2D_RecA[grep('Halomonadacea',NMDS2D_RecA$family_gtdb),]

##
NMDS2D_RecAgg<-ggplot(NMDS2D_RecA, aes(x=NMDS_1, y=NMDS_2)) + 
  geom_point(aes(fill=class_gtdb, color=class_gtdb,shape=class_gtdb), size=2)+
  theme_bw()+
  scale_color_manual(name='GTDB classification',labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria','test'),values=c("black","black","black",'pink'))+
  scale_fill_manual(name='GTDB classification',labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria','test'),values=c("#4BD0E4","#FF8F3F","#BB0044",'#BB0044'))+
  scale_shape_manual(name='GTDB classification',labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria','test'),values=c(21, 22, 24,24))+
  xlab("NMDS1") + ylab("NMDS2")+
  geom_point(data= NMDS2D_RecA[grep('Pseudomonas',NMDS2D_RecA$genus_gtdb),],
             aes(x=NMDS2D_RecA$NMDS_1[grep('Pseudomonas',NMDS2D_RecA$genus_gtdb)], y=NMDS2D_RecA$NMDS_2[grep('Pseudomonas',NMDS2D_RecA$genus_gtdb)]),
             fill= '#BB0044', color ='pink',pch=24,size=2)+
  geom_point(data= NMDS2D_RecA[grep('Halomonas',NMDS2D_RecA$genus_gtdb),],
             aes(x=NMDS2D_RecA$NMDS_1[grep('Halomonas',NMDS2D_RecA$genus_gtdb)], y=NMDS2D_RecA$NMDS_2[grep('Halomonas',NMDS2D_RecA$genus_gtdb)],
                 fill= NMDS2D_RecA$class_gtdb[grep('Halomonas',NMDS2D_RecA$genus_gtdb)]),color ='yellow',pch=24,size=2,
             show.legend = FALSE)

NMDS2D_RecAgg
ggsave(file="/home/tslaird/leveau_lab/cluster_search/Figures/RecA_NMDS.svg", plot=NMDS2D_RecAgg, width=5, height=3, dpi=600 )








NMDS2D_RecAgg<- ggplot(NMDS2D_RecA, aes(x = NMDS_1, y =  NMDS_2, shape = class_gtdb, 
                                        fill = class_gtdb, color=class_gtdb))+ 
  geom_point(size = 2,stroke = 0.3) +
  geom_point(data = Pseudomonas, aes(fill = "Pseudomonas", shape = "Pseudomonas",color="Pseudomonas"),stroke = 0.3) +
  geom_point(data = Halomonadacea, aes(fill = "Halomonadaceae", shape = "Halomonadaceae",color="Halomonadaceae"),stroke = 0.3) +
  scale_shape_manual(name='GTDB classification', labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria',"Halomonadaceae","Pseudomonas"),
                     values = c(21, 21, 21,21,21) ) +
  scale_fill_manual(name='GTDB classification', labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria',"Halomonadaceae","Pseudomonas"),
                    values =c("#4BD0E4","#FF8F3F","#BB0044",'#BB0044','#BB0044') )+
  scale_color_manual(name='GTDB classification', labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria',"Halomonadaceae","Pseudomonas"),
                     values =c('black',"black","black","green",'pink') )+
  theme_bw()
NMDS2D_RecAgg
ggsave(file="/Users/tslaird/Box/leveau_lab/iac_database/IacD_analysis_101219/RecA_NMDS2D_ggplot_103019.svg", plot=NMDS2D_RecAgg, width=8, height=6)
ggsave(file="/Users/tslaird/Box/leveau_lab/iac_database/IacD_analysis_101219/RecA_NMDS2D_ggplot_103019.pdf", plot=NMDS2D_RecAgg, width=5, height=3, dpi=600 )

test<- NMDS2D_RecA[grep('Streptomyces',NMDS2D_RecA$genus_gtdb),]



pal<-c("#427df6",
       "#00a94b",
       "#c20049")


p <- plotly::plot_ly()
NMDS2D_RecA_plot <- plotly::plot_ly()
NMDS2D_RecA_plot<-add_markers(NMDS2D_RecA_plot,data=NMDS2D_RecA, x = NMDS2D_RecA$NMDS_1, y = NMDS2D_RecA$NMDS_2,color = NMDS2D_RecA$class_gtdb,
                         colors=pal,
                         marker=list(line=list(color='black', width=0.5),size=8,opacity = 1 ),
                         hoverinfo = 'text',  text = rownames(NMDS2D_RecA))
NMDS2D_RecA_plot<-add_markers(NMDS2D_RecA_plot, x = NMDS2D_RecA$NMDS_1[grep('Pseudomonas',NMDS2D_RecA$genus_gtdb)], y = NMDS2D_RecA$NMDS_2[grep('Pseudomonas',NMDS2D_RecA$genus_gtdb)], data = NMDS2D_RecA,inherit = TRUE,
                         color = NMDS2D_RecA$class_gtdb[grep('Pseudomonas',NMDS2D_RecA$genus_gtdb)],colors=pal, text=rownames(NMDS2D_RecA[grep('Pseudomonas',NMDS2D_RecA$genus_gtdb),]),
                         marker=list(line=list(color='orange', width=0.8),size=8,opacity = 1 ))
NMDS2D_RecA_plot

NMDS2D_RecAgg<-ggplot(NMDS2D_RecA, aes(x=NMDS2D_RecA$NMDS_1, y=NMDS2D_RecA$NMDS_2)) + 
  geom_point(aes(fill=class_gtdb),color="black",pch=21,size=2 )+
  theme_bw()+
  geom_point(data= NMDS2D_RecA[grep('Pseudomonas',NMDS2D_RecA$genus_gtdb),],
             aes(x=NMDS2D_RecA$NMDS_1[grep('Pseudomonas',NMDS2D_RecA$genus_gtdb)], y=NMDS2D_RecA$NMDS_2[grep('Pseudomonas',NMDS2D_RecA$genus_gtdb)],
                 fill= NMDS2D_RecA$class_gtdb[grep('Pseudomonas',NMDS2D_RecA$genus_gtdb)]),color ='pink',pch=21,size=2,
             show.legend = FALSE)+
  #stat_ellipse(aes(color=NMDS2D_RecA$class_gtdb, values=c("#427df6","#00a94b","#c20049")),type = "norm")+
  scale_fill_manual(name = 'GTDB_class',values=c("#4BD0E4","#FF8F3F","#BB0044"))+
  xlab("NMDS1") + ylab("NMDS2")
NMDS2D_RecAgg

ggsave(file="/Users/tslaird/Box/leveau_lab/iac_database/IacD_analysis/RecA/RecA_NMDS2D_ggplot_09182019.svg", plot=NMDS2D_RecAgg, width=8, height=6)
ggsave(file="/Users/tslaird/Box/leveau_lab/iac_database/IacD_analysis/RecA/RecA_NMDS2D_ggplot_09182019.png", plot=NMDS2D_RecAgg,  width=5, height=3, dpi=600)


# recA 






# comparing IacD distance to cluster distance
library(stringr)
library(spaa)
distmat_flat<-dist2list(as.dist(distmat))
distmat_flat$col<-(str_match(distmat_flat$col, '(?<=!!_).+(?=_WP)'))


#IacB----

distmat_IacB <-read.table('/Users/tslaird/Box/leveau_lab/iac_database/IacB_analysis_101419/IacB_seqs_101419.fa.trim.msa.uncorrected.distmat', fill = TRUE, skip = 7, sep ='\t')
distmat_IacB<-distmat_IacB[-c(1,ncol(distmat_IacB)-1)]
distmat_IacB[distmat_IacB==""]<- NA
rownames(distmat_IacB)<- distmat_IacB[1:nrow(distmat_IacB),ncol(distmat_IacB)]
colnames(distmat_IacB)<- distmat_IacB[1:nrow(distmat_IacB),ncol(distmat_IacB)]
distmat_IacB<-distmat_IacB[-c(ncol(distmat_IacB))]
distmat_IacB_tri<-as.dist(t(distmat_IacB))

library(stringr)
library(randomcoloR)
gtdb_tax<-read.table('/Users/tslaird/Box/leveau_lab/iac_database/iac_positive_gtdb_taxonomy.tsv',sep='\t', header = T, fill = T)
accessions<-unlist(lapply(rownames(distmat_IacB), function(x) str_extract(x,'GCF_\\d+_\\d')))
taxonomy<-lapply(accessions, function(x) gtdb_tax[gtdb_tax$assembly_x== gsub( "(?<=\\d)_(?=\\d)","\\.",x, perl = T),]  )
taxonomy<-do.call('rbind', taxonomy)


NMDS_IacB_2D<-metaMDS(distmat_IacB_tri,k =2,try=1,trymax=1)
NMDS_IacB2D = data.frame(NMDS_1 = NMDS_IacB_2D$points[,1], NMDS_2 = NMDS_IacB_2D$points[,2]  )
NMDS_IacB2D<-cbind.data.frame(NMDS_IacB2D,taxonomy)

NMDS_IacB2D$ABICDE<- ifelse(NMDS_IacB2D$assembly_x %in% gtdb[gtdb$synteny_dir_nhbr =='JArJBrJIrJCrJDrJEr','assembly_x'],'Yes','No')

NMDS_IacB2Dgg_synteny<- ggplot(NMDS_IacB2D, aes(x = NMDS_1, y =  NMDS_2, shape = class_gtdb, 
                                      fill = class_gtdb, color=ABICDE, size=ABICDE ))+  
  geom_point(stroke=0.7) +
  scale_shape_manual(name='GTDB classification', labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria',"ABICDE"),
                     values = c(21, 24, 22)) +
  scale_fill_manual(name='GTDB classification', labels=c('Actinobacteria','Alphaproteobacteria','Gammaproteobacteria',"ABICDE"),
                    values =c("#4BD0E4","#FF8F3F","#BB0044") )+
  scale_color_manual(name='ABICDE synteny', labels=c('No','Yes'),
                     values =c("black", 'white'))+
  scale_size_manual(name='ABICDE synteny', labels=c('No','Yes'),
                    values =c(2,2))+
  theme_dark()

NMDS_IacB2Dgg_synteny

# phylogenetic tree IacD
library(ggtree)
IacD_tree=read.tree('/Users/tslaird/Box/leveau_lab/iac_database/IacD_analysis_101219/IacD_seqs_101219.fa.trim.msa.fasttree.tree')
accessions_tree<-unlist(lapply(IacD_tree$tip.label, function(x) str_extract(x,'GCF_\\d+_\\d')))
taxonomy_tree<-lapply(accessions, function(x) gtdb_tax[gtdb_tax$assembly_x== gsub( "(?<=\\d)_(?=\\d)","\\.",x, perl = T),]  )
taxonomy_tree<-do.call('rbind', taxonomy_tree)
rownames(taxonomy_tree) <- taxonomy_tree$tip.label


p<-ggtree(IacD_tree, layout = 'circular')

p1<-  gheatmap(IacD_tree, as.character(taxonomy_tree[,'class_gtdb', drop=F]), offset=.8, width=.1)

df <- data.frame("category" = rep(c("a", "b", "c", "d"), 3))

set.seed(102)
sim_tree <- rtree(48, rooted = FALSE)
ggsim <- ggtree(sim_tree)

df <- data.frame("category" = rep(c("a", "b", "c", "d"), 3))
rownames(taxonomy_tree) = IacD_tree$tip.label

gheatmap(p, taxonomy_tree[, "class_gtdb", drop = FALSE], width = 0.6, color = NA,colnames = FALSE) +
  scale_fill_manual(breaks = c("c__Gammaproteobacteria", "c__Alphaproteobacteria", "c__Actinobacteria"),
                    values = c("#4BD0E4","#FF8F3F","#BB0044"))
                    

