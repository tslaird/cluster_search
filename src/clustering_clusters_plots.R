library(dendextend)


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

#start----
library(kmer)
library(stringi)
library(gtools)
library(ecodist)
library(stringr)
#get GTDB taxonomy order
cluster_df<-read.csv('/home/tslaird/leveau_lab/cluster_search/FARM_output/iac_positive_all_data.tsv', sep='\t', stringsAsFactors = T)

GTDB_taxonomy<-cluster_df[,colnames(cluster_df) %in% c('domain_gtdb','phylum_gtdb','class_gtdb','order_gtdb','family_gtdb','genus_gtdb','species_gtdb')]
GTDB_taxonomy<-GTDB_taxonomy[GTDB_taxonomy$domain_gtdb!='',]
colnames(GTDB_taxonomy)<-c('Domain','Phylum','Class','Order','Family','Genus','Species')
GTDB_taxonomy$Domain<- factor(GTDB_taxonomy$Domain,levels(GTDB_taxonomy$Domain)[order(table(GTDB_taxonomy$Domain), decreasing = F)])
GTDB_taxonomy$Phylum<- factor(GTDB_taxonomy$Phylum, levels(GTDB_taxonomy$Phylum)[order(table(GTDB_taxonomy$Phylum), decreasing = F)])
GTDB_taxonomy$Class<- factor(GTDB_taxonomy$Class, levels(GTDB_taxonomy$Class)[order(table(GTDB_taxonomy$Class), decreasing = F)])
GTDB_taxonomy$Order<- factor(GTDB_taxonomy$Order,levels(GTDB_taxonomy$Order)[order(table(GTDB_taxonomy$Order), decreasing = F)])
#specially make Burkholderia come after
#GTDB_taxonomy$Order<- factor(GTDB_taxonomy$Order,c("Streptomycetales","Actinomycetales","Azospirillales",
#  "Propionibacteriales","Rhizobiales","Acetobacterales" ,"Rhodobacterales","Sphingomonadales",
#  "Mycobacteriales" , "Burkholderiales" ,"Xanthomonadales", "Enterobacterales", "Pseudomonadales"))
#
GTDB_taxonomy$Family<- factor(GTDB_taxonomy$Family,levels(GTDB_taxonomy$Family)[order(table(GTDB_taxonomy$Family), decreasing = F)])
GTDB_taxonomy$Genus<- factor(GTDB_taxonomy$Genus,levels(GTDB_taxonomy$Genus)[order(table(GTDB_taxonomy$Genus), decreasing = F)])
GTDB_taxonomy$Species<- factor(GTDB_taxonomy$Species,levels(GTDB_taxonomy$Species)[order(table(GTDB_taxonomy$Species), decreasing = F)])

GTDB_taxonomy$Domain_freq <-as.numeric(GTDB_taxonomy$Domain)
GTDB_taxonomy$Phylum_freq <-as.numeric(GTDB_taxonomy$Phylum)
GTDB_taxonomy$Class_freq <-as.numeric(GTDB_taxonomy$Class)
GTDB_taxonomy$Order_freq <-as.numeric(GTDB_taxonomy$Order)
GTDB_taxonomy$Family_freq <-as.numeric(GTDB_taxonomy$Family)
GTDB_taxonomy$Genus_freq <-as.numeric(GTDB_taxonomy$Genus)
GTDB_taxonomy$Species_freq <-as.numeric(GTDB_taxonomy$Species)

GTDB_taxonomy<-GTDB_taxonomy[ order( GTDB_taxonomy$Domain_freq,
                                     GTDB_taxonomy$Phylum_freq,
                                     GTDB_taxonomy$Class_freq,
                                     GTDB_taxonomy$Order_freq,
                                     GTDB_taxonomy$Family_freq,
                                     GTDB_taxonomy$Genus_freq,
                                     GTDB_taxonomy$Species_freq,
                                     decreasing = T),]


GTDB_taxonomy<-GTDB_taxonomy[,c(1:6)]
GTDB_taxonomy<-data.frame(lapply(GTDB_taxonomy, function(x) gsub('-','_',x)))

gtdb<- cluster_df[cluster_df$gtdb_tax !='' & !sapply(cluster_df$nhbrhood_prot_ids,function(x) grepl(x,'nan')),]

synteny_dist<-adist(unique(gtdb$synteny_alphabet_nhb), costs =  list(ins=1, del=10, sub=5) )

combos<-permutations(8, 2, c('A','B','C','D','E','H','I','X'), repeats = T)[1:63,]
combos_str<-paste0(combos[,1],combos[,2])
combos_counts<-t(sapply(unique(gtdb$synteny_alphabet_nhb), function(y) sapply(combos_str, function(x) str_count(y,x) + str_count(stri_reverse(y),x) )))
synteny_dist_kmer<-distance(combos_counts, method='jaccard')
synteny_dist<-as.matrix(synteny_dist_kmer)

methods<-c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")
correlation<-c()
for(i in seq(1:length(methods))){
  hc<- hclust(as.dist(synteny_dist),method=methods[i])
  hc_coph<-cophenetic(hc)
  correlation[i]<- cor(as.dist(synteny_dist_kmer),hc_coph)
}
results<-cbind.data.frame(methods,correlation)
print(results[order(-results$correlation),])
colnames(synteny_dist)<-unique(gtdb$synteny_alphabet_nhb)
rownames(synteny_dist)<-unique(gtdb$synteny_alphabet_nhb)
hclust<-hclust(as.dist(synteny_dist), method='average')
plot(hclust(as.dist(synteny_dist)))

crosstab<-as.data.frame.matrix(table(gtdb$synteny_alphabet_nhbr,gtdb$family_gtdb))[-1]
colnames(crosstab)<-gsub('f__','',colnames(crosstab))
#crosstab_rel<-t(apply(crosstab,1,function(x) x/sum(x)))
#hclust<-hclust(as.dist(synteny_dist))
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

svglite::svglite('cluster_of_clusters_pheatmap_only.svg',width = 10, height=11)
pdf('clusters_pheatmap_bw.pdf',width = 10, height=15)
ph<-pheatmap(x, cluster_rows = hclust,cluster_cols = FALSE, cellwidth = 8, cellheight = 5,
         border_color = 'grey', labels_row = hclust$labels,fontsize_row = 6.5,
         color = c('white','black'))
ph$gtable$grobs[[1]]$gp <- gpar(lex=3.5)
ph
dev.off()

color_map<-data.frame('gene'=c('A','B','C','D','E','H','I','X'),
                      'color'=c('red','orange','yellow','green','dodgerblue1','purple','violet','grey'))

clusters<-list()
for( i in hclust$labels[hclust$order]){
  names=unlist(str_extract_all(i,'[AaBbCcDdEeHhIiXx]'))
  names_upper=toupper(names)
  starts<-c()
  ends<-c()
  strands<-c()
  cols<-c()
  for(n in c(1:(length(names)))){
    starts<-c(starts,n*5)
    ends<-c(ends,(n+1)*5)
    strands<-c(strands,ifelse(grepl("^[[:upper:]]+$",names[n]),1,-1 ) )
    cols<-c(cols, toString(color_map[color_map$gene== names_upper[n],]$color) )
  }
  clusters[[i]]<-dna_seg(data.frame(name=names_upper, start=starts, end=ends,
                                          strand=strands, fill=cols, col='black')) 
}


#function----
arrow_plot<- function(input,start=1, left=15,right=25,total_length=100, arrow.width=3,horizontal_x=0, close.gene.setback=5,font_cex=1, center_on=NULL){
  #lines( c(min(as.numeric(input$start)) ,max(as.numeric(input$end))) ,c(0,0) )
  if(!is.null(center_on)){
    setback<- input[input$name==center_on,]$start[1]
    input$start<-input$start-setback
    input$end<-input$end-setback
  }

  apply(input, 1, function(x){ 
    mean_point<-((as.numeric(x['start'])+as.numeric(x['end']))/2)
    span<-as.numeric(x['end'])-as.numeric(x['start'])
    if(x['strand']==-1){
      polygon(c(as.numeric(x['start'])+(span*0.2),x['start'],as.numeric(x['start'])+(span*0.2),x['end'],x['end'],as.numeric(x['start'])+(span*0.2)),
              c(-arrow.width+horizontal_x,horizontal_x,arrow.width+horizontal_x,arrow.width+horizontal_x,-arrow.width+horizontal_x,-arrow.width+horizontal_x),col=x['fill'])}
    else{
      polygon(c(x['start'],x['start'],as.numeric(x['end'])-(span*0.2),x['end'],as.numeric(x['end'])-(span*0.2),x['start']),
              c(-arrow.width+horizontal_x,arrow.width+horizontal_x,arrow.width+horizontal_x,horizontal_x,-arrow.width+horizontal_x,-arrow.width+horizontal_x),col=x['fill'])}
    for(i in seq(1:nrow(input))){
      len<-input[i,]$end-input[i,]$start
      mean_point<-((as.numeric(input[i,]$start)+as.numeric(input[i,]$end))/2)
      #mean_point_ahead<-((as.numeric(input[i+1,]$start)+as.numeric(input[i+1,]$end))/2)
      #mean_point_behind<-((as.numeric(input[i-1,]$start)+as.numeric(input[i-1,]$end))/2)
      #segments(c(mean_point),c(0),c(mean_point),c(arrow.width+2))
      #segments(c(mean_point),c(arrow.width+2),c(mean_point),c(arrow.width+5))
      text(c(mean_point),c((arrow.width-0.5*arrow.width)+horizontal_x),paste(input[i,]$name),cex=font_cex,pos=1,srt=0,offset = 0)
    }
    
   })
  }


svglite::svglite('cluster_of_clusters_arrows.svg', width=10, height =22)
plot(NA, xlim=c(-100,100), ylim=c(-10,1000),axes=F, xlab=NA, ylab=NA)
for( i in rev(c(1:length(clusters)))){
  arrow_plot(rev(clusters)[[i]], horizontal_x = i*10, center_on = 'D', font_cex = 0.4)
}
dev.off()
