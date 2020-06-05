library(plotly)

hdist<-0.16
p <- plot_ly(
  type = "sankey",
  orientation = "h",
  arrangement='floating',
  
  node = list(
    label = nodes$name,
    color = 'darkblue',
    pad = 15,
    #x= c(0.1,(0.1+hdist),(0.1+hdist),(0.1+2*hdist)),
    #y= c(0.3,0.3,0.5,0.3),
    thickness = 20,
    line = list(
      color = "black",
      width = 0.5
    )
  ),
  
  link = list(
    source = links$source,
    target = links$target,
    value =  links$value
  )
) %>% 
  layout(
    title = "Basic Sankey Diagram",
    font = list(
      size = 8
    )
  )
p %>% layout(autosize = F)

links

#----
library(dplyr)
#GTDB_taxonomy<-read.csv('/Users/tslaird/Box/leveau_lab/iac_database/GTDB_table_iac_positive.tsv', sep='\t', header = T )
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

categories<-unique(unlist(GTDB_taxonomy))

library(networkD3)
library(reshape2)
library(htmlwidgets)

links<-data.frame()
for(i in c(1:(ncol(GTDB_taxonomy)))){
  if(i<ncol(GTDB_taxonomy)){
    links<-rbind.data.frame(links, melt(table(GTDB_taxonomy[,i],GTDB_taxonomy[,i+1])))
  }
}

colnames(links)<-c('source','target','value')


nodes<-data.frame('name' = unique(as.character(unlist(GTDB_taxonomy))))
node_counts<-sapply(nodes$name, function(x) sum(GTDB_taxonomy==as.character(x)))
nodes<-cbind.data.frame(nodes,node_counts)
nodes$combined_name<-paste0(nodes$name,' (',nodes$node_counts,')')

links$source<- sapply(links$source, function(x) which(nodes$name==as.character(x))-1)
links$target<- sapply(links$target, function(x) which(nodes$name==as.character(x))-1)
links<-links[links$value>0,]


library(colourvalues)
color_palettes()
color_vals<-colourvalues::colour_values(nodes$node_counts, palette = 'magma',summary=T)
color_vals<-toString(shQuote(color_vals))
color_vals<-gsub("'",'"',color_vals)

color_vals<-rep('#001f59',nrow(nodes))
color_vals<-toString(shQuote(color_vals))
color_vals<-gsub("'",'"',color_vals)


names(links) = c("source", "target", "value")
sn<-sankeyNetwork(Links = links, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name",
              colourScale = JS(paste0("d3.scaleOrdinal().range([",color_vals,"]);")),
              fontFamily = 'Arial Bold',nodePadding = 12,
              height = 1550, width=1800,
              fontSize= 14, nodeWidth = 50, sinksRight = F, iterations = 10000)
sn

for(i in seq(1,100)){
  rect(ncol(df)+3,sbr[i],ncol(df)+4,sbr[i]+1.1,col=as.character(col[i,1]),border = NA)
}
rect(ncol(df)+3,max(sbr)+1.1,ncol(df)+4,max(sbr)+2.2,col='white',border = NA)

text(ncol(df)+5,1,as.character(min(df)),cex=font.size,offset=0,pos=4)
text(ncol(df)+5,max(sbr)+1.1,as.character(max(df)),cex=font.size,offset=0,pos=4)

saveNetwork(sn, file='sankey_taxonomy_101419.html')

base::print(type)
svg('testsn.svg')
sn
dev.off()

onRender(
  sn,
  '
  function(el, x) {
  d3.selectAll(".node text").attr("text-anchor", "begin").attr("x", 20);
  }
  '
)


#make a csv----
GTDB_table<-as.data.frame(apply(GTDB_taxonomy,2,function(x) sapply(x, function(y) nodes$combined_name[which(nodes$name == y)]  )),stringsAsFactors = FALSE)

GTDB_table<-GTDB_table[!duplicated(GTDB_table), ]

GTDB_table$Domain[duplicated(GTDB_table$Domain)]<-''
GTDB_table$Phylum[duplicated(GTDB_table$Phylum)]<-''
GTDB_table$Class[duplicated(GTDB_table$Class)]<-''
GTDB_table$Order[duplicated(GTDB_table$Order)]<-''
GTDB_table$Family[duplicated(GTDB_table$Family)]<-''
GTDB_table$Genus[duplicated(GTDB_table$Genus)]<-''

write.table(GTDB_table,file='/Users/tslaird/Box/leveau_lab/iac_database/Figures/GTDB_table.tsv', sep='\t', row.names = F)
####


##Sankey for metadata----
library(readxl)
metadata<-read_excel('/Users/tslaird/Box/leveau_lab/iac_database/metadata_mapping_09132019.xlsx')
metadata<-data.frame( metadata[metadata$A != 'Unspecified',], stringsAsFactors=T)

metadata$A<- factor(metadata$A,levels(factor(metadata$A))[order(table(metadata$A), decreasing = F)])
metadata$B<- factor(metadata$B,levels(factor(metadata$B))[order(table(metadata$B), decreasing = F)])
metadata$C<- factor(metadata$C,levels(factor(metadata$C))[order(table(metadata$C), decreasing = F)])
metadata$class_gtdb<- factor(metadata$class_gtdb,levels(factor(metadata$class_gtdb))[order(table(metadata$class_gtdb), decreasing = F)])

metadata$A_freq <-as.numeric(metadata$A)
metadata$B_freq <-as.numeric(metadata$B)
metadata$class_gtdb_freq <-as.numeric(metadata$class_gtdb)


metadata<-metadata[ order( metadata$A_freq,
                            metadata$B_freq,
                           metadata$class_gtdb_freq,
                           decreasing = T), c('A','B','class_gtdb')]

nodes<-data.frame('name' = unique(as.character(unlist(metadata))))


links<-data.frame()
links<-rbind.data.frame(links, melt(table(metadata$A,metadata$B)))
links<-rbind.data.frame(links, melt(table(metadata$B,metadata$class_gtdb)))
colnames(links)<-c('source','target','value')
links$source<- sapply(links$source, function(x) which(nodes$name==as.character(x))-1)
links$target<- sapply(links$target, function(x) which(nodes$name==as.character(x))-1)
links<-links[links$value>0,]
sn<-sankeyNetwork(Links = links, Nodes = nodes,
                  Source = "source", Target = "target",
                  Value = "value", NodeID = "name",nodePadding = 10, fontFamily = "Arial Bold",
                  height = 1000, width=1000,
                  fontSize= 12, nodeWidth = 50 ,sinksRight = F, iterations = 0)
sn

library(networkD3)



z<-data.frame(table(metadata$A, metadata$species_gtdb))
z<-dcast(z, Var2 ~ Var1)
z<-z[((z$Environmental>0) & (z$Host>0)),]

#as a table of species 
metadata<-read_excel('/Users/tslaird/Box/leveau_lab/iac_database/metadata_mapping_09132019.xlsx')
metadata<-data.frame( metadata[metadata$A != 'Unspecified',], stringsAsFactors=T)
metadata$A<- factor(metadata$A,levels(factor(metadata$A))[order(table(metadata$A), decreasing = F)])
metadata$B<- factor(metadata$B,levels(factor(metadata$B))[order(table(metadata$B), decreasing = F)])
metadata$C<- factor(metadata$C,levels(factor(metadata$C))[order(table(metadata$C), decreasing = F)])
metadata$class_gtdb<- factor(metadata$class_gtdb,levels(factor(metadata$class_gtdb))[order(table(metadata$class_gtdb), decreasing = F)])

metadata$A_freq <-as.numeric(metadata$A)
metadata$B_freq <-as.numeric(metadata$B)
metadata$C_freq <-as.numeric(metadata$C)
metadata$class_gtdb_freq <-as.numeric(metadata$class_gtdb)


metadata<-metadata[ order( metadata$A_freq,
                           metadata$B_freq,
                           metadata$C_freq,
                           metadata$class_gtdb_freq,
                           decreasing = T), c('A','B','C','class_gtdb')]


library(reshape2)

z<-data.frame(table(metadata$class_gtdb,metadata$A,metadata$B))

metadata_table<-data.frame()
for(cat1 in unique(metadata$A)){
  cat2s<-as.character(unique(metadata[metadata$A == cat1,]$B))
  df1<-metadata[metadata$A == cat1,]
  cat1_df<-data.frame(table(as.character( metadata[metadata$A == cat1,]$A), metadata[metadata$A == cat1,]$class_gtdb))
  cat1_dcast<-dcast(cat1_df, Var1 ~Var2)
  cat1_dcast$category<-"1"
  metadata_table<-rbind.data.frame(metadata_table, cat1_dcast)
  for(cat2 in cat2s){
    cat3s<-as.character(unique(metadata[metadata$B == cat2,]$C))
    df2<-metadata[metadata$B == cat2,]
    cat2_df<-data.frame(table(as.character( df1[df1$B == cat2,]$B), df1[df1$B == cat2,]$class_gtdb))
    cat2_dcast<-dcast(cat2_df, Var1 ~ Var2)
    cat2_dcast$category<-"2"
    metadata_table<-rbind.data.frame(metadata_table, cat2_dcast)
    for(cat3 in cat3s){
      if(cat3 !=0){
      cat3_df<-data.frame(table(as.character( df2[df2$C == cat3,]$C), df2[df2$C == cat3,]$class_gtdb))
      cat3_dcast<-dcast(cat3_df, Var1 ~ Var2)
      cat3_dcast$category<-"3"
      metadata_table<-rbind.data.frame(metadata_table, cat3_dcast)
      }
    }
  }}
    
write.csv(metadata_table, file = '/Users/tslaird/Box/leveau_lab/iac_database/metadata_table_101719.csv',
          row.names = F)

zz<- data.frame(table( metadata[metadata$B=='Plant',]$genus_gtdb))

zz<-data.frame(table(metadata$species_gtdb)

 
