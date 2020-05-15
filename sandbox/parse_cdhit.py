#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 17:03:36 2020

@author: tslaird
"""
import scipy.spatial
import re
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from scipy import stats
import numpy as np
import statsmodels.stats.multitest as multi
import math
import networkx
import time
import concurrent.futures



def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

def parse_cdhit(cluster_file):
    with open(cluster_file) as file:
        cdhit_clstr=file.read()
    all_genomes= set(re.findall("(?<=!!)GCF_.+?(?=!!)",cdhit_clstr))
    all_clusters= cdhit_clstr.split(">Cluster ")[1:]
    rep_names=[]
    rep_full_names=[]
    clstr_tally=[]
    for cluster in all_clusters:
        #print(cluster)
        split_cluster= cluster.split("!!")
        rep_names.append(split_cluster[7]+"_"+split_cluster[9])
        #rep_names.append(split_cluster[6]+"_"+split_cluster[8])
        rep_full_name = cluster.split(", >")[1].split("... ")[0]
        rep_full_names.append(rep_full_name)
        genome_tally=[]
        for genome in all_genomes:
            genome_tally.append(cluster.count("!!"+genome+"!!"))
            #genome_tally.append(cluster.count(">"+genome+"!!"))
        clstr_tally.append(genome_tally)
    cluster_df=pd.DataFrame(data=clstr_tally, index=["Cluster_"+str(i) for i in range(0,len(all_clusters))], columns= all_genomes)
    #make edge list
    edge_list=[]
    for i in cluster_df.index:
        for j in cluster_df.columns:
            if cluster_df.loc[i,j]>0:
                edge_list.append([i,j])
    edge_list_df=pd.DataFrame(edge_list, columns=['cluster','genome'])
    #make adjacency matrix for genes
    
    #cluster_df['freq'] = (cluster_df > 0).sum(axis=1)/len(all_genomes)
    #calculate distance matrix based on presence absence
    dist_list=[]
    for i in cluster_df.columns:
        for j in cluster_df:
            dist_list.append(scipy.spatial.distance.jaccard( cluster_df[i], cluster_df[j]) )
        #scipy.spatial.distance.jaccard(u, v, w=None)[source]
    data_array=np.array(dist_list)
    data_array.shape=(len(cluster_df.columns),len(cluster_df.columns))
    hclust=scipy.cluster.hierarchy.average(np.triu(data_array))
    tree=scipy.cluster.hierarchy.to_tree(hclust)
    newick=getNewick(tree, "", data_array, cluster_df.columns[scipy.cluster.hierarchy.leaves_list(hclust)])
    cluster_df['rep_name'] = rep_names
    cluster_df['rep_full_name']= rep_full_name
    cluster_df = cluster_df[ ['rep_name'] +['rep_full_name']+ [ col for col in cluster_df.columns if col not in ['rep_name','rep_full_name' ] ]]
    return([cluster_df,edge_list_df,newick])
    

cd_hit_pg=parse_cdhit('/home/tslaird/leveau_lab/cluster_search/all_Pput_proteins_cdhit_60.clstr')

cd_hit_pg[0].to_csv('/home/tslaird/leveau_lab/cluster_search/all_Pput_proteins_cdhit_60.matrix')
cd_hit_pg[1].to_csv('/home/tslaird/leveau_lab/cluster_search/all_Pput_proteins_cdhit_60.edges.tsv',sep='\t',index=False, header=False)

with open('/home/tslaird/leveau_lab/cluster_search/all_Pput_proteins_cdhit_60.newick','w') as file:
    file.write(cd_hit_pg[2])


adj=cd_hit_pg[0].iloc[:,3:37].dot(cd_hit_pg[0].iloc[:,3:37].transpose())
adj.to_csv('/home/tslaird/leveau_lab/cluster_search/all_Pput_proteins_cdhit_60.adjacency_matrix.csv')

cd_hit_pg[0].iloc[:,3:34]
mat=cd_hit_pg[0].iloc[:,2:37]
mat[mat>1] =1


df_tuples=[row for row in mat.itertuples()]

start=time.time()
results=[]
pval=[]
phi_values=[]
jaccard=[]
all_contingency_tables=[]
contingency_tables=[]
A=df_tuples[4529]
for row in mat.itertuples():
    combo=tuple(zip(A[1::],row[1::]))
    pospos=combo.count((1,1))
    posneg=combo.count((1,0))
    negpos=combo.count((0,1))
    negneg=combo.count((0,0))
    contingency_table=tuple([[pospos,posneg],[negpos,negneg]]) 
    #jaccard.append(scipy.spatial.distance.jaccard(A[1::],row[1::]))
    try:
        phi= ((pospos*negneg)-(posneg*negpos))/math.sqrt((pospos+posneg)*(negpos+negneg)*(pospos+negpos)*(posneg+negneg))
        
    except:
        phi=0
    phi_values.append(phi)
    if contingency_table in all_contingency_tables:
        output=pval[all_contingency_tables.index(contingency_table)]
    else:
        all_contingency_tables.append(contingency_table)   
        fe_test=scipy.stats.fisher_exact(contingency_table, alternative='greater')
        pval.append(fe_test[1])
        output=fe_test[1]
    results.append(output)
    
print(time.time()-start)   

p_adj=multi.multipletests(results,method='fdr_bh')
pval_df1=pd.DataFrame([mat.index,cd_hit_pg[0].iloc[:,0],results,p_adj[1]]).transpose()


#map version

df_tuples=[row for row in mat.itertuples()]
results=[]
pval=[]
phi_values=[]
jaccard=[]
all_contingency_tables=[]
contingency_tables=[]
A=df_tuples[4529]

def fetest(A,B):
    combo=tuple(zip(A[1::],B[1::]))
    pospos=combo.count((1,1))
    posneg=combo.count((1,0))
    negpos=combo.count((0,1))
    negneg=combo.count((0,0))
    contingency_table=tuple([[pospos,posneg],[negpos,negneg]]) 
    #jaccard.append(scipy.spatial.distance.jaccard(A[1::],row[1::]))
    try:
        phi= ((pospos*negneg)-(posneg*negpos))/math.sqrt((pospos+posneg)*(negpos+negneg)*(pospos+negpos)*(posneg+negneg))
        
    except:
        phi=0
    phi_values.append(phi)
    if contingency_table in all_contingency_tables:
        output=pval[all_contingency_tables.index(contingency_table)]
    else:
        all_contingency_tables.append(contingency_table)   
        fe_test=scipy.stats.fisher_exact(contingency_table, alternative='greater')
        pval.append(fe_test[1])
        output=fe_test[1]
    results.append(output)
    return(output)


C=df_tuples*3
D=list(itertools.chain.from_iterable(itertools.repeat(x, 3) for x in df_tuples))
x=pd.DataFrame([C,D])
start=time.time()
mapout=map(fetest, C,D) 
out=list(mapout)
print(time.time()-start)   









#concurrent futures version #### one to use
df_tuples=list(mat.itertuples(name=None))

results=[]
pval=[]
phi_values=[]
jaccard=[]
all_contingency_tables={}
contingency_tables=[]

def iter_fetest(row):
    results=[]
    for row2 in mat.itertuples(name=None):
        combo=tuple(zip(row[1::],row2[1::]))
        pospos=combo.count((1,1))
        posneg=combo.count((1,0))
        negpos=combo.count((0,1))
        negneg=combo.count((0,0))
        contingency_table=tuple([tuple([pospos,posneg]),tuple([negpos,negneg])]) 
        #jaccard.append(scipy.spatial.distance.jaccard(A[1::],row[1::]))
        #try:
        #    phi= ((pospos*negneg)-(posneg*negpos))/math.sqrt((pospos+posneg)*(negpos+negneg)*(pospos+negpos)*(posneg+negneg))     
        #except:
        #    phi=0
        #phi_values.append(phi)
        if contingency_table in all_contingency_tables:
            output=all_contingency_tables[contingency_table]
        else:
            fe_test=scipy.stats.fisher_exact(contingency_table, alternative='greater')
            all_contingency_tables[contingency_table]=fe_test[1]  
            output=fe_test[1]
        results.append(output)
    return(results)

start=time.time()
out=[]
with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
    for i in executor.map(iter_fetest, df_tuples[0:500]):
        out.append(i)
        pass    
print(time.time()-start)   

Z=np.array(out)
Z_indices=np.triu_indices(len(Z),-1)
Z[Z_indices]=0
#Z_tri=np.tril(Z)
#np.fill_diagonal(Z_tri, 0)
#Z_masked=np.ma.masked_equal(Z_tri,0)
#Z_data=ma.getdata(Z_masked)
Z_data=Z[Z_masked>0]
p_adj=multi.multipletests(Z_data,method='fdr_bh')[1]


Z_adj = np.zeros(Z.shape)
indices = np.tril_indices(len(Z),-1)
Z_adj[indices] = p_adj






###############################




#vectorized version


def co_occur(array_1,array_2):
    combo=tuple(zip(array_1,array_2))
    pospos=combo.count((1,1))
    posneg=combo.count((1,0))
    negpos=combo.count((0,1))
    negneg=combo.count((0,0))
    contingency_table=tuple([[pospos,posneg],[negpos,negneg]]) 
    fe_test=scipy.stats.fisher_exact(contingency_table, alternative='greater')
    return(fe_test[1])
                

mat_t=mat.transpose()
all_data=pd.Series([i for i in mat_t])
A=mat_t['Cluster_1'].values
co_occur(A,np.array([[A,A]])) 



def haversine(a, b):
    c=a+b
    B=c-a
    A=a
    return tuple(zip(A,B))

z=haversine([1,0,1], np.array([[0,1,0],[0,1,1]] )  )




def myfunc(a, b):
    c=a+b
    
    return c


vfunc = np.vectorize(myfunc)
vfunc([[1, 2], [3, 4]], [2,3])


numbers1 = np.array([[0,1,0],[0,1,1]])
numbers2 = [1,0,1]
  
result = map(lambda x, y: x - y, numbers1, numbers2) 
print(list(result)) 



scipy.stats.pearsonr(mat.loc['Cluster_17525'],mat.loc[i])


x=adj.iloc[4000:7000,4000:7000].stack().reset_index()
x.to_csv('/home/tslaird/leveau_lab/cluster_search/all_Pput_proteins_cdhit_60.adjacency_edges.tsv',sep='\t', index=False, header=False)




import networkx as nx
import matplotlib
from networkx.algorithms.community import greedy_modularity_communities
G = nx.karate_club_graph()
G = nx.fast_gnp_random_graph(10000,0.4)
c = list(greedy_modularity_communities(G))
sorted(c[0])
nx.draw(G)




   
#cluster_df.to_csv(input_file+".matrix.csv")
sum(cluster_
    df['freq'] >0.8)/ len(cluster_df['freq'])
core_threshold=params.threshold
rep_input_file =str(input.repfile)
#with open(rep_input_file) as file:
#   rep_proteins=file.read()
core_rep_proteins=list(cluster_df[cluster_df['freq']>0.8]['rep_full_name'])
rep_proteins_split=rep_proteins.split('\n\n')[:-1]
all_prots=[]
all_seqs=[]
for i in rep_proteins_split:
    all_prots.append(i.split("\n",1)[0])
    all_seqs.append(i.split("\n",1)[1])
from itertools import compress
bool_index=[i in core_rep_proteins for i in all_prots]
core_fasta="\n\n".join(list(compress(rep_proteins_split, bool_index)))+"\n"
#with open(rep_input_file+".core_proteins"+".fasta","w+") as outfile:
#    outfile.write(core_fasta)
    
  start=time.time()

results=[]
pval=[]
phi_values=[]
jaccard=[]
all_contingency_tables=[]
contingency_tables=[]
A=mat.loc['Cluster_4217']
for i in mat.index[1:1000]:
    contingency_table=np.array(pd.crosstab(A,mat.loc[i],dropna = False))[::-1,::-1]
    #combo=tuple(zip(A,mat.loc[i]))
    #pospos=combo.count((1,1))
    #posneg=combo.count((1,0))
    #negpos=combo.count((0,1))
    #negneg=combo.count((0,0))
    #contingency_table=tuple([[pospos,posneg],[negpos,negneg]]) 
    #jaccard.append(scipy.spatial.distance.jaccard(mat.loc['Cluster_17525'],mat.loc[i]))
    #try:
    #    phi= ((pospos*negneg)-(posneg*negpos))/math.sqrt((pospos+posneg)*(negpos+negneg)*(pospos+negpos)*(posneg+negneg))
        
    #except:
    #    phi=0
    #phi_values.append(phi)
    contingency_table_ph=tuple(contingency_table.tolist())
    if contingency_table_ph in all_contingency_tables:
        output=pval[all_contingency_tables.index(contingency_table_ph)]
    else:
        all_contingency_tables.append(contingency_table_ph)   
        fe_test=scipy.stats.fisher_exact(contingency_table, alternative='greater')
        pval.append(fe_test[1])
        output=fe_test[1]
    results.append(output)
    
print(time.time()-start)     