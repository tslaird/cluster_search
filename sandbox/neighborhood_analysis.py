#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 11:46:47 2020

@author: tslaird
"""
import pandas as pd
from collections import Counter
import glob
import concurrent.futures
from tqdm import tqdm
from scipy import stats
import re

#from neighborhoods
cluster_df=pd.read_pickle('/home/tslaird/leveau_lab/cluster_search/FARM_output/iac_positive_all_data_gtdb_20_20.pickle')
cluster_df["nhbrhood_prot_name"]
all_proteins=list([prot for clustr in cluster_df["nhbrhood_prot_name"] for prot in clustr])


master_tally=[]
for i in cluster_df["nhbrhood_prot_name"]:
    i= pd.Series(i)
    master_tally.append((i.value_counts()))

freq = pd.concat(master_tally,axis=1)
freq['combined_freq'] = (freq > 0).sum(axis=1)/len(cluster_df["nhbrhood_prot_name"])



#from index files

def get_annotations(index_file):
    index_prot = pd.read_csv(index_file, sep = "!!" ,header = None, engine='python')
    index_prot.columns = ["filename","assembly","accession","locus_tag","old_locus_tag","name","biosample","protein_name","coordinates","protein_id","pseudogene","protein_seq"]
    annotation_counts= index_prot['protein_name']
    #print("fetched annotations from: "+str(index_file))
    return((annotation_counts.value_counts()))

def get_annotations_fasta(fasta_file):
    prot_re = re.compile('^>[\s|\S]+?(?=\n\s*\n)', re.MULTILINE)
    with open(fasta_file, 'r') as file:
        file_contents = file.read()
        proteins = prot_re.findall(file_contents)
        out_prots = []
        for p in proteins:
            p1 = re.sub('NULL$','NULL!!',p)
            p2 = re.sub('PSEUDOGENE$','PSEUDOGENE!!',p1)
            p3 = p2.replace('\n','!!',1)
            p4 = p3.replace('\n','??')
            out_prots.append(p4)
    df=pd.DataFrame([i.split("!!") for i in out_prots], columns = ["filename","assembly","accession","locus_tag","old_locus_tag","name","biosample","protein_name","coordinates","protein_id","pseudogene","protein_seq"])
    annotation_counts= df['protein_name']
    annotation_counts=annotation_counts.rename(str(fasta_file),axis='columns')
    #print("fetched annotations from: "+str(index_file))
    print(annotation_counts)
    return((annotation_counts.value_counts()))


index_files=glob.glob("*.indexprot")[0:100]
whole_genome_tally=[]
with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
    for i in executor.map(get_annotations, index_files):
        whole_genome_tally.append(i)
        pass
    

def run(f, my_iter):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(f, my_iter), total=len(my_iter)))
    pass
    return results

index_files=glob.glob("/home/tslaird/leveau_lab/cluster_search/index_files/*")
fasta_files= glob.glob("/home/tslaird/leveau_lab/cluster_search/fasta_files/*.fa")
cf_results=run(get_annotations_fasta, fasta_files)
freq_df = pd.concat(cf_results,axis=1)
freq_df['combined_freq'] = (freq_df > 0).sum(axis=1)/len(index_files)
freq_df=freq_df.fillna(0)
freq_df['combined_freq'].sort_values(ascending=False)[1:200]
freq_df.to_csv('/home/tslaird/leveau_lab/cluster_search/test_annotation_matrix.csv')


pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)


#random sample


GTDB_meta = pd.read_csv('/home/tslaird/leveau_lab/cluster_search/gtdb/bac120_metadata_r89.tsv', sep='\t', low_memory = False)
rs_assembly_summary = pd.read_csv('/home/tslaird/leveau_lab/cluster_search/assembly_summary.txt', sep='\t', low_memory = False, skiprows=[0])
GTDB_notation_rs_acc = ["RS_"+ i for i in rs_assembly_summary['# assembly_accession']]
GTDB_meta[GTDB_meta[]]                                                              

GTDB_dict=dict(zip(GTDB_meta.accession, GTDB_meta.gtdb_taxonomy))
rs_assembly_summary['gtdb_notation']=["RS_"+ i for i in rs_assembly_summary['# assembly_accession']]
rs_assembly_summary['gtdb_taxonomy']=rs_assembly_summary['gtdb_notation'].map(GTDB_dict)
rs_assembly_summary[['domain_gtdb','phylum_gtdb','class_gtdb','order_gtdb','family_gtdb','genus_gtdb','species_gtdb']]=rs_assembly_summary['gtdb_taxonomy'].str.split(";", expand = True)
rs_assembly_summary.groupby('phylum_gtdb').size()
                   
def compare_counts(cluster_df,level):
    df=cluster_df.drop_duplicates(subset='assembly')
    iac_pos_counts=df.groupby(level).size()
    RefSeq_counts=rs_assembly_summary.groupby(level).size()
    comparison_df=pd.merge(iac_pos_counts.to_frame(name='iac_positive'),RefSeq_counts.to_frame(name='RefSeq_GTDB'), on= level )
    comparison_df=comparison_df.reset_index()
    comparison_df['pct']=comparison_df['iac_positive']/comparison_df['RefSeq_GTDB']
    comparison_df['pct_of_iac_positive_genomes']=comparison_df['iac_positive']/sum(comparison_df['iac_positive'])
    return(comparison_df)
    
x=compare_counts(cluster_df,'genus_gtdb')



def iac_pvn(cluster_df,taxonomic_name, level):
    pos=cluster_df[cluster_df[level+'_gtdb'] == taxonomic_name]['assembly']
    rs_all=rs_assembly_summary[rs_assembly_summary[level+'_gtdb'] == taxonomic_name ]['# assembly_accession']
    rs_all_name=rs_assembly_summary[rs_assembly_summary[level+'_gtdb'] == taxonomic_name]['organism_name']
    rs_all_level=rs_assembly_summary[rs_assembly_summary[level+'_gtdb'] == taxonomic_name]['assembly_level']
    rs_all_filepath=rs_assembly_summary[rs_assembly_summary[level+'_gtdb'] == taxonomic_name]['ftp_path']
    result=[i in list(pos) for i in list(rs_all)]
    pvn_out = pd.DataFrame({'all_refseq_gtdb':rs_all,'name':rs_all_name,'level':rs_all_level, 'iac_positive':result,'file':rs_all_filepath})
    return(pvn_out)

a_ven_pvn=iac_pvn(cluster_df,'s__Acinetobacter baumannii', 'species') 







cluster_df=pd.read_pickle('/home/tslaird/leveau_lab/cluster_search/FARM_output/iac_positive_all_data_gtdb_20_20.pickle')
cluster_df.sample(frac = 1.0).groupby('genus_gtdb').head(1)

complete_genome=cluster_df[cluster_df['complete_genome']==1]
cluster_df.sample(n=100)



scipy.stats.pearsonr()


#p putida test



