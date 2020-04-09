import pandas as pd
import urllib.request
import concurrent.futures
import glob
import re
import sys
import pandas as pd
import pycurl
import psutil
import hashlib
from io import BytesIO
import xml.etree.ElementTree as ET
import numpy as np
import multiprocessing as mp
import os
import ast


rule all:
    input:
         "results/IacD_analysis/IacD_seqs.fa",
         "results/IacD_analysis/IacD_seqs.fa.trim.msa",
         "results/IacD_analysis/IacD_seqs.fa.trim.msa.fasttree.tree",
         "results/IacD_analysis/IacD_seqs.fa.trim.msa.uncorrected.distmat",
         "results/IacD_analysis/IacD_seqs.fa.trim.msa.jukes_cantor.distmat",
         "results/IacD_analysis/IacD_seqs.fa.trim.msa.kimura_protein.distmat",
         "results/IacD_analysis/RecA_seqs.fa",
         "results/IacD_analysis/RecA_seqs.fa.trim.msa",
         "results/IacD_analysis/RecA_seqs.fa.trim.msa.fasttree.tree",
         "results/IacD_analysis/RecA_seqs.fa.trim.msa.uncorrected.distmat",
         "results/IacD_analysis/RecA_seqs.fa.trim.msa.jukes_cantor.distmat",
         "results/IacD_analysis/RecA_seqs.fa.trim.msa.kimura_protein.distmat",
         "pangenome_analysis/complete_gtdb_genomes-proteins_combined.fasta",
         "pangenome_analysis/complete_gtdb_genomes-proteins_combined.fasta.cdhit.clstr",
         "pangenome_analysis/complete_gtdb_genomes-proteins_combined.fasta.cdhit",
         "pangenome_analysis/complete_gtdb_genomes-proteins_combined.fasta.cdhit.clstr.matrix.csv",
         "pangenome_analysis/complete_gtdb_genomes-proteins_combined.fasta.cdhit.core_proteins.fasta"


rule fetch_IacD:
    input: "results/iac_positive_all_data_gtdb.tsv"
    output: "results/IacD_analysis/IacD_seqs.fa"
    run:
        gtdb=pd.read_csv("results/iac_positive_all_data_gtdb.tsv", sep ='\t')
        protein_seqs=[]
        protein_info=[]
        for i in range(len(gtdb['accession'])):
            #get index of query gene matching iacD
            indices = [i for i, x in enumerate(ast.literal_eval(gtdb['nhbrhood_hit_list'].iloc[i])) if x == "IacD"]
            #index = ast.literal_eval(tree_df['query_list'].iloc[i]).index("IacA")
            for n in range(len(indices)):
                index = indices[n]
                protein_seqs.append(ast.literal_eval( re.sub('nan',"'absent'",gtdb['nhbrhood_prot_seq'].iloc[i]))[index])
                protein_info.append(gtdb['assembly'].iloc[i]+"_"+ast.literal_eval(gtdb['nhbrhood_hit_list'].iloc[i])[index]+"_n"+str(n+1)+"_"+str(gtdb['species_gtdb'].iloc[i])+"_"+ast.literal_eval(re.sub('nan',"'absent'",gtdb['nhbrhood_locus_tags'].iloc[i]))[index]+"_"+ast.literal_eval(re.sub('nan',"'absent'",gtdb['nhbrhood_prot_ids'].iloc[i]))[index])
        protein_seqs=[re.sub('\?\?','\n',i) for i in protein_seqs]
        protein_info=[">"+re.sub(' ','_',i) for i in protein_info]
        protein_info=[re.sub('\.','_',i) for i in protein_info]
        fasta_out=''
        for i,j in zip(protein_info,protein_seqs):
            if re.search('IacD_n1_s__Sphingomonas',i):
                    pass
            else:
                fasta_out+=i+'\n'+j+'\n\n'
        with open("results/IacD_analysis/IacD_seqs.fa",'w+') as output:
            output.writelines(fasta_out)

rule tree_iacD:
    input: "results/IacD_analysis/IacD_seqs.fa"
    output:
        "results/IacD_analysis/IacD_seqs.fa.trim.msa",
        "results/IacD_analysis/IacD_seqs.fa.trim.msa.fasttree.tree"
    shell:
        '''muscle -in {input} -out {input}.msa -maxiters 2 && trimal -in {input}.msa -automated1 > {input}.trim.msa && fasttree {input}.trim.msa > {input}.trim.msa.fasttree.tree'''

rule cluster_IacD:
    input:"results/IacD_analysis/IacD_seqs.fa.trim.msa"
    conda: "envs/emboss.yml"
    output:
        "results/IacD_analysis/IacD_seqs.fa.trim.msa.uncorrected.distmat",
        "results/IacD_analysis/IacD_seqs.fa.trim.msa.jukes_cantor.distmat",
        "results/IacD_analysis/IacD_seqs.fa.trim.msa.kimura_protein.distmat"
    shell:
        """
        distmat -protmethod 0 {input} -outfile {input}.uncorrected.distmat
        distmat -protmethod 1 {input} -outfile {input}.jukes_cantor.distmat
        distmat -protmethod 2 {input} -outfile {input}.kimura_protein.distmat
        """

rule fetch_RecA:
    input: "results/iac_positive_all_data_gtdb.pickle"
    output: "results/IacD_analysis/RecA_seqs.fa"
    run:
        import multiprocessing as mp
        def levenshtein(s, t):
            rows = len(s)+1
            cols = len(t)+1
            dist = [[0 for x in range(cols)] for x in range(rows)]
            # source prefixes can be transformed into empty strings
            # by deletions:
            for i in range(1, rows):
                dist[i][0] = i
            # target prefixes can be created from an empty source string
            # by inserting the characters
            for i in range(1, cols):
                dist[0][i] = i

            for col in range(1, cols):
                for row in range(1, rows):
                    if s[row-1] == t[col-1]:
                        cost = 0
                    else:
                        cost = 1
                    dist[row][col] = min(dist[row-1][col] + 1,      # deletion
                                         dist[row][col-1] + 1,      # insertion
                                         dist[row-1][col-1] + cost) # substitution
            #for r in range(rows):
                #print(dist[r])
            return(dist[row][col])
        def fetch_recA(filename):
            all_seqs=[]
            prot_df = pd.read_csv(filename, sep = "!!" ,header = None, engine='python')
            prot_df.columns = ["filename","assembly","accession","locus_tag","old_locus_tag","name","biosample","protein_name","coordinates","protein_id","pseudogene","protein_seq"]
            # add length
            prot_df['length']=[len(re.sub('\n','',str(i))) for i in prot_df['protein_seq']]
            #filter by name and size
            #prot_df_filtered= prot_df[(prot_df['protein_name'].str.contains('DNA recombination\/repair protein RecA|recombinase A|recombinase RecA')) & (prot_df['length']>=200) & (prot_df['length']<=400) ]
            prot_df_filtered= prot_df[(prot_df['protein_name'].str.contains('DNA_recombination\/repair_protein_RecA|recombinase_A|recombinase_RecA')) & (prot_df['length']>=200) ]
            if len(prot_df_filtered)==1:
                for index, row in prot_df_filtered.iterrows():
                    sequence= str(row['protein_seq'])
                    sequence= re.sub('\?\?','\n',sequence)
                    fasta_entry='>'+'_'.join([str(row['assembly']),str(row['accession']),str(row['locus_tag']),str(row['name']),str(row['biosample']),str(row['protein_name']),str(row['protein_id'])])+'\n'+str(sequence)
                    fasta_entry=re.sub(" |;|:|-|,|\[|\]|\.|\*|\(|\)|\/","_",str(fasta_entry))
                    fasta_entry=re.sub("__|___","_",fasta_entry)
                    all_seqs.append(fasta_entry)
                    print('fetched recA from: ' + filename)
            if len(prot_df_filtered)>1:
                distances=[]
                entries=[]
                Ec_recA='MAIDENKQKALAAALGQIEKQFGKGSIMRLGEDRSMDVETISTGSLSLDIALGAGGLPMGRIVEIYGPESSGKTTLTLQVIAAAQREGKTCAFIDAEHALDPIYARKLGVDIDNLLCSQPDTGEQALEICDALARSGAVDVIVVDSVAALTPKAEIEGEIGDSHMGLAARMMSQAMRKLAGNLKQSNTLLIFINQIRMKIGVMFGNPETTTGGNALKFYASVRLDIRRIGAVKEGENVVGSETRVKVVKNKIAAPFKQAEFQILYGEGINFYGELVDLGVKEKLIEKAGAWYSYKGEKIGQGKANATAWLKDNPETAKEIEKKVRELLLSNPNSTPDFSVDDSEGVAETNEDF'
                for index, row in prot_df_filtered.iterrows():
                    sequence= str(row['protein_seq'])
                    sequence= re.sub('\?\?','\n',sequence)
                    sequence_str= re.sub('\n','',sequence)
                    distances.append(levenshtein(sequence_str, Ec_recA))
                    fasta_entry='>'+'_'.join([str(row['assembly']),str(row['accession']),str(row['locus_tag']),str(row['name']),str(row['biosample']),str(row['protein_name']),str(row['protein_id'])])+'\n'+str(sequence)
                    fasta_entry=re.sub(" |;|:|-|,|\[|\]|\.|\*|\(|\)|\/","_",str(fasta_entry))
                    fasta_entry=re.sub("__|___","_",fasta_entry)
                    entries.append(fasta_entry)
                print(entries[distances.index(min(distances))])
                all_seqs.append(entries[distances.index(min(distances))])
                print('fetched recA from: ' + filename)
            return(all_seqs)

        iac_positive_all_data_gtdb= pd.read_pickle("results/iac_positive_all_data_gtdb.pickle")
        inputs_fetch_recA = [re.sub('.gbff','_proteins.fa.indexprot', i) for i in set(list("index_files/"+iac_positive_all_data_gtdb['filename']))]
        print(inputs_fetch_recA)
        all_RecA_list=[]
        print(all_RecA_list)
        for i in inputs_fetch_recA:
            all_RecA_list.append(fetch_recA(i)[0])
        all_RecA_list=[x for x in all_RecA_list if x]
        all_RecA_str= '\n\n'.join(all_RecA_list) + '\n'
        with open('results/IacD_analysis/RecA_seqs.fa','w+') as file:
            file.writelines(all_RecA_str)

rule tree_RecA:
    input: "results/IacD_analysis/RecA_seqs.fa"
    output:
        "results/IacD_analysis/RecA_seqs.fa.trim.msa",
        "results/IacD_analysis/RecA_seqs.fa.trim.msa.fasttree.tree"
    shell:
        '''muscle -in {input} -out {input}.msa -maxiters 2 && trimal -in {input}.msa -automated1 > {input}.trim.msa && fasttree {input}.trim.msa > {input}.trim.msa.fasttree.tree'''

rule cluster_RecA:
    input:"results/IacD_analysis/RecA_seqs.fa.trim.msa"
    conda: "envs/emboss.yml"
    output:
        "results/IacD_analysis/RecA_seqs.fa.trim.msa.uncorrected.distmat",
        "results/IacD_analysis/RecA_seqs.fa.trim.msa.jukes_cantor.distmat",
        "results/IacD_analysis/RecA_seqs.fa.trim.msa.kimura_protein.distmat"
    shell:
        """
        distmat -protmethod 0 {input} -outfile {input}.uncorrected.distmat
        distmat -protmethod 1 {input} -outfile {input}.jukes_cantor.distmat
        distmat -protmethod 2 {input} -outfile {input}.kimura_protein.distmat
        """


rule make_iac_pan_genome_data:
    input: "results/iac_positive_all_data_gtdb.pickle"
    output: "pangenome_analysis/complete_gtdb_genomes-proteins_combined.fasta"
    run:
        iac_positive_all_data =pd.read_pickle(str(input))
        complete_gtdb_genomes=iac_positive_all_data[iac_positive_all_data['complete_genome']==1]['filename']
        genome_list = ["fasta_files/"+re.sub('.gbff','_proteins.fa', i) for i in list(set(complete_gtdb_genomes))]
        with open("pangenome_analysis/complete_gtdb_genomes-proteins_combined.fasta", 'w') as outfile:
            files = genome_list
            for f in files:
                    with open(f) as infile:
                            outfile.write(infile.read())

rule cluster_proteins_cdhit:
    input: rules.make_iac_pan_genome_data.output
    output:
        clusterfile="pangenome_analysis/complete_gtdb_genomes-proteins_combined.fasta.cdhit.clstr",
        repfile="pangenome_analysis/complete_gtdb_genomes-proteins_combined.fasta.cdhit"
    shell:
        '''cd-hit -i {input} -o {input}.cdhit -d 0 -n 4 -c 0.6'''

rule parse_cdhit:
    input:
        clusterfile=rules.cluster_proteins_cdhit.output.clusterfile,
        repfile=rules.cluster_proteins_cdhit.output.repfile
    output:
        "pangenome_analysis/complete_gtdb_genomes-proteins_combined.fasta.cdhit.clstr.matrix.csv",
        "pangenome_analysis/complete_gtdb_genomes-proteins_combined.fasta.cdhit.core_proteins.fasta"
    params:
        threshold=0.8
    run:
        input_file = str(input.clusterfile)
        with open(input_file) as file:
            cdhit_clstr=file.read()
        all_genomes= set(re.findall("(?<=!!)GCF_.+?(?=!!)",cdhit_clstr))
        all_clusters= cdhit_clstr.split(">Cluster ")[1:]
        rep_names=[]
        rep_full_names=[]
        clstr_tally=[]
        for cluster in all_clusters:
                split_cluster= cluster.split("!!")
                rep_names.append(split_cluster[7]+"_"+split_cluster[9])
                rep_full_name = cluster.split(", ")[1].split("... ")[0]
                rep_full_names.append(rep_full_name)
                genome_tally=[]
                for genome in all_genomes:
                    genome_tally.append(cluster.count("!!"+genome+"!!"))
                clstr_tally.append(genome_tally)
        clustr_df=pd.DataFrame(data=clstr_tally, index=["Cluster_"+str(i) for i in range(0,len(all_clusters))], columns= all_genomes)
        clustr_df['freq'] = (clustr_df > 0).sum(axis=1)/len(all_genomes)
        clustr_df['rep_name'] = rep_names
        clustr_df['rep_full_name']=rep_full_names
        clustr_df.to_csv(input_file+".matrix.csv")
        sum(clustr_df['freq'] >0.8)/ len(clustr_df['freq'])
        core_threshold=params.threshold
        rep_input_file =str(input.repfile)
        with open(rep_input_file) as file:
            rep_proteins=file.read()
        core_rep_proteins=list(clustr_df[clustr_df['freq']>0.8]['rep_full_name'])
        rep_proteins_split=rep_proteins.split('\n\n')[:-1]
        all_prots=[]
        all_seqs=[]
        for i in rep_proteins_split:
            all_prots.append(i.split("\n",1)[0])
            all_seqs.append(i.split("\n",1)[1])
        from itertools import compress
        bool_index=[i in core_rep_proteins for i in all_prots]
        core_fasta="\n\n".join(list(compress(rep_proteins_split, bool_index)))+"\n"
        with open(rep_input_file+".core_proteins"+".fasta","w+") as outfile:
            outfile.write(core_fasta)
