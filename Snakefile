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

if os.path.exists("assembly_summary.txt"):
    print("assembly summary already downloaded")
else:
    url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt'
    urllib.request.urlretrieve(url,"assembly_summary.txt")

if os.path.exists('assembly_summary_filtered.txt'):
    assembly_summary_filtered= pd.read_csv("assembly_summary_filtered.txt", sep='\t')
else:
    assembly_summary = pd.read_csv('assembly_summary.txt', sep = '\t', low_memory = False, skiprows=[0])
    assembly_summary_filtered = assembly_summary[((assembly_summary['assembly_level'] == 'Complete Genome') & (assembly_summary['version_status']=='latest') & assembly_summary['organism_name'].str.contains('Pseudomonas putida'))]
    assembly_summary_filtered.to_csv("assembly_summary_filtered.txt", sep='\t')
#assembly_summary_filtered = assembly_summary[(assembly_summary['version_status']=='latest')]
all_paths = assembly_summary_filtered['ftp_path']+ "/"+[i.split('/')[-1] for i in assembly_summary_filtered['ftp_path']] + \
    "_genomic.gbff.gz"
all_paths.to_csv("ftp_paths.txt", index=False)
sample_names_df=pd.DataFrame([re.search("GCF_.+genomic(?=.gbff.gz)?",i.split('/')[-1])[0] for i in all_paths])
sample_names_df.to_csv("sample_names.txt", sep='\t', index=False, header=False)
sample_names=[re.search("GCF_.+genomic(?=.gbff.gz)?",i.split('/')[-1])[0] for i in all_paths]



chunk_size=5000
fasta_names= ["fasta_files/"+i+"_proteins.fa" for i in sample_names]
fasta_file_chunks= [fasta_names[i * chunk_size:(i + 1) * chunk_size] for i in range((len(fasta_names) + chunk_size - 1) // chunk_size )]
combined_fasta_chunks_index=list(range(0,len(fasta_file_chunks)))




download_dict= dict(zip(sample_names, all_paths))

rule all:
    input:
         expand("gbff_files/{sample}.gbff.gz", sample=sample_names),
         expand("gbff_files/{sample}.gbff.gz_md5checksums.txt", sample= sample_names),
         expand("gbff_files_unzipped/{sample}.gbff", sample=sample_names),
         expand("fasta_files/{sample}_proteins.fa", sample=sample_names),
         expand("fasta_files_combined/all_proteins_combined_{db_id}.fa", db_id= combined_fasta_chunks_index),
         expand("fasta_files_combined/all_proteins_combined_{db_id}.fa.phr", db_id=combined_fasta_chunks_index),
         expand("fasta_files_combined/all_proteins_combined_{db_id}.fa.pin", db_id=combined_fasta_chunks_index),
         expand("fasta_files_combined/all_proteins_combined_{db_id}.fa.psq", db_id=combined_fasta_chunks_index),
         # "fasta_files_combined/all_proteins_combined_master.pal",
         # "results/blast_out",
         # "results/blast_output_table.txt"

rule download_gbff_files:
    output: "gbff_files/{sample}.gbff.gz" , "gbff_files/{sample}.gbff.gz_md5checksums.txt"
    params:
        # dynamically generate the download link directly from the dictionary
        download_link = lambda wildcards: download_dict[wildcards.sample]
    run:
        path = params["download_link"]
        name = 'gbff_files/'+path.split('/')[-1]
        md5path= path.rsplit('/',1)[0]+"/md5checksums.txt"
        md5name= 'gbff_files/'+path.split('/')[-1]+"_md5checksums.txt"
        print(md5path)
        print(md5name)
        retries=20
        while(retries > 0):
            try:
                urllib.request.urlretrieve(path,name)
                urllib.request.urlretrieve(md5path,md5name)
                with open(name, 'rb') as file_to_check:
                    data = file_to_check.read()
                    md5_returned = str(hashlib.md5(data).hexdigest())
                with open(md5name,'r') as md5_file:
                    md5_file_text=md5_file.read()
                    md5_sum=str(re.findall('\w+(?=\s\s.+gbff.gz)',md5_file_text)[0])
                if md5_returned==md5_sum:
                    print("Fetched " + name)
                    break
            except:
                print("Retrying download from " + path)
                retries = retries - 1
                continue

rule unzip_gbff_files:
    input:
        "gbff_files/{sample}.gbff.gz"
    output:
        "gbff_files_unzipped/{sample}.gbff"
    run:
        import psutil
        import gzip
        import os
        file = str(input)
        output_name = 'gbff_files_unzipped/'+file.split('/')[1].split('.gz')[0]
        if os.path.exists(output_name):
            print(file)
            print(file + " already unzipped")
        else:
            f = gzip.open(file, 'rb')
            file_content = f.read()
            f.close()
            print('Unzipping '+ output_name)
            output = open(output_name, 'wb')
            output.write(file_content)
            output.close()

rule parse_gbff_files:
    input:
        "gbff_files_unzipped/{sample}.gbff"
    output:
        "fasta_files/{sample}_proteins.fa"
    run:
        import re
        import textwrap
        import multiprocessing as mp
        import concurrent.futures
        import glob
        import os
        locus_re= re.compile('(?<=LOCUS\s{7})\w+')
        definition_re = re.compile('(?<=DEFINITION\s{2})[\s\S]+?(?=\nACC)')
        definition_sub = re.compile('\n|(\s{12})')
        locus_tag_re = re.compile('(?<=gene\W{12})[\S\s]+?locus_tag[\S\s]+?(?="\n)')
        locus_tag_sub = re.compile('[\s\S]+"')
        old_locus_tag_re = re.compile('(?<=gene\W{12})[\S\s]+?old_locus_tag[\S\s]+?(?="\n)')
        product_re = re.compile('(?<=product=")[\S\s]+?(?=")')
        coords_re = re.compile('(?<=(tRNA|rRNA|gene|ncRN)(\s|A){12})\S+')
        protein_id_re = re.compile('(?<=protein_id=")[\s\S]+?(?=")')
        protein_re = re.compile('(?<=translation=")[\S\s]+?(?=")')
        biosample_re = re.compile('(?<=BioSample:\s)\w+')
        assembly_re = re.compile('(?<=Assembly:\s)\S+')
        features_re = re.compile('(?<=\n)\W{5}(?=gene\W{12})')
        product_sub = re.compile('\n|\s{19}')
        protein_out_sub = re.compile('\n|(\s+)')
        separate_re = re.compile('//\n')
        gbff_file= str(input)
        with open(gbff_file) as file:
            file_text=file.read()
        loci = separate_re.split(file_text)
        all_proteins = []
        for i in loci:
            if locus_re.search(i):
                locus = locus_re.search(i).group(0)
                definition = definition_re.search(i).group(0)
                definition = definition_sub.sub('',definition)
                biosample =  biosample_re.search(i)
                if biosample:
                    biosample_out = biosample.group(0)
                else:
                    biosample_out = 'NULL'
                if assembly_re.search(i):
                    assembly = assembly_re.search(i).group(0)
                else:
                    assembly = re.search('GCF_\d+?(?=\.)',gbff_file).group(0)
                features = features_re.split(i)
                for f in features[1:]:
                     locus_tag = locus_tag_re.search(f).group(0)
                     locus_tag = locus_tag_sub.sub('', locus_tag)
                     if old_locus_tag_re.search(f):
                         old_locus_tag = old_locus_tag_re.search(f).group(0)
                         old_locus_tag = locus_tag_sub.sub('', old_locus_tag)
                     else:
                         old_locus_tag = 'NULL'
                     if product_re.search(f):
                         product = product_re.search(f).group(0)
                         product = product_sub.sub('',product)
                     else:
                         product = 'NULL'
                     coords = coords_re.search(f).group(0)
                     protein_id = protein_id_re.search(f)
                     if protein_id:
                         protein_id_out = protein_id.group(0)
                     else:
                         protein_id_out = 'NULL'
                     protein = protein_re.search(f)
                     if protein:
                         protein_out = protein.group(0)
                         protein_out = protein_out_sub.sub('', protein_out)
                         protein_out = textwrap.fill(protein_out, 70)
                         almost_all_items = ">%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s" % (gbff_file.split("/")[-1],assembly,locus,locus_tag,old_locus_tag,definition,biosample_out,product,coords, protein_id_out)
                         almost_all_items = re.sub('\n|\s{4,}','',almost_all_items)
                         all_items = "%s\n%s\n" % (almost_all_items, protein_out)
                     else:
                         almost_all_items = ">%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s" % (gbff_file.split("/")[-1],assembly,locus,locus_tag,old_locus_tag,definition,biosample_out,product,coords, protein_id_out)
                         almost_all_items = re.sub('\n|\s{4,}','',almost_all_items)
                         all_items = "%s\n" % (almost_all_items)
                     #all_proteins.add(">"+assembly+'•'+locus+'•'+locus_tag+'•'+biosample+'•'+product+'•'+coords+'•'+protein_id+'\n'+protein+'\n')
                     all_proteins.append(all_items)
        result= '\n'.join(all_proteins)+'\n'
        outname = [re.sub('.gbk|.gbff','_proteins.fa', gbff_file)]
        outname = ''.join(outname).split('/')[-1]
        with open("fasta_files/"+outname,'w+') as output:
            output.writelines(result)
        print("extracted proteins for: " + gbff_file)



combined_fasta_dict= dict(zip(combined_fasta_chunks_index,fasta_file_chunks))

def get_combined_fasta_samples(wildcards):
    return [i for i in combined_fasta_dict[int(wildcards.db_id)]]

rule combine_fasta:
    input: get_combined_fasta_samples
    # input: expand("fasta_files/{sample}_proteins.fa", sample=sample_names)
    output: "fasta_files_combined/all_proteins_combined_{db_id}.fa"
    params:
        file_list = lambda wildcards: combined_fasta_dict[int(wildcards.db_id)]
    run:
        import glob
        import concurrent.futures
        import psutil
        import subprocess
        number_of_cpus = psutil.cpu_count()
        # file_chunks= [fasta_files[i * chunk_size:(i + 1) * chunk_size] for i in range((len(fasta_files) + chunk_size - 1) // chunk_size )]
        # fasta_files = glob.glob('fasta_files/*.fa')
        # file_chunks = [fasta_files[i * chunk_size:(i + 1) * chunk_size] for i in range((len(fasta_files) + chunk_size - 1) // chunk_size )]
        # os.mkdir("fasta_files_combined")
        with open('fasta_files_combined/all_proteins_combined_'+wildcards.db_id+".fa", 'w') as outfile:
            files = params.file_list
            for f in files:
                    with open(f) as infile:
                            outfile.write(infile.read())

rule make_blastdb:
    input:
        file= rules.combine_fasta.output
    output:
        "fasta_files_combined/all_proteins_combined_{db_id}.fa.phr",
        "fasta_files_combined/all_proteins_combined_{db_id}.fa.pin",
        "fasta_files_combined/all_proteins_combined_{db_id}.fa.psq",
    shell:''' makeblastdb -in {input.file} -dbtype prot'''
