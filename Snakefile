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

combined_fasta_chunks=list(range(0,len(fasta_file_chunks)))

download_dict= dict(zip(sample_names, all_paths))

rule all:
    input:
         expand("gbff_files/{sample}.gbff.gz", sample=sample_names),
         expand("gbff_files_unzipped/{sample}.gbff", sample=sample_names)

rule download_gbff_files:
    output: "gbff_files/{sample}.gbff.gz"
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
