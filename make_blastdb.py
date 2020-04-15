#!/usr/bin/env python3

import glob
import os
import concurrent.futures
import psutil
import subprocess

number_of_cpus = psutil.cpu_count()
chunk_size = 5000
#locates all files in the fasta_files directory ending in .fa
files = glob.glob('fasta_files/*.fa')
file_chunks = [files[i * chunk_size:(i + 1) * chunk_size] for i in range((len(files) + chunk_size - 1) // chunk_size )]

def combine_files(file_list_index):
    with open('fasta_files_combined/fasta_files_combined_'+str(file_list_index), 'w') as outfile:
        files = file_chunks[file_list_index]
        print("Writing "+'fasta_files_combined/fasta_files_combined_'+str(file_list_index))
        for f in files:
                with open(f) as infile:
                        outfile.write(infile.read())

if not os.path.exists("fasta_files_combined"):
    os.mkdir("fasta_files_combined")

print("Making chunks of "+chunk_size+" fasta files")
with concurrent.futures.ProcessPoolExecutor(max_workers=number_of_cpus) as executor:
    for _ in executor.map(combine_files, range(len(file_chunks))):
        pass

###########
#locates all files in the all_proteins_combined directory
combined_files_list = glob.glob('fasta_files_combined/fasta_files_combined_*')
def make_blastdb(file_list_index):
    combined_file = 'fasta_files_combined/fasta_files_combined_'+str(file_list_index)
    print("Making blastdb for: "+combined_file)
    subprocess.call(['makeblastdb','-in',combined_file,'-dbtype','prot'])

with concurrent.futures.ProcessPoolExecutor(max_workers=number_of_cpus*2) as executor:
    for _ in executor.map(make_blastdb, range(len(combined_files_list))):
        pass
print("Finished making blast databases")
###########
#make alias db
print("Making alias database")
combined_files_str = " ".join(combined_files_list)
subprocess.call(['blastdb_aliastool','-dblist',combined_files_str,'-dbtype','prot','-out','fasta_files/all_proteins_combined_master','-title','all_proteins_combined_master'])
