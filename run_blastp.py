#!/usr/bin/env python3
import psutil
import os
import subprocess

number_of_cpus = psutil.cpu_count()

print("Initiating blastp search")
#blast the db using blastp

if not os.path.exists("results/"):
    os.mkdir("results/")

subprocess.call(['blastp','-db','fasta_files_combined/fasta_files_combined_master','-query','iac_proteins.fasta','-out','results/blast_out.txt','-outfmt','6 qseqid qgi qacc sseqid sallseqid sgi sallgi sacc sallacc qstart qend sstart send qseq sseq evalue bitscore score length qlen slen pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus','-max_target_seqs','50000','-num_threads',str(number_of_cpus),'-evalue','0.1'])
