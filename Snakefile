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



chunk_size=5
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
         "fasta_files_combined/all_proteins_combined_master.pal",
         "results/blast_out",
         "results/blast_output_table.txt",
         "results/iac_positive_accessions.txt",
         "results/iac_positive_df.tsv",
         "results/iac_positive_df.pickle",
         directory("index_files"),
         "results/iac_positive_all_data.xlsx",
         "results/iac_positive_all_data_gtdb.xlsx",
         "results/iac_positive_all_data_gtdb.tsv"

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
        pseudogene_re = re.compile('(?<=\/)pseudo(?=\n)')
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
                     if pseudogene_re.search(f):
                         pseudogene = "PSEUDOGENE"
                     else:
                         pseudogene = "NULL"
                     if protein_id:
                         protein_id_out = protein_id.group(0)
                     else:
                         protein_id_out = 'NULL'
                     protein = protein_re.search(f)
                     if protein:
                         protein_out = protein.group(0)
                         protein_out = protein_out_sub.sub('', protein_out)
                         protein_out = textwrap.fill(protein_out, 70)
                         almost_all_items = ">%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s" % (gbff_file.split("/")[-1],assembly,locus,locus_tag,old_locus_tag,definition,biosample_out,product,coords, protein_id_out, pseudogene)
                         almost_all_items = re.sub('\n|\s{4,}','',almost_all_items)
                         all_items = "%s\n%s\n" % (almost_all_items, protein_out)
                     else:
                         almost_all_items = ">%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s" % (gbff_file.split("/")[-1],assembly,locus,locus_tag,old_locus_tag,definition,biosample_out,product,coords, protein_id_out, pseudogene)
                         almost_all_items = re.sub('\n|\s{4,}','',almost_all_items)
                         all_items = "%s\n" % (almost_all_items)
                     #all_proteins.add(">"+assembly+'•'+locus+'•'+locus_tag+'•'+biosample+'•'+product+'•'+coords+'•'+protein_id+'\n'+protein+'\n')
                     all_items= re.sub(' |,','_',all_items)
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

rule make_master_db:
    input:
        psqfiles=expand("fasta_files_combined/all_proteins_combined_{db_id}.fa.psq", db_id = combined_fasta_chunks_index ),
        fastafiles=expand("fasta_files_combined/all_proteins_combined_{db_id}.fa", db_id = combined_fasta_chunks_index )
    output: "fasta_files_combined/all_proteins_combined_master.pal"
    # params:
    #     combined_files_str= " ".join(glob.glob('fasta_files_combined/all_proteins_combined_*.fa'))
    # wildcard_constraints: db_id= '\d+'
    run:
        combined_files_list = input.fastafiles
        combined_files_str = " ".join(combined_files_list)
        subprocess.call(['blastdb_aliastool','-dblist',combined_files_str,'-dbtype','prot','-out','fasta_files_combined/all_proteins_combined_master','-title','all_proteins_combined_master'])

rule blastp:
    input: rules.make_master_db.output
    output:"results/blast_out"
    params:
        blastp_threads = psutil.cpu_count()
    # run:
    #     subprocess.call(['blastp','-db','fasta_files_combined/all_proteins_combined_master','-query','iac_proteins.fasta','-out','blast_out','-outfmt','6 qseqid qgi qacc sseqid sallseqid sgi sallgi sacc sallacc qstart qend sstart send qseq sseq evalue bitscore score length qlen slen pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus','-max_target_seqs','50000','-num_threads','4','-evalue','0.1'])
    shell:
        '''blastp -db fasta_files_combined/all_proteins_combined_master -query iac_proteins.fasta -out results/blast_out -outfmt "6 qseqid qgi qacc sseqid sallseqid sgi sallgi sacc sallacc qstart qend sstart send qseq sseq evalue bitscore score length qlen slen pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus" -max_target_seqs 50000 -num_threads {params.blastp_threads} -evalue 0.1'''

rule add_header_to_blast_out:
    input: rules.blastp.output
    output: "results/blast_output_table.txt"
    shell:'''echo -e "qseqid\tqgi\tqacc\tsseqid\tsallseqid\tsgi\tsallgi\tsacc\tsallacc\tqstart\tqend\tsstart\tsend\tqseq\tsseq\tevalue\tbitscore\tscore\tlength\tqlen\tslen\tpident\tnident\tmismatch\tpositive\tgapopen\tgaps\tppos\tframes\tqframe\tsframe\tbtop\tstaxids\tsscinames\tscomnames\tsblastnames\tsskingdoms\tstitle\tsalltitles\tsstrand\tqcovs\tqcovhsp\tqcovus" | cat - results/blast_out > results/blast_output_table.txt'''

rule download_gtdb:
    input: "results/blast_output_table.txt"
    output: 'gtdb/bac120_metadata_r89.tsv'
    run:
        print("downloading gtdb data")
        urllib.request.urlretrieve("https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/bac120_metadata_r89.tsv",'gtdb/bac120_metadata_r89.tsv')


rule parse_blastp:
    input: 'gtdb/bac120_metadata_r89.tsv'
    threads: 1
    output:
        "results/iac_positive_accessions.txt",
        "results/iac_positive_df.tsv",
        "results/iac_positive_df.pickle"
    run:
        import re
        import sys
        import pandas as pd
        import pycurl
        from io import BytesIO
        import time
        import xml.etree.ElementTree as ET
        import numpy as np
        import multiprocessing as mp
        import glob
        import os
        blastout = pd.read_csv("results/blast_output_table.txt", sep='\t', low_memory = False)
        print('Editing blastp output file')
        blastout[["filename","assembly","accession","locus_tag","old_locus_tag","name","biosample","protein_name","coordinates","protein_id","pseudogene"]] = blastout['stitle'].str.split("!!", expand = True)
        blastout['direction']= [-1 if re.match('complement',str(c)) else 1 for c in blastout['coordinates']]
        blastout['join_status'] = [1 if re.match('join',str(c)) else 0 for c in blastout['coordinates']]
        #extracts start and end coordinates from the coordinates column
        blastout['start_coord'] = [re.search('\d+?(?=\.\.(\d|\>))',str(c)).group(0) for c in blastout['coordinates'] ]
        blastout['start_coord'] = [re.sub('join|complement|join\(|complement|>|<|\)|\(|,\d+',"",str(c)) for c in blastout['start_coord'] ]
        blastout['start_coord'] = blastout['start_coord'].astype(int)
        blastout['end_coord'] = [re.search('(?<=\.(\.|\>))\d+',str(c)).group(0) for c in blastout['coordinates'] ]
        blastout['end_coord'] = [re.sub('>|<|\)|\(',"",c) for c in blastout['end_coord'] ]
        blastout['end_coord'] = blastout['end_coord'].astype(int)
        blastout['middle_coord'] = (  blastout['start_coord'] + blastout['end_coord'] )/2
        print('Filtering blastp output file')
        blastout_filtered = blastout[(blastout["qseqid"].isin(["IacA","IacB","IacC","IacD","IacE","IacH","IacI"])) &
                (blastout['evalue'] <= 1e-10) & (blastout['bitscore']>=50)]
        #print raw blast out statistics:
        iac_proteins = ['IacA','IacB','IacC','IacD','IacE','IacF','IacG','IacR','IacH','IacI']
        raw_blast_stats=[]
        for i in iac_proteins:
            genome_count = len(list(set(blastout[blastout['qseqid']== i]['assembly'] )))
            total_count = len(blastout[blastout['qseqid']== i])
            raw_blast_stats.append('Genomes with homologue to '+i+" "+str(genome_count)+"\n")
            raw_blast_stats.append('Total homologue count for '+i+" "+str(total_count)+"\n")
        with open("results/raw_blast_stats.txt", 'w') as f:
            f.writelines(raw_blast_stats)
        filtered_blast_stats=[]
        for i in iac_proteins:
            genome_count = len(list(set(blastout_filtered[blastout_filtered['qseqid']== i]['assembly'] )))
            total_count = len(blastout_filtered[blastout_filtered['qseqid']== i])
            filtered_blast_stats.append('Genomes with homologue to '+i+" "+str(genome_count)+"\n")
            filtered_blast_stats.append('Total homologue count for '+i+" "+str(total_count)+"\n")
        with open("results/filtered_blast_stats.txt", 'w') as f:
            f.writelines(filtered_blast_stats)
        print('Defining function')
        def parseblastout2(accession):
            print("parsing proteins from: "+ accession )
            genome_match= blastout_filtered[blastout_filtered['accession'] ==accession]
            #if len(genome_match) > 4:
            #the following filters the genomes if they have 6 of 7 of the iacABCDEH or I
            #this filter takes into account having multiples of a hit and filters based on presence or absence
            if sum(gene in list(genome_match['qseqid']) for gene in ["IacA","IacB","IacC","IacD","IacE","IacH","IacI"]) >= 6:
                genome_match_sorted = genome_match.sort_values(by='start_coord')
                middle_coords = genome_match_sorted['middle_coord']
                middle_coords_np = np.array(middle_coords)
                genome_match_sorted['groups']=list(np.cumsum([0] + list(1*(middle_coords_np[1:] - middle_coords_np[0:-1] >= 10000))) + 1)
                list_of_clusters = []
                for cluster_number, df in genome_match_sorted.groupby('groups'):
                    if len(df) >= 1:
                        hit_list = list(df['locus_tag'])
                        old_locus_hit_list = list(df['old_locus_tag'])
                        protein_name_list = list(df['protein_name'])
                        protein_id_list = list(df['protein_id'])
                        query_list = list(df['qseqid'])
                        coord_list = list(zip(df['start_coord'], df['end_coord'],df['direction'],df['qseqid']))
                        if sum(df['direction']) < 0:
                            df['actual_start_tmp'] = df['start_coord']
                            df['start_coord']= df['end_coord'] * -1
                            df['end_coord']= df['actual_start_tmp']* -1
                            df['direction'] = df['direction'] *-1
                            df.sort_values(by='start_coord',inplace=True)
                        order = [ ("| " + gene['qseqid'] + " 〉") if gene['direction'] == 1 else ("〈 " + gene['qseqid'] + " |") for index,gene in df.iterrows()  ]
                        order_pident = [ ("| " + gene['qseqid']+':'+ str(round(gene['pident'],2)) + " 〉") if gene['direction'] == 1 else ("〈 " +  gene['qseqid']+':'+ str(round(gene['pident'],2)) + " |") for index,gene in df.iterrows()  ]
                        dist = list(np.array(df['start_coord'][1:]) -  np.array(df['end_coord'][:-1]))
                        dist = ["-"+str(d)+"-" for d in dist]
                        #obtains "| A 〉-23-| B 〉-23-| C 〉"
                        synteny_dir_dist=''.join(sum(zip(order, dist+[0]), ())[:-1])
                        synteny_dir_dist = re.sub("Iac" ,"", synteny_dir_dist)
                        #obtains "| A 〉| B 〉| C 〉"
                        synteny_dir =''.join(order)
                        synteny_dir = re.sub("Iac" ,"", synteny_dir)
                        #obtains "| A:23.23 〉| B:23.23〉| C:23.23 〉"
                        synteny_dir_pident =''.join(order_pident)
                        synteny_dir_pident = re.sub("Iac" ,"", synteny_dir_pident)
                        #obtains "A-B-C"
                        synteny= re.sub("\n" ,"-", df['qseqid'].to_string(index=False))
                        synteny= re.sub("Iac" ,"", synteny)
                        filename = list(set(df['filename']))[0]
                        biosample = list(set(df['biosample']))[0]
                        name = list(set(df['name']))[0]
                        assembly= re.sub("\{|\}|\'","", str(set(df['assembly'])) )
                        accession = re.sub("\{|\}|\'","", str(set(df['accession'])) )
                        #sgi = re.sub("\{|\}|\'","", str(set(df['sgi'])) )
                        cluster_len= max(df['end_coord']) - min(df['start_coord'])
                        number_of_hits = len(df)
                        pseudogene_list= list(df['pseudogene'])
                        list_of_clusters.append([accession, filename, biosample, number_of_hits , cluster_len, synteny, synteny_dir_dist, synteny_dir_pident, synteny_dir, assembly, accession, name, hit_list, old_locus_hit_list, protein_name_list, protein_id_list, pseudogene_list, query_list, coord_list, cluster_number])
                return(list_of_clusters)
        print('\nParsing filtered blastp output file')
        parse_blastp_input= list(set(blastout_filtered['accession']))
        print(parse_blastp_input)
        result_list=[]
        with concurrent.futures.ThreadPoolExecutor(max_workers=None) as executor:
            for i in executor.map(parseblastout2, parse_blastp_input):
                result_list.append(i)
                pass
        # pool_parse_blastp = mp.Pool(mp.cpu_count()-2)
        # pool_parse_blastp_output = pool_parse_blastp.map(parseblastout2, [i for i in parse_blastp_input])
        # pool_parse_blastp.close()
        print('\nObtaining iac cluster positive data frame')
        iac_positive = list(filter(None, result_list))
        iac_positive_flat = [item for sublist in iac_positive for item in sublist]
        iac_positive_df= pd.DataFrame(iac_positive_flat, columns=('genome_acc','filename','biosample', 'hits', 'cluster_length','synteny','synteny_dir_dist','synteny_dir_pident','synteny_dir','assembly','accession','name','hit_list','old_locus_hit_list','protein_name_list','protein_id_list','pseudogene_list','query_list','coord_list','cluster_number'))
        iac_positive_df=iac_positive_df[[len(set(i))>=6 for i in iac_positive_df['synteny'].str.findall('\w')]]
        iac_positive_df['contig'] = iac_positive_df['name'].str.contains('supercont|ctg|node|contig|scaffold|contigs',case=False)
        iac_positive_df['complete_genome']= iac_positive_df['name'].str.contains('complete',case=False)
        iac_positive_df['plasmid'] = iac_positive_df['name'].str.contains('plasmid')
        iac_positive_df['has_overlap'] = iac_positive_df['synteny_dir_dist'].str.contains("--")
        iac_positive_df['duplicated']= iac_positive_df.duplicated(subset="filename")
        print('Writing positive accessions to file')
        assembly_list = "\n".join(list(set(iac_positive_df['assembly']))) + "\n"
        with open("results/iac_positive_accessions.txt",'w+') as output:
            output.writelines(assembly_list)
        print('Fetching genome metadata')
        all_accs=list(set(iac_positive_df['accession']))
        all_biosamples = list(iac_positive_df['biosample'])
        print(all_biosamples)
        batch=150
        id_batches = [all_biosamples[i * batch:(i + 1) * batch] for i in range((len(all_biosamples) + batch - 1) // batch )]
        id_batch_str=[]
        for i in id_batches:
            id_list=re.sub('\n',',',str(i))
            id_list=re.sub(' ','',id_list)
            id_list=re.sub("NA\',|\'",'',id_list)
            id_list=re.sub("\[|\]","'",id_list)
            id_batch_str.append(id_list)
        isolation_source_dict={}
        isolation_source_dict["NA"] = "NA"
        env_biome_dict={}
        env_biome_dict["NA"] = "NA"
        env_feature_dict={}
        env_feature_dict["NA"] = "NA"
        host_dict={}
        host_dict["NA"] = "NA"
        host_sciname_dict={}
        host_sciname_dict["NA"] = "NA"
        strain_dict={}
        strain_dict["NA"] = "NA"
        all_meta_dict={}
        all_meta_dict["NA"]="NA"
        for i in id_batch_str:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=" + i + "&api_key=52553dfd9c090cfba1c3b28a45d8a648fd09"
            print(url)
            buffer = BytesIO()
            c = pycurl.Curl()
            c.setopt(c.URL, url)
            c.setopt(c.WRITEDATA, buffer)
            c.perform()
            c.close()
            body = buffer.getvalue()
            out=body.decode('iso-8859-1')
            root = ET.fromstring(out)
            for i in root.findall('./BioSample'):
                bs_text = ET.tostring(i).decode('ISO-8859-1')
                bs_acc=re.findall('(?<=\saccession\=")\w+',bs_text)[0]
                isolation_source= i.findall('./Attributes/Attribute[@attribute_name="isolation_source"]')
                if isolation_source:
                    isolation_source_dict[bs_acc] = isolation_source[0].text
                else:
                    isolation_source_dict[bs_acc] = "NA"
                env_biome= i.findall('./Attributes/Attribute[@attribute_name="env_biome"]')
                if env_biome:
                    env_biome_dict[bs_acc] = env_biome[0].text
                else:
                    env_biome_dict[bs_acc] = "NA"
                env_feature= i.findall('./Attributes/Attribute[@attribute_name="env_feature"]')
                if env_feature:
                    env_feature_dict[bs_acc] = env_feature[0].text
                else:
                    env_feature_dict[bs_acc] = "NA"
                host= i.findall('./Attributes/Attribute[@attribute_name="host"]')
                if host:
                    host_dict[bs_acc] = host[0].text
                else:
                    host_dict[bs_acc] = "NA"
                host_sciname= i.findall('./Attributes/Attribute[@attribute_name="host scientific name"]')
                if host_sciname:
                    host_sciname_dict[bs_acc] = host_sciname[0].text
                else:
                    host_sciname_dict[bs_acc] = "NA"
                strain= i.findall('./Attributes/Attribute[@attribute_name="strain"]')
                if strain:
                    strain_dict[bs_acc] = strain[0].text
                else:
                    strain_dict[bs_acc] = "NA"
                all_metadata = i.findall('./Attributes/Attribute[@attribute_name]')
                if all_metadata:
                    all_metadata_list=[]
                    for m in all_metadata:
                        all_metadata_list.append( (m.attrib['attribute_name'], m.text ))
                    all_meta_dict[bs_acc]= all_metadata_list
                else:
                    all_meta_dict[bs_acc] = "NA"
        print(isolation_source_dict)
        iac_positive_df['isolation_src'] = iac_positive_df['biosample'].map(isolation_source_dict)
        iac_positive_df['env_biome'] = iac_positive_df['biosample'].map(env_biome_dict)
        iac_positive_df['env_feature'] = iac_positive_df['biosample'].map(env_feature_dict)
        iac_positive_df['host'] = iac_positive_df['biosample'].map(host_dict)
        iac_positive_df['host_sciname'] = iac_positive_df['biosample'].map(host_sciname_dict)
        iac_positive_df['strain'] = iac_positive_df['biosample'].map(strain_dict)
        iac_positive_df['all_biosample_metadata'] = iac_positive_df['biosample'].map(all_meta_dict)
        all_assemblies= ['RS_'+ i for i in list(set(iac_positive_df['assembly']))]
        gtdb_metadata = pd.read_csv('gtdb/bac120_metadata_r89.tsv', sep='\t')
        gtdb_matches= gtdb_metadata[gtdb_metadata['accession'].isin(all_assemblies)]
        gtdb_matches.loc[:,'accession'] = gtdb_matches['accession'].str.replace('RS_','')
        gtdb_dict = gtdb_matches[['accession','gtdb_taxonomy']].to_dict()
        gtdb_dict = dict(zip(gtdb_matches['accession'], gtdb_matches['gtdb_taxonomy']))
        ncbi_dict = dict(zip(gtdb_matches['accession'], gtdb_matches['ncbi_taxonomy']))
        iac_positive_df['gtdb_tax'] = iac_positive_df['assembly'].map(gtdb_dict)
        iac_positive_df['ncbi_tax'] = iac_positive_df['assembly'].map(ncbi_dict)
        iac_positive_df['same_taxonomy'] = iac_positive_df['gtdb_tax'] == iac_positive_df['ncbi_tax']
        iac_positive_df[['domain_gtdb','phylum_gtdb','class_gtdb','order_gtdb','family_gtdb','genus_gtdb','species_gtdb']]=iac_positive_df['gtdb_tax'].str.split(";", expand = True)
        print('\nWriting iac cluster positive data frame to file')
        iac_positive_df.to_csv("results/iac_positive_df.tsv",sep = '\t', index = False)
        iac_positive_df.to_pickle("results/iac_positive_df.pickle")

rule index_files:
    input: rules.parse_blastp.output
    output: directory("index_files")
    threads: 1
    run:
        def make_indexprot(file, replace = False):
            #header_re = re.compile('^>.+', re.MULTILINE)
            prot_re = re.compile('^>[\s|\S]+?(?=\n\s*\n)', re.MULTILINE)
            print("attempting indexing for " + file)
            if glob.glob('fasta_files/'+file):
                filename = glob.glob('fasta_files/'+ file)[0]
                #outfilename = 'index_files/'+ file + ".index"
                outfilename = 'index_files/'+ file + ".indexprot"
                if os.path.exists(outfilename) and replace ==False:
                    print(outfilename + " already indexed")
                else:
                    with open(filename, 'r') as file:
                        file_contents = file.read()
                        proteins = prot_re.findall(file_contents)
                        out_prots = []
                        for p in proteins:
                            p1 = re.sub('NULL$','NULL!!',p)
                            p2 = re.sub('PSEUDOGENE$','PSEUDOGENE!!',p1)
                            p3 = p2.replace('\n','!!',1)
                            p4 = p3.replace('\n','??')
                            out_prots.append(p4)
                    with open(outfilename, 'w+') as outfile:
                        outfile.write('\n'.join(out_prots) + '\n')
                    print('made index for: ' + filename)
            else:
                print('No match for' + file)
        iac_positive_df=pd.read_pickle("results/iac_positive_df.pickle")
        inputs_indexprot = [re.sub('.gbff','_proteins.fa', i) for i in list(set(iac_positive_df['filename']))]
        # pool_index_prot = mp.Pool(mp.cpu_count()-2)
        # pool_index_prot.map(make_indexprot, inputs_indexprot)
        # pool_index_prot.close()
        os.mkdir("index_files")
        with concurrent.futures.ThreadPoolExecutor(max_workers=None) as executor:
            for i in executor.map(make_indexprot, inputs_indexprot):
                pass
        # with open("index_files/done") as file:
        #     file.write("done indexing")

rule fetch_neighborhoods:
    input: directory("index_files")
    output: "results/iac_positive_all_data.xlsx", "results/iac_positive_all_data_gtdb.xlsx", "results/iac_positive_all_data_gtdb.tsv"
    threads: 1
    run:
        def fetchneighborhood2(index):
            cluster = iac_positive_df.iloc[index,:]
            acc= cluster['accession']
            assembly = re.sub('.gbff','_proteins.fa.indexprot', cluster['filename'])
            print(acc + assembly)
            #make the genome database from the .fa.index file
            assembly_index_file = glob.glob('index_files/'+ assembly)[0]
            print(assembly_index_file)
            db = pd.read_csv(assembly_index_file, sep = "!!" ,header = None, engine='python')
            #db.columns = ["filename","assembly","accession","locus_tag","old_locus_tag","name","biosample","protein_name","coordinates","protein_id"]
            db.columns = ["filename","assembly","accession","locus_tag","old_locus_tag","name","biosample","protein_name","coordinates","protein_id","pseudogene","protein_seq"]
            db['direction']= [-1 if re.match('complement',c) else 1 for c in db['coordinates']]
            db['start_coord'] = [re.search('\d+?(?=\.\.(\d|\>))',str(c)).group(0) for c in db['coordinates'] ]
            db['start_coord'] = [re.sub('complement|>|<|\)|\(',"",c) for c in db['start_coord'] ]
            db['start_coord'] = db['start_coord'].astype(int)
            db['end_coord'] = [re.search('(?<=\.(\.|\>))\d+',str(c)).group(0) for c in db['coordinates'] ]
            db['end_coord'] = [re.sub('>|<|\)|\(',"",c) for c in db['end_coord'] ]
            db['end_coord'] = db['end_coord'].astype(int)
            features_upstream = 0
            features_downstream = 0
            hit_list = cluster['hit_list']
            query_list = cluster['query_list']
            cluster_number = cluster['cluster_number']
            hit_dict = dict(zip(hit_list,query_list))
            genome = db.loc[db['accession'] == acc].copy()
            start = genome[genome['locus_tag'] == hit_list[0]].index.values.astype(int)[0] - features_upstream
            stop = genome[genome['locus_tag'] == hit_list[-1]].index.values.astype(int)[0] + features_downstream
            neighborhood = genome.loc[start:stop,].copy()
            neighborhood['query_match'] = neighborhood['locus_tag'].map(hit_dict)
            coord_list = list(zip(neighborhood['start_coord'], neighborhood['end_coord'],neighborhood['direction'],neighborhood['query_match']))
            #function to find GC content of cluster vs genome
            gbff_str=str(db['filename'][0][1:])
            with open("gbff_files_unzipped/"+gbff_str) as file:
                gbff_file = file.read()
            genome_seq = "".join(re.findall("(?<=ORIGIN)[\s+\S+]+(?=\/\/)",gbff_file))
            Gg=genome_seq.count("g")
            Gc=genome_seq.count("c")
            Ga=genome_seq.count("a")
            Gt=genome_seq.count("t")
            genomeGC= (Gg+Gc)/(Gg+Gc+Ga+Gt)
            start=min(coord_list)[0]
            end=max(coord_list)[1]
            regex_str=acc+"[\s+\S+]+\/\/"
            all_cluster_fasta = re.findall(regex_str,gbff_file)[0]
            all_cluster_fasta = re.findall("(?<=ORIGIN)[\s+\S+]+(?=\/\/)",all_cluster_fasta)[0]
            all_cluster_fasta = re.sub(" |\d|\n","",all_cluster_fasta)
            cluster_seq = all_cluster_fasta[start-1:end-1]
            g=cluster_seq.count("g")
            c=cluster_seq.count("c")
            a=cluster_seq.count("a")
            t=cluster_seq.count("t")
            clusterGC = (g+c)/(g+c+a+t)
            diffGC = abs(clusterGC - genomeGC)
            ####
            if sum( neighborhood[neighborhood['query_match'].notnull()]['direction'] ) < 0:
                    neighborhood['actual_start_tmp'] = neighborhood['start_coord']
                    neighborhood['start_coord']= neighborhood['end_coord'] * -1
                    neighborhood['end_coord']= neighborhood['actual_start_tmp']* -1
                    neighborhood['direction'] = neighborhood['direction'] *-1
                    neighborhood = neighborhood.sort_values(by='start_coord')
            neighborhood['query_match'] = neighborhood['query_match'].replace(np.nan,"x")
            nhbrhood_hit_list= list(neighborhood['query_match'])
            nhbrhood_locus_tags= list(neighborhood['locus_tag'])
            nhbrhood_old_locus_tags= list(neighborhood['old_locus_tag'])
            nhbrhood_prot_ids= list(neighborhood['protein_id'])
            nhbrhood_prot_name= list(neighborhood['protein_name'])
            nhbrhood_prot_seq= list(neighborhood['protein_seq'])
            order = [ ("| " + gene['query_match'] + " 〉") if gene['direction'] == 1 else ("〈 " + gene['query_match'] + " |") for index,gene in neighborhood.iterrows()  ]
            dist = list(np.array(neighborhood['start_coord'][1:]) -  np.array(neighborhood['end_coord'][:-1]))
            dist = ["-"+str(d)+"-" for d in dist]
            adj_coord_list = list(zip(neighborhood['start_coord'], neighborhood['end_coord'],neighborhood['direction'],neighborhood['query_match']))
            if min(neighborhood['start_coord']) <0:
                tare_value = abs(min(neighborhood['start_coord']))
                tared_adj_coord_list = list(zip([v + tare_value for v in neighborhood['start_coord']], [v + tare_value for v in neighborhood['end_coord']],neighborhood['direction'],neighborhood['query_match']))
            else:
                tare_value = min(neighborhood['start_coord'])
                tared_adj_coord_list = list(zip([v - tare_value for v in neighborhood['start_coord']], [v - tare_value for v in neighborhood['end_coord']],neighborhood['direction'],neighborhood['query_match']))
            # making an ITOL compatible string
            gene_color_dict={ 'IacA':'#ff5969',
                                'IacB':'#2db34e',
                                'IacC':'#fb77e0',
                                'IacD':'#00bc7e',
                                'IacE':'#8d006e',
                                'IacF':'#cfdd63',
                                'IacG':'#0060d0',
                                'IacR':'#bb7b00',
                                'IacH':'#7c2c29',
                                'IacI':'#f1d17a',
                                'x':'#d1d1d1'}
            max_len = tared_adj_coord_list[-1][1]
            itol_diagram=[]
            for g in tared_adj_coord_list:
                gene_string=[]
                gene_length=g[1]-g[0]
                if g[2] > 0:
                    gene_string.append('RE')
                    gene_string.append(str(g[0]))
                    gene_string.append(str(g[1]-(0.1*gene_length)))
                    #gene_string.append('#34b4eb')
                    gene_string.append(gene_color_dict[g[3]])
                    gene_string.append(str(g[3]))
                    gene_string.append(',')
                    gene_string.append('TR')
                    gene_string.append(str(g[1]-(0.1*gene_length)))
                    gene_string.append(str(g[1]))
                    #gene_string.append('#34b4eb')
                    gene_string.append(gene_color_dict[g[3]])
                    gene_string.append('')
                else:
                    gene_string.append('TL')
                    gene_string.append(str(g[0]))
                    gene_string.append(str(g[0]+(0.1*gene_length)))
                    #gene_string.append('#34b4eb')
                    gene_string.append(gene_color_dict[g[3]])
                    gene_string.append('')
                    gene_string.append(',')
                    gene_string.append('RE')
                    gene_string.append(str(g[0]+(0.1*gene_length)))
                    gene_string.append(str(g[1]))
                    #gene_string.append('#34b4eb')
                    gene_string.append(gene_color_dict[g[3]])
                    gene_string.append(str(g[3]))
                itol_gene='|'.join(gene_string)
                itol_diagram.append(itol_gene)

            itol_diagram_joined=",".join(map(str,itol_diagram))
            itol_diagram_string= str(max_len)+','+itol_diagram_joined
            itol_diagram_string = re.sub(',\|',',',itol_diagram_string)
            #obtains "| A 〉-23-| B 〉-23-| C 〉"
            synteny_dir_dist=''.join(sum(zip(order, dist+[0]), ())[:-1])
            synteny_dir_dist = re.sub("Iac" ,"", synteny_dir_dist)
            #obtains "| A 〉| B 〉| C 〉"
            synteny_dir =''.join(order)
            synteny_dir = re.sub("Iac" ,"", synteny_dir)
            #obtains "| A:23.23 〉| B:23.23〉| C:23.23 〉"
            #synteny_dir_pident =''.join(order_pident)
            #synteny_dir_pident = re.sub("Iac" ,"", synteny_dir_pident)
            #obtains "A-B-C"
            synteny= re.sub("\n" ,"-", neighborhood['query_match'].to_string(index=False))
            synteny= re.sub("Iac| " ,"", synteny)
            cluster_len= max(neighborhood['end_coord']) - min(neighborhood ['start_coord'])
            assembly= re.sub("\{|\}|\'|>","", str(set(neighborhood['assembly'])) )
            accession = re.sub("\{|\}|\'","", str(set(neighborhood['accession'])) )
            title= re.sub("\{|\}|\'", "",str(set(neighborhood['name'])) )
            print(assembly_index_file + " successfully used")
            return([accession, assembly, title, len(neighborhood), cluster_len, synteny, synteny_dir_dist, synteny_dir, cluster_number,coord_list,adj_coord_list,tared_adj_coord_list,itol_diagram_string, nhbrhood_hit_list,nhbrhood_locus_tags,nhbrhood_old_locus_tags,nhbrhood_prot_ids,nhbrhood_prot_name,nhbrhood_prot_seq, clusterGC, genomeGC,diffGC, cluster_seq])
        iac_positive_df=pd.read_pickle("results/iac_positive_df.pickle")
        inputs_fetchneighborhood = list(range(0,len(iac_positive_df)))
        # pool_fetchneighborhood= mp.Pool(mp.cpu_count()-2)
        # pool_fetchneighborhood_output = pool_fetchneighborhood.map(fetchneighborhood2, inputs_fetchneighborhood)
        # pool_fetchneighborhood.close()
        outputs_fetchneighborhood=[]
        with concurrent.futures.ThreadPoolExecutor(max_workers=None) as executor:
            for i in executor.map(fetchneighborhood2, inputs_fetchneighborhood):
                outputs_fetchneighborhood.append(i)
                pass
        iac_pos_neighborhoods = list(filter(None, outputs_fetchneighborhood))
        iac_pos_neighborhoods_df= pd.DataFrame(iac_pos_neighborhoods, columns=('accession', 'assembly', 'title', 'feature_count_nhbr', 'cluster_len_nhbr', 'synteny_nhbr', 'synteny_dir_dist_nhbr', 'synteny_dir_nhbr','cluster_number','coord_list','adj_coord_list','tared_adj_coord_list','itol_cluster_string', 'nhbrhood_hit_list','nhbrhood_locus_tags','nhbrhood_old_locus_tags','nhbrhood_prot_ids','nhbrhood_prot_name','nhbrhood_prot_seq', 'clusterGC', 'genomeGC','diffGC','cluster_seq'))
        iac_pos_neighborhoods_df.to_csv("results/iac_pos_neighborhoods.tsv",sep = '\t', index = False)
        iac_pos_neighborhoods_df.to_pickle("results/iac_pos_neighborhoods.pickle")
        # merge the neighborhoods and iac_positive dataframes
        print('\nMerging all data and writing it to file')
        iac_positive_df = pd.read_pickle("results/iac_positive_df.pickle")
        iac_pos_neighborhoods_df = pd.read_pickle("results/iac_pos_neighborhoods.pickle")
        iac_positive_all = pd.merge(iac_pos_neighborhoods_df,iac_positive_df, on = ["accession","cluster_number","assembly"])
        iac_positive_all = iac_positive_all.rename(columns={"coord_list_y": "coord_list"}).drop(columns=['coord_list_x'])
        iac_positive_all_gtdb= iac_positive_all[iac_positive_all['gtdb_tax'].notnull()]
        iac_positive_all.to_pickle("results/iac_positive_all_data.pickle")
        iac_positive_all.to_csv("results/iac_positive_all_data.tsv", sep = '\t', index = False)
        iac_positive_all.to_excel("results/iac_positive_all_data.xlsx", index = False)
        iac_positive_all_gtdb.to_pickle("results/iac_positive_all_data_gtdb.pickle")
        iac_positive_all_gtdb.to_csv("results/iac_positive_all_data_gtdb.tsv", sep = '\t', index = False)
        iac_positive_all_gtdb.to_excel("results/iac_positive_all_data_gtdb.xlsx", index = False)
