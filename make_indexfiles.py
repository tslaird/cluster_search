
import glob
import pandas as pd
import re
import concurrent.futures
import os


def make_indexprot(file, replace = True):
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
                    p2 = p1.replace('\n','!!',1)
                    p3 = p2.replace('\n','??')
                    #substitution needed for database created 070919
                    p4 = re.sub("(?<!!)!{1}(?!!)","!!",p3)
                    out_prots.append(p4)
            with open(outfilename, 'w+') as outfile:
                outfile.write('\n'.join(out_prots) + '\n')
            print('made index for: ' + filename)
    else:
        print('No match for' + file)

iac_positive_df=pd.read_pickle("results_041320/iac_positive_df.pickle")
inputs_indexprot = [re.sub('.gbff','_proteins.fa', i) for i in list(set(iac_positive_df['filename']))]
# pool_index_prot = mp.Pool(mp.cpu_count()-2)
# pool_index_prot.map(make_indexprot, inputs_indexprot)
# pool_index_prot.close()
if not os.path.exists("index_files"):
    os.mkdir("index_files")

with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
    for i in executor.map(make_indexprot, inputs_indexprot):
        pass

def fetchneighborhood2(index):
    cluster = iac_positive_df.iloc[index,:]
    acc= cluster['accession']
    assembly = re.sub('.gbff','_proteins.fa.indexprot', cluster['filename'])
    print(acc +" "+ assembly)
    #make the genome database from the .fa.index file
    assembly_index_file = glob.glob('index_files/'+ assembly)[0]
    #print(assembly_index_file)
    db = pd.read_csv(assembly_index_file, sep = "!!" ,header = None, index_col=False, engine='python', names=["assembly","accession","locus_tag","old_locus_tag","name","biosample","protein_name","coordinates","protein_id","protein_seq"])
    #db.columns = ["filename","assembly","accession","locus_tag","old_locus_tag","name","biosample","protein_name","coordinates","protein_id"]
    #db.columns = ["assembly","accession","locus_tag","old_locus_tag","name","biosample","protein_name","coordinates","protein_id","protein_seq"]
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
    gbff_filename="gbff_files_unzipped/"+ cluster['filename']
    print(gbff_filename)
    with open(gbff_filename) as file:
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
    print(diffGC)
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
    return([accession, assembly, title, len(neighborhood), cluster_len, synteny, synteny_dir_dist, synteny_dir, cluster_number,coord_list,adj_coord_list,tared_adj_coord_list,itol_diagram_string, nhbrhood_hit_list,nhbrhood_locus_tags,nhbrhood_old_locus_tags,nhbrhood_prot_ids,nhbrhood_prot_name,nhbrhood_prot_seq, clusterGC, genomeGC ,diffGC, cluster_seq])



iac_positive_df=pd.read_pickle("results_041320/iac_positive_df.pickle")
inputs_fetchneighborhood = list(range(0,len(iac_positive_df)))
# pool_fetchneighborhood= mp.Pool(mp.cpu_count()-2)
# pool_fetchneighborhood_output = pool_fetchneighborhood.map(fetchneighborhood2, inputs_fetchneighborhood)
# pool_fetchneighborhood.close()
outputs_fetchneighborhood=[]
with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
    for i in executor.map(fetchneighborhood2, inputs_fetchneighborhood):
        outputs_fetchneighborhood.append(i)
        pass

iac_pos_neighborhoods = list(filter(None, outputs_fetchneighborhood))
iac_pos_neighborhoods_df= pd.DataFrame(iac_pos_neighborhoods, columns=('accession', 'assembly', 'title', 'feature_count_nhbr', 'cluster_len_nhbr', 'synteny_nhbr', 'synteny_dir_dist_nhbr', 'synteny_dir_nhbr','cluster_number','coord_list','adj_coord_list','tared_adj_coord_list','itol_cluster_string', 'nhbrhood_hit_list','nhbrhood_locus_tags','nhbrhood_old_locus_tags','nhbrhood_prot_ids','nhbrhood_prot_name','nhbrhood_prot_seq', 'clusterGC', 'genomeGC','diffGC','cluster_seq'))
iac_pos_neighborhoods_df.to_csv("results_041320/iac_pos_neighborhoods.tsv",sep = '\t', index = False)
iac_pos_neighborhoods_df.to_pickle("results_041320/iac_pos_neighborhoods.pickle")
# merge the neighborhoods and iac_positive dataframes
print('\nMerging all data and writing it to file')
iac_positive_df = pd.read_pickle("results_041320/iac_positive_df.pickle")
iac_pos_neighborhoods_df = pd.read_pickle("results_041320/iac_pos_neighborhoods.pickle")
iac_positive_all = pd.merge(iac_pos_neighborhoods_df,iac_positive_df, on = ["accession","cluster_number","assembly"])
iac_positive_all = iac_positive_all.rename(columns={"coord_list_y": "coord_list"}).drop(columns=['coord_list_x'])
iac_positive_all_gtdb= iac_positive_all[iac_positive_all['gtdb_tax'].notnull()]
iac_positive_all.to_pickle("results_041320/iac_positive_all_data.pickle")
iac_positive_all.to_csv("results_041320/iac_positive_all_data.tsv", sep = '\t', index = False)
iac_positive_all.to_excel("results_041320/iac_positive_all_data.xlsx", index = False)
iac_positive_all_gtdb.to_pickle("results_041320/iac_positive_all_data_gtdb.pickle")
iac_positive_all_gtdb.to_csv("results_041320/iac_positive_all_data_gtdb.tsv", sep = '\t', index = False)
iac_positive_all_gtdb.to_excel("results_041320/iac_positive_all_data_gtdb.xlsx", index = False)
