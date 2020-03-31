def clusterGC()
iac_positive_df=pd.read_pickle("results/iac_positive_df.pickle")
inputs_fetchGC = list(range(0,len(iac_positive_df)))
cluster = iac_positive_df.iloc[index,:]
acc= cluster['accession']
assembly = cluster['filename']
with open("gbff_files_unzipped/"+assembly) as file:
    gbff_file = file.read()
all_fasta = re.findall( "(?<=ORIGIN)[\s+\S+]+(?=\/\/)",gbff_file)
for i in all_fasta:
    g=i.count("g")
    c=i.count("c")
    a=i.count("a")
    t=i.count("t")
    GCcontent= (g+c)/(g+c+a+t)
    print(GCcontent)
