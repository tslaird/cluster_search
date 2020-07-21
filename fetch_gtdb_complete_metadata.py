import re
import sys
import pandas as pd
import pycurl
from io import BytesIO
import time
import xml.etree.ElementTree as ET
import numpy as np
import concurrent.futures
import urllib.request
import glob
import os
from sourmash import MinHash
import itertools
import scipy
from scipy import spatial

ncbi_api_key="52553dfd9c090cfba1c3b28a45d8a648fd09"
df=pd.read_excel('/home/tslaird/leveau_lab/cluster_search/GTDB_RS_complete_condensed.xlsx')
all_accs=list(set(df['nucleotide_id']))
#fetch biosample id from accession
# this works better than getting the metadata from the biosample name
# there were errors
batch=150
accs_batches = [all_accs[i * batch:(i + 1) * batch] for i in range((len(all_accs) + batch - 1) // batch )]
accs_batch_str=[]
gi_dict={}
biosample_ids_dict={}
for i in accs_batches:
    id_list=re.sub('\n',',',str(i))
    id_list=re.sub(' ','',id_list)
    id_list=re.sub("NA\',|\'",'',id_list)
    id_list=re.sub("\[|\]","'",id_list)
    id_list=re.sub(",","&id=",id_list)
    id_list=re.sub("'","",id_list)
    accs_batch_str.append(id_list)

for i in accs_batch_str:
    accs_list = re.findall("(?<=\=)\w+|^\w+",i)
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=biosample&id=" + i + "&api_key="+ncbi_api_key
    buffer = BytesIO()
    c = pycurl.Curl()
    c.setopt(c.URL, url)
    c.setopt(c.WRITEDATA, buffer)
    c.perform()
    c.close()
    body = buffer.getvalue()
    out=body.decode('iso-8859-1')
    root = ET.fromstring(out)
    for j in range(len(root.findall('./LinkSet'))):
        sgi= root.findall('./LinkSet')[j].findall('./IdList/Id')
        bsid= root.findall('./LinkSet')[j].findall('./LinkSetDb/Link/Id')
        accs_list[j]
        if bsid:
            biosample_ids_dict[str(accs_list[j])] = bsid[0].text
        else:
            biosample_ids_dict[str(accs_list[j])] = "NA"
        gi_dict[accs_list[j]] = sgi[0].text

df['nucleotide_id'] = df['nucleotide_id'].str.replace('\.\d+','')
df['biosample_id'] = df['nucleotide_id'].map(biosample_ids_dict)
df['gi'] = df['nucleotide_id'].map(gi_dict)
all_ids=list(df['biosample_id'])
batch=150
id_batches = [all_ids[i * batch:(i + 1) * batch] for i in range((len(all_ids) + batch - 1) // batch )]
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
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=" + i + "&api_key="+ncbi_api_key
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
        bs_id=re.findall('(?<=\sid\=")\w+',bs_text)[0]
        isolation_source= i.findall('./Attributes/Attribute[@attribute_name="isolation_source"]')
        if isolation_source:
            isolation_source_dict[bs_id] = isolation_source[0].text
        else:
            isolation_source_dict[bs_id] = "NA"
        env_biome= i.findall('./Attributes/Attribute[@attribute_name="env_biome"]')
        if env_biome:
            env_biome_dict[bs_id] = env_biome[0].text
        else:
            env_biome_dict[bs_id] = "NA"
        env_feature= i.findall('./Attributes/Attribute[@attribute_name="env_feature"]')
        if env_feature:
            env_feature_dict[bs_id] = env_feature[0].text
        else:
            env_feature_dict[bs_id] = "NA"
        host= i.findall('./Attributes/Attribute[@attribute_name="host"]')
        if host:
            host_dict[bs_id] = host[0].text
        else:
            host_dict[bs_id] = "NA"
        host_sciname= i.findall('./Attributes/Attribute[@attribute_name="host scientific name"]')
        if host_sciname:
            host_sciname_dict[bs_id] = host_sciname[0].text
        else:
            host_sciname_dict[bs_id] = "NA"
        strain= i.findall('./Attributes/Attribute[@attribute_name="strain"]')
        if strain:
            strain_dict[bs_id] = strain[0].text
        else:
            strain_dict[bs_id] = "NA"
        all_metadata = i.findall('./Attributes/Attribute[@attribute_name]')
        if all_metadata:
            all_metadata_list=[]
            for m in all_metadata:
                all_metadata_list.append( (m.attrib['attribute_name'], m.text ))
            all_meta_dict[bs_id]= all_metadata_list
        else:
            all_meta_dict[bs_id] = "NA"
#print(isolation_source_dict)
df['isolation_src'] = df['biosample_id'].map(isolation_source_dict)
df['env_biome'] = df['biosample_id'].map(env_biome_dict)
df['env_feature'] = df['biosample_id'].map(env_feature_dict)
df['host'] = df['biosample_id'].map(host_dict)
df['host_sciname'] = df['biosample_id'].map(host_sciname_dict)
df['strain'] = df['biosample_id'].map(strain_dict)
df['all_biosample_metadata'] = df['biosample_id'].map(all_meta_dict)
df['Main_category'] = ''
df['Subcategory_1'] = ''
df['Subcategory_2'] = ''



n = 2000  #chunk row size
list_df = [df[i:i+n] for i in range(0,df.shape[0],n)]

for i in range(0, len(list_df)):
    list_df[i].to_excel('/home/tslaird/leveau_lab/cluster_search/GTDB_RS_complete_condensed_metadata_'+str(i)+'.xlsx', index=False)
