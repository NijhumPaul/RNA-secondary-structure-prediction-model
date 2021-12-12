#!/usr/bin/env python
# coding: utf-8

# In[26]:


#percentage_of_parameters: Scrap the UniProt website 
#                         and calculates the parameters of each accession ID(bp, kpbp, hp, kphp, rc, kprc)
#protein_list_input_file: Read a CSV file, extract column 'Accession', 
#                         return parameter value of all accession IDs in a .csv file


# In[27]:


import re
from bs4 import BeautifulSoup
import requests
import pandas as pd
import urllib.request
import re
import logging


# In[28]:


from platform import python_version
print(python_version())


# In[29]:


'''All errors will be logged in a file named Error_logFile.txt
It will include each error's type, name and when the error occured''' 

logging.basicConfig(filename="Error_logFile.txt",
                   filemode='a',
                   format='%(asctime)s %(levelname)s-%(message)s',
                    datefmt='%Y-%m-%d ' '%H:%M:%S')


# In[30]:


'''
 This function scrap the UniProt website and calculates the 6 parameters for each accession ID (accession, kpbp, hp, kphp, rc, kprc)
  bp-> Betasheet percentage
  kpbp-> known partial betasheet percentage
  hp-> helix percentage
  kphp-> known partial helix percentage
  rc-> random coil percentage
  kprc-> known partial random coil percentage
  Input is accession ID, i.e. P35579
  Output is a tuple of 7 items (accession, bp, kpbp, hp, kphp, rc, kprc) of each accession ID'''

def percentage_of_parameters(accession):
    try:
        url_old="https://www.uniprot.org/uniprot/accession" 
        url_new=re.sub("accession", accession, url_old)
        headers = {
            'Access-Control-Allow-Origin': '*',
            'Access-Control-Allow-Methods': 'GET',
            'Access-Control-Allow-Headers': 'Content-Type',
            'Access-Control-Max-Age': '3600',
            'User-Agent': 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:52.0) Gecko/20100101 Firefox/52.0'
            }
        req = requests.get(url_new, headers)
        soup = BeautifulSoup(req.content, 'html.parser')
        #sec_structure = soup.find("div", {"id": "secondarystructure"})
        sec_structure = soup.find_all('text')[1]
        total_length = int(sec_structure.string)
        sec_table = soup.find_all('table', class_='featureTable', id="secstructure_section")
        each_feature_length_list = []
        feature_names_list = []
        for features in sec_table:
            rows = features.find_all('tr', class_='feature_row')
            for row in rows:
                feature_key = str(row.select_one('td').text)
                fts = re.sub(".*>", "", feature_key).strip()
                feature_names_list.append(str(fts))
                length = row.find('td', class_='featlength numeric').text
                each_feature_length_list.append(int(length))
        sum_all_features_length = sum(each_feature_length_list)
        df = pd.DataFrame(list(zip(feature_names_list, each_feature_length_list)),
                    columns =['feature_names_list', 'each_feature_length_list'])
        beta = df.loc[df['feature_names_list'] == 'Beta strandi', 'each_feature_length_list'].sum()
        bp = (beta/sum(each_feature_length_list))*100
        kpbp = (beta/total_length)*100
        alpha_helix = df.loc[df['feature_names_list'] == 'Helixi', 'each_feature_length_list'].sum()
        hp = (alpha_helix/sum(each_feature_length_list))*100
        kphp = (alpha_helix/total_length)*100
        random_coil = df.loc[df['feature_names_list'] == 'Turni', 'each_feature_length_list'].sum()
        rp = (random_coil/sum(each_feature_length_list))*100
        kprp = (random_coil/total_length)*100
        return accession, bp, kpbp, hp, kphp, rp, kprp
        
    except Exception as e:
        logging.error(e)
        bp = (float("nan"))
        kpbp = (float("nan"))
        hp = (float("nan"))
        kphp = (float("nan"))
        rp = (float("nan"))
        kprp = (float("nan"))
        return accession, bp, kpbp, hp, kphp, rp, kprp


# In[34]:


'''This function reads a CSV file and 
  extracts the column: Accession
  Returns a dataframe and a csv file containing all accession IDs and their corresponding parameters
  '''
def protein_list_input_file(file_name):
    try:
        parameters_list = []
        file = pd.read_csv(file_name)
        accession_clm = list(file["Accession"])
        for accesion_id in accession_clm:
            result = percentage_of_parameters(accesion_id)
            parameters_list.append(result)
        parameters_in_dataframe = pd.DataFrame(parameters_list, columns=['Accession', 'bp', 'kpbp', 'hp', 'kphp', 'rp', 'kprp'])
        #parameters_in_dataframe.to_csv("parameters_in_dataframe.csv")
        #print(parameters_in_dataframe)
        return parameters_in_dataframe.to_csv("Parameters_in_dataframe.csv", na_rep='NA')

    except Exception as e:
        logging.error(e)
        print(e)


# In[35]:


#protein_list_input_file("./MC3_BB6.csv")


# In[ ]:




