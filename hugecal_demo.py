"""
Streamlit Cheat Sheet
App to summarise streamlit docs v0.81.0 for quick reference
There is also an accompanying png version
https://github.com/daniellewisDL/streamlit-cheat-sheet
v0.71.0 November 2020 Daniel Lewis and Austin Chen
"""

import streamlit as st
from pathlib import Path
import altair as alt
import base64
import functools
import requests
import graphql.utilities as graphql
import itertools
import pandas as pd
import re
import urllib.parse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
# Initial page config

st.set_page_config(
     page_title='HUGE Calculator',
     layout="wide",
     initial_sidebar_state="expanded",
)

gene_input = st.sidebar.text_input("Enter a gene","SLC30A8").upper()
phenotype_input = "T2D"
prior_input = 0.3696
isGwasSignificant = False
isExomewideSignificant = False
def main():
    cs_sidebar()
    #get_genes()
    cs_body()
    #bf_commonvariation()
    
    return None

# Thanks to streamlitopedia for the following code snippet

def img_to_bytes(img_path):
    img_bytes = Path(img_path).read_bytes()
    encoded = base64.b64encode(img_bytes).decode()
    return encoded

# sidebarAssociations

def cs_sidebar():
    st.sidebar.header('HuGE Calculator')
    st.sidebar.markdown('''<small>https://hugeamp.org/hugecalculator.html.</small>''', unsafe_allow_html=True)
    #Dropdown menus
    
    phenotype_input = st.sidebar.selectbox("Select a phenotype", (get_phenotypes()))
    
    return None
    
def get_genes():
    data = requests.get(f'https://bioindex.hugeamp.org/api/bio/match/gene?q=s').json()['data']
    # build a dataframe
    dataframe = pd.DataFrame(data)
    st.write("Here's our first attempt at using data to create a table:")
    st.write(dataframe)
    return data
    
def get_52kassociations():
    data =requests.get(f'https://bioindex.hugeamp.org/api/bio/query/gene-associations-52k?q=slc30a8&fmt=row').json()['data']
    dataframe_52kassoc = pd.DataFrame(data)
    return dataframe_52kassoc
    


def get_phenotypes():
    dataframe_52kassoc = get_52kassociations()
    return dataframe_52kassoc["phenotype"].tolist();


def graph_client(q, concat=False):
    resp = requests.post(f'https://bioindex.hugeamp.org/api/bio/query', data=q)
    if resp.status_code != 200:
        st.write(resp)
        raise RuntimeError(resp.json()['detail'])
        #concatenate all the results together into a single frame
    if concat:
        return pd.DataFrame(itertools.chain.from_iterable(resp.json()['data'].values()))
            #values are in insertion order of the query
#        st.write(pd.DataFrame())
    return {k: pd.DataFrame(rs) for k, rs in resp.json()['data'].items()}
    
def get_Category(bayes_factor):
    category = ""
   
    if bayes_factor <=1:
        category = "No"
    if bayes_factor >1 and bayes_factor < 3:
        category = "Anecdotal"
    if bayes_factor >=3 and bayes_factor <10:
        category = "Moderate"
    if bayes_factor >=10 and bayes_factor <30:
        category = "Strong"
    if bayes_factor >=30 and bayes_factor <100:
        category = "Very Strong"
    if bayes_factor >= 100 and bayes_factor < 350:
        category = "Extreme"
    if bayes_factor >= 350:
        category = "Compelling"
    return category
   


def getEGLData(dataset,trait):
    dataset = f'"{dataset}"'
    trait = f'"{trait}"'
    egl_data =requests.get("""https://kp4cd.org/egldata/dataset?dataset=mccarthy&trait=t2d""").json()['data']
    return pd.DataFrame(egl_data)

def isGwasSignificant(df):
    for index, row in df.iterrows():
        if row['phenotype'] == phenotype_input:
            if row['pValue'] <= 5e-8:
                isGwasSignificant = True
                return isGwasSignificant
                
  

def bf_commonvariation():
    gene = f'"{gene_input}"'
    phenotype = f'"{phenotype_input}"'
    firstBF = 1
    secondBF = 1
    thirdBF = 1
    commonBF = 1
    firstEvidence = ""
    secondEvidence = ""
    thirdEvidence = ""
    q = "{Associations(locus:" + gene + ",phenotype:" + phenotype + ") {pValue,phenotype}}"
    assoc_data = graph_client(q)
    assoc_data_df = assoc_data["Associations"]
    egl_data_df = getEGLData("mccarthy",trait="t2d")
    
    if isGwasSignificant(assoc_data_df):
        firstEvidence = "Yes (P-value <= 5e-8)"
        firstBF = 3
        for index,row in egl_data_df.iterrows():
            if row['gene'] == gene_input:
            
                if(row.genetic == "1C"):
                    secondEvidence = "1C"
                    secondBF = 117
                if(row.genetic == "2C"):
                    secondEvidence = "2C"
                    secondBF = 5
                if(row.genomic == "2R"):
                    thirdEvidence = "2R"
                    thirdBF = 5
                if(row.genomic == "3R"):
                    thirdEvidence = "3R"
                    thirdBF = 2.2
                    
    else:
        firstEvidence = "Not GWAS Significant"
       
    
    commonBF = firstBF * secondBF * thirdBF
    commonCategory = get_Category(commonBF)
 
#        st.write(f"{firstBF}" + " * " +  f"{secondBF}" + " * " + f"{thirdBF}" )
    d = {'GWAS Evidence': [firstEvidence,firstBF ],
            'Coding Evidence':[secondEvidence, secondBF],
            'Regulatory Evidence':[thirdEvidence, thirdBF],
            'Bayes Factor':["",commonBF],
            'Category':[commonCategory,""]
          
            }
   
    common_bayes_score = pd.DataFrame(d)
    return common_bayes_score
        
        
def isExomewideSignificant(df):
    for index, row in df.iterrows():
        if row['phenotype'] == phenotype_input:
            if row['pValue'] <= 2.5e-6:
                isExomewideSignificant = True
                return isExomewideSignificant
    
    
def rare_variation_strength(prior,beta,stdErr):
    w = prior
    v = stdErr**2
    f1 = v/(v+w)
    sqrt_f1 =  math.sqrt(f1)
    f2 = w * math.pow(beta,2)
    f3 =  2 * v * (v + w)
    f4 = f2/f3
    bayes_factor = sqrt_f1 * math.exp(f4)
    
    return bayes_factor


def find_most_significantmask(df):
    for index, row in df.iterrows():
        if row['phenotype'] == phenotype_input:
            masks = row['masks']
            newlist = sorted(masks, key=lambda k: k['pValue'])
            #st.write(pd.DataFrame(newlist))
            return newlist[0]
                   
        
        
def bf_rarevariation(prior):
    gene = f'"{gene_input}"'
    q = "{GeneAssociations52k(gene:" + gene + ") {pValue,phenotype,beta,dataset,gene, masks {beta, mask, pValue, stdErr, combinedAF}}}"
    assoc_52kdata = graph_client(q)
    assoc_52kdata_df = assoc_52kdata["GeneAssociations52k"]
    #st.write(assoc_52kdata_df)
    rareBF = 1
    rareEvidence = ""
    
    most_siginificant_mask = find_most_significantmask(assoc_52kdata_df)
    stdErr = most_siginificant_mask['stdErr']
    beta = most_siginificant_mask['beta']
   
    if isExomewideSignificant(assoc_52kdata_df):
        rareBF = 1650
        rareEvidence = "Exome wide significant"
    else:
        rareBF = rare_variation_strength(prior_input,beta,stdErr)
        rareEvidence = "Not Exome wide significant"
    
    rareCategory = get_Category(rareBF)
    d = {'Exome Evidence': [rareEvidence,rareBF ],
            'Bayes Factor':["",rareBF],
            'Category':[rareCategory,""]
            }
   
    rare_bayes_score = pd.DataFrame(d)
    return rare_bayes_score


    
    
    
    
#1. input for Prior inside my_expander
#2. with change in prior, bayes factor for rare variation should change
#3. draw a plot using this prior and ppr using plotly

        
##########################
# Main body of Huge Cal
##########################


def cs_body():
    gene = f'"{gene_input}"'
    phenotype = f'"{phenotype_input}"'
    q = "{Associations(locus:" + gene + ",phenotype:" + phenotype + ") {pValue,phenotype}}"
    assoc_data = graph_client(q)
    assoc_data_df = assoc_data["Associations"]
    
    q = "{GeneAssociations52k(gene:" + gene + ") {pValue,phenotype,beta,dataset,gene, masks {beta, mask, pValue, stdErr, combinedAF}}}"
    assoc_52kdata = graph_client(q)
    assoc_52kdata_df = assoc_52kdata["GeneAssociations52k"]
    
    
    st.header("Human Genetic Evidence for" + " "+ gene_input + " in " + phenotype_input)
    bf_combinedvariation = bf_rarevariation(0.3696) * bf_commonvariation()["Bayes Factor"][1]
    # Combined Evidence
    with st.beta_container():
        st.subheader("Combined Evidence")
        st.write(pd.DataFrame({
            'Bayes factor for Common variation':[bf_commonvariation()["Bayes Factor"][1]],
            'Bayes Factor for Rare Variation':[bf_rarevariation(0.3696)["Bayes Factor"][1]],
            'Combined Bayes Factor':[bf_combinedvariation["Bayes Factor"][1]],
            'Combined Evidence':[get_Category(bf_combinedvariation["Bayes Factor"][1])]
                              }))
        combined_expander = st.beta_expander("See explaination")
        #prior_input = combined_expander.text_input('change prior',0.3696)
        
    col1, col2 = st.beta_columns(2)
    
    #Common Evidence
    col1.subheader('Common Evidence')
    if isGwasSignificant(assoc_data_df):
        col1.code("* " + f'{gene_input}' +  " is GWAS siginificant" + "\n* " + bf_commonvariation()["Category"][0] + " Evidence" + "\n* Bayes Score = " + str(bf_commonvariation()["Bayes Factor"][1]) )
    else:
        col1.code("* " + f'{gene_input}' +  " is not GWAS siginificant" + "\n* "  + bf_commonvariation()["Category"][0] + " Evidence" + "\n* Bayes Score = " + str(bf_commonvariation()["Bayes Factor"][1]) )
    commonvariation_expander = col1.beta_expander("See explaination")
    commonvariation_expander.code(bf_commonvariation())
    
   
   
    
    df = pd.DataFrame(
    {
        'prior': [0, 0.001, 0.002, 0.003, 0.004],
        'ppr': [10,20,30,40,50],
    },
    columns=['prior', 'ppr']
)
   

#    st.write(df)
    chart = alt.Chart(df).mark_line().encode(
    x=alt.X('prior:N'),
    y=alt.Y('ppr:Q'),
    color=alt.Color("name:N")).properties(title="Posterior probability vs Prior")
    
#    #Rare Variation
#    col2.subheader('Rare Evidence')
#    if isExomewideSignificant(assoc_52kdata_df):
#        col2.code('''Bayes factor = 1 (Exome significant)''')
#
#    else:
##        col2.st.text_input("Enter Prior",0.3696)
#        my_expander2 = col2.beta_expander("See explaination")
#        my_expander2.altair_chart(chart, use_container_width=True)

    #Common Evidence
    col2.subheader('Rare Evidence')
    if isExomewideSignificant(assoc_52kdata_df):
        col2.code("* " + f'{gene_input}' +  " is Exome wide siginificant" + "\n* " + bf_rarevariation(0.3696)["Category"][0] + " Evidence" + "\n* Bayes Score = " + str(bf_rarevariation(0.3696)["Bayes Factor"][1]) )
    else:
        col2.code("* " + f'{gene_input}' +  " is not GWAS siginificant" + "\n* "  + bf_rarevariation(0.3696)["Category"][0] + " Evidence" + "\n* Bayes Score = " + str(bf_rarevariation(0.3696)["Bayes Factor"][1]) )
    rarevariation_expander = col2.beta_expander("See explaination")
    rarevariation_expander.code(bf_rarevariation(0.3696))
    

# Run main()

if __name__ == '__main__':
    main()

