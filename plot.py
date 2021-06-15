import streamlit as st
import pandas as pd
import functools
import numpy as np
import requests
import graphql.utilities as gl

def query(q, concat=False):
    resp = requests.post(f'https://bioindex.hugeamp.org/api/bio/query', data=q)
    if resp.status_code != 200:
        st.write(resp)
        raise RuntimeError(resp.json()['detail'])
        #concatenate all the results together into a single frame
    if concat:
        return pd.DataFrame(itertools.chain.from_iterable(resp.json()['data'].values()))
            #values are in insertion order of the query
#        st.write(pd.DataFrame())
    for k, rs in resp.json()['data'].items():
#        st.write(pd.DataFrame(rs))
        x = {k:pd.DataFrame(rs)}
    
    return {k: pd.DataFrame(rs) for k, rs in resp.json()['data'].items()}
        
#Create a text element and let the reader know the data is loading.
data_load_state = st.text('Loading data...')
q = """{GlobalEnrichment(phenotype: "ldl") {annotation, tissue, SNPs, expectedSNPs}}"""
 # Load 10,000 rows of data into the dataframe.
data = query(q)
# Notify the reader that the data was successfully loaded.
data_load_state.text('Loading data...using st cache!')
