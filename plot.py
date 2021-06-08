import streamlit as st
import pandas as pd
import functools
import numpy as np
import requests
import graphql.utilities as graphql


DATE_COLUMN = 'date/time'
DATA_URL = ('https://s3-us-west-2.amazonaws.com/'
         'streamlit-demo-data/uber-raw-data-sep14.csv.gz')



def load_data():
    data = requests.get(f'https://bioindex.hugeamp.org/api/bio/query/genes?q=slc30a8&fmt=row').json()['data']
    # build a dataframe
    dataframe = pd.DataFrame(data)
    st.write("Here's our first attempt at using data to create a table:")
    st.write(dataframe)
    return data



# Create a text element and let the reader know the data is loading.
data_load_state = st.text('Loading data...')
# Load 10,000 rows of data into the dataframe.
data = load_data()
# Notify the reader that the data was successfully loaded.
data_load_state.text('Loading data...using st cache!')
