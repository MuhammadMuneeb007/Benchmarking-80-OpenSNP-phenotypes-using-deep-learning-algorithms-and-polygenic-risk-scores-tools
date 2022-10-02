import pandas as pd
import numpy as np
import os

# Main directory
maindirec = "opensnp_datadump.current"
# Phenotype FILE
pf = "phenotypes_202208291325.csv"
data = pd.read_csv(maindirec+os.sep+pf,sep=";")

# Iterate each column and print the unique values and number of Nan.

# PhenotypeName
phenotypename = []

# Number of nan
numbernan = []

# unique values 
unique = []


columns = data.columns

# First 5 columns doe not contain the phenotype.
columns = columns[5:]
for col in columns:
 data[col] = data[col].astype(str)
 phenotypename.append(col)
 numbernan.append(len(data)-data[col].str.contains(r'-').sum())
 unique.append(data[col].nunique())

 
analysis = pd.DataFrame()
analysis['Phenotype Name'] = phenotypename
analysis['Number of People for which data is available'] = numbernan
analysis['Number of Unique Values'] = unique
analysis.to_html("Analysis1.html")
