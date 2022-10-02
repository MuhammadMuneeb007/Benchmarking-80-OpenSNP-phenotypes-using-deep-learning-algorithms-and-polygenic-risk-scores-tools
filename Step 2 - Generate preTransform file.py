import pandas as pd
import numpy as np
import os

# Main directory
maindirec = "opensnp_datadump.current"
# Phenotype FILE
pf = "phenotypes_202208291325.csv"
data = pd.read_csv(maindirec+os.sep+pf,sep=";")

# The first thing is to replace the comma and semi-colon from the data to parse the data properly.
data.replace(';','', regex=True, inplace=True)
data.replace(',','', regex=True, inplace=True)
data.replace('\t','', regex=True, inplace=True)


# Iterate each column and print the unique values and number of Nan.

# PhenotypeName
phenotypename = []

# Number of nan
numbernan = []

# unique values 
unique = []


columns = data.columns
# First 5 columns does not contain the phenotype.
columns = columns[5:]
count=0
valuesdataframe = pd.DataFrame()
phenovalue = []
phenovaluecount = []
for col in columns:
 data[col] = data[col].astype(str)
 phenotypename.append(col)
 numbernan.append(len(data)-data[col].str.contains(r'-').sum())
 unique.append(data[col].nunique())
 print(data[col].value_counts())
 phenovalue.append("Phenotype")
 phenovaluecount.append(col)
 phenovalue.extend(data[col].value_counts().tolist())
 phenovaluecount.extend(data[col].value_counts().index.tolist())
 #print(col)
 #print(data[col].value_counts().tolist())
 #print(data[col].value_counts().index.tolist())
 #if count==20:
 # exit(0)
 #count=count+1



 
valuesdataframe['phenovaluecount'] =  phenovalue
valuesdataframe['phenovalue'] = phenovaluecount
valuesdataframe.to_csv("preTransform.csv")

