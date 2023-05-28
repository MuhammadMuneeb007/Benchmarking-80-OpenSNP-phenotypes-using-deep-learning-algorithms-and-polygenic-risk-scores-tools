import pandas as pd
import numpy as np
import os


def Punctuation(string):

    # punctuation marks
    punctuations = '''!()-[]{};:'"\,<>./?@#$%^&*_~'''

    # traverse the given string and if any punctuation
    # marks occur replace it with null
    for x in string.lower():
        if x in punctuations:
            string = string.replace(x, "")

    # Print string without punctuation
    return string

# Create a final phenotype structure.
finalphenotypename = []
numberofcases = []
numberofcontrols = [] 
actualstringsforcasesandcontrols = []
totalnumberofsamples = []

# Main directory
maindirec = "opensnp_datadump.current"

# Phenotype FILE
pf = "phenotypes_202208291325.csv"
data = pd.read_csv(maindirec+os.sep+pf,sep=";")
data.replace(';','', regex=True, inplace=True)
data.replace(',','', regex=True, inplace=True)
data.replace('\t','', regex=True, inplace=True)

# Post pre-processing
postpheno = pd.read_csv("postTransform.csv",index_col=None)
postpheno = postpheno[['phenovaluecount','prephenovalue','postphenovalue']]
postpheno = postpheno[postpheno.postphenovalue != 'x']
#postpheno.dropna(subset=['prephenovalue'], inplace=True)
postpheno.to_csv("temp.csv")



columns = data.columns
# First 5 columns does not contain the phenotype.
columns = columns[5:]
count = 0
valuesdataframe = pd.DataFrame()
phenovalue = []
phenovaluecount = []
onlyphenopostpheno =  postpheno[postpheno.phenovaluecount == 'Phenotype']

# Read each phenotype name.
for col in columns:
 data[col] = data[col].astype(str)
 
 if col in onlyphenopostpheno['prephenovalue'].unique():
  
  temp = postpheno.head(data[col].nunique()+3)
  temp = temp.iloc[2:,:]
  actualtemp = pd.DataFrame()
  actualtemp['Name'] = data['user_id'].values
  actualtemp['File']= data['genotype_filename'].values
  actualtemp['dateofbith']= data['date_of_birth'].values
  actualtemp['Phenotype'] = data[col].values
  old = temp['prephenovalue'].values
  new = temp['postphenovalue'].values
  
  for loop in range(0,len(temp['prephenovalue'].values)):
   actualtemp['Phenotype'] = actualtemp['Phenotype'].replace(old[loop],new[loop])
 
   
   
  actualtemp = actualtemp[actualtemp.Phenotype != '-']
  actualtemp = actualtemp[actualtemp.Phenotype != 'U']
  actualtemp = actualtemp[actualtemp.Phenotype != '']
  actualtemp = actualtemp[actualtemp.Phenotype != 'nan']
  actualtemp = actualtemp[actualtemp['Phenotype'].notna()]
  ttp = Punctuation(col)
  ttp = ttp.replace(" ","")
  if not os.path.isdir(ttp): 
   os.mkdir(ttp)
  actualtemp.to_csv(ttp+os.sep+ttp+".csv",index=False)
  casesandcontrols = actualtemp['Phenotype'].value_counts()
  finalphenotypename.append(col)
  numberofcases.append(actualtemp.Phenotype.str.contains(actualtemp['Phenotype'].unique()[0]).sum())
  numberofcontrols.append(actualtemp.Phenotype.str.contains(actualtemp['Phenotype'].unique()[1]).sum())
  actualstringsforcasesandcontrols.append(actualtemp['Phenotype'].unique())
  totalnumberofsamples.append(int(actualtemp.Phenotype.str.contains(actualtemp['Phenotype'].unique()[0]).sum())+int(actualtemp.Phenotype.str.contains(actualtemp['Phenotype'].unique()[1]).sum()))
  
  finaldataframe = pd.DataFrame()
  finaldataframe['Phenotype Name'] = finalphenotypename
  finaldataframe['Number of Samples of Class 1'] = numberofcases
  finaldataframe['Number of Samples of Class 2'] = numberofcontrols
  finaldataframe['Total Number of Samples'] = totalnumberofsamples
  finaldataframe['Actual Classes'] = actualstringsforcasesandcontrols
  finaldataframe.to_html("Analysis2.html")
  import pdfkit
  pdfkit.from_file('Analysis2.html', 'Analysis2.pdf')
  postpheno = postpheno.iloc[data[col].nunique()+1: , :]



