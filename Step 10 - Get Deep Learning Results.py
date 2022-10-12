import pandas as pd
import numpy as np
import os
from os.path import exists
import re
import pandas as pd
import numpy as np
import os
from os.path import exists
import re

def sorted_nicely( l ):
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)


total = 0
for loop in pd.read_csv("allphenotypesname2.txt",header=None)[0].values:
 count = 0
 merge = pd.DataFrame()
 for loop2 in range(1,6):
  f = os.listdir("./"+loop+os.sep+str(loop2)+os.sep)
  for ff in f:
   if "pv_" in ff:
    if exists("./"+loop+os.sep+str(loop2)+os.sep+ff+os.sep+"Results3.csv"):
     count=count+1
 if count==35:
  total = total+1
  print(loop) 
 else:
  pass
  print(count,loop)


for loop in pd.read_csv("allphenotypesname2.txt",header=None)[0].values:
 count = 0
 
 for loop2 in range(1,6):
  f = os.listdir("./"+loop+os.sep+str(loop2)+os.sep)
  for ff in f:
   if "pv_" in ff:
    if exists("./"+loop+os.sep+str(loop2)+os.sep+ff+os.sep+"Results3.csv"):
     count=count+1
 if count==35:
  for loop2 in range(1,6):
   merge = pd.DataFrame()
   f = os.listdir("./"+loop+os.sep+str(loop2)+os.sep)
   for ff in f:
    if "pv_" in ff:
     data = pd.read_csv("./"+loop+os.sep+str(loop2)+os.sep+ff+os.sep+"Results3.csv")
     merge = pd.concat([merge, data],axis=1)
     del data
   merge.to_csv("./"+loop+os.sep+str(loop2)+os.sep+"Results2.csv",sep ="\t",index=False)
   #print(merge.head())

 else:
  print(count,loop)

  
  
hu = pd.DataFrame()
disease = []
auc = []
SNP = []
STD = []
ALGO = []

for loop in pd.read_csv("allphenotypesname2.txt",header=None)[0].values:
 count=0
 for loop2 in range(1,6):
   if exists("./"+loop+os.sep+str(loop2)+os.sep+"Results2.csv"):
    count=count+1
 #print(loop,count)
 if count==5:
   shape = pd.read_csv("./"+loop+os.sep+str(loop2)+os.sep+"Results2.csv",sep="\t").shape
   average = np.zeros((shape[0], shape[1]))
   #try:
   all = []
   for loop2 in range(1,6):
    data = pd.read_csv("./"+loop+os.sep+str(loop2)+os.sep+"Results2.csv",sep="\t")
    data.columns = data.columns.str.replace(r"SNPs:", "")
    x = sorted_nicely(data.columns)
    
    data = data[list(x)]
    
    for col in data.columns:
     data[["Train"+col, 'Test'+col]] = data[col].str.split('-', expand=True)
     data['Test'+col] = pd.to_numeric(data['Test'+col], errors='coerce')
     del data[col]
     del data["Train"+col]
    all.append(data.values) 
    average = average + data.values
    #data = data.reindex_axis(sorted(data.columns, key=lambda x: float(x[1:])), axis=1)
   all = np.std((all[0],all[1],all[2],all[3],all[4]), axis=0, ddof=1)
   average = average/5
   result = np.where(average == np.amax(average))
   # Find the minimum standard deviation parameters
   
   aa =1000
   ind=0
   row = result[0]
   col = result[1]
   minrow=100
   mincol=100
   for xx in range(0,len(row)):
    #print(all[row[xx]][col[xx]])
    if aa>all[row[xx]][col[xx]]:
     aa=all[row[xx]][col[xx]]
     minrow =row[xx]
     mincol =col[xx]
   maximum = average[minrow][mincol]
   std = all[minrow][mincol]
   #print(maximum,std,"ML_"+str(minrow+1))
   disease.append(loop)
   auc.append(np.amax(average))
   STD.append(std)
   #print(result)
   x = np.array(x)
   SNP.append(x[mincol])
   #print(SNP)
   ALGO.append("DL_"+str(minrow+1))
   #print(ALGO)
   #except:
   #pass

print(len(disease),len(auc),len(STD),len(ALGO),len(SNP))
hu['Phenotype'] = disease
hu['Test AUC 5 Iterations Average'] = auc
hu['Standard Deviation'] = STD
hu['Deep learning algorithm index'] = ALGO
hu['Number of SNPs'] = SNP
hu.to_html("Deeplearningbasedbechmarking.html")
hu.to_csv("Deeplearningbasedbechmarking.csv",index=False,sep=",")

