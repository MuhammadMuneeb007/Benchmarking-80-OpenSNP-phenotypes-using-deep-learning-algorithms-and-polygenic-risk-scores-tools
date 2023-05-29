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

import sys
metric = sys.argv[1]
tool = sys.argv[3]

total = 0
f = []
foo = []
foo2 = []

for loop in pd.read_csv("allphenotypesname2.txt",header=None)[0].values:
 count = 0
 merge = pd.DataFrame()
 for loop2 in range(1,6):
  if exists("./"+loop+os.sep+str(loop2)+os.sep+"Results_"+tool+"_"+metric+".txt"):
   tempd = pd.read_csv("./"+loop+os.sep+str(loop2)+os.sep+"Results_"+tool+"_"+metric+".txt")
   if len(tempd)<675:
    foo.append(loop)
    foo2.append(loop2)
   count=count+1
  else:
   print(loop,loop2)
 if count!=5:
  pass
  #f.append(loop) 
 else:
  #print(count,loop)
  pass
#data = pd.DataFrame()
#data['A'] = foo
#data.to_csv("temppheno.txt",index=False,header=False)

#data['A'] = foo2
#data.to_csv("tempiteration.txt",index=False,header=False)
#exit(0)

data = pd.DataFrame()
data[0] = f
#data[0].to_csv("allphenotypesname3.txt",header=False,index=False) 
#exit(0)
import sys

hu = pd.DataFrame()
disease = []
auc = []
SNP = []
STD = []
ALGO = []
ss = 0
for loop in pd.read_csv("allphenotypesname2.txt",header=None)[0].values:
  count=0
  for loop2 in range(1,6):
   if exists("./"+loop+os.sep+str(loop2)+os.sep+"Results_"+tool+"_"+metric+".txt"):
    count=count+1
  #continue
  if count==5:
   shape = pd.read_csv("./"+loop+os.sep+str(loop2)+os.sep+"Results_"+tool+"_"+metric+".txt",sep=",").shape
   average = np.zeros((shape[0], 1))
  try:
   all = []
   
   for loop2 in range(1,6):
    data = pd.read_csv("./"+loop+os.sep+str(loop2)+os.sep+"Results_"+tool+"_"+metric+".txt",sep=",")
    temp = pd.DataFrame()
    temp[sys.argv[2]] = data[sys.argv[2]].values
    data = temp
    #data.columns = data.columns.str.replace(r"SNPs:", "")
    #data.columns = data.columns.str.replace(r"SNPs:", "")    
    x = sorted_nicely(data.columns)
    data = data[list(x)]
    all.append(data.values) 
    average = average + data.values
   ss=ss+1
   #continue
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
  except:
   pass

print(len(disease),len(auc),len(STD),len(ALGO),len(SNP))
hu['Phenotype'] = disease
hu['Test '+metric+' 5 Iterations Average'] = auc
hu['Test '+metric+' 5 Iterations Average'] = hu['Test '+metric+' 5 Iterations Average']*100

hu['Standard Deviation'] = STD
hu['Standard Deviation'] = hu['Standard Deviation']*100

hu['Best parameters index'] = ALGO
#hu['Number of SNPs'] = SNP
hu.to_html("PRS_"+sys.argv[2]+"_"+sys.argv[1]+"basedbechmarking.html")
hu.to_csv("PRS_"+sys.argv[2]+"_"+sys.argv[1]+"basedbechmarking.csv",index=False,sep=",")

