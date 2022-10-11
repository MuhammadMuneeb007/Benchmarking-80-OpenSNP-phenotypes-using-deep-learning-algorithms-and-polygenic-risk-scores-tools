import pandas as pd
import numpy as np
import os
import shutil
import re 
import glob


finalphenotypename = []
numberofcases = []
numberofcontrols = [] 
actualstringsforcasesandcontrols = []
totalnumberofsamples = []



def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)







allfiles = os.listdir("./")
allfiles = sorted_nicely(allfiles)
count = 1
for loop in allfiles:
 if ".sh" not in loop and "PRSice" not in loop and "plink" not in loop and "gtool" not in loop and "snptest" not in loop and ".err" not in loop and ".out" not in loop and ".txt" not in loop and "23andme" not in loop and ".csv" not in loop and ".py" not in loop and ".html" not in loop and "opensnp_datadump.current" not in loop:
  if "final.QC.bed" in os.listdir(loop+"/"):
   #print("\n")
   #print(count,loop,len(os.listdir(loop+"/")))
   #print(count,loop)
   #print(os.listdir(loop+"/"))
   #print(loop)
   
   try:
    phenotypefile2 = pd.read_csv(loop+os.sep+"Phenotype_process2.csv",sep="\t")
    phenotypefile1 = pd.read_csv(loop+os.sep+"Phenotype_process1.csv",sep=",")
    phenotypefile2['Phenotype'] = phenotypefile2['Phenotype'].astype(str)
    phenotypefile1['Phenotype'] = phenotypefile1['Phenotype'].astype(str)
    totalnumberofsamples.append(int(phenotypefile2.Phenotype.str.contains(phenotypefile2['Phenotype'].unique()[0]).sum())+int(phenotypefile2.Phenotype.str.contains(phenotypefile2['Phenotype'].unique()[1]).sum()))
   
    finalphenotypename.append(loop)
    numberofcases.append(phenotypefile2.Phenotype.str.contains(phenotypefile2['Phenotype'].unique()[0]).sum())
    numberofcontrols.append(phenotypefile2.Phenotype.str.contains(phenotypefile2['Phenotype'].unique()[1]).sum())
    actualstringsforcasesandcontrols.append(phenotypefile1['Phenotype'].unique())
    
   except:
    continue
   print(loop,phenotypefile2['Phenotype'].unique())

   count=count+1
   

print(count)
finaldataframe = pd.DataFrame()
finaldataframe['Phenotype Name'] = finalphenotypename
finaldataframe['Number of Samples of Class 1'] = numberofcases
finaldataframe['Number of Samples of Class 2'] = numberofcontrols
finaldataframe['Total Number of Samples'] = totalnumberofsamples
finaldataframe['Actual Classes'] = actualstringsforcasesandcontrols
finaldataframe.to_html("analysis3.html")
