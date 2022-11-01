import pandas as pd
import numpy as np
import os
import sys
np.seterr(divide='ignore',invalid='ignore')

def modify(path1):
 pheno = pd.read_csv(path1+os.sep+"train.fam",sep="\s+",header=None,names=["a","b","c","d","e","f"])
 pheno.f[pheno['f']==1]=0
 pheno.f[pheno['f']==2]=1
 
 sample = pd.read_csv(path1+os.sep+"train.samples",sep="\s+")
 
 data = []
 data.append('B')
 data.extend(pheno["f"].values)
 sample['pheno'] = data
 sample.to_csv(path1+os.sep+"train.sample",sep="\t",index= False)

def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))

def makecommon(x,y):
 x['ID'] = x['CHR'].astype(str)+":"+x['BP'].astype(str)+"_"+x['A2']+"_"+x['A1']
 x['ID'] = x['ID'].astype(str)
 y['alternate_ids'] = y['alternate_ids'].astype(str)
 x= x.drop_duplicates(subset='ID', keep="last")
 y = y.drop_duplicates(subset='alternate_ids', keep="last")
 m = common(list(x['ID'].values), list(y['alternate_ids'].values))
 x =x[x['ID'].isin(m)]
 y =y[y['alternate_ids'].isin(m)]
 del x['ID']
 return x, y

import math
from scipy.stats import norm
from collections import Counter


def gwasmaker(path1):

 snpteststats = pd.read_csv(path1+os.sep+"train.sum",sep=" ",low_memory=False,skiprows = 10)
 snpteststats = snpteststats.head(len(snpteststats)-1)
 plinkstats = pd.read_csv(path1+os.sep+"train.assoc.fisher",sep="\s+",low_memory=False)
 plinkstats, snpteststats = makecommon(plinkstats,snpteststats)
 plinkstats, snpteststats = makecommon(plinkstats,snpteststats)

 # Find common SNPs between plink's GWAS and SNPTEST's GWAS.
 gwasstats = pd.DataFrame()
 gwasstats['CHR'] = plinkstats['CHR'].values
 gwasstats['BP'] = plinkstats['BP'].values
 gwasstats['SNP'] = plinkstats['SNP'].values
 gwasstats['A1'] = snpteststats['alleleA'].values
 gwasstats['A2'] = snpteststats['alleleB'].values
 gwasstats['N'] = snpteststats['all_total'].values
 gwasstats['SE'] = 1
 gwasstats['P'] = plinkstats['P'].values
 gwasstats['OR'] = plinkstats['OR'].values
 gwasstats['INFO'] = 1
 gwasstats['MAF'] = snpteststats['all_maf'].values
 gwasstats.to_csv(path1+os.sep+"plinkData.txt",index=False, sep="\t")
 gwasstats.to_csv(path1+os.sep+"plinkData.txt.gz",index=False, sep="\t",compression='gzip')
 print("DONE")
 '''
 # If you want to use LDpred-2, consider the following tool. LDpred-2 requires standard error calculation to find the heritability.
 a = plinkstats['OR'].astype(float)
 a =a.values
 b = plinkstats['P'].astype(float)
 b =b.values
 newvalues = []
 count = 0
 for loop in range(0,len(a)):
  try:
   if a[loop]==0:
    newvalues.append(-1)
   else:
    newvalues.append(float(abs(math.log(a[loop])/norm.ppf(b[loop]/2))))
  except Exception as e:
   if "math domain error" in e:
    newvalues.append(-1)
    count= count+1
    print(count,e)
   else:
    pass
 gwasstats = pd.DataFrame()
 gwasstats['CHR'] = plinkstats['CHR'].values
 gwasstats['BP'] = plinkstats['BP'].values
 gwasstats['SNP'] = plinkstats['SNP'].values
 gwasstats['A1'] = snpteststats['alleleA'].values
 gwasstats['A2'] = snpteststats['alleleB'].values
 gwasstats['N'] = snpteststats['all_total'].values
 gwasstats['SE'] = newvalues
 gwasstats['P'] = plinkstats['P'].values
 gwasstats['OR'] = plinkstats['OR'].values
 gwasstats['INFO'] = 1
 gwasstats['MAF'] = snpteststats['all_maf'].values
 gwasstats.replace([np.inf, -np.inf], np.nan, inplace=True)
 gwasstats.dropna(inplace=True)
 gwasstats =gwasstats[gwasstats['SE']!=-1]
 gwasstats.to_csv(path1+os.sep+"Data.txt",index=False, sep="\t")
 gwasstats.to_csv(path1+os.sep+"Data.txt.gz",index=False, sep="\t",compression='gzip')
 '''

pheno = [sys.argv[1].replace('\r', '')]
iteration = [sys.argv[2].replace('\r', '')]

# For each phenotype and iteration generate the GWAS summary statistic file.
for loop in pheno:
 for loop2 in iteration:
  try:
   path1 = "./"+loop+os.sep+str(loop2)+os.sep+"train"+os.sep
   
   # Generate GWAS using plink.
   os.system("./plink --bfile "+path1+"train --recode vcf --out "+path1+"train")
   # Generate vcf file for the training data, so GWAS can be generated from the SNPTEST.
   os.system("bcftools convert "+path1+"train.vcf -g   "+path1+"train")
   
   # Change the phenotype values. Plink requires phenotype to be 1 and 2. SNPTEST considers phenotype to be 0 and 1.
   modify(path1)
   os.system("./snptest -data "+path1+"train.gen.gz "+path1+"train.sample -o "+path1+"train.sum -frequentist 1 -method score -pheno pheno")
   # Merge GWAS.
   gwasmaker(path1)
  except:
   print("Error for (phenotype,iteration)",loop,loop2)
   pass


