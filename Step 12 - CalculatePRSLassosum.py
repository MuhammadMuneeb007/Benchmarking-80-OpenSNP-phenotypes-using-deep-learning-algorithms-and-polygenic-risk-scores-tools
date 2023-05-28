from cgi import test
import pandas as pd
import os
import glob
import numpy as np
import shutil
import sys
import numpy as np
import scipy as sp
from scipy import stats
import os
import pandas as pd
import numpy as np
import math
from sklearn.metrics import roc_auc_score
import numpy as np
import scipy.stats as stats
import math
import numpy as np
from scipy.stats import norm
import statistics
from sklearn.metrics import matthews_corrcoef

def Average(lst):
    return sum(lst) / len(lst)
def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def calculateMCClasso(direc):
  best = pd.read_csv(direc+os.sep+"result"+os.sep+"test.txt",sep="\s+",header=None)
  pheno = pd.read_csv(direc+os.sep+"test/test.fam",sep="\s+",header=None,names=["a","b","c","d","e","f"])
  pheno.f[pheno['f']==1]=0
  pheno.f[pheno['f']==2]=1

  if all(item == 0 for item in best[0].values) or all(item == 1 for item in best[0].values):
   best = best.fillna(0)
   best[0] = NormalizeData(best[0].values)
   best = best.fillna(0)
   best[0].values[best[0] >=0.5] = 1
   best[0].values[best[0] < 0.5] = 0
   print(best[0].values)
   print(pheno['f'].values)
   return matthews_corrcoef(best[0].values,pheno['f'].values)
  else:
   best = best.fillna(0)
   best[0] = NormalizeData(best[0].values)
   best = best.fillna(0)
   best[0].values[best[0] >=0.5] = 1
   best[0].values[best[0] < 0.5] = 0
   print(list(best[0].values))
   print(list(pheno['f'].values))
   auc = roc_auc_score(best[0].values,pheno['f'].values)
   mcc = matthews_corrcoef(best[0].values,pheno['f'].values)
   return mcc

def calculateAUClasso(direc):
  best = pd.read_csv(direc+os.sep+"result"+os.sep+"test.txt",sep="\s+",header=None)
  pheno = pd.read_csv(direc+os.sep+"test/test.fam",sep="\s+",header=None,names=["a","b","c","d","e","f"])
  pheno.f[pheno['f']==1]=0
  pheno.f[pheno['f']==2]=1  
  
  if all(item == 0 for item in best[0].values) or all(item == 1 for item in best[0].values):
   best[0].values[best[0] >=0.5] = 1
   best[0].values[best[0] < 0.5] = 0
   return 0.5
  else:
   best[0] = best[0].fillna(0)
   best[0] = NormalizeData(best[0].values)
   best =  best.fillna(0)
   best[0].values[best[0] >=0.5] = 1
   best[0].values[best[0] < 0.5] = 0
   auc = roc_auc_score(best[0].values,pheno['f'].values)
   return auc

def calculateAUCPRScice(path,test):
  best = pd.read_csv(path+os.sep+"result"+os.sep+"PRScice_PRS.best",sep="\s+")
  pheno = pd.read_csv(test+os.sep+"test.fam",sep="\s+",header=None,names=["a","b","c","d","e","f"])
  pheno.f[pheno['f']==1]=0
  pheno.f[pheno['f']==2]=1
  if all(item == 0 for item in best["PRS"].values) or all(item == 1 for item in best["PRS"].values):
   best["PRS"].values[best["PRS"] >=0.5] = 1
   best["PRS"].values[best["PRS"] < 0.5] = 0
   return 0
  else:
   best["PRS"] = NormalizeData(np.array(best["PRS"].values))
   best["PRS"].values[best["PRS"] >=0.5] = 1
   best["PRS"].values[best["PRS"] < 0.5] = 0
   roc_auc_score(best["PRS"].values, pheno['f'].values)

  return roc_auc_score(best["PRS"].values, pheno['f'].values)

def calculateMCCPRScice(path,test):
  best = pd.read_csv(path+os.sep+"result"+os.sep+"PRScice_PRS.best",sep="\s+")
  pheno = pd.read_csv(test+os.sep+"test.fam",sep="\s+",header=None,names=["a","b","c","d","e","f"])
  pheno.f[pheno['f']==1]=0
  pheno.f[pheno['f']==2]=1
  if all(item == 0 for item in best["PRS"].values) or all(item == 1 for item in best["PRS"].values):
   best["PRS"].values[best["PRS"] >=0.5] = 1
   best["PRS"].values[best["PRS"] < 0.5] = 0
   return 0
  else:
   best["PRS"] = NormalizeData(np.array(best["PRS"].values))
   best["PRS"].values[best["PRS"] >=0.5] = 1
   best["PRS"].values[best["PRS"] < 0.5] = 0
   roc_auc_score(best["PRS"].values, pheno['f'].values)
  return   int(matthews_corrcoef(best["PRS"].values, pheno['f'].values))




def calculateMCCPlink(direc):
  pheno = pd.read_csv(direc+os.sep+"/test/YRI.pheno",sep="\s+")
  pheno.phenotype[pheno['phenotype']==1]=0
  pheno.phenotype[pheno['phenotype']==2]=1
  files = os.listdir(direc+os.sep+"files")
  maxxacc=0
  maxauc = 0
  profile = ""
  temp = []
  for loop in files:
    if ".profile" in loop:      
      best = pd.read_csv(direc+os.sep+"files"+os.sep+loop,sep="\s+")
      #plt.hist(best["PRS"].values)
      #plt.show()

      #best["PRS"].values[best["PRS"] >= sum(best["PRS"].values)/len(best["PRS"].values)] = 1
      #best["PRS"].values[best["PRS"] < sum(best["PRS"].values)/len(best["PRS"].values)] = 0
      best["SCORE"] = NormalizeData(best["SCORE"].values)

      best["SCORE"].values[best["SCORE"] >=0.5] = 1
      best["SCORE"].values[best["SCORE"] < 0.5] = 0

      if maxauc<matthews_corrcoef(best["SCORE"].values, pheno['phenotype'].values):
        maxauc = matthews_corrcoef(best["SCORE"].values, pheno['phenotype'].values)
        temp = NormalizeData(pd.read_csv(direc+os.sep+"files"+os.sep+loop,sep="\s+")["SCORE"].values)
        profile = loop
        a= best["SCORE"].values
        b = pheno['phenotype'].values
        accuracy = len([a[i] for i in range(0, len(a)) if a[i] == b[i]]) / len(a)
        maxxacc = accuracy
  return  maxauc

def calculateAUCPlink(direc):
  pheno = pd.read_csv(direc+os.sep+"/test/YRI.pheno",sep="\s+")
  pheno.phenotype[pheno['phenotype']==1]=0
  pheno.phenotype[pheno['phenotype']==2]=1

  files = os.listdir(direc+os.sep+"files")
  #print(pheno.head())

  maxxacc=0
  maxauc = 0
  profile = ""
  temp = []
  for loop in files:

    if ".profile" in loop:      
      best = pd.read_csv(direc+os.sep+"files"+os.sep+loop,sep="\s+")
      #plt.hist(best["PRS"].values)
      #plt.show()
    
      #best["PRS"].values[best["PRS"] >= sum(best["PRS"].values)/len(best["PRS"].values)] = 1
      #best["PRS"].values[best["PRS"] < sum(best["PRS"].values)/len(best["PRS"].values)] = 0
      best["SCORE"] = NormalizeData(best["SCORE"].values)
    
      best["SCORE"].values[best["SCORE"] >=0.5] = 1
      best["SCORE"].values[best["SCORE"] < 0.5] = 0
      
      if maxauc<roc_auc_score(best["SCORE"].values, pheno['phenotype'].values):
        maxauc = roc_auc_score(best["SCORE"].values, pheno['phenotype'].values)
        temp = NormalizeData(pd.read_csv(direc+os.sep+"files"+os.sep+loop,sep="\s+")["SCORE"].values)
        profile = loop
        a= best["SCORE"].values
        b = pheno['phenotype'].values
        accuracy = len([a[i] for i in range(0, len(a)) if a[i] == b[i]]) / len(a)
        maxxacc = accuracy

  return  maxauc




def changeIDS(direct):
  data = pd.read_csv(direct, sep="\s+",index_col=False)
  data['FID'] = data['FID'].str.split('_').str[0]
  data['IID'] = data['IID'].str.split('_').str[0]
  data.to_csv(direct,sep="\t",index=False)




def gwasquality(traindirec,fdirec):
 os.system("gunzip -c ./"+traindirec+os.sep+"Data.txt.gz | awk 'NR==1 || ($11 > 0.01) && ($10 > 0.8) {print}' | gzip  > ./"+fdirec+os.sep+"Data.gz") 
 os.system("gunzip -c ./"+fdirec+os.sep+"Data.gz | awk '{seen[$3]++; if(seen[$3]==1){ print}}' | gzip -> ./"+fdirec+os.sep+"Data.nodup.gz")
 os.system("gunzip -c ./"+fdirec+os.sep+"Data.nodup.gz | awk '!( ($4==\"A\" && $5==\"T\") || ($4==\"T\" && $5==\"A\") || ($4==\"G\" && $5==\"C\") || ($4==\"C\" && $5==\"G\")) {print}' | gzip > ./"+fdirec+os.sep+"Data.QC.gz")

def gwasqualityplink(traindirec,fdirec):
 os.system("gunzip -c ./"+traindirec+os.sep+"plinkData.txt.gz | awk 'NR==1 || ($11 > 0.01) && ($10 > 0.8) {print}' | gzip  > ./"+fdirec+os.sep+"Data.gz") 
 os.system("gunzip -c ./"+fdirec+os.sep+"Data.gz | awk '{seen[$3]++; if(seen[$3]==1){ print}}' | gzip -> ./"+fdirec+os.sep+"Data.nodup.gz")
 os.system("gunzip -c ./"+fdirec+os.sep+"Data.nodup.gz | awk '!( ($4==\"A\" && $5==\"T\") || ($4==\"T\" && $5==\"A\") || ($4==\"G\" && $5==\"C\") || ($4==\"C\" && $5==\"G\")) {print}' | gzip > ./"+fdirec+os.sep+"Data.QC.gz")




 
def makephenofile():
  data = pd.read_csv(fdirec+os.sep+"test.QC.fam",sep="\s+",header=None)
  new = pd.DataFrame()
  new['FID'] = data[0].values
  new['IID'] = data[1].values
  new['phenotype'] = data[5].values
  new.to_csv(testdirec+os.sep+"YRI.pheno",sep="\t",index=False)



def testquality(testdirec,fdirec,x1,y1,z1,c1,c2,c3):
 os.system("./plink --bfile ./"+testdirec+os.sep+"test --maf 0.01 --hwe 1e-6 --geno 0.01 --mind 0.7 --write-snplist --make-just-fam --allow-no-sex -out ./"+fdirec+os.sep+"test.QC") 
 os.system("./plink --bfile ./"+testdirec+os.sep+"test --keep ./"+testdirec+os.sep+"test.fam  --allow-no-sex --indep-pairwise "+str(x1)+" "+str(y1)+" "+str(z1)+" --out ./"+fdirec+os.sep+"test.QC")
 os.system("./plink --bfile ./"+testdirec+os.sep+"test --extract ./"+fdirec+os.sep+"test.QC.prune.in --keep ./"+testdirec+os.sep+"test.fam --het --out ./"+fdirec+os.sep+"test.QC")
 os.system("Rscript QCtarget2.R "+direc+" 1")
 os.system("./plink --bfile ./"+testdirec+os.sep+"test --make-bed --allow-no-sex  --out ./"+fdirec+os.sep+"test.QC --extract ./"+fdirec+os.sep+"test.QC.snplist --exclude ./"+fdirec+os.sep+"test.mismatch --a1-allele ./"+fdirec+os.sep+"test.a1")
 os.system("Rscript QCtarget2.R  "+direc+"  3")
 os.system("./plink --bfile ./"+testdirec+os.sep+"test --clump-p1 "+str(c2)+" --clump-r2 "+str(c3)+" --clump-kb "+str(c1)+" --clump ./"+fdirec+os.sep+"Data.QC.Transformed --clump-snp-field SNP --clump-field P --out ./"+fdirec+os.sep+"test")
 os.system("awk 'NR!=1{print $3}' ./"+fdirec+os.sep+"test.clumped >  ./"+fdirec+os.sep+"test.valid.snp")
 os.system("awk '{print $3,$8}' ./"+fdirec+os.sep+"Data.QC.Transformed >  ./"+fdirec+os.sep+"SNP.pvalue")
 os.system("echo \"0.001 0 0.001\" > ./"+fdirec+os.sep+"range_list")
 os.system("echo \"0.05 0 0.05\" >> ./"+fdirec+os.sep+"range_list")
 os.system("echo \"0.1 0 0.1\" >> ./"+fdirec+os.sep+"range_list")
 os.system("echo \"0.2 0 0.2\" >> ./"+fdirec+os.sep+"range_list")
 os.system("echo \"0.3 0 0.3\" >> ./"+fdirec+os.sep+"range_list")
 os.system("echo \"0.4 0 0.4\" >> ./"+fdirec+os.sep+"range_list")
 os.system("echo \"0.5 0 0.5\" >> ./"+fdirec+os.sep+"range_list")
 os.system("./plink --bfile ./"+testdirec+os.sep+"test --score ./"+fdirec+os.sep+"Data.QC.Transformed 3 4 9 header --q-score-range ./"+fdirec+os.sep+"range_list ./"+fdirec+os.sep+"SNP.pvalue --extract ./"+fdirec+os.sep+"test.valid.snp --out ./"+fdirec+os.sep+"test")
 os.system("./plink --bfile ./"+fdirec+os.sep+"test.QC --indep-pairwise "+str(x1)+" "+str(y1)+" "+str(z1)+" --out ./"+fdirec+os.sep+"test")
 numberofpca = len(pd.read_csv(testdirec+os.sep+"test.fam",sep="\s+",header=None))

 print("Number of PCA",numberofpca)

 if numberofpca>6:
  os.system("./plink --bfile ./"+fdirec+os.sep+"test.QC --extract ./"+fdirec+os.sep+"test.prune.in --pca 6 --out ./"+fdirec+os.sep+"test")
 else:
  os.system("./plink --bfile ./"+fdirec+os.sep+"test.QC --extract ./"+fdirec+os.sep+"test.prune.in --pca "+str(numberofpca)+" --out ./"+fdirec+os.sep+"test")
  modifypca(fdirec,numberofpca)

 makephenofile()

def modifypca(fdirec,numberofpca):
 data = pd.read_csv(fdirec+os.sep+"test.eigenvec",header=None,sep="\s+")
 del data[numberofpca+1] 
 data.to_csv(fdirec+os.sep+"test.eigenvec",index=False,header=False,sep="\t")




def changechromosomeforlasso():
 # clean Lasso result file.
 try:
  os.remove(direc+os.sep+"result"+os.sep+"test.txt")
 except:
  pass
 data = pd.read_csv(direc+os.sep+"files"+os.sep+"Data.QC.gz",compression='gzip',sep="\t")
 data =data[data['CHR'].isin(list(range(1,23)))]
 print(data.head())
 data.to_csv(direc+os.sep+"files"+os.sep+"Data.QC.gz",index=False, sep="\t",compression='gzip')


def ldpred(direc):
  # change GWAS headers for ldpred-2
  data = pd.read_csv(direc+os.sep+"files"+os.sep+"Data.QC.gz",compression='gzip',sep="\t")
  data = data.rename(columns={'CHR': 'chr', 'SNP': 'rsid', 'BP': 'pos', 'A1': 'a0', 'A2': 'a1'})
  data =data[data['chr'].isin(list(range(1,22)))]
  #print(data.head())
  
  data.to_csv(direc+os.sep+"files"+os.sep+"Data.QC.gz",index=False, sep="\t",compression='gzip')
  os.system("./plink  -bfile  ./"+fdirec+os.sep+"test.QC  --chr 1-22 --make-bed --out ./"+fdirec+os.sep+"test.QC2")
  try:
    os.system("rm -rf "+direc+os.sep+"files/test.QC2.bk")
    os.remove(direc+os.sep+"files/test.QC2.bk")
    
  except:
    pass
  os.system("Rscript QCtarget.R  "+direc+" 5")

 

# Define directories and paths
direc = sys.argv[1].replace('\r', '')+ os.sep + sys.argv[2].replace('\r', '')

traindirec = direc+os.sep+"train"+os.sep
testdirec =direc+os.sep+"test"+os.sep
result = direc+os.sep+"result"+os.sep
fdirec = direc+os.sep+"files"+os.sep

if not os.path.isdir(direc+os.sep+"files"):
  os.mkdir(direc+os.sep+"files")
if not os.path.isdir(direc+os.sep+"result"):
  os.mkdir(direc+os.sep+"result")  
 
# Pruning parameters
x = [200,500,1000]
y = [50,100,150]
z = [0.1,0.3,0.5]

# Clumping parameters
clump1 = list(range(200,1001,200))
clump2 = [1]
clump3 = np.linspace(0.1,1,5)




# Variables to store the final AUC for all tools.
plinkAUC = []
plinkMCC = []

lassoAUC = []
lassoMCC = []

prsiceAUC = []
prsiceMCC = []

ldpredAUC1 = []
ldpredAUC2 = []


def lasso(direc):
 os.system("Rscript QCtarget2.R  "+direc+" 6")
 lassoAUC.append(calculateAUClasso(direc))
 lassoMCC.append(calculateMCClasso(direc))


def PRSice(fdirec,testdirec):
  #os.system("Rscript PRSice.R --prsice PRSice --base "+fdirec+"/Data.QC.gz --target "+fdirec+"/test.QC --thread 1 --print-snp --stat OR --binary-target T --interval 5e-01 --upper 1.0 --clump-r2 "+str(c3)+" --clump-kb "+ str(c1) +" --out "+result+os.sep+"PRScice_PRS")
  #print("Rscript PRSice.R --prsice PRSice --base "+fdirec+"/Data.QC.gz --target "+fdirec+"/test.QC --thread 1 --print-snp --stat OR --binary-target T --interval 5e-01 --upper 1.0  --out "+result+os.sep+"PRScice_PRS")
  #os.system("Rscript PRSice.R --prsice PRSice --base "+fdirec+"/Data.QC.gz --target "+fdirec+"/test.QC  --thread 1 --print-snp --stat OR --binary-target T --interval 5e-01 --upper 1.0  --out "+result+os.sep+"PRScice_PRS")
  #os.system("Rscript PRSice.R --prsice PRSice --base "+traindirec+"/train.assoc.fisher --target "+fdirec+"/test.QC --thread 1 --print-snp --stat OR --binary-target T --interval 5e-01 --upper 1.0  --out "+result+os.sep+"PRScice_PRS")

  #prsiceAUC.append(calculateAUCPRScice(direc,testdirec))
  #prsiceMCC.append(calculateMCCPRScice(direc,testdirec))
  #print(prsiceAUC)
  pass

def plink(direc):
 os.system("Rscript QCtarget2.R  "+direc+" 4")
 #plinkAUC.append(calculateAUCPlink(direc))
 #plinkMCC.append(calculateMCCPlink(direc))

 #print(plinkAUC)

def gwasquality(traindirec,fdirec):
 os.system("gunzip -c ./"+traindirec+os.sep+"Data.txt.gz | awk 'NR==1 || ($11 > 0.01) && ($10 > 0.8) {print}' | gzip  > ./"+fdirec+os.sep+"Data.gz") 
 os.system("gunzip -c ./"+fdirec+os.sep+"Data.gz | awk '{seen[$3]++; if(seen[$3]==1){ print}}' | gzip -> ./"+fdirec+os.sep+"Data.nodup.gz")
 os.system("gunzip -c ./"+fdirec+os.sep+"Data.nodup.gz | awk '!( ($4==\"A\" && $5==\"T\") || ($4==\"T\" && $5==\"A\") || ($4==\"G\" && $5==\"C\") || ($4==\"C\" && $5==\"G\")) {print}' | gzip > ./"+fdirec+os.sep+"Data.QC.gz")

def gwasqualityplink(traindirec,fdirec):
 os.system("gunzip -c ./"+traindirec+os.sep+"plinkData.txt.gz | awk 'NR==1 || ($11 > 0.01) && ($10 > 0.8) {print}' | gzip  > ./"+fdirec+os.sep+"Data.gz") 
 os.system("gunzip -c ./"+fdirec+os.sep+"Data.gz | awk '{seen[$3]++; if(seen[$3]==1){ print}}' | gzip -> ./"+fdirec+os.sep+"Data.nodup.gz")
 os.system("gunzip -c ./"+fdirec+os.sep+"Data.nodup.gz | awk '!( ($4==\"A\" && $5==\"T\") || ($4==\"T\" && $5==\"A\") || ($4==\"G\" && $5==\"C\") || ($4==\"C\" && $5==\"G\")) {print}' | gzip > ./"+fdirec+os.sep+"Data.QC.gz")

def calculateAUCldpred1(path,test):
  best = pd.DataFrame()
  best["PRS"] = pd.read_csv(path+os.sep+"result"+os.sep+"autoldpredtest.txt",sep="\s+")["x"].values
  pheno = pd.read_csv(test+os.sep+"test.fam",sep="\s+",header=None,names=["a","b","c","d","e","f"])
  pheno.f[pheno['f']==1]=0
  pheno.f[pheno['f']==2]=1
  temp = NormalizeData(best["PRS"].values)
  best["PRS"] = NormalizeData(np.array(best["PRS"].values))
  best["PRS"].values[best["PRS"] >=0.5] = 1
  best["PRS"].values[best["PRS"] < 0.5] = 0

  #print("AUC",roc_auc_score(np.array(best["PRS"].values), np.array(pheno['f'].values)))

  ldpredAUC1.append(roc_auc_score(best["PRS"].values, pheno['f'].values))


def calculateAUCldpred3(path,test):
  best = pd.DataFrame()
  best["PRS"] = pd.read_csv(path+os.sep+"result"+os.sep+"infldpredtest.txt",sep="\s+")["x"].values
  pheno = pd.read_csv(test+os.sep+"test.fam",sep="\s+",header=None,names=["a","b","c","d","e","f"])
  pheno.f[pheno['f']==1]=0
  pheno.f[pheno['f']==2]=1
  temp = NormalizeData(best["PRS"].values)
  best["PRS"] = NormalizeData(np.array(best["PRS"].values))
  best["PRS"].values[best["PRS"] >=0.5] = 1
  best["PRS"].values[best["PRS"] < 0.5] = 0

  ldpredAUC2.append(roc_auc_score(best["PRS"].values, pheno['f'].values))
import sys,os

for x1 in x:
 for x2 in y:
  for x3 in z:
   for c1 in clump1:
    for c2 in clump2:
     for c3 in clump3:

      rs = pd.DataFrame()
      rs2 = pd.DataFrame()
      gwasqualityplink(traindirec,fdirec) 
      testquality(testdirec,fdirec,x1,x2,x3,c1,c2,c3)
      plink(direc)

      changechromosomeforlasso()

      try:
       os.remove("Rplots.pdf") 
      except:
       pass
      lasso(direc)
      rs['LassosumAUC'] = lassoAUC
      rs2['LassosumMCC'] = lassoMCC
      
      print("Working!")

      rs.to_csv(direc+os.sep+"Results_Lassosum_AUC.txt",index=False)
      rs2.to_csv(direc+os.sep+"Results_Lassosum_MCC.txt",index=False)
      print(rs)
      print(rs2)



rs = pd.DataFrame()
rs['LassosumAUC'] = lassoAUC
rs2 = pd.DataFrame()
rs2['LassosumMCC'] = lassoMCC


rs.to_csv(direc+os.sep+"Results_Lassosum_AUC.txt",index=False)
rs2.to_csv(direc+os.sep+"Results_Lassosum_MCC.txt",index=False)

