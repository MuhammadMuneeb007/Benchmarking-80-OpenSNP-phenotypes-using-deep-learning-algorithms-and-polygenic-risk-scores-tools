import pandas as pd
import numpy as np
import os
import shutil
import re 
import glob

def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def _23andmetoBED_BIM_FAM(input_file):
    files = []
    allfiles = os.listdir("./"+input_file+"/")
    print(allfiles)
    personname = [] 
    sexinfo = "0"
    for loop in allfiles:
      if ".23andme" in loop:
        data = loop.split(".")[0].split("_")
        personname.append(data[0])
        print(data)
        if data[5] =="XX":
          sexinfo = "2"
        elif data[5] =="XY":
          sexinfo = "1"
        else:
          sexinfo = "0" 
        
        os.system("./plink --23file ./"+input_file+os.sep+loop+" --snps-only --make-bed  --out ./"+input_file+os.sep+data[0])
        if os.path.exists("./"+input_file+os.sep+data[0]+".fam"):
          data2 = pd.read_csv("./"+input_file+os.sep+data[0]+".fam",header=None, sep="\s+")
          data2[0] = data[0]
          data2[1] = data[0]
          data2[4] = sexinfo    
          data2.to_csv("./"+input_file+os.sep+data[0]+".fam",sep="\t",header=False,index=False)
    
    allfiles =  os.listdir("./"+input_file+"/")
    count=0
    files=[]
    for loop in allfiles:
        if ".txt" in loop and ".bed" not in loop and ".fam" not in loop and ".bim" not in loop:
          print(loop)
          x = loop.split("_")[0]
          x = x + ".fam"
          me = os.path.exists("./"+input_file+os.sep+x)
          if me==True:
            x = x.split(".")[0]
            x = "./"+input_file+os.sep+x +".bed " + "./"+input_file+os.sep+x + ".bim " +  "./"+input_file+os.sep+x + ".fam"	
            files.append(x)
          else:
            count=count+1
          print(count," People removed due to missing fam file")
    with open("./"+input_file+os.sep+"All.txt", "w") as filehandle:
      for listitem in files:
        filehandle.write('%s\n' % listitem)
    
    os.system("./plink --merge-list ./"+input_file+os.sep+"/All.txt --make-bed --out 23andmetoBED")
    
    if os.path.exists("23andmetoBED.bed"):
      exit(0)
    else:
      allfiles =  os.listdir("./"+input_file+"/")
      count=0
      for loop in allfiles:
        if ".bed" in loop:
          x = loop
          x = x.split(".")[0]
    
          command = "./plink --bfile ./"+input_file+os.sep+x+" --exclude 23andmetoBED-merge.missnp --make-bed --out ./"+input_file+os.sep + x
          os.system(command)
      allfiles =  os.listdir("./"+input_file+"/")
      files=[]
      
      for loop in allfiles:
          if ".txt" in loop and ".bed" not in loop and ".fam" not in loop and ".bim" not in loop:
            print(loop)
            x = loop.split("_")[0]
            x = x + ".fam"
            me = os.path.exists("./"+input_file+os.sep+x)
            if me==True:
              x = x.split(".")[0]
              x = "./"+input_file+os.sep+x +".bed " + "./"+input_file+os.sep+x + ".bim " +  "./"+input_file+os.sep+x + ".fam"	
              files.append(x)
            else:
              count=count+1
            print(count," People removed due to missing fam file")
      with open("./"+input_file+os.sep+"All.txt", "w") as filehandle:
        for listitem in files:
          filehandle.write('%s\n' % listitem)      
      os.system("./plink --merge-list ./"+input_file+os.sep+"/All.txt --make-bed --out ./"+input_file.split("/")[0]+"/23andmetoBED")    

def ancestryto23andme(source,filename,dest):
 data = pd.read_csv(source+os.sep+filename,sep="\t",comment='#')
 new = pd.DataFrame()
 new['Rsid'] = data['rsid'].values
 new['Chromosome'] = data['chromosome'].values
 new['position'] = data['position'].values
 new['genotype'] = data['allele2']+ data['allele1']
 new['Chromosome'] = new['Chromosome'].replace(23, 'X')
 new['Chromosome'] = new['Chromosome'].replace(24, 'Y')
 new['Chromosome'] = new['Chromosome'].replace(25, 'XY')
 new['Chromosome'] = new['Chromosome'].replace(26, 'MT')
 filename = filename.replace("ancestry","23andme")
 new.to_csv(loop+os.sep+"Data"+os.sep+filename, sep="\t",index=False,header=False)
 
originalpath = "./opensnp_datadump.current"
originalfiles = os.listdir(originalpath)
allfiles = os.listdir("./")
allfiles = sorted_nicely(allfiles)
actualfile = [] 
filename = [] 
actualpheno = []
data2 = pd.DataFrame()
x = 0

import sys
allfiles=[str(sys.argv[1])]
allfiles = pd.read_csv("phenotypes.txt",header=None)[0].values


for loop in allfiles:

 originalpath = "./opensnp_datadump.current"
 originalfiles = os.listdir(originalpath)
 allfiles = os.listdir("./")
 allfiles = sorted_nicely(allfiles)
 actualfile = [] 
 filename = [] 
 actualpheno = []
 data2 = pd.DataFrame()
 x = 0

 if ".txt" not in loop and "23andme" not in loop and ".csv" not in loop and ".py" not in loop and ".html" not in loop and "opensnp_datadump.current" not in loop:
  data = pd.read_csv(loop+os.sep+loop+".csv")
  tempfiles = data['Name'].values
  # Files in each phenotype
  for index, row in data.iterrows():
   strings_with_substring = [string for string in originalfiles if str(row['Name']) in string]
   for subfiles in strings_with_substring:
    if "ftdna-illumina" not in subfiles:
     tempuser = subfiles.split("_")
     tempuser = tempuser[0].replace("user","")
     tempfile = subfiles.split("_")
     tempfile = tempfile[1].replace("file","")
     try:
      if str(tempuser)==str(row['Name']) and str(tempfile)==str(row['File'].split(".")[-1]):
       if not os.path.isdir(loop+os.sep+"Data"): 
        os.mkdir(loop+os.sep+"Data")
       if "ancestry" in subfiles:
        ancestryto23andme(originalpath,subfiles,loop+os.sep+"Data"+os.sep)
       else:
        shutil.copy(originalpath+os.sep+subfiles, loop+os.sep+"Data"+os.sep)
       actualfile.append(tempuser)
       filename.append(subfiles)
       actualpheno.append(row['Phenotype'])
     except:
      pass
  try:
   print("Total number of people",data.shape[0])
   data2['Username'] =  actualfile
   data2['File name'] = filename
   data2['Phenotype'] = actualpheno
   data2.to_csv(loop+os.sep+"Phenotype_process1.csv")
   _23andmetoBED_BIM_FAM(loop+os.sep+"Data")
   shutil.rmtree(loop+os.sep+"Data/")
  except:
   pass
  os.system("./plink --bfile ./"+loop+os.sep+"23andmetoBED --maf 0.01 --hwe 1e-6 --geno 0.01 --mind 0.7 --allow-no-sex --make-bed -out ./"+loop+os.sep+"final.QC")
  #os.remove(loop+os.sep+"23andmetoBED.*")
  phenotype = pd.read_csv(loop+os.sep+loop+".csv")
  phenotype.sort_values('Name',inplace=True)
  values =  phenotype['Phenotype'].unique()

  phenotype['Phenotype'] = phenotype['Phenotype'].replace(values[0],1)
  phenotype['Phenotype'] = phenotype['Phenotype'].replace(values[1],2)
  famfile = pd.read_csv(loop+os.sep+"final.QC.fam",header=None,sep="\s+")
  famfile[0] = famfile[0].str.replace('user','').astype(int)
  phenotype = phenotype.loc[phenotype['Name'].isin(list(famfile[0].values))]
  phenotype = phenotype.drop_duplicates(subset=['Name'], keep='first')

  famfile[0] = famfile[1].values
  famfile[5] = phenotype['Phenotype'].values

  famfile.to_csv(loop+os.sep+"final.QC.fam",header=False,index=False,sep="\t")
  phenotype.to_csv(loop+os.sep+"Phenotype_process2.csv",index=False,sep="\t")
 


