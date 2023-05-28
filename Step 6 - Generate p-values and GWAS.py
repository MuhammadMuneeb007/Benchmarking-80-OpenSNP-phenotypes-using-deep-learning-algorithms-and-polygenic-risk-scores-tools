from email.header import Header
import pandas as pd
import numpy as np
import os
import shutil
import re 
import glob
import sys
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedShuffleSplit


phenotypename = sys.argv[1]
iterations = ["1","2","3","4","5"]


# Train/test split.
famfile = pd.read_csv(phenotypename+os.sep+"final.QC.fam",header=None,sep="\s+",names=["a","b","c","d","e","f"])
sss = StratifiedShuffleSplit(n_splits=5, test_size=0.25)
X = famfile['a'].values
y = famfile['f'].values
sss.get_n_splits(X,y)
funk=0
for train_index, test_index in sss.split(famfile,y):
    iteration = iterations[funk]
    funk = funk+1
    workingdirec = phenotypename+ os.sep+ iteration
    testdirec = phenotypename +os.sep+ iteration +os.sep+ "test"
    traindirec = phenotypename +os.sep+ iteration +os.sep+ "train"

    if not os.path.isdir(workingdirec): 
        os.mkdir(workingdirec)
    
    if not os.path.isdir(testdirec): 
        os.mkdir(testdirec)
    if not os.path.isdir(traindirec): 
        os.mkdir(traindirec)
    if not os.path.isdir(phenotypename+os.sep+iteration): 
        os.mkdir(phenotypename+os.sep+iteration)

    print("TRAIN:", train_index, "TEST:", test_index)
    usertrain, usertest = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index] 

    usertrain = famfile[famfile['a'].isin(usertrain)]
    usertest = famfile[famfile['a'].isin(usertest)]


    print("Train",list(y_train).count(2))
    print("Train",list(y_train).count(1))

    print("Test",list(y_test).count(2))
    print("Test",list(y_test).count(1))



    usertest.to_csv(workingdirec+os.sep+"testname.txt",index=False,header=False,sep="\t")
    usertrain.to_csv(workingdirec+os.sep+"trainname.txt",index=False,header=False,sep="\t")


    genofilename = phenotypename+os.sep+"final.QC"

    os.system("./plink --bfile "+genofilename+" --keep "+workingdirec+os.sep+"testname.txt --make-bed --out ./"+testdirec+os.sep+"/test")
    os.system("./plink --bfile "+genofilename+" --keep "+workingdirec+os.sep+"trainname.txt --make-bed --out ./"+traindirec+os.sep+"train")
    os.system("./plink --bfile ./"+traindirec+os.sep+"train --allow-no-sex --fisher --out ./"+traindirec+os.sep+"train")

    gwasfile = pd.read_csv(traindirec+os.sep+"train.assoc.fisher",sep="\s+")
    print(gwasfile.head())

    numberofsnps = [50,100,200,500,1000,5000,10000,5]

    pvalues = []
    x=0
    for me in range(0,6): 
        gwasfile['P']=pd.to_numeric(gwasfile['P'],errors='coerce')
        l = np.linspace(gwasfile['P'].min(),gwasfile['P'].max(),100000)
        for loop in l:
            if len(pvalues)==7:
                x=7
                break
            subpvalues = gwasfile[gwasfile['P']<=float(loop)]
            count = len(subpvalues)
            if count>numberofsnps[x]:
                pvalues.append(loop)
                print(pvalues,x)
                x=x+1
                if x==6:
                    break

    for loop in pvalues:
        pdirectory = workingdirec+os.sep+"pv_"+str(loop)
        if not os.path.isdir(pdirectory): 
            os.mkdir(pdirectory)
        subpvalues = gwasfile[gwasfile['P']<=float(loop)]
        subpvalues['SNP'].to_csv(pdirectory+os.sep+str(loop)+".txt",index=False,header=False)
        os.system("./plink --bfile ./"+workingdirec+os.sep+"test/test   --extract ./"+pdirectory+os.sep+str(loop)+".txt --recodeA --out "+pdirectory+os.sep+"ptest")
        os.system("./plink --bfile ./"+workingdirec+os.sep+"train/train   --extract ./"+pdirectory+os.sep+str(loop)+".txt --recodeA --out "+pdirectory+os.sep+"ptrain")
    




print("DONE")



