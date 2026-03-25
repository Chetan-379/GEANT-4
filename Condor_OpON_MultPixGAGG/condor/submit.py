import os,sys,re,fileinput,string,shutil
from datetime import date
min_count=0
max_count=100
count=0;
case=["Correlated", "Uncorrelated"]

for i in range(max_count):
        for j in case:                
                condorSubmit = "condor_submit Submit_condor/submitCondor_OpON_%s_%s"%(j,i)
                print(condorSubmit)
                os.system(condorSubmit)
