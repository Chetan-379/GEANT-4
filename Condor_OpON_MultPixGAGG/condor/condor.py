import os,sys,re,fileinput,string,shutil
from datetime import date
count=0
min_count=0
max_count = 100
case=["Correlated", "Uncorrelated"]
                                                                     
count1=0;
for i in range (max_count):
        for j in case:
                condorSubmit = "Submit_condor/submitCondor_OpON_%s_%s"%(j,i)        
                fname = "/afs/cern.ch/user/c/cagrawal/public/run_G4/GEANT-4/OpON_MultPixGAGG/condor/OpON_%s_100kEvts_%s.root"%(j,i)
                executable = "test_OpON_%s"%(j)
                shutil.copyfile("proto_condor_submit",condorSubmit)
                for line in fileinput.FileInput(condorSubmit, inplace=1):
                        line=line.replace("JobID", str(i))
                        line=line.replace("outfile", fname)
                        line=line.replace("execute", executable)
                        line=line.replace("case", j)
                        print(line.rstrip())
                        
                print(condorSubmit)        
count1+=1
            
