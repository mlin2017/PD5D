#!/usr/bin/env python3

import pandas as pd
import os
lsdir = os.listdir("/home/jap0606/PD5D_Repo/PD5D/scRNA_and_scATAC-seq_Markdowns")
lsdir.sort()
stepper=1
for i in lsdir:
    if re.search("batch[0-9]+_ASAP_snRNA-seq_[0-9]+",i):
        temppath="/n/scratch3/users/j/jap0606/FullIntegration/"+i+"/matrix.mtx"
        tempframe=pd.read_csv(temppath,sep=' ',skiprows=3)
        unicells = list(dict.fromkeys(tempframe.iloc[:2]))
        for cell in unicells.iloc:
            cellframe = tempframe[tempframe.iloc[:2] == cell]
            
            
 = dataframe.loc[dataframe['Percentage'] > 80]
        if stepper == 1:
            finalframe=tempfile
        else:
            finalframe = pd.concat([finalframe, tempfile])     

finalframe.to_csv('/n/scratch3/users/j/jap0606/FullIntegration/matrix.mtx2', sep=' ', index=False, header=True, mode='a')
