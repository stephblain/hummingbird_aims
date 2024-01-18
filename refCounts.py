# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 13:16:14 2024

@author: Steph
"""

#get a list of ancestry informative markers (AIMs) from total and alt allele counts
#for two species

#four command line arguments: (1) minimum allele frequency difference
#(2) project ID, (3 & 4) species 1 & 2 allele count files locations


import pandas as pd
import os
# import sys

freq="0.5"
project="HummerRefPanel"
sp1File="C:/Users/Steph/GitHub_data/hummers/counts.test.ruby"
sp2File="C:/Users/Steph/GitHub_data/hummers/counts.test.black"

os.chdir("C:/Users/Steph/GitHub_data/hummers")



# freq=sys.argv[1]
# project=sys.argv[2]
# sp1File=sys.argv[3]
# sp2File=sys.argv[4]


print(f"minimum frequency difference : {freq}")
print(f"species 1 count file : {sp1File}")
print(f"species 2 count file : {sp2File}")

#read in count files and rename columns
sp1Counts=pd.read_csv(sp1File,header=None,sep=" ")
sp1Counts.columns=["chr","pos","ref","alt","All","Alt"]
sp2Counts=pd.read_csv(sp2File,header=None,sep=" ")
sp2Counts.columns=["chr","pos","ref","alt","All","Alt"]

#get and save absolute value of the difference between species counts
spDiff=abs((sp1Counts["Alt"]/sp1Counts["All"])-(sp2Counts["Alt"]/sp2Counts["All"]))

#select chromosomes and positions
aims=sp1Counts[spDiff>float(freq)][["chr","pos","ref","alt"]]

print(f"number of AIMs above cutoff : {len(aims)}")

freq2=str(int(float(freq)*100)) #format frequency for output

aims.to_csv(project+".f"+freq2+".aims",header=None,sep="\t",
            index=False)

#make a positions file for selecting aims from vcfs - original scaffold names
aims[["chr","pos"]].to_csv(project+".f"+freq2+".aims.pos",header=None,sep="\t",
            index=False)

#make a version of aims file with no underscores
aims2=aims.copy()
aims2["chr"]=aims2["chr"].str.replace("_","-",regex=False).str.replace(".","-",regex=False)

aims2.to_csv(project+".f"+freq2+".noUnderscores.aims",header=None,sep="\t",
            index=False)

#and make the mod file for it
aims2["chr"]=aims2["chr"]+"_"+aims2["pos"].astype(str)
aims2.to_csv(project+".f"+freq2+".noUnderscores.aims.mod",header=None,sep="\t",
            index=False)



#join chr and pos columns for matching among dataframes
aims["chrpos"]=aims["chr"]+"_"+aims["pos"].astype(str)
sp1Counts["chrpos"]=sp1Counts["chr"]+"_"+sp1Counts["pos"].astype(str)
sp2Counts["chrpos"]=sp2Counts["chr"]+"_"+sp2Counts["pos"].astype(str)

#only keep counts at aims
sp1aims=sp1Counts.loc[sp1Counts["chrpos"].isin(aims["chrpos"])]
sp2aims=sp2Counts.loc[sp2Counts["chrpos"].isin(aims["chrpos"])]

#estimate ref count from total and alt counts
par=aims.assign(Refsp1=(sp1aims["All"]-sp1aims["Alt"]),
                Altsp1=sp1aims["Alt"],
                Refsp2=(sp2aims["All"]-sp2aims["Alt"]),
                Altsp2=sp2aims["Alt"])

#select relevant columns in correct order
par=par[["chr","pos","Refsp1","Altsp1","Refsp2","Altsp2"]]

par["chr"]=par["chr"].str.replace("_","-",regex=False).str.replace(".","-",regex=False)

par.to_csv(project+".f"+freq2+".parhmm",header=None,sep="\t",
            index=False)