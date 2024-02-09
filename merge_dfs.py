#!/usr/bin/env python
# coding: utf-8

#Combine input data for AncestryHMM counts file

import pandas as pd
import sys

#collect all command line arguments
projectName=sys.argv[1] #"HummerHybrids"
f_tag=sys.argv[2] #"f50"
aimsFile=sys.argv[3] #"HummerRefPanel.f50.noUnderscores.aims"
aims_cmsFile=sys.argv[4] #"HummerRefPanel.f50.noUnderscores.aims.cMs"
parhmmFile=sys.argv[5] #"HummerRefPanel.f50.parhmm"
hybridCountsFile=sys.argv[6] #"../HummerHybrids.f50.readCounts"
indListFile=sys.argv[7] #"../HummerHybrids.f50.readCounts.ind"

#read in files
aims=pd.read_csv(aimsFile,header=None,sep="\t")
aims_cms=pd.read_csv(aims_cmsFile,header=None,sep="\t")
parhmm=pd.read_csv(parhmmFile,header=None,sep="\t")
hybCounts=pd.read_csv(hybridCountsFile,header=None,sep="\t")
indList=pd.read_csv(indListFile,header=None,sep="\t")

#rename columns; all df's need a "chr" and "pos" column
aims.columns=["chr","pos","Ref","Alt"]
aims_cms.columns=["chr","pos","Ref","Alt","cMs"]
parhmm.columns=["chr","pos","sp1ref","sp1alt","sp2ref","sp2alt"]
hybCols=['chr','pos'] #first two cols are chromosome and position
#duplicate individual names (ind1, ind2 > ind1, ind1, ind2, ind2)
hybCols.extend([val for val in list(indList[0]) for _ in (0, 1)])
hybCounts.columns=hybCols

#left join so that everything aligns with aims file
aims=aims.merge(parhmm,how="left").merge(aims_cms,how="left").merge(hybCounts,how="left")

aims=aims.drop(["Ref","Alt"],axis=1)

#print("first few lines of aims + parents + cms")
#aims.head(10)

print("NAs per column")
aimsNAs=aims.isna().sum()
print(aimsNAs)

indList['ploidy']="2" #declare all individuals to be diploid

outFile="ancestryHMM_counts_"+projectName+"_"+f_tag
aims.to_csv(outFile,index=False,header=False,sep="\t",float_format='%.15f')
print("counts output printed to: "+outFile)

outFile2="ancestryHMM_samples_"+projectName+"_"+f_tag
indList.to_csv(outFile2,index=False,header=False,sep="\t")
print("samples output printed to: "+outFile2)

