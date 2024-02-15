#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import sys


projectName=sys.argv[1] #"HummerHybrids"
f_tag=sys.argv[2] #"f50"
aimsFile=sys.argv[3] #"HummerRefPanel.f50.noUnderscores.aims"
aims_cmsFile=sys.argv[4] #"HummerRefPanel.f50.noUnderscores.aims.cMs"
parhmmFile=sys.argv[5] #"HummerRefPanel.f50.parhmm"
hybridCountsFile=sys.argv[6] #"../HummerHybrids.f50.hmmFormat"
indsFile=sys.argv[7] #"../HummerHybrids.f50.hmmFormat.ind2"



############################################################
#read in data files and rename columns
############################################################

aims=pd.read_csv(aimsFile,header=None,sep="\t")
aims_cms=pd.read_csv(aims_cmsFile,header=None,sep="\t")
parhmm=pd.read_csv(parhmmFile,header=None,sep="\t")
inds=pd.read_csv(indsFile,header=None,sep="\t")
counts=pd.read_csv(hybridCountsFile,header=None,sep="\t")


aims.columns=["chr","pos","Ref","Alt"]
aims_cms.columns=["chr","pos","Ref","Alt","cMs"]
parhmm.columns=["chr","pos","sp1ref","sp1alt","sp2ref","sp2alt"]

hybCols=['chr','pos','Ref','Alt'] #first two cols are chromosome and position
#duplicate individual names (ind1, ind2 > ind1, ind1, ind2, ind2)
hybCols.extend([val for val in list(inds[0]) for _ in (0, 1)])
counts.columns=hybCols
#counts.head(10)


############################################################
#merge all dataframes, drop ref and alt nucleotide columns
############################################################

aims=aims.merge(parhmm,how="left").merge(aims_cms,how="left").merge(counts,how="left")

aims=aims.drop(["Ref","Alt"],axis=1)
aims.head(10)


#check that the number of rows in the output make sense
print("Number of rows in output")
print(aims.shape[0])

print("Number of rows in cMs input")
print(aims_cms.shape[0])

#aims['chr'].unique()


inds['ploidy']="2" #declare all individuals to be diploid


############################################################
#write output files
############################################################


outFile="ancestryHMM_counts_"+projectName+"_"+f_tag
aims.to_csv(outFile,index=False,header=False,sep="\t",float_format='%.15f')
print("counts output printed to: "+outFile)

outFile2="ancestryHMM_samples_"+projectName+"_"+f_tag
inds.to_csv(outFile2,index=False,header=False,sep="\t")
print("samples output printed to: "+outFile2)

