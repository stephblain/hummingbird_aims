#!/usr/bin/env python3


############
#estimates the distance in cM between bp based on a recomb map
#inputs:
	##"aimsFile" = a tab-delimited file with chromosome, position, reference allele and alternative allele
	##"recomb_M_bp" = average recombination rate in Morgans per base pair, 0.0000000369 M/bp in flycatchers (Kawakami et al. 2017)
############

print("start py script")

import pandas as pd
import numpy as np
import sys

aimsFile=sys.argv[1]
recomb_M_bp_str=sys.argv[2]

recomb_M_bp=float(recomb_M_bp_str)

aims=pd.read_csv(aimsFile,header=None,sep="\t")
aims.columns=["chr","pos","Ref","Alt"]

print("recomb rate Morgans per bp:",recomb_M_bp)

#all scaffolds, including mtDNA
print("outputting scaffolds:",aims.chr.unique())


for chr1 in aims.chr.unique():

	aimsChr=aims[aims['chr']==chr1]
	aimsChr=aimsChr.reset_index(drop=True)
	print("start scaffold:",chr1)

	for i in range(1,len(aimsChr)): #start at 1 to skip first position

		dist_bp=aimsChr.iloc[i]["pos"]-aimsChr.iloc[i-1]["pos"] #get distance in bp from previous site
		dist_cM=dist_bp*recomb_M_bp*100 #get distance in cM
		aimsChr.loc[i,"cM"]=dist_cM# #add value to dataframe

    aimsChr.loc[0,"cM"]=0

	outFile=aimsFile+"_"+chr1+".cMs"
	aimsChr.to_csv(outFile,index=False,header=False,sep="\t",float_format='%.15f')

print("outputted scaffolds:",aims.chr.unique())