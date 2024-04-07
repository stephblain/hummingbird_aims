#!/bin/bash

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=run_hmm_f50
#SBATCH --time=24:00:00
#SBATCH --ntasks=6
#SBATCH --mem=360G
#SBATCH --output=error_hmm_f50.%j

source activate AncestryEnv
source activate py2 #this one after so bamutils runs
module load GCC/10.3.0 OpenMPI/4.1.1 Armadillo/10.7.5 GSL/2.7

ancestry_hmm -a 2 0.5 0.5 -p 0 -1000000 0.5 -p 1 -3000 0.5 -s ancestryHMM_samples_HummerHybrids_f50 -i ancestryHMM_counts_HummerHybrids_f50_autosomes