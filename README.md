# Hummingbird AIMs Pipeline

Pipeline for identifying and calling ancestry informative markers between blackchinned and rubythroated hummingbirds and their hybrids

All sequences were aligned to the Anna's hummingbird reference genome

## Main pipeline scripts

### 01_refPanel_Snakefile.py

Processes reference panel from fastq's to output divergent markers between the species and the genetic distances between those markers

When running in the cluster, put this in it's own folder and rename it "Snakefile"

Run with snakemake - see file header for modules to load in TAMU cluster

### 02_hybridCounts_Snakefile.py

Processes hummingbird hybrid sequences from bams to read counts, then calls script to combine hybrid read counts with reference panel markers

When running in the cluster, put this in it's own folder and rename it "Snakefile"

Run with snakemake - see file header for modules to load in TAMU cluster

### 03_run_hmm.sh

Run AncestryHMM

This is a slurm job submission script formatted for the TAMU cluster

### 04_combine_posteriors.sh

AncestryHMM outputs a separate file for each individual, this script processes those outputs to make:

(1) a combined file with posterior probabilities

(2) a combined file with called states for sites where probability > 0.9

## Scripts called by others in the pipeline

### aims_bp_to_cms.py

Estimate genetic distance in centimorgans between ancestry informative markers, based on physical distance and recombination rate

Called by 01_refPanel_Snakefile.py

### get_posterior_states.R

Call states from posterior probabilities outputted by AncestryHMM

Called by 04_combine_posteriors.sh

### get_refPanel_genomes.sh

Access publicly available reference genome and reference panel high coverage sequences from SRA and NCBI

Do this before running 01_refPanel_Snakefile.py

### HummerHybrids_check.R

Just a cute little script to run on the states file output from 04_combine_posteriors.sh to make sure the whole process didn't got terribly horribly wrong

### merge_dfs.py

Python script to merge various inputs (reference panel allele counts, genetic distances in cMs, hybrid read counts) for AncestryHMM

Called by 02_hybridCounts_Snakefile.py

### refCounts.py

Count alleles in the reference panel to determine divergent markers

Called by 01_refPanel_Snakefile.py

### pseudoref_Snakefile

Another Snakefile - this one to make a pseudoref from fastq'script

Did not end up being part of this pipeline, but could be part of a future pipeline pending new genomic resources