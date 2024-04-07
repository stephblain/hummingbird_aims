#Snakemake pipeline to get read counts at previously defined ancestry informative markers from a set of bam files
#bams should be in a folder called "hybrid_bams" and be named in the format SAMPLE.combo.bam, with a corresponding .bai
#make another folder called "hybrid_counts" for intermediate files to live in

#Before running, load modules for various programs (Bcftools, Vcftools, snakemake)
#Ordered by dependencies for TAMU Grace cluster:
##module load GCC/11.2.0 BCFtools/1.14 VCFtools/0.1.16 BEDTools/2.30.0 OpenMPI/4.1.1 snakemake/6.10.0 Python/3.9.6 #pandas installed in python with pip

#test as: snakemake -np
#run as:  nohup snakemake --profile . &
#before running, make config.yaml file in directory

#config.yaml example:
##cluster: "sbatch --time=24:00:00 --mem-per-cpu=80G
##						--nodes=1 --ntasks=1 --output=errors/err_{rule}_{wildcards}.%j"
##jobs: 18


import pandas as pd

samples_df = pd.read_table('hummer_samples.tsv').set_index("sample", drop=False)
SAMPLES= list(samples_df['sample'])


#path to aims list, make sure there are no underscores in the scaffold names, for Hummers I also removed periods because I'm paranoid
#4 tab-delimited columns: chromosome, position, ref allele, alt allele
#in same folder, have aims.mod file with tab-delimited: chromosome_position, position, ref allele, alt allele
# and aims.mod.bed file with tab-delimited: chromosome, position, position, ref allele, alt allele
AIMS_LIST= "/scratch/user/sblain/hummer_aims/ref_panel/HummerRefPanel.f50.noUnderscores.aims"
#positions file (chr and pos) with original chromosome names, including underscores
AIMS_POS= "/scratch/user/sblain/hummer_aims/ref_panel/HummerRefPanel.f50.aims.pos"
#path to ancestry infer
A_INFER="/scratch/user/sblain/tools/ancestryinfer"
#tab-deliminted file with old scaffold name and new scaffold name on each line, to use for removing underscores from chromosome names
RENAME_SCAFFOLDS="/scratch/user/sblain/hummer_aims/ref_panel/rename_hummer_scaffolds.txt"
#allele counts for the parents
#6 tab-deliminted columns: chromosome, position, species 1 ref allele count, species 1 alt allele count, species 2 ref allele count, species 2 alt allele count
PARHMM="/scratch/user/sblain/hummer_aims/ref_panel/HummerRefPanel.f50.parhmm"

#these need to have underscores replaced with dashes, as in AIMs file
ZCHROM="NC-044274-1" #if >1 scaffold, structure as "scaffold1|scaffold2"
WCHROM="NC-044276-1"

#full path to indexed reference genome
REF="/scratch/user/sblain/hummer_aims/identifyAIMs/hummer_refs/bCalAnn1_v1.fasta"

PROJECT="HummerHybrids" #label to attach to files
F_TAG="f50"  #keep this format, f is "frequency", the number is FREQUENCY*100, where FREQUENCY varies between 0 and 1 (probably 0.3-0.9) - same as in AIMS_LIST

rule all:
	input:
		expand("{project}.{f_tag}.hmmFormat", project=PROJECT, f_tag=F_TAG)


rule sort_bam:
	input:
		bam=expand("hybrid_bams/{{sample}}.combo.bam", sample=SAMPLES),
		bai=expand("hybrid_bams/{{sample}}.combo.bam.bai", sample=SAMPLES)
	output:
		bam=expand("hybrid_bams/{{sample}}.combo.sort.bam", sample=SAMPLES)
	shell:
		"samtools sort {input.bam} -o {output.bam}"

rule index_sorted_bams:
	input:
		bam=expand("hybrid_bams/{{sample}}.combo.sort.bam", sample=SAMPLES)
	output:
		bai=expand("hybrid_bams/{{sample}}.combo.sort.bam.bai", sample=SAMPLES)
	shell:
		"samtools index {input.bam}"

rule bcftools_call:
	input:
		bam=expand("hybrid_bams/{{sample}}.combo.sort.bam", sample=SAMPLES),
		bai=expand("hybrid_bams/{{sample}}.combo.sort.bam.bai", sample=SAMPLES),
		ref=REF
	threads: 4
	output:
		vcf=expand("hybrid_counts/{{sample}}.vcf.gz", sample=SAMPLES)
	shell:
		"bcftools mpileup -f {input.ref} {input.bam} | bcftools call -m -Oz -o {output.vcf}"
		

rule index_vcf:
	input:
		vcf=expand("hybrid_counts/{{sample}}.vcf.gz", sample=SAMPLES)
	output:
		expand("hybrid_counts/{{sample}}.vcf.gz.tbi", sample=SAMPLES)
	shell:
		"tabix {input}"

#clean up bams
rule remove_intermediate_files_bam:
	input:
		bam1=expand("hybrid_bams/{{sample}}.combo.bam", sample=SAMPLES),
		bai1=expand("hybrid_bams/{{sample}}.combo.bam.bai", sample=SAMPLES),
		bam2=expand("hybrid_bams/{{sample}}.combo.sort.bam", sample=SAMPLES),
		bai2=expand("hybrid_bams/{{sample}}.combo.sort.bam.bai", sample=SAMPLES),
		#require vcf to exist before removal
		vcf=expand("hybrid_counts/{{sample}}.vcf.gz.tbi", sample=SAMPLES)
	params:
		sample=expand("{{sample}}", sample=SAMPLES)
	output:
		expand("errors/rm_bams_{{sample}}.txt", sample=SAMPLES)
	run:
		shell("rm {input.bam1} {input.bai1} {input.bam2} {input.bai2}")
		shell("echo {params.sample} 'bams removed' >> {output}")

#select AIMs and remove underscores from scaffold names
rule vcf_aims:
	input:
		vcf=expand("hybrid_counts/{{sample}}.vcf.gz", sample=SAMPLES),
		tbi=expand("hybrid_counts/{{sample}}.vcf.gz.tbi", sample=SAMPLES),
		rm_bams=expand("errors/rm_bams_{{sample}}.txt", sample=SAMPLES)
	params:
		aims=AIMS_POS,
		scaf=RENAME_SCAFFOLDS
	threads: 4
	output:
		expand("hybrid_counts/{{sample}}.noUnderscores.vcf", sample=SAMPLES)
	shell:
		"bcftools view -R {params.aims} {input.vcf} | bcftools annotate --rename-chrs {params.scaf} > {output}"

#bcftools view -R /scratch/user/sblain/hummer_aims/ref_panel/HummerRefPanel.f50.aims.pos hybrid_counts/AE18K01.vcf.gz
#bcftools annotate --rename-chrs {params.scaf} > {output}"

#clean up vcfs etc
rule remove_intermediate_files_vcf:
	input:
		vcf=expand("hybrid_counts/{{sample}}.vcf.gz", sample=SAMPLES),
		tbi=expand("hybrid_counts/{{sample}}.vcf.gz.tbi", sample=SAMPLES),
		vcf2=expand("hybrid_counts/{{sample}}.noUnderscores.vcf", sample=SAMPLES)
	params:
		sample=expand("{{sample}}", sample=SAMPLES)
	output:
		expand("errors/rm_vcf_{{sample}}.txt", sample=SAMPLES)
	run:
		shell("rm {input.vcf} {input.tbi}")
		shell("echo {params.sample} 'vcfs removed' >> {output}")

#runs Schumer lab perl script to get counts
#this script uses DP4 values to get counts
#format: 0,0,0,0 - first two ref (reverse and forward), second alt (reverse and forward)
#so DP4=0,2,3,1 translates to a count of 2 4 and DP4=5,1,0,0 translates to a count of 6 0
rule vcf_to_counts:
	input:
		vcf=expand("hybrid_counts/{{sample}}.noUnderscores.vcf", sample=SAMPLES),
		aims=expand("{aims_list}", aims_list=AIMS_LIST),
		rm_vcf=expand("errors/rm_vcf_{{sample}}.txt", sample=SAMPLES)
	params:
		script=expand("{a_infer}/vcf_to_counts_non-colinear.pl", a_infer=A_INFER),
		path=expand("{a_infer}", a_infer=A_INFER)
	threads: 4
	output:
		expand("hybrid_counts/{{sample}}.noUnderscores.vcf_counts", sample=SAMPLES)
	shell:
		"perl {params.script} {input.vcf} {input.aims} {params.path}"

rule counts_to_bed:
	input:
		expand("hybrid_counts/{{sample}}.noUnderscores.vcf_counts", sample=SAMPLES)
	output:
		expand("hybrid_counts/{{sample}}.noUnderscores.vcf_counts.bed", sample=SAMPLES)
	shell:
		"sed 's/_/\t/g' {input} | awk '{{print $1,$2,$2,$3,$4,$5,$6}}' OFS='\t' > {output}"

#find overlap (and then difference) between AIMs bed file and individual bed file
#then make one big happy bed file
rule bedtools_intersections:
	input:
		aims=expand("{aims_list}.mod.bed", aims_list=AIMS_LIST),
		bed=expand("hybrid_counts/{{sample}}.noUnderscores.vcf_counts.bed", sample=SAMPLES)
	output:
		bed=expand("hybrid_counts/{{sample}}.noUnderscores.vcf_counts.bed_combined", sample=SAMPLES)
	run:
		shell("bedtools intersect -a {input.aims} -b {input.bed} -wb -f 1 > {output.bed}")
		shell("bedtools intersect -a {input.aims} -b {input.bed} -wb -v -f 1 | awk '{{print $1,$2,$3,$4,$5,$1,$2,$3,$4,$5}}' OFS='\t' |  sed 's/\s*$/\t0\t0/' >> {output.bed}")

rule reformat_counts:
	input:
		aims=expand("{aims_list}.mod.bed", aims_list=AIMS_LIST),
		bed=expand("hybrid_counts/{{sample}}.noUnderscores.vcf_counts.bed_combined", sample=SAMPLES)
	params:
		script=expand("{a_infer}/merge_files_using_two_columns_sharing_values_stdout.pl", a_infer=A_INFER)
	output:
		counts=expand("hybrid_counts/{{sample}}.noUnderscores.vcf_counts.hmmCounts_{f_tag}_{project}", sample=SAMPLES,f_tag=F_TAG,project=PROJECT)
	shell:
		"perl {params.script} {input.aims} 0 1 {input.bed} 0 1 | cut -f 1-5 --complement | cut -f 11,12 > {output.counts}"

rule format_hmm:
	input:
		aims=expand("{aims_list}.pos", aims_list=AIMS_LIST),
		counts=expand("hybrid_counts/{sample}.noUnderscores.vcf_counts.hmmCounts_{f_tag}_{project}", sample=SAMPLES,f_tag=F_TAG,project=PROJECT)
	params:
		sample=expand("{sample}", sample=SAMPLES)
	output:
		hmmFormat=expand("{project}.{f_tag}.readCounts", project=PROJECT, f_tag=F_TAG),
		indList=expand("{project}.{f_tag}.readCounts.ind", project=PROJECT, f_tag=F_TAG)
	run:
		shell("paste -d '\t' {input.aims} {input.counts} > {output.hmmFormat}")
		shell("echo {params.sample} | perl -pe 's/ /\n/g' > {output.indList}")

rule merge_hmm_inputs:
	input:
		aimsFile=AIMS_LIST,
		aims_cmsFile=expand("{aims_list}.cMs", aims_list=AIMS_LIST),
		parhmmFile=PARHMM,
		hybridCountsFile=expand("{project}.{f_tag}.readCounts", project=PROJECT, f_tag=F_TAG),
		indsFile=expand("{project}.{f_tag}.readCounts.ind", project=PROJECT, f_tag=F_TAG)
	params:
		projectName=PROJECT,
		f_tag=F_TAG
	output:
		counts=expand("ancestryHMM_counts_.{project}.{f_tag}", project=PROJECT, f_tag=F_TAG),
		samples=expand("ancestryHMM_samples_.{project}.{f_tag}", project=PROJECT, f_tag=F_TAG)
	shell:
		"python scripts/merge_hmm_dfs.py {params.projectName} {params.f_tag} {input.aimsFile} {input.aims_cmsFile} {input.parhmmFile} {input.hybridCountsFile} {input.indsFile}"

rule remove_sex_chroms:
	input:
		expand("ancestryHMM_counts_{project}_{f_tag}", project=PROJECT, f_tag=F_TAG)
	params:
		Zchrom=ZCHROM,
		Wchrom=WCHROM
    output:
		autosomes=expand("ancestryHMM_counts_{project}_{f_tag}_autosomes", project=PROJECT, f_tag=F_TAG),
		Zchrom=expand("ancestryHMM_counts_{project}_{f_tag}_Z", project=PROJECT, f_tag=F_TAG)
	run:
		shell(grep -v -E {params.Zchrom} {input} | grep -v -E {params.Wchrom} > {output.autosomes})
		shell(grep -E {params.Zchrom} {input} > {output.Zchrom})

rule remove_intermediate_files_counts:
	input:
		counts1=expand("hybrid_counts/{{sample}}.noUnderscores.vcf_counts", sample=SAMPLES),
		bed1=expand("hybrid_counts/{{sample}}.noUnderscores.vcf_counts.bed", sample=SAMPLES),
		bed2=expand("hybrid_counts/{{sample}}.noUnderscores.vcf_counts.bed_combined", sample=SAMPLES),
		counts2=expand("hybrid_counts/{{sample}}.noUnderscores.vcf_counts.hmmCounts_{f_tag}_{project}", sample=SAMPLES,f_tag=F_TAG,project=PROJECT),
		#make sure formatted inputs exist before deleting other files
		hmmFormat=expand("{project}.{f_tag}.hmmFormat", project=PROJECT, f_tag=F_TAG),
		indList=expand("{project}.{f_tag}.hmmFormat.ind", project=PROJECT, f_tag=F_TAG)
	params:
		sample=expand("{{sample}}", sample=SAMPLES)
	output:
		expand("errors/rm_counts_{{sample}}.txt", sample=SAMPLES)
	run:
		shell("rm {input.counts1} {input.bed1} {input.bed2} {input.counts2}")
		shell("echo {params.sample} 'counts removed' >> {output}")
		
