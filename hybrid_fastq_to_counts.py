# import pandas as pd

# samples_df = pd.read_table('hummer_samples.tsv').set_index("sample", drop=False)
# SAMPLES= list(samples_df['sample'])

SAMPLES=["AE18K01", "AE18K02", "AE18K03", "AE18K04", "AE18K05", "AE18K06", "AE18K07", "AE18K08", "AE18K11", "AE20K09", "AE20K16", "AE24K02", "AE24K03", "AE24K04", "AE24K05", "AE24K06", "AE24K07", "AE24K08", "AE24K09", "AE24K16", "AE24K17", "AE24K18", "AE24K19", "AE24K20", "AE24K21", "AE24K22", "AE24K23", "AE24K24", "AE24K25", "AE24K26", "AE25K02", "AE25K03", "AE25K04", "AE25K05", "AE25K06", "AE25K07", "AE25K09", "AE25K10", "AE25K13", "AE25K14", "AE25K15", "AE25K16", "AE25K17", "AE25K18", "AE25K19", "AE25K32", "AE25K33", "AE25K34", "AE25K35", "AE25K36", "AE25K37", "AE25K38", "AE25K40", "AE25K41", "AE25K42", "AE26K01", "AE26K02", "AE26K03", "AE26K04", "AE26K05", "AE26K06", "AE26K07", "AE26K08", "AE26K09", "AE26K10", "AE26K11", "AE26K12", "AE26K13", "AE26K14", "AE26K15", "AE26K16", "AE26K17", "AE26K18", "AE26K21", "AE26K24", "AE26K31", "AE26K32", "AE26K33", "AE26K34", "AE26K35", "AE26K36", "AE27K01", "AE27K02", "AE27K03", "AE27K04", "AE27K05", "AE27K09", "AE27K10", "AE27K11", "AE27K13", "AE27K15", "AE27K16", "AE27K17", "AE27K18", "AE27K23", "AF21K01", "AF21K02", "AF21K03", "AF21K04", "AF21K05", "AF21K06", "AF22K01", "AF22K02", "AF22K03", "AF22K04", "AF22K05", "AF22K06", "AF22K07", "AF22K08", "AF22K09", "AF22K10", "AF22K11", "AF23K01", "AF23K01A", "AF23K02", "AF23K03", "AF23K04", "AF23K05", "AF23K06", "AF23K07", "AF23K08", "AF23K09", "AF23K10", "AF23K11", "AF23K12", "AF23K13", "AF24K01", "AF24K05", "AF24K11", "AF24K12", "AF24K18", "AF25K01", "AF25K03", "AF25K17", "AF25K20", "AF25K23", "AF26K25", "AF26K26", "AF26K27", "AG01K01", "AG01K02", "AG01K03", "AG01K04", "BE01K01", "BE01K02", "BE01K03", "BE01K04", "BE02K01", "BE02K02", "BE02K03", "BE03K01", "BE04K01", "BE06K01", "BE06K02", "BE06K03", "BE06K04", "BE06K05", "BE06K06", "BE09K01", "BE10K01", "BE10K02", "BE10K03", "BE10K04", "BE10K05", "BE10K06", "BE10K07", "BE10K08", "BE10K09", "BE11K02", "BE11K03", "BE11K05", "BE11K06", "BE11K07", "BE11K08", "BE11K09", "BE11K10", "BE11K11", "BE11K12", "BE11K13", "BE11K14", "BE11K15", "BE11K16", "BE11K17", "BE11K18", "BE11K19", "BE11K20", "BE11K21", "BE11K22", "BE11K23", "BE13K01", "BE13K03", "BE13K04", "BE28K02", "BE28K03", "BE29K01", "BE29K03", "BF04K01", "BF05K01", "BF16K01", "BF18K01", "BF19K01", "BF19K02", "BF23K01", "BF24K01", "BF24K02", "BF25K01", "BF27K01", "BF27K02", "BF28K01", "BF29K01", "BF29K02", "BG20K01", "CE12K01", "CE12K02", "CE19K01", "CE21K01", "CE21K02", "CE22K01", "CE30K02", "CE31K01", "CF17K01", "CF17K02", "CF17K03", "CF17K04", "CF17K05", "CF17K06", "CF17K07", "CF17K08", "CF17K09", "CF17K10", "CF20K01", "CF22K01", "CF23K01", "CF23K02", "CF23K03", "CF23K04", "CF23K05", "CF23K06", "CF23K07", "CF23K08", "CF23K09", "CF23K10", "CF23K11", "CF23K12", "CF23K13", "CF23K14", "CF23K15", "CF23K16", "CF23K17", "CF23K18", "DE14K01", "DE14K02", "DE14K03", "DE14K04", "DE14K05", "DE14K06", "DE14K07", "DE14K08", "DE14K09", "DE14K10", "DE14K11", "DE14K12", "DE14K13", "DE14K14", "DE14K15", "DE14K16", "DE14K17", "DE14K18", "DE14K19", "DE14K20", "DE14K21", "DE14K22", "DE14K23", "DE14K24", "DE14K25", "DE14K26", "DE14K27", "DE14K28", "DE14K29", "DE14K30", "DE15K01", "DE15K02", "DE15K03", "DE15K04", "DE15K05", "DE16K01", "DE16K02", "DE16K03", "DE16K04", "DE16K05", "DE16K06", "DE16K07", "DE16K08", "DE16K09", "DE16K10", "DE16K11", "DE25K03", "DE25K04", "DE28K01", "DE28K02", "DE28K03", "DE28K04", "DE28K05", "DE28K06", "DE28K07", "DE28K08", "DE28K09", "DE28K10", "DE28K11", "DE28K12", "DE28K13", "DE28K14", "DE28K15", "DE28K16", "DE28K17", "DE28K18", "DE28K19", "DE28K20", "DE28K21", "DE29K02", "DE29K05", "DF01K01", "DF06K01", "DF06K02", "DF06K03", "DF06K04", "DF06K05", "DF06K06", "DF08K09", "DF09K01", "DF09K02", "DF09K03", "DF10K01", "DG13K01", "DG13K02", "SRR12247318", "SRR12247319", "SRR12247320", "SRR12247322", "SRR12247323", "SRR12247324", "SRR12247325", "SRR12247326", "SRR12247327", "SRR12247328", "SRR12247329", "SRR12247340", "SRR12247342", "SRR12247344", "SRR12247345", "SRR12247346", "SRR12247347", "SRR12247348", "TC23787", "TC24315", "TC24316"]

#path to aims list, make sure there are no underscores in the scaffold names, for Hummers I also removed periods because I'm paranoid
AIMS_LIST= "/scratch/user/sblain/hummer_aims/ref_panel/HummerRefPanel.f50.noUnderscores.aims"
A_INFER="/scratch/user/sblain/tools/ancestryinfer" #path to ancestry infer
RENAME_SCAFFOLDS="/scratch/user/sblain/hummer_aims/ref_panel/rename_hummer_scaffolds.txt" #tab-deliminted file with old scaffold name and new scaffold name on each line

#full path to indexed reference genome
REF="/scratch/user/sblain/hummer_aims/identifyAIMs/hummer_refs/bCalAnn1_v1.fasta"

PROJECT="HummerRefPanel" #label to attach to files

#Before running, load modules for various programs (Bcftools, Vcftools, snakemake)
#Ordered by dependencies for TAMU Grace cluster:
##module load GCC/11.2.0 BCFtools/1.14 VCFtools/0.1.16 OpenMPI/4.1.1 snakemake/6.10.0

#test as: snakemake -np
#run as:  nohup snakemake --profile . &
#before running, make config.yaml file in directory

#config.yaml example:
##cluster: "sbatch --time=24:00:00 --mem-per-cpu=80G
##						--nodes=1 --ntasks=1 --output=errors/err_{rule}_{wildcards}.%j"
##jobs: 18

rule all:
	input:
		expand("hybrid_counts/{sample}.vcf.gz.tbi", sample=SAMPLES)


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


#select AIMs and remove underscores from scaffold names
rule vcf_aims:
	input:
		vcf=expand("hybrid_counts/{{sample}}.vcf.gz", sample=SAMPLES),
		tbi=expand("hybrid_counts/{{sample}}.vcf.gz.tbi", sample=SAMPLES)
	params:
		aims=AIMS_LIST
        scaf=RENAME_SCAFFOLDS
	output:
		expand("hybrid_counts/{{sample}}.noUnderscores.vcf", sample=SAMPLES)
	shell:
		"bcftools view -R {params.aims} {input.vcf} | bcftools annotate --rename-chrs {params.scaf} > {output}"

#clean up bams, vcfs, etc
rule remove_intermediate_files:
	input:
		bam=expand("hybrid_bams/{{sample}}.combo.sort.bam", sample=SAMPLES),
		bai=expand("hybrid_bams/{{sample}}.combo.sort.bam.bai", sample=SAMPLES),
		vcf=expand("hybrid_counts/{{sample}}.vcf.gz", sample=SAMPLES),
		tbi=expand("hybrid_counts/{{sample}}.vcf.gz.tbi", sample=SAMPLES)
	output:
		expand("errors/{{sample}}_rm.txt", sample=SAMPLES)
	run:
		shell("rm {input.bam} {input.bai} {input.bai} {input.vcf} {input.tbi}")
		shell("echo 'bams and vcfs removed' > {output}")

#runs Schumer lab perl script to get counts
rule vcf_to_counts:
	input:
		vcf=expand("hybrid_counts/{{sample}}.noUnderscores.vcf", sample=SAMPLES),
		tbi=expand("hybrid_counts/{{sample}}.vcf.gz.tbi", sample=SAMPLES),
		aims=expand("{a_infer}", aims_list=AIMS_LIST)
		rm_out=expand("errors/{{sample}}_rm.txt", sample=SAMPLES)
	output:
		expand("hybrid_counts/{{sample}}.vcf.gz_counts", sample=SAMPLES)
	params:
		script=expand("{a_infer}/vcf_to_counts_non-colinear.pl", a_infer=A_INFER),
		path=expand("{a_infer}", a_infer=A_INFER)
	shell:
		"perl {params.script} {input.vcf} {input.aims} {params.path}"

rule counts_to_bed:
	input:
		expand("hybrid_counts/{{sample}}.vcf.gz_counts", sample=SAMPLES)
	output:
		expand("hybrid_counts/{{sample}}.vcf.gz_counts.bed", sample=SAMPLES)
	shell:
		"cat {input} | perl -p -e 's/_/\t/g' | awk -v OFS="\t" \$1=\$1\"\\t\"\$2\ > {output}"

