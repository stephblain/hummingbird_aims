#Fastqs to a list of AIMs based on allele frequency differences between a reference panel


#list of samples to run
SAMPLES = ["rubythroated_SRR12247318", "blackchinned_SRR12247328", "blackchinned_SRR12247340", "blackchinned_SRR12247342", "blackchinned_SRR12247344", "blackchinned_SRR12247345", "blackchinned_SRR12247346", "blackchinned_SRR12247347", "blackchinned_SRR12247348", "rubythroated_SRR12247319", "rubythroated_SRR12247320", "rubythroated_SRR12247322", "rubythroated_SRR12247325", "rubythroated_SRR12247326", "rubythroated_SRR12247327", "rubythroated_SRR12247324", "blackchinned_SRR12247329", "rubythroated_SRR12247323"]

#list of scaffolds to include in final vcf
SCAFFOLDS = ["NC_044244.1", "NC_044245.1", "NC_044246.1", "NC_044247.1", "NC_044248.1", "NC_044249.1", "NC_044250.1", "NC_044251.1", "NC_044252.1", "NC_044253.1", "NC_044254.1", "NC_044255.1", "NC_044256.1", "NC_044257.1", "NC_044258.1", "NC_044259.1", "NC_044260.1", "NC_044261.1", "NC_044262.1", "NC_044263.1", "NC_044264.1", "NC_044265.1", "NC_044266.1", "NC_044267.1", "NC_044268.1", "NC_044269.1", "NC_044270.1", "NC_044271.1", "NC_044272.1", "NC_044273.1", "NC_044274.1", "NC_044276.1", "NC_044275.1", "NW_022045418.1", "NW_022045419.1", "NW_022045420.1", "NW_022045421.1", "NW_022045422.1", "NW_022045423.1", "NW_022045424.1", "NW_022045425.1", "NW_022045426.1", "NW_022045427.1", "NW_022045428.1", "NW_022045429.1", "NW_022045430.1", "NW_022045431.1", "NW_022045432.1", "NW_022045433.1", "NW_022045434.1", "NW_022045435.1", "NW_022045436.1", "NW_022045437.1", "NW_022045438.1", "NW_022045439.1", "NW_022045440.1", "NW_022045441.1", "NW_022045442.1", "NW_022045443.1", "NW_022045444.1", "NW_022045445.1", "NW_022045446.1", "NW_022045447.1", "NW_022045448.1", "NW_022045449.1", "NW_022045450.1", "NW_022045451.1", "NW_022045452.1", "NW_022045453.1", "NW_022045454.1", "NW_022045455.1", "NW_022045456.1", "NW_022045457.1", "NW_022045458.1", "NW_022045459.1", "NW_022045460.1", "NW_022045461.1", "NW_022045462.1", "NW_022045463.1", "NW_022045464.1", "NW_022045465.1", "NW_022045466.1", "NW_022045467.1", "NW_022045468.1", "NW_022045469.1", "NW_022045470.1", "NW_022045471.1", "NW_022045472.1", "NW_022045473.1", "NW_022045474.1", "NW_022045475.1", "NW_022045476.1", "NW_022045477.1", "NW_022045478.1", "NW_022045479.1", "NW_022045480.1", "NW_022045481.1", "NW_022045482.1", "NW_022045483.1", "NW_022045484.1", "NW_022045485.1", "NW_022045486.1", "NW_022045487.1", "NW_022045488.1", "NW_022045489.1", "NW_022045490.1", "NW_022045491.1", "NW_022045492.1", "NW_022045493.1", "NW_022045494.1", "NW_022045495.1", "NW_022045496.1", "NW_022045497.1", "NW_022045498.1", "NW_022045499.1", "NW_022045500.1", "NW_022045501.1", "NW_022045502.1", "NW_022045503.1", "NW_022045504.1", "NW_022045505.1", "NW_022045506.1", "NW_022045507.1", "NW_022045508.1", "NW_022045509.1", "NW_022045510.1", "NW_022045511.1", "NW_022045512.1", "NW_022045513.1", "NW_022045514.1", "NW_022045515.1", "NW_022045516.1", "NW_022045517.1", "NW_022045518.1", "NW_022045519.1", "NW_022045520.1", "NW_022045521.1", "NW_022045522.1", "NW_022045523.1", "NW_022045524.1", "NW_022045525.1", "NW_022045526.1", "NW_022045527.1", "NW_022045528.1", "NW_022045529.1", "NW_022045530.1", "NW_022045531.1", "NW_022045532.1", "NW_022045533.1", "NW_022045534.1", "NW_022045535.1", "NW_022045536.1", "NW_022045537.1", "NW_022045538.1", "NW_022045539.1", "NW_022045540.1", "NW_022045541.1", "NW_022045542.1", "NW_022045543.1"]
SCAFFOLDS_AIMS = ["NC-044244-1", "NC-044245-1", "NC-044246-1", "NC-044247-1", "NC-044248-1", "NC-044249-1", "NC-044250-1", "NC-044251-1", "NC-044252-1", "NC-044253-1", "NC-044254-1", "NC-044255-1", "NC-044256-1", "NC-044257-1", "NC-044258-1", "NC-044259-1", "NC-044260-1", "NC-044261-1", "NC-044262-1", "NC-044263-1", "NC-044264-1", "NC-044265-1", "NC-044266-1", "NC-044267-1", "NC-044268-1", "NC-044269-1", "NC-044270-1", "NC-044271-1", "NC-044272-1", "NC-044273-1", "NC-044274-1", "NC-044275-1", "NC-044276-1", "NW-022045418-1", "NW-022045420-1", "NW-022045422-1", "NW-022045423-1", "NW-022045425-1", "NW-022045426-1", "NW-022045427-1", "NW-022045428-1", "NW-022045430-1", "NW-022045431-1", "NW-022045434-1", "NW-022045438-1", "NW-022045439-1", "NW-022045440-1", "NW-022045441-1", "NW-022045442-1", "NW-022045443-1", "NW-022045444-1", "NW-022045445-1", "NW-022045446-1", "NW-022045447-1", "NW-022045448-1", "NW-022045449-1", "NW-022045450-1", "NW-022045451-1", "NW-022045452-1", "NW-022045453-1", "NW-022045455-1", "NW-022045456-1", "NW-022045457-1", "NW-022045458-1", "NW-022045461-1", "NW-022045462-1", "NW-022045465-1", "NW-022045469-1", "NW-022045470-1", "NW-022045472-1", "NW-022045473-1", "NW-022045474-1", "NW-022045475-1", "NW-022045479-1", "NW-022045481-1", "NW-022045482-1", "NW-022045483-1", "NW-022045484-1", "NW-022045485-1", "NW-022045486-1", "NW-022045487-1", "NW-022045488-1", "NW-022045489-1", "NW-022045490-1", "NW-022045493-1", "NW-022045495-1", "NW-022045496-1", "NW-022045497-1", "NW-022045499-1", "NW-022045500-1", "NW-022045501-1", "NW-022045502-1", "NW-022045503-1", "NW-022045506-1", "NW-022045508-1", "NW-022045511-1", "NW-022045512-1", "NW-022045513-1", "NW-022045514-1", "NW-022045517-1", "NW-022045518-1", "NW-022045519-1", "NW-022045522-1", "NW-022045524-1", "NW-022045525-1", "NW-022045526-1", "NW-022045527-1", "NW-022045528-1", "NW-022045529-1", "NW-022045530-1", "NW-022045532-1", "NW-022045533-1", "NW-022045534-1", "NW-022045535-1", "NW-022045536-1", "NW-022045537-1", "NW-022045539-1", "NW-022045540-1", "NW-022045541-1"]

#full path to indexed reference genome
REF="/scratch/user/sblain/hummer_aims/identifyAIMs/hummer_refs/bCalAnn1_v1.fasta"

PROJECT="HummerRefPanel" #label to attach to files

#path to lists of samples (one per row) divided by species
SP1_SAMPLES = "/scratch/user/sblain/hummer_aims/ref_panel/variantCalling/inds_rubythroated"
SP2_SAMPLES = "/scratch/user/sblain/hummer_aims/ref_panel/variantCalling/inds_blackchinned"

#species labels - make sure these correspond to SP_SAMPLES lists
SP1_NAME="rubythroated"
SP2_NAME="blackchinned"

FREQUENCY="0.5" #minimum allele frequency difference to designate an AIM
F_TAG="f50" #keep this format, f is "frequency", the number is FREQUENCY*100

RECOMB_M_BP=0.0000000369 #recombination rate in morgans per bp, Flycatchers are 3.69 cM/Mb (Kawakami et al. 2017)


#Fastqs to bams pipeline adapted from: https://github.com/kdelmore/delmore_lab/blob/master/SWTH/ref_panel/01_trim_align.sh

#Before running, load modules for various programs (Trim Galore, FastQC, Bcftools, Vcftools, snakemake, picard)
#Ordered by dependencies for TAMU Grace cluster:
##module load GCCcore/11.2.0 Trim_Galore/0.6.7 FastQC/0.11.9-Java-11
##module load GCC/11.2.0 BCFtools/1.14 VCFtools/0.1.16
##module load OpenMPI/4.1.1 BWA/0.7.17 picard/2.25.1-Java-11
##module load snakemake/6.10.0 Python/3.9.6 #pandas installed in python with pip


#test as: snakemake -np
#run as:  nohup snakemake --profile . &
#before running, make config.yaml file in directory

#config.yaml example:
##cluster: "sbatch --time=24:00:00 --mem-per-cpu=80G
##						--nodes=1 --ntasks=1 --output=errors/err_{rule}_{wildcards}.%j"
##jobs: 18
#possible solution to weird parallelization: https://stackoverflow.com/questions/68541814/snakemake-waiting-to-finish-all-parallel-jobs-before-starting-next-parallel-job

################################
# Trim and align fastqs
################################

rule all:
	input:
		expand("{project}.{f_tag}.noUnderscores.aims.cMs", project=PROJECT, f_tag=F_TAG)

rule trim_reads:
	input:
		fq1=expand("{{sample}}_1.fastq", sample=SAMPLES),
		fq2=expand("{{sample}}_2.fastq", sample=SAMPLES)
	output:
		val1=expand("{{sample}}_1_val_1.fq", sample=SAMPLES),
		val2=expand("{{sample}}_2_val_2.fq", sample=SAMPLES),
		fq1=expand("{{sample}}_1_unpaired_1.fq", sample=SAMPLES),
		fq2=expand("{{sample}}_2_unpaired_2.fq", sample=SAMPLES)
	shell:
		"trim_galore --paired -fastqc --clip_R1 15 --clip_R2 15 --three_prime_clip_R1 5 --three_prime_clip_R2 5 --retain_unpaired {input.fq1} {input.fq2}"

rule bwa_align:
	input:
		val1=expand("{{sample}}_1_val_1.fq", sample=SAMPLES),
		val2=expand("{{sample}}_2_val_2.fq", sample=SAMPLES),
		fq1=expand("{{sample}}_1_unpaired_1.fq", sample=SAMPLES),
		fq2=expand("{{sample}}_2_unpaired_2.fq", sample=SAMPLES),
		ref=REF
	output:
		sam=expand("{{sample}}.sam", sample=SAMPLES),
		log=expand("{{sample}}.bwape.log", sample=SAMPLES),
		sam1=expand("{{sample}}_1_unpaired.sam", sample=SAMPLES),
		log1=expand("{{sample}}.bwase1.log", sample=SAMPLES),
		sam2=expand("{{sample}}_2_unpaired.sam", sample=SAMPLES),
		log2=expand("{{sample}}.bwase2.log", sample=SAMPLES)
	run:
		shell("bwa mem -M -t 2 {input.ref} {input.val1} {input.val2} > {output.sam} 2> {output.log}")
		shell("bwa mem -M {input.ref} {input.fq1} > {output.sam1} 2> {output.log1}")
		shell("bwa mem -M {input.ref} {input.fq2} > {output.sam2} 2> {output.log2}")

rule sam_to_bam:
	input:
		sam=expand("{{sample}}.sam", sample=SAMPLES),
		sam1=expand("{{sample}}_1_unpaired.sam", sample=SAMPLES),
		sam2=expand("{{sample}}_2_unpaired.sam", sample=SAMPLES)
	output:
		bam=expand("{{sample}}.bam", sample=SAMPLES),
		log=expand("{{sample}}.sampe.log", sample=SAMPLES),
		bam1=expand("{{sample}}_1_unpaired.bam", sample=SAMPLES),
		log1=expand("{{sample}}.samse1.log", sample=SAMPLES),
		bam2=expand("{{sample}}_2_unpaired.bam", sample=SAMPLES),
		log2=expand("{{sample}}.samse2.log", sample=SAMPLES)
	run:
		shell("samtools view -Sb {input.sam} > {output.bam} 2> {output.log}")
		shell("samtools view -Sb {input.sam1} > {output.bam1} 2> {output.log1}")
		shell("samtools view -Sb {input.sam2} > {output.bam2} 2> {output.log2}")


################################
# Sorting, cleaning, etc. with picard
################################


rule clean_bam:
	input:
		bam=expand("{{sample}}.bam", sample=SAMPLES),
		bam1=expand("{{sample}}_1_unpaired.bam", sample=SAMPLES),
		bam2=expand("{{sample}}_2_unpaired.bam", sample=SAMPLES)
	output:
		bam=expand("{{sample}}.clean.bam", sample=SAMPLES),
		log=expand("{{sample}}.cleansam.log", sample=SAMPLES),
		bam1=expand("{{sample}}_1_unpaired.clean.bam", sample=SAMPLES),
		log1=expand("{{sample}}.cleansam1.log", sample=SAMPLES),
		bam2=expand("{{sample}}_2_unpaired.clean.bam", sample=SAMPLES),
		log2=expand("{{sample}}.cleansam2.log", sample=SAMPLES)
	run:
		shell("java -jar $EBROOTPICARD/picard.jar CleanSam INPUT={input.bam} OUTPUT={output.bam}  2> {output.log}")
		shell("java -jar $EBROOTPICARD/picard.jar CleanSam INPUT={input.bam1} OUTPUT={output.bam1}  2> {output.log1}")
		shell("java -jar $EBROOTPICARD/picard.jar CleanSam INPUT={input.bam2} OUTPUT={output.bam2}  2> {output.log2}")
		
rule sort_sam:
	input:
		bam=expand("{{sample}}.clean.bam", sample=SAMPLES),
		bam1=expand("{{sample}}_1_unpaired.clean.bam", sample=SAMPLES),
		bam2=expand("{{sample}}_2_unpaired.clean.bam", sample=SAMPLES)
	output:
		bam=expand("{{sample}}.sort.bam", sample=SAMPLES),
		log=expand("{{sample}}.sortsam.log", sample=SAMPLES),
		bam1=expand("{{sample}}_1_unpaired.sort.bam", sample=SAMPLES),
		log1=expand("{{sample}}.sortsam1.log", sample=SAMPLES),
		bam2=expand("{{sample}}_2_unpaired.sort.bam", sample=SAMPLES),
		log2=expand("{{sample}}.sortsam2.log", sample=SAMPLES)
	run:
		shell("java -jar $EBROOTPICARD/picard.jar SortSam SO=coordinate INPUT={input.bam} OUTPUT={output.bam}  2> {output.log}")
		shell("java -jar $EBROOTPICARD/picard.jar SortSam SO=coordinate INPUT={input.bam1} OUTPUT={output.bam1}  2> {output.log1}")
		shell("java -jar $EBROOTPICARD/picard.jar SortSam SO=coordinate INPUT={input.bam2} OUTPUT={output.bam2}  2> {output.log2}")

rule add_replace_read_groups:
	input:
		bam=expand("{{sample}}.sort.bam", sample=SAMPLES),
		bam1=expand("{{sample}}_1_unpaired.sort.bam", sample=SAMPLES),
		bam2=expand("{{sample}}_2_unpaired.sort.bam", sample=SAMPLES)
	params:
		project=PROJECT,
		sampleid=expand("{{sample}}", sample=SAMPLES)
	output:
		bam=expand("{{sample}}.sortrg.bam", sample=SAMPLES),
		log=expand("{{sample}}.addRG.log", sample=SAMPLES),
		bam1=expand("{{sample}}_1_unpaired.sortrg.bam", sample=SAMPLES),
		log1=expand("{{sample}}.addRG1.log", sample=SAMPLES),
		bam2=expand("{{sample}}_2_unpaired.sortrg.bam", sample=SAMPLES),
		log2=expand("{{sample}}.addRG2.log", sample=SAMPLES)
	run:
		shell("java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I={input.bam} O={output.bam} MAX_RECORDS_IN_RAM=5000000 SORT_ORDER=coordinate RGID={params.sampleid} RGLB={params.project} RGPL=ILLUMINA RGPU={params.project} RGSM={params.sampleid} CREATE_INDEX=True 2> {output.log}")
		shell("java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I={input.bam1} O={output.bam1} SORT_ORDER=coordinate RGID={params.sampleid} RGLB={params.project} RGPL=ILLUMINA RGPU={params.project} RGSM={params.sampleid} CREATE_INDEX=True 2> {output.log1}")
		shell("java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I={input.bam2} O={output.bam2} SORT_ORDER=coordinate RGID={params.sampleid} RGLB={params.project} RGPL=ILLUMINA RGPU={params.project} RGSM={params.sampleid} CREATE_INDEX=True 2> {output.log2}")

rule mark_duplicates:
	input:
		bam=expand("{{sample}}.sortrg.bam", sample=SAMPLES)
	output:
		bam=expand("{{sample}}.duprem.bam", sample=SAMPLES),
		log=expand("{{sample}}.duprem.log", sample=SAMPLES)	
	shell:
		"java -Xmx7g -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT={input.bam} MAX_RECORDS_IN_RAM=5000000 OUTPUT={output.bam} M={output.log} REMOVE_DUPLICATES=true 2> {output.log}"

################################
# Merge, index, and sort final bam
################################

rule merge_bams:
	input:
		bam=expand("{{sample}}.duprem.bam", sample=SAMPLES),
		bam1=expand("{{sample}}_1_unpaired.sortrg.bam", sample=SAMPLES),
		bam2=expand("{{sample}}_2_unpaired.sortrg.bam", sample=SAMPLES)
	output:
		bam=expand("{{sample}}.combo.bam", sample=SAMPLES),
		log=expand("{{sample}}.sammerge.log", sample=SAMPLES)
	shell:
		"samtools merge {output.bam} {input.bam} {input.bam1} {input.bam2} > {output.log}"


rule index_bams:
	input:
		bam=expand("{{sample}}.combo.bam", sample=SAMPLES)
	output:
		bai=expand("{{sample}}.combo.bam.bai", sample=SAMPLES)
	shell:
		"samtools index {input.bam}"
		
rule sort_bam:
	input:
		bam=expand("{{sample}}.combo.bam", sample=SAMPLES),
		bai=expand("{{sample}}.combo.bam.bai", sample=SAMPLES)
	output:
		bam=expand("{{sample}}.combo.sort.bam", sample=SAMPLES)		
	shell:
		"samtools sort {input.bam} -o {output.bam}"
		

rule index_sorted_bams:
	input:
		bam=expand("{{sample}}.combo.sort.bam", sample=SAMPLES)
	output:
		bai=expand("{{sample}}.combo.sort.bam.bai", sample=SAMPLES)
	shell:
		"samtools index {input.bam}"


################################
# haplotype calling with bcftools
################################

#see example in d stats folder
rule bcftools_call:
	input:
		bam=expand("{sample}.combo.sort.bam", sample=SAMPLES),
		bai=expand("{sample}.combo.sort.bam.bai", sample=SAMPLES),
		ref=REF
	params:
		scaf=expand("{{scaffold}}", scaffold=SCAFFOLDS)
	threads: 4
	output:
		bcf=expand("{project}.{{scaffold}}.bcf", scaffold=SCAFFOLDS, project=PROJECT)
	shell:
		"bcftools mpileup -r {params.scaf} -f {input.ref} {input.bam} | bcftools call -mO b -o {output.bcf}"
		
#concatenate bcf files
#keep only biallelic SNPs
#filter out sites with > 25% missing data or a minor allele frequency < 0.05
rule bcftools_filter:
	input:
		bcf=expand("{project}.{scaffold}.bcf", scaffold=SCAFFOLDS, project=PROJECT)
	output:
		vcf=expand("{project}.filtered.vcf.gz",project=PROJECT)
	shell:
		"bcftools concat {input.bcf} | bcftools view -m2 -M2 -v snps | bcftools filter -e 'F_MISSING > 0.25 || MAF <= 0.05' -O z -o {output.vcf}"

rule index_vcf:
	input:
		expand("{project}.filtered.vcf.gz",project=PROJECT)
	output:
		expand("{project}.filtered.vcf.gz.tbi",project=PROJECT)
	shell:
		"tabix {input}"

################################
# species allele counts
################################

rule allele_counts:
	input:
		vcf=expand("{project}.filtered.vcf.gz",project=PROJECT),
		tbi=expand("{project}.filtered.vcf.gz.tbi",project=PROJECT),
		sp1=SP1_SAMPLES,
		sp2=SP2_SAMPLES
	output:
		counts1=expand("{project}.{sp1_name}.counts",project=PROJECT,sp1_name=SP1_NAME),
		counts2=expand("{project}.{sp2_name}.counts",project=PROJECT,sp2_name=SP2_NAME)
	run:
		shell("bcftools view -S {input.sp1} {input.vcf} | bcftools query -f '%CHROM %POS %REF %ALT %AN %AC\n' > {output.counts1}")
		shell("bcftools view -S {input.sp2} {input.vcf} | bcftools query -f '%CHROM %POS %REF %ALT %AN %AC\n' > {output.counts2}")

#select aims based on a minimum frequency difference between species
#output lists of aims (with and without underscores in scaffolds) and parental subspecies input for ancestryHMM
rule select_aims:
	input:
		counts1=expand("{project}.{sp1_name}.counts",project=PROJECT,sp1_name=SP1_NAME),
		counts2=expand("{project}.{sp2_name}.counts",project=PROJECT,sp2_name=SP2_NAME)
	params:
		freq=FREQUENCY,
		project=PROJECT
	output:
		expand("{project}.{f_tag}.aims",project=PROJECT,f_tag=F_TAG),
		expand("{project}.{f_tag}.aims.pos",project=PROJECT,f_tag=F_TAG),
		expand("{project}.{f_tag}.noUnderscores.aims",project=PROJECT,f_tag=F_TAG)
	shell:
		"python scripts/get_ref_counts.py {params.freq} {params.project} {input.counts1} {input.counts2}"

rule bp_to_cMs:
	input:
		aims=expand("{project}.{f_tag}.noUnderscores.aims",project=PROJECT,f_tag=F_TAG)
	params:
		recombRate=RECOMB_M_BP
	output:
		cMs=expand("{project}.{f_tag}.noUnderscores.aims_{{scaffold}}.cMs", scaffold=SCAFFOLDS_AIMS, project=PROJECT, f_tag=F_TAG)
	shell:
		"python scripts/aims_bp_to_cms.py {input.aims} {params.recombRate}"

rule combine_cMs:
	input:
		expand("{project}.{f_tag}.noUnderscores.aims_{scaffold}.cMs", scaffold=SCAFFOLDS_AIMS, project=PROJECT, f_tag=F_TAG)
	output:
		expand("{project}.{f_tag}.noUnderscores.aims.cMs", project=PROJECT, f_tag=F_TAG)
	run:
		shell("cat {input} > {output}")
		shell("rm {input}")
#