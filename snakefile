#  module load  GCC/11.2.0  OpenMPI/4.1.1 snakemake/6.10.0 BWA/0.7.17 picard/2.25.1-Java-11 BCFtools/1.14


################################
# Index reference genome
################################

rule bwa_index:
    input:
        ref="bCalAnn1_v1.fasta"
    output:
        idx=multiext("bCalAnn1_v1.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
	shell:
		"bwa index {input.ref} > {output.idx}"

rule picard_dict:
	input:
        ref="bCalAnn1_v1.fasta"
	output:
		dict="bCalAnn1_v1.dict"
	shell:
		"java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R={input.ref} O={output.dict}"

rule samtools_index:
	input:
        ref="bCalAnn1_v1.fasta"
	output:
		fai="bCalAnn1_v1.fasta.fai"
	shell:
		"samtools faidx {input.ref} > {output.fai}"


################################
# Align, sort and index fastqs
################################

rule bwa_mem:
	input:
		ref="bCalAnn1_v1.fasta",
		fq1="{sample}_1.fastq",
		fq2="{sample}_2.fastq"
	output:
		sam="{sample}.sam"
	shell:
		"bwa mem -M {input.ref} {input.fq1} {input.fq2} > {output.sam}"
		
#fill in info on corresponding paired end reads
	
rule samtools_fixmate:
	input:
		sam="{sample}.sam"
	output:
		bam="{sample}.bam"
	shell:
		"samtools fixmate -O bam {input.sam} {output.bam}"

rule samtools_sort:
	input:
		bam="{sample}.bam"
	output:
		sortBam="{sample}.sorted.bam"
	shell:
		"samtools sort {input.bam} -o {output.sortBam}"
		
rule samtools_index_2:
	input:
        bam="{sample}.sorted.bam"
	output:
		bai="{sample}.sorted.bai"
	shell:
		"samtools index {input.bam} > {output.bai}"

################################
# Filter to qual > 30 and index
################################

rule samtools_filter:
	input:
        bam="{sample}.sorted.bam"
	output:
		q30Bam="{sample}.sorted.q30.bam"
	shell:
		"samtools view -b -q 30 {input.bam} > {output.q30Bam}"

rule samtools_index_3:
	input:
        bam="{sample}.sorted.q30.bam"
	output:
		bai="{sample}.sorted.q30.bai"
	shell:
		"samtools index {input.bam} > {output.bai}"
		
################################
# Call variants
################################
	
#generates genotype likelihoods at each position
rule bcftools_mpileup:
	input:
		bam="{sample}.sorted.q30.bam"
		ref="bCalAnn1_v1.fasta"
	output:
		pileup="{sample}.sorted.q30.pileup"
	shell:
		"bcftools mpileup -o {output.pileup} -f {input.ref} {input.bam}"

#actually call variants
rule bcftools_call:
	input:
		pileup="{sample}.sorted.q30.pileup"
	output:
		vcf="{sample}.sorted.q30.vcf.gz"
	shell:
		"bcftools call -mO z -o {output.vcf} {input.pileup}"

#index called variants
rule bcftools_index:
	input:
		vcf="{sample}.sorted.q30.vcf.gz"
	output:
		csi="{sample}.sorted.q30.vcf.gz.csi"
	shell:
		"samtools index {input.vcf} > {output.csi}"

################################
# filter variants
################################

#1- no indels
#2- homozygous only (will also filter to only biallelic for this dataset)
#3- other assorted thresholds (quality scoring and depth)
rule bcftools_index:
	input:
		vcfIn="{sample}.sorted.q30.vcf.gz"
	output:
		vcfOut="{sample}.sorted.q30.filtered.vcf.gz"
	shell:
		"bcftools view --exclude-types indels {input.vcfIn} | bcftools view --genotype hom | bcftools view -e 'QUAL < 30 || DP < 7 || DP > 60 || MQ < 40 ' -o {output.vcfOut}"

################################
# get a consensus sequence (pseudoref)
################################

rule bcftools_consensus:
	input:
		vcf="{sample}.sorted.q30.filtered.vcf.gz"
		ref="bCalAnn1_v1.fasta"
	output:
		fasta="{sample}_consensus.fasta"
	shell:
		cat {input.ref} | bcftools consensus {input.vcf} > {output.fasta}

