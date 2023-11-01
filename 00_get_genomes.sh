module load GCC/11.3.0  OpenMPI/4.1.4  SRA-Toolkit/3.0.3

#get blackchinned
prefetch SRR12247329
cd SRR12247329
fasterq-dump --split-files SRR12247329.sra
mv SRR12247329_1.fastq ../identifyAIMs/blackchinned_SRR12247329_1.fastq
mv SRR12247329_2.fastq ../identifyAIMs/blackchinned_SRR12247329_2.fastq
rm -r SRR12247329

#get rubythroated
prefetch SRR12247324
fasterq-dump --split-files SRR12247324/SRR12247324.sra
mv SRR12247324_1.fastq identifyAIMs/rubythroated_SRR12247324_1.fastq
mv SRR12247324_2.fastq identifyAIMs/rubythroated_SRR12247324_2.fastq
rm -r SRR12247324

#get Anna's hummingbird reference
#downloaded from https://www.ncbi.nlm.nih.gov/assembly/GCA_003957555.2
# on 2023-06-13
tar -xvf genome_assemblies_genome_fasta.tar
mv ncbi-genomes-2023-06-13/GCF_003957555.1_bCalAnn1_v1.p_genomic.fna.gz identifyAIMs/bCalAnn1_v1.fasta.gz
rm -r ncbi-genomes-2023-06-13/




#get blackchinned ref panel

cd /scratch/user/sblain/hummer_aims/ref_panel

for SRR in SRR12247328 SRR12247340 SRR12247342 SRR12247344 SRR12247345 SRR12247346 SRR12247347 SRR12247348
do
	prefetch $SRR
	fasterq-dump --split-files $SRR/$SRR.sra
	mv "$SRR"_1.fastq blackchinned_"$SRR"_1.fastq
	mv "$SRR"_2.fastq blackchinned_"$SRR"_2.fastq
	rm -r $SRR
done


#get rubythroated ref panel

cd /scratch/user/sblain/hummer_aims/ref_panel

for SRR in SRR12247318 SRR12247319 SRR12247320 SRR12247322 SRR12247323 SRR12247325 SRR12247326 SRR12247327

do
	prefetch $SRR
	fasterq-dump --split-files $SRR/$SRR.sra
	mv "$SRR"_1.fastq rubythroated_"$SRR"_1.fastq
	mv "$SRR"_2.fastq rubythroated_"$SRR"_2.fastq
	rm -r $SRR
done

