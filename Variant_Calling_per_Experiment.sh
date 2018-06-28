#!/bin/bash

#sample list - same as before but will feasibly need to be edited to accommodate arbitrary input lists of samples
samples=("cellranger210_count_24933_5149STDY7274846_Danio_rerio_GRCz10_91" "cellranger210_count_24933_5149STDY7274847_Danio_rerio_GRCz10_91" "cellranger210_count_24933_5149STDY7274848_Danio_rerio_GRCz10_91")

#reference file used for variant calling - must match the reference file used in the initial alignment phase of scRNA-seq
#pipeline
reference_genome = "Danio_rerio.GRCz10.dna.toplevel.fa"

#Creates directory to eventually store the product of merging all cell-associated BAMs from a given experiment
mkdir out/merged

for SAMPLE in "${samples[@]}"
do
  	#separate experiment directory 
  	mkdir out/merged/${SAMPLE}
  
  	#merge all cell-associate BAMs from a given experiment in an output BAM
  	samtools merge out/merged/${SAMPLE}/sorted_QC_${SAMPLE}.bam out/${SAMPLE}/*.bam
	
	#indexing merged BAM - used for mpileup
	samtools index out/merged/${SAMPLE}/sorted_QC_${SAMPLE}.bam

	#use mpileup whilst ignoring read group - mpileup will usually consider transcripts on a sample by sample basis which
	#in our case is a cell by cell basis.
	#Therefor we want to treat all transcripts as if they are from the same sample for this phase of variant calling.
	bcftools mpileup --ignore-RG --skip-indels -f ${reference_genome} \
		out/merged/${SAMPLE}/sorted_QC_${SAMPLE}.bam \
	       	| bcftools call -mv -Oz -o out/merged/${SAMPLE}/unfiltered_calls_${SAMPLE}.vcf.gz
  	
	#filtering protocol
	bcftools filter -e"%QUAL<50 || (RPB<0.5 && %QUAL<70) || (AC<5 && %QUAL<70) || IDV <=3 || IDV/%MAX(DP)<=0.3" -Oz -o out/merged/${SAMPLE}/filtered_calls_${SAMPLE}.vcf.gz

	echo "Filtered Calls for ${SAMPLE} located at out/merged/${SAMPLE}/filtered_calls_${SAMPLE}.vcf.gz"
	echo "There are"
	bcftools view -H out/merged/${SAMPLE}/filtered_calls_${SAMPLE}.vcf.gz | wc -l
	echo "calls for ${SAMPLE}"
done
