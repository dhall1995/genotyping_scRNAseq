#!/bin/bash

samples=("cellranger210_count_24933_5149STDY7274846_Danio_rerio_GRCz10_91" "cellranger210_count_24933_5149STDY7274847_Danio_rerio_GRCz10_91" "cellranger210_count_24933_5149STDY7274848_Danio_rerio_GRCz10_91")
reference_genome = "Danio_rerio.GRCz10.dna.toplevel.fa"

for SAMPLE in "${samples[@]}"
do
  samtools merge out/merged/${SAMPLE}/sorted_QC_${SAMPLE}.bam out/${SAMPLE}/*.bam

	samtools index out/merged/${SAMPLE}/sorted_QC_${SAMPLE}.bam

	bcftools mpileup --ignore-RG --skip-indels -f ${reference_genome} \
		out/merged/${SAMPLE}/sorted_QC_${SAMPLE}.bam \
	       	| bcftools call -mv -Oz -o out/merged/${SAMPLE}/unfiltered_calls_${SAMPLE}.vcf.gz
  bcftools filter -e"%QUAL<50 || (RPB<0.5 && %QUAL<70) || (AC<5 && %QUAL<70) || IDV <=3 || IDV/%MAX(DP)<=0.3" -Oz -o out/merged/${SAMPLE}/filtered_calls_${SAMPLE}.vcf.gz

	echo "Filtered Calls for ${SAMPLE} located at out/merged/${SAMPLE}/filtered_calls_${SAMPLE}.vcf.gz"
	echo "There are"
	bcftools view -H out/merged/${SAMPLE}/filtered_calls_${SAMPLE}.vcf.gz | wc -l
	echo "calls for ${SAMPLE}"
done
