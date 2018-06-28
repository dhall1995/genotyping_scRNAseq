#!/bin/bash
set -e
samples=("cellranger210_count_24933_5149STDY7274847_Danio_rerio_GRCz10_91" "cellranger210_count_24933_5149STDY7274846_Danio_rerio_GRCz10_91" "cellranger210_count_24933_5149STDY7274848_Danio_rerio_GRCz10_91")
reference_genome = "../../Danio_rerio.GRCz10.dna.toplevel.fa"

cd out

#create storage directory for VCF outputs
mkdir all_vcfs
for file in "${samples[@]}"
do
  #within sample directory
  cd ${file}
  
  #indexing the separate cell-associated BAMs ready for variant calling
  for BAM in *.bam
  do
    samtools index ${BAM}
  done
  
  #now within the out directory
  cd ..

  echo ${file}
  #storage directory for each experiment
  mkdir all_vcfs/${file}
  
  #now within the sample directory
  cd ${file}

  #for loop runs over indexed files in case indexing of a particular BAM has failed
  for name in *.bai
  do
    #recover barcode bam file name
    newname=${name%.bai}
    echo ${newname}
    echo 'about to start mpileup'
    
    #variant calling using variant sites as reference. Save output in the form {barcode}.vcf.gz
    bcftools mpileup -E -Oz -o ${newname%-1.bam}.vcf.gz --skip-indels -R ../merged/${SAMPLE}/filtered_calls_${SAMPLE}.vcf.gz -f ${reference_genome} ${newname}

    #index the compressed mpileup output ready for variant calling
    bcftools index ${newname%-1.bam}.vcf.gz

    #variant calling with ouput to barcode.vcf.gz. Script assumes barcodes are of the form e.g AAACCTGAGAAGCCCA-1
    bcftools call -m -R ../merged/${SAMPLE}/filtered_calls_${SAMPLE}.vcf.gz \
    -Oz -o ${newname%-1.bam}-1.vcf.gz ${newname%-1.bam}.vcf.gz

    #remove indexed vcf file
    rm ${newname%.bam}.vcf.gz*
    
    #note down those BAM files which have have variant calling completed successfully.
    mv ${newname} ${newname}.done
  done

  #check if there are any non-completed BAM files - i.e. cells for which variant calling failed. Note them down
  #in a csv
  for second_file in *
  do
    if [[ ${second_file} == *.bam ]]
    then
      echo ${second_file%.bam} >> uncalled_barcodes_${index}.tsv
    fi
  done

  #return completed BAM files back to normal file names
  for called_barcode in *.done
  do
    mv ${called_barcode} ${called_barcode%.done}
  done

  #now in out folder
  cd ..
  
  #move completed cell VCFs to the relevant output folder
  mv ${file}/*.vcf.gz all_vcfs/${file}/
done

#now within VCF output file
cd all_vcfs

for file in "${samples[@]}"
do
  #within experiment specific VCF output file
  cd ${file}
  #within each VCF the sample name will now be the same sample name as the entire experiment. For each barcode VCF 
  #we therefore change the sample name to agree with the barcode as we are now treating each cell as a separate sample
  for vcf_file in *.vcf
  do
    bcftools index ${vcf_file}
    gunzip ${vcf_file}
  	sed -i "s/${file}/${vcf_file%.vcf.gz}/g" ${vcf_file%.gz}
  	bcftools view -Oz -o ${vcf_file} ${file%.gz}
  	rm ${vcf_file%.gz}
  done
  
  #merge all barcode VCFs from a given experiment together
  bcftools merge -Oz -o merged_calls_${file}.vcf.gz
  
  #filter out variants present in <5% of cells in that experiment. Output to filtered_merged_calls_[experiment].vcf.gz
  #This filtered vcf is the output file to 
  bcftools filter -Oz -o filtered_merged_calls_${file}.vcf.gz -e "F_MISSING > 0.95" merged_calls_${file}.vcf.gz
  
  #return to VCF output file
  cd ..
  echo "filtered merged calls for ${file} are located at out/all_vcfs/${file}/filtered_merged_calls_${file}.vcf.gz"
done

cd ../..
