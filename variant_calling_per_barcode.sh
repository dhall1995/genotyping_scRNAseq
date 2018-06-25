#!/bin/bash
set -e
samples=("cellranger210_count_24933_5149STDY7274847_Danio_rerio_GRCz10_91" "cellranger210_count_24933_5149STDY7274846_Danio_rerio_GRCz10_91" "cellranger210_count_24933_5149STDY7274848_Danio_rerio_GRCz10_91")
reference_genome = "../../Danio_rerio.GRCz10.dna.toplevel.fa"

cd out
for file in "${samples[@]}"
do
  cd ${file}
  for BAM in *.bam
  do
    samtools index ${BAM}
  done
  cd ..

  echo ${file}
  mkdir all_vcfs/${file}
  cd ${file}

  for name in *.bai
  do
    newname=${name%.bai}
    echo ${newname}
    echo 'about to start mpileup'
    bcftools mpileup -E -Oz -o ${newname%-1.bam}.vcf.gz --skip-indels -R ../merged/${SAMPLE}/filtered_calls_${SAMPLE}.vcf.gz -f ${reference_genome} ${newname}

    bcftools index ${newname%-1.bam}.vcf.gz

    bcftools call -m -R ../merged/${SAMPLE}/filtered_calls_${SAMPLE}.vcf.gz \
    -Oz -o ${newname%1.bam}${index}.vcf.gz ${newname%-1.bam}.vcf.gz

    rm ${newname%-1.bam}.vcf.gz*
    mv ${newname} ${newname}.done
  done

  for second_file in *
  do
    if [[ ${second_file} == *.bam ]]
    then
      echo ${second_file%.bam} >> uncalled_barcodes_${index}.tsv
    fi
  done

  for called_barcode in *.done
  do
    mv ${called_barcode} ${called_barcode%.done}
  done

  cd ..
  mv ${file}/*.vcf.gz all_vcfs/${file}/
done

cd all_vcfs
for file in "${samples[@]}"
do
  cd ${file}
  for vcf_file in *.vcf
  do
    bcftools index ${vcf_file}
    gunzip ${vcf_file}
  	sed -i "s/${file}/${vcf_file%.vcf.gz}/g" ${vcf_file%.gz}
  	bcftools view -Oz -o ${vcf_file} ${file%.gz}
  	rm ${vcf_file%.gz}
  done
  
  bcftools merge -Oz -o merged_calls_${file}.vcf.gz
  bcftools filter -Oz -o filtered_merged_calls_${file}.vcf.gz -e "F_MISSING > 0.95" merged_calls_${file}.vcf.gz
  cd ..
  echo "filtered merged calls for ${file} are located at out/all_vcfs/${file}/filtered_merged_calls_${file}.vcf.gz"
done

cd ../..
