#!/bin/bash
set -e

samples=("cellranger210_count_24933_5149STDY7274846_Danio_rerio_GRCz10_91" "cellranger210_count_24933_5149STDY7274847_Danio_rerio_GRCz10_91" "cellranger210_count_24933_5149STDY7274848_Danio_rerio_GRCz10_91")

mkdir out
mkdir out/merged

touch mark1

mkdir temp_barcodes

for SAMPLE in "${samples[@]}"
do
echo ${SAMPLE}
done

for SAMPLE in "${samples[@]}"
do
mkdir out/${SAMPLE}
  experiment = ${SAMPLE#cellranger210_count_24933_5149STDY72748}
  experiment = ${experiment%_Danio_rerio_GRCz10_91}
  split -l 500 ${SAMPLE}/cell_ranger_barcodes${experiment}.tsv temp_barcodes/barcodes_


  for BARCODES in temp_barcodes/barcodes_*
  do
    echo ${BARCODES}
    samtools view -@ 4 -h $SAMPLE/possorted_genome_bam.bam | perl -nle 'use strict; use autodie; our %h; BEGIN{open(my$fh,q(<),shift@ARGV); my$od=shift@ARGV; $od//=q(); while(<$fh>){chomp; open(my$f2,"| samtools view -u - |bamstreamingmarkduplicates level=0 tag=CB | samtools view -b - > $od/$_.bam");$h{$_}=$f2; }close $fh}  if(/^@/){foreach my$fh (values %h){print {$fh} $_ }}elsif(m{\tCB:Z:(\S+)\b}){ my$fh=$h{$1}||next; print {$fh} $_;} END{close $_ foreach values %h; warn "closed BAMs\n"}' ${BARCODES} out/${SAMPLE}
  done

rm temp_barcodes/*
touch mark2${SAMPLE}


done

touch mark3

rm -r temp_barcodes

echo "Barcode BAMs created"
