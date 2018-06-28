#!/bin/bash
set -e

#Directories in which the samples BAM files are located. Note that BAM files are assumed to be named
#possorted_genome_bam.bam within each of these directories. In my case, each directory held results fomr one Cell Ranger run
#and was associated with a single experiment consisting of 3 donors each. 
samples=("cellranger210_count_24933_5149STDY7274846_Danio_rerio_GRCz10_91" "cellranger210_count_24933_5149STDY7274847_Danio_rerio_GRCz10_91" "cellranger210_count_24933_5149STDY7274848_Danio_rerio_GRCz10_91")

#script assumes there is no dirrectory named out within the current directory. Results are placed in the out folder once 
#it is created.
mkdir out
#Creates directory to eventually store the product of merging all cell-associated BAMs from a given experiment
mkdir out/merged

touch mark1

#Barcode names are held in a temporary location and processed 500 at a time in order to avoid too many files being open at once.
#This does slow down the code somewhat for very large lists of BAMs but does avoid crashes
mkdir temp_barcodes


#for SAMPLE in "${samples[@]}"
#do
#echo ${SAMPLE}
#done

for SAMPLE in "${samples[@]}"
do

#Create output file for a given sample to store cell-associated BAM files for each cell in that sample
mkdir out/${SAMPLE}
  #Isolate the experiment number - this section is mostly cosmetic and not necessarily needed for arbitrary sample names.
  #in this case it tallies with how my barcode lists were named but this process could be streamlined for arbitrary sample names
  #and barcode list file names
  experiment = ${SAMPLE#cellranger210_count_24933_5149STDY72748}
  experiment = ${experiment%_Danio_rerio_GRCz10_91}
  
  #Split barcode list into 500 barcode chunks - store in temp_barcodes
  split -l 500 ${SAMPLE}/cell_ranger_barcodes${experiment}.tsv temp_barcodes/barcodes_


  for BARCODES in temp_barcodes/barcodes_*
  do
    echo ${BARCODES}
    #actual splitting script
    samtools view -@ 4 -h $SAMPLE/possorted_genome_bam.bam | perl -nle 'use strict; use autodie; our %h; BEGIN{open(my$fh,q(<),shift@ARGV); my$od=shift@ARGV; $od//=q(); while(<$fh>){chomp; open(my$f2,"| samtools view -u - |bamstreamingmarkduplicates level=0 tag=CB | samtools view -b - > $od/$_.bam");$h{$_}=$f2; }close $fh}  if(/^@/){foreach my$fh (values %h){print {$fh} $_ }}elsif(m{\tCB:Z:(\S+)\b}){ my$fh=$h{$1}||next; print {$fh} $_;} END{close $_ foreach values %h; warn "closed BAMs\n"}' ${BARCODES} out/${SAMPLE}
  done

#clear temp_barcodes file ready fo the next sample
rm temp_barcodes/*
touch mark2${SAMPLE}


done

touch mark3

rm -r temp_barcodes

echo "Barcode BAMs created"
