Selection of scripts to genotype single cells within a BAM file created from a scRNA-seq run using barcoded cells
where the donor.id of each barcode is unknown. For a given number of experiments the scripts assume that for each experiment
we have a directory containing:
1) a BAM file named possorted_genome_bam.bam
2) a .tsv file containing a list of barcodes of interest

For my project at the sanger I had 3 cellranger runs: 
    "cellranger210_count_24933_5149STDY7274846_Danio_rerio_GRCz10_91"
    "cellranger210_count_24933_5149STDY7274847_Danio_rerio_GRCz10_91" 
    "cellranger210_count_24933_5149STDY7274848_Danio_rerio_GRCz10_91"

I therefore had three directories, one for each cellranger run. Within each directory I had:

1) a .bam file named possorted_genome_bam.bam

2) a list of all barcodes of interest called

e.g. cellranger210_count_24933_5149STDY7274846_Danio_rerio_GRCz10_91_barcodes.tsv
    
All shell scripts are intended to be run in the directory directly above the experiment directories. 

- To separate each BAM file into one BAM file for each barcode of interest, we first run ./Extract_barcodes.sh
- We then run Variant_Calling_per_Experiment.sh. This merges the barcode BAMs for each experiment into a QC-ed BAM
  file for each experiment and then performs de novo variant calling on this merged BAM and filters out poor quality SNPs
- We then run variant_calling_per_barcode.sh. This Uses our vcf file generated for each experiment as a reference file in
  order to genotype each barcode BAM at each variant site. We then filter out SNPs which are missing in too many cells 
  and finally produce a compressed vcf for each experiment in which each barcoded cell BAM file is treated as a different
  sample.
  
  
  
The R script is currently untested as I have only run those commands directly from within the rstudio terminal. However, 
if the script works it should output some descriptive plots about the sparsity of the SNP matrix generated from each experiment
and also attempt to run probabilistic pca on the SNP matrix, assigning each barcode to a likely cluster. The number of donors
can be input directly into that script and the ppca algorithm will try and come up with the most likely partitioning of the 
cells into donors given the SNP profiles of each cell. I have only tested this so far with small numbers of donors and I expect
that for n>4 donors the clustering may become unreliable. 


The scripts require a quite large amount of storage capacity to work since it generates a lot of BAM files. If storage is low
then it is safe to delete the possorted_genome_bam.bam files after running Extract_barcodes.sh. Similarly, if they are unneeded
then it is safe to delete any merged BAM ***after*** variant calling has been performed. 
    
