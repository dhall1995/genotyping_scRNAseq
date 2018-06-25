#!/usr/bin/env Rscript

library("optparse")

mywd <- getwd()
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="merged vcf file", metavar="character"),
	make_option(c("-G", "--Gene-build"), type="character", default="GRCz10",
              help="Gene Build", metavar="character")
  make_option(c("-N", "--Num-donors"), type="numeric", default ="4",
              help="Number of donors", metavar="numeric")
  make_option(c("-o", "--Output-file"), type="character", default = mywd,
              help="Output directory", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

library(VariantAnnotation)
library(snpStats)
library(pcaMethods)
library(mclust)

vcf <- readVcf(opt$file, opt$Gene-build)

#only pick out single nucleotide variations
vcf <- vcf[isSNV(vcf)]
#create snp matrix
snpmatrix <- genotypeToSnpMatrix(vcf)

#coerce genotype matrix into numeric form, center and normalise it then run ppca
p <- as(snpmatrix$genotypes , “numeric”)
p <- prep(p, scale = “none”, center = TRUE)
resPPCA <- pca(p, method = “ppca”, center = TRUE, nPcs = 3)
ppca_pc1_pc2 <- scores(resPPCA)[,1:2]

#collect data on names, statistics for cells and statistics for snps
names <- rownames(snpmatrix$genotypes)
cell_stats <- row.summary("snpmatrix$genotypes")
snp_stats <- col.summary("snpmatrix$genotypes")

#save histograms for call rate for cells and snps
cell_stats_plot <- hist((1-cell_stats$Call.rate)*100, main = "Percentage of Missing SNPs per Cell", xlab = "Percent Missing (%)", ylab = "Frequency")
dev.copy(png, paste(opt$Output-file,"/cell_stats_plot.png")
dev.off()
snp_stats_plot <- hist((1-snp_stats$Call.rate)*100, main = "Percentage of Missing Cells per SNP", xlab = "Percent Missing (%)", ylab = "Frequency")
dev.copy(png, paste(opt$Output-file,"/snp_stats_plot.png")
dev.off()

#perform optimal clustering and save the BIC result and classifcation plot
df <- data.frame(sample.id = names, PC1 = ppca_pc1_pc2[,1], PC2 = ppca_pc1_pc2[,2])
mclust_pc1_pc2 <- Mclust(ppca_pc1_pc2)
classification_plot_optimal <- plot.Mclust(mclust_pc1_pc2)[2]
dev.copy(png, paste(opt$Output-file,"/optimal_clustering.png")
dev.off()
BIC_optimal <- plot.Mclust(mclust_pc1_pc2)[1]
dev.copy(png, paste(opt$Output-file,"/Mclust_BIC_results.png",)
dev.off()

#perform clustering with specified number of donors and save the plot
mclust_pc1_pc2_component <- Mclust(ppca_pc1_pc2, G=opt$Num-donors)
classification_plot_component <- plot.Mclust(mclust_pc1_pc2_component)[2]
dev.copy(png, paste(opt$Output-file, "/", as.character(opt$Num-donors),"_component_clustering.png")
dev.off()

write.csv(mclust_pc1_pc2_component$classification, paste(opt$Output-file,"/clusters.csv"))
