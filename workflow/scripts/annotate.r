#!/usr/bin/Rscript

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ChIPseeker")

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 2){
  bedfile <- args[1]
  output <- args[2]
} else {
  stop("Requires [bedfile] [output]", call.=FALSE)
}

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peaks <- readPeakFile(bedfile)

annotations <- annotatePeak(peaks, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)

write.table(as.data.frame(annotations@anno), file = output, quote=FALSE, sep = "\t", row.names=FALSE)