#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0){
  stop("Requires [bedfile]", call.=FALSE)
} else {
  bedfile = args[1]
}

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peaks <- readPeakFile(bedfile)

annotations <- annotatePeak(peaks, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)

write.table(annotations@anno, file = "", quote=FALSE, sep = "\t", row.names=FALSE)