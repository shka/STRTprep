#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
idx <- ifelse(is.na(args[1]), 0, as.numeric(args[1]))
path_samples <- ifelse(is.na(args[2]), 'out/byGene/samples.csv', args[2])
path_reads <- ifelse(is.na(args[3]), 'out/byGene/reads.RData', args[3])
path_nreads <- ifelse(is.na(args[4]), 'out/byGene/nreads.RData', args[4])
path_output <- ifelse(is.na(args[5]), 'out/byGene/sources0.RData', args[5])

samples <- read.table(path_samples, header=T, sep=',', quote='', check.names=F)
tmp.classes <- samples[, sprintf('CLASS.%d', idx)]
tmp.blocks <- samples [, sprintf('BLOCK.%d', idx)]
names(tmp.classes) <- names(tmp.blocks) <- sprintf('%s.%s|%s', samples[, 'LIBRARY'], samples[, 'WELL'], samples[, 'NAME'])
classes <- tmp.classes[which(!is.na(tmp.classes))]
blocks <- tmp.blocks[which(!is.na(tmp.classes))]

load(path_reads)
tmp.reads <- reads[, names(classes)]
tmp.reads2 <- tmp.reads[which(rowSums(tmp.reads) > 0), ]

if(regexpr('/byTFE/', path_reads)) {
  tmp <- apply(tmp.reads2, 1, function(row) { length(which(row>0)) > 2 })
  tmp.reads3 <- tmp.reads2[which(tmp), ]
  reads <- tmp.reads3[order(rownames(tmp.reads3)), order(colnames(tmp.reads3))]
} else
  reads <- tmp.reads2[order(rownames(tmp.reads2)), order(colnames(tmp.reads2))]

load(path_nreads)
nreads <- nreads[rownames(reads), colnames(reads)]

save(reads, nreads, classes, blocks, file=path_output, compress='gzip')
