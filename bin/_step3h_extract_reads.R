args <- commandArgs(trailingOnly=T)
idx <- ifelse(is.na(args[1]), 0, as.numeric(args[1]))
path_samples <- ifelse(is.na(args[2]), 'out/byGene/samples.csv', args[2])
path_reads <- ifelse(is.na(args[3]), 'out/byGene/reads.RData', args[3])
path_output <- ifelse(is.na(args[4]), 'out/byGene/reads0.txt.gz', args[4])

samples <- read.table(path_samples, header=T, sep=',', quote='', check.names=F)
tmp.classes <- samples[, sprintf('CLASS.%d', idx)]
tmp.blocks <- samples [, sprintf('BLOCK.%d', idx)]
names(tmp.classes) <- names(tmp.blocks) <- sprintf('%s.%s|%s', samples[, 'LIBRARY'], samples[, 'WELL'], samples[, 'NAME'])
classes <- tmp.classes[which(!is.na(tmp.classes))]
blocks <- tmp.blocks[which(!is.na(tmp.classes))]

load(path_reads)
tmp.reads <- reads[, names(classes)]
tmp.reads2 <- tmp.reads[which(rowSums(tmp.reads) > 0), ]
reads <- tmp.reads2[order(rownames(tmp.reads2)), order(colnames(tmp.reads2))]

save(reads, file=path_output, compress='gzip')
