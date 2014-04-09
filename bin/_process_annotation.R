args <- commandArgs(trailingOnly=T)

load(args[1])
ann <- read.table(args[2], header=T, sep="\t", row.names=1, quote='')
ann.fluctuated <- ann[rownames(reads.uniq.detected.normalized.fluctuated), ]
tmp <- cbind(ann.fluctuated, reads.uniq.detected.normalized.fluctuated)

gz <- gzfile(args[3], 'w')
write.table(tmp, file=gz, quote=F, sep="\t", col.names=NA)
close(gz)
