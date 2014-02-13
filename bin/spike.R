tbl <- read.table("out/stat/spike.txt.gz", sep="\t", header=F)
tmp <- table(tbl[, c('V2', 'V3')])
write.table(tmp[, order(colSums(tmp), decreasing=T)], "out/stat/spike_density.txt", sep="\t", quote=F, col.names=NA, row.names=T)
write.table(table(tbl[, c('V2', 'V4')]), "out/stat/spike_position.txt", sep="\t", quote=F, col.names=NA, row.names=T)

plot(tmp[21, ], tmp[35, ], log='xy'); lines(c(1, 30000), c(1, 30000))
