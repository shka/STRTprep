reads <- read.table("out/byGene/reads.txt.gz", header=T, sep="\t", quote='', row.names=1, check.names=F)

reads.spike <- reads[which(substr(rownames(reads), 1, 10) == 'RNA_SPIKE_'), ]
reads.spike.colSums <- colSums(reads.spike)
reads.scale <- 1/reads.spike.colSums
nreads <- as.matrix(reads*rep(reads.scale, each=nrow(reads)))

save(nreads, file="out/byGene/nreads.RData", compress='gzip')

gz <- gzfile("out/byGene/nreads.txt.gz", 'w')
write.table(nreads, file=gz, quote=F, sep="\t", row.names=T, col.names=NA)
close(gz)

###

### Local Variables:
### eval: (setq tmppath (file-name-directory (buffer-file-name)))
### eval: (add-to-list 'exec-path (concat tmppath "../.homebrew/bin"))
### eval: (setq default-directory (concat tmppath ".."))
### End:
