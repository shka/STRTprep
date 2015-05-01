args <- commandArgs(trailingOnly=T)
readsfile <- args[1]
nreadsfile <- args[2]

load(readsfile)

reads.spike <- reads[which(substr(rownames(reads), 1, 10) == 'RNA_SPIKE_'), ]
reads.spike.colSums <- colSums(reads.spike)
reads.scale <- 1/reads.spike.colSums
nreads <- as.matrix(reads*rep(reads.scale, each=nrow(reads)))

gz <- gzfile(nreadsfile, 'w')
write.table(nreads, file=gz, quote=F, sep="\t", row.names=T, col.names=NA)
close(gz)

###

### Local Variables:
### eval: (setq tmppath (file-name-directory (buffer-file-name)))
### eval: (add-to-list 'exec-path (concat tmppath "../.homebrew/bin"))
### eval: (setq default-directory (concat tmppath ".."))
### End:
