source('bin/_library_expression.R')

normalize <- function(reads, wells.exclude) {
    reads.spike <- extract_spike(reads)
    reads.spike.colSums <- colSums(reads.spike)
    reads.spike.wells <- sapply(colnames(reads.spike), function(n) { strsplit(n, '\\.')[[1]][2] })
    reads.scale <- min(reads.spike.colSums[which(!is.element(reads.spike.wells, strsplit(wells.exclude, ',')[[1]]))])/reads.spike.colSums
    reads*rep(reads.scale, each=nrow(reads))
}

args <- commandArgs(trailingOnly=T)
reads.uniq.file <- args[1]
wells.failure <- args[2]
wells.exception <- args[3]
    
reads.uniq.success <- exclude_wells(as.matrix(read.table(reads.uniq.file, header=T, sep="\t", row.names=1)), wells.failure)
reads.uniq.success.file <- sub('.txt', '.RData', sub('.uniq', '.uniq.success', reads.uniq.file))
gz <- gzfile(reads.uniq.success.file, 'wb')
save(reads.uniq.success, file=gz)
close(gz)

nreads.uniq.success <- normalize(reads.uniq.success, wells.exception)
nreads.uniq.success.file <- sub('.reads', '.nreads', reads.uniq.success.file)
gz <- gzfile(nreads.uniq.success.file, 'wb')
save(nreads.uniq.success, file=gz)
close(gz)
