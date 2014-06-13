library(MASS)

exclude_wells <- function(reads, wells) {
    if(wells != '') {
        reads.wells <- sapply(colnames(reads), function(n) { strsplit(n, '\\.')[[1]][2] })
        tmp <- reads[, which(!is.element(reads.wells, strsplit(wells, ',')[[1]]))]
        tmp[which(rowSums(tmp)>0), ]
    } else
        reads
}

extract_spike <- function(reads) {
    reads[which(substr(rownames(reads), 1, 10) == 'RNA_SPIKE_'), ]
}

