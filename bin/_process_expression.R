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

normalize <- function(reads, wells.exclude) {
    reads.spike <- extract_spike(reads)
    reads.spike.colSums <- colSums(reads.spike)
    reads.spike.wells <- sapply(colnames(reads.spike), function(n) { strsplit(n, '\\.')[[1]][2] })
    reads.scale <- min(reads.spike.colSums[which(!is.element(reads.spike.wells, strsplit(wells.exclude, ',')[[1]]))])/reads.spike.colSums
    reads*rep(reads.scale, each=nrow(reads))
}

rowCV2s <- function(reads, reads.rowMeans=NA) {
    if(!is.vector(reads.rowMeans)) reads.rowMeans <- rowMeans(reads)
    reads.rowVars <- apply(reads, 1, var)
    reads.rowVars/(reads.rowMeans^2)
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
rm(reads.uniq.success)

nreads.uniq.success.spike <- exclude_wells(extract_spike(nreads.uniq.success), wells.exception)
nreads.uniq.success.spike.rowMeans <- rowMeans(nreads.uniq.success.spike)
nreads.uniq.success.spike.rowCV2s <-
    rowCV2s(nreads.uniq.success.spike, nreads.uniq.success.spike.rowMeans)

errormodel.uniq.success <- glm(nreads.uniq.success.spike.rowCV2s ~ I(1/nreads.uniq.success.spike.rowMeans), family=Gamma(link='identity'))
gz <- gzfile(sub('.nreads', '.errormodel', nreads.uniq.success.file), 'wb')
save(errormodel.uniq.success, file=gz)
close(gz)

nreads <- exclude_wells(nreads.uniq.success, wells.exception)
nreads.rowMeans <- rowMeans(nreads)
nreads.rowCV2s <- rowCV2s(nreads, nreads.rowMeans)

fit.x <- exp(seq(log(min(nreads.rowMeans)), log(max(nreads.rowMeans)), length.out=100))
fit.y <- errormodel.uniq.success$coefficients[2]/fit.x+errormodel.uniq.success$coefficients[1]
fit.a <- gamma.shape(errormodel.uniq.success, it.lim=1000)$alpha
fit.s <- fit.y/fit.a
fit.y.upper <- qgamma(.99, shape=fit.a, scale=fit.s)
fit.y.lower <- qgamma(.01, shape=fit.a, scale=fit.s)
nreads.y <- errormodel.uniq.success$coefficients[2]/nreads.rowMeans+errormodel.uniq.success$coefficients[1]
nreads.s <- nreads.y/fit.a
nreads.p.adj <- p.adjust(1-pgamma(nreads.rowCV2s, shape=fit.a, scale=nreads.s), 'fdr')
nreads.isFluctuated <- which(nreads.p.adj<0.01)

nreads.uniq.success.fluctuated <- nreads[nreads.isFluctuated, ]
nreads.uniq.success.fluctuated.file <- sub('.success', '.success.fluctuated', nreads.uniq.success.file)
gz <- gzfile(nreads.uniq.success.fluctuated.file, 'wb')
save(nreads.uniq.success.fluctuated, file=gz)
close(gz)

cex <- .5
options(scipen=3)
tryCatch ({
    png(sub('.RData.gz', '.png', nreads.uniq.success.fluctuated.file),
        width=3.34, height=3.34, units='in', res=600, pointsize=8)
    par(mar=c(3.2, 3.2, 1, 1), mgp=c(1.9, 0.8, 0))
    plot(nreads.rowMeans, nreads.rowCV2s, cex=cex, col='gray', log='xy',
         xlab='Normalized spike-in read count', ylab=expression(CV^2))
    lines(fit.x, fit.y, col='blue')
    lines(fit.x, fit.y.upper, col='blue', lty=2)
    lines(fit.x, fit.y.lower, col='blue', lty=2)
    points(nreads.rowMeans[nreads.isFluctuated],
           nreads.rowCV2s[nreads.isFluctuated], cex=cex, col='magenta')
    points(nreads.uniq.success.spike.rowMeans,
           nreads.uniq.success.spike.rowCV2s, cex=cex, col='black')
    legend(x=min(nreads.rowMeans), y=min(nreads.rowCV2s), yjust=0,
           legend=c(
               sprintf('%d detected spike-ins, %d samples',
                       length(nreads.uniq.success.spike.rowMeans),
                       ncol(nreads.uniq.success.spike)),
               sprintf("Gamma fit; %1.3f+%1.3f/x, AIC:%1.3f",
                       errormodel.uniq.success$coefficients[1],
                       errormodel.uniq.success$coefficients[2],
                       errormodel.uniq.success$aic),
               'Gamma fit; 99% confidence',
               sprintf('%d detected regions', length(nreads.rowMeans)),
               sprintf('%d fluctuated regions; FDR<0.01',
                       length(nreads.isFluctuated))),
           lty=c(NA, 1, 2, NA, NA), pch=c(1, NA, NA, 1, 1),
           col=c('black', 'blue', 'blue', 'gray', 'red'), cex=6/8)
}, finally={ dev.off() })
options(scipen=0)
