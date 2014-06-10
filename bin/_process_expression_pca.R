library(maptools)

args <- commandArgs(trailingOnly=T)
nreads.uniq.success.fluctuated.file <- args[1]
lib <- sub('out/exp/', '', strsplit(nreads.uniq.success.fluctuated.file, '\\.')[[1]][1])
nreads.uniq.success.file <- sub('.fluctuated', '', nreads.uniq.success.fluctuated.file)
pca.uniq.success.fluctuated.file <- sub('nreads.', 'pca.', nreads.uniq.success.fluctuated.file)

process_prcomp <- function(dat) {
    n <- nrow(dat)
    r <- prcomp(dat)
    pcs <- sum(r$sdev^2 >= 1)
    r$loadings <- t(r$sdev*t(r$rotation))[, 1:pcs, drop=F]
    r$contributions <- rowSums(r$loadings^2)
    r$eigen_values <- r$sdev[1:pcs]^2
    r$denominator <- sum(r$sdev^2)
    r$proportions <- r$eigen_values/r$denominator*100
    r$cumulative_proportions <- cumsum(r$proportions)
    r$scores <- r$x * sqrt(n/(n-1))
    invisible(r)
}

load(nreads.uniq.success.fluctuated.file)
load(nreads.uniq.success.file)
nreads.uniq.success.colOffsets <- apply(nreads.uniq.success[, colnames(nreads.uniq.success.fluctuated)], 2, function(cols) { min(cols[which(cols > 0)])/2 })
pca.uniq.success.fluctuated <- process_prcomp(t(log10(nreads.uniq.success.fluctuated+runif(length(nreads.uniq.success.fluctuated), min=0, max=rep(nreads.uniq.success.colOffsets, each=nrow(nreads.uniq.success.fluctuated))))))

gz <- gzfile(pca.uniq.success.fluctuated.file, 'wb')
save(pca.uniq.success.fluctuated, file=gz)
close(gz)

plot_pca <- function(pcx, pcy) {
    x <- pca.uniq.success.fluctuated$scores[, pcx]
    y <- pca.uniq.success.fluctuated$scores[, pcy]
    plot(x, y, xlab=sprintf('PC%d', pcx), ylab=sprintf('PC%d', pcy), col='gray', cex=1/3)
    pointLabel(x, y, labels=sub(sprintf('%s.', lib), '', rownames(pca.uniq.success.fluctuated$scores)), adj=c(.5, .5))
}

pdf(sub('.RData.gz', '.pdf', pca.uniq.success.fluctuated.file), width=6.69/2, height=6.69/2, pointsize=6)
layout(matrix(c(1, 4, 2, 3), 2, 2, byrow=T), c(1, 1), c(1, 1), T)
par(mar=c(3, 3, 1, 1), mgp=c(2, 1, 0)); plot_pca(1, 2)
par(mar=c(3, 3, 1, 1), mgp=c(2, 1, 0)); plot_pca(1, 3)
par(mar=c(3, 3, 1, 1), mgp=c(2, 1, 0)); plot_pca(2, 3)
par(mar=c(3, 3, 1, 1), mgp=c(2, 1, 0))
plot(0, 0, col='white', axes=F, xlab='', ylab='')
text(0, 0, lib)
dev.off()
