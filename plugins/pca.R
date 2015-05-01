#!/usr/bin/env Rscript

source('Rlib/STRTprepGateway.R')
gw <- STRTprepGateway$new('pca')

ex <- gw$getExpressions()
nreads <- ex$getNormalizedReads()
nreads.fluctuated <- ex$significant()$getNormalizedReads()

prefix <- gw$outputPrefix
if(nrow(nreads.fluctuated) == 0) {
  message("No significant ones.")
  system(sprintf("rm -f %s.pdf", prefix))
  quit(save='no')
}

##

pca_by_spearman <- function(dat) {
    nr <- nrow(dat)
    nc <- ncol(dat)
    vname <- colnames(dat)
    heikin <- colMeans(dat)
    bunsan <- apply(dat, 2, var)
    sd <- sqrt(bunsan)
    r <- cor(dat, method='spearman') # use='pairwise.complete.obs' - too slow
    result <- eigen(r)
    eva <- result$values
    eve <- result$vectors
    contr <- eva/nc*100
    cum.contr <- cumsum(contr)
    fl <- t(sqrt(eva)*t(eve))
    fs <- scale(dat)%*%eve*sqrt(nr/(nr-1))
    names(heikin) <- names(bunsan) <- names(sd) <- rownames(r) <- colnames(r) <- rownames(fl) <- colnames(dat)
    names(eva) <- names(contr) <- names(cum.contr) <- colnames(fl) <- colnames(fs) <- paste("PC", 1:nc, sep="")
    return(structure(list(mean=heikin, variance=bunsan, standard.deviation=sd, r=r, factor.loadings=fl, eva=eva, contribution=contr, cum.contribution=cum.contr, factor.scores=fs), class="pca"))
}

tmp.nreads <- nreads.fluctuated+min(nreads[which(nreads>0)])
pca <- pca_by_spearman(t(tmp.nreads))

pdf(sprintf('%s.pdf', prefix), width=6.69, height=6.69)
plot(pca$factor.score[, 1], pca$factor.score[, 2], col='white', xlab='PC1', ylab='PC2')
text(pca$factor.score[, 1], pca$factor.score[, 2], labels=colnames(tmp.nreads))
dev.off()

save(pca, file=sprintf('%s.RData', prefix), compress='gzip')
