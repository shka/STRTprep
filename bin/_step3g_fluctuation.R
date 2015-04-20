source('bin/_fluctuation.R')

draw_fluctuation <- function(fluctuation, sig=0.01) {
    ord <- order(fluctuation$score, decreasing=T)
    score <- fluctuation$score[ord]
    p.adj <- -log10(fluctuation$p.adj[ord])
    prop <- (1:length(score))/length(score)
    
    pos <- which(score>0)
    score.pos <- score[pos]
    p.adj.pos <- p.adj[pos]
    prop.pos <- prop[pos]
    
    isFinite <- which(is.finite(p.adj.pos))
    score.pos.isFinite <- score.pos[isFinite]
    p.adj.pos.isFinite <- p.adj.pos[isFinite]
    
    par(mar=c(4, 4, .5, 4), mgp=c(2.5, 1, 0))
    
    plot(range(score.pos), range(p.adj.pos.isFinite), log='x', col='white', xlab='Fluctuation score', ylab=expression(-log[10]('adjusted p-value')), axes=F)
    points(score.pos.isFinite, p.adj.pos.isFinite, cex=.25)
    axis(1)
    axis(2)
    
    par(new=T)
    
    plot(range(score.pos), range(prop.pos), log='x', col='white', xlab='', ylab='', axes=F)
    points(score.pos, prop.pos, cex=.25, col='red')
    axis(4, col='red')
    mtext(sprintf('Proportion (%d regions)', length(score)), side=4, line=2.5, col='red')
    
    tmp <- score[length(which(p.adj > -log10(sig)))]
    abline(v=tmp, col='blue', lty=3)
    mtext(sprintf(' score=%.2f', tmp), line=-1, at=tmp, adj=0, col='blue')
    mtext(sprintf(' adj.p=%s', sig), line=-2, at=tmp, adj=0, col='blue')
    mtext(sprintf(' %d regions', length(which(score>=tmp))), line=-3, at=tmp, adj=0, col='blue')
}

###

args <- commandArgs(trailingOnly=T)
p.fluctuation <- as.numeric(args[1])
dir <- args[2]

load(sprintf("%s/nreads.RData", dir))

errormodels <- estimate_errormodels(nreads)
save(errormodels, file=sprintf('%s/errormodels.RData', dir), compress='gzip')

errormodel <- merge_errormodels(errormodels)
save(errormodel, file=sprintf('%s/errormodel.RData', dir), compress='gzip')

fluctuation <- estimate_fluctuation(errormodel)
save(fluctuation, file=sprintf('%s/fluctuation.RData', dir), compress='gzip')

pdf(sprintf('%s/fig_fluctuation.pdf', dir), width=2.26, height=2.26, pointsize=8)
draw_fluctuation(fluctuation, sig=p.fluctuation)
dev.off()

tmp <- data.frame(Gene=names(fluctuation$p.adj), pvalue=fluctuation$p.adj)
gz <- gzfile(sprintf('%s/fluctuation.txt.gz', dir), 'w')
write.table(tmp, file=gz, quote=F, sep="\t", row.names=F, col.names=T)
close(gz)

###

row.fluctuated <- names(fluctuation$p.adj)[which(fluctuation$p.adj < p.fluctuation)]
row.spikes <- rownames(extract_spikein_reads(nreads))
nreads.fluctuated <- nreads[union(row.fluctuated, row.spikes), ]

###

library(Heatplus)
distfun <- function(x) as.dist((1-cor(t(x), method='spearman'))/2)
clustfun <- function(d) hclust(d, method='ward.D2')
tmp.nreads.pre <- nreads.fluctuated+min(nreads[which(nreads>0)])
tmp.nreads <- tmp.nreads.pre[setdiff(rownames(tmp.nreads.pre), rownames(extract_spikein_reads(tmp.nreads.pre))), ]
hm <- annHeatmap2(log10(tmp.nreads), scale='none', col=heat.colors, legend=2,
                  dendrogram=list(clustfun=clustfun, distfun=distfun, lwd=.5),
                  label=list(
                      Col=list(nrow=max(nchar(colnames(tmp.nreads)))/2.25),
                      Row=list(nrow=max(nchar(rownames(tmp.nreads)))/2.25)))
pdf(sprintf('%s/fig_heatmap_global.pdf', dir), width=6.69, height=6.69)
plot(hm)
dev.off()

###

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

tmp <- pca_by_spearman(t(tmp.nreads))
pdf(sprintf('%s/fig_pca.pdf', dir), width=6.69, height=6.69, pointsize=6)
plot(tmp$factor.score[, 1], tmp$factor.score[, 2], col='white', xlab='PC1', ylab='PC2')
text(tmp$factor.score[, 1], tmp$factor.score[, 2], labels=unlist(sapply(colnames(tmp.nreads), function(colname) { strsplit(colname, '\\|')[[1]][2] })))
dev.off()

###

### Local Variables:
### eval: (setq tmppath (file-name-directory (buffer-file-name)))
### eval: (add-to-list 'exec-path (concat tmppath "../.homebrew/bin"))
### eval: (setq default-directory (concat tmppath ".."))
### End:
