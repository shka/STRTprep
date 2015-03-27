library(MASS)

parse_libid <- function(colname) unlist(strsplit(colname, '\\.'))[1]

rowCV2s <- function(reads) {
    reads.rowVars <- apply(reads, 1, var)
    reads.rowVars/(rowMeans(reads)^2)
}

extract_spikein_reads <- function(reads)
    reads[which(substr(rownames(reads), 1, 10) == 'RNA_SPIKE_'), ]

estimate_errormodels <- function(nreads) {
    libids <- unique(sapply(colnames(nreads), parse_libid))
    tmp <- lapply(libids, function(libid) {
        cat(sprintf('estimate_errormodels : %s\n', libid))
        targets <- sapply(colnames(nreads),
                          function(colname) { parse_libid(colname) == libid })
        nreads.tmp <- nreads[, which(targets)]
        nreads.spike <- extract_spikein_reads(nreads.tmp)
        x.spike <- 1/rowMeans(nreads.spike)
        y.spike <- rowCV2s(nreads.spike)
        model <- glm(y.spike ~ x.spike, family=Gamma(link='identity'))
        shape <- gamma.shape(model, it.lim=1000)$alpha
        scales <- (model$coefficients[1]+model$coefficients[2]/rowMeans(nreads.tmp))/shape
        responses <- rowCV2s(nreads.tmp)
        list(model=model, shape=shape, scales=scales, responses=responses)
    })
    names(tmp) <- libids
    tmp
}

create_scale_matrix <- function(errormodels) {
    tmp <- data.frame(errormodels[[1]]$scales)
    for(i in 2:length(errormodels))
        tmp <- cbind(tmp, data.frame(errormodels[[i]]$scales))
    colnames(tmp) <- names(errormodels)
    as.matrix(tmp)
}

create_response_matrix <- function(errormodels) {
    tmp <- data.frame(errormodels[[1]]$responses)
    for(i in 2:length(errormodels))
        tmp <- cbind(tmp, data.frame(errormodels[[i]]$responses))
    colnames(tmp) <- names(errormodels)
    as.matrix(tmp)
}

create_shape_vector <- function(errormodels) {
    tmp <- errormodels[[1]]$shape
    for(i in 2:length(errormodels))
        tmp <- c(tmp, errormodels[[i]]$shape)
    names(tmp) <- names(errormodels)
    tmp
}

merge_errormodels <- function(errormodels) {
    tmp.scales <- create_scale_matrix(errormodels)
    scales <- tmp.scales[, 1]
    tmp.responses <- create_response_matrix(errormodels)
    tmp.responses.scaled <- tmp.responses*(rep(scales, ncol(tmp.scales))/tmp.scales)
    responses <- rowSums(tmp.responses.scaled)
    tmp.shapes <- create_shape_vector(errormodels)
    shape <- sum(tmp.shapes)
    list(responses=responses, shape=shape, scales=scales)
}

estimate_fluctuation <- function(errormodel) {
    p.adj <- p.adjust(pgamma(q=errormodel$responses, shape=errormodel$shape, scale=errormodel$scale, lower.tail=F), 'fdr')
    score <- errormodel$responses/(errormodel$scale*errormodel$shape)
    list(p.adj=p.adj, score=score)
}

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

load("out/cg/nreads.RData")

errormodels <- estimate_errormodels(nreads)
save(errormodels, file='out/cg/errormodels.RData', compress='gzip')

errormodel <- merge_errormodels(errormodels)
save(errormodel, file='out/cg/errormodel.RData', compress='gzip')

fluctuation <- estimate_fluctuation(errormodel)
save(fluctuation, file='out/cg/fluctuation.RData', compress='gzip')

pdf('out/cg/fig_fluctuation.pdf', width=2.26, height=2.26, pointsize=8)
draw_fluctuation(fluctuation, sig=.25)
dev.off()

tmp <- data.frame(Gene=names(fluctuation$p.adj), pvalue=fluctuation$p.adj)
gz <- gzfile('out/cg/fluctuation.txt.gz', 'w')
write.table(tmp, file=gz, quote=F, sep="\t", row.names=F, col.names=T)
close(gz)

row.fluctuated <- names(fluctuation$p.adj)[which(fluctuation$p.adj < 0.25)]
row.spikes <- rownames(extract_spikein_reads(nreads))
nreads.fluctuated <- nreads[union(row.fluctuated, row.spikes), ]
save(nreads.fluctuated, file='out/cg/nreads.fluctuated.RData', compress='gzip')

###

library(Heatplus)
distfun <- function(x) as.dist((1-cor(t(x), method='spearman'))/2)
clustfun <- function(d) hclust(d, method='ward.D2')
tmp.nreads.pre <- nreads.fluctuated+min(nreads.fluctuated[which(nreads.fluctuated>0)])
tmp.nreads <- tmp.nreads.pre[setdiff(rownames(tmp.nreads.pre), rownames(extract_spikein_reads(tmp.nreads.pre))), ]
hm <- annHeatmap2(log10(tmp.nreads), scale='none', dendrogram=list(clustfun=clustfun, distfun=distfun, lwd=.5), col=function(n) g2r.colors(n, min.tinge=0), labels=list(nrow=10))
pdf('out/cg/fig_heatmap.pdf', width=6.69, height=6.69, pointsize=6)
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
pdf('out/cg/fig_pca.pdf', width=6.69, height=6.69, pointsize=6)
plot(tmp$factor.score[, 1], tmp$factor.score[, 2], col='white', xlab='PC1', ylab='PC2')
text(tmp$factor.score[, 1], tmp$factor.score[, 2], labels=unlist(sapply(colnames(tmp.nreads), function(colname) { strsplit(colname, '\\|')[[1]][2] })))
dev.off()

###

### Local Variables:
### eval: (setq tmppath (file-name-directory (buffer-file-name)))
### eval: (add-to-list 'exec-path (concat tmppath "../.homebrew/bin"))
### eval: (setq default-directory (concat tmppath ".."))
### End:
