#!/usr/bin/env Rscript

source('Rlib/STRTprepGateway.R')
gw <- STRTprepGateway$new()

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

fluctuations <- gw$getFluctuations()
p.fluctuation <- gw$configuration$DEFAULT$FLUCTUATION

draw_fluctuation(fluctuations, p.fluctuation)
