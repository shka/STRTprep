#!/usr/bin/env Rscript

draw_fluctuation <- function(fluctuations, sig=0.01, qt='byGene') {
    ord <- order(fluctuations$fluctuationScore, decreasing=T)
    score <- fluctuations$fluctuationScore[ord]
    p.adj <- -log10(fluctuations$fluctuation[ord])
    prop <- (1:length(score))/length(score)

    unit <- ifelse(qt == 'byGene', 'genes', 'TFEs')
    
    pos <- which(score>0)
    score.pos <- score[pos]
    p.adj.pos <- p.adj[pos]
    prop.pos <- prop[pos]
    
    isFinite <- which(is.finite(p.adj.pos))
    score.pos.isFinite <- score.pos[isFinite]
    p.adj.pos.isFinite <- p.adj.pos[isFinite]
    
    par(mar=c(4, 4, .5, 4), mgp=c(2.5, 1, 0))
    
    plot(range(score.pos), range(p.adj.pos.isFinite),
         log='x', col='white', axes=F,
         xlab='Fluctuation score',
         ylab=expression(-log[10]('adjusted p-value')))
    points(score.pos.isFinite, p.adj.pos.isFinite, cex=.25)
    axis(1)
    axis(2)
    
    par(new=T)
    
    plot(range(score.pos), range(prop.pos),
         log='x', col='white', xlab='', ylab='', axes=F)
    points(score.pos, prop.pos, cex=.25, col='red')
    axis(4, col='red')
    mtext(sprintf('Proportion (%d %s)', length(score), unit),
          side=4, line=2.5, col='red')

    tmp.sig <- length(which(p.adj > -log10(sig)))
    if(tmp.sig > 0) {
      tmp <- score[tmp.sig]
      abline(v=tmp, col='blue', lty=3)
      mtext(sprintf(' score=%.2f', tmp), line=-1, at=tmp, adj=0, col='blue')
      mtext(sprintf(' adj.p=%s', sig), line=-2, at=tmp, adj=0, col='blue')
      mtext(sprintf(' %d %s', length(which(score>=tmp)), unit),
            line=-3, at=tmp, adj=0, col='blue')
    }
}

source('Rlib/STRTprepGateway.R')
gw <- STRTprepGateway$new('flucuationSummary')

pdf(sprintf('%s.pdf', gw$outputPrefix), width=3.34, height=3.34)
draw_fluctuation(gw$getExpressions()$getFluctuations(),
                 sig=gw$configuration$DEFAULT$FLUCTUATION,
                 qt=gw$quantificationType)
dev.off()
