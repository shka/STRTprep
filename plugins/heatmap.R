#!/usr/bin/env Rscript

source('Rlib/STRTprepGateway.R')
gw <- STRTprepGateway$new('heatmap', requiredPackages=c('Heatplus'))

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

library(Heatplus)

distfun <- function(x) as.dist((1-cor(t(x), method='spearman'))/2)
clustfun <- function(d) hclust(d, method='ward.D2')

tmp.nreads <- nreads.fluctuated+min(nreads[which(nreads>0)])

pdf(sprintf('%s.pdf', prefix), width=6.69, height=6.69)

if(ex$comparisonClass == 'global') {
  hp <- annHeatmap2(log10(tmp.nreads), scale='none', col=heat.colors, legend=2,
                    dendrogram=list(clustfun=clustfun, distfun=distfun, lwd=.5))
  plot(hp)
} else {
  samples <- gw$getSamples()
  classes <- samples$getClasses()
  blocks <- samples$getBlocks()
  classes.uniq <- sort(unique(classes))
  blocks.uniq <- sort(unique(blocks))
  ann <- data.frame(class1=classes == 1, class2=classes==2)
  rownames(ann) <- samples$getNames()
  if(length(classes.uniq) > 2) {
    for(n in seq(3, length(classes.uniq))) {
      cls = classes.uniq[n]
      ann[, sprintf('class%d', cls)] <- classes == cls
    }
  }
  for(n in seq(1, length(blocks.uniq))) {
    blk = blocks.uniq[n]
    ann[, sprintf('block%d', blk)] <- blocks == blk
  }
  hp <- annHeatmap2(log10(tmp.nreads),
                    dendrogram=list(
                        clustfun=clustfun, distfun=distfun, lwd=.5),
                    annotation=list(inclRef=F, Col=list(data=ann)),
                    col=heat.colors, legend=2, scale='none')
  plot(hp)
}

dev.off()

save(hp, file=sprintf('%s.RData', prefix), compress='gzip')
