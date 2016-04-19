#!/usr/bin/env Rscript

source('Rlib/STRTprepHelper.R', chdir=T)
helper <- STRTprepHelper$newPlugin(
  name='transcripts',
  required_packages=c('beeswarm', 'RColorBrewer'))

library(beeswarm)
library(RColorBrewer)

annotations <- helper$samples$annotations
options <- helper$options

if(helper$comparison == 'global') {
  annotations[, 'CLASS'] <-
    apply(annotations, 1,
          function(r) paste(r[options$CLASSES], collapse='.'))
} else {
  annotations[, 'CLASS'] <-
    factor(annotations[, 'CLASS'], labels=options$LABELS)
}

pdf(sprintf("%s.pdf", helper$output_prefix), pointsize=6)
ppar <- par()
if(is.null(options$COLOR)) {
  color_indexes <- factor(rep(1, nrow(annotations)))
  color_palette <- c('black')
  width_legend <- 0
} else {
  color_indexes <- factor(annotations[, options$COLOR])
  color_palette <- c(brewer.pal(9, 'Set1'), brewer.pal(7, 'Set2'))
  width_legend <- max(strwidth(levels(color_indexes), units='inches'))
}
width_label <- max(strwidth(annotations[, 'CLASS'], units='inches'))
dev.off()

tmp <- data.frame(VALUE=annotations[, 'MAPPED/SPIKEIN'],
                  CLASS=annotations[, 'CLASS'])

pdf(sprintf("%s.pdf", helper$output_prefix),
    width=(.05+1+2.5+(length(unique(tmp[, 2]))+1)*2+1+1+.05)*ppar$csi+width_legend,
    height=2.23+width_label, pointsize=6)
par(mai=c(width_label+ppar$mgp[2]*ppar$csi,
          (1+2.5)*ppar$csi,
          .1,
          width_legend+(1+1)*ppar$csi)+.05, las=3, mgp=c(2.5, 1, 0))
transcripts <- boxplot(VALUE ~ CLASS, data=tmp, log='y', cex=2/3,
                       ylab='Relative amount of poly(A)+ transcripts')
beeswarm(VALUE ~ CLASS, data=tmp,
         pwcol=color_palette[color_indexes], add=T, cex=2/3, spacing=2/3)
if(width_legend > 0)
  legend(par()$usr[2], 10^par()$usr[4], yjust=1, legend=levels(color_indexes),
         xpd=TRUE, ncol=1, bty='n', col=color_palette, pch=1)
dev.off()

transcripts$values <- tmp$VALUE
transcripts$classes <- tmp$CLASS
save(transcripts,
     file=sprintf('%s.RData', helper$output_prefix), compress='gzip')

###

### Local Variables:
### eval: (setq tmppath (file-name-directory (buffer-file-name)))
### eval: (add-to-list 'exec-path (concat tmppath "../.homebrew/bin"))
### eval: (setq default-directory (concat tmppath ".."))
### End:
