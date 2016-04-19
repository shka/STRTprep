#!/usr/bin/env Rscript

source('Rlib/STRTprepHelper.R', chdir=T)
helper <- STRTprepHelper$newPlugin(
  name='heatmap_diffexp',
  required_packages=c('renozao/pkgmaker@develop', 'renozao/NMF'))
options <- helper$options
if(helper$comparison != 'global') {
  annotations <-
    helper$samples$annotations[, union(options$ANNOTATIONS, 'CLASS')]
  if(is.null(options$LABELS)) {
    annotations[, 'CLASS'] <- helper$samples$annotations[, 'CLASS']
  } else {
    annotations[, 'CLASS'] <-
      factor(helper$samples$annotations[, 'CLASS'], labels=options$LABELS)
  }
} else {
  annotations <- helper$samples$annotations[, options$ANNOTATIONS]
}

nreads <- helper$expressions$offset$significant$normalized_levels

distfun <- function(x) as.dist((1-cor(t(x), method='spearman'))/2)
clustfun <- function(d) hclust(d, method='ward.D2')
scaleBlackRedYellow <-
  function(n) colorRampPalette(c('black', 'red', 'yellow', 'white'),
                               space = "rgb")(n)

if(nrow(nreads) > 3) {
  library(NMF)
  nmf.options(grid.patch=T)
  heatmap <- aheatmap(
    log10(nreads),
    hclustfun=clustfun, distfun=distfun,
    ## color=scaleBlackRedYellow(100), breaks=0,
    scale='row',
    annCol=annotations, layout='dlmL|dalm',
    fontsize=6, filename=sprintf('%s.pdf', helper$output_prefix))
  save(heatmap,
       file=sprintf('%s.RData', helper$output_prefix), compress='gzip')
} else
  warning('Skipped - significant features were less than 3.')

###

### Local Variables:
### eval: (setq tmppath (file-name-directory (buffer-file-name)))
### eval: (add-to-list 'exec-path (concat tmppath "../.homebrew/bin"))
### eval: (setq default-directory (concat tmppath ".."))
### End:
