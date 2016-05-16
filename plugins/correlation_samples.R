#!/usr/bin/env Rscript

source('Rlib/STRTprepHelper.R', chdir=T)
helper <- STRTprepHelper$newPlugin(
  name='correlation_samples',
  required_packages=c('renozao/pkgmaker@develop', 'renozao/NMF'))
if(length(helper$options$ANNOTATIONS) == 1) {
  annotations <- list(helper$samples$annotations[, helper$options$ANNOTATIONS])
  names(annotations) <- helper$options$ANNOTATIONS
} else {
  annotations <- helper$samples$annotations[, helper$options$ANNOTATIONS]
}

solarized <- c('#dc322f', '#859900', '#268bd2', '#b58900', '#cb4b16', '#6c71c4', '#d33682', '#2aa198', '#002b36')

colors <- list()
for (i in 1:ncol(annotations)) {
  col <- colnames(annotations)[i]
  if (is.factor(annotations[, col])) {
    annotations[, col] <- factor(as.character(annotations[, col]))
    colors[[col]] <- rep(solarized, length=length(levels(annotations[, col])))
  } else
    colors[[col]] <- rep(solarized, length=ncol(annotations))[i]
}

nreads <- helper$expressions$fluctuated$mask$normalized_levels

if(nrow(nreads) > 3) {
  correlations <- cor(nreads, use='pairwise.complete.obs', method='spearman')
  library(NMF)
  nmf.options(grid.patch=T)
  correlation_samples <- aheatmap(
    correlations,
    annCol=annotations, annColors=colors,
    breaks=0, Rowv=F, revC=T, layout='lmL|dalm', fontsize=6, 
    filename=sprintf('%s.pdf', helper$output_prefix))
  save(correlation_samples,
       file=sprintf('%s.RData', helper$output_prefix), compress='gzip')
} else
  warning('Skipped - significant features were less than 3.')

###

### Local Variables:
### eval: (setq tmppath (file-name-directory (buffer-file-name)))
### eval: (add-to-list 'exec-path (concat tmppath "../.homebrew/bin"))
### eval: (setq default-directory (concat tmppath ".."))
### End:
