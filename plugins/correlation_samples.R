#!/usr/bin/env Rscript

source('Rlib/STRTprepHelper.R', chdir=T)
helper <- STRTprepHelper$new(name='correlation_samples',
                             required_packages=c('renozao/pkgmaker@develop',
                                                 'renozao/NMF'))
annotations <- helper$samples$annotations[, helper$options$ANNOTATIONS]
nreads <- helper$expressions$significant$mask$normalized_levels

if(nrow(nreads) > 3) {
  correlations <- cor(nreads, use='pairwise.complete.obs', method='spearman')
  library(NMF)
  nmf.options(grid.patch=T)
  correlation_samples <- aheatmap(
    correlations,
    breaks=0, Rowv=F, revC=T, annCol=annotations, layout='lmL|dalm', fontsize=6, 
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
