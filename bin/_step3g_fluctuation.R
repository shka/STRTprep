#!/usr/bin/env Rscript

source('Rlib/fluctuation.R')

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

tmp <- data.frame(Gene=names(fluctuation$p.adj), pvalue=fluctuation$p.adj, score=fluctuation$score)
gz <- gzfile(sprintf('%s/fluctuation.txt.gz', dir), 'w')
write.table(tmp, file=gz, quote=F, sep="\t", row.names=F, col.names=T)
close(gz)

###

### Local Variables:
### eval: (setq tmppath (file-name-directory (buffer-file-name)))
### eval: (add-to-list 'exec-path (concat tmppath "../.homebrew/bin"))
### eval: (setq default-directory (concat tmppath ".."))
### End:
