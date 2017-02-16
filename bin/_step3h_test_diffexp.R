#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
idx <- ifelse(is.na(args[1]), 0, as.numeric(args[1]))
path_sources <-
  ifelse(is.na(args[2]), sprintf('out/byGene/sources%d.RData', idx), args[2])
path_diffexp <-
  ifelse(is.na(args[3]), sprintf('out/byGene/diffexp%d.txt.gz', idx), args[3])
dir <- ifelse(is.na(args[4]), 'out/byGene', args[4])

load(path_sources)

##

library(SAMstrt)

test_twoClasses <- function(reads, classes, blocks, dir, idx) {
  samfit <- SAMseq(reads, y=sprintf('%dBlock%d', classes, blocks),
                   random.seed=1, resp.type='Two class unpaired', nperms=1000)
  save(samfit, file=sprintf('%s/samfit%d.RData', dir, idx), compress='gzip')
  tmp.sig <- samr.compute.siggenes.table(samfit$samr.obj, samfit$del,
                                         samfit$samr.obj, samfit$delta.table,
                                         all.genes=TRUE)
  tmp.pv <-
    samr.pvalues.from.perms(samfit$samr.obj$tt, samfit$samr.obj$ttstar)
  tmp2 <- rbind(tmp.sig$genes.up, tmp.sig$genes.lo)
  rowidxs <- as.numeric(tmp2[, 1])-1
  data.frame(Gene=rownames(reads)[rowidxs],
             Score=as.numeric(tmp2[, 'Score(d)']),
             pvalue=p.adjust(tmp.pv[rowidxs], method='BH'),
             qvalue=as.numeric(tmp2[, 'q-value(%)'])/100)
}

test_pairedClasses <- function(reads, classes, blocks, dir, idx) {
  samfit <- SAMseq(reads, y=sprintf('%d', classes),
                   random.seed=1, resp.type='Two class paired', nperms=1000)
  save(samfit, file=sprintf('%s/samfit%d.RData', dir, idx), compress='gzip')
  tmp.sig <- samr.compute.siggenes.table(samfit$samr.obj, samfit$del,
                                         samfit$samr.obj, samfit$delta.table,
                                         all.genes=TRUE)
  tmp.pv <-
    samr.pvalues.from.perms(samfit$samr.obj$tt, samfit$samr.obj$ttstar)
  tmp2 <- rbind(tmp.sig$genes.up, tmp.sig$genes.lo)
  rowidxs <- as.numeric(tmp2[, 1])-1
  data.frame(Gene=rownames(reads)[rowidxs],
             Score=as.numeric(tmp2[, 'Score(d)']),
             pvalue=p.adjust(tmp.pv[rowidxs], method='BH'),
             qvalue=as.numeric(tmp2[, 'q-value(%)'])/100)
}

test_multiClasses <- function(reads, classes, blocks, dir, idx) {
  samfit <- SAMseq(reads, y=sprintf('%dBlock%d', classes, blocks),
                   random.seed=1, resp.type='Multiclass', nperms=1000)
  save(samfit, file=sprintf('%s/samfit%d.RData', dir, idx), compress='gzip')
  tmp.sig <- samr.compute.siggenes.table(samfit$samr.obj, samfit$del,
                                         samfit$samr.obj, samfit$delta.table,
                                         all.genes=TRUE)
  tmp.pv <-
    samr.pvalues.from.perms(samfit$samr.obj$tt, samfit$samr.obj$ttstar)
  tmp2 <- tmp.sig$genes.up
  rowidxs <- as.numeric(tmp2[, 1])-1
  data.frame(Gene=rownames(reads)[rowidxs],
             Score=as.numeric(tmp2[, 'Score(d)']),
             pvalue=p.adjust(tmp.pv[rowidxs], method='BH'),
             qvalue=as.numeric(tmp2[, 'q-value(%)'])/100)
}

classes.uniq <- as.numeric(unique(classes))

if(setequal(classes.uniq, c(1, 2))) {
  diffexp <- test_twoClasses(reads, classes, blocks, dir, idx)
} else if(setequal(classes.uniq[which(classes.uniq<0)], -classes.uniq[which(classes.uniq>0)])) {
  diffexp <- test_pairedClasses(reads, classes, blocks, dir, idx)
} else {
  diffexp <- test_multiClasses(reads, classes, blocks, dir, idx)
}

gz <- gzfile(path_diffexp, 'w')
write.table(diffexp, file=gz, quote=F, sep="\t", row.names=F, col.names=T)
close(gz)

##

source('Rlib/fluctuation.R')

errormodels <- estimate_errormodels(nreads, blocks=blocks)
errormodel <- merge_errormodels(errormodels)
fluctuation <- estimate_fluctuation(errormodel)
save(errormodels, errormodel, fluctuation,
     file=sprintf('%s/fluctuation%d.RData.RData', dir, idx), compress='gzip')
tmp <- data.frame(Gene=names(fluctuation$p.adj), pvalue=fluctuation$p.adj, score=fluctuation$score)
gz <- gzfile(sprintf('%s/fluctuation%d.txt.gz', dir, idx), 'w')
write.table(tmp, file=gz, quote=F, sep="\t", row.names=F, col.names=T)
close(gz)

###

### Local Variables:
### eval: (setq tmppath (file-name-directory (buffer-file-name)))
### eval: (add-to-list 'exec-path (concat tmppath "../.homebrew/bin"))
### eval: (setq default-directory (concat tmppath ".."))
### End:
