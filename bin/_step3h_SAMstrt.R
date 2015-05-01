args <- commandArgs(trailingOnly=T)
idx <- ifelse(is.na(args[1]), 0, as.numeric(args[1]))
path_samples <- ifelse(is.na(args[2]), 'out/byGene/samples.csv', args[2])
path_reads <- ifelse(is.na(args[3]), 'out/byGene/reads.RData', args[3])
path_nreads <- ifelse(is.na(args[3]), 'out/byGene/nreads.RData', args[4])
path_diffexp <- ifelse(is.na(args[5]), sprintf('out/byGene/diffexp%d.txt.gz', idx), args[5])
q.diffexp <- ifelse(is.na(args[6]), 0.05, as.numeric(args[6]))
p.fluctuation <- ifelse(is.na(args[7]), 0.05, as.numeric(args[7]))
dir <- ifelse(is.na(args[8]), 'out/byGene', args[8])

samples <- read.table(path_samples, header=T, sep=',', quote='', check.names=F)
tmp.classes <- samples[, sprintf('CLASS.%d', idx)]
tmp.blocks <- samples [, sprintf('BLOCK.%d', idx)]
names(tmp.classes) <- names(tmp.blocks) <- sprintf('%s.%s|%s', samples[, 'LIBRARY'], samples[, 'WELL'], samples[, 'NAME'])
classes <- tmp.classes[which(!is.na(tmp.classes))]
blocks <- tmp.blocks[which(!is.na(tmp.classes))]

load(path_reads)
tmp.reads <- reads[, names(classes)]
reads <- tmp.reads[which(rowSums(tmp.reads) > 0), ]

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
               pvalue=tmp.pv[rowidxs],
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
               pvalue=tmp.pv[rowidxs],
               qvalue=as.numeric(tmp2[, 'q-value(%)'])/100)
}

if(setequal(unique(classes), c(1, 2))) {
    diffexp <- test_twoClasses(reads, classes, blocks, dir, idx)
} else {
    diffexp <- test_multiClasses(reads, classes, blocks, dir, idx)
}

gz <- gzfile(path_diffexp, 'w')
write.table(diffexp, file=gz, quote=F, sep="\t", row.names=F, col.names=T)
close(gz)

###

source('Rlib/fluctuation.R')

load(path_nreads)
nreads.org <- nreads
nreads <- nreads[rownames(reads), colnames(reads)]
errormodels <- estimate_errormodels(nreads)
save(errormodels, file=sprintf('%s/errormodels_diffexp%d.RData', dir, idx), compress='gzip')
errormodel <- merge_errormodels(errormodels)
save(errormodel, file=sprintf('%s/errormodel_diffexp%d.RData', dir, idx), compress='gzip')
fluctuation <- estimate_fluctuation(errormodel)
save(fluctuation, file=sprintf('%s/fluctuation_diffexp%d.RData', dir, idx), compress='gzip')
tmp <- data.frame(Gene=names(fluctuation$p.adj), pvalue=fluctuation$p.adj, score=fluctuation$score)
gz <- gzfile(sprintf('%s/fluctuation_diffexp%d.txt.gz', dir, idx), 'w')
write.table(tmp, file=gz, quote=F, sep="\t", row.names=F, col.names=T)
close(gz)

###

### Local Variables:
### eval: (setq tmppath (file-name-directory (buffer-file-name)))
### eval: (add-to-list 'exec-path (concat tmppath "../.homebrew/bin"))
### eval: (setq default-directory (concat tmppath ".."))
### End:
