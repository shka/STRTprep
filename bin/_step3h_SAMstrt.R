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

source('bin/_fluctuation.R')

load(path_nreads)
nreads.org <- nreads
nreads <- nreads[rownames(reads), colnames(reads)]
errormodels <- estimate_errormodels(nreads)
save(errormodels, file=sprintf('%s/errormodels_diffexp%d.RData', dir, idx), compress='gzip')
errormodel <- merge_errormodels(errormodels)
save(errormodel, file=sprintf('%s/errormodel_diffexp%d.RData', dir, idx), compress='gzip')
fluctuation <- estimate_fluctuation(errormodel)
save(fluctuation, file=sprintf('%s/fluctuation_diffexp%d.RData', dir, idx), compress='gzip')
tmp <- data.frame(Gene=names(fluctuation$p.adj), pvalue=fluctuation$p.adj)
gz <- gzfile(sprintf('%s/fluctuation_diffexp%d.txt.gz', dir, idx), 'w')
write.table(tmp, file=gz, quote=F, sep="\t", row.names=F, col.names=T)
close(gz)

classes.uniq <- sort(unique(classes))
blocks.uniq <- sort(unique(blocks))
ann <- data.frame(class1=classes == 1, class2=classes==2)
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

library(Heatplus)
distfun <- function(x) as.dist((1-cor(t(x), method='spearman'))/2)
clustfun <- function(d) hclust(d, method='ward.D2')
row.diffexp <- intersect(diffexp[which(diffexp[, 'qvalue'] < q.diffexp), 'Gene'], names(fluctuation$p.adj)[which(fluctuation$p.adj < p.fluctuation)])
nreads.diffexp <- nreads[row.diffexp, ]
tmp.nreads.pre <- nreads.diffexp+min(nreads.org[which(nreads.org>0)])
tmp.nreads <- tmp.nreads.pre[setdiff(rownames(tmp.nreads.pre), rownames(extract_spikein_reads(tmp.nreads.pre))), ]
if(nrow(tmp.nreads)>1 & ncol(tmp.nreads)>1) {
    hm <- annHeatmap2(log10(tmp.nreads),
                      dendrogram=list(
                          clustfun=clustfun, distfun=distfun, lwd=.5),
                      annotation=list(inclRef=F, Col=list(data=ann)),
                      label=list(
                          Col=list(nrow=max(nchar(colnames(tmp.nreads)))/2.25),
                          Row=list(nrow=max(nchar(rownames(tmp.nreads)))/2.25)),
                      col=heat.colors, legend=2, scale='none')
    pdf(sprintf('%s/fig_heatmap_diffexp%d.pdf', dir, idx), width=6.69, height=6.69)
    plot(hm)
    dev.off()
}

###

### Local Variables:
### eval: (setq tmppath (file-name-directory (buffer-file-name)))
### eval: (add-to-list 'exec-path (concat tmppath "../.homebrew/bin"))
### eval: (setq default-directory (concat tmppath ".."))
### End:
