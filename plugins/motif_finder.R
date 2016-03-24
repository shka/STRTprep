#!/usr/bin/env Rscript

source('Rlib/STRTprepHelper.R', chdir=T)
helper <- STRTprepHelper$new(name='motif_finder',
                             required_packages=c('rGADEM',
                                                 'snow',
                                                 'TFBSTools'))
if(helper$quantification != 'byTFE') {
  warning('Skipped - this plugin is only for TFE-based quantification.')
  q(save='no')
}

nreads <- helper$expressions$offset$significant$normalized_levels
if(nrow(nreads) < 20) {
  warning('Skipped - fluctuated TFEs were less than 20.')
  q(save='no')
}

library(parallel)
library(rGADEM)
library(snow)
library(TFBSTools)

## preparation of promoter sequences

fai <- read.table(sprintf('%s.fa.fai', helper$reference_genome), header=F, sep="\t", quote='', row.names=1, check.names=F)

bed.all <- helper$expressions$table[, c('Chr', 'Peak', 'Peak', 'Gene', 'Peak', 'Str', 'diffexpScore')]
bed.all[, 2] <- bed.all[, 2]-ifelse(bed.all[, 6] == '+', 2001, 500)
bed.all[, 3] <- bed.all[, 3]+ifelse(bed.all[, 6] == '+', 499, 2000)
bed.all[, 4] <- rownames(bed.all)
bed.all[, 5] <- 0
bed.all <- bed.all[which(bed.all[, 2] > 0 & bed.all[, 3] <= fai[as.character(bed.all[, 1]), 1]), ]

bed <- bed.all[which(is.element(rownames(bed.all), rownames(nreads))), ]
bed.pos <- bed[which(bed[, 7] > 0), ]
bed.neg <- bed[which(bed[, 7] < 0), ]

create_promoter_sequences <- function(bed, suffix) {
  gz <- gzfile(sprintf('%s_%s.bed.gz', helper$output_prefix, suffix), 'w')
  write.table(bed[, 1:6], file=gz, quote=F, sep="\t", row.names=F, col.names=F)
  close(gz)
  system(sprintf("bedtools getfasta -s -name -fi %s.fa -bed %s_%s.bed.gz -fo - | fold -w 50 | pigz -c > %s_%s.fa.gz",
                 helper$reference_genome,
                 helper$output_prefix, suffix,
                 helper$output_prefix, suffix))
}

create_promoter_sequences(bed.all, 'all')
create_promoter_sequences(bed.pos, 'pos')
if(nrow(bed.neg) > 0)
  create_promoter_sequences(bed.neg, 'neg')

system(sprintf("bedtools shuffle -i %s_all.bed.gz -g %s.fa.fai -seed 1 | pigz -c > %s_ref.bed.gz",
               helper$output_prefix,
               helper$reference_genome,
               helper$output_prefix))
system(sprintf("bedtools getfasta -s -name -fi %s.fa -bed %s_ref.bed.gz -fo - | fold -w 50 | pigz -c > %s_ref.fa.gz",
               helper$reference_genome,
               helper$output_prefix,
               helper$output_prefix))
  
## de-novo motif extraction

execute_GADEM <- function(bed, suffix) {
  pSeqs <- readDNAStringSet(sprintf('%s_%s.fa.gz', helper$output_prefix, suffix), 'fasta')
  gadem <- GADEM(pSeqs, verbose=1, seed=1)
  save(gadem, file=sprintf('%s_%s.gadem.RData', helper$output_prefix, suffix), compress='gzip')
  for(i in 1:length(gadem@motifList)) {
    motif <- gadem@motifList[[i]]
    pdf(sprintf('%s_%s.motif%d.pwm.pdf', helper$output_prefix, suffix, i),
        width=6.69, height=2.23, pointsize=6)
    plot(motif, main=sprintf('%s.motif%d', suffix, i))
    dev.off()
    write.table(t(motif@pwm), file=sprintf('%s_%s.motif%d.pwm.txt', helper$output_prefix, suffix, i), quote=F, sep="\t", row.names=F, col.names=F)
    tmp <- data.frame(TFE=sapply(motif@alignList, function(align) rownames(bed)[align@seqID]),
                      position=sapply(motif@alignList, function(align) align@pos),
                      strand=sapply(motif@alignList, function(align) align@strand),
                      pvalue=sapply(motif@alignList, function(align) align@pval))
    write.csv(tmp, file=sprintf('%s_%s.motif%d.alignment.csv', helper$output_prefix, suffix, i), quote=F, row.names=F)
  }
  gadem
}

gadem.pos <- execute_GADEM(bed.pos, 'pos')
if(nrow(bed.neg) > 0)
  gadem.neg <- execute_GADEM(bed.neg, 'neg')

## overall motif occurrence

pwms <- PWMatrixList()
for(i in 1:length(gadem.pos@motifList)) {
  pwms[[i]] <- PWMatrix(ID=sprintf('pos.motif%d', i),
                       profileMatrix=gadem.pos@motifList[[i]]@pwm)
}
if(nrow(bed.neg) > 0)
  for(i in 1:length(gadem.neg@motifList)) {
    pwms[[i+length(gadem.pos@motifList)]] <-
      PWMatrix(ID=sprintf('neg.motif%d', i),
               profileMatrix=gadem.neg@motifList[[i]]@pwm)
  }

pSearchSeq <- function(suffix) {
  pSeqs <- readDNAStringSet(sprintf('%s_%s.fa.gz', helper$output_prefix, suffix), 'fasta')
  tmp.sites <-
    parLapply(clusters, 0:(detectCores()-1),
              function(i, pwms, pSeqs) {
                targets <- which(i == (1:length(pSeqs)) %% detectCores())
                searchSeq(pwms, pSeqs[targets])
              }, pwms=pwms, pSeqs=pSeqs)
  sites <- tmp.sites[[1]]
  for(i in 2:length(tmp.sites)) sites <- c(sites, tmp.sites[[i]])
  rm(tmp.sites)
  save(sites, file=sprintf('%s_%s.sites.RData', helper$output_prefix, suffix), compress='gzip')
}

clusters <- makeCluster(detectCores())
clusterExport(clusters, c('searchSeq', 'detectCores'))
pSearchSeq('all')
pSearchSeq('ref')
stopCluster(clusters)

##

motifCountByRanges <- function(suffix) {
  load(sprintf('%s_%s.sites.RData', helper$output_prefix, suffix))
  tmp.seqs <- unlist(sapply(sites, function(site) rep(site@seqname, length(site@views))))
  tmp.ranges <- factor(floor(unlist(sapply(sites, function(site) site@views@ranges@start))/50), levels=0:49, labels=sprintf("%d", 0:49*50-2000))
  tmp.motifs <- unlist(sapply(sites, function(site) rep(site@pattern@ID, length(site@views))))
  data.frame(seq=tmp.seqs, range=tmp.ranges, motif=tmp.motifs)
}

mc <- motifCountByRanges('all')
mc.ref <- motifCountByRanges('ref')

##

plot_posPref <- function(tmp.motif, tmp.tc, suffix) {
  tmp <- as.matrix(table(mc[which(mc[, 'motif'] == tmp.motif), c('seq', 'range')]))
  tmp.mc.all <- apply(tmp, 2, function(col) length(which(col>0)))
  tmp.n.all <- rep(nrow(bed.all), length(tmp.mc.all))
  tmp.ref <- as.matrix(table(mc.ref[which(mc.ref[, 'motif'] == tmp.motif), c('seq', 'range')]))
  tmp.mc.ref <- apply(tmp.ref, 2, function(col) length(which(col>0)))
  tmp.n.ref <- rep(nrow(bed.all), length(tmp.mc.ref))
  tmp.mc.target <- apply(tmp[intersect(rownames(tmp), tmp.tc), ], 2, function(col) length(which(col>0)))
  tmp.n.target <- rep(length(tmp.tc), length(tmp.mc.target))
  
  tmp.all.ave <- qhyper(0.5, tmp.mc.all, tmp.n.all-tmp.mc.all, tmp.n.target)
  tmp.all.upper <- qhyper(0.05, tmp.mc.all, tmp.n.all-tmp.mc.all, tmp.n.target, lower.tail=F)
  tmp.all.lower <- qhyper(0.95, tmp.mc.all, tmp.n.all-tmp.mc.all, tmp.n.target, lower.tail=F)
  tmp.ref.ave <- qhyper(0.5, tmp.mc.ref, tmp.n.ref-tmp.mc.ref, tmp.n.target)
  tmp.ref.upper <- qhyper(0.05, tmp.mc.ref, tmp.n.ref-tmp.mc.ref, tmp.n.target, lower.tail=F)
  tmp.ref.lower <- qhyper(0.95, tmp.mc.ref, tmp.n.ref-tmp.mc.ref, tmp.n.target, lower.tail=F)
  
  tmp.ymax <- max(tmp.mc.target, tmp.all.upper, tmp.ref.upper)
  pdf(sprintf('%s_%s_%s.posPref.pdf', helper$output_prefix, tmp.motif, suffix), width=2.23, height=1.6725, pointsize=6)
  par(mar=c(3.5, 3.5, 2.5, 1), mgp=c(2.25, 1, 0))
  plot(0:49+.5, tmp.ref.ave, xlim=c(0, 50), ylim=c(0, tmp.ymax), type='l', col='#CCCCCC', xlab='Position', axes=F, ylab='TFEs', main=sprintf('%s in %d %s vs %d bg TFEs', tmp.motif, length(tmp.tc), suffix, nrow(bed.all)))
  polygon(c(0:49, 49:0)+.5, c(tmp.ref.upper, rev(tmp.ref.lower)), border=F, col='#CCCCCC88')
  lines(0:49+.5, tmp.all.ave, col='#6c71c4cc')
  polygon(c(0:49, 49:0)+.5, c(tmp.all.upper, rev(tmp.all.lower)), border=F, col='#6c71c444')
  abline(v=40, lty=2, col='darkgray')
  lines(0:49+.5, tmp.mc.target, col='#dc322f')
  axis(1, at=0:10*5, labels=0:10*250-2000)
  axis(2)
  dev.off()
}

for(i in 1:length(gadem.pos@motifList))
  plot_posPref(sprintf('pos.motif%d', i), as.character(bed.pos[, 4]), 'pos')
if(nrow(bed.neg) > 0) {
  for(i in 1:length(gadem.pos@motifList))
    plot_posPref(sprintf('pos.motif%d', i), as.character(bed.neg[, 4]), 'neg')
  for(i in 1:length(gadem.neg@motifList))
    plot_posPref(sprintf('neg.motif%d', i), as.character(bed.pos[, 4]), 'pos')
  for(i in 1:length(gadem.neg@motifList))
    plot_posPref(sprintf('neg.motif%d', i), as.character(bed.neg[, 4]), 'neg')
}

##

## tmp.motif_counts <- mc[, c(1, 2)]
## motif_counts <- as.matrix(table(tmp.motif_counts$seq, tmp.motif_counts$motif))
## nreads <- nreads[as.character(bed[, 4]),]
## motif_zscores <- scale(motif_counts)[rownames(nreads)[which(is.element(rownames(nreads), rownames(motif_counts)))], ]
