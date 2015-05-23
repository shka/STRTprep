#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
libid <- ifelse(is.na(args[1]), 'GBM', args[1])
th <- ifelse(is.na(args[2]), 50, as.numeric(args[2]))
if(is.na(args[3])) {
  seqs <- c('tmp/GBM.0.step1b.stat', 'tmp/GBM.1.step1b.stat')
} else {
  seqs <- args[3:length(args)]
}

mat <- NA
for(j in 1:length(seqs)) {
  tmp <- read.table(seqs[j])
  tmp2 <- data.frame(x=1:(max(tmp[, 2])+1), y=0)
  for(i in 1:nrow(tmp)) tmp2[tmp[i, 2]+1, 'y'] <- tmp[i, 1]
  tmp3 <- data.frame(x=tmp2[, 'x']-1, y=(sum(tmp2[, 'y'])-cumsum(c(0, tmp2[, 'y'][1:(nrow(tmp2)-1)])))/sum(tmp2[, 'y']))[1:max(tmp[, 2]), ]
  colnames(tmp3) <- c('x', sprintf('y%d', j))
  if(! is.data.frame(mat))
    mat <- tmp3
  else
    mat <- merge(mat, tmp3, all=T)
}

pdf(sprintf('out/fig.%s.qualifiedReads.pdf', libid), width=3.34, height=3.34)
par(mar=c(3.5, 3.5, 2.5, .5)+.1, mgp=c(2.5, 1, 0))
plot(mat[, 1], mat[, 2]*100, type='l', ylim=range(mat[, 2:3]*100, na.rm=T), 
     xlab='Low quality base position', ylab='Qualified reads (%)', col='white',
     main=libid)
abline(v=th, col='red')
for(i in 2:ncol(mat))
  lines(mat[, 1], mat[, i]*100, lty=i-1)
legend(0, min(mat[,2:ncol(mat)], na.rm=T)*100, xjust=0, yjust=0,
       legend=c(sprintf("Fastq file %d", 1:length(seqs)), 'Length threshold'),
       col=c(rep('black', length(seqs)), 'red'),
       lty=c(1:length(seqs), 1),
       bty='n')
dev.off()
