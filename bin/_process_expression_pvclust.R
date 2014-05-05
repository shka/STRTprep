library(pvclust)
library(snow)

args <- commandArgs(trailingOnly=T)
nreads.uniq.success.fluctuated.file <- args[1]
cores <- args[2]
lib <- sub('out/exp/', '', strsplit(nreads.uniq.success.fluctuated.file, '\\.')[[1]][1])
pvclust.uniq.success.fluctuated.file <- sub('nreads.', 'pvclust.', nreads.uniq.success.fluctuated.file)

load(nreads.uniq.success.fluctuated.file)
nreads <- as.matrix(nreads.uniq.success.fluctuated)
nreads[which(nreads == 0)] <- NA
cl <- makeCluster(cores, type='MPI')
pvclust.uniq.success.fluctuated <- parPvclust(cl, log10(nreads), method.dist='correlation', method.hclust='average', nboot=1000)
stopCluster(cl)

gz <- gzfile(pvclust.uniq.success.fluctuated.file, 'wb')
save(pvclust.uniq.success.fluctuated, file=gz)
close(gz)

## load(pvclust.uniq.success.fluctuated.file)

pvclust.uniq.success.fluctuated$hclust$labels <- sub(sprintf('%s.', lib), '', pvclust.uniq.success.fluctuated$hclust$labels)
pdf(sub('.RData.gz', '.pdf', pvclust.uniq.success.fluctuated.file), width=6.69/2, height=6.69/3, pointsize=4)
plot(pvclust.uniq.success.fluctuated, sub=lib, xlab='Distance: correlation, cluster method: average')
pvrect(pvclust.uniq.success.fluctuated)
dev.off()
