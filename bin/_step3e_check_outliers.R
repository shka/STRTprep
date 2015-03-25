draw_outliers <- function(x, y, ylab) {
    tmp.names <- unique(x)
    tmp.stats <- boxplot(y ~ x, plot=F)
    tmp.stats.mean <- colMeans(tmp.stats$stats)
    tmp.stats.mean.stats <- boxplot(tmp.stats.mean, plot=F)
    for(i in which(is.element(tmp.stats.mean, tmp.stats.mean.stats$out))) {
        tmp.names[i] <- sprintf("* %s", tmp.names[i])
    }
    par(las=3, mar=c(6, 4, 2, 1), mgp=c(2.5, 1, 0))
    tmp <- boxplot(y ~ x, axes=F, ylab=ylab)
    axis(1, at=1:length(tmp.names), labels=tmp.names, lty=0)
    axis(2)
    tmp.stats
}

extract_outliers <- function(x, y, stat) {
    which(is.element(x, stat$names) & is.element(y, stat$out))
}

###

samples.all <- read.table('tmp/cg/samples.txt', header=T, sep="\t", quote='', check.names=F)
libwellids <- sprintf('%s.%s', samples.all[, 'LIBRARY'], samples.all[, 'WELL'])
samples <- samples.all[which(!is.na(samples.all[, 'NAME'])), ]

###

draw_outliers_spikeinReads <- function(samples) {
    draw_outliers(samples[, 'LIBRARY'],
                  log10(samples[, 'SPIKEIN_READS']),
                  expression(log[10]('Spike-in reads')))
}

extract_outliers_spikeinReads <- function(samples, stat) {
    tmp <- extract_outliers(samples[, 'LIBRARY'],
                            log10(samples[, 'SPIKEIN_READS']),
                            stat)
    sprintf("%s.%s", samples[tmp, 'LIBRARY'], samples[tmp, 'WELL'])
}

pdf('out/cg/fig_outliers_spikeinReads.pdf', width=1.13, height=2.26, pointsize=6)
tmp <- draw_outliers_spikeinReads(samples)
dev.off()

samples.all[, 'SPIKEIN_READS.OUTLIER'] <-
    is.element(libwellids, extract_outliers_spikeinReads(samples, tmp))

###

draw_outliers_mappedPerSpikein <- function(samples) {
    draw_outliers(samples[, 'LIBRARY'],
                  log10(samples[, 'MAPPED/SPIKEIN']),
                  'Mapped reads / Spike-in reads')
}

extract_outliers_mappedPerSpikein <- function(samples, stat) {
    tmp <- extract_outliers(samples[, 'LIBRARY'],
                            log10(samples[, 'MAPPED/SPIKEIN']),
                            stat)
    sprintf("%s.%s", samples[tmp, 'LIBRARY'], samples[tmp, 'WELL'])
}

pdf('out/cg/fig_outliers_mappedPerSpikein.pdf', width=1.13, height=2.26, pointsize=6)
tmp <- draw_outliers_mappedPerSpikein(samples)
dev.off()

samples.all[, 'MAPPED/SPIKEIN.OUTLIER'] <-
    is.element(libwellids, extract_outliers_mappedPerSpikein(samples, tmp))

###

draw_outliers_spikeinCapture <- function(samples) {
    draw_outliers(samples[, 'LIBRARY'],
                  samples[, 'SPIKEIN_5END_RATE'],
                  "Spike-in 5'-end capture rate")
}

extract_outliers_spikeinCapture <- function(samples, stat) {
    tmp <- extract_outliers(samples[, 'LIBRARY'],
                            samples[, 'SPIKEIN_5END_RATE'],
                            stat)
    sprintf("%s.%s", samples[tmp, 'LIBRARY'], samples[tmp, 'WELL'])
}

pdf('out/cg/fig_outliers_spikeinCapture.pdf', width=1.13, height=2.26, pointsize=6)
tmp <- draw_outliers_spikeinCapture(samples)
dev.off()

samples.all[, 'SPIKEIN_5END_RATE.OUTLIER'] <-
    is.element(libwellids, extract_outliers_spikeinCapture(samples, tmp))

###

draw_outliers_codingCapture <- function(samples) {
    draw_outliers(samples[, 'LIBRARY'],
                  samples[, 'CODING_5END_RATE'],
                  "Spike-in 5'-end capture rate")
}

extract_outliers_codingCapture <- function(samples, stat) {
    tmp <- extract_outliers(samples[, 'LIBRARY'],
                            samples[, 'CODING_5END_RATE'],
                            stat)
    sprintf("%s.%s", samples[tmp, 'LIBRARY'], samples[tmp, 'WELL'])
}

pdf('out/cg/fig_outliers_codingCapture.pdf', width=1.13, height=2.26, pointsize=6)
tmp <- draw_outliers_codingCapture(samples)
dev.off()

samples.all[, 'CODING_5END_RATE.OUTLIER'] <-
    is.element(libwellids, extract_outliers_codingCapture(samples, tmp))

###

write.table(samples.all, 'out/cg/samples_all.txt', quote=F, sep="\t", row.names=F, col.names=T)

samples <- samples.all[which(!is.na(samples.all[, 'NAME'])
                             & !samples.all[, 'SPIKEIN_READS.OUTLIER']
                             & !samples.all[, 'MAPPED/SPIKEIN.OUTLIER']
                             & !samples.all[, 'SPIKEIN_5END_RATE.OUTLIER']
                             & !samples.all[, 'CODING_5END_RATE.OUTLIER']), ]
write.table(samples, 'out/cg/samples.txt', quote=F, sep="\t", row.names=F, col.names=T)

###

### Local Variables:
### eval: (setq tmppath (file-name-directory (buffer-file-name)))
### eval: (add-to-list 'exec-path (concat tmppath "../.homebrew/bin"))
### eval: (setq default-directory (concat tmppath ".."))
### End:
