draw_outliers <- function(x, y, ylab) {
    tmp.stats <- boxplot(y ~ x, plot=F)
    par(las=3, mar=c(6, 4, 2, 1), mgp=c(2.5, 1, 0))
    tmp <- boxplot(y ~ x, axes=F, ylab=ylab)
    axis(1, at=1:length(tmp.stats$names), labels=tmp.stats$names, lty=0)
    axis(2)
    tmp.stats
}

extract_outliers <- function(stat, lowerOnly=T) {
  names(stat$out)[which(lowerOnly == FALSE
                  | stat$out < stat$stats[3, stat$group])]
}

###

samples.all <- read.table('tmp/byGene/samples.csv', header=T, sep=',', quote='', check.names=F)
if(!is.element('FORCE_APPROVAL', colnames(samples.all)))
    samples.all[, 'FORCE_APPROVAL'] <- FALSE 
libwellids <- sprintf('%s.%s', samples.all[, 'LIBRARY'], samples.all[, 'WELL'])
rownames(samples.all) <- libwellids
samples <- samples.all[which(!is.na(samples.all[, 'NAME'])), ]

###

draw_outliers_spikeinReads <- function(samples) {
    y <- log10(as.numeric(samples[, 'SPIKEIN_READS']))
    names(y) <- rownames(samples)
    draw_outliers(samples[, 'LIBRARY'], y,
                  expression(log[10]('Spike-in reads')))
}

pdf('out/byGene/fig_outliers_spikeinReads.pdf', width=1.13, height=2.26, pointsize=6)
tmp <- draw_outliers_spikeinReads(samples)
dev.off()

samples.all[, 'SPIKEIN_READS.OUTLIER'] <-
    is.element(libwellids, extract_outliers(tmp))

###

draw_outliers_mappedPerSpikein <- function(samples) {
    y <- log10(as.numeric(samples[, 'MAPPED/SPIKEIN']))
    names(y) <- rownames(samples)
    draw_outliers(samples[, 'LIBRARY'], y,
                  expression(log[10]('Mapped / Spike-in reads')))
}

pdf('out/byGene/fig_outliers_mappedPerSpikein.pdf', width=1.13, height=2.26, pointsize=6)
tmp <- draw_outliers_mappedPerSpikein(samples)
dev.off()

samples.all[, 'MAPPED/SPIKEIN.OUTLIER'] <-
    is.element(libwellids, extract_outliers(tmp))

###

draw_outliers_spikeinCapture <- function(samples) {
    y <- as.numeric(samples[, 'SPIKEIN_5END_RATE'])
    names(y) <- rownames(samples)
    draw_outliers(samples[, 'LIBRARY'], y,
                  "Spike-in 5'-end capture rate")
}

pdf('out/byGene/fig_outliers_spikeinCapture.pdf', width=1.13, height=2.26, pointsize=6)
tmp <- draw_outliers_spikeinCapture(samples)
dev.off()

samples.all[, 'SPIKEIN_5END_RATE.OUTLIER'] <-
  is.element(libwellids, extract_outliers(tmp))

###

draw_outliers_codingCapture <- function(samples) {
    y <- as.numeric(samples[, 'CODING_5END_RATE'])
    names(y) <- rownames(samples)
    draw_outliers(samples[, 'LIBRARY'], y,
                  "Coding gene 5'-end capture rate")
}

pdf('out/byGene/fig_outliers_codingCapture.pdf', width=1.13, height=2.26, pointsize=6)
tmp <- draw_outliers_codingCapture(samples)
dev.off()

samples.all[, 'CODING_5END_RATE.OUTLIER'] <-
    is.element(libwellids, extract_outliers(tmp))

###

samples.all <- samples.all[order(samples.all[, 'LIBRARY'], samples.all[, 'WELL']), ]
write.table(samples.all, 'out/byGene/samples_all.csv', quote=F, sep=',', row.names=F, col.names=T)

###

samples <-
    samples.all[which(!is.na(samples.all[, 'NAME'])
                      & (samples.all[, 'FORCE_APPROVAL']
                         | (!samples.all[, 'FORCE_APPROVAL']
                            & !samples.all[, 'SPIKEIN_READS.OUTLIER']
                            & !samples.all[, 'MAPPED/SPIKEIN.OUTLIER']
                            & !samples.all[, 'SPIKEIN_5END_RATE.OUTLIER']
                            & !samples.all[, 'CODING_5END_RATE.OUTLIER']))), ]
write.table(samples, 'out/byGene/samples.csv', quote=F, sep=',', row.names=F, col.names=T)

###

### Local Variables:
### eval: (setq tmppath (file-name-directory (buffer-file-name)))
### eval: (add-to-list 'exec-path (concat tmppath "../.homebrew/bin"))
### eval: (setq default-directory (concat tmppath ".."))
### End:
