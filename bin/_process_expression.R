extract_detected <- function(reads.uniq.file) {
    reads.uniq.head <- read.table(reads.uniq.file, nrow=10)
    reads.uniq.classes <- c('character', rep('numeric', ncol(reads.uniq.head)))
    reads.uniq <- read.table(reads.uniq.file, check.names=F, comment.char='',
                             colClasses=reads.uniq.classes)
    
    reads.uniq.detected <-
        reads.uniq[which(apply(reads.uniq, 1, function(rows) { max(rows)>4 })), ]
    reads.uniq.detected.file <- sub('.txt', '.detected.RData', reads.uniq.file)
    gz <- gzfile(reads.uniq.detected.file, 'wb')
    save(reads.uniq.detected, file=gz)
    close(gz)

    reads.uniq.detected.file
}

normalize <- function(reads.uniq.detected.file) {
    load(reads.uniq.detected.file)
    
    reads.uniq.detected.spikein <- reads.uniq.detected[which(substr(rownames(reads.uniq.detected), 1, 10) == 'RNA_SPIKE_'), ]
    reads.uniq.detected.spikein.colSums <- colSums(reads.uniq.detected.spikein)
    reads.uniq.detected.spikein.colSums.noOutlier <- reads.uniq.detected.spikein.colSums[which(mean(reads.uniq.detected.spikein.colSums)-3*sd(reads.uniq.detected.spikein.colSums) <= reads.uniq.detected.spikein.colSums & reads.uniq.detected.spikein.colSums <= mean(reads.uniq.detected.spikein.colSums)+3*sd(reads.uniq.detected.spikein.colSums))]
    reads.uniq.detected.spikein.scale <- mean(reads.uniq.detected.spikein.colSums.noOutlier)/reads.uniq.detected.spikein.colSums
    reads.uniq.detected.normalized <- reads.uniq.detected*rep(reads.uniq.detected.spikein.scale, each=nrow(reads.uniq.detected))

    reads.uniq.detected.normalized.file <- sub('.detected', '.detected.normalized', reads.uniq.detected.file)
    gz <- gzfile(reads.uniq.detected.normalized.file, 'wb')
    save(reads.uniq.detected.normalized, file=gz)
    close(gz)

    reads.uniq.detected.normalized.file
}

plot_overdispersion2poisson <- function(reads.uniq.detected.normalized.file) {
    load(reads.uniq.detected.normalized.file)

    obs <- reads.uniq.detected.normalized[which(substr(rownames(reads.uniq.detected.normalized), 1, 10) == 'RNA_SPIKE_'), ]
    obs.rowMeans <- rowMeans(obs)
    
    x <- rep(as.numeric(factor(obs.rowMeans)), times=ncol(obs))
    y.obs <- as.matrix(obs)
    y.exp.1 <- t(sapply(obs.rowMeans, function(x) { rpois(ncol(obs), x) }))
    
    xlim <- range(x)
    ylim <- range(y.obs)
    cex <- 1/3

    pdf(sub('.reads.uniq.detected.normalized.RData.gz', '.fig.overdisperson2poisson.pdf', sub('out/exp', 'out/stat', reads.uniq.detected.normalized.file)),
        width=3.34, height=3.34/2, pointsize=8)
    par(mar=c(3, 3, 1, 1), mgp=c(1.9, 1, 0))
    plot(x-.15, y.obs, xlim=xlim, ylim=ylim, cex=cex, col='black',
         xlab=sprintf('%d detected spike-ins', length(obs.rowMeans)), 
         ylab='Normalized spike-in reads', axes=F)
    axis(2)
    par(new=T)
    plot(x+.15, y.exp.1, xlim=xlim, ylim=ylim, cex=cex, xlab='', ylab='',
         bty='n', axes=F, col='red', pch=2)
    legend(xlim[1], ylim[2], col=c('black', 'red'), pch=1:3, cex=6/8,
           legend=c(sprintf('Observed (n=%d)', ncol(obs)),
               'Poisson'))
    dev.off()
}

fit <- function(reads.uniq.detected.normalized.file) {
    load(reads.uniq.detected.normalized.file)

    obs.spike <- reads.uniq.detected.normalized[which(substr(rownames(reads.uniq.detected.normalized), 1, 10) == 'RNA_SPIKE_'), ]
    obs.spike.rowMeans <- rowMeans(obs.spike)
    obs.spike.rowVars <- apply(obs.spike, 1, var)
    obs.spike.invrowCV2s <- 1/(sqrt(obs.spike.rowVars)/obs.spike.rowMeans)^2

    x.spike <- log10(obs.spike.rowMeans)
    y.spike <- log10(obs.spike.invrowCV2s)
    errormodel.uniq.detected.normalized <- lm(y.spike ~ x.spike)

    errormodel.uniq.detected.normalized.file <- sub('reads.', 'errormodel.', reads.uniq.detected.normalized.file)
    gz <- gzfile(errormodel.uniq.detected.normalized.file, 'wb')
    save(errormodel.uniq.detected.normalized, file=gz)
    close(gz)

    errormodel.uniq.detected.normalized.file
}

plot_errormodel <- function(errormodel.uniq.detected.normalized.file) {
    load(errormodel.uniq.detected.normalized.file)

    errormodel <- errormodel.uniq.detected.normalized
    x.spike <- errormodel$model$x.spike
    y.spike <- errormodel$model$y.spike
    x.fit <- sort(errormodel$model$x.spike)
    errormodel.predicted <- predict(errormodel, newdata=data.frame(x.spike=x.fit), interval='prediction', level=.99)

    options(scipen=3)
    tryCatch ({
        pdf(sub('.errormodel.uniq.detected.normalized.RData.gz', '.fig.errormodel.pdf', sub('out/exp', 'out/stat', errormodel.uniq.detected.normalized.file)),
            width=3.34, height=3.34, pointsize=8)
        par(mar=c(3.2, 3.2, 1, 1), mgp=c(1.9, 0.8, 0))
        ## observation
        plot(x.spike, y.spike, cex=.5, col='black', xlab=expression(paste(log[10], '(normalized spike-in read count)')), ylab=expression(log[10](1/CV^2)))
        ## Poisson
        points(x.fit, sapply(x.fit, function(l) { log10(1/(sqrt(10^l)/10^l)^2) }), type='b', lty=1, col='red', pch=2, cex=.5)
        ## Regression line
        lines(x.fit, errormodel.predicted[, 1], col='blue', type='b', pch=3, cex=.5)
        lines(x.fit, errormodel.predicted[, 2], col='blue', lty=2)
        lines(x.fit, errormodel.predicted[, 3], col='blue', lty=2)
        ##
        legend(x=max(x.spike), y=min(y.spike), xjust=1, yjust=0, legend=c(sprintf('%d detected spike-ins', length(x.spike)), 'Poisson', sprintf("Linear regression; %1.5f+%1.5f*x", errormodel$coefficients[1], errormodel$coefficients[2]), 'Linear reguression; 99% confidence'), lty=c(NA, 2, 1, 3), pch=c(1, 2, 3, NA), col=c('black', 'red', 'blue', 'blue'), cex=6/8)
    }, finally={ dev.off() })
    options(scipen=0)
}

extract_fluctuated <- function(reads.uniq.detected.normalized.file, errormodel.uniq.detected.normalized.file) {
    load(reads.uniq.detected.normalized.file)
    load(errormodel.uniq.detected.normalized.file)

    obs <- reads.uniq.detected.normalized
    obs.rowMeans <- rowMeans(obs)
    obs.rowVars <- apply(obs, 1, var)
    obs.invrowCV2s <- 1/(sqrt(obs.rowVars)/obs.rowMeans)^2
    x <- log10(obs.rowMeans)
    y <- log10(obs.invrowCV2s)

    errormodel <- errormodel.uniq.detected.normalized
    errormodel.residuals.mean <- mean(residuals(errormodel))
    errormodel.residuals.sd <- sd(residuals(errormodel))
    errormodel.pred.all.residuals <- y-(errormodel$coefficients[1]+errormodel$coefficients[2]*x)
    errormodel.pred.all.residuals.p <- pnorm(errormodel.pred.all.residuals, mean=errormodel.residuals.mean, sd=errormodel.residuals.sd)
    reads.uniq.detected.normalized.fluctuated <- obs[which(errormodel.pred.all.residuals.p<0.01/length(x)), ]

    reads.uniq.detected.normalized.fluctuated.file <- sub('.normalized', '.normalized.fluctuated', reads.uniq.detected.normalized.file)
    gz <- gzfile(reads.uniq.detected.normalized.fluctuated.file, 'wb')
    save(reads.uniq.detected.normalized.fluctuated, file=gz)
    close(gz)
    
    reads.uniq.detected.normalized.fluctuated.file
}

plot_fluctuated <- function(reads.uniq.detected.normalized.file, reads.uniq.detected.normalized.fluctuated.file, errormodel.uniq.detected.normalized.file) {
    load(reads.uniq.detected.normalized.file)
    load(reads.uniq.detected.normalized.fluctuated.file)
    load(errormodel.uniq.detected.normalized.file)

    obs <- reads.uniq.detected.normalized
    obs.rowMeans <- rowMeans(obs)
    obs.rowVars <- apply(obs, 1, var)
    obs.invrowCV2s <- 1/(sqrt(obs.rowVars)/obs.rowMeans)^2
    
    x <- sort(log10(obs.rowMeans))
    y.obs <- log10(obs.invrowCV2s[names(x)])

    xlim <- range(x)
    ylim <- range(y.obs)
    cex <- .5
    x.fit <- seq(xlim[1], xlim[2], .1)

    fluctuated <- reads.uniq.detected.normalized.fluctuated
    fluctuated.rowMeans <- rowMeans(fluctuated)
    fluctuated.rowVars <- apply(fluctuated, 1, var)
    fluctuated.invrowCV2s <- 1/(sqrt(fluctuated.rowVars)/fluctuated.rowMeans)^2

    x.fluctuated <- log10(fluctuated.rowMeans)
    y.fluctuated <- log10(fluctuated.invrowCV2s)

    errormodel <- errormodel.uniq.detected.normalized
    x.spike <- errormodel$model$x.spike
    y.spike <- errormodel$model$y.spike
    errormodel.pred.fit <- predict(errormodel, newdata=data.frame(x.spike=x.fit), interval='prediction', level=.99)
    
    options(scipen=3)
    tryCatch ({
        png(sub('.errormodel.uniq.detected.normalized.RData.gz', '.fig.fluctuated.png', sub('out/exp', 'out/stat', errormodel.uniq.detected.normalized.file)),
            width=3.34, height=3.34, units='in', res=600, pointsize=8)
        par(mar=c(3.2, 3.2, 1, 1), mgp=c(1.9, 0.8, 0))
        ## observation
        plot(x, y.obs, cex=cex, col='gray', xlab=expression(paste(log[10], '(normalized spike-in read count)')), ylab=expression(log[10](1/CV^2)))
        points(x.spike, y.spike, cex=cex, col='black')
        points(x.fluctuated, y.fluctuated, cex=cex, col='magenta')
        ## Poisson
        lines(x.fit, sapply(x.fit, function(l) { log10(1/((sqrt(10^l)/10^l)^2)) }), type='l', col='red')
        ## Linear regression
        lines(x.fit, errormodel.pred.fit[, 1], col='blue')
        lines(x.fit, errormodel.pred.fit[, 2], col='blue', lty=2)
        lines(x.fit, errormodel.pred.fit[, 3], col='blue', lty=2)
        ## ##
        legend(x=min(x.fit), y=max(y.obs), legend=c(sprintf('Observed (%d regions)', length(obs.rowMeans)), sprintf('%d detected spike-ins', length(x.spike)), sprintf('Fluctuated (%d regions; P<0.01 w/ Bonferroni corr.)', nrow(fluctuated)), 'Poisson', sprintf("Linear regression; %1.5f+%1.5f*x", errormodel$coefficients[1], errormodel$coefficients[2]), 'Linear reguression; 99% confidence'), lty=c(NA, NA, NA, 2, 1, 3), pch=c(1, 1, 1, 2, 3, NA), col=c('gray', 'black', 'magenta', 'red', 'blue', 'blue'), cex=6/8)
    }, finally={ dev.off() })
    options(scipen=0)
}

args <- commandArgs(trailingOnly=T)
reads.uniq.detected.file <- extract_detected(args[1])
reads.uniq.detected.normalized.file <- normalize(reads.uniq.detected.file)
plot_overdispersion2poisson(reads.uniq.detected.normalized.file)
errormodel.uniq.detected.normalized.file <- fit(reads.uniq.detected.normalized.file)
plot_errormodel(errormodel.uniq.detected.normalized.file)
reads.uniq.detected.normalized.fluctuated.file <- extract_fluctuated(reads.uniq.detected.normalized.file, errormodel.uniq.detected.normalized.file)
plot_fluctuated(reads.uniq.detected.normalized.file, reads.uniq.detected.normalized.fluctuated.file, errormodel.uniq.detected.normalized.file)
