library(MASS)

parse_libid <- function(colname) unlist(strsplit(colname, '\\.'))[1]

rowCV2s <- function(reads) {
    reads.rowVars <- apply(reads, 1, var)
    reads.rowVars/(rowMeans(reads)^2)
}

extract_spikein_reads <- function(reads)
    reads[which(substr(rownames(reads), 1, 10) == 'RNA_SPIKE_'), ]

estimate_errormodels <- function(nreads) {
    libids <- unique(sapply(colnames(nreads), parse_libid))
    tmp <- lapply(libids, function(libid) {
        cat(sprintf('estimate_errormodels : %s\n', libid))
        targets <- sapply(colnames(nreads),
                          function(colname) { parse_libid(colname) == libid })
        nreads.tmp <- nreads[, which(targets)]
        nreads.spike <- extract_spikein_reads(nreads.tmp)
        x.spike <- 1/rowMeans(nreads.spike)
        y.spike <- rowCV2s(nreads.spike)
        model <- glm(y.spike ~ x.spike, family=Gamma(link='identity'))
        shape <- gamma.shape(model, it.lim=1000)$alpha
        scales <- (model$coefficients[1]+model$coefficients[2]/rowMeans(nreads.tmp))/shape
        responses <- rowCV2s(nreads.tmp)
        list(model=model, shape=shape, scales=scales, responses=responses)
    })
    names(tmp) <- libids
    tmp
}

create_scale_matrix <- function(errormodels) {
    tmp <- data.frame(errormodels[[1]]$scales)
    if(length(errormodels) > 1) {
        for(i in 2:length(errormodels))
            tmp <- cbind(tmp, data.frame(errormodels[[i]]$scales))
    }
    colnames(tmp) <- names(errormodels)
    as.matrix(tmp)
}

create_response_matrix <- function(errormodels) {
    tmp <- data.frame(errormodels[[1]]$responses)
    if(length(errormodels) > 1) {
        for(i in 2:length(errormodels))
            tmp <- cbind(tmp, data.frame(errormodels[[i]]$responses))
    }
    colnames(tmp) <- names(errormodels)
    as.matrix(tmp)
}

create_shape_vector <- function(errormodels) {
    tmp <- c(errormodels[[1]]$shape)
    if(length(errormodels) > 1) {
        for(i in 2:length(errormodels))
            tmp <- c(tmp, errormodels[[i]]$shape)
    }
    names(tmp) <- names(errormodels)
    tmp
}

merge_errormodels <- function(errormodels) {
    tmp.scales <- create_scale_matrix(errormodels)
    scales <- tmp.scales[, 1]
    tmp.responses <- create_response_matrix(errormodels)
    tmp.responses.scaled <- tmp.responses*(rep(scales, ncol(tmp.scales))/tmp.scales)
    responses <- rowSums(tmp.responses.scaled)
    tmp.shapes <- create_shape_vector(errormodels)
    shape <- sum(tmp.shapes)
    list(responses=responses, shape=shape, scales=scales)
}

estimate_fluctuation <- function(errormodel) {
    p.adj <- p.adjust(pgamma(q=errormodel$responses, shape=errormodel$shape, scale=errormodel$scale, lower.tail=F), 'fdr')
    score <- errormodel$responses/(errormodel$scale*errormodel$shape)
    list(p.adj=p.adj, score=score)
}
