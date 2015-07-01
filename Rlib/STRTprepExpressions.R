library(R6)

STRTprepExpressions <- R6Class(
  "STRTprepExpressions",
  
  private = list(
    content = NA,
    helper = NA,
    samples = NA,
    
    columns_of_level = function() {
      c(private$columns_of_normalized_level(),
        private$columns_of_raw_level())
    },
    
    columns_of_normalized_level = function()
      which(sapply(colnames(private$content),
                   function(n) {
                     tmp <- strsplit(n, '\\|')[[1]]
                     tmp[1] == 'N' & is.element(tmp[3], private$samples$names)
                   })),
    
    columns_of_raw_level = function()
      which(sapply(colnames(private$content),
                   function(n) {
                     tmp <- strsplit(n, '\\|')[[1]]
                     tmp[1] == 'R' & is.element(tmp[3], private$samples$names)
                   }))),
  
  public = list(
    initialize = function(helper) {
      private$samples <- helper$samples
      private$content <- read.table(
        helper$expressions_path,
        header=T, check.names=F, sep=',', quote='', row.names=1)
                                        #
      tmp <- colnames(private$content)
      tmp2 <- helper$comparison
      if(tmp2 != 'global') {
        tmp[which(tmp == sprintf('diffexpScore.%s', tmp2))] <- 'diffexpScore'
        tmp[which(tmp == sprintf('pvalue.%s', tmp2))] <- 'pvalue'
        tmp[which(tmp == sprintf('qvalue.%s', tmp2))] <- 'qvalue'
      }
      tmp[which(tmp == sprintf('fluctuationScore.%s', tmp2))] <-
        'fluctuationScore'
      tmp[which(tmp == sprintf('fluctuation.%s', tmp2))] <- 'fluctuation'
      colnames(private$content) <- tmp
                                        #
      tmp <- private$columns_of_raw_level()
      private$content <-
        private$content[which(rowSums(private$content[, tmp]) > 0), ]
    }),
  
  ## getFluctuations = function() {
  ##   if(!is.data.frame(private$fluctuations)) {
  ##     tmp <- private$content[, c('fluctuation', 'fluctuationScore')]
  ##     private$fluctuations <- tmp[which(!is.na(tmp[, 1])), ]
  ##   }
  ##   private$fluctuations
  ## },
  
  ## getNormalizedReadsOfSignificantWithOffset =
  ##     function(excludeSpikeins=TRUE) {
  ##       nreads <- self$getNormalizedReads(excludeSpikeins)
  ##       nreads.fluctuated <-
  ##         self$significant()$getNormalizedReads(excludeSpikeins)
  ##       if(nrow(nreads.fluctuated) == 0)
  ##         nreads.fluctuated
  ##       else
  ##         nreads.fluctuated + min(nreads[which(nreads > 0)])
  ##     }
  
  
  active = list(
    mask = function() {
      for(i in private$columns_of_level()) {
        tmp <- private$content[, i]
        tmp[which(tmp == 0)] <- NA
        private$content[, i] <- tmp
      }
      self
    },

    normalized_levels = function() {
      tmp <- as.matrix(private$content[, private$columns_of_normalized_level()])
      colnames(tmp) <-
        sapply(colnames(tmp), function(n) strsplit(n, '\\|')[[1]][3])
      tmp
    },

    offset = function() {
      tmp <- as.matrix(private$content[, private$columns_of_normalized_level()])
      ofs <- min(tmp[which(tmp > 0)])/2
      for(i in private$columns_of_normalized_level())
        private$content[, i] <- private$content[, i]+ofs
      self
    },
  
    significant = function() {
      tmp <- private$content
      options <- helper$options
      if(helper$comparison != 'global') {
        targets1 <- 1:nrow(tmp)
        if(!is.null(options$DIFFEXPP))
          targets1 <- which(tmp[, 'pvalue'] < options$DIFFEXPP)
        targets2 <- 1:nrow(tmp)
        if(!is.null(options$DIFFEXPQ))
          targets2 <- which(tmp[, 'qvalue'] < options$DIFFEXPQ)
        targets3 <- 1:nrow(tmp)
        if(!is.null(options$FLUCTUATIONP))
          targets3 <- which(tmp[, 'fluctuation'] < options$FLUCTUATIONP)
        targets <- intersect(targets1, intersect(targets2, targets2))
      } else {
        targets <- 1:nrow(tmp)
        if(!is.null(options$FLUCTUATIONP))
          targets <- which(tmp[, 'fluctuation'] < options$FLUCTUATIONP)
      }
      private$content <- tmp[targets, ]
      self
    }))
  
