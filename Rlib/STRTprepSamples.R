library(R6)

STRTprepSamples <- R6Class(
  "STRTprepSamples",
  
  private = list(
    content = NA,
    helper = NA
  ),
  
  public = list(
    initialize = function(helper) {
      private$helper = helper
      tmp <- read.table(
        helper$samples_path, header=T, check.names=F, sep=',', quote='', colClasses=list(WELL="character"))
      rownames(tmp) <- tmp[, which(colnames(tmp) == 'NAME')]
      comparison <- helper$comparison
      if(comparison != 'global') {
        class_name <- sprintf('CLASS.%s', comparison)
        block_name <- sprintf('BLOCK.%s', comparison)
        classification <- tmp[, class_name]
        tmp <- tmp[which(!is.na(classification)), ]
        tmp2 <- colnames(tmp)
        tmp2[which(tmp2 == class_name)] <- 'CLASS'
        tmp2[which(tmp2 == block_name)] <- 'BLOCK'
        colnames(tmp) <- tmp2
      }
      private$content <- tmp
    }),
  
  active = list(
    annotations = function() {
      tmp <- colnames(private$content)
      private$content[, which(substr(tmp, 1, 6) != 'CLASS.'
                              & substr(tmp, 1, 6) != 'BLOCK.')]
    },
    
    names = function() rownames(private$content)))
