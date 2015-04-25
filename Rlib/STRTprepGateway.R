library(R6)
library(yaml)

STRTprepGateway <-
  R6Class(
      "STRTprepGateway",
      
      public = list(
          
          quantificationType = NA,
          expressionsPath = NA,
          comparisonClass = NA,
          samplesPath = NA,
          configurationPath = NA,
          configuration = NA,
          
          initialize = function() {
            args <- commandArgs(trailingOnly=T)
            self$quantificationType <- ifelse(is.na(args[1]), 'byGene', args[1])
            self$comparisonClass <- ifelse(is.na(args[2]), '0', args[2])
            self$expressionsPath <-
              ifelse(is.na(args[3]),
                     sprintf("out/%s/diffexp.csv", self$quantificationType),
                     args[3])
            self$samplesPath <-
              ifelse(is.na(args[4]),
                     sprintf("out/%s/samples.csv", self$quantificationType),
                     args[4])
            self$configurationPath <-
              ifelse(is.na(args[5]), 'conf.yaml', args[5])
            self$configuration <- yaml.load_file(self$configurationPath)
          },

          getExpressions = function() {
            if(!is.data.frame(private$expressions))
              private$expressions <-
                read.table(self$expressionsPath, header=T, check.names=F,
                           sep=',', quote='', row.names=1)
            private$expressions
          },
          
          getFluctuations = function() {
            if(!is.data.frame(private$fluctuations)) {
              e <- self$getExpressions()
              tmp <- 
                e[, c(sprintf('fluctuation.%s', self$comparisonClass),
                      sprintf('fluctuationScore.%s', self$comparisonClass))]
              private$fluctuations <- tmp[which(!is.na(tmp[, 1])), ]
              colnames(private$fluctuations) <- c('pvalue', 'score')
            }
            private$fluctuations
          }
          ),

      private = list(
          expressions = NA,
          fluctuations = NA
          )
      )

## nreads <- NA
## getNreads <- function() {
##     if(is.na(nreads)) {
##         expression <- getExpression()
##         tmp <- 
##             expression[, which(substr(colnames(expression), 1, 2) == 'N|')]
##         colnames(tmp) <-
##             sapply(colnames(tmp), function(n) { strsplit(n, '|')[3] })
##         nreads <<- tmp
##     }
##     return nreads
## }
