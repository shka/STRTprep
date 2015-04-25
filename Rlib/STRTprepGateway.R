library(R6)
library(yaml)

STRTprepGateway <-
  R6Class(
      "STRTprepGateway",
      
      public = list(
          
          analysisType = NA,
          quantificationType = NA,
          expressionsPath = NA,
          samplesPath = NA,
          configurationPath = NA,
          configuration = NA,
          
          initialize = function() {
            args <- commandArgs(trailingOnly=T)
            self$analysisType <- ifelse(is.na(args[1]), 'fluctuation', args[1])
            self$quantificationType <- ifelse(is.na(args[2]), 'byGene', args[2])
            self$expressionsPath <-
              ifelse(is.na(args[3]),
                     sprintf("out/%s/%s.csv",
                             self$quantificationType, self$analysisType),
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
            if(is.na(private$expressions))
              private$expressions <-
                read.table(self$expressionsPath, header=T, check.names=F,
                           sep=',', quote='', row.names=1)
            private$expressions
          },
          
          getFluctuations = function() {
            if(is.na(private$fluctuations)) {
              e <- self$getExpressions()
              private$fluctuations <- e[, 'fluctuation']
              names(private$fluctuations) <- rownames(e)
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
