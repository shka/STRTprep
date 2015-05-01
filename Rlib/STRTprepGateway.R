library(R6)
library(yaml)

## TODO: private fields, as many as possible
## TODO: field access '$...' by active bindings rather than '$get...()'

STRTprepSamples <-
  R6Class(
      "STRTprepSamples",

      public = list(
          comparisonClass = NA,
          configuration = NA,
          content = NA,
          path = NA,
          
          initialize = function(samplesPath, comparisonClass, configuration) {
            self$path <- samplesPath
            self$comparisonClass <- comparisonClass
            self$configuration <- configuration
            tmp <-
              read.table(self$path, header=T, check.names=F, sep=',', quote='')
            rownames(tmp) <- tmp[, which(colnames(tmp) == 'NAME')]
            if(comparisonClass != 'global') {
              className <- sprintf('CLASS.%s', comparisonClass)
              blockName <- sprintf('BLOCK.%s', comparisonClass)
              classification <- tmp[, className]
              tmp <- tmp[which(!is.na(classification)), ]
              tmp2 <- colnames(tmp)
              tmp2[which(tmp2 == className)] <- 'CLASS'
              tmp2[which(tmp2 == blockName)] <- 'BLOCK'
              colnames(tmp) <- tmp2
            }
            self$content <- tmp
          },

          getNames = function() rownames(self$content),

          getClasses = function() self$content[, 'CLASS'],

          getBlocks = function() self$content[, 'BLOCK']))

STRTprepExpressions <-
  R6Class(
      "STRTprepExpressions",

      public = list(
          comparisonClass = NA,
          configuration = NA,
          content = NA,
          path = NA,
          samples = NA,

          initialize = function(expressionsPath, samples) {
            self$path <- expressionsPath
            self$samples <- samples
            self$configuration <- samples$configuration
            self$comparisonClass <- samples$comparisonClass
            e <- read.table(self$path, header=T, check.names=F,
                            sep=',', quote='', row.names=1)
            tmp <- private$extractNormalizedReads(e)
            tmp <- tmp[, which(is.element(colnames(tmp), samples$getNames()))]
            tmp2 <- colnames(e)
            cls <- self$comparisonClass
            if(cls != 'global') {
              tmp2[which(tmp2 == sprintf('diffexpScore.%s', cls))] <-
                'diffexpScore'
              tmp2[which(tmp2 == sprintf('pvalue.%s', cls))] <- 'pvalue'
              tmp2[which(tmp2 == sprintf('qvalue.%s', cls))] <- 'qvalue'
            }
            tmp2[which(tmp2 == sprintf('fluctuationScore.%s', cls))] <-
              'fluctuationScore'
            tmp2[which(tmp2 == sprintf('fluctuation.%s', cls))] <- 'fluctuation'
            colnames(e) <- tmp2
            self$content <- e[which(rowSums(tmp) > 0), ]
          },

          getFluctuations = function() {
            if(!is.data.frame(private$fluctuations)) {
              tmp <- self$content[, c('fluctuation', 'fluctuationScore')]
              private$fluctuations <- tmp[which(!is.na(tmp[, 1])), ]
            }
            private$fluctuations
          },

          getNormalizedReads = function(excludeSpikeins=TRUE) {
            tmp <- private$extractNormalizedReads(self$content)
            if(excludeSpikeins == TRUE)
              tmp <- tmp[which(substr(rownames(tmp), 1, 10) != 'RNA_SPIKE_'), ]
            as.matrix(tmp[,  which(is.element(colnames(tmp),
                                              self$samples$getNames()))])
          },

          significant = function() {
            tmp <- STRTprepExpressions$new(self$path, self$samples)
            if(self$comparisonClass != 'global')
              tmp$content <-
                self$content[which(self$content[, 'fluctuation'] <
                                       self$configuration$DEFAULT$FLUCTUATION
                                   & self$content[, 'qvalue'] <
                                       self$configuration$DEFAULT$DIFFEXP), ]
            else
              tmp$content <- 
                self$content[which(self$content[, 'fluctuation'] <
                                       self$configuration$DEFAULT$FLUCTUATION),]
            tmp
          }),

      private = list(
          fluctuations = NA,
          
          extractNormalizedReads = function(e) {
            tmp <- e[, which(substr(colnames(e), 1, 2) == 'N|')]
            colnames(tmp) <-
              sapply(colnames(tmp),
                     function(n) { strsplit(n, '\\|')[[1]][3] })
            tmp
          }))

STRTprepGateway <-
  R6Class(
      "STRTprepGateway",
      
      public = list(
          quantificationType = NA,
          expressionsPath = NA,
          comparisonClass = NA,
          samplesPath = NA,
          configurationPath = NA,
          outputPrefix = NA,
          configuration = NA,
          
          initialize = function(
              pluginName,
              requiredPackages=c(),
              quantificationType='byGene',
              comparisonClass='global') {
            args <- commandArgs(trailingOnly=T)
            self$quantificationType <-
              ifelse(is.na(args[1]), quantificationType, args[1])
            self$comparisonClass <-
              ifelse(is.na(args[2]), comparisonClass, args[2])
            self$expressionsPath <-
              ifelse(is.na(args[3]),
                     sprintf("out/%s/diffexp.csv", self$quantificationType),
                     args[3])
            self$samplesPath <-
              ifelse(is.na(args[4]),
                     sprintf("out/byGene/samples.csv", self$quantificationType),
                     args[4])
            self$configurationPath <-
              ifelse(is.na(args[5]), 'conf.yaml', args[5])
            self$outputPrefix <-
              ifelse(is.na(args[6]),
                     sprintf('out/%s/plugin_%s_%s',
                             self$quantificationType,
                             pluginName,
                             self$comparisonClass),
                     args[6])
            self$configuration <- yaml.load_file(self$configurationPath)

            sapply(requiredPackages,
                   function(requiredPackage) {
                     if(!is.element(requiredPackage,
                                    rownames(installed.packages()))) {
                       library(BiocInstaller)
                       biocLite(requiredPackage, ask=F)
                     }
                   })
          },

          getExpressions = function() {
            if(!is.data.frame(private$expressions))
              private$expressions <-
                STRTprepExpressions$new(self$expressionsPath, self$getSamples())
            private$expressions
          },

          getSamples = function() {
            if(!is.data.frame(private$samples))
              private$samples <-
                STRTprepSamples$new(self$samplesPath,
                                    self$comparisonClass,
                                    self$configuration)
            private$samples
          }),

      private = list(
          expressions = NA,
          samples = NA))
