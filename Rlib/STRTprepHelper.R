source('STRTprepSamples.R')
source('STRTprepExpressions.R')

library(R6)
library(yaml)

STRTprepHelper <- R6Class(
  "STRTprepHelper",
  
  private = list(
    base_path = NA,
    configuration = NA,
    quantification = NA),
  
  public = list(
    comparison = NA,
    expressions_path = NA,
    options = NA,
    output_prefix = NA,
    samples_path = NA,
    
    initialize = function(name = 'dummy',
                          quantification = 'byGene',
                          comparison = 'global',
                          base_path = '.',
                          required_packages = c()) {
      args <- commandArgs(trailingOnly = T)
      private$quantification <- ifelse(is.na(args[1]), quantification, args[1])
      self$comparison <- ifelse(is.na(args[2]), comparison, args[2])
      private$base_path <- base_path
      self$expressions_path <- ifelse(
        is.na(args[3]),
        sprintf("%s/out/%s/diffexp.csv", base_path, private$quantification),
        args[3])
      self$samples_path <- ifelse(
        is.na(args[4]),
        sprintf("%s/out/byGene/samples.csv", base_path),
        args[4])
      private$configuration <- yaml.load_file(ifelse(
        is.na(args[5]), sprintf("%s/src/conf.yaml", base_path), args[5]))
      self$output_prefix <- ifelse(
        is.na(args[6]),
        sprintf('%s/out/%s/plugin_%s_%s',
                base_path, private$quantification, name, self$comparison),
        args[6])
      sapply(
        required_packages,
        function(package) {
          if(regexpr('/', package) > 0) {
            if(!is.element(strsplit(package, '[/@]')[[1]][2],
                           rownames(installed.packages()))) {
              library(devtools)
              install_github(package)
            }
          } else if(!is.element(package, rownames(installed.packages()))) {
            library(BiocInstaller)
            biocLite(package, ask=F)
          }})
      self$options <-
        private$configuration[['PLUGINS']][[name]][[self$comparison]]
    }),
  
  active = list(
    expressions = function() STRTprepExpressions$new(self),
    samples = function() STRTprepSamples$new(self)))
