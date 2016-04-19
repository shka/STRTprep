source('STRTprepSamples.R')
source('STRTprepExpressions.R')

library(R6)
library(yaml)

STRTprepHelper <- R6Class(
  "STRTprepHelper",
  
  private = list(
    base_path = NA,
    configuration = NA),
  
  public = list(
    comparison = NA,
    expressions_path = NA,
    options = NA,
    output_prefix = NA,
    quantification = NA,
    samples_path = NA,
    
    initialize =
      function(name = 'dummy',
               quantification = 'byGene',
               comparison = 'global',
               base_path = '.',
               expressions_path =
                 sprintf("%s/out/%s/diffexp.csv", base_path, quantification),
               samples_path = sprintf("%s/out/byGene/samples.csv", base_path),
               configuration_path = sprintf("%s/src/conf.yaml", base_path),
               output_prefix =
                 sprintf('%s/out/%s/plugin_%s_%s',
                         base_path, quantification, name, comparison),
               required_packages = c()) {
        self$quantification <- quantification
        self$comparison <- comparison
        private$base_path <- base_path
        self$expressions_path <- expressions_path
        self$samples_path <- samples_path
        private$configuration <- yaml.load_file(configuration_path)
        self$output_prefix <- output_prefix
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
              biocLite(
                package, ask=F, lib.loc=.libPaths()[1], lib=.libPaths()[1])
            }})
        self$options <-
          private$configuration[['PLUGINS']][[name]][[self$comparison]]
      }
  ),
  
  active = list(
    expressions = function() STRTprepExpressions$new(self),
    samples = function() STRTprepSamples$new(self),
    reference_genome =
      function() private$configuration$PREPROCESS$GENOMESPIKERIBO
  )
)

STRTprepHelper$newPlugin <- function(name, required_packages=c()) {
  args <- commandArgs(trailingOnly = T)
  quantification <- ifelse(is.na(args[1]), 'byGene', args[1])
  comparison <- ifelse(is.na(args[2]), 'global', args[2])
  STRTprepHelper$new(
    name=name,
    quantification=quantification,
    comparison=comparison,
    expressions_path=ifelse(is.na(args[3]),
                            sprintf("./out/%s/diffexp.csv", quantification),
                            args[3]),
    samples_path=
      ifelse(is.na(args[4]), './out/byGene/samples.csv', args[4]),
    configuration_path=ifelse(is.na(args[5]), './src/conf.yaml', args[5]),
    output_prefix=ifelse(is.na(args[6]),
                         sprintf('./out/%s/plugin_%s_%s',
                                 quantification, name, comparison),
                         args[6]),
    required_packages=required_packages
  )
}
