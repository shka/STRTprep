#!/usr/bin/env Rscript

source('Rlib/STRTprepHelper.R', chdir=T)
helper <- STRTprepHelper$new(name='simple_gsea',
                             required_packages=c('renozao/pkgmaker@develop',
                                                 'renozao/NMF',
                                                 'XML'))

if(helper$quantification != 'byGene') {
  warning('Skipped - this plugin is only for gene-based quantification.')
  q(save='no')
}
options <- helper$options
prefix <- helper$output_prefix

nreads <- helper$expressions$offset$significant$normalized_levels
if(nrow(nreads) < 5) {
  warning('Skipped - fluctuated TFEs were less than 5.')
  q(save='no')
}
significant <- rownames(nreads)

background <- rownames(helper$expressions$table)
background <- background[which(substr(background, 1, 10) != 'RNA_SPIKE_'
                               & substr(background, 1, 5) != 'RIBO_')]

if(helper$comparison == 'global') {
  warning('Skipped - this plugin is not for global comparison.')
  q(save='no')
}

annotations <-
  helper$samples$annotations[, union(options$ANNOTATIONS, 'CLASS')]
if(is.null(options$LABELS)) {
  annotations[, 'CLASS'] <- helper$samples$annotations[, 'CLASS']
} else {
  annotations[, 'CLASS'] <-
    factor(helper$samples$annotations[, 'CLASS'], labels=options$LABELS)
}

library(NMF)
library(parallel)
library(XML)

msigdb <- xmlRoot(xmlTreeParse(options$MSIGDB))
xmlAttrs(msigdb)

geneSets <- xmlSApply(msigdb, function(geneset) strsplit(xmlAttrs(geneset)['MEMBERS_SYMBOLIZED'], ','))
names(geneSets) <- xmlSApply(msigdb, function(geneset) paste(xmlAttrs(geneset)[c('CATEGORY_CODE', 'SUB_CATEGORY_CODE', 'STANDARD_NAME')], collapse='/'))

test_gse <- function(id) {
  tmp <- strsplit(id, '/')[[1]]
  geneSet <- geneSets[[id]]
  notSignificant <- setdiff(background, significant)
  notGeneSet <- setdiff(background, geneSet)
  stat <- c(length(intersect(geneSet, significant)),
            length(intersect(geneSet, notSignificant)),
            length(intersect(notGeneSet, significant)),
            length(intersect(notGeneSet, notSignificant)))
  test <- chisq.test(matrix(stat, ncol=2, byrow=2))
  stat[5] <- ifelse(is.nan(test$p.value), NA, test$p.value)
  ## stat <- c(tmp, stat, test$expected)
  stat <- c(stat, test$expected)
  stat
}
tmp <- mclapply(names(geneSets), test_gse, mc.cores=detectCores())

tbl <- data.frame(matrix(unlist(tmp), ncol=9, byrow=T), stringsAsFactors=F)
colnames(tbl) <- c('OBS1', 'OBS2', 'OBS3', 'OBS4', 'CHISQP', 'EXP1', 'EXP2', 'EXP3', 'EXP4')
tbl[, 'CATEGORY'] <- sapply(names(geneSets)[1:length(tmp)], function(id) strsplit(id, '/')[[1]][1])
tbl[, 'SUBCATEGORY'] <- sapply(names(geneSets)[1:length(tmp)], function(id) strsplit(id, '/')[[1]][2])
tbl[, 'GENESET'] <- sapply(names(geneSets)[1:length(tmp)], function(id) strsplit(id, '/')[[1]][3])
tbl[, 'ENRICHMENTP'] <- p.adjust(as.numeric(as.character(tbl[, 'CHISQP'])), method='BH')
tbl[, 'RESULT'] <- ifelse((tbl[, 'ENRICHMENTP'] < options$ENRICHMENTP
                           & tbl[, 'OBS1'] > tbl[, 'EXP1']
                           & tbl[, 'EXP1'] > 5), 'OR', 'NA')
rownames(tbl) <- names(geneSets)[1:nrow(tbl)]
write.csv(tbl[, c('CATEGORY', 'SUBCATEGORY', 'GENESET', 'RESULT', 'ENRICHMENTP',
                  'OBS1', 'OBS2', 'OBS3', 'OBS4',
                  'EXP1', 'EXP2', 'EXP3', 'EXP4',
                  'CHISQP')],
          sprintf('%s.csv', prefix), quote=F)

##

sigs <- length(which(tbl[, 'RESULT'] == 'OR'))
sprintf('... %d out of %d gene-sets were enriched.', sigs, nrow(tbl))

nmf.options(grid.patch=T)
distfun <- function(x) as.dist((1-cor(t(x), method='spearman'))/2)
clustfun <- function(d) hclust(d, method='ward.D2')

plot_heatmaps <- function(id) {
  stat <- tbl[id, ]
  if(stat['RESULT'] == 'OR') {
    dir <- sprintf('%s.figures/%s/%s', prefix, stat['CATEGORY'], stat['SUBCATEGORY'])
    if(!dir.exists(dir))
      dir.create(dir, recursive=T)
    geneSet <- geneSets[[id]]
    targets <- intersect(geneSet, significant)
    heatmap <- aheatmap(
      log10(nreads[targets, ]),
      hclustfun=clustfun, distfun=distfun,
      scale='row', Colv=order(annotations[, 'CLASS']),
      annCol=annotations, layout='dlmL|alm', main=id,
      fontsize=6, filename=sprintf('%s/%s.pdf', dir, stat['GENESET']))
    save(heatmap,
         file=sprintf('%s/%s.RData', dir, stat['GENESET']), compress='gzip')
  }
}

unlink(sprintf('%s.figures', prefix), recursive=T, force=T)
tmp2 <- mclapply(rownames(tbl)[1:nrow(tbl)], plot_heatmaps, mc.cores=detectCores())
