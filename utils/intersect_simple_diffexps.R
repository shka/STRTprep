#!/usr/bin/env Rscript

source('Rlib/STRTprepHelper.R', chdir=T)
args <- commandArgs(trailingOnly = T)
table1 <- STRTprepHelper$new(
  name='simple_diffexp',
  quantification=args[1],
  comparison=args[2])$expressions$significant$table[, 1:5]
table2 <- STRTprepHelper$new(
  name='simple_diffexp',
  quantification=args[1],
  comparison=args[3])$expressions$significant$table[, 1:5]
targets <- intersect(rownames(table1), rownames(table2))
tmp <- cbind(table1[targets, ], table2[targets, ])
colnames(tmp) <- c(sprintf("%s.%s", colnames(table1), args[2]),
                   sprintf("%s.%s", colnames(table2), args[3]))
write.csv(tmp, quote=F)
