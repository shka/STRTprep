#!/usr/bin/env Rscript

source('Rlib/STRTprepHelper.R', chdir=T)
helper <- STRTprepHelper$new(name='simple_diffexp')
table <- helper$expressions$significant$table
if(nrow(table) > 0) {
  write.csv(table, file=sprintf('%s.csv', helper$output_prefix), quote=F)
} else
  warning('Skipped - no significant features.')
  
