library(xlsx)

args <- commandArgs(trailingOnly=T)

wb.template.file <- args[1]
wb.output.file <- args[2]
spike.file <- args[3]
demultiplex.file <- args[4]
alignment.file <- args[5]
annotation.file <- args[6]
sequence.files <- args[7]

wb <- loadWorkbook(wb.template.file)
sheets <- getSheets(wb)

tbl.spike <- read.table(spike.file, header=T, sep="\t", quote='', check.names=F, fill=T)
addDataFrame(tbl.spike, sheets[['spike']], col.names=F, row.names=F, startRow=2, startColumn=1); 

tbl.demultiplex <- read.table(demultiplex.file, header=T, sep="\t", quote='', check.names=F, fill=T)
addDataFrame(tbl.demultiplex, sheets[['demultiplex']], col.names=F, row.names=F, startRow=2, startColumn=1); 

tbl.alignment <- read.table(alignment.file, header=T, sep="\t", quote='', check.names=F, fill=T)
addDataFrame(tbl.alignment, sheets[['alignment']], col.names=F, row.names=F, startRow=2, startColumn=1); 

tbl.annotation <- read.table(annotation.file, header=T, sep="\t", quote='', check.names=F, fill=T)
addDataFrame(tbl.annotation, sheets[['annotation']], col.names=F, row.names=F, startRow=2, startColumn=1); 

sequences.sheet <- sheets[['summary']]
sequences.rows <- getRows(sequences.sheet, rowIndex=15)
sequences.cells <- getCells(sequences.rows, colIndex=2)
print(getCellValue(sequences.cells[['15.2']]))
setCellValue(sequences.cells[['15.2']], sequence.files)

saveWorkbook(wb, wb.output.file)
