#!/usr/bin/env Rscript

makefname = function(name,dir) { if (dir=="") { name } else { sprintf("%s/%s",dir,name) } }

args = commandArgs(trailingOnly=TRUE)

# input: file with LFCS, rundir (optional)

rundir = ""
if (length(args)==2) { rundir = args[2] }

#data = read.table("temp_LFCs.txt",sep='\t',head=T)
data = read.table(makefname(args[1],rundir),sep='\t',head=T)
vals = as.matrix(data[,3:length(colnames(data))])
N = length(colnames(vals))

library(corrplot)

png(makefname("conditions_corrplot.png",rundir),width=20*N+300,height=20*N+300)
corrplot(cor(vals))
dev.off()


