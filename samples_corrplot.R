#!/usr/bin/env Rscript

makefname = function(name,dir) { if (dir=="") { name } else { sprintf("%s/%s",dir,name) } }

args = commandArgs(trailingOnly=TRUE)

# input: file with insertion counts, rundir (optional)

rundir = ""
if (length(args)==2) { rundir = args[2] }

#data = read.table("temp_counts_TTR.txt",sep='\t',head=T)
data = read.table(makefname(args[1],rundir),sep='\t',head=T)
vals = as.matrix(data[,4:length(colnames(data))])
N = length(colnames(vals))

library(corrplot)

fname = "samples_corrplot.png"
if (rundir!="") { fname = sprintf("%s/%s",rundir,fname) }
png(fname,width=300+20*N,height=300+20*N)
#corrplot(cor(vals)) # among TAsites (not good)

TAsites = aggregate ( coord~ORF,data=data,length)
colnames(TAsites) = c("ORF","sites")
temp = aggregate(vals,by=list(data$ORF),FUN=mean)
corrplot(cor(temp[TAsites$sites>=10,2:length(colnames(temp))]))

dev.off()


