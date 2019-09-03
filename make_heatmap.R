#!/usr/bin/env Rscript

makefname = function(name,dir) { if (dir=="") { name } else { sprintf("%s/%s",dir,name) } }

args = commandArgs(trailingOnly=TRUE)

# input: temp_LFCs.txt, rundir (optional)

rundir = ""
if (length(args)==2) { rundir = args[2] }

#data = read.table("temp_LFCs.txt",head=T,sep='\t')
data = read.table(makefname(args[1],rundir),head=T,sep='\t')
lfcs = as.matrix(data[,3:length(colnames(data))])
labels = paste(data$Rv,'/',data$Gene,sep="")
rownames(lfcs) = labels

library(corrplot)
library(RColorBrewer)
library(gplots)

# Red for down regulated, green for upregulated. 
redgreen <- function(n) { c( hsv(h=0/6, v=seq(1,0,length=n/2) ), hsv(h=2/6, v=seq(0,1,length=n/2) ), ) }
colors <- colorRampPalette(c("red", "white", "green"))(n = 1000)

C = length(colnames(lfcs))
R = length(rownames(lfcs))
W = 300+C*30
H = 300+R*15

fname = "heatmap.png"
if (rundir!="") { fname = sprintf("%s/%s",rundir,fname) }
png(fname,width=W,height=H)
#par(oma=c(10,1,1,10)) # b,l,t,r; units="lines of space"
heatmap.2(lfcs,col=colors,margin=c(12,12),lwid=c(1,8),lhei=c(1,8),trace="none",cexCol=1.4,cexRow=1.4) # make sure white=0
dev.off()

