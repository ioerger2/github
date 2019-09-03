#!/usr/bin/env Rscript

makefname = function(name,dir) { if (dir=="") { name } else { sprintf("%s/%s",dir,name) } }

args = commandArgs(trailingOnly=TRUE)

K = as.integer(args[2]) # gene clusters
D = as.integer(args[3]) # varimax dimensions (number of condition groups)
rundir = ""
if (length(args)==4) { rundir = args[4] }

# input: temp_LFCs.txt
data = read.table(makefname(args[1],rundir),head=T,sep='\t')
lfcs = data[,3:length(colnames(data))]
labels = paste(data$Rv,'/',data$Gene,sep="")
rownames(lfcs) = labels

R = length(rownames(lfcs)) # genes
C = length(colnames(lfcs)) # conditions

library(corrplot)

fname = makefname("lfcs_boxplot.png",rundir)
png(fname,width=300+15*C,height=500)
boxplot(lfcs,las=2,ylab="log2 fold change of genes (relative to the median)")
dev.off()

km = kmeans(lfcs,K)
clusters = km$cluster

write.table(data.frame(lfcs,clust=clusters),makefname("temp_clust.txt",rundir),sep='\t',quote=F)

pca = prcomp(t(lfcs))

library(MASS)
mat = as.matrix(lfcs)
ldapca = lda(clusters~.,lfcs)
X = mat %*% ldapca$scaling[,1]
Y = mat %*% ldapca$scaling[,2]
fname = makefname("pca_genes.png",rundir)
png(fname,width=1000,height=1000)
plot(X,Y,pch=20)
text(X,Y,label=labels,adj=c(0.5,1.5),col=clusters)
dev.off()

actpca = prcomp(lfcs,center=TRUE,scale=TRUE)
fname = makefname("pca_conditions.png",rundir)
png(fname,width=500,height=500)
plot(actpca$rotation[,1:2],pch=20)#,xlim=c(-1,1),ylim=c(-1,1))
text(actpca$rotation[,1:2],label=colnames(lfcs),adj=c(0.5,1.5),cex.lab=1.5)
dev.off()

hc = hclust(dist(lfcs),method="ward.D2")
fname = makefname("hclust_genes.png",rundir)
png(fname,width=300+15*R,height=600)
plot(hc)
K2 = max(2,min(K,max(hc$height)))
rect.hclust(hc,k=K2,border='red')
dev.off()

hc = hclust(dist(t(lfcs)),method="ward.D2")
fname = makefname("hclust_conditions.png",rundir)
png(fname,width=300+15*C,height=600)
plot(hc)
D2 = max(2,min(D,C-1))
rect.hclust(hc,k=D2,border='red')
dev.off()

fname = makefname("cluster_opt.png",rundir)
png(fname)
library(factoextra)
fviz_nbclust(lfcs, kmeans, method = "wss",k.max=30) # or silhouette or gap_stat
#fviz_nbclust ( lfcs, kmeans, method = "wss",k.max=30)+geom_hline ( yintercept=52 )+geom_vline(xintercept=10)
dev.off()

####################################################
# required factoextra and corrplot libraries

var = get_pca_var(actpca) 
fname = makefname("condition_PCs.png",rundir)
png(fname,width=1000,height=1000)
corrplot(var$cor,main="Principle Components",mar=c(1,1,1,1))
dev.off()

S = diag(actpca$sdev,D,D)
rawLoadings = actpca$rotation[,1:D] %*% S
vmax = varimax(rawLoadings)
rotatedLoadings = vmax$loadings
rotatedScores = scale(actpca$x[,1:D]) %*% vmax$rotmat

# https://stats.stackexchange.com/questions/59213/how-to-compute-varimax-rotated-principal-components-in-r
# scores = U.S = X.V = actpca$x = scale(lfcs) %*% actpca$rotation
combined_transform = (actpca$rotation[,1:D] %*% S) %*% vmax$rotmat
# vmax$loadings is same as combined_transform, but with abs(vals)<0.1 blanked out
write.table(round(combined_transform,6),makefname("varimax_loadings.txt",rundir),sep='\t',quote=F)

fname = makefname("varimax.png",rundir)
png(fname,width=800,height=800)
corrplot(rotatedLoadings,main="Varimax loadings",mar=c(1,1,1,1)) 
dev.off()

scores = rotatedScores
colnames(scores) = paste("score",seq(1:D),sep="")
squares = scores*scores
ss = apply(squares,1,sum)
cos2 = squares/ss
best = as.matrix(apply(cos2,1,max))
assoc = cos2*sign(scores)
colnames(assoc) = paste("assoc",seq(1:D),sep="")

# # generate null distribution for cos2 to closest axis
# 
# library(MASS)
# SAMPLES = 10000
# sample = mvrnorm(n=SAMPLES,mu=rep(0,D),Sigma=diag(D))
# # no need to apply Varimax rotation
# squares = sample*sample
# ss = apply(squares,1,sum)
# Xcos2 = squares/ss
# Xbest = as.matrix(apply(Xcos2,1,max))
# #distn = ecdf(best)
# pvals = apply(best,1,function (x) { length(Xbest[Xbest>=x]) })
# pvals = pvals/SAMPLES
# padj = p.adjust(pvals,method="BH")
 
#res = data.frame(data,scores,assoc,best,pvals,padj)
res = data.frame(data,scores,assoc)
write.table(format(res,digits=3,scientific=F),makefname("temp_scores.txt",rundir),sep='\t',quote=F,row.names=F)
