#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# args: <input_deseq_file> <output_file <df>
# it also creates Rplots.pdf

df = as.numeric(args[3])

options(width=200)

# https://cran.r-project.org/web/packages/locfdr/vignettes/locfdr-example.pdf
# https://biodatascience.github.io/compbio/test/localfdr.html

require(locfdr)

#data(lfdrsim)
#zex = lfdrsim[,2]
#w = locfdr(zex)
#print(summary(zex))
#hist(zex,breaks=100)

###############################

data = read.table(args[1],head=T,sep='\t') # "deseq.SatS_KO_vs_MsmegWT"
# data already has genes as rownames
print(dim(data))

Tstats = data[,4]
Pvals = data[,5]
print("setting NAs to 0")
print(length(Tstats[is.na(Tstats)==T]))
Tstats[is.na(Tstats)==T] = 0

hist(Tstats,breaks=100)
hist(Pvals,breaks=100)

a = sort(Tstats)
qqnorm(a)
qqline(a)

require(MASS)
Tfit = fitdistr(Tstats,"t")
cat("Tfit:\n")
print(Tfit)
m = Tfit$estimate["m"]
s = Tfit$estimate["s"]
dfmle = Tfit$estimate["df"]

hist(Tstats,freq=F,breaks=100,main="T-distribution fitted to Tstats")
sample = m+s*rt(df=dfmle,n=2000)
dens = density(sample,bw=1.0)
lines(dens)

Z = qnorm(pt(Tstats,df)) # Z-transformation; use user's df arg, or dfmle from Tfit?
hist(Z,breaks=100)

Zfit = fitdistr(Z,"t")
cat("Zfit:\n")
print(Zfit)

a = sort(Z)
qqnorm(a)
qqline(a)

####################################

Zfdr = locfdr(Z,nulltype=0) # default=1 (empirical); 0 assumes null stats are ~N(0,1)
print(attributes(Zfdr))
print(Zfdr$Efdr)
cat("tail cutoffs of Z (fdr<0.2):\n")
print(Zfdr$z.2)
print(Zfdr$call)
plot(sort(Zfdr$fdr))

print(Zfdr$mat) # 120 bins
print(sum(Zfdr$fdr<0.2)) # 35/120 in the 2 tails

data$Z = Z
data$locfdr = Zfdr$fdr

####################################
# FDR control - optimize cutoffs

alpha = 0.05

temp = data
temp = temp[order(temp$locfdr),]

EFDR = NULL
for (i in 1:nrow(temp))
{
  efdr = mean(temp$locfdr[1:i])
  EFDR = rbind(EFDR,efdr)
}
temp$efdr = EFDR

print(head(temp,n=20))

data$efdr = temp[rownames(data),"efdr"]

Zhits = temp[which(temp$efdr<alpha),"Z"]
print(dim(Zhits))
a = max(Zhits[Zhits<0])
b = min(Zhits[Zhits>0])
na = sum(Zhits <= a)
nb = sum(Zhits >= b)
cat(sprintf("implied cutoffs for FDR<%s: Z<=%s (%s), Z>=%s (%s)\n",alpha,round(a,3),na,round(b,3),nb))

####################################

write.table(data,args[2],sep='\t',quote=F)

