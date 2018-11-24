#!/usr/bin/env Rscript

####################################
# process args, read inputs

args = commandArgs(trailingOnly=TRUE)
options(width=200)
require(reshape)
require(pscl)

if (length(args)<5) { 
  print("usage: Rscript ZINB.R <prot_table> <combined_wigs_TTR> <samples-metadata> <output_fname> <comma_separated_list_of_conditions> [-no_sat_adjust] [-gene ORF_id]")
  quit("no") }

prot_table = args[1]
wigs = args[2]
metadata = args[3]
outf = args[4]

if ("-no_sat_adjust" %in% args) { print("found") }

genes = read.table(prot_table,sep='\t',quote="",stringsAsFactors=F)

data = read.table(wigs,sep='\t',head=T)
file = file(wigs,"r")
lines = readLines(file)
close(file)
filenames = c()
for (i in 1:length(lines)) {
  if (substr(lines[i],1,6)=='#File:') { filenames = c(filenames,substr(lines[i],8,nchar(lines[i]))) } }
colnames(data) = c("coord",filenames,"ORF")

metadata = read.table(metadata,head=T,sep='\t',stringsAsFactors=F)
groupnames = strsplit(args[5],",")[[1]]
ngroups = length(groupnames)
for (group in groupnames) { 
  if (!(group %in% metadata$Condition)) { print(paste("error:",group,"not found in samples metadata")); quit("no") } }

sat_adjust = T
if ("-no_sat_adjust" %in% args) { sat_adjust = F }
cat("saturation adjustment:",sat_adjust,"\n")

one_gene = NULL
if ("-gene" %in% args) {
  i = which(args %in% "-gene")
  one_gene = args[i+1] }

#####################################
# utility functions

# returns (coord, cnt, cond, dataset) for each observation
# depends on genes and data as globals

get_melted_counts = function(rv,groupnames) 
{
  start = genes[genes$V9==rv,2]
  end = genes[genes$V9==rv,3]
  strand = genes[genes$V9==rv,4]
  if (strand=="+") { end = end-3 } # ignore TA in stop codon, like resampling
  else { start = start+3 } 
  sites = data[data$coord>=start & data$coord<=end,] 
  if (length(sites$coord)==0) { return(NULL) }
  rows = c()
  for (group in groupnames) {
    for (j in 1:length(metadata$Filename)) {
      if (metadata[j,"Condition"]==group) {
        dataset = metadata[j,"Filename"]
        obs = sites[,c('coord',dataset)] # cols of melted are: coord, variable(=filename), value(=count)
        colnames(obs) = c('coord','cnt')
        obs$cond=group
        obs$dataset=dataset
        rows = rbind(rows,obs) } } }
  counts = data.frame(rows)
  return(counts) 
}

get_params = function(mod,counts) 
{
  factors = as.data.frame(colnames(counts)) # use colnames(mod.model)?
  colnames(factors) = c('cond')
  params = cbind(factors, 
   Counts = predict(mod, newdata = factors, type = "count"),
   Zeros = predict(mod, newdata = factors, type = "zero"))
  return(params) 
}

ZINB_signif = function(melted,sat_adjust) 
{
  require(pscl)
  melted$cnt = as.integer(melted$cnt) # note: ZINB requires that counts are integers
  melted$dataset = as.factor(melted$dataset)
  mod1 = tryCatch(
     { if (sat_adjust) { zeroinfl(cnt~0+cond+offset(log(NZmean))|0+cond+offset(logitZperc),data=melted,dist="negbin") }
                  else { zeroinfl(cnt~0+cond,data=melted,dist="negbin") } },
     error=function(err) { return(NULL) } )
print(summary(mod1))
  mod0 = tryCatch( # null model, independent of conditions
     { if (sat_adjust) { zeroinfl(cnt~1+offset(log(NZmean))|1+offset(logitZperc),data=melted,dist="negbin") }
                  else { zeroinfl(cnt~1,data=melted,dist="negbin") } },
     error=function(err) { return(NULL) } )
  if (is.null(mod1) | is.null(mod0)) { return(1) }
  df1 = attr(logLik(mod1),"df"); df0 = attr(logLik(mod0),"df") # should be (2*ngroups+1)-3
  pval = pchisq(2*(logLik(mod1)-logLik(mod0)),df=df1-df0,lower.tail=F) # alternatively, could use lrtest()
  # this gives same answer, but I would need to extract the Pvalue...
  #require(lmtest)
  #print(lrtest(mod1,mod0))
  return(pval) 
}

# compute stats like mean and NZmean for each condition

get_stats = function(melted) 
{
  tot = NULL 
  nonzeros = NULL 
  obs = NULL 
  for (cond in groupnames) {
    temp = melted[melted$cond==cond,]
    tot = c(tot,sum(temp$cnt))
    nonzeros = c(nonzeros,sum(temp$cnt!=0))
    obs = c(obs,length(temp$cnt))
  }
  stats = data.frame(tot,nonzeros,obs)
  rownames(stats) = groupnames
  stats$mean = round(stats$tot/stats$obs,1)
  stats$NZmean = round(stats$tot/max(1,stats$nonzeros),1)
  stats$NZperc = round(stats$nonzeros/stats$obs,3)
  return(stats)
}

# analyze a particular gene, given its Rv number

gene_variability = function(rv) 
{
  gene = genes[genes$V9==rv,8]
  if (! rv %in% genes$V9) { cat("gene not found\n"); quit("no") }

  melted = get_melted_counts(rv,groupnames)
  nTA = length(table(melted$coord))

  vals = c(Rv=rv,gene=gene,nTA=nTA)
  cat("\n-----------------------------------------\n")
  print(vals)
  flush.console()
  if (nTA<=1) { return(vals) }
  if (sum(melted$cnt)==0) { return(vals) } # skip gene if completely empty in all conds

  stats = get_stats(melted)
  cat("stats:\n")
  print(stats)

  means = apply(t(stats$mean),2,mean)
  names(means) = paste0("mean_",groupnames)
  vals = c(vals,means)

  NZmeans = t(stats$NZmean)
  NZpercs = t(stats$NZperc)
  names(NZmeans) = paste0("NZmean_",groupnames)
  names(NZpercs) = paste0("NZperc_",groupnames) 
  vals = c(vals,NZmeans,NZpercs)

  # do ANOVA

  anova = aov(cnt~cond,data=melted) 
  ANOVA_pval = summary(anova)[[1]][["Pr(>F)"]][1]
  print(paste("ANOVA_pval:",ANOVA_pval))
  vals = c(vals,ANOVA_pval=ANOVA_pval)
  #print(summary(anova))

  # add global NZmean and Zperc for each dataset to melted (for offsets in ZINB model, i.e. saturation adjustment)

  Zperc = NULL
  for (j in 2:(length(colnames(data))-1)) {
    Zperc[[colnames(data)[j]]] = sum(data[,j]==0)/length(data[,j]) }
  melted$Zperc = Zperc[melted$dataset]
  melted$logitZperc = log(melted$Zperc/(1-melted$Zperc))
  globalNZmean = NULL # for each dataset
  for (j in 2:(length(colnames(data))-1)) {
    cts = data[,j]; a = sum(cts); b = sum(cts>0)
    globalNZmean[[colnames(data)[j]]] = a/b }
  melted$NZmean = globalNZmean[melted$dataset]

  # run ZINB

  ZINB_pval = ZINB_signif(melted,sat_adjust)
  print(paste("ZINB_pval:",ZINB_pval))
  vals = c(vals,ZINB_pval=ZINB_pval)

  print(vals)
  return(vals)
}

###########################################
# main

if (length(one_gene)!=0) { 
  if (one_gene %in% genes$V8) { one_gene = genes[genes$V8==one_gene,9] } # convert gene name to ORFid (Rv)
  gene_variability(one_gene) # will check if ORFid exists
  quit("no") } 

# this is the driver that loops over all genes

res = NULL
for (rv in genes[,9]) { 
  x = gene_variability(rv)
  if (length(x)<4+3*ngroups) { x = c(x,rep(0,3+3*ngroups-length(x)),1,1) }
  res = rbind(res,x)
  # what a pain: rbind only uses col names from first gene, which might have had all 0's and thus no pval
  if ("ZINB_pval" %in% names(x)) { colnames(res) = names(x) } 
}

res = as.data.frame(res,stringsAsFactors=F)
ANOVA_pvals = res$ANOVA_pval
ANOVA_qvals = p.adjust(ANOVA_pvals,method="BH")
ZINB_pvals = res$ZINB_pval
ZINB_qvals = p.adjust(ZINB_pvals,method="BH")
res2 = res[,1:(length(colnames(res))-2)]
res2$ANOVA_pval = ANOVA_pvals
res2$ANOVA_qval = ANOVA_qvals
res2$ZINB_pval = ZINB_pvals
res2$ZINB_qval = ZINB_qvals
print(paste("writing",outf))
cat(c("# command-line: Rscript ZINB2.R",args,"\n"),file=outf)
write.table(res2,outf,append=T,sep='\t',quote=F,row.names=F)
