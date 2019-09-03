import sys,os,math,random
import numpy,scipy
import statsmodels.stats.multitest
import statsmodels.stats.multicomp
sys.path.append("/pacific/home/ioerger/transit/src/")
import pytransit.norm_tools as norm_tools
import pytransit.tnseq_tools as tnseq_tools
sys.path.append("/pacific/home/ioerger/genomics")
from utils import read_genome,read_genes,hash_genes
from Spreadsheet import *

DIR = "/pacific/home/ioerger/TRASH/macrophages" # location of R scripts

EOL = "\n"

def quantiles(X,Q):
  temp,n = sorted(X),len(X)
  return [temp[min(n-1,int(n*q))] for q in Q]

def max_diff_of_means(countsvec):
  means = [numpy.mean(x) for x in countsvec]
  a,b = min(means),max(means)
  return b-a

def get_range_null_distribution(countsvec,iters=10000): 
  sizes = [x.size for x in countsvec]
  start,ranges = 0,[]
  for n in sizes:
    ranges.append((start,start+n))
    start += n
  allcounts,mdms = [],[]
  for x in countsvec: allcounts += x.tolist()
  for i in range(iters):
    random.shuffle(allcounts)
    vec = [allcounts[x[0]:x[1]] for x in ranges]
    mdm = max_diff_of_means(vec)
    mdms.append(mdm)
  mdms.sort(reverse=True)
  return mdms

def multiway_resampling(countsvec): # returns a p-value
  samples = 10000
  distn = get_range_null_distribution(countsvec,iters=samples)
  mdm = max_diff_of_means(countsvec)
  i = 0
  while i<samples and mdm<=distn[i]: i += 1 # a binary search would be more efficient
  return i/float(samples)

def adaptive_multiway_resampling(countsvec,iters=10000): # returns a p-value
  obsmdm = max_diff_of_means(countsvec)
  sizes = [x.size for x in countsvec]
  start,ranges = 0,[]
  for n in sizes:
    ranges.append((start,start+n))
    start += n
  allcounts,cnt = [],0
  for x in countsvec: allcounts += x.tolist()
  checkpoint = 100
  mdms,BINSIZE = {},10 # for histograms
  pval = None
  for i in range(iters):
    if i==checkpoint:
      if cnt>10: pval = cnt/float(i); break
      checkpoint *= 10
    random.shuffle(allcounts)
    vec = [allcounts[x[0]:x[1]] for x in ranges]
    mdm = max_diff_of_means(vec)
    if mdm>=obsmdm: cnt += 1

    bin = int(mdm/BINSIZE)
    if bin not in mdms: mdms[bin] = 0
    mdms[bin] += 1
  print obsmdm
  a,b = min(mdms),max(mdms)  
  for i in range(40):
    print i*BINSIZE,mdms.get(i,0)

  if pval==None: pval = cnt/float(iters)
  return pval

#######################################

def read_combined_wig(fname):
  sites,counts,files = [],[],[]
  for line in open(fname):
    if line.startswith("#File:"):
      files.append(line.split()[1])
    if line[0]=='#': continue
    w = line.split()
    w = w[:len(files)+1] # chop off suffix, like "Rv gene"
    w = [int(float(x)) for x in w]
    sites.append(w[0])
    counts.append(w[1:])
  return (numpy.array(counts),numpy.array(sites),files)

def makefname(name,dir):
  if dir=="": return name
  else: return "%s/%s" % (dir,name)

def pipeline(ref,annot,combined_wig,samples_metadata_fname,conditions_metadata_fname,K=10,D=5,debatch=False,stdnorm=False,rundir=""):

  Samples = Spreadsheet(makefname(samples_metadata_fname,rundir))
  Conditions = Spreadsheet(makefname(conditions_metadata_fname,rundir))

  genome = read_genome(makefname(ref,rundir))
  genes = read_genes(makefname(annot,rundir))
#  temp = [] ###
#  for (start,end,rv,gene,strand) in genes:
#    if rv not in "Rv0181c Rv3842c Rv0877".split(): temp.append((start,end,rv,gene,strand))
#    #if rv=="Rvnr01" or int(rv[2:6]) not in range(2930,2941): temp.append((start,end,rv,gene,strand))
#    else: print 'excluding',rv
#  genes = temp
  ghash = hash_genes(genes,genome)

  datasets = Samples.getcol("Filename") # filenames

  #(data, sites) = tnseq_tools.get_data(datasets) # data: rows are datasets, cols are TA sites
  (data, sites, filenamesInCombWig) = read_combined_wig(makefname(combined_wig,rundir)) # data: rows are TA sites, cols are datasets

  # "datasets" are columns (indexes) in combined wig (data)
  # SampleIndexes are list of rows in Samples to be analyzed (contains Filenames and Conditions), i.e. which datasets (row indexes in data) represent each condition?
  # wigindexes are corresponding columns in data (from original order, matched by Filename)

  SampleIndexes,wigindexes = [],[]
  for cond in Conditions.keys:
    for row in range(Samples.nrows):
      if Samples.get(row,"Condition")==cond:
        fname = Samples.get(row,"Filename")
        if fname not in filenamesInCombWig: 
          print "error: filename '%s' listed in samples metadata not found in combined wig file" % fname; sys.exit(0)
        SampleIndexes.append(row)
        wigindexes.append(filenamesInCombWig.index(fname))

  # reordering data matrix (select columns, order by cond) for normalization; now should be parallel to SampleIndexes
  data = data[:,wigindexes]
  neworder = [filenamesInCombWig[i] for i in wigindexes]
  filenamesInData = neworder

  # foreach cond, list of indexes into data
  cond2datasets = {}
  for cond in Conditions.keys:
    cond2datasets[cond] = []
    for row in range(Samples.nrows):
      if Samples.get(row,"Condition")==cond:
        fname = Samples.get(row,"Filename")
        cond2datasets[cond].append(filenamesInData.index(fname))
    if len(cond2datasets[cond])==0:
      print "error: no samples found in metadata for condition %s" % cond; sys.exit(0)

  (normed,factors) = norm_tools.normalize_data(data.transpose(), method='TTR') 
  Nds,Nsites = normed.shape # transposed: rows are datasets, cols are TA sites

  file = open(makefname("temp_counts_TTR.txt",rundir),"w")
  vals = "coord ORF gene".split()+[Samples.get(r,"Id") for r in SampleIndexes] # in order of Conditions
  file.write('\t'.join(vals)+EOL)
  for i,co in enumerate(sites):
    rv,gene = "igr","igr"
    if co in ghash:
      annot = ghash[co]
      rv,gene = annot[2],annot[3]
    vals = [str(co),rv,gene]+["%0.1f" % (int(x)) for x in list(normed[:,i])]
    file.write('\t'.join([str(x) for x in vals])+EOL)
  file.close()

  os.system("Rscript %s/samples_corrplot.R temp_counts_TTR.txt %s > /dev/null" % (DIR,rundir))

  # write table of stats (saturation,NZmean)
  file = open(makefname("temp_stats.txt",rundir),"w")
  file.write("dataset\tdensity\tmean_ct\tNZmean\tNZmedian\tmax_ct\ttotal_cts\tskewness\tkurtosis\n")
  for i in range(data.shape[1]):
    density, meanrd, nzmeanrd, nzmedianrd, maxrd, totalrd, skew, kurtosis = tnseq_tools.get_data_stats(data[:,i])
    vals = [filenamesInData[i], "%0.2f" % density, "%0.1f" % meanrd, "%0.1f" % nzmeanrd, "%d" % nzmedianrd, maxrd, totalrd, "%0.1f" % skew, "%0.1f" % kurtosis]
    file.write('\t'.join([str(x) for x in vals])+EOL)
  file.close()

  ################################

  # new version, hashes on Rv

  siteshash = {}
  for i,TA in enumerate(sites): siteshash[TA] = i

  TAsites = {}
  for g,(start,end,Rv,gene,strand) in enumerate(genes):
    siteindexes = []
    for i in range(start,end): # end+1?
      co = i+1
      if co in siteshash: siteindexes.append(siteshash[co]) 
    TAsites[Rv] = siteindexes



  if False: # print means for each gene for each dataset
    print '\t'.join(filenamesInData)
    print '\t'.join([Samples.get(x,"Id") for x in SampleIndexes])
    for (start,end,Rv,gene,strand) in genes:
      if len(TAsites[Rv])>0: # skip genes with no TA sites
        means = []
        for j in range(Nds):
          obs = normed[j,TAsites[Rv]]
          means.append(numpy.mean(obs))
      print '\t'.join([Rv,gene,str(len(TAsites[Rv]))]+["%0.1f" % x for x in means])
    #sys.exit(0)


  Counts = {} # sub-arrays of (datasets X sites) for genes X conds, where #TAs>0
  for (start,end,Rv,gene,strand) in genes:
    siteindexes = TAsites[Rv]
    local = normed[:,siteindexes]
    counts = []
    for j in range(Conditions.nrows):
      cond = Conditions.get(j,"Condition")
      obs = local[cond2datasets[cond],:]
      counts.append(obs) # sub-matrices
    Counts[Rv] = counts
  
  Means = {} 
  for (start,end,Rv,gene,strand) in genes:
   if len(TAsites[Rv])>0: # skip genes with no TA sites
    means = []
    for j in range(Conditions.nrows):
      cond = Conditions.get(j,"Condition")
      obs = Counts[Rv][j]
      means.append(numpy.mean(obs))
    Means[Rv] = means
  
  file = open(makefname("temp_gene_means.txt",rundir),"w")
  file.write('\t'.join("ORF Gene TAs".split()+Conditions.getcol('Condition'))+EOL)
  for (start,end,Rv,gene,strand) in genes:
    if Rv not in Means: continue # skip genes with no TA sites
    means,sites = Means[Rv],TAsites[Rv]
    file.write('\t'.join([Rv,gene,str(len(sites))]+["%0.1f" % x for x in means])+EOL)
  file.close()
  
  ################################
  # do ANOVA to identify genes with significant variability
  # it is much faster to do this in python than R (and don't have to create the melted file)
  
  pvals,Rvs, = [],[] # lists
  for (start,end,Rv,gene,strand) in genes:
   if Rv in Means:
    countsvec = Counts[Rv]
    countsvec = [x.flatten() for x in countsvec]
    stat,pval = scipy.stats.f_oneway(*countsvec)
    pvals.append(pval)
    Rvs.append(Rv)
  
  pvals = numpy.array(pvals)
  mask = numpy.isfinite(pvals)
  qvals = numpy.full(pvals.shape,numpy.nan)
  qvals[mask] = statsmodels.stats.multitest.fdrcorrection(pvals[mask])[1] # BH, alpha=0.05

  p,q = {},{}
  for i,rv in enumerate(Rvs):
     p[rv],q[rv] = pvals[i],qvals[i]
  pvals,qvals = p,q, # hashes on Rv

  # post-hoc analysis to identify "high" and "low" subgroups
  Sortedmeans,Divisions = {},{} # sortedmeans: list of (mean,cond); divisions: list of (lowsubset,highsubset) or None
  for (start,end,Rv,gene,strand) in genes:
    if Rv in qvals and qvals[Rv]<0.05:
      #print Rv,gene
      countsvec = Counts[Rv]
      countsvec = [x.flatten() for x in countsvec]
      allcounts,allconds = numpy.array([]),numpy.array([])
      for j in range(len(Conditions.keys)):
        allcounts = numpy.append(allcounts,countsvec[j])
        for k in range(len(countsvec[j])):
          allconds = numpy.append(allconds,Conditions.keys[j])
      mc = statsmodels.stats.multicomp.MultiComparison(allcounts,allconds)
      tuk = mc.tukeyhsd() # alpha=0.1
      #print tuk
      reject,n = {},0
      groups = tuk.groupsunique.tolist()
      for j,group1 in enumerate(groups):
        for k,group2 in enumerate(groups):
          if j<k: reject[(group1,group2)] = reject[(group2,group1)] = tuk.reject[n]; n += 1
      sortedmeans = [(numpy.mean(countsvec[k]),Conditions.keys[k]) for k in range(len(countsvec))] # list of (mean,condname)
      sortedmeans = sorted(sortedmeans)
      Sortedmeans[Rv] = sortedmeans
      sortedgroups = [x[1] for x in sortedmeans]
      #for (m,c) in sortedmeans: print c,"%0.1f" % m
      candidate_divisions = [] # lists of condition names in 2 subgroups (lower and higher, but not necessarily in that order)
      for j in range(1,len(sortedmeans)):
        lowsubset,highsubset = sortedgroups[:j],sortedgroups[j:]
        alldistinct = True
        for group1 in lowsubset:
          for group2 in highsubset: alldistinct = alldistinct and reject[(group1,group2)]; # print group1,group2,reject[(group1,group2)]
        diff = sortedmeans[j][0]-sortedmeans[j-1][0] # could be negative if conds not sorted
        vals = lowsubset+['<<<--->>>']+highsubset+["%s" % alldistinct,"%0.1f" % diff]
        #print '\t'.join(vals)
        if diff<0: diff,lowsubset,highsubset = -diff,highsubset,lowsubset
        if alldistinct==True: candidate_divisions.append((diff,lowsubset,highsubset))
      if len(candidate_divisions)>0: 
        candidate_divisions.sort(reverse=True)
        (diff,lowsubset,highsubset) = candidate_divisions[0]
        Divisions[Rv] = (lowsubset,highsubset)

#The tukeyhsd of statsmodels doesn't return P value.
#So, if you want to know P value, calculate from these outputted value or use R.
#After saving the output into a variable res, you can get p-values by applying 
# psturng(np.abs(res.meandiffs / res.std_pairs), len(res.groupsunique), res.df_total), where psturng comes from from statsmodels.stats.libqsturng import psturng

  # save pvals as temp_anova.txt
  file = open(makefname("temp_anova.txt",rundir),"w")
  vals = "Rv Gene TAs".split()+Conditions.keys+"pval padj".split()
  file.write('\t'.join(vals)+EOL)  
  for (start,end,Rv,gene,strand) in genes:
   if Rv in Means:
    m,s = numpy.mean(Means[Rv]),numpy.std(Means[Rv])
    scv = s/m
    vals = [Rv,gene,str(len(TAsites[Rv]))]+["%0.1f" % x for x in Means[Rv]]+["%f" % x for x in [pvals[Rv],qvals[Rv]]] # could append SCV, but beware of how multi_TnSeq.py counts hits in temp_anova.txt
    file.write('\t'.join(vals)+EOL)
  file.close()

  # write sorted means and high/low subgroups into temp_subgroups.txt
  file = open(makefname("temp_divisions.txt",rundir),"w")
  vals = "Rv Gene TAs".split()+Conditions.keys
  file.write('\t'.join(vals)+EOL)  
  for (start,end,Rv,gene,strand) in genes:
   if Rv in qvals and qvals[Rv]<0.05:
    vals = [Rv,gene,str(len(TAsites[Rv]))]+["%0.1f" % x for x in Means[Rv]]
    vals += [x[1] for x in Sortedmeans[Rv]]+["%0.1f" % x[0] for x in Sortedmeans[Rv]]
    if Rv in Divisions:
      (lowsubset,highsubset) = Divisions[Rv]
      vals += [','.join(lowsubset),','.join(highsubset)]
    file.write('\t'.join(vals)+EOL)
  file.close()


  if False: # print gene means for each condition, and batch means
    for (start,end,Rv,gene,strand) in genes:
     vals = []
     if Rv in Means:
      batches = {}
      for r in range(Conditions.nrows):
        batch = Conditions.get(r,"Batch")
        if batch not in batches: batches[batch] = []
        batches[batch].append(Means[Rv][r])
        vals.append(Means[Rv][r]) ###
      for b in "CC1 CC2 CC3 KO".split(): vals.append(numpy.mean(batches[b]))
      print '\t'.join([Rv,gene]+["%0.1f" % x for x in vals]) ###
  ###sys.exit(0)



  # remove batch effects (subtract batch means from mean gene counts; adjust each to global mean)  
  if debatch and "Batch" in Conditions.headers:
    print "<BR>correcting means for batch effects..."
    for (start,end,Rv,gene,strand) in genes:
     if Rv in Means:
      batches = {}
      for r in range(Conditions.nrows):
        batch = Conditions.get(r,"Batch")
        if batch not in batches: batches[batch] = []
        batches[batch].append(Means[Rv][r])

      globalmean = numpy.mean(Means[Rv])
      keys = sorted(batches.keys())
      batchmeans = {}
      for b in keys: batchmeans[b] = numpy.mean(batches[b])

      for r in range(Conditions.nrows):
        batch = Conditions.get(r,"Batch")
        delta = globalmean-batchmeans[batch] # correction for each batch
        Means[Rv][r] = max(0,Means[Rv][r]+delta)
  
  if False: # print gene means for each condition (with batch corrections)
    vals = ['ORF','gene']+Conditions.getcol("Condition")
    print '\t'.join(vals)
    for (start,end,Rv,gene,strand) in genes:
       if Rv not in Means: continue
       vals = [Means[Rv][r] for r in range(Conditions.nrows)]
       print '\t'.join([Rv,gene]+["%0.1f" % x for x in vals]) 
    sys.exit(0)

  ################################
  # calc LFCs

  # consider normalizing by reference conditions?
  # consider only calculating for non-reference conditions?  
  
  lfcs = {}
  for (start,end,Rv,gene,strand) in genes:
   if Rv in qvals and qvals[Rv]<0.05:
    lfcvec = []
    for j in range(Conditions.nrows): # could remove means for reference conditions
      a = Means[Rv][j]
      b = numpy.median(Means[Rv])
      PC = 5
      lfc = math.log((a+5)/float(b+5),2)
      lfcvec.append(lfc)
    lfcs[Rv] = lfcvec

  if stdnorm:
    print "<BR>applying standard normalization to LFCs"
    for i,cond in enumerate(Conditions.keys):
      col = [row[i] for row in lfcs.values()]
      m,s = numpy.mean(col),numpy.std(col) # consider using quantiles(col,[0,0.25,0.5,0.75,1.0])
      for Rv in lfcs.keys():
        lfcs[Rv][i] = (lfcs[Rv][i]-m)/s

  file = open(makefname("temp_LFCs.txt",rundir),"w")
  vals = "Rv Gene".split()+Conditions.keys # or non-ref-conds
  file.write('\t'.join(vals)+EOL)
  cnt = 0
  for (start,end,Rv,gene,strand) in genes:
   if Rv in qvals and qvals[Rv]<0.05:
    vals = [Rv,gene]+["%0.3f" % x for x in lfcs[Rv]]
    file.write('\t'.join([str(x) for x in vals])+EOL)
    cnt += 1
  file.close()
  if cnt==0: print "error: no significantly varying genes found by ANOVA"; return

  os.system("Rscript %s/conditions_corrplot.R temp_LFCs.txt %s > /dev/null" % (DIR,rundir))

  os.system("Rscript %s/make_heatmap.R temp_LFCs.txt %s > /dev/null" % (DIR,rundir))

  os.system("Rscript %s/clustering.R temp_LFCs.txt %s %s %s > /dev/null" % (DIR,K,D,rundir))


#######################################

if __name__=="__main__":
  if len(sys.argv)<6:
    print "usage: python pipeline.py <ref.fna> <ref.prot_table> <combined_wig> <samples_metadata> <conditions_metadata> [-K n] [-D n] [-debatch] [-stdnorm]"
    sys.exit(-1)
  K,D = 10,5
  debatch,stdnorm = False,False
  if "-K" in sys.argv: K = int(sys.argv[sys.argv.index("-K")+1])
  if "-D" in sys.argv: D = int(sys.argv[sys.argv.index("-D")+1])
  if "-debatch" in sys.argv: debatch = True
  if "-stdnorm" in sys.argv: stdnorm = True
  print "K=%s, D=%s" % (K,D)
  pipeline(*sys.argv[1:6],K=K,D=D,debatch=debatch,stdnorm=stdnorm)

