import sys

class Spreadsheet:
  def __init__(self,filename): # ,key="Id"
    self.keys,self.data,self.headers = [],[],None
    self.rowhash,self.colhash = {},{}
    for line in open(filename,'rU'): # also handle mac files with x0d vs x0a
      if len(line)==0 or line[0]=='#': continue
      w = line.rstrip().split('\t') 
      if self.headers==None: 
        self.headers = w
        for i,h in enumerate(self.headers): 
          if h in self.colhash: print "error: '%s' appears multiple times in first row of '%s' - column headers must be unique" % (h,filename); sys.exit(-1)
          self.colhash[h] = i
        continue
      self.data.append(w)
      key = w[0]
      if key=="": continue
      if key in self.rowhash: print "error: '%s' appears multiple times in first column of '%s' - keys must be unique" % (key,filename); sys.exit(-1)
      self.rowhash[key] = len(self.keys)
      self.keys.append(key) # check if unique?
    self.ncols = max([len(row) for row in self.data])
    self.nrows = len(self.keys)
  def getrow(self,r):
    if r in self.rowhash: r = self.rowhash[r] 
    # check that it is an integer (in range)?
    return self.data[r]
  def getcol(self,c):
    if c in self.colhash: c = self.colhash[c] # otherwise, assume it is an integer
    return [r[c] for r in self.data] # what if not all same length?
  def get(self,r,c):
    if c in self.colhash: c = self.colhash[c] # otherwise, assume it is an integer
    if r in self.rowhash: r = self.rowhash[r]
    return self.data[r][c]

if __name__=="__main__":
  dat = Spreadsheet(sys.argv[1])

  print dat.nrows,dat.ncols
  #for x in dat.keys: print x
  #print dat.getrow("Tcell1")
  #print dat.get("Tcell2","Condition")
  #print dat.get("IFNg3","Filename")
  #print dat.getcol("Condition")

  for r in range(dat.nrows): # could transpose
    for c in range(dat.ncols):
      print dat.get(r,c),
    print
   
