import sys

# note this parse requires full sentences to be on 1 line
# case-sensitive
# '#' is comment char (anywhere on line)
# no quotes
# all chars besides ( and ) are treated as parts of tokens
# top level is assumed to be only 1 Sexpr (atom or list), not multiple exprss automatically enclosed in parens

def tokenize(s):
  toks = []
  i,n = 0,len(s)
  while i<n:
    if s[i]==' ': i += 1; continue # skip whitespace
    elif s[i]=='(' or s[i]==')': toks.append(s[i]); i += 1; continue
    else:
      j = i
      while j<n and s[j] not in " ()": j += 1 # scan for end of token
      toks.append(s[i:j])
      i = j
  return toks

class Sexpr:
  def __init__(self,toks,i):
    self.atom = None
    self.list = []
    self.next = i
    # parse expr
    n = len(toks)
    if i>=n: return # empty list
    if toks[i]=='(':
      i += 1
      while i<n and toks[i]!=')':
        expr = Sexpr(toks,i)
        self.list.append(expr)
        i = expr.next
      if i==n: raise Exception("syntax error: list not closed: %s" % ' '.join(toks))
      self.next = i+1; return
    elif toks[i]==')': self.next = i; return # this is not my close paren; parent might be interested in it
    else: self.atom = toks[i]; self.next = i+1; return
  
  def toString(self):
    if self.atom!=None: return self.atom
    else: 
      x = ' '.join([x.toString() for x in self.list])
      return "(%s)" % x

  # assume expr is a clause

  def toDIMACS(self):
    if self.atom!=None: return self.atom
    lits = ""
    for arg in self.list[1:]:
      if arg.atom!=None: lits += " "+arg.atom
      else: lits += " -"+arg.list[1].atom # assume it is a neg lit
    return lits.strip()

##########################
# convert to CNF
# 1. eliminate <-> and xor
# 2. eliminate implications
# 3. push negations inward (using DeMorgan's Law if necessary)
# 4. distribute and over or
# returns an atom, a clause (1-level with 'or' oper), or a 2-level CNF expression (with 'and' over 'or') which will get broken into multiple clauses by caller
# also collapses nested or's: (or (or P Q) (or R S)) -> (or P Q R S)

def convCNF(expr):
  expr = eliminate_biconditionals(expr)
  expr = eliminate_xors(expr)
  expr = eliminate_implications(expr)
  expr = push_negs_inward(expr)
  expr = distribute(expr)
  expr = or_wrapper(expr) 
  return expr

# note: don't do this if top-level oper is 'and' (it will be handling by caller by breaking into clauses)

def or_wrapper(expr):
  if expr.list==None or (is_disjunction(expr)==False and is_conjunction(expr)==False):
    return make_Sexpr_from_list(['or',expr])
  return expr

# make new expr from list of strings or Sexpr objects (because python does not have polymorphic constructors)

def make_Sexpr_from_list(lst): 
  expr = Sexpr(tokenize(""),0) # make an empty expr
  expr.list = [x if isinstance(x,Sexpr) else Sexpr(tokenize(x),0) for x in lst]
  return expr

def negate(expr): return make_Sexpr_from_list(['not',expr])

# assume there is always at least an oper, and implications always have 2 args (could do error checking)

def eliminate_implications(expr):
  if expr.atom!=None: return expr
  oper = expr.list[0]
  if oper.atom=="->" or oper.atom=="implies":
    if len(expr.list)!=3: raise Exception("implication can only have 2 arguments")
    a,b = expr.list[1],expr.list[2]
    newlist = ['or',negate(a),b]
    newexpr = make_Sexpr_from_list(newlist)
    return eliminate_implications(newexpr)
  else:
    return make_Sexpr_from_list([eliminate_implications(x) for x in expr.list])

def eliminate_biconditionals(expr):
  if expr.atom!=None: return expr
  oper = expr.list[0]
  if oper.atom=="<->":
    if len(expr.list)!=3: raise Exception("biconditional can only have 2 arguments")
    a,b = expr.list[1],expr.list[2]
    newlist1 = ['->',a,b]
    newlist2 = ['->',b,a]
    newexpr1 = make_Sexpr_from_list(newlist1) 
    newexpr2 = make_Sexpr_from_list(newlist2) 
    newlist = ["and",newexpr1,newexpr2]
    newexpr = make_Sexpr_from_list(newlist) 
    return eliminate_biconditionals(newexpr)
  else:
    return make_Sexpr_from_list([eliminate_biconditionals(x) for x in expr.list])

def eliminate_xors(expr):
  if expr.atom!=None: return expr
  oper = expr.list[0]
  if oper.atom=="xor":
    if len(expr.list)!=3: raise Exception("xor can only have 2 arguments")
    a,b = expr.list[1],expr.list[2]
    newlist1 = ['or',a,negate(b)]
    newlist2 = ['or',b,negate(a)]
    newexpr1 = make_Sexpr_from_list(newlist1) 
    newexpr2 = make_Sexpr_from_list(newlist2) 
    newlist = ["and",newexpr1,newexpr2]
    newexpr = make_Sexpr_from_list(newlist) 
    return eliminate_xors(newexpr)
  else:
    return make_Sexpr_from_list([eliminate_xors(x) for x in expr.list])

# assume all lists have >=1 element (oper), and negations have 2 elements (1 subexpr)

def push_negs_inward(expr):
  if expr.atom!=None: return expr
  if expr.list[0].atom=='not': # assume each list has at least 1 item, the oper, which is an atom
    subexpr = expr.list[1]
    if subexpr.atom!=None: return expr # (not X) is fine
    oper = subexpr.list[0]
    if oper.atom=="and" or oper.atom=="or":
      args = subexpr.list[1:]
      args = [negate(x) for x in args]
      args = [push_negs_inward(x) for x in args]
      flipop = "and" if oper.atom=="or" else "or"
      return make_Sexpr_from_list([flipop]+args)
    elif oper.atom=="not": return push_negs_inward(subexpr.list[1]) # double-neg elimination
    else: raise Exception("error: don't know how to push negation inward over the following operator: %s" % oper.atom)
  else: return make_Sexpr_from_list([push_negs_inward(x) for x in expr.list]) 

def is_literal(expr): return expr.atom!=None or expr.list[0].atom=="not"

def is_conjunction(expr): return expr.atom==None and expr.list[0].atom=="and"

def is_disjunction(expr): return expr.atom==None and expr.list[0].atom=="or"

# assume negs have already been pushed inward
# invariant: always returns at most 2-level expr: 
#   atom (prop), not (neg lit), or (clause, 1-level), and (cnd, 2-level)
#

def distribute(expr):
  #print("distributing: "+expr.toString())
  if is_literal(expr): return expr
  oper = expr.list[0].atom
  #if oper not in "and or": raise Exception("error: unexpected operator %s in distribute()" % oper) # what about caps, AND and OR?
  if not (is_conjunction(expr) or is_disjunction(expr)):
    raise Exception("error: unexpected operator %s in distribute()" % oper) # what about caps, AND and OR?
  args = [distribute(x) for x in expr.list[1:]] # bottom-up

  if oper=="and":
    absorbed = []
    for arg in args:
      #if is_literal(arg) or arg.list[0].atom=="or": absorbed.append(arg)
      if is_conjunction(arg): absorbed += arg.list[1:]
      else: absorbed.append(arg)
    return make_Sexpr_from_list(["and"]+absorbed)

  elif oper=="or":
    # if one of the args is an 'and', distribute its terms over all the other args, and recurse; if another 'or', absorb its args
    # special case: if there is only one arg, drop the 'or'
    if len(args)==1: return args[0]
    absorbed = []
    for i in range(len(args)):
      arg = args[i]
      if is_conjunction(arg):
        terms = arg.list[1:] # in conjunct
        others = args[:i]+args[i+1:] # disjuncts
        newargs = []
        for other in others: 
          for term in terms:
            disjoined = make_Sexpr_from_list(["or",term]+others)
            newargs.append(disjoined)
        return distribute(make_Sexpr_from_list(["and"]+newargs)) # after processing this arg, use recursion to handle any others
      if is_disjunction(arg): absorbed += arg.list[1:]
      else: absorbed.append(arg)
    return make_Sexpr_from_list(["or"]+absorbed) # if no args are conjunctions

def validate_PropLog(expr):
  if expr.atom!=None: return True
  if len(expr.list)==0: return False
  oper = expr.list[0]
  if oper.atom==None: return False
  oper = oper.atom
  if oper not in ["<->","->","implies","and","or","not","xor"]: return False 
  args = expr.list[1:]
  for arg in args:
    if validate_PropLog(arg)==False: return False
  if (oper=="<->" or oper=="->" or oper=="xor") and len(args)!=2: return False
  if oper=="not" and len(args)!=1: return False
  # 'and' and 'or' can have any number of args (0 or more)
  return True

##########################

DIMACS = False

if __name__=="__main__":
 if len(sys.argv)<2: sys.stderr.write("usage: python convCNF.py <propositional_kb> [-DIMACS]\n"); sys.exit(0)
 if "-DIMACS" in sys.argv: DIMACS = True

 for line in open(sys.argv[1]):
  line = line.rstrip()
  print("# "+line)
  if line.startswith('#') or len(line)==0: continue
  if '#' in line: line = line[:line.index('#')]
  toks = tokenize(line)
  expr = Sexpr(toks,0) # I could check this: validate_PropLog(expr)

  cnf = convCNF(expr)

  if is_conjunction(cnf): clauses = cnf.list[1:]
  else: clauses = [cnf]
  clauses = [or_wrapper(c) for c in clauses] # make sure clauses have an 'or' at the top level; (should not have 'and's at top level)

  print("")
  if DIMACS:
    for clause in clauses: print(clause.toDIMACS())
  else:
    for clause in clauses: print(clause.toString())
  print("")

