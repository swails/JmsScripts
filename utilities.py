#####################################################################################
#                                                                                   #
#  This is a module that contains functions generally useful for doing all sorts of #
#  things. Whenever I write a quick utility/script that manipulates stuff, I will   #
#  add generally useful functions here. This is expected to grow, and all functions #
#  return -1 in error (a protocol I am adopting). If -1 is a valid result, then     #
#  'none' will be returned in string form.                                          #
#                                                                                   #
#          Last updated: 04/23/2010                                                 #
#                                                                                   #
#####################################################################################



# Definition of included utilities. Add descriptions to whatever function(s) you add!

################################################################################
#
# which: searches the path for a given executable
# linecount: count number of lines in a file
# average: returns the average of a given list
# stdev: returns the standard deviation of a given list and average (a non-float
#        value for average will make the function call 'average' to compute it
# digit: returns the given digit (power of 10, so it starts 0 for 1's place) of
#        a given number
# resnum: returns the number of residues in a given amber topology file
# getallresinfo: prints all values for a given prmtop FLAG in a given topology
# fileexists: alias for os.path.exists (didn't know about that at first)
# getresindex: returns the rank in alphabetical order of a given amino acid
# getresdecmp: returns the amino acid decomposition of a given prmtop in an
#              alphabeitcal ordered array
# minmax: returns the maximum and minimum of a data set
# add_log: adds logs of two numbers and stores it as the log of the sum
#
################################################################################

from chemistry.amber.readparm import AmberParm
import os, math

def which(program):
   def is_exe(fpath):
      return os.path.exists(fpath) and os.access(fpath, os.X_OK)

   fpath, fname = os.path.split(program)
   if fpath:
      if is_exe(program):
         return program
   else:
      for path in os.environ["PATH"].split(os.pathsep):
         exe_file = os.path.join(path, program)
         if is_exe(exe_file):
            return exe_file
   return None

def linecount(fname):
   if type(fname).__name__ != 'str':
      raise TypeError('linecount: Expected a string!')
   fl = open(fname, 'r')
   for i, l in enumerate(fl): pass
   return i+1

def average(list):
   from math import fsum

   return fsum(list) / float(len(list))

def stdev(list,avg):
# it must be passed 2 arguments. pass it a string if you want stdev to calculate
# the average first.
   import math

   if len(list) < 2:
      return 0
   try:
      ave = float(avg)
   except ValueError:
      ave = average(list)

   sum = 0

   for x in range(len(list)):

      sum = sum + (list[x] - ave) ** 2

   return math.sqrt( sum / len(list) )

def digit(num, powten):

   import math
   tmp = math.floor(num / 10 ** powten)
   return int(tmp % 10)

def round(number, decimals):
# YOU DO NOT NEED TO EVER USE THIS FUNCTION. Just use the format utility
# built into python itself: '{0:.decimal}'.format(floating_number)

   import math
   try:
      num = float(number)
   except ValueError:
      return -1

   try:
      dec = int(decimals)
   except ValueError:
      return -1

   test = (num * 10 ** decimals) - math.floor(num * 10 ** decimals)
   test = math.floor(test * 10)
   if test < 5:
      toreturn = math.floor(num * 10 ** decimals) / 10 ** decimals
      return toreturn
   else:
      toreturn = math.ceil(num * 10 ** decimals) / 10 ** decimals
      return toreturn


def resnum(topfile):

   parm = AmberParm(topfile)
   return parm.ptr("NRES")

def natom(topfile):
   
   parm = AmberParm(topfile)
   return parm.ptr("NATOM")

def getresinfo(res, topname, flag):

   parm = AmberParm(topname)
   return parm.parm_data[flag][res-1] # and simply return the residue of interest

def getallresinfo(topname, flag):
   parm = AmberParm(topname)
   return parm.parm_data[flag]

def fileexists(file):
   if os.path.exists(file): return 0
   else: return -1

def getresindex(resname):
# return the amino acid's rank in terms of alpabetical order
   if resname == 'ALA':
      return 1
   if resname == 'ARG':
      return 2
   if resname == 'ASN':
      return 3
   if resname in 'ASP ASH AS4':
      return 4
   if resname in 'CYS CYX CYM':
      return 5
   if resname in 'GLU GLH GL4':
      return 6
   if resname == 'GLN':
      return 7
   if resname == 'GLY':
      return 8
   if resname in 'HIP HID HIE HIS':
      return 9
   if resname == 'ILE':
      return 10
   if resname == 'LEU':
      return 11
   if resname in 'LYN LYS':
      return 12
   if resname == 'MET':
      return 13
   if resname == 'PHE':
      return 14
   if resname == 'PRO':
      return 15
   if resname == 'SER':
      return 16
   if resname == 'THR':
      return 17
   if resname == 'TRP':
      return 18
   if resname == 'TYR':
      return 19
   if resname == 'VAL':
      return 20
   if resname == 'WAT':
      return 21

   return 22

def getresdecmp_prmtop(file):

   import sys

   residue_decomp=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
# residue indices are : ALA, ARG, ASN, ASP, CYS, GLU, GLN, GLY
#                       HIS, ILE, LEU, LYS, MET, PHE, PRO, SER
#                       THR, TRP, VAL, WAT, UNK
   unkres = []

   try:
      prmtop = open(file,'r')
   except IOError:
      print >> sys.stderr, "Error: Topology file " + file +            \
               " does not exist!"
      sys.exit()
   
   lines = prmtop.readlines()
   prmtop.close()
   residues = getallresinfo('RESIDUE_LABEL')
   
   for x in range(len(residues)):
      index = getresindex(residues[x]) - 1
      residue_decomp[index] = residue_decomp[index] + 1
      if index == 21:
         add = True
         for y in range(len(unkres)):
            if unkres[y] == residues[x]:
               add = False
         if add:
            unkres.append(residues[x])
   
   print unkres 

   return residue_decomp

def getresdecmp_pdb(file):
   import sys
   
   residue_decomp=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
# residue indices are : ALA, ARG, ASN, ASP, CYS, GLU, GLN, GLY
#                       HIS, ILE, LEU, LYS, MET, PHE, PRO, SER
#                       THR, TRP, VAL, WAT, UNK, POS, NEG

   try:
      pdb = open(file,'r')
   except IOError:
      print >> sys.stderr, "Error: PDB file " + file + "doesn't exist!"

   lines = pdb.readlines()
   pdb.close()
   residues = []
   curres = 0
   resindex = 4
   foundindex = False

   for x in range(len(lines)):
      if not lines[x].startswith("ATOM") and not lines[x].startswith(   \
         "HETAT"):
         continue
      else:
         while not foundindex:
            try:
               test = int(lines[x].split()[resindex])
            except ValueError:
               resindex += 1
         break

   for x in range(len(lines)):
      if not lines[x].startswith("ATOM") and not lines[x].startswith(   \
         "HETAT"):
         continue
      words = lines[x].split()
      if curres == int(words[resindex]):
         continue
      curres = int(words[4])
      residues.append(words[3])
      
   for x in range(len(residues)):
      index = getresindex(residues[x]) - 1
      residue_decomp[index] = residue_decomp[index] + 1
      if residues[x] in "ASP GLU":
         residue_decomp[21] += 1
      elif residues[x] in "LYS ARG HIP":
         residue_decomp[20] += 1

   return residue_decomp

def minmax(list):
   return [min(list),max(list)]

def add_log(a, b):
   """ Adds log(a) + log(b) and returns log(a+b) """
   # Make sure a is the largest and b is the smallest
   a, b = max(a,b), min(a,b)
   return a + math.log(1 + math.exp(b-a))
