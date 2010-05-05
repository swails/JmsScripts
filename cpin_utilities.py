# This python file contains functions useful to cpinutil.py

########################################################################
#  List of functions and a brief description of their purpose
#
#  fileexists: determines whether a passed file is found in the current 
#              directory
#
#  fileexists_noprint: determines whether a passed file is found in the
#                      current directory, but does not print an error 
#                      message
#
#  getallresinfo: Takes a given flag and returns every value in given 
#                 prmtop corresponding to that flag
#
#  printusage: prints the usage statement for the script
#
#  getbondiset: Gets the radius set used to build topology file
#
########################################################################

def fileexists(file):
   import sys

   try:
      f = open(file,'r')
   except IOError:
      print >> sys.stderr, 'Error: Specified file (' + file +          \
               ') does not exist!'
      return -1
   f.close()
   return 0

def fileexists_noprint(file):
   try:
      f = open(file,'r')
   except IOError:
      return -1
   f.close()
   return 0

def getallresinfo(prmtop, flag):
   import sys

   top = open(prmtop,'r')
   toplines = top.readlines()
   top.close()
   items = []

   for x in range(len(toplines)):
      if flag in toplines[x]:
         x = x + 1
         while not 'FORMAT' in toplines[x]:
            x = x + 1
         x = x + 1
         while not 'FLAG' in toplines[x]:
            words = toplines[x].split()
            for y in range(len(words)):
               items.append(words[y])
            x = x + 1
         break

   return items

def printusage():
   import sys

   print "cpinutil.py -p <prmtop>                  \\"
   print "            -igb <2 or 5>                \\"
   print "            -resname <resname list>      \\"
   print "            -notresname <resname list>   \\"
   print "            -resnum <residue numbers>    \\"
   print "            -notresnum <residue numbers> \\"
   print "            -minpKa <value>              \\"
   print "            -maxpKa <value>              \\"
   print "            -states <list of states>     \\"
   print "            -system <system name>        \\"
   print "           [--ignore-warnings]" 

   sys.exit()

def getbondiset(topname):

   topfile = open(topname, 'r')
   toplines = topfile.readlines()
   topfile.close()

   for x in range(len(toplines)):

      if '%FLAG RADIUS_SET' in toplines[x]:
         x = x + 1
         while not toplines[x].startswith('%FORMAT'):
            x = x + 1
         return toplines[x+1].strip()

   return 'NO RADII SET SPECIFIED'

