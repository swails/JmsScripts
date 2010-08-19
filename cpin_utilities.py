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
#  printusage: prints the usage statement for the script
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

