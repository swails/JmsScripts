#!/usr/bin/env python

import sys

def usage():
   print >> sys.stderr, 'Usage: splitNMR_PDB.py <PDB>'
   sys.exit()


fileno = 1

if len(sys.argv) != 2:
   usage()

file = open(sys.argv[1],'r')

wfile = open('%s.%s' % (sys.argv[1], fileno), 'w')

for line in file:
   if line.strip() == "ENDMDL":
      wfile.close()
      fileno += 1
      wfile = open('%s.%s' % (sys.argv[1], fileno), 'w')
      continue
   elif not "MODEL" in line:
      wfile.write(line)

file.close()
wfile.close()

print >> sys.stdout, 'Processed %s models from %s' % (fileno, sys.argv[1])
