#!/usr/bin/env python

from __future__ import division
from optparse import OptionParser
import sys
from math import log

parser = OptionParser(epilog='This adjusts the reference pKa based on input')
parser.add_option('--ore', dest='ore', type='float', default=None,
                  help='Original Reference Energy', metavar='FLOAT')
parser.add_option('--target-pka', dest='tpka', type='float', default=None,
                  help='Target pKa (what you want ideally)', metavar='FLOAT')
parser.add_option('--observed-pka', dest='titpka', type='float', default=None,
                  help='Observed pKa (what you actually got)', metavar='FLOAT')
parser.add_option('--deprotonating', dest='deprot', action='store_true',
                  default=False, help='Is this residue deprotonating from the '
                  'first to the second state? Default is No.')

opt, arg = parser.parse_args()

if arg:
   parser.print_help()
   print 'Too many arguments!'
   sys.exit(1)

if opt.ore is None or opt.tpka is None or opt.titpka is None:
   print 'I need --ore, --target-pka, and --observed-pka!'
   parser.print_help()
   sys.exit(1)

KB = 0.00199
TEMP = 300
DEPROTONATING = opt.deprot

def nre(ore, tpka, titpka):
   global KB, TEMP, DEPROTONATING
   if DEPROTONATING:
      return ore + KB * log(10) * TEMP * (titpka - tpka)
   else:
      return ore - KB * log(10) * TEMP * (titpka - tpka)

def newti(nre, tpka):
   global DEPROTONATING
   if DEPROTONATING:
      return nre - KB * log(10) * TEMP * tpka
   else:
      return nre + KB * log(10) * TEMP * tpka


mynre = nre(opt.ore, opt.tpka, opt.titpka)
print 'The new pKa-adjusted reference energy is  %.8f' % mynre
print 'The new TI-equivalent reference energy is %.8f' % newti(mynre, opt.tpka)
