#!/usr/bin/env python

from optparse import OptionParser
import re
import sys

# This allows for a name mapping
NAME_DICT = {
             'AS4' : 'ASP',
             'GL4' : 'GLU',
             'HIP' : 'HIS'
            }

def resname(r):
   """ Return a string with the 'normal' residue name """
   global NAME_DICT
   try:
      return NAME_DICT[r]
   except KeyError:
      return str(r)

class TitratableResidue(object):
   """ A titratable residue with data extracted from a cpout file """
   def __init__(self, name, num, offset, pred, frac_prot, transitions):
      self.name = name
      self.num = int(num)
      self.offset = float(offset)
      self.pred = float(pred)
      self.frac_prot = float(frac_prot)
      self.transitions = int(transitions)

   def __str__(self):
      return ("%s" * 6) % (self.name, self.num, self.offset, self.pred,
                           self.frac_prot, self.transitions)

class pHRecord(object):
   """ A pH Record from a given output file """
   startre = re.compile(r'Solvent pH is *([-+]?\d+\.\d+)')
   resre = re.compile(r'([A-Z4]+) *(\d+) *: Offset *([-+]?\d+(?:\.\d*)?|\.\d+|-*[Ii]nf) *Pred *([-+]?\d+(?:\.\d*)?|\.\d+|-*[Ii]nf) *Frac Prot *(\d+(?:\.\d*)?) *Transitions *(\d+)')

   def __init__(self, lines):
      """ Open the file and parse out the information """
      self.residues = []
      for line in lines:
         rematch = self.startre.match(line)
         if rematch:
            self.pH = float(rematch.groups()[0])
            continue
         rematch = self.resre.match(line)
         if rematch:
            self.residues.append(TitratableResidue(*rematch.groups()))
            continue

   def write_header(self, dest):
      """ Writes a header for a gnuplot data file """
      dest.write('#%7s' % 'pH')
      for r in self.residues:
         dest.write(' %8s' % (resname(r.name) + str(r.num)))
      dest.write('\n')

   def write_line(self, dest, deprot=True):
      dest.write('%8.3f' % self.pH)
      if deprot:
         for r in self.residues:
            dest.write(' %8.3f' % (1-r.frac_prot))
      else:
         for r in self.residues:
            dest.write(' %8.3f' % r.frac_prot)
      dest.write('\n')

def main():
   """ Parses the main file and writes the results """
   parser = OptionParser()
   parser.add_option('-i', '--input-file', dest='infile', metavar='FILE',
                     default=None, help='Input file with titration data')
   parser.add_option('-o', '--output-file', dest='outfile', metavar='FILE',
                     default=None, help='Output file suitable for graphing')
   parser.add_option('-p', '--frac-prot', dest='frac_prot', action='store_true',
                     default=False, help='Print out pH vs. fraction protonated.'
                     ' NOT default behavior')
   parser.add_option('-d', '--frac-deprot', dest='frac_prot',
                     action='store_false', help='Print out pH vs. fraction '
                     'protonated. Default behavior')
   
   opt, arg = parser.parse_args()

   if opt.infile is None or opt.outfile is None:
      print 'Error: Must specify input/output file!'
      sys.exit(1)

   if opt.frac_prot:
      print 'Printing fraction protonated'
   else:
      print 'Printing fraction deprotonated'

   lines = []
   ph_groups = []
   for line in open(opt.infile, 'r'):
      if line.startswith('Average'):
         continue
      elif not line.strip():
         if lines:
            ph_groups.append(pHRecord(lines))
            lines = []
         continue
      lines.append(line)

   # Now print out the file
   f = open(opt.outfile, 'w')
   ph_groups[0].write_header(f)
   for g in ph_groups:
      g.write_line(f, deprot=not opt.frac_prot)

   print 'Done!'

if __name__ == '__main__':
   main()
