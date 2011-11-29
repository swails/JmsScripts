#!/usr/bin/env python

import csv, sys, re, math
from optparse import OptionParser

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Residue(object):
   """ A holder for all of the data for a residue """
   def __init__(self, name='', number=0, offset=0.0, frac_prot=0.0, 
                pred=0.0, transitions=0):
      self.name = name
      self.number = int(number)
      self.offset = float(offset)
      self.frac_prot = float(frac_prot)
      self.pred = float(pred)
      self.transitions = int(transitions)

   def __eq__(self, other):
      return self.name == other.name and self.number == other.number
   def __gt__(self, other):
      return self.number > other.number
   def __lt__(self, other):
      return self.number < other.number
   def __ne__(self, other):
      return not self == other
   def __ge__(self, other):
      return self.number >= other.number
   def __le__(self, other):
      return self.number <= other.number

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class PhSet(object):
   """ A set of pH values """
   def __init__(self, pH=7.0):
      self.pH = float(pH)
      self.residue_list = []

   def add_residue(self, **kwargs):
      """ Adds a residue to the list of residues """
      self.residue_list.append(Residue(**kwargs))

   def __gt__(self, other):
      return self.pH > other.pH
   def __eq__(self, other):
      return self.pH == other.pH
   def __ge__(self, other):
      return self > other or self == other
   def __lt__(self, other):
      return not self >= other
   def __le__(self, other):
      return not self > other
   def __ne__(self, other):
      return self.pH != other.pH
   def __contains__(self, element):
      for res in self.residue_list:
         if res == element: return True
      return False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def set_list(lines):
   """ Converts lines of input into a list of PhSet's """
   start_of_set = \
          re.compile(r'Solvent pH is *([-+]?\d+(?:\.\d*)?|\.\d+)')
   residue_line = \
          re.compile(r'([A-Z4]+) *(\d+) *: Offset *([-+]?\d+(?:\.\d*)?|\.\d+|-*[Ii]nf) *Pred *([-+]?\d+(?:\.\d*)?|\.\d+|-*[Ii]nf) *Frac Prot *(\d+(?:\.\d*)?) *Transitions *(\d+)')
   end_of_set = \
          re.compile(r'Average total molecular protonation: *(\d+(?:\.\d*))')

   ph_setlist = []
   i = 0
   while i < len(lines):
      start = start_of_set.match(lines[i])
      if start:
         my_set = PhSet(start.groups()[0])
         end = end_of_set.match(lines[i])
         while i < len(lines) and not end:
            my_info = residue_line.match(lines[i])
            if my_info:
               info = my_info.groups()
               my_set.add_residue(name=info[0], number=info[1], offset=info[2],
                                  pred=info[3], frac_prot=info[4],
                                  transitions=info[5])
            else:
               end = end_of_set.match(lines[i])
               if end:
                  my_set.avg_prot = float(end.groups()[0])
                  break
            i += 1
         ph_setlist.append(my_set)
      i += 1

   return ph_setlist
            
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def _get_residues(ph_setlist):
   # Get all of the residues in the ph_setlist
   residues = []
   for phset in ph_setlist:
      if not residues: residues = phset.residue_list[:]
      for residue in phset.residue_list:
         if residue in residues: continue
         residues.append(residue)
   # Sort the residue list
   residues.sort()

   return residues
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def _individual_info(lines, outfile, data_type, data_desc, hill=False):
   # Write the fraction protonated

   ph_setlist = set_list(lines)
   ph_setlist.sort()

   residues = _get_residues(ph_setlist)

   # Write out the first couple rows (just headers)
   firstline = []
   secondline = []
   for res in residues:
      firstline += ['%s %d' % (res.name, res.number), '', '']
      secondline += ['pH', data_desc, '']
   outfile.writerow(firstline)
   outfile.writerow(secondline)

   # Now write out each row
   for phset in ph_setlist:
      line = []
      for res in residues:
         if not res in phset.residue_list:
            line += ['-----', '', '']
            continue
         myres = phset.residue_list[phset.residue_list.index(res)]
         if not hill:
            line += [phset.pH, getattr(myres, data_type), '']
         else:
            if myres.frac_prot == 1:
               line += [phset.pH, '-inf', '']
            elif myres.frac_prot == 0:
               line += [phset.pH, 'inf', '']
            else:
               line += [phset.pH, 
                        math.log10((1-myres.frac_prot)/myres.frac_prot), '']
      outfile.writerow(line)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def pka_info(lines, outfile):
   _individual_info(lines, outfile, 'pred', 'predicted pKa')
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def frac_info(lines, outfile):
   _individual_info(lines, outfile, 'frac_prot', 'Fraction Protonated')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def hill_info(lines, outfile):
   _individual_info(lines, outfile, 'frac_prot', 'Hill plot', hill=True)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def trans_info(lines, outfile):
   _individual_info(lines, outfile, 'transitions', 'Transitions')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def offset_info(lines, outfile):
   _individual_info(lines, outfile, 'offset', 'Offset From pH')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def full_info(lines, outfile):
   
   ph_setlist = set_list(lines)
   ph_setlist.sort()

   # Write out the first couple rows (just headers)
   firstline = ['Residue', '']
   secondline = ['', '']
   for phset in ph_setlist:
      firstline += ['', 'pH %.2f' % phset.pH, '', '']
      secondline += ['pKa', 'Frac Prot', 'Transitions', '']
   outfile.writerow(firstline)
   outfile.writerow(secondline)

   residues = _get_residues(ph_setlist)

   # Now write out all of the rows
   res_rows = [['%4s %4d' % (res.name, res.number), ''] for res in residues]
   for i, res in enumerate(residues):
      for phset in ph_setlist:
         if not res in phset:
            res_rows[i] += ['-----', '-----', '-----', '']
            continue
         res_item = phset.residue_list[phset.residue_list.index(res)]
         res_rows[i] += [res_item.pred, res_item.frac_prot, res_item.transitions, '']

   for row in res_rows: outfile.writerow(row)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def main():
   parser = OptionParser(usage='%prog [Options] stat_file1 [stat_file2 ...]')
   parser.add_option('-o', '--output', dest='output', default='results.csv',
                help='Final CSV file with total statistics. [Default %default]')
   parser.add_option('-i', '--info', dest='info', default=None,
                   help='What info do you want printed out?. Available ' +
                   'options are [pKa, frac_prot, transitions, offset, hill]. ' +
                   'Default: full description')
   opt, args = parser.parse_args()

   if not opt.output.endswith('.csv'): opt.output += '.csv'

   if opt.info == None:
      run_method = full_info
   elif opt.info.lower() == 'pka':
      run_method = pka_info
   elif opt.info.lower().startswith('frac'):
      run_method = frac_info
   elif opt.info.lower().startswith('trans'):
      run_method = trans_info
   elif opt.info.lower() == 'offset':
      run_method = offset_info
   elif opt.info.lower() == 'hill':
      run_method = hill_info
   else:
      run_method = full_info

   outfile = csv.writer(open(opt.output, 'wb'))
   lines = []
   if not args:
      lines = sys.stdin.readlines()
   else:
      for arg in args:
         tmp = open(arg, 'r')
         lines += tmp.readlines()
         tmp.close()

   for i in range(len(lines)):
      lines[i] = lines[i].strip()

   run_method(lines, outfile)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':
   main()
   print 'Done!'
