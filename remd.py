"""
This file provides classes to interact with replica exchange files and provide
statistics for the simulation
"""

import re
import numpy as np

class Replica(object):
   """ Just a container """

class RemdError(Exception):
   """ Error when reading remd files """

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class RemLog(object):
   """ Replica exchange log file """

   numexch_desc_re = re.compile(r'# numexchg is *(\d+)')

   #================================================

   def __init__(self, fname):
      """ Loads a replica exchange log file """
      self.fname = fname
      self.reps = [] # list of all replicas
      self.numexchg = 0
      self.file = open(fname, 'r')
      self.numexchg = self._get_numexchg()
      self.values = []
      self._get_replicas()
      self._parse()

   #================================================

   def _get_numexchg(self):
      rawline = self.file.readline()
      while rawline:
         rematch = self.numexch_desc_re.match(rawline)
         if rematch:
            return int(rematch.groups()[0])
         rawline = self.file.readline()
      raise RemdError('Could not find numexchg!')

   #================================================

   def _get_replicas(self):
      " Gets replica information from the first block -- must be inherited! "
      raise RemdError('_get_replicas: Virtual method only!')

   #================================================

   def _parse(self):
      """ Parses the file -- must be inherited! """
      raise RemdError('_parse: Virtual method only!')

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class TempRemLog(RemLog):
   """ Replica exchange log file """

   line_re = re.compile(r' *(\d+) *([-+]?\d\.\d+) *(\d+\.\d+) *([-+]?\d+\.\d+) *(\d+\.\d+) *(\d+\.\d+) *(\d\.\d+) *([+-]?\d+)')

   #================================================

   def _get_replicas(self):
      " Gets all of the replica information from the first block of repinfo "
      rawline = self.file.readline()
      while rawline:
         rematch = self.line_re.match(rawline)
         if not rematch:
            rawline = self.file.readline()
            continue
         # Now we have our first block of replicas!
         while rematch:
            rep = Replica()
            self.reps.append(rep)
            # Each replica will have an index array to trace its trajectory
            # through replica space. Indexing starts from 0
            rep.index = [0 for i in range(self.numexchg)]
            rep.potene = np.zeros(self.numexchg)
            rep.success_rate = np.zeros(self.numexchg)
            rep.old_temp = np.zeros(self.numexchg)
            rep.new_temp = np.zeros(self.numexchg)
            rep.index[0] = int(rematch.groups()[0]) - 1
            rep.potene[0] = float(rematch.groups()[3])
            rep.success_rate[0] = float(rematch.groups()[6])
            rep.old_temp[0] = float(rematch.groups()[4])
            rep.new_temp[0] = float(rematch.groups()[5])
            self.values.append(rep.old_temp[0])
            rematch = self.line_re.match(self.file.readline())
         # Finished our first block of replicas. Bail out now
         break
      # Now sort our temperatures and adjust our index
      self.values.sort()
      for rep in self.reps:
         rep.index[0] = self.values.index(rep.old_temp[0])

   #================================================

   def _parse(self):
      """ Parses the rem.log file and loads the data arrays """
      rawline = self.file.readline()
      num_record = 1
      while rawline:
         rematch = self.line_re.match(rawline)
         if not rematch:
            rawline = self.file.readline()
            continue
         while rematch:
            repnum,j1,j2,potene,old_temp,new_temp,success_rate,j4 = \
                                 rematch.groups()
            repnum = int(repnum) - 1
            potene = float(potene)
            old_temp = float(old_temp)
            new_temp = float(new_temp)
            success_rate = float(success_rate)
            self.reps[repnum].potene[num_record] = potene
            self.reps[repnum].old_temp[num_record] = old_temp
            self.reps[repnum].new_temp[num_record] = new_temp
            self.reps[repnum].success_rate[num_record] = success_rate
            self.reps[repnum].index[num_record] = self.values.index(old_temp)
            rematch = self.line_re.match(self.file.readline())
         num_record += 1
         rawline = self.file.readline()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class HRemLog(RemLog):
   """ A class for H-REMD log file """

   line_re = re.compile(r' *(\d+) *(\d+) *([-+]?\d+\.\d+) *([-+]?\d+\.\d+) *([-+]?\d+\.\d+) *([-+]?\d+\.\d+) *([-+]?\d+\.\d+) *([TF]) *([-+]?\d+\.\d+)')

   def _get_replicas(self):
      " Gets all of the replica information from the first block of repinfo "
      rawline = self.file.readline()
      while rawline:
         rematch = self.line_re.match(rawline)
         if not rematch:
            rawline = self.file.readline()
            continue
         # Now we have our first block of replicas
         while rematch:
            rep = Replica()
            self.reps.append(rep)
            # Each replica will have an index array to trace its trajectory
            # through replica space. Indexing starts from 0
            rep.index = [0 for i in range(self.numexchg)]
            rep.neighbor_index = [0 for i in range(self.numexchg)]
            rep.temp = np.zeros(self.numexchg)
            rep.potene1 = np.zeros(self.numexchg)
            rep.potene2 = np.zeros(self.numexchg)
            rep.left_fe = np.zeros(self.numexchg)
            rep.right_fe = np.zeros(self.numexchg)
            rep.success = ['F' for i in range(self.numexchg)]
            rep.success_rate = np.zeros(self.numexchg)
            rep.index[0] = int(rematch.groups()[0]) - 1
            rep.neighbor_index[0] = int(rematch.groups()[1]) - 1
            rep.temp[0] = float(rematch.groups()[2])
            rep.potene1[0] = float(rematch.groups()[3])
            rep.potene2[0] = float(rematch.groups()[4])
            rep.left_fe[0] = float(rematch.groups()[5])
            rep.right_fe[0] = float(rematch.groups()[6])
            rep.success[0] = rematch.groups()[7]
            rep.success_rate[0] = float(rematch.groups()[8])
            rematch = self.line_re.match(self.file.readline())
         # Finished with our first block of replicas. Bail out now
         break

   #================================================

   def _parse(self):
      """ Parses the rem.log file and loads the data arrays """
      rawline = self.file.readline()
      num_record = 1
      while rawline:
         rematch = self.line_re.match(rawline)
         if not rematch:
            rawline = self.file.readline()
            continue
         while rematch:
            n, m, t, p1, p2, lfe, rfe, suc, suc_rate = rematch.groups()
            n = int(n) - 1
            m = int(m) - 1
            self.reps[n].index[num_record] = n
            self.reps[n].neighbor_index[num_record] = m
            self.reps[n].temp[num_record] = float(t)
            self.reps[n].potene1[num_record] = float(p1)
            self.reps[n].potene2[num_record] = float(p2)
            self.reps[n].left_fe[num_record] = float(lfe)
            self.reps[n].right_fe[num_record] = float(rfe)
            self.reps[n].success[num_record] = suc
            self.reps[n].success_rate[num_record] = float(suc_rate)
            rematch = self.line_re.match(self.file.readline())
         num_record += 1
         rawline = self.file.readline()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class pHRemLog(RemLog):
   """ A class for a rem.log file in pH exchange """

   line_re = re.compile(r' *(\d+) *(\d+) *([-+]?\d+\.\d+) *(\d+\.\d+) *([-+]?\d+\.\d+)')

   #================================================

   def _get_replicas(self):
      " Gets all of the replica information from the first block of repinfo "
      rawline = self.file.readline()
      while rawline:
         rematch = self.line_re.match(rawline)
         if not rematch:
            rawline = self.file.readline()
            continue
         # Now we have our first block of replicas!
         while rematch:
            rep = Replica()
            self.reps.append(rep)
            # Each replica will have an index array to trace its trajectory
            # through replica space. Indexing starts from 0
            rep.index = [0 for i in range(self.numexchg)]
            rep.prot_cnt = np.zeros(self.numexchg)
            rep.success_rate = np.zeros(self.numexchg)
            rep.old_pH = np.zeros(self.numexchg)
            rep.new_pH = np.zeros(self.numexchg)
            rep.index[0] = int(rematch.groups()[0]) - 1
            rep.prot_cnt[0] = float(rematch.groups()[1])
            rep.old_pH[0] = float(rematch.groups()[2])
            rep.new_pH[0] = float(rematch.groups()[3])
            rep.success_rate[0] = float(rematch.groups()[4])
            self.values.append(rep.old_pH[0])
            rematch = self.line_re.match(self.file.readline())
         # Finished our first block of replicas. Bail out now
         break
      # Now sort our temperatures and adjust our index
      self.values.sort()
      for rep in self.reps:
         rep.index[0] = self.values.index(rep.old_pH[0])

   #================================================

   def _parse(self):
      """ Parses the rem.log file and loads the data arrays """
      rawline = self.file.readline()
      num_record = 1
      while rawline:
         rematch = self.line_re.match(rawline)
         if not rematch:
            rawline = self.file.readline()
            continue
         while rematch:
            repnum, prot_cnt, old_pH, new_pH, success_rate = rematch.groups()
            repnum = int(repnum) - 1
            self.reps[repnum].prot_cnt[num_record] = int(prot_cnt)
            self.reps[repnum].old_pH[num_record] = float(old_pH)
            self.reps[repnum].new_pH[num_record] = float(new_pH)
            self.reps[repnum].success_rate[num_record] = float(success_rate)
            self.reps[repnum].index[num_record] = self.values.index(float(old_pH))
            rematch = self.line_re.match(self.file.readline())
         num_record += 1
         rawline = self.file.readline()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

if __name__ == '__main__':
   from optparse import OptionParser
   import sys

   parser = OptionParser()
   parser.add_option('-l', '--remlog', dest='input', metavar='FILE',
                     help='Input rem.log file.', default=None)
   parser.add_option('-o', '--output', dest='output', metavar='FILE',
                     help='Output file', default=None)
   parser.add_option('-t', '--type', dest='type', default='TEMP',
                     help='Type of REM log file (TEMP/ph)')
   
   opt, args = parser.parse_args()

   if args or not opt.input:
      parser.print_help()
      sys.exit(1)

   if opt.type.upper() == 'TEMPERATURE'[:len(opt.type)] and len(opt.type) >= 4:
      opt.type = 'TEMP'
   elif opt.type.upper() == 'PH':
      opt.type = 'pH'

   if not opt.output: output = sys.stdout
   else: output = open(opt.output, 'w')

   if opt.type == 'TEMP':
      remlog = TempRemLog(opt.input)
   elif opt.type == 'pH':
      remlog = pHRemLog(opt.input)

   output.write('NUMEXCHG = %d\n' % remlog.numexchg)

   output.write('There are %d replicas\n\n' % len(remlog.reps))

   output.write("The final success rates are:\n")
   output.write("\n".join(['Val %d: %6.3f' % (rep.index[remlog.numexchg-1]+1,
                rep.success_rate[remlog.numexchg-1]) for rep in remlog.reps]))
   
   output.write('\n\n')
   output.write('The final resting positions are:\n')
   
   if opt.type == 'TEMP':
      output.write('\n'.join(["Rep %d: %.3fK"%(i+1,rep.new_temp[remlog.numexchg-1])
                   for i,rep in enumerate(remlog.reps)]))
   if opt.type == 'pH':
      output.write('\n'.join(["Rep %d: pH %.3f"%(i+1,rep.new_pH[remlog.numexchg-1]) 
                   for i,rep in enumerate(remlog.reps)]))
   output.write('\n')
