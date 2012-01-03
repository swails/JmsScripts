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

class TempRemLog(object):
   """ Replica exchange log file """

   templine_re = re.compile(r'\d+ *([-+]?\d\.\d{2}) *(\d+\.\d{2}) *([-+]?\d+\.\d{2}) *(\d+\.\d{2}) *(\d+\.\d{2}) *(\d\.\d{2}) *([+-]?\d+)')
   numexch_desc_re = re.compile(r'#numexchg is *(\d+)')

   #================================================

   def __init__(self, fname):
      """ Loads a replica exchange log file """
      # Properties corresponding to each replica
      self.props = [0 for i in range(self.numexchg)]

      self.fname = fname          # Name of the REMD file
      self.reps = []              # List of all replicas
      self.numexchg = 0           # How many exchanges we do
      self.file = open(fname, 'r')# File object for rem.log
      self._get_numexchg()        # Parse the numexchg out and return
      self._get_replicas()        # Set self.reps from first section
      self._parse()               # Get all replica information

   #================================================

   def _get_numexchg(self):
      """ Get the number of exchanges """
      for line in self.file:
         rematch = self.numexchg_desc_re.match(line)
         if rematch:
            self.numexchg = int(rematch.groups()[0])
            return
      raise RemdError('Could not find #numexchg. Bad remlog file [%s]' %
                      self.fname)

   #================================================

   def _get_replicas(self):
      " Gets all of the replica information from the first block of repinfo "
      rawline = self.file.readline()
      while rawline:
         rematch = self.templine_re.match(rawline)
         if not rematch:
            rawline = self.file.readline()
            continue
         # Now we have our first block of replicas!
         while rematch:
            rep = Replica()
            self.reps.append(rep)
            # Each replica will have an index array to trace its trajectory
            # through replica space. Indexing starts from 0
            rep.index = np.zeros(self.numexchg)
            rep.potene = np.zeros(self.numexchg)
            rep.success_rate = np.zeros(self.numexchg)
            rep.temp = np.zeros(self.numexchg)
            rep.index[0] = int(rematch.groups()[0]) - 1
            rep.potene[0] = float(rematch.groups()[3])
            rep.success_rate[0] = float(rematch.groups()[6])
            rep.temp[0] = float(rematch.groups()[4])
            self.values[rep.index[0]] = rep.temp[0]
            rematch = self.templine_re.match(self.file.readline())
         # Finished our first block of replicas. Bail out now
         break
      # Now sort our temperatures and adjust our index
      self.values.sort()
      for rep in self.reps:
         rep.index[0] = self.values.index(rep.temp[0])

   #================================================

   def _parse(self):
      """ Parses the rem.log file and loads the data arrays """
      rawline = self.file.readline()
      num_record = 1
      while rawline:
         rematch = self.templine_re.match(rawline)
         if not rematch:
            rawline = self.file.readline()
            continue
         while rematch:
            repnum,j1,j2,potene,temp,j3,success_rate,j4 = rematch.groups()
            repnum = int(repnum)
            potene = float(potene)
            temp = float(temp)
            success_rate = float(success_rate)
            self.reps[repnum].potene[num_record] = potene
            self.reps[repnum].temp[num_record] = temp
            self.reps[repnum].success_rate[num_record] = success_rate
            self.reps[repnum].index[num_record] = self.values.index(temp)
            rematch = self.templine_re.match(self.file.readline())
         num_record += 1

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

