"""
Read the mdout file from an H-REMD simulation (not general)
"""
from __future__ import division
import numpy as np
import re
from remd import Replica

class HRemMdout(object):
   """ Mdout class for H-REMD compiled with -DVERBOSE_REMD """
   def __init__(self, file_list):
      self.reps = [Replica() for j in range((len(file_list)))]
      self.numexchg = self._find_numexchg(open(file_list[0]))
      for i, f in enumerate(file_list):
         self._parse_file(open(f, 'r'), i)

   def _find_numexchg(self, f):
      """ Finds the number of exchange attempts """
      myre = re.compile(r'.*numexchg= *(\d+)')
      for line in f:
         rematch = myre.match(line)
         if rematch:
            return int(rematch.groups()[0])

   def _parse_file(self, f, i):
      """ Parses the H-REMD mdout file with -DVERBOSE_REMD """
      self.reps[i].index = [i for j in range((self.numexchg))]
      self.reps[i].neighbor_index = [0 for j in range((self.numexchg))]
      self.reps[i].temp = np.zeros(self.numexchg)
      self.reps[i].potene1 = np.zeros(self.numexchg)
      self.reps[i].potene2 = np.zeros(self.numexchg)
      self.reps[i].left_fe = np.zeros(self.numexchg)
      self.reps[i].right_fe = np.zeros(self.numexchg)
      self.reps[i].success = ['F' for j in range((self.numexchg))]

      ene1_re = re.compile(r'My Eptot_1: *([-+]?\d*\.\d+)')
      ene2_re = re.compile(r'My Eptot_2: *([-+]?\d*\.\d+)')
      dir_re = re.compile(r'Jumping (Left|Right)')
      
      rawline = f.readline()
      cnt = 0
      while rawline:
         rematch = ene1_re.match(rawline)
         if rematch:
            self.reps[i].potene1[cnt] = float(rematch.groups()[0])
            rematch = ene2_re.match(f.readline())
            self.reps[i].potene2[cnt] = float(rematch.groups()[0])
            # Eat the next 2 lines
            f.readline()
            f.readline()
            rematch = dir_re.match(f.readline())
            if rematch.groups()[0] == 'Left':
               if i == 0:
                  self.reps[i].neighbor_index[cnt] = len(self.reps) - 1
               else:
                  self.reps[i].neighbor_index[cnt] = i - 1
            else:
               if i == len(self.reps) - 1:
                  self.reps[i].neighbor_index[cnt] = 0
               else:
                  self.reps[i].neighbor_index[cnt] = i + 1

            cnt += 1
         rawline = f.readline()
