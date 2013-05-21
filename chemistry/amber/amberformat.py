"""
This is a generalization of the readparm.AmberParm class to handle similar
Amber-style files with %FLAG/%FORMAT tags
"""
from __future__ import division

from array import array
from chemistry.amber.constants import AMBER_ELECTROSTATIC
from chemistry.exceptions import AmberParmError, AmberFormatWarning, FlagError
import datetime
from math import ceil
import os
from warnings import warn

class FortranFormat(object):
   """ Handles fortran formats """

   #===================================================

   def __init__(self, format_string, strip_strings=True):
      """
      Sets the format string and determines how we will read and write strings
      using this format
      """
      self.format = format_string

      # Define a function that processes all arguments prior to adding them to
      # the returned list. By default, do nothing, but this allows us to
      # optionally strip whitespace from strings.
      self.process_method = lambda x: x

      if 'a' in format_string.lower():
         # replace our write() method with write_string to force left-justify
         self.type, self.write = str, self.write_string
         try:
            self.nitems, self.itemlen = format_string.split('a')
            self.nitems, self.itemlen = int(self.nitems), int(self.itemlen)
         except ValueError:
            self.nitems, self.itemlen = 1, 80
         self.fmt = '%s'
         # See if we want to strip the strings
         if strip_strings: self.process_method = lambda x: x.strip()

      elif 'I' in format_string.upper():
         self.type = int
         if len(format_string.upper().split('I')[0]) == 0:
            self.nitems = 1
         else:
            self.nitems = int(format_string.upper().split('I')[0])
         self.itemlen = int(format_string.upper().split('I')[1])
         self.fmt = '%%%dd' % self.itemlen

      elif 'E' in format_string.upper():
         self.type = float
         format_parts = format_string.upper().split('E')
         if len(format_parts[0]) == 0:
            self.nitems = 1
         else:
            self.nitems = int(format_parts[0])
         self.itemlen = int(format_parts[1].split('.')[0])
         self.num_decimals = int(format_parts[1].split('.')[1])
         self.fmt = '%%%s.%sE' % (self.itemlen, self.num_decimals)

      elif 'F' in format_string.upper():
         self.type = float
         format_parts = format_string.upper().split('F')
         if len(format_parts[0].strip()) == 0:
            self.nitems = 1
         else:
            self.nitems = int(format_parts[0])
         self.itemlen = int(format_parts[1].split('.')[0])
         self.num_decimals = int(format_parts[1].split('.')[1])
         self.fmt = '%%%s.%sF' % (self.itemlen, self.num_decimals)

      else:
         # replace our write() method with write_string to force left-justify
         self.type, self.write = str, self.write_string
         warn('Unrecognized format "%s". Assuming string.' % format_string,
              AmberFormatWarning)
         self.fmt = '%s'
         self.nitems, self.itemlen = 1, 80
         # See if we want to strip the strings
         if strip_strings: self.process_method = lambda x: x.strip()

   #===================================================

   def __str__(self):
      return self.format

   #===================================================

   def write(self, items, dest):
      """ Writes a list/tuple of data (or a single item) """
      if isinstance(items, list) or isinstance(items, tuple) or \
            isinstance(items, array):
         mod = self.nitems - 1
         for i, item in enumerate(items):
            dest.write(self.fmt % item)
            if i % self.nitems == mod:
               dest.write('\n')
         if i % self.nitems != mod:
            dest.write('\n')
      else:
         dest.write(self.fmt % item)
         dest.write('\n')

   #===================================================

   def write_string(self, items, dest):
      """ Writes a list/tuple of strings """
      if isinstance(items, list) or isinstance(items, tuple) or \
            isinstance(items, array):
         mod = self.nitems - 1
         for i, item in enumerate(items):
            dest.write((self.fmt % item).ljust(self.itemlen))
            if i % self.nitems == mod:
               dest.write('\n')
         if i % self.nitems != mod:
            dest.write('\n')
      else:
         dest.write((self.fmt % item).ljust(self.itemlen))
         dest.write('\n')

   #===================================================

   def read(self, line):
      """ Reads the line and returns the converted data """
      line = line.rstrip()
      ret = [0 for i in range(int(ceil(len(line) / self.itemlen)))]
      start, end = 0, self.itemlen
      for i in range(len(ret)):
         ret[i] = self.process_method(self.type(line[start:end]))
         start += self.itemlen
         end += self.itemlen
      return ret

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AmberFormat(object):
   """ 
   Generalization of the AmberParm class without some of the assumptions made
   about Amber topology files specifically
   """
   
   #===================================================

   def __init__(self, fname=None):
      """ Constructor.  Read a file if given """
      self.parm_data = {}
      self.formats = {}
      self.parm_comments = {}
      self.flag_list = []
      self.version = None
      self.prm_name = fname
      self.charge_flag = 'CHARGE'
      self.valid = False

      if fname is not None:
         self.rdparm(fname)

   #===================================================

   def rdparm(self, fname):
      """ Parses the Amber format file """

      self.prm_name = fname
      current_flag = ''
      gathering_data = False
      self.version = None # reset all top infor each time rdparm is called
      self.formats = {}
      self.parm_data = {}
      self.parm_comments = {}
      self.flag_list = []

      # Open up the file and read the data into memory
      prm = open(self.prm_name, 'r')
      self.valid = False

      for line in prm:

         if line[0:8] == '%VERSION':
            self.version = line.strip()

         elif line[0:5] == '%FLAG':
            current_flag = line[6:].strip()
            self.formats[current_flag] = ''
            self.parm_data[current_flag] = []
            self.parm_comments[current_flag] = []
            self.flag_list.append(current_flag)
            gathering_data = False

         elif line[0:8] == '%COMMENT':
            self.parm_comments[current_flag].append(line[9:].strip())

         elif line[0:7] == '%FORMAT':
            fmt = FortranFormat(line[8:len(line.strip())-1])
            self.formats[current_flag] = fmt
            gathering_data = True

         elif gathering_data:
            self.parm_data[current_flag].extend(fmt.read(line))

      # convert charges to fraction-electrons
      try:
         for i, chg in enumerate(self.parm_data[self.charge_flag]):
            self.parm_data[self.charge_flag][i] = chg / AMBER_ELECTROSTATIC
      except KeyError:
         pass
      self.valid = True
      return
      
   #===================================================

   def set_version(self, verstring):
      """ Sets the version string """
      self.version = '%%VERSION %s' % verstring.strip()

   #===================================================

   def writeParm(self, name, overwrite=False):
      """
      Writes the current data in parm_data into a new topology file with
      the given name
      """
      # global variable(s)
      global AMBER_ELECTROSTATIC

      # make sure we want to write the new prmtop file
      if not overwrite and os.path.exists(name):
         raise AmberParmError('%s exists. Not overwriting' % name)

      # now that we know we will write the new prmtop file, open the new file
      new_prm = open(name, 'w')

      # get current time to put into new prmtop file if we had a %VERSION
      now = datetime.datetime.now()
      datestring = "DATE = %02d/%02d/%02d  %02d:%02d:%02d" % ( now.month, 
                   now.day, now.year % 100, now.hour, now.minute, now.second)
      for i in range(len(self.version)): # replace the date in version string
         if self.version[i:i+2] == "DA":
            self.version = self.version[:i] + datestring
            break

      # convert charges back to amber charges...
      if self.charge_flag in self.parm_data.keys():
         for i in range(len(self.parm_data[self.charge_flag])):
            self.parm_data[self.charge_flag][i] *= AMBER_ELECTROSTATIC

      # write version to top of prmtop file
      new_prm.write('%s\n' % self.version)

      # write data back to prmtop file, inserting blank line if it's an empty field
      for i in range(len(self.flag_list)):
         flag = self.flag_list[i]
         new_prm.write('%%FLAG %s\n' % flag)
         # Insert any comments before the %FORMAT specifier
         for comment in self.parm_comments[flag]:
            new_prm.write('%%COMMENT %s\n' % comment)
         new_prm.write('%%FORMAT(%s)\n' % self.formats[flag])
         if len(self.parm_data[flag]) == 0: # empty field...
            new_prm.write('\n')
            continue
         self.formats[flag].write(self.parm_data[flag], new_prm)

      new_prm.close() # close new prmtop

      if self.charge_flag in self.parm_data.keys():
         # Convert charges back to electron-units
         for i in range(len(self.parm_data[self.charge_flag])):
            self.parm_data[self.charge_flag][i] /= AMBER_ELECTROSTATIC

   #===================================================

   def addFlag(self, flag_name, flag_format, data=None, num_items=-1,
               comments=[]):
      """
      Adds a new flag with the given flag name and Fortran format string and
      initializes the array with the values given, or as an array of 0s
      of length num_items
      """
      self.flag_list.append(flag_name.upper())
      self.formats[flag_name.upper()] = FortranFormat(flag_format)
      if data is not None:
         self.parm_data[flag_name.upper()] = list(data)
      else:
         if num_items < 0:
            raise FlagError("If you do not supply prmtop data, num_items must "
                            "be non-negative!")
         self.parm_data[flag_name.upper()] = [0 for i in range(num_items)]
      if comments:
         if isinstance(comments, str):
            comments = [comments]
         elif isinstance(comments, tuple):
            comments = list(comments)
         elif isinstance(comments, list):
            pass
         else:
            raise TypeError('Comments must be string, list, or tuple')
         self.parm_comments[flag_name.upper()] = comments
      else:
         self.parm_comments[flag_name.upper()] = []

   #===================================================

   def deleteFlag(self, flag_name):
      """ Removes a flag from the topology file """
      flag_name = flag_name.upper()
      if not flag_name in self.flag_list:
         return # already gone
      del self.flag_list[self.flag_list.index(flag_name)]
      del self.parm_comments[flag_name]
      del self.formats[flag_name]
      del self.parm_data[flag_name]

