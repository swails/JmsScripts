"""
This is a generalization of the readparm.AmberParm class to handle similar
Amber-style files with %FLAG/%FORMAT tags
"""
import os, datetime
from chemistry.amber.readparm import (AmberParm, _parseFormat, 
                                      AMBER_ELECTROSTATIC)
from chemistry.exceptions import AmberParmError

class AmberFormat(AmberParm):
   """ 
   Generalization of the AmberParm class without some of the assumptions made
   about Amber topology files specifically
   """
   
   def __init__(self, fname = None):
      """ Constructor.  Read a file if given """
      self.parm_data = {}
      self.formats = {}
      self.parm_comments = {}
      self.flag_list = []
      self.version = ''
      self.prm_name = ''
      self.charge_flag = 'CHARGE'

      if fname is not None:
         self.rdparm(fname)

   def rdparm(self, fname):
      """ Parses the Amber format file """

      self.prm_name = fname
      current_flag = ''
      dat_type = ''
      gathering_data = False
      number_items_perline = 0
      size_item = 0
      self.version = '' # reset all topology information each time rdparm is called
      self.formats = {}
      self.parm_data = {}
      self.parm_comments = {}
      self.flag_list = []

      try: # open up the prmtop file, catching if it doesn't exist
         prmtop = open(self.prm_name, 'r')
      except IOError:
         self.exists = False
         self.valid = False
         return

      self.exists = True
      prmlines = prmtop.readlines() # load all lines into memory
      prmtop.close() # close the file now

      for i in range(len(prmlines)):

         if prmlines[i][0:8] == '%VERSION':
            self.version = prmlines[i].strip()

         elif prmlines[i][0:5] == '%FLAG':
            current_flag = prmlines[i][6:].strip()
            self.formats[current_flag] = ''
            self.parm_data[current_flag] = []
            self.parm_comments[current_flag] = []
            self.flag_list.append(current_flag)
            gathering_data = False

         elif prmlines[i][0:8] == '%COMMENT':
            self.parm_comments[current_flag].append(prmlines[i][9:].strip())

         elif prmlines[i][0:7] == '%FORMAT':
            self.formats[current_flag] = prmlines[i][8:len(prmlines[i].strip())-1]
            number_items_perline, size_item, dat_type, junk = _parseFormat(self.formats[current_flag])
            gathering_data = True

         elif gathering_data:
            position_inline = 0
            for j in range(number_items_perline):
               data_holder = prmlines[i][position_inline:position_inline+size_item].strip()
               if len(data_holder) == 0: # we've reached the end...
                  break # break out of the loop
               if dat_type == 'int': # if int, put held data item in as integer
                  self.parm_data[current_flag].append(int(data_holder))
               elif dat_type == 'dec': # if dec, put held data item in as float
                  self.parm_data[current_flag].append(float(data_holder))
               else: # otherwise, it stays a string
                  self.parm_data[current_flag].append(data_holder)
               position_inline += size_item # move to the next item

      # eliminate multiplicative constant on charges to reduce to fraction-e charges
      try:
         for i in range(len(self.parm_data[self.charge_flag])):
            self.parm_data[self.charge_flag][i] /= AMBER_ELECTROSTATIC
      except KeyError:
         return
      
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
      if self.version != '':
         now = datetime.now()
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
         number_items_perline, size_item, dat_type, decnum = _parseFormat(self.formats[flag])
         line = ''
         num_items = 0
         if len(self.parm_data[flag]) == 0: # empty field...
            new_prm.write('\n')
            continue
         for j in range(len(self.parm_data[flag])): # write data in new_prm
            if dat_type == 'dec' and 'E' in self.formats[flag].upper():
               line += ('%%%s.%sE' % (size_item, decnum)) % self.parm_data[flag][j] 
            elif dat_type == 'dec' and 'F' in self.formats[flag].upper():
               line += ('%%%s.%sF' % (size_item, decnum)) % self.parm_data[flag][j] 
            elif dat_type == 'int':
               line += ('%%%sd' % size_item) % self.parm_data[flag][j] 
            else:
               line += ('%s' % self.parm_data[flag][j]).ljust(size_item)

            num_items += 1
            if num_items == number_items_perline: # flush line to prmtop
               new_prm.write(line + '\n')
               line = ''
               num_items = 0
         if len(line.strip()) > 0:
            new_prm.write(line + '\n') # flush what's left to the prmtop

      new_prm.close() # close new prmtop

      if self.charge_flag in self.parm_data.keys():
         # Convert charges back to electron-units
         for i in range(len(self.parm_data[self.charge_flag])):
            self.parm_data[self.charge_flag][i] /= AMBER_ELECTROSTATIC
