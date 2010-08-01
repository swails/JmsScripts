from sys import stderr, stdout
from datetime import datetime

# Global constants
AMBER_ELECTROSTATIC = 18.2223

class amberParm:

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def __init__(self, prm_name): # set up necessary variables

      self.formats = {}        # List of formats associated with each section
      self.parm_data = {}      # Master dictionary with all of the prmtop data
      self.flag_list = []      # List of all of the flags in the prmtop
      self.version = ''        # version of the prmtop file
      self.prm_name = prm_name # name of the prmtop file
      self.overwrite = False   # whether writeparm will overwrite prm_name
      self.prm_exists = False  # whether or not prm_name exists

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def parseFormat(self, format_string):  # parse a format statement and send back details

      if 'a' in format_string: # this is a string
         format_parts = format_string.split('a')
         return int(format_parts[0]), int(format_parts[1]), 'str'

      elif 'I' in format_string: # this is an integer
         format_parts = format_string.split('I')
         return int(format_parts[0]), int(format_parts[1]), 'int'

      elif 'E' in format_string: # this is a floating point decimal
         format_parts = format_string.split('E')
         decimal_parts = format_parts[1].split('.')
         return int(format_parts[0]), int(decimal_parts[0]), 'dec'

      else:
         print >> stderr, 'Error: Unrecognized format "{0}"!'.format(format_string)
         return 1, 80, 'str'

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def rdparm(self):   # read topology file and load all data into arrays/dictionaries

      # necessary variables
      current_flag = ''
      dat_type = ''
      gathering_data = False
      number_items_perline = 0
      size_item = 0
      self.version = '' # reset version each time we read the topology file

      try: # open up the prmtop file, catching if it doesn't exist
         prmtop = open(self.prm_name, 'r')
      except IOError:
         print >> stderr, 'Error: %s does not exist!' % self.prm_name
         return

      self.prm_exists = True
      prmlines = prmtop.readlines() # load all lines into memory
      prmtop.close() # close the file now

      for i in range(len(prmlines)):

         if prmlines[i][0:8] == '%VERSION':
            if self.version != '':
               print >> stderr, 'Warning: %VERSION string defined multiple times.'
            self.version = prmlines[i].strip()

         elif prmlines[i][0:5] == '%FLAG':
            current_flag = prmlines[i][6:].strip()
            self.formats[current_flag] = ''
            self.parm_data[current_flag] = []
            self.flag_list.append(current_flag)
            gathering_data = False

         elif prmlines[i][0:7] == '%FORMAT':
            self.formats[current_flag] = prmlines[i][8:len(prmlines[i].strip())-1]
            number_items_perline, size_item, dat_type = self.parseFormat(self.formats[current_flag])
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
      for i in range(len(self.parm_data["CHARGE"])):
         self.parm_data["CHARGE"][i] /= AMBER_ELECTROSTATIC

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def writeParm(self, name):   # write a new prmtop with the current prmtop data
      # make sure we want to write the new prmtop file
      if not self.overwrite and name == self.prm_name:
         print >> stderr, 'Error: Object\'s overwrite set to False! Will not overwrite original prmtop.'
         return
      elif self.version == '':
         print >> stderr, 'Error: Version string missing! Load prmtop data before writing new prmtop.'
         return

      # now that we know we will write the new prmtop file, open the new file
      new_prm = open(name, 'w')

      # get current time to put into new prmtop file
      now = datetime.now()
      datestring = "DATE = {0:02d}/{1:02d}/{2:02d}  {3:02d}:{4:02d}:{5:02d}".format(
                  now.month, now.day, now.year % 100, now.hour, now.minute, now.second)
      for i in range(len(self.version)): # replace the date in version string
         if self.version[i:i+2] == "DA":
            self.version = self.version[:i] + datestring
            break

      # convert charges back to amber charges...
      for i in range(len(self.parm_data["CHARGE"])):
         self.parm_data["CHARGE"][i] *= AMBER_ELECTROSTATIC

      # write version to top of prmtop file
      new_prm.write('{0}\n'.format(self.version))

      # write data back to prmtop file, inserting blank line if it's an empty field
      for i in range(len(self.flag_list)):
         flag = self.flag_list[i]
         new_prm.write('%FLAG {0}\n'.format(flag))
         new_prm.write('%FORMAT({0})\n'.format(self.formats[flag]))
         number_items_perline, size_item, dat_type = self.parseFormat(self.formats[flag])
         if dat_type == 'dec':
            decnum = int(self.formats[flag].split('E')[1].split('.')[1])
         line = ''
         num_items = 0
         if len(self.parm_data[flag]) == 0: # empty field...
            new_prm.write('\n')
            continue
         for j in range(len(self.parm_data[flag])): # write data in new_prm
            if dat_type == 'dec':
               line += '{0:{1}.{2}e}'.format(self.parm_data[flag][j],size_item,decnum)
            elif dat_type == 'int':
               line += '{0:{1}d}'.format(self.parm_data[flag][j],size_item)
            else:
               line += '{0}'.format(self.parm_data[flag][j]).ljust(size_item)

            num_items += 1
            if num_items == number_items_perline: # flush line to prmtop
               new_prm.write(line + '\n')
               line = ''
               num_items = 0
         new_prm.write(line + '\n') # flush what's left to the prmtop

      new_prm.close() # close new prmtop

      # eliminate multiplicative constant on charges to reduce to fraction-e charges
      for i in range(len(self.parm_data["CHARGE"])):
         self.parm_data["CHARGE"][i] /= AMBER_ELECTROSTATIC


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
