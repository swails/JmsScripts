"""
This module contains code useful for parsing energy terms out of an mdout file
and exposing the data (numpy arrays) for various analyses. It will also figure
out what kind of simulation is being run based on the input variables

It does NOT do:

   o  Store decomp data when idecomp != 0
"""

from optparse import OptionParser
import numpy as np
import re

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class AmberMdout(object):
   """ Class for Mdout object """
   # If we don't know how many terms we have, just make it a million, since
   # that's (should be!) more data than we could possibly have. We can reshape
   # it afterwards
   UNKNOWN = 1000000

   #================================================

   def __init__(self, filename):
      self.filename = filename # File name of mdout file
      self.data = {}           # Dict that matches mdout term names to arrays
      self.properties = {}     # List of input variables with their values
      self.get_properties()    # Get properties of output file
      keys = self.properties.keys()
      # Figure out if we're a minimization or an MD
      if 'imin' in keys:
         self.is_md = self.properties['imin'] == 0 # This is MD iff imin == 0
         self.is_min = not self.is_md
      else:
         self.is_min = True
         self.is_md = False
      
      # Figure out if we're a restart or not
      if 'irest' in keys:
         self.is_restart = self.properties['irest'] != 0
      else:
         self.is_restart = False

      # Determine how many steps we've done
      if (self.is_min and 'maxcyc' in keys)or(self.is_md and 'nstlim' in keys):
         if self.is_min: self.num_steps = self.properties['maxcyc']
         if self.is_md: self.num_steps = self.properties['nstlim']
         # If maxcyc and nstlim are in properties, then ntpr HAS to be
         self.num_terms = self.num_steps / self.properties['ntpr'] + 1
         # For restart, we don't have that extra term at the beginning
         if self.is_restart: self.num_terms -= 1
      else:
         self.num_steps = AmberMdout.UNKNOWN
         self.num_terms = AmberMdout.UNKNOWN
      self.get_data()
      
   #================================================

   def get_properties(self):
      """ Searches through and finds properties """
      fl = open(self.filename, 'r')
      rawline = fl.readline()
      get_prop = re.compile(r"""(\w+ *= *[\w\.'"]+),*""")
      amber_comment = re.compile(r'^\|')
      while rawline:
         # Get up to where we want to go
         if rawline[:35] != '   2.  CONTROL  DATA  FOR  THE  RUN':
            rawline = fl.readline()
            continue
         # Skip over comments
         if amber_comment.match(rawline): 
            rawline = fl.readline()
            continue
         # Now we start reading in the 
         while rawline:
            # If we've reached the end, stop
            if rawline[:40] == '   3.  ATOMIC COORDINATES AND VELOCITIES':
               break
            # Ignore amber comments
            if amber_comment.match(rawline):
               rawline = fl.readline()
               continue
            items = get_prop.findall(rawline)
            # Now get all of the individual items and add them to the prop dict
            for item in items:
               var, val = item.split('=')
               var, val = var.strip(), val.strip()
               # See if we can make it an int, then a float, or just a string
               try:
                  self.properties[var] = int(val)
               except ValueError:
                  # Now try float
                  try:
                     self.properties[var] = float(val)
                  except ValueError:
                     # Finally, it's just a string
                     self.properties[var] = str(val)
               # end except
            # Get the next line
            rawline = fl.readline()
         # We only reach here if we broke out of the while loop above
         break
      # end while rawline
      # Close out the file now
      fl.close()

   #================================================

   def get_data(self):
      """ Extracts the data from the mdout file """
      energy_fields = re.compile(r'([A-Z\(\) 1-9\-]+(?:tot){0,1} *= *\-*\d+[.\d]*)')
      if self.is_md:
         start_of_record = energy_fields
         ignore_record = re.compile(r'^      A V E R A G E S   O V E R|^      R M S  F L U C T U A T I O N S')
      elif self.is_min:
         start_of_record = re.compile(r'^   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER')
         ignore_record = re.compile(r'^                    FINAL RESULTS')
      fl = open(self.filename, 'r')
      rawline = fl.readline()
      first_record_done = False
      num_record = 0
      ignore_this_record = False

      while rawline:
         # First find the Results section
         if rawline[:14] != '   4.  RESULTS':
            rawline = fl.readline()
            continue
         # Now that we've found our results section, start parsing it
         while rawline and rawline[:14] != '   5.  TIMINGS':
            # See if we hit a line that says we need to ignore the next group
            # of terms. We construct it like this since there are blank lines
            # between the ignore marker and the terms, so once the re.match
            # becomes true (and toggles ignore_this_record to True), we don't
            # want future misses of this re.match to turn this to False. All
            # is set well for the next record by setting ignore_this_record to
            # False at the end of the while loop
            ignore_this_record = ignore_this_record or \
                                 bool(ignore_record.match(rawline))
            # See if we found the start of the record
            items = start_of_record.findall(rawline)
            # If items is blank, read the next line and continue
            if not items:
               rawline = fl.readline()
               continue
            # If we are doing a minimization, the format is different -- the
            # starting line is by itself, and their values are right below them
            if self.is_min:
               terms = items[0].split()
               if not first_record_done:
                  for term in terms: 
                     # NAME is an atom name -- skip over this one
                     if term == 'NAME': continue
                     self.data[term] = np.zeros(self.num_terms)
               term_vals = fl.readline().split()
               for i, term in enumerate(terms):
                  # Skip over NAME again (see above)
                  if term == 'NAME': continue
#                 print num_record, len(self.data[term]), 'vs.', i, len(term_vals)
                  if not ignore_this_record:
                     self.data[term][num_record] = float(term_vals[i])
               # Eat the next line (it's blank)
               fl.readline()
               # The next line has stuff on it
               rawline = fl.readline()
               items = energy_fields.findall(rawline)
            # Now if we've found a record, cycle through it
            while items:
               # Load the items into their respective places in the data dict
               for item in items:
                  term, term_val = item.split('=')
                  term, term_val = term.strip(), term_val.strip()
                  # If we're on our first one, we need to create the data array
                  if not first_record_done:
                     self.data[term] = np.zeros(self.num_terms)
                  if not ignore_this_record:
                     self.data[term][num_record] = float(term_val)
               # Now we're done with the terms, get the next line
               rawline = fl.readline()
               items = energy_fields.findall(rawline)
            # Now that we're (definitely) done with our first record...
            first_record_done = True
            # We do not want to off-hand ignore the next record
            if not ignore_this_record: num_record += 1
            ignore_this_record = False
            rawline = fl.readline()
         # end while rawline and rawline[:14] != '   5.  TIMINGS':

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

if __name__ == '__main__':
   """ Test our code """
   parser = OptionParser()
   parser.add_option('--mdout', dest='mdout', help='Input mdout to analyze',
                     default=None)
   opt, args = parser.parse_args()
   if not opt.mdout:
      parser.print_help()
      quit()

   my_mdout = AmberMdout(opt.mdout)
   longest_option = max([len(it) for it in my_mdout.properties.keys()])
   print 'Testing get_properties() (Check with mdout file)'
   print '-'*longest_option + '------------'
   for prop in my_mdout.properties.keys():
      print '  %%%ds | %%s' % longest_option % (prop, my_mdout.properties[prop])

   print ''
   if my_mdout.is_md:
      print 'This mdout file is a MOLECULAR DYNAMICS calculation'
   if my_mdout.is_min:
      print 'This mdout file is a MINIMIZATION calculation'
   if my_mdout.is_restart:
      print 'This mdout file is from a RESTARTED simulation'
   if not my_mdout.is_restart:
      print 'This mdout file is NOT a restarted simulation'

   print 'Testing get_data()'
   print 'Found data keys:', ', '.join(my_mdout.data.keys())
   print 'Averages of each are:'
   keys = my_mdout.data.keys()
   keys.sort()
   for key in keys:
      print '   %10s : %15.4f' % (key, np.average(my_mdout.data[key]))
   print '\nDone testing.'
