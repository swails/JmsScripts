#!/usr/bin/env python
from __future__ import division

"""
This script will calculate a pKa as a function of time fitting to a titration
curve
"""

import os, sys, re, math

# Make sure we have numpy and scipy
try:
   import numpy as np
   from scipy import optimize
   from scipy import __version__ as sp_version
except ImportError:
   print >> sys.stderr, 'Error: running_titration_curve.py depends on both ' + \
                        'numpy and scipy. Quitting'
   sys.exit(1)

# Make sure we have a version of scipy with curve_fit
if not hasattr(optimize, 'curve_fit'):
   print >> sys.stderr, ('Error: Your version of scipy [%s] does not have ' + \
         'curve_fit. You need at least scipy version 0.9.0') % sp_version
   sys.exit(1)

#-------------------------------------------------------------------------------

class pHStatFile(object):
   """ This loads a running pH statistics file """
   # Regular expressions for extracting info from the running pH stat file
   cumulativere = re.compile('========================== CUMULATIVE')
   chunkre = re.compile('============================ CHUNK')
   my_re = cumulativere # change this to chunkre if you want chunks
   solvphre = re.compile(r'Solvent pH is *([-+]?\d+\.\d+)')
   reslinere = re.compile(r'([A-Z4]+) *(\d+) *: Offset *([-+]?\d+(?:\.\d*)?|\.\d+|-*[Ii]nf) *Pred *([-+]?\d+(?:\.\d*)?|\.\d+|-*[Ii]nf) *Frac Prot *(\d+(?:\.\d*)?) *Transitions *(\d+)')
   
   def __init__(self, infile):
      """
      Constructor for pHStatFile:
         .  infile should be a "file" type object
      """
      self.f = infile
      # Search the file for the pH, then rewind the file
      self.pH = pHStatFile.get_pH(self)
      self.f.seek(0)
      # Get the list of titrating residues, then rewind the file
      self.list_of_residues = pHStatFile.titrating_residues(self)
      self.f.seek(0)

   def get_pH(self):
      """ Searches the beginning of the file for the pH """
      rawline = self.f.readline()
      while rawline:
         rematch = self.solvphre.match(rawline)
         if rematch:
            return float(rematch.groups()[0])
         rawline = self.f.readline()

   def titrating_residues(self):
      """
      Finds out which residues are titrating (i.e., have statistics printed in
      this file)
      """
      list_of_residues = []
      resname,resnum,offset,pred,frac_prot,trans = self.get_next_residue()
      while not '%s_%d' % (resname, resnum) in list_of_residues:
         list_of_residues.append('%s_%d' % (resname, resnum))
         resname,resnum,offset,pred,frac_prot,trans = self.get_next_residue()
      return list_of_residues

   def get_next_residue(self):
      """
      This command gets the next residue from this file.
      Returns a tuple: (resname, resnum, frac_prot)
      """
      rawline = self.f.readline()
      while rawline:
         rematch = self.reslinere.match(rawline)
         if rematch:
            return (rematch.groups()[0], int(rematch.groups()[1]),
                    float(rematch.groups()[2]), float(rematch.groups()[3]),
                    float(rematch.groups()[4]), int(rematch.groups()[5]))
         # If we make it to a blank line, keep skipping forward until we hit
         # another CUMULATIVE record
         elif not rawline.strip():
            rematch2 = self.my_re.match(rawline)
            while not rematch2:
               rawline = self.f.readline()
               # Check if we hit EOF
               if not rawline: return None
               rematch2 = self.my_re.match(rawline)
            # end while not rematch2
         rawline = self.f.readline()
      # If we hit here, we are out of lines, or something
      return None

#-------------------------------------------------------------------------------

def curve_with_hillcoef(ph, pka, hillcoef):
   """ Callable function with a variable Hill coefficient """
#  return hillcoef * ph - pka
   return 1/(1+10**(hillcoef*(pka-ph)))

#-------------------------------------------------------------------------------

def curve_no_hillcoef(ph, pka):
   """ Callable function with a fixed Hill coefficient of 1 """
#  return ph - pka
   return 1/(10**(pka-ph)+1)

#-------------------------------------------------------------------------------

def get_avg_pka(predicted_pkas):
   """ Gets average pKas, omitting the infinities """
   # Use a stupid (but effective in this case) test for infinity
   sm, np = 0, 0
   for val in predicted_pkas:
      if abs(val) > 99999999.0: continue
      sm += val
      np += 1
   if np == 0: return 1
   return sm / np

#-------------------------------------------------------------------------------

def main(file_list, outname, fit_func, starting_guess, chunk, hill):
   """
   This function is the main driver. It fits the data to the given fit_func (it
   should be one of the Callable functions defined above).
   
   Variable explanations:
      file_list:      List of input files with running titration data
      outname:        File prefix to dump all of the statistics to
      fit_func:       The function we're fitting to
      starting_guess: The starting guess for the parameters

   All error checking should be done on this input before calling main, or
   suffer the exceptions! Output files are named "outname_RES_NUM.dat"
   """
   
   class StopWhile(Exception): pass

   # See if we want to analyze chunks
   if chunk: pHStatFile.my_re = pHStatFile.chunkre

   xdata = np.zeros(len(file_list))
   ydata = np.zeros(len(file_list))
   # Convert the file_list into a list of pHStatFile objects if it's not yet
   if type(file_list[0]).__name__ == 'str':
      tmp = [pHStatFile(open(fname, 'r')) for fname in file_list]
      file_list = tmp
      del tmp
   # Build the list of output files
   output_files = {}
   for resid in file_list[0].list_of_residues:
      output_files[resid] = open('%s_%s.dat' % (outname, resid), 'w', 0)
  
   # Generate the x-data (the pHs). This never changes
   for i, frec in enumerate(file_list): xdata[i] = frec.pH

   # Now loop through all of our data
   numres = 0      # Number of residues we've looped through so far
   numframes = 0   # Number of frames we've looped through so far

   # This is the easiest way to bust out of an infinite loop -- Engulf the whole
   # thing in a big try-except, and catch a specialized exception.
   try:
      while True:
         numres += 1
         # If we've looped through all of our residues, then we know we've hit
         # the next frame, so update our counters accordingly
         if numres % len(output_files) == 0:
            numframes += 1
            numres = 1
         # Zero out the y-data, because we're about to fill it up
         ydata = np.zeros(len(file_list))           # fraction protonated
         offset = np.zeros(len(file_list))          # Offset for pKa
         pred = np.zeros(len(file_list))            # Predicted pKas
         trans = [0 for i in range(len(file_list))] # num of transitions
         # Loop through all of the files and get our next residue -- they should
         # be synchronized, so this should pull the same residue from each file
         for i, frec in enumerate(file_list):
            stuff = frec.get_next_residue()
            # If we got nothing bust out of the loop
            if not stuff:
               raise StopWhile
            resname,resnum,offset[i],pred[i],ydata[i],trans[i] = stuff
            ydata[i] = 1-ydata[i] # Get fraction DEprotonated
            # Make the y-data into a hill-plottable form
         if fit_func:
            # If we're doing a hill plot, adjust our starting guess to be
            # relatively close -- hill will start as 1, and pKa will start
            # as the average of pKa values (not including infinity)
            if hill:
               starting_guess = (get_avg_pka(pred), 1)
            try:
               params, cov = optimize.curve_fit(fit_func, xdata,
                                                ydata, starting_guess)
            except (RuntimeError, ValueError):
               # If we can't fit the data (expected at the beginning) just go on
               continue
            line = '%d ' % numframes
            try:
               for i, param in enumerate(params):
                  try:
#                    line += '%.4f %.4f ' % (param, math.sqrt(cov[i][i]))
                     line += '%.4f ' % param
                  except TypeError:
#                    line += '%.4f %.4f ' % (param, cov)
                     line += '%.4f ' % param
            except ValueError:
               continue
         else:
            # Average all of the predicted pKas, ignoring values whose offset is
            # >= 3 pH units
            runsum = runsum2 = numpts = 0
            for i in range(len(file_list)):
               if abs(offset[i]) < 3:
                  runsum += pred[i]
                  runsum2 += pred[i] * pred[i]
                  numpts += 1

            if numpts == 0: continue
            avg = runsum / numpts
            stdev = math.sqrt(abs(runsum2/numpts - avg*avg))
            line = '%d %.4f %.4f' % (numframes, avg, stdev)
            
         # Now write out the data as: Frame # pKa1 std.dev. [hill.coef. std.dev.]
         # but only write out if we actually got a pKa this time around
         ofile = output_files['%s_%d' % (resname, resnum)]
         ofile.write(line + os.linesep)

   except StopWhile: pass
   
if __name__ == '__main__':
   """ Main program """
   from optparse import OptionParser, OptionGroup
   from time import time
   
   start = time()
   usage = '%prog [Options] <pH_data1> <pH_data2> ... <pH_dataN>'
   epilog = ('This program will generate running pKa values with error bars ' +
             'taken from the quality of the fit.')

   parser = OptionParser(usage=usage, epilog=epilog)
   group = OptionGroup(parser, 'Fitting Options',
                       'These options control how the data are fit')
   group.add_option('-i', '--hill-coefficient', dest='hill', default=False,
                    action='store_true',
                    help='Fit allowing the Hill coefficient to change.' +
                    'Default behavior is to fix the Hill coefficient to 1.')
   group.add_option('-n', '--no-hill-coefficient', dest='hill',
                    action='store_false',
                    help='Fix the Hill coefficient to 1. Default behavior.')
   group.add_option('-a', '--average', dest='avg', default=False,
                    action='store_true',
                    help='Do simple averaging to get pKa and error bars.')
   parser.add_option_group(group)
   group = OptionGroup(parser, 'Data Options', 'These options control the ' +
                       'data that gets processed')
   group.add_option('-c', '--chunk', dest='chunk', default=False,
                    action='store_true', help='Analyze independent chunks ' +
                    'rather than the running pKa. NOT default behavior')
   group.add_option('-u', '--cumulative', dest='chunk', action='store_false',
                    help='Calculate running pKas. Default behavior')
   parser.add_option_group(group)
   group = OptionGroup(parser, 'Output Options', 'These options control output')
   group.add_option('-o', '--output', dest='outprefix', metavar='FILE_PREFIX',
                    default='running_pkas',
                    help='Prefix to apply to output files. The output files ' +
                    'will be named FILE_PREFIX_[resname]_[resnum].dat where ' +
                    'resname is the residue name of each titratable residue ' +
                    'and resnum is the residue number of that residue. ' +
                    'Default (%default)')
   parser.add_option_group(group)

   opt, args = parser.parse_args()

   # Check that we have input data files
   if not args:
      print >> sys.stderr, 'Error: Missing pH data files'
      parser.print_help()
      sys.exit(1)
   # Check that all files exist
   for fname in args:
      if not os.path.exists(fname):
         print >> sys.stderr, 'Error: File [%s] cannot be found!' % fname
         sys.exit(1)
   # Select our fitting function
   fit_func = curve_no_hillcoef
   starting_guess = 1
   if opt.hill:
      fit_func = curve_with_hillcoef
      starting_guess = (1,1)
   if opt.avg:
      fit_func = None
   # Now call our main function
   main(args, opt.outprefix, fit_func, starting_guess, opt.chunk, opt.hill)
   print 'Finished: I took %.3f minutes' % ((time()-start)/60)
