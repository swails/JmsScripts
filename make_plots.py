#! /usr/bin/env python

"""
This program will collect whatever statistics you want from given mdout files
and/or trajectory files and plots them.
"""

# Imports
import os
import numpy as np
from mdcrd import AmberTraj, RmsdData, TrajError, RMSError
from mdout import AmberMdout, MdoutError

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class PlotError(Exception): pass
class ImplementError(Exception):
   def __init__(self, msg='Not implemented yet!'): self.msg = msg
   def __str__(self): return self.msg

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def dump_1d_2d_data(fname, data, xcol=None, do_xcol=False, overwrite=False):
   """ 
   Dumps 1-D or 2-D data into a file.  If data is 2-D already, then the first
   dimension is assumed to be the X-axis and the second is assumed to be the
   Y-axis.  If it's 1-D, then the X-column can be (optionally) set to the set
   of natural numbers (1, 2, 3, ...) up to the number of elements in the data
   array. The output file is dumped in gnuplot-ready format
   """
   # Force do_xcol if an xcol was provided
   if not do_xcol and xcol: do_xcol = True
   # If we want to do the x-column, get the default
   if len(data.shape) == 1:
      if do_xcol:
         if not xcol:
            xcol = np.empty(data.shape[0])
            for i in range(len(data.shape[0])): xcol[i] = i+1
         # Stack them vertically
         data = np.vstack((xcol, data))

   if os.path.exists(fname) and not overwrite:
      raise PlotError('%s exists!' % fname)

   f = open(fname, 'w')

   for i in range(data.shape[1]):
      f.write('%s %s' % (data[0][i], data[1][i]) + os.linesep)
   f.close()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def gnuplot_2D_script(fname, datanames, xlabel='X values', ylabel='Y values',
                      key=None, lw=2, with_='lines', outfile_name=None,
                      overwrite=False, multiplot=False, size='1280,1024',
                      xrange=None, yrange=None, title='Plot'):
   """
   Writes a gnuplot script to plot given data:
      o  fname     : file name for gnuplot script
      o  datanames : all of the data files to plot (string or list)
      o  xlabel    : title of the x-axis
      o  ylabel    : title of the y-axis
      o  key       : location of the key (top right, bottom left, etc.), default
                     no key
      o  lw        : line width
      o  with_     : lines, points, linespoints?
      o  outfile_name : gif image to save (file name of that image). Forced gif
      o  overwrite : allow overwriting of existing files
      o  multiplot : if multiple datanames, put each in a separate plot?
      o  size      : how big to make the term
   """
   # Make sure our files don't exist
   if not overwrite:
      if os.path.exists(fname):
         raise PlotError('%s exists!' % fname)
      if outfile_name and os.path.exists(outfile_name):
         raise PlotError('%s exists!' % fname)
   # Open the file
   f = open(fname, 'w')
   if outfile_name:
      if not outfile_name.endswith('.gif'): outfile_name += '.gif'
      f.write('set term gif size \'%s\' enhanced font \'courier,20\'\n' % size)
      f.write('set output "%s"\n' % outfile_name)
   else:
      f.write('set term x11 size \'%s\' enhanced font \'courier,20\'\n' % size)
   # Do we have a key?
   if key:
      f.write('set key %s\n' % key)
   else:
      f.write('unset key\n')
   # Set the x label and y label and title and ranges
   f.write('set xl "%s"\n' % xlabel)
   f.write('set yl "%s"\n' % ylabel)
   f.write('set title "%s\'\n' % title)

   if xrange: 
      if not type(xrange).__name__ in ('list','tuple') and len(xrange) != 2:
         raise TypeError('xrange must be a list or tuple with 2 elements!')
      f.write('set xr %s\n' % list(xrange))
   if yrange: 
      if not type(yrange).__name__ in ('list','tuple') and len(yrange) != 2:
         raise TypeError('yrange must be a list or tuple with 2 elements!')
      f.write('set yr %s\n' % list(yrange))
   
   if type(datanames).__name__ == 'str':
      datanames = [datanames]

   if len(datanames) > 1 and multiplot:
      # Do multiplot
      raise ImplementError('Multiplot not implemented in gnuplot_2D_script yet')
   else:
      add_str = "plot "
      for i,dataname in enumerate(datanames):
         last = i == len(datanames) - 1
         add_str += "'%s' with %s" % (dataname, with_)
         if last: add_str += '\n'
         else: add_str += ', '
      f.write(add_str)
   f.close()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

