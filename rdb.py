"""
This module contains code useful for reading, writing, and manipulating RDB
database tables (i.e., tab-delimited data fields with a header and format line)
"""
from __future__ import division, with_statement, print_function

import array
from collections import OrderedDict
from functools import wraps
from math import sqrt
import re
from sys import stderr
import warnings

try:
   from numpy import sqrt
   import numpy as np
   _HAS_NUMPY = True
except ImportError:
   np = None
   _HAS_NUMPY = False

def checksize(fcn):
   @wraps(fcn)
   def newfcn(self, other):
      if len(self) != len(other):
         raise SizeError('Vector sizes incompatible')
      return fcn(self, other)
   return newfcn

def check_vector_or_scalar(fcn):
   @wraps(fcn)
   def newfcn(self, other):
      scalar = False
      if len(other) == 1: # Scalar
         scalar = True
      elif len(self) != len(other):
         raise SizeError('Vector sizes incompatible')
      return fcn(self, other)
   return newfcn

class Vector(array.array):
   # Make vector addition, subtraction work the way it's expected to
   @checksize
   def __add__(self, other):
      return Vector(self.typecode,
                    [self[i] + other[i] for i in range(len(self))])

   @checksize
   def __sub__(self, other):
      return Vector(self.typecode,
                    [self[i] - other[i] for i in range(len(self))])

   @check_vector_or_scalar
   def __mul__(self, other):
      global scalar
      if scalar: return Vector([i*other for i in self])
      return Vector(self.typecode,
                    [self[i] * other[i] for i in range(len(self))])

   @check_vector_or_scalar
   def __div__(self, other):
      global scalar
      if scalar: return Vector([i/other for i in self])
      return Vector(self.typecode,
                    [self[i] / other[i] for i in range(len(self))])

   @checksize
   def __iadd__(self, other):
      for i in range(len(self)):
         self[i] += other[i]
      return self

   @checksize
   def __isub__(self, other):
      for i in range(len(self)):
         self[i] -= other[i]
      return self

   @check_vector_or_scalar
   def __idiv__(self, other):
      global scalar
      if scalar:
         for i in range(len(self)): self[i] /= other
      else:
         for i in range(len(self)): self[i] /= other[i]
      return self

   @check_vector_or_scalar
   def __imul__(self, other):
      global scalar
      if scalar:
         for i in range(len(self)): self[i] *= other
      else:
         for i in range(len(self)): self[i] *= other[i]
      return self

   def mean(self):
      """ Return the mean of this list """
      return sum(self) / len(self)

   def std(self):
      """ Return the standard deviation of this list """
      sum = sum2 = 0
      for v in self:
         sum += v
         sum2 += v * v
      avg = sum / len(self)
      return sqrt(abs(sum2 / len(self) - avg * avg))

# Python 3 support
try:
   unicode
   str = unicode
except NameError:
   pass # Nothing to do

# Exception hierarchy
class RDBError(Exception): pass
class SizeError(RDBError): pass

# Decorators
def _ensurefile(mode, clsmeth=False):
   """
   Makes sure that the first object passed to the function is an open file
   with the requested mode
   """
   def wrapper(fcn):
      @wraps(fcn)
      def newfcn(self, f, *args, **kwargs):
         if not hasattr(f, 'write') or not hasattr(f, 'read'):
            # Pretend it's a string
            with open(f, mode) as _f:
               return fcn(self, _f, *args, **kwargs)
         return fcn(self, f, *args, **kwargs)
      return newfcn
   return wrapper

class RDB(OrderedDict):
   """ RDB data """

   _fmtre = re.compile(r'%(?:[0-9-\.]*)([sgfdi])')
   def __init__(self, *args, **kwargs):
      super(RDB, self).__init__(*args, **kwargs)
      self.formats = []
      self._types = []
      self._format = {str : '%s', float : '%g', int : '%d'}

   @property
   def format(self):
      return self._format

   @format.setter
   def format(self, fmt):
      rematch = RDB._fmtre.match(fmt)
      if rematch is None:
         raise RDBError('Invalid format statement!')
      typecode = rematch.groups()[0]
      if typecode == 's': self._format[str] = fmt
      if typecode in ('i', 'd'): self._format[int] = fmt
      if typecode in ('f', 'g'): self._format[float] = fmt

   @_ensurefile('w')
   def write_to_file(self, f):
      """ Writes the RDB table to a file """
      if not self._types:
         # Set our types
         self._set_types()
      print('\t'.join(self), file=f)
      print('\t'.join(self.formats), file=f)
      mylen = len(self[self.keys()[0]])
      for key in self:
         if len(self[key]) != mylen:
            raise RDBError('data length mismatch')
      for i in range(mylen):
         print('\t'.join(
            [self._format[self._types[j]] % (self[key][i]) for j, key in
                          enumerate(self)]),
            file=f
         )

   def _set_types(self):
      """ Determine the type of each of our data sets """
      global _HAS_NUMPY
      for key in self:
         if all([isinstance(val, int) for val in self[key]]):
            self._types.append(int)
         elif all([isinstance(val, float) for val in self[key]]):
            self._types.append(float)
         elif _HAS_NUMPY:
            if isinstance(self[key], np.ndarray):
               # Now see if it is any of the float types
               if str(self[key].dtype).startswith('float'):
                  self._types.append(float)
               elif str(self[key].dtype).startswith('int'):
                  self._types.append(int)
               else:
                  self._types.append(str)
            else:
               self._types.append(str)
         else:
            self._types.append(str)
         
   @classmethod
   def load_from_dict(cls, data):
      """ Load this from another dict class """
      inst = cls()
      for key in data: inst[key] = data[key]

   def check(self, reporter=stderr):
      """ Checks that everything is valid. Return 0 for OK; 1 for not """
      retval = 0
      same = all([len(self[key]) == len(self[self.keys()[0]]) for key in self])
      if not same:
         reporter.write('Length mismatch:\n\t')
         reporter.write('\n\t'.join(['%s: %d' % (key, len(self[key])) 
                                     for key in self]))
         reporter.write('\n')
         retval = 1
      return retval

   @classmethod
   @_ensurefile('r', clsmeth=True)
   def load_from_file(cls, f, use_numpy=True):
      """ Loads the data structure from an RDB file """
      global _HAS_NUMPY
      if use_numpy and not _HAS_NUMPY:
         warnings.warn('numpy could not be found. Using array.')
         use_numpy = False
      inst = RDB()
      for key in f.readline()[:-1].split('\t'):
         inst[key] = []
      fmtline = f.readline().strip()
      if fmtline:
         inst.formats = fmtline.split('\t')
         if len(inst.formats) != len(inst):
            raise RDBError('badly-formed formats line. wrong number of formats')
      line = f.readline()[:-1]
      # determine initial types
      data = line.split('\t')
      for i, d in enumerate(data):
         try:
            inst[inst.keys()[i]].append(int(d))
            inst._types.append(int)
         except ValueError:
            # Not integer, try float
            try:
               inst[inst.keys()[i]].append(float(d))
               inst._types.append(float)
            except ValueError:
               # Not integer, store as unicode string
               inst[inst.keys()[i]].append(str(d))
               inst._types.append(str)
      for line in f:
         data = line[:-1].split('\t')
         for i, d in enumerate(data):
            try:
               inst[inst.keys()[i]].append(inst._types[i](d))
            except ValueError:
               try:
                  if inst._types[i] is int:
                     inst[inst.keys()[i]].append(float(d))
                     inst._types[i] = float
                  else:
                     # Originally a float, now make it a string
                     inst[inst.keys()[i]].append(str(d))
                     inst._types[i] = str
               except ValueError:
                  # Only hit here if original type was int and current is string
                  inst[inst.keys()[i]].append(str(d))
                  inst._types[i] = str
      if use_numpy:
         # Convert data to numpy arrays
         for i, key in enumerate(inst):
            if inst._types[i] is int:
               inst[key] = np.asarray(inst[key], dtype=np.int)
            elif inst._types[i] is float:
               inst[key] = np.asarray(inst[key], dtype=np.float)
      else:
         # Convert data to Vector arrays
         for i, key in enumerate(inst):
            if inst._types[i] is int:
               inst[key] = Vector('i', inst[key])
            elif inst._types[i] is float:
               inst[key] = Vector('d', inst[key])
      return inst
