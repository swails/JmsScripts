"""
This module contains code useful for reading, writing, and manipulating RDB
database tables (i.e., tab-delimited data fields with a header and format line)
"""
from __future__ import division, with_statement, print_function

from collections import OrderedDict
from functools import wraps
import re

# Python 3 support
try:
   unicode
except NameError:
   unicode = str

# Exception hierarchy
class RDBError(Exception): pass

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
#              if clsmeth:
                  return fcn(self, _f, *args, **kwargs)
#              return fcn(_f, *args, **kwargs)
      return newfcn
   return wrapper

class RDB(OrderedDict):
   """ RDB data """

   _fmtre = re.compile(r'%(?:[0-9-\.]*)([sgfdi])')
   def __init__(self, *args, **kwargs):
      super(RDB, self).__init__(*args, **kwargs)
      self.formats = []
      self._types = []
      self._format = {unicode : '%s', float : '%g', int : '%d'}

   @property
   def format(self):
      return self._format

   @format.setter
   def format(self, fmt):
      rematch = RDB._fmtre.match(fmt)
      if rematch is None:
         raise RDBError('Invalid format statement!')
      typecode = rematch.groups()[0]
      if typecode == 's': self._format[unicode] = fmt
      if typecode in ('i', 'd'): self._format[int] = fmt
      if typecode in ('f', 'g'): self._format[float] = fmt

   @_ensurefile('w')
   def write_to_file(self, f):
      """ Writes the RDB table to a file """
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

   @classmethod
   def load_from_dict(cls, data):
      """ Load this from another dict class """
      inst = cls()
      for key in data: inst[key] = data[key]

   @classmethod
   @_ensurefile('r', clsmeth=True)
   def load_from_file(cls, f):
      """ Loads the data structure from an RDB file """
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
               inst[inst.keys()[i]].append(unicode(d))
               inst._types.append(unicode)
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
                     inst[inst.keys()[i]].append(unicode(d))
                     inst._types[i] = unicode
               except ValueError:
                  # Only hit here if original type was int and current is string
                  inst[inst.keys()[i]].append(unicode(d))
                  inst._types[i] = unicode
      return inst
