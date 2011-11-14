#!/usr/bin/env python
""" 
This program will combine a column from 2 different files, read in an operator,
and then output a column with the combination (col1 <operator> col2). It will
also read in an arbitrary delimiter
"""
from optparse import OptionParser
import os, sys

class StopLoop(Exception):
   """ Stops the loop when a line is broken """
   pass

def add(n1, n2): return n1 + n2
def sub(n1, n2): return n1 - n2
def mul(n1, n2): return n1 * n2
def div(n1, n2): return n1 / n2

parser = OptionParser()
parser.add_option('-f', '--file1', dest='file1', default=None,
              help='File containing first data column. Must specify.')
parser.add_option('-c', '--column1', dest='col1', default=None, type='int',
              help='Column with data from first file. Must specify.')
parser.add_option('-F', '--file2', dest='file2', default=None,
              help='File containing second data column. Must specify.')
parser.add_option('-C', '--column2', dest='col2', default=None, type='int',
              help='Column with data from second file. Default COL1')
parser.add_option('-d', '--delimiter', dest='delim', default=' ',
              help='Delimiter character to split fields. Default <space>')
parser.add_option('-o', '--operator', dest='op', default='+',
              help='Operator to combine the fields. Allowed + - / *, Default +')

(opt, arg) = parser.parse_args()

if arg:
   sys.stderr.write('Unrecognized arguments %s' % arg + os.linesep)
   parser.print_help()
   sys.exit(1)

if not (opt.file1 and opt.file2 and opt.col1):
   sys.stderr.write('Missing crucial arguments!' + os.linesep)
   parser.print_help()
   sys.exit(1)

if not opt.col2: opt.col2 = opt.col1

# Define the operation based on opt.op

if opt.op == '+': oper = add
elif opt.op == '-': oper = sub
elif opt.op == '/': oper = div
elif opt.op == '*': oper = mul
else:
   sys.stderr.write('Bad operator %s. Must be +,-,*, or /.'%opt.op + os.linesep)
   parser.print_help()

fl1 = open(opt.file1, 'r')
fl2 = open(opt.file2, 'r')

line1 = fl1.readline()
line2 = fl2.readline()

try:
   while line1 or line2:
      words1 = line1.split(opt.delim)
      d1 = False
      try: d1 = float(words1[opt.col1-1])
      except (ValueError, IndexError):
         while d1 == False:
            if not line1: raise StopLoop()
            line1 = fl1.readline()
            words1 = line1.split(opt.delim)
            try: d1 = float(words1[opt.col1-1])
            except (ValueError, IndexError): pass
      words2 = line2.split(opt.delim)
      d2 = False
      try: d2 = float(words2[opt.col2-1])
      except (ValueError, IndexError):
         while d2 == False:
            if not line2: raise StopLoop()
            line2 = fl2.readline()
            words2 = line2.split(opt.delim)
            try: d2 = float(words2[opt.col2-1])
            except (ValueError, IndexError): pass
      
      print '%s' % oper(d1, d2)
   
      line1 = fl1.readline()
      line2 = fl2.readline()
except StopLoop:
   pass
