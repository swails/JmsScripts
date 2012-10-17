#!/usr/bin/env python

from tkFileDialog import askopenfilenames
from Tkinter import Tk, BOTH
from mdoutanalyzer import __version__, __author__, __date__
from mdoutanalyzer.mdout import AmberMdout
from mdoutanalyzer.toplevel_app import MdoutAnalyzerApp
from optparse import OptionParser

verstring = """
   %%prog : An AMBER MD output file parser and graphing utility

                              Version %s
                             %s

   Written by %s
""" % (__version__, __date__, __author__)

parser = OptionParser(usage='%prog [mdout1] [mdout2] ... [mdoutN]',
                      version=verstring)

opt, arg = parser.parse_args()

root = Tk()
root.title('Mdout Analyzer')
if not arg:
   arg = askopenfilenames(title='Select Mdout File(s)', parent=root,
                          filetypes=[('Mdout Files', '*.mdout'),
                                     ('All Files', '*')])

if not arg:
   print ('No mdout files chosen!')

for f in arg:
   try:
      mdout += AmberMdout(f)
   except NameError:
      mdout = AmberMdout(f)

app = MdoutAnalyzerApp(root, mdout)
app.pack(fill=BOTH)
# Update idle tasks here to make sure the whole app is filled in before making
# our window non-resizable.  In some instances, this can chop off the second
# frame
app.update_idletasks()
# Now make our window non-resizable
root.resizable(False, False)
# Enter our mainloop
root.mainloop()
