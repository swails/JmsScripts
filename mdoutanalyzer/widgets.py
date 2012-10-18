"""
List of useful widgets we use here.
"""
from __future__ import division

from csv import writer
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mdoutanalyzer.windows import TextWindow
import numpy as np
from tkFileDialog import asksaveasfilename
from tkMessageBox import showerror
from Tkinter import *

class InputEntryWindow(Toplevel):
   """ Widget that takes a single text input """
   
   def __init__(self, master, var, label, nchar=50):
      Toplevel.__init__(self, master)
      self.entry = Entry(self, textvariable=var, width=nchar)
      self.label = Label(self, text=label)
      self.button = Button(self, text='OK', command=self.destroy)
      self.entry.pack()
      self.label.pack()
      self.button.pack(fill=X)
      self.resizable(False, False)
      self.grab_set()

class DataButton(Checkbutton):
   """ Button for each data set parsed from the mdout file """
   
   def __init__(self, master, var, label):
      Checkbutton.__init__(self, master, indicatoron=False, text=label,
                           variable=var, width=20, height=2)

class _AnaButton(Button):
   """ General base class for analysis buttons """

   def __init__(self, master, datasets, activelist, graph_props, text):
      Button.__init__(self, master, text=text, width=20,
                      height=2, command=self.execute)
      self.datasets = datasets
      self.graph_props = graph_props
      self.activelist = activelist
      self.keylist = self.datasets.keys()
   
class GraphButton(_AnaButton):
   """ This is the button to graph all of the selected data sets """

   def execute(self):
      """ Graphs the data sets """
      if sum([v.get() for v in self.activelist]) == 0:
         showerror('No Data Sets!', 'No data sets chosen!', parent=self)
         return

      plt.clf()
      plt.cla()
      # Try to get our x data from the Time
      try:
         xdata = self.datasets['TIME(PS)'].copy()
      except KeyError:
         xdata = np.arange(1, len(self.datasets[self.keylist[0]])+1)
      props = {'linestyle' : '',
               'marker' : '',
               'linewidth' : self.graph_props.linewidth()}
      if self.graph_props.lines():
         props['linestyle'] = '-'
      if self.graph_props.points():
         props['marker'] = 'o'
      
      # Set the graph properties
      plt.xlabel(self.graph_props.xlabel())
      plt.ylabel(self.graph_props.ylabel())
      plt.title(self.graph_props.title())
      plt.grid(self.graph_props.gridlines())
      for i, a in enumerate(self.activelist):
         if not a.get(): continue
         # plot me
         if self.graph_props.legend():
            label = self.keylist[i]
         else:
            label = '_nolegend_'
         plt.plot(xdata, self.datasets[self.keylist[i]].copy(), label=label,
                  color=self.graph_props.get_next_color(), **props)
      # Deiconify the root
      if self.graph_props.legend():
         plt.legend(loc=0)
      self.graph_props.reset_colors()
      plt.show()

class SaveButton(_AnaButton):
   """ For saving the data to a file """

   def execute(self):
      """ Graphs the data sets """
      fname = asksaveasfilename(parent=self, defaultextension='.dat',
                                filetypes=[('Data File', '*.dat'),
                                           ('CSV File', '*.csv'),
                                           ('All Files', '*')])
      # Bail out if we cancelled
      if not str(fname.strip()):
         return

      xdata = np.arange(1, len(self.datasets[self.keylist[0]])+1)
      actives, keys = [xdata], ['Frame']
      for i, val in enumerate(self.activelist):
         if not val.get(): continue
         keys.append(self.keylist[i])
         actives.append(self.datasets[self.keylist[i]])
      
      f = open(fname, 'w')
      # Detect csv or not
      if fname.endswith('.csv'):
         csvwriter = writer(f)
         # Header
         csvwriter.writerow(keys)
         for i in range(len(self.datasets[keys[1]])):
            csvwriter.writerow([v[i] for v in actives])
      else:
         f.write('#' + ''.join(['%16s' % n for n in keys]) + '\n')
         for i in range(len(self.datasets[keys[1]])):
            f.write(' ' + ''.join(['%16.4f' % v[i] for v in actives]) + '\n')
      f.close()

class StatButton(_AnaButton):
   """ Prints average and standard deviation of each data set """

   def execute(self):

      report_str = '%20s%20s%20s\n' % ('Data Set', 'Avg.', 'Std. Dev.')
      report_str += '-'*60 + '\n'
      for i, val in enumerate(self.activelist):
         if not val.get(): continue
         dset = self.datasets[self.keylist[i]]
         report_str += '%20s%20.4f%20.8f\n' % (self.keylist[i],
                       np.sum(dset) / len(dset), dset.std())
      
      window = TextWindow(self.master, width=62, height=i+5)
      window.write(report_str)
      window.grab_set()

class HistButton(_AnaButton):
   """ Histograms the data """

   def execute(self):
      """ Graphs the histograms of the data sets """
      if sum([v.get() for v in self.activelist]) == 0:
         showerror('No Data Sets!', 'No data sets chosen!', parent=self)
         return

      plt.clf()
      plt.cla()
      props = {'linestyle' : '',
               'marker' : '',
               'linewidth' : self.graph_props.linewidth()}
      if self.graph_props.lines():
         props['linestyle'] = '-'
      if self.graph_props.points():
         props['marker'] = 'o'
      
      # Set the graph properties
      plt.xlabel(self.graph_props.xlabel())
      plt.ylabel(self.graph_props.ylabel())
      plt.title(self.graph_props.title())
      plt.grid(self.graph_props.gridlines())
      for i, a in enumerate(self.activelist):
         if not a.get(): continue
         # plot me
         if self.graph_props.legend():
            label = self.keylist[i]
         else:
            label = '_nolegend_'
         dset = self.datasets[self.keylist[i]].copy()
         bw = self.graph_props.binwidth()
         if bw == 0:
            bw = 3.5 * dset.std() / len(dset) ** (1/3)
         if bw > 0:
            nbins = int((np.max(dset) - np.min(dset)) / bw)
         else:
            nbins = -int(bw)
         hist, bin_edges = np.histogram(dset, nbins,
                                        density=self.graph_props.normalize())
         plt.plot(bin_edges[:len(hist)], hist, label=label,
                  color=self.graph_props.get_next_color(), **props)
      # Deiconify the root
      if self.graph_props.legend():
         plt.legend(loc=0)
      self.graph_props.reset_colors()
      plt.show()

class AutoCorrButton(_AnaButton):
   """
   Does the normed autocorrelation function (with plotting) of the data set
   """

   def execute(self):
      if sum([v.get() for v in self.activelist]) == 0:
         showerror('No Data Sets!', 'No data sets chosen!', parent=self)
         return

      plt.clf()
      plt.cla()
      try:
         xdata = self.datasets['TIME(PS)'].copy()
      except KeyError:
         xdata = np.arange(1, len(self.datasets[self.keylist[0]])+1)
      props = {'linestyle' : '',
               'marker' : '',
               'linewidth' : self.graph_props.linewidth()}
      if self.graph_props.lines():
         props['linestyle'] = '-'
      if self.graph_props.points():
         props['marker'] = 'o'
      
      # Set the graph properties
      plt.xlabel(self.graph_props.xlabel())
      plt.ylabel(self.graph_props.ylabel())
      plt.title(self.graph_props.title())
      plt.grid(self.graph_props.gridlines())
      for i, a in enumerate(self.activelist):
         if not a.get(): continue
         # plot me
         if self.graph_props.legend():
            label = self.keylist[i]
         else:
            label = '_nolegend_'
         dset = self.datasets[self.keylist[i]].copy()
         dset -= dset.sum() / len(dset)
         dset /= dset.std()
         dset2 = dset.copy() / len(dset)
         acor = np.correlate(dset, dset2, 'full')
         acor = acor[len(acor)//2:]
         plt.plot(xdata, acor, label=label,
                  color=self.graph_props.get_next_color(), **props)
      # Deiconify the root
      if self.graph_props.legend():
         plt.legend(loc=0)
      self.graph_props.reset_colors()
      plt.show()
      
if __name__ == '__main__':
   root = Tk()
   myvar = StringVar()
   ient = InputEntryWindow(root, myvar, 'Test!')
   root.mainloop()

   print 'myvar == ', myvar.get()
