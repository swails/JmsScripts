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
from tkMessageBox import showerror, showwarning
from Tkinter import *
try:
   from scipy.stats import gaussian_kde
except ImportError:
   gaussian_kde = None

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

class LabeledEntry(Frame):
   """ This is a labeled entry widget """

   def __init__(self, master, text, *args, **kwargs):
      Frame.__init__(self, master)
      self.entry = Entry(self, *args, **kwargs)
      self.label = Label(self, text=text)
      self.entry.pack(fill=BOTH)
      self.label.pack(fill=BOTH)

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
      nexcl = self.graph_props.nexcl()
      if nexcl < 0 or nexcl > len(self.datasets[self.keylist[0]]):
         showerror('Bad Exclusions!',
                   'Number of excluded points must be greater than '
                   'zero and less than the number of data points!',
                   parent=self)
         return
      # Try to get our x data from the Time
      try:
         xdata = self.datasets['TIME(PS)'].copy()[nexcl:]
      except KeyError:
         xdata = np.arange(nexcl+1, len(self.datasets[self.keylist[0]])+1)
      
      if not self.graph_props.use_time():
         xdata = np.arange(nexcl+1, len(self.datasets[self.keylist[0]])+1)
      
      # Set the graph properties
      plt.xlabel(self.graph_props.xlabel())
      plt.ylabel(self.graph_props.ylabel())
      plt.title(self.graph_props.title())
      plt.grid(self.graph_props.gridlines())
      for i, a in enumerate(self.activelist):
         if not a.get(): continue
         # plot me
         props = self.graph_props.graph_options()
         if self.graph_props.legend():
            label = self.keylist[i]
         else:
            label = '_nolegend_'
         plt.plot(xdata, self.datasets[self.keylist[i]].copy()[nexcl:], 
                  label=label, **props)
      # Show the legend or not
      if self.graph_props.legend():
         plt.legend(loc=0)
      self.graph_props.reset_props()
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

      nexcl = self.graph_props.nexcl()
      if nexcl < 0 or nexcl > len(self.datasets[self.keylist[0]]):
         showerror('Bad Exclusions!',
                   'Number of excluded points must be greater than '
                   'zero and less than the number of data points!',
                   parent=self)
         return

      xdata = np.arange(nexcl+1, len(self.datasets[self.keylist[0]])+1)
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
         for i in range(nexcl, len(self.datasets[keys[1]])):
            f.write(' ' + ''.join(['%16.4f' % v[i] for v in actives]) + '\n')
      f.close()

class StatButton(_AnaButton):
   """ Prints average and standard deviation of each data set """

   def execute(self):

      report_str = '%20s%20s%20s\n' % ('Data Set', 'Avg.', 'Std. Dev.')
      report_str += '-'*60 + '\n'

      nexcl = self.graph_props.nexcl()
      if nexcl < 0 or nexcl > len(self.datasets[self.keylist[0]]):
         showerror('Bad Exclusions!',
                   'Number of excluded points must be greater than '
                   'zero and less than the number of data points!',
                   parent=self)
         return

      for i, val in enumerate(self.activelist):
         if not val.get(): continue
         dset = self.datasets[self.keylist[i]][nexcl:]
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
      
      nexcl = self.graph_props.nexcl()
      if nexcl < 0 or nexcl > len(self.datasets[self.keylist[0]]):
         showerror('Bad Exclusions!',
                   'Number of excluded points must be greater than '
                   'zero and less than the number of data points!',
                   parent=self)
         return
      
      if self.graph_props.use_kde() and gaussian_kde is None:
         showwarning('No scipy!', 'You must have scipy in order to use a '
                     'kernel density estimate to smooth the histograms!',
                     parent=self)
         self.graph_props.noscipy()
         
      # Set the graph properties
      plt.xlabel(self.graph_props.xlabel())
      plt.ylabel(self.graph_props.ylabel())
      plt.title(self.graph_props.title())
      plt.grid(self.graph_props.gridlines())
      for i, a in enumerate(self.activelist):
         if not a.get(): continue
         dset = self.datasets[self.keylist[i]].copy()[nexcl:]
         props = self.graph_props.graph_options()
         if self.graph_props.legend():
            label = self.keylist[i]
         else:
            label = '_nolegend_'
         # Plot either with a KDE or not
         if self.graph_props.use_kde():
            # Use a kernel density estimate for the histogramming to provide a
            # smooth curve
            kde = gaussian_kde(dset)
            # Use 200 points to characterize the surface. Go 1/100th of the range
            # out past either side of the max and min
            kmin = dset.min() - (dset.max() - dset.min()) / 100
            kmax = dset.max() + (dset.max() - dset.min()) / 100
            xdata = np.arange(kmin, kmax+0.000000001, (kmax-kmin)/200)
            ydata = np.asarray([kde.evaluate(x) for x in xdata])
            plt.plot(xdata, ydata, label=label, **props)
         else:
            # No KDE -- straight-out histogramming
            bw = self.graph_props.binwidth()
            if bw == 0:
               bw = 3.5 * dset.std() / len(dset) ** (1/3)
            if bw > 0:
               nbins = int((np.max(dset) - np.min(dset)) / bw)
            else:
               nbins = -int(bw)
            hist, bin_edges = np.histogram(dset, nbins,
                                           density=self.graph_props.normalize())
            plt.plot(bin_edges[:len(hist)], hist, label=label, **props)
      # Plot our function
      if self.graph_props.legend():
         plt.legend(loc=0)
      self.graph_props.reset_props()
      plt.show()

class AutoCorrButton(_AnaButton):
   """
   Does the normed autocorrelation function (with plotting) of the data set
   """

   def execute(self):
      if sum([v.get() for v in self.activelist]) == 0:
         showerror('No Data Sets!', 'No data sets chosen!', parent=self)
         return

      nexcl = self.graph_props.nexcl()
      if nexcl < 0 or nexcl > len(self.datasets[self.keylist[0]]):
         showerror('Bad Exclusions!',
                   'Number of excluded points must be greater than '
                   'zero and less than the number of data points!',
                   parent=self)
         return

      plt.clf()
      plt.cla()

      xdata = np.arange(nexcl+1, len(self.datasets[self.keylist[0]])+1)
      
      # Set the graph properties
      plt.xlabel(self.graph_props.xlabel())
      plt.ylabel(self.graph_props.ylabel())
      plt.title(self.graph_props.title())
      plt.grid(self.graph_props.gridlines())
      for i, a in enumerate(self.activelist):
         if not a.get(): continue
         # plot me
         props = self.graph_props.graph_options()
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
         plt.plot(xdata, acor, label=label, **props)
      # Deiconify the root
      if self.graph_props.legend():
         plt.legend(loc=0)
      self.graph_props.reset_props()
      plt.show()
      
if __name__ == '__main__':
   root = Tk()
   myvar = StringVar()
   ient = InputEntryWindow(root, myvar, 'Test!')
   root.mainloop()

   print 'myvar == ', myvar.get()
