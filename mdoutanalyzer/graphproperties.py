"""
Data structure controlling the appearance of the graph and a GUI control panel
exposing those options
"""

from Tkinter import *
from mdoutanalyzer.widgets import LabeledEntry

class GraphProperties(object):
   """ Defines the graph properties """

   def __init__(self):
      self._lines = BooleanVar()
      self._points = BooleanVar()
      self._xlabel = StringVar()
      self._ylabel = StringVar()
      self._title = StringVar()
      self._linewidth = DoubleVar()
      self._linestyle = StringVar()
      self._pointstyle = StringVar()
      self._legend = BooleanVar()
      self._bin_width = DoubleVar()
      self._normalize = BooleanVar()
      self._use_kde = BooleanVar()
      self._gridlines = BooleanVar()
      self._nexcl = IntVar()
      self._use_time = BooleanVar()
      # Set defaults
      self._lines.set(True)
      self._points.set(False)
      self._legend.set(True)
      self._normalize.set(True)
      self._use_kde.set(False)
      self._gridlines.set(True)
      self._use_time.set(True)
      self._bin_width.set(0.0)
      self._nexcl.set(0)
      self._linewidth.set(1.0)
      self._linestyle.set('-')
      self._pointstyle.set('o')
      self._colorschemes = ['k', 'b', 'r', 'g', 'c', 'm', 'y']
      self._used_colors = 0
   
   def lines(self):
      return self._lines.get()

   def points(self):
      return self._points.get()

   def xlabel(self):
      return self._xlabel.get()

   def ylabel(self):
      return self._ylabel.get()

   def title(self):
      return self._title.get()

   def legend(self):
      return self._legend.get()

   def get_next_color(self):
      i = self._used_colors
      self._used_colors += 1
      return self._colorschemes[i % len(self._colorschemes)]

   def reset_colors(self):
      self._used_colors = 0

   def linewidth(self):
      return self._linewidth.get()

   def binwidth(self):
      return self._bin_width.get()

   def normalize(self):
      return self._normalize.get()

   def gridlines(self):
      return bool(self._gridlines.get())

   def use_time(self):
      return self._use_time.get()

   def linestyle(self):
      if self.lines():
         return self._linestyle.get()
      else:
         return ''

   def use_kde(self):
      return self._use_kde.get()

   def pointstyle(self):
      if self.points():
         return self._pointstyle.get()
      else:
         return ''

   def nexcl(self):
      return self._nexcl.get()

   def graph_options(self):
      """ Returns a dict of the graph options """
      return {'linestyle' : self.linestyle(),
              'marker' : self.pointstyle(),
              'linewidth' : self.linewidth()
             }
   
   def noscipy(self):
      """ Call this if we have no scipy -- turns off scipy-only options """
      self._use_kde.set(False)


class GraphControlWindow(Toplevel):
   """ Frame controlling the graph options """

   def __init__(self, master, graph_props):
      Toplevel.__init__(self, master)
      self.resizable(False, False)
      self.graph_props = graph_props
      # Draw the checkbuttons
      self.lines = Checkbutton(self, text='Use Lines', onvalue=True,
                       offvalue=False, variable=self.graph_props._lines)
      self.points = Checkbutton(self, text='Use Points', onvalue=True,
                       offvalue=False, variable=self.graph_props._points)
      self.legend = Checkbutton(self, text='Show Legend', onvalue=True,
                       offvalue=False, variable=self.graph_props._legend)
      self.gridl = Checkbutton(self, text='Show Grid Lines', onvalue=True,
                       offvalue=False, variable=self.graph_props._gridlines)
      self.norm = Checkbutton(self, text='Normalize Histograms', onvalue=True,
                       offvalue=False, variable=self.graph_props._normalize)
      self.time = Checkbutton(self, text='Use Time as X-Coord', onvalue=True,
                       offvalue=False, variable=self.graph_props._use_time)
      self.kde = Checkbutton(self, text='Use Gaussian KDE For Histograms',
                       onvalue=True, offvalue=False,
                       variable=self.graph_props._use_kde)
      # Now grid the entries
      self.lines.grid(column=0, row=0, padx=2, pady=2, sticky=N+S+E+W)
      self.points.grid(column=1, row=0, padx=2, pady=2, sticky=N+S+E+W)
      self.legend.grid(column=0, row=1, padx=2, pady=2, sticky=N+S+E+W)
      self.gridl.grid(column=1, row=1, padx=2, pady=2, sticky=N+S+E+W)
      self.norm.grid(column=0, row=2, padx=2, pady=4, sticky=N+S+E+W)
      self.time.grid(column=1, row=2, padx=2, pady=4, sticky=N+S+E+W)
      self.kde.grid(column=1, row=3, padx=2, pady=4, columnspan=2, sticky=N+S)
      # Draw the entries
      self.ti = LabeledEntry(self, 'Plot Title', width=80,
                             textvariable=self.graph_props._title)
      self.xl = LabeledEntry(self, 'X-axis Label', width=80,
                             textvariable=self.graph_props._xlabel)
      self.yl = LabeledEntry(self, 'Y-axis Label', width=80,
                             textvariable=self.graph_props._ylabel)
      self.lw = LabeledEntry(self, 'Line width', width=40,
                             textvariable=self.graph_props._linewidth)
      self.ls = LabeledEntry(self, 'Line Style', width=40,
                             textvariable=self.graph_props._linestyle)
      self.ps = LabeledEntry(self, 'Point Style', width=40,
                             textvariable=self.graph_props._pointstyle)
      self.bins = LabeledEntry(self, 'Bin Width or -Number of Bins', width=40,
                               textvariable=self.graph_props._bin_width)
      self.nexcl = LabeledEntry(self,
                                'Points to skip from beginning of data sets',
                                width=40, textvariable=self.graph_props._nexcl)
      # Now grid the entry buttons
      self.ti.grid(column=0, row=4, columnspan=2, pady=4, sticky=N+S+E+W)
      self.xl.grid(column=0, row=5, columnspan=2, pady=4, sticky=N+S+E+W)
      self.yl.grid(column=0, row=6, columnspan=2, pady=4, sticky=N+S+E+W)
      self.lw.grid(column=0, row=7, pady=4, padx=2, sticky=N+S+E+W)
      self.ls.grid(column=1, row=7, pady=4, padx=2, sticky=N+S+E+W)
      self.ps.grid(column=0, row=8, pady=4, padx=2, sticky=N+S+E+W)
      self.bins.grid(column=1, row=8, pady=4, padx=2, sticky=N+S+E+W)
      self.nexcl.grid(column=0, row=9, columnspan=2, pady=4, sticky=N+S+E+W)
