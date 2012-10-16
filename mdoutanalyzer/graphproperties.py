from Tkinter import BooleanVar, IntVar, StringVar, DoubleVar
class GraphProperties(object):
   """ Defines the graph properties """

   def __init__(self):
      self._lines = BooleanVar()
      self._points = BooleanVar()
      self._xlabel = StringVar()
      self._ylabel = StringVar()
      self._title = StringVar()
      self._linewidth = IntVar()
      self._legend = BooleanVar()
      self._bin_width = DoubleVar()
      self._normalize = BooleanVar()
      self._gridlines = BooleanVar()
      self._lines.set(True)
      self._points.set(False)
      self._legend.set(True)
      self._normalize.set(True)
      self._gridlines.set(True)
      self._bin_width.set(0.0)
      self._linewidth.set(1)
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
