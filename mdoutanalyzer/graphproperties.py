from Tkinter import BooleanVar
class GraphProperties(object):
   """ Defines the graph properties """

   def __init__(self):
      self._lines = BooleanVar()
      self._points = BooleanVar()
      self._lines.set(True)
      self._points.set(True)
   
   def lines(self):
      return self._lines.get()

   def points(self):
      return self._points.get()
