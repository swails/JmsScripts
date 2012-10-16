from Tkinter import *
from tkFileDialog import askopenfilenames
from mdout import AmberMdout
from mdoutanalyzer.widgets import InputEntryWindow

class FileMenu(Menu):
   """ The main file menu """
   def __init__(self, master, mdout):
      Menu.__init__(self, master, tearoff=0)
      self.mdout = mdout
      self.add_command(label='Add Mdout File(s)', command=self._add_mdout)
      self.add_separator()
      self.add_command(label='Quit', command=master.master.destroy)

   def _add_mdout(self):
      """ Get open filenames """
      fnames = askopenfilenames(title='Select Mdout File(s)', parent=self,
                                filetypes=[('Mdout Files', '*.mdout'),
                                           ('All Files', '*')])
      for f in fnames:
         self.mdout += AmberMdout(f)

class GraphPropMenu(Menu):
   """ Controls the graph properties """
   def __init__(self, master, graph_props):
      Menu.__init__(self, master, tearoff=0)
      self.graph_props = graph_props
      self.add_checkbutton(onvalue=True, offvalue=False, label='Grid Lines',
                           variable=self.graph_props._gridlines)
      self.add_checkbutton(onvalue=True, offvalue=False, label='Use Points',
                           variable=self.graph_props._points)
      self.add_checkbutton(onvalue=True, offvalue=False, label='Use Lines',
                           variable=self.graph_props._lines)
      self.add_checkbutton(onvalue=True, offvalue=False, label='Show Legend',
                           variable=self.graph_props._legend)
      self.add_checkbutton(onvalue=True, offvalue=False,
                           label='Normalize Histograms',
                           variable=self.graph_props._normalize)
      self.add_separator()
      self.add_command(command=lambda: InputEntryWindow(self,
                               self.graph_props._xlabel, 'X-axis Label'),
                       label='Set X-axis Label')
      self.add_command(command=lambda: InputEntryWindow(self,
                               self.graph_props._ylabel, 'Y-axis Label'),
                       label='Set Y-axis Label')
      self.add_command(command=lambda: InputEntryWindow(self,
                               self.graph_props._title, 'Graph Title'),
                       label='Graph Title')
      self.add_command(command=lambda: InputEntryWindow(self,
                               self.graph_props._bin_width,
                              'Histogram Bin Width (negative means # of bins)'),
                       label='Histogram Bins')

class HelpMenu(Menu):
   """ About and Help """
   def __init__(self, master):
      Menu.__init__(self, master, tearoff=0)
      # Holders
      self.add_command(label='About', command=lambda: 0)
      self.add_separator()
      self.add_command(label='Help', command=lambda: 0)

class MainMenu(Menu):
   """ The main file menu for the main App """

   def __init__(self, master, graph_props, mdout):
      """ Sets up the main menu """
      Menu.__init__(self, master, tearoff=0, relief=FLAT,
                    activeborderwidth=0)
      self.graph_props = graph_props
      self.mdout = mdout
      self.file_menu = FileMenu(self, self.mdout)
      self.graph_prop_menu = GraphPropMenu(self, self.graph_props)
      self.help = HelpMenu(self)
      # Add the menus
      self.add_cascade(label='File', underline=0, menu=self.file_menu)
      self.add_cascade(label='Graph Options', underline=0,
                       menu=self.graph_prop_menu)
      self.add_cascade(label='Help', underline=0, menu=self.help)

if __name__ == '__main__':
   from mdoutanalyzer.graphproperties import GraphProperties
   # Testing part

   root = Tk()
   root.config(menu=MainMenu(root, GraphProperties(), object()))
   root.mainloop()
