from Tkinter import *
from mdoutanalyzer.widgets import (DataButton, GraphButton, SaveButton,
                                   StatButton, HistButton, AutoCorrButton,
                                   RunningAvgButton)
from mdoutanalyzer.graphproperties import GraphProperties
from mdoutanalyzer.menus import MainMenu

class MdoutAnalyzerApp(Frame):
   """ Main application """
   
   NCOLS  = 4
   NCOLS2 = 3

   def __init__(self, root, mdout):
      """ Set up the main app """
      self.mdout = mdout
      self.graph_props = GraphProperties()
      Frame.__init__(self, root)
      # Create the main file menu
      self.main_menu = MainMenu(root, self.graph_props, self.mdout)
      root.config(menu=self.main_menu)
      # Create a frame/canvas where we will put all of the Data buttons
      self.data_canv = Frame(self)
      self.data_active = [IntVar() for key in self.mdout.data]
      self.button_list = [DataButton(self.data_canv, self.data_active[i], key)
                          for i, key in enumerate(self.mdout.data)]
      # Label this
      Label(self.data_canv, text='Select Data to Analyze').grid(row=0, column=0,
            columnspan=self.NCOLS, sticky=N+S+E+W)
            
      # Grid the buttons in the canvas
      for i, b in enumerate(self.button_list):
         b.grid(row=1+i//self.NCOLS, column=i%self.NCOLS, sticky=N+S+E+W,
                padx=2, pady=2)
      self.data_canv.grid(column=0, row=0, sticky=N+S+E+W)

      # What do we want to do with the data?
      self.action_canv = Frame(self)
      self.acts = []
      self.acts.append(GraphButton(self.action_canv, self.mdout.data,
                                   self.data_active, self.graph_props,
                                   'Graph Them!'))
      self.acts.append(SaveButton(self.action_canv, self.mdout.data,
                                  self.data_active, self.graph_props,
                                  'Save to File'))
      self.acts.append(StatButton(self.action_canv, self.mdout.data,
                                  self.data_active, self.graph_props,
                                  'Show Statistics'))
      self.acts.append(HistButton(self.action_canv, self.mdout.data,
                                  self.data_active, self.graph_props,
                                  'Histogram Data'))
      self.acts.append(AutoCorrButton(self.action_canv, self.mdout.data,
                                      self.data_active, self.graph_props,
                                      'Autocorrelation'))
      self.acts.append(RunningAvgButton(self.action_canv, self.mdout.data,
                                        self.data_active, self.graph_props,
                                        'Running Average'))

      Label(self.action_canv, 
            text='What do you want to do with your data?').grid(
            row=0, column=0, columnspan=self.NCOLS2, sticky=N+S+E+W)

      # Grid the buttons in the canvas
      for i, b in enumerate(self.acts):
         b.grid(row=1+i//self.NCOLS2, column=i%self.NCOLS2, sticky=N+S+E+W,
                padx=2, pady=2)
      self.action_canv.grid(column=0, row=1, sticky=N+S+E+W)
