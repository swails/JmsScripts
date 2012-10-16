from Tkinter import *
from mdoutanalyzer.menus import MainMenu
from mdoutanalyzer.graphproperties import GraphProperties

class MdoutAnalyzerApp(Frame):
   """ Main application """

   def __init__(self, root, mdout):
      """ Set up the main app """
      self.mdout = mdout
      self.graph_properties = GraphProperties
      Frame.__init__(self, root)
      # Create the main file menu
      self.main_menu = MainMenu(root, self.graph_properties, self.mdout)
      
