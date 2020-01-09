#!/usr/bin/env python

from tkinter import *
from tkinter.ttk import Button, Scrollbar
from subprocess import Popen, PIPE

class TextWindow(Frame):
   """ A basic, scrollable text window that can either be editable or not """
   def __init__(self, master):
      """ Make a scrollable, resizeable text """
      super().__init__(master)
      # Add a horizontal and vertical scroller
      self.hscroller = Scrollbar(self, orient=HORIZONTAL)
      self.vscroller = Scrollbar(self, orient=VERTICAL)
      # Make the text box
      self.text = Text(self, width=90, state=DISABLED, height=15, wrap=NONE,
                       xscrollcommand=self.hscroller.set,
                       yscrollcommand=self.vscroller.set)
      # Pack everything in there, nice and tight. Let the text expand and the
      # scroll bars lengthen, but do not let the scroll bars thicken.
      self.text.grid(column=0, row=0, sticky=N+S+E+W)
      self.hscroller.grid(column=0, row=1, sticky=N+S+E+W)
      self.vscroller.grid(column=1, row=0, sticky=N+S+E+W)
      self.columnconfigure(0, weight=1)
      self.columnconfigure(1, weight=0)
      self.rowconfigure(0, weight=1)
      self.rowconfigure(1, weight=0)
      # Now make the scroll bars actually work
      self.hscroller.configure(command=self.text.xview)
      self.vscroller.configure(command=self.text.yview)

   def write(self, s):
      """
      Writes 's' to the window, such that it will emulate a file.  We have to
      change the state to ACTIVE in order to add text, but then change it back
      to the original state afterwards
      """
      self.text.configure(state=NORMAL)
      self.text.insert(END, s)
      self.text.configure(state=DISABLED)

   def clear(self, event=None):
      """ Clears all text from this window """
      self.text.config(state=NORMAL)
      self.text.delete('0.0', END)
      self.text.config(state=DISABLED)

class RefreshButton(Button):
   """ Allows refreshing of the fortune """
   def __init__(self, master, textwindow):
      super().__init__(text='Speak again!', command=self.update)
      self.textwindow = textwindow

   def update(self):
      """ Updates the text """
      self.textwindow.clear()
      self.textwindow.write('\n\n\n')
      self.textwindow.write(self.get_fortune())

   def get_fortune(self):
      """ Reads a fortune """
      try:
         process = Popen(['fortune'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
      except OSError:
         return '\tInstall "fortune" for Wanda to give your fortune'
      out, err = process.communicate(b'')
      if process.wait() != 0:
         raise RuntimeError('fortune failed!')
      return '\t' + '\n\t'.join(out.decode('utf-8').split('\n'))

def main():
    root = Tk()
    root.title('Wanda!')

    textbox = TextWindow(root)
    button = RefreshButton(root, textbox)

    textbox.pack(fill=BOTH, expand=1)
    button.pack(fill=BOTH, expand=0)
    button.update()

    root.mainloop()

if __name__ == '__main__':
    main()
