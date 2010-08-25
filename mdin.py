# This module will create a sander/pmemd input

class mdin:

   import cntrl

   cntrl_vars = cntrl.cntrl

   program = ''

   def __init__(self, program = 'sander'):
      if program
