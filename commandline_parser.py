"""
This module contains classes and such that are responsible for parsing
command-line arguments for MMPBSA.py.  All of the files specified for use
in MMPBSA.py will be assigned as attributes to the returned class.

Last updated: 4/17/2011

                          GPL LICENSE INFO                             

  Copyright (C) 2009  Dwight McGee, Billy Miller III, and Jason Swails

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
   
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330,
  Boston, MA 02111-1307, USA.
"""

__version__ = 1.0

class CommandlineError(Exception):
   def __init__(self, msg='Command-line error'):
      self.msg=msg
   def __str__(self):
      return self.msg

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class Option:
   """ Base Option class that has normal attributes. It has a flag and has a
       single option that can follow it. It also stores the help string that
       will be printed upon request.
   """

   __name__ = 'Option'

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __init__(self, flag, variable, help='', default=None, optional=True):
      """ Sets up an option:
            flag: the CL flag that tags this option

            variable: the name of the variable this option is to be stored to
            help:     the message that's displayed when help is asked for
            default:  the default value this option will take on
      """
      self.flag = flag
      self.variable = variable
      self.help = help
      self.default = default
      self.optional = optional

      self.is_set = False # if this variable is set on CL or not

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Parse(self, option_list, arg_number):
      """ Parses this argument from the command-line, assuming arg_number
          is the argument that corresponds to the flag in question. It will
          add that attribute to option_list and return the location of the
          next command-line argument that should be parsed.

            option_list: where to add the value of this option
            arg_number:  which argument to start at
      """
      from sys import argv as args

      if len(args) < arg_number + 1: # We expect the item to follow
         raise CommandlineError('Error parsing %s. Not enough arguments!' % 
                                self.flag) 
      
      if args[arg_number] != self.flag:
         raise CommandlineError('expecting flag %s, but got %s' % 
                                self.flag, args[arg_number]) 
      
      self.is_set = True
      arg_number += 1

      # Make sure we didn't follow with another flag if we weren't expecting one
      if args[arg_number][0] == '-':
         raise CommandlineError('No variable following flag %s!' %
                                self.flag) 

      setattr(option_list, self.variable, args[arg_number])

      return arg_number + 1

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def UseDefault(self, option_list):
      """ This method is invoked to add the default option to option_list """
      if self.is_set:
         return 0 # no need for a value

      if self.default == None and not self.optional: 
         raise CommandlineError('%s has no default, and no value was given!' %
                                 self.flag) 

      setattr(option_list, self.variable, self.default)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class ArrayOption(Option):
   """ This is a command-line option that takes an arbitrary number of arguments
       following it and loads the result in an array
   """

   __name__ = 'ArrayOption'

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __init__(self, num_entries=-1, **kwargs):
      """ Adds a max number of entries to the base Option __init__ method """
      Option.__init__(self, **kwargs)  # call base class's __init__
      # if negative, take an arbitrary number of arguments
      self.num_entries = num_entries
      self.entries_left = num_entries 

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Parse(self, option_list, arg_number):
      """ Parses this argument from the command-line, assuming arg_number is the
          argument that corresponds to the flag in question. It will add that
          attribute to option_list and return the location of the next command-
          line argument that should be parsed.
      """
      from sys import argv as args

      if len(args) < arg_number + 1:
         raise CommandlineError('Error parsing %s. Not enough arguments!' % 
                                 self.flag) 

      if args[arg_number] != self.flag:
         raise CommandlineError('expecting flag %s, but got %s' % 
                                self.flag, args[arg_number]) 
      
      arg_number += 1
      arry = []

      # Parse through arguments until we've either reached our maximum number
      # of allowed arguments (self.entries_left), we hit a new flag, or we reach
      # the end of the argument list. If self.entries_left < 0, then this means
      # that we take an arbitrary number of arguments (it'll never reach 0). We
      # allow comma- and whitespace-delimited lists (any combo of them)

      while arg_number < len(args) and args[arg_number][0] != '-' and \
            self.entries_left != 0:
         if ',' in args[arg_number]:
            comma_args = args[arg_number].split(',')
            for i in range(len(comma_args)):
               comma_args[i] = comma_args[i].strip()
               if len(comma_args[i]) > 0:
                  arry.append(comma_args[i])
                  self.entries_left -= 1
               if self.entries_left == 0:
                  if i != len(comma_args) - 1:
                     raise CommandlineError('Too many options given for %s' %
                                             self.flag)
                  break
         else:
            arry.append(args[arg_number])
            self.entries_left -= 1

         arg_number += 1

      # Make sure we've found enough entries if we put a limit on it

      if self.entries_left > 0:
         raise CommandlineError('Expected %s more entries following %s!'%
                                self.entries_left, self.flag) 
      elif self.entries_left == 0:
         arg_number += 1 # add one here if we stopped before hitting next flag

      setattr(option_list, self.variable, arry)

      self.is_set = True

      return arg_number  # this is already at the next flag or at the end

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SetOption(Option):
   """ This class simply sets an attribute to true or false. It will always
       set the *opposite* of what the default is if it's found. 
   """

   __name__ = 'SetOption'

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Parse(self, option_list, arg_number):
      """ Parses this option from the command-line and returns the value of the
          next argument to parse.
      """
      from sys import argv as args

      if len(args) < arg_number:
         raise CommandlineError('Error parsing %s. Not enough arguments!' % 
                                self.flag) 

      self.is_set = True

      setattr(option_list, self.variable, not self.default)

      return arg_number + 1

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class OptionList:
   """ Just a container to hold the command-line options """
   pass

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class OptionParser:
   """ This is the parser class that organizes all of the Options and ultimately
       parses the command-line, filling the option list. It automatically 
       formats a help message based on descriptions you give. There are 3 kinds
       of options you can choose from: a standard option whose variable
       immediately follows the flag, an option that sets an attribute as true or
       false (so no argument follows the flag), or an option that takes an array
       of a specified or arbitrary number of arguments. You use it as follows:

       parser = OptionParser()

       parser.addOption('-f', my_var, help='This is my variable')
       parser.addOption('-list', my_list, help='Many options', num_entries=-1)
       parser.addOption('-noplot', plot, help='Disable plotting', default=True,
                        num_entries=0)
       parser.SetHelp(['-helpme', '-H', '-h', '--help'])

       parser.Parse()
   """

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __init__(self, stdout=None, stderr=None):
      """ Sets up the option parser """
      import sys
      from sys import argv as args
      from os.path import split

      self.options = {}
      self.option_keys_ordered = []
      self.help_flag = ['--help', '-h', '--h']
      self.program = split(args[0])[1]
      self.parsed_options = OptionList()
      if not stdout: self.stdout = sys.stdout
      else: self.stdout = stdout
      if not stderr: self.stderr = sys.stderr
      else: self.stderr = stderr

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def SetHelp(self, flags):
      """ Redefines the help flags """
      if type(flags).__name__ == 'list':
         self.help_flag = flags
      else:
         self.help_flag = [flags.strip()]

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
   
   def addHelp(self, flag):
      """ Adds a new flag to have the parser recognize a call for help """
      self.help_flag.append(flag)

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def addOption(self, flag, variable, help='', default=None, num_entries=1,
                 optional=True):
      """ Adds an option to the list of parsable options. It figures out
          what kind of option to assign based on how many entries are
          asked for. If > 1 or < 0, it adds an ArrayOption, if == 1, it adds
          an Option, if == 0, it adds a SetOption
      """

      # Guard against adding the same flag multiple times

      if flag in self.option_keys_ordered:
         raise CommandlineError('Flag %s added multiple times!' % flag) 
      
      self.option_keys_ordered.append(flag)

      if num_entries == 1:
         self.options[flag] = Option(flag, variable, help, default, optional)
      elif num_entries == 0:
         self.options[flag] = SetOption(flag, variable, help, default, optional)
      else:
         # Now it's an ArrayOption. Find out if we have to make our
         # default into an array or if it already is. If it already is,
         # pass it as-is. If it's not, typecast it to a string and
         # split it -- return that.
         if default != None:
            if type(default).__name__ == 'list':
               pass # nothing to do
            else:
               default = split(str(default))

         self.options[flag] = ArrayOption(num_entries=num_entries, flag=flag, 
                                variable=variable, help=help, default=default, 
                                optional=optional)

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Parse(self):
      """ Parses the command-line arguments, fills up, and returns the options,
          which is basically just a namespace where all of the CL options are
          dumped
      """
      from sys import argv as args

      i = 1
      j = 1 # This is a protection against an infinite loop!

      # the Parse routines of each Option class SHOULD increment i. We have
      # A j-counter running just in case that will jump out if it's cycled
      # through the command-line arguments more than 100x over, which should
      # never even come close to happening (in fact it should stop before 
      # len(args))

      if len(args) > 1 and args[1].lower() in self.help_flag:
         return self.print_help(self.stdout)

      if len(args) > 1 and args[1].lower() == '--version':
         return stdout.write(__version__)

      while i < len(args):
         
         try:
            i = self.options[args[i]].Parse(self.parsed_options, i)
         except KeyError:
            raise CommandlineError('Unknown flag %s!' % args[i]) 

         if j > len(args) * 100: # this is being very generous
            raise CommandlineError("I'm scared! I've been in the Parse loop " +
                  " long! I'm quitting now before I get in trouble.") 
      
         j += 1

      # Now go through and call each Option's "UseDefault" method. This should
      # simply return after doing nothing if it's already been set (thereby not
      # overwriting anything on the command-line with its default)

      for flag in self.option_keys_ordered:
         self.options[flag].UseDefault(self.parsed_options)
      
      return self.parsed_options

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def print_help(self, file=None):
      """ Prints the help message, and works on properly formatting it. 
          It writes it to a file as well, which will do its best to comply
          with both python2 and python3. Expect this to be stdout or stderr...
      """
      import sys
      if not file: file = self.stdout

      # Maximum Line length
      max_line_len = 72

      # Tab stop sizes for the [Options]
      leading_tabstop = 8
      second_tabstop = 25

      # Make strings of spaces 1 character shorter than the tabstops. This lets
      # us use these in "%s begin string" to get the desired spacing
      tab1 = '%%%ds' % (leading_tabstop - 1) % ' '
      tab2 = '%%%ds' % (second_tabstop - 1) % ' '

      # Number of spaces separating longest CL option and its variable name
      flag_option_sep_space = 3

      # Find the longest CL option that is not a SetOption. This includes the
      # help
      longest_noset_option = 0
      for key in self.option_keys_ordered:
         if self.options[key].__name__ != 'SetOption':
            longest_noset_option = max(longest_noset_option, len(key))
      
      # Now it's time to start writing. Write the header

      file.write('Usage: %s   [Options]\n\nOptions:\n%s ' %
                 (self.program, tab1))
      first_tag = True
      for msg in self.help_flag:
         if first_tag:
            file.write('%s' % msg)
            first_tag = False
         else:
            file.write(', %s' % msg)

      file.write('\n%s show this help message and exit\n\n' % tab2)

      for flag in self.option_keys_ordered:
         num_spaces_needed = flag_option_sep_space + longest_noset_option - \
                             len(flag)
         spaces = '%%%ds' % num_spaces_needed % ' '
         if self.options[flag].__name__ == 'SetOption':
            file.write('%s %s\n' % (tab1, flag))
         elif self.options[flag].__name__ == 'Option':
            file.write('%s %s%s%s\n' % (tab1, flag, spaces, 
                                        self.options[flag].variable))
         elif self.options[flag].__name__ == 'ArrayOption':
            if self.options[flag].num_entries == -1:
               file.write('%s %s%s%s1,%s2,...,%sN\n' % (tab1, flag, spaces,
                    self.options[flag].variable, self.options[flag].variable,
                    self.options[flag].variable))
            else:
               file.write('%s %s%s%s1,' % (tab1, flag, spaces, 
                          self.options[flag].variable))
               for i in range(1,self.options[flag].num_entries):
                  file.write('%s%d,' % (self.options[flag].variable, i+1))
               file.write('\n')

         file.write(self._format(flag, tab2, max_line_len) + '\n')

      sys.exit(1)

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def _format(self, flag, tab, max_line_len):
      """ Formats the description statement to make sure it fits in the line """

      help_msg = self.options[flag].help.strip()
      if len(help_msg) == 0:
         return ''

      if len(help_msg) + len(tab) < max_line_len:
         return '%s %s\n' % (tab, help_msg)

      # Now is the tricky part -- the description is sufficiently long that
      # we have to break it up into several lines.
      return_string = ''
      return_line = '%s ' % tab

      help_words = help_msg.split()
      for word in help_words:
         if len(return_line) + len(word) > max_line_len:
            return_string += return_line + '\n'
            return_line = '%s %s ' % (tab, word)
         else:
            return_line += '%s ' % word
      
      # Now make sure that the return_line is flushed to return_string
      if len(return_line.strip()) > 0:
         return_string += return_line + '\n'

      return return_string

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
