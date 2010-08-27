#!/usr/bin/env python

import signal
import sys
import os
from optparse import OptionParser
from mdin import mdin

files = {}

# General strategy for dealing with a situation where a later
# option forces a revision of an earlier one: Go back to the
# earlier option, which should be in its own function. At the
# end of executing that function, go back via a GOTO or Python's
# equivalent of the same (as messy as those are). To facilitate
# such GOTOs, a label identifying where the execution flow
# arrived from should be passed to the "earlier" function. If
# the label is blank, the function should not attempt to go
# anywhere, but should rather continue normal execution.

#======================================================================
def get_filenames():
   
   # parse_command_line
   #
   # Take an array of arguments (excepting the program name)
   # and use it to determine the names of two key entities:
   # the script file (*.bash) and the config file (*.mdin).
      
   # By default, an OptionParser will do the following:
   # If called with an invalid argument, print thus:
   #   usage: program [options]
   # If called with -h or --help, print the full help strings
   
   # Simulation input file is accessible from -m or --mdin=
   # Script name is accessible from -s or --script=
   parser = OptionParser()
   
   parser.add_option("-m", "--mdin", action="store", dest="mdin_name",
                     help="mdin file name")
   parser.add_option("-s", "--script", action="store", dest="scriptname",
                     help="script file name")
   
   (cl_options, cl_arguments) = parser.parse_args()
      
   if (not cl_options.scriptname):
      sys.stdout.write("Script file name: ")
      cl_options.scriptname = sys.stdin.readline().strip()
   
   if (not cl_options.mdin_name):
      sys.stdout.write("Amber mdin file name: ")
      cl_options.mdin_name = sys.stdin.readline().strip()
   
   return (cl_options.scriptname, cl_options.mdin_name)
#======================================================================



#======================================================================   
def interrupt_handler(signal, frame):
   
   # interrupt_handler
   #
   # Print out a message if the script has been interrupted (by
   # pressing Ctrl-C, for example) instead of finishing in the
   # usual manner.
   
   sys.stdout.write("\n")
   sys.stdout.write("\n")
   sys.stdout.write("           *********************\n")
   sys.stdout.write("           *                   *\n")
   sys.stdout.write("           *    Interrupted    *\n")
   sys.stdout.write("           *                   *\n")
   sys.stdout.write("           *********************\n")
   sys.stdout.write("\n")
   exit(0)
#======================================================================



#======================================================================
def close_files(exec_name):
   
   # Close the script and mdin files
   script.close()
   mdin.close()
   
   # Make the script executable (and, more generally, readable by everyone
   # and writable by its owner).
   os.chmod(exec_name, 0755)
#======================================================================
   


#======================================================================
def finished():
   
   # finished
   #
   # Print out a message confirming that the script has
   # run to completion, and exit.
   
   sys.stdout.write("\n")
   sys.stdout.write("           o----------------o\n")
   sys.stdout.write("           |                |\n")
   sys.stdout.write("           |    Finished    |\n")
   sys.stdout.write("           |                |\n")
   sys.stdout.write("           o----------------o\n")
   sys.stdout.write("\n")
   exit(0)
#======================================================================



#======================================================================
def get_title():
   
   # get_title
   # Get the simulation title for printing to the mdin file
   
   sys.stdout.write("Simulation title: ")
   title = sys.stdin.readline()
   mdin.write(title)
#======================================================================



#======================================================================
def query_yes_no(question, default="yes"):
   
   # query_yes_no
   # Ask a yes/no question via raw_input() and return their answer.
   #
   # Adapted from http://code.activestate.com/recipes/577058-query-yesno/
   # (BPR 16 April 2010)
   #
   # "question" is a string that is presented to the user.
   # "default" is the presumed answer if the user just hits <Enter>.
   # It must be "yes" (the default), "no" or None (meaning
   # an answer is required of the user).
   
   # The "answer" return value is one of "yes" or "no".
   valid = {"yes":True, "y":True, "ye":True,
            "no":False,   "n":False}
   if default == None:
      prompt = " [y/n] "
   elif default == "yes":
      prompt = " [Y/n] "
   elif default == "no":
      prompt = " [y/N] "
   else:
      raise ValueError("invalid default answer: '%s'" % default)
   
   while 1:
      sys.stdout.write(question + prompt)
      choice = raw_input().lower()
      if default is not None and choice == '':
         return valid[default]
      elif choice in valid.keys():
         return valid[choice]
      else:
         sys.stdout.write("Please respond with 'yes' or 'no' "\
                          "(or 'y' or 'n').\n")
#======================================================================



#======================================================================
def query_parallel():
   
   # query_parallel
   # Determine whether the simulation is to be run in parallel or not
   
   do_parallel = query_yes_no("Run sander in parallel?", None)
   
   if do_parallel:
      parallel = True
      setup_parallel()
   else:
      parallel = False
      files["sander_exe"] = "sander"
#======================================================================
   


#======================================================================
def setup_parallel():
   
   # setup_parallel
   # Set up the MPI executor
   
   files["sander_exe"] = "sander.MPI"
   
   sys.stdout.write("MPI executable name: ")
   files["mpiexec"] = sys.stdin.readline().rstrip('\r\n')
   
   sys.stdout.write("Number of CPUs to use: ")
   num_cpus = sys.stdin.readline().rstrip('\r\n')
#======================================================================
   




# Now, start the main execution

# First, we invoke the interrupt handler.
signal.signal(signal.SIGINT, interrupt_handler)


(scriptname, mdin_name) = get_filenames()

# Attempt to open the script file. If the open process generates
# an IO Error (for example, a file by that name already exists and
# is not writable by the current user, or the containing directory
# is not writable by the current user), then keep asking for an
# alternative until the system is happy.
while 1:
   try:
      script = open(scriptname, 'w')
      break
   except IOError:
      print "Could not open file", scriptname, "for writing."
      sys.stdout.write("Alternative script file name: ")
      scriptname = sys.stdin.readline().rstrip('\r\n')

# Print out the shebang.
script.write("#!/bin/sh\n\n")

# Do the same for the mdin file.
while 1:
   try:
      mdin = open(mdin_name, 'w')
      break
   except IOError:
      print "Could not open file", mdin_name, "for writing."
      sys.stdout.write("Alternative mdin file name: ")
      mdin_name = sys.stdin.readline().rstrip('\r\n')

query_parallel()

get_title()


mdin_contents = mdin("sander")



mdin_contents.write(mdin_name)


close_files(scriptname)

finished()
