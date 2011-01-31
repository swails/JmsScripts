#!/usr/bin/env python

import sys, getatms

# list the variables to be used
residuelist = []
residuenamelist = []
firstatom = []
statelist = []
indexedarray = []
tempholder = []
holder = []

if len(sys.argv) < 3:
   print "Usage: VMD_pH.py <cpin> <cpout1> <cpout2> ... <cpoutN>"
   sys.exit()

# open the CPIN file for reading only
try:
   cpin = open(sys.argv[1],'r')
except IOError:
   print >> sys.stderr, "Error: cpin file, " + sys.argv[1] + " not found!"
   sys.exit()

#parse the cpin file.
for line in cpin:
  for line_index in range(len(line)):
    if line[line_index:line_index+7] == 'Residue':
      residuenamelist.append(line[line_index+9:line_index+12])
      num_toadd = ''
      line_index = line_index + 12
      while ( line_index < len(line) and line[line_index] != '\'' ):
        if line[line_index] != ' ' and line[line_index] != ':' and line[line_index] != '\n':
          num_toadd = num_toadd + line[line_index]
        line_index = line_index + 1
      residuelist.append(int(num_toadd))

    if line[line_index:line_index+10] == 'FIRST_ATOM': 
      line_index = line_index + 10
      num_toadd = ''
      while ( line_index < len(line) and line[line_index] != ',' ):
        if line[line_index] != ' ' and line[line_index] != '=':
          num_toadd = num_toadd + line[line_index]
        line_index = line_index + 1
      if num_toadd != '0':
        firstatom.append(int(num_toadd))
      
    # end for line_index		
# end parse of cpin file
cpin.close()
if len(residuenamelist) == 0: # check for CPIN data
   print "Error: invalid CPIN file! No residues found."
   sys.exit()

linenum = 0
numframes = 0
foundrst = 0
# loop through all of the cpout files
for x in range(2,len(sys.argv)):
  cpout = open(sys.argv[x],'r')
  isfirst = 1
  for line in cpout: # process each full record, since those accompany snapshots
    if line[0:10] == 'Solvent pH':
      numframes = numframes + 1 # this is a new frame
      foundrst = 1 # we have now found a full record
      continue   # done for this line, go to next step in loop
    if foundrst == 0: #if we have not yet found a full record
      continue # skip to next line
    elif foundrst < 4: # the first four lines have no residue data
      foundrst = foundrst + 1 # so move on to the next line
    elif foundrst < 4 + len(residuelist): # while we're still looking at residues
      if isfirst == 0: # The first full record corresponds to initial state, not to any frame\
        statelist.append(int(line[21]))  # only add to this if it's not initial state
      foundrst = foundrst + 1  # add one to this counter so we know when we've scanned everything
    else:
      foundrst = 0 # This is only hit if it's the first full record, so reset these values
      isfirst = 0 
numberframes = len(statelist)/len(residuenamelist)
for x in range(numberframes + 1):
  indexedarray.append(0) # initialize array for pairlist, but only the pointer part

indexpointer = numberframes + 1

for framenum in range(numberframes): # loop through each frame to find protons to remove
  indexedarray[framenum] = indexpointer # put the pointer into its location for each frame
  for resnum in range(len(residuelist)): # scan through all titratable residues to find hydrogens to disappear
    statelistindex = framenum * len(residuelist) + resnum 
    tempholder = tempholder + getatms.GetAtoms(residuenamelist[resnum], firstatom[resnum], statelist[statelistindex])
    holder = holder + tempholder
    indexpointer = indexpointer + len(tempholder)
    tempholder = []
   
indexedarray[numberframes] = indexpointer # set the last index pointer
tempholder = [] # empty out tempholder
indexedarray = indexedarray + holder # add the holder array to the indexed array
holder = []    # empty out holder

# now it's time to write the tcl script
outputfile = open('pHVMDscript.tcl', 'w')

for x in range(numberframes):
  line = "atomselect top \"index "
  for atom in range(indexedarray[x+1]-indexedarray[x]):
    line = line + str(indexedarray[indexedarray[x] + atom])
    line = line + " "
  line = line + "\" frame "
  line = line + str(x)
  outputfile.write(line)
  line = "atomselect"
  line = line + str(x)
  line = line + " moveby {99 99 99}" # move them far away to disappear
  outputfile.write("\n")
  outputfile.write(line)
  outputfile.write("\n")
