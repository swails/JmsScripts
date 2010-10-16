#!/usr/bin/env python

# This python utility will take a prmtop and some optional arguments
# and return a CPIN file to be used with sander for constant pH
# simulations.

import sys, string
try:
   import cpin_data
   from cpin_utilities import *
   from readparm import amberParm
except ImportError:
   print >> sys.stdout, 'Error: You must put cpin_data.py, cpin_utilities.py, and ' + \
            'readparm.py in a directory referenced by PYTHONPATH!'
   sys.exit()

############################# Begin user-defined modules #########################
def addOn(line, toadd):
   import sys

   if len(line) + len(toadd) > 80:
      print >> sys.stdout, line
      return ' ' + toadd
   else:
      return line + toadd
#############################  End  user-defined modules #########################

if len(sys.argv) < 3:
   printusage()

if sys.argv[1].startswith('-h') or sys.argv[1].startswith('--h'):
   printusage()

# list all of the residues that can currently be titrated:
titratable = cpin_data.TITRATABLE
exp_pkas   = cpin_data.EXP_PKAS

# set some default variables
solvent_ions = ['WAT','Na+','Br-','Cl-','Cs+','F-','I-','K+','Li+','Mg+','Rb+','CIO','IB','MG2']
titrate_residues = titratable
igb = ' '
ignore_warn = False
prmtop = 'prmtop'
prmtop_residues = []
prmtop_indices = []
residue_states = []
system_name = "Unknown"
maxpKa = 100
minpKa = -100
resnums = []
notresnums = []
resname_spec = False
nresname_spec = False
resnum_spec = False
nresnum_spec = False
titrated_residue_names = []
titrated_residue_nums = []


# read command-line arguments
try:
   for x in range(1,len(sys.argv)):
      if sys.argv[x] == '-p':
         prmtop = sys.argv[x+1]
      elif sys.argv[x] == '-igb':
         igb = sys.argv[x+1]
      elif sys.argv[x] == '-system':
         system_name = sys.argv[x+1]
      elif sys.argv[x] == '-minpKa':
         minpKa = sys.argv[x+1]
      elif sys.argv[x] == '-maxpKa':
         maxpKa = sys.argv[x+1]
      elif sys.argv[x] == '--ignore-warnings':
         ignore_warn = True

      # Get arbitrary num of resnames to add to the list
      elif sys.argv[x] == '-resname':           # resname flag
         if nresname_spec:
            print >> sys.stderr, 'Error: You cannot specify -notresname and -resname ' + \
                  'simultaneously!'
            sys.exit()
         resname_spec = True
         titrate_residues = ''
         x += 1
         while x < len(sys.argv) and not sys.argv[x].startswith("-"):
            dataholder = sys.argv[x].split(",")
            for y in range(len(dataholder)):
               if len(dataholder[y].strip()) != 0:
                  titrate_residues += dataholder[y].strip() + ' '
            x += 1
         # while loop may have been triggered by a flag; rewind to read flag
         x -= 1

      # Get arbitrary num of resnames to remove from the list
      elif sys.argv[x] == '-notresname':        # notresname flag
         if resname_spec:
            print >> sys.stderr, 'Error: You cannot specify -notresname and -resname ' + \
                  'simultaneously!'
            sys.exit()
         nresname_spec = True
         x += 1
         while x < len(sys.argv) and not sys.argv[x].startswith("-"):
            dataholder = sys.argv[x].split(",")
            for y in range(len(dataholder)):
               titrate_residues = " %s " % titrate_residues # add space in so we match only whole residue names
               if len(dataholder[y].strip()) != 0:
                  titrate_residues = string.replace(titrate_residues, " %s " % dataholder[y].strip(), ' ')
            x += 1
         # while loop may have been triggered by a flag; rewind to read flag
         x -= 1

      # Get arbitrary num of resnames to remove from the list
      elif sys.argv[x] == '-resnum':        # resnum flag
         if nresnum_spec:
            print >> sys.stderr, 'Error: You cannot specify -notresnum and -resnum ' + \
                  'simultaneously!'
            sys.exit()
         resnum_spec = True
         x += 1
         while x < len(sys.argv) and not sys.argv[x].startswith("-"):
            dataholder = sys.argv[x].split(",")
            for y in range(len(dataholder)):
               if len(dataholder[y].strip()) != 0:
                  resnums.append(dataholder[y].strip())
            x += 1
         # while loop may have been triggered by a flag; rewind to read flag
         x -= 1

      # Get arbitrary num of resnames to remove from the list
      elif sys.argv[x] == '-notresnum':        # notresnum flag
         if resnum_spec:
            print >> sys.stderr, 'Error: You cannot specify -notresnum and -resnum ' + \
                  'simultaneously!'
            sys.exit()
         nresnum_spec = True
         x += 1
         while x < len(sys.argv) and not sys.argv[x].startswith("-"):
            dataholder = sys.argv[x].split(",")
            for y in range(len(dataholder)):
               if len(dataholder[y].strip()) != 0:
                  notresnums.append(dataholder[y].strip())
            x += 1
         # while loop may have been triggered by a flag; rewind to read flag
         x -= 1

      # Get arbitrary num of resnames to remove from the list
      elif sys.argv[x] == '-states':        # states flag
         x += 1
         while x < len(sys.argv) and not sys.argv[x].startswith("-"):
            dataholder = sys.argv[x].split(",")
            for y in range(len(dataholder)):
               if len(dataholder[y].strip()) != 0:
                  residue_states.append(dataholder[y].strip())
            x += 1
         # while loop may have been triggered by a flag; rewind to read flag
         x -= 1
except IndexError:
   print >> sys.stderr, 'Error: Command line error!'
   printusage()

prmtop_object = amberParm(prmtop)
if not prmtop_object.exists:
   sys.exit()
if igb == ' ':
   if not ignore_warn:
      print >> sys.stderr, 'Warning: igb value not specified. igb = 5 used as default. ' + \
            'Be sure to use this value for your constant pH simulations!'
   igb = '5'

try:
   prmtop_residues = prmtop_object.parm_data["RESIDUE_LABEL"]
   prmtop_indices = prmtop_object.parm_data['RESIDUE_POINTER']
except KeyError:
   'Error: Could not find RESIDUE_LABEL or RESIDUE_POINTER in %s!' % prmtop
   sys.exit()

try:
   radius_set = prmtop_object.parm_data["RADIUS_SET"][0].strip()
except KeyError:
   radius_set = "NO RADIUS SET SPECIFIED"

if not ignore_warn:
   if radius_set != 'H(N)-modified Bondi radii (mbondi2)':
      print >> sys.stderr, 'Error: mbondi2 radius set should be used for constant pH simulations. All reference energies were calculated'
      print >> sys.stderr, '       using these radii! Set --ignore-warnings to ignore this warning.'
      sys.exit()

# Turn resnum and notresnum into integers. Turn igb into integer. Turn maxpKa and
# minpKa into floats.

try:
   if len(resnums) > 0:
      for x in range(len(resnums)):
         resnums[x] = int(resnums[x])
   elif len(notresnums) > 0:
      for x in range(len(notresnums)):
         notresnums[x] = int(notresnums[x])
   if len(residue_states) > 0:
      for x in range(len(residue_states)):
         residue_states[x] = int(residue_states[x])
   igb = int(igb)
   maxpKa = float(maxpKa)
   minpKa = float(minpKa)
except ValueError:
   print >> sys.stderr, 'Error: Make sure you have the right data types for the right arguments!'
   printusage()

# Filter residue names out of titrate_residues based on maxpKa and minpKa

titrate_name_holder = titratable.split()
titrate_pka_holder = exp_pkas.split()

for x in range(len(titrate_name_holder)):
   titrate_name_holder[x] = titrate_name_holder[x].strip()
   titrate_pka_holder[x] = float(titrate_pka_holder[x].strip())

   # Now tack a space onto the beginning and end of titrate_residues so we can match whole residue names
   titrate_residues = " %s " % titrate_residues.strip()

   if titrate_pka_holder[x] < minpKa or titrate_pka_holder[x] > maxpKa:
      if " %s " % titrate_name_holder[x] in titrate_residues:
         titrate_residues = string.replace(titrate_residues, titrate_name_holder[x], '').strip()


# Now tack a space onto the beginning and end of titrate_residues so we can match whole residue names
titrate_residues = " %s " % titrate_residues.strip()

# Start compiling the list of titrated_residues
# First start by only pulling residue numbers if specified. If not, look for all
# matching residues in the prmtop. Then filter out the "notresnums". 

# add in all specified residues

if len(resnums) != 0:
   for x in range(len(resnums)):
      if resnums[x] > len(prmtop_residues):
         print >> sys.stderr, 'Error: You specified a residue greater than the number of '
         print >> sys.stderr, '       residues in the topology file ' + prmtop
         sys.exit()

      if not (" %s " % prmtop_residues[resnums[x]-1]) in titrate_residues:
         print >> sys.stderr, 'Error: Residue ' + str(resnums[x]) + ' is a ' + prmtop_residues[resnums[x]-1] + \
               '. This is not in the list of titratable residues!'
         sys.exit()
      else:
         titrated_residue_nums.append(resnums[x])

# otherwise, add in all matching residues

else:
   for x in range(len(prmtop_residues)):
      if " %s " % prmtop_residues[x] in titrate_residues:
         titrated_residue_nums.append(x+1)

# Check to make sure that at least one residue is going to be titrated

if len(titrated_residue_nums) == 0:
   print >> sys.stderr, 'Error: No titratable residues conform to your criteria!'
   sys.exit()

# Now filter out the notresnums

for x in range(len(notresnums)):
   for y in range(len(titrated_residue_nums)):
      if notresnums[x] == titrated_residue_nums[y]:
         titrated_residue_nums.pop(y)
         break

# We should now have all of the residue numbers that we want to titrate

# Get the names of the residues:

for x in range(len(titrated_residue_nums)):
   if len(titrated_residue_names) == 0:
      titrated_residue_names.append(prmtop_residues[titrated_residue_nums[x]-1])

   else:
      present = False
      for y in range(len(titrated_residue_names)):
         if titrated_residue_names[y] == prmtop_residues[titrated_residue_nums[x]-1]:
            present = True

      if not present:
         titrated_residue_names.append(prmtop_residues[titrated_residue_nums[x]-1])

# check that if -states was specified, that as many states were specified as residues
# are being titrated. Otherwise, fill it up with zeroes

if len(residue_states) != 0 and len(residue_states) != len(titrated_residue_nums):
   print >> sys.stderr, 'Error: You specified ' + str(len(residue_states)) + ' states ' + \
            'for ' + str(len(titrated_residue_nums)) + ' titrated residues!'
   sys.exit()

elif len(residue_states) == 0:
   for x in range(len(titrated_residue_nums)):
      residue_states.append(0)

# now it's time to write the cnstph namelist in the cpin file. This will be written
# to stdout the same way the last one did.

print >> sys.stdout, "&CNSTPH"


# load all of the necessary data in the same order as it is in titrated_residue_names
resdata = []

for x in range(len(titrated_residue_names)):
   resdata.append(cpin_data.getData(titrated_residue_names[x], igb))


# First print out the charge data

line = " CHRGDAT="

for x in range(len(resdata)):
   curdata = resdata[x]

   if curdata == -1:
      print >> sys.stderr, 'Error: Couldn\'t get data from cpin_data.py!'
      sys.exit()

   for y in range(len(curdata)):
      for z in range(2,len(curdata[y])):
         toadd = str(curdata[y][z]) + ','
         line = addOn(line, toadd)

print >> sys.stdout, line

# Next print out protcnt

line = ' PROTCNT='

for x in range(len(resdata)):
   curdata = resdata[x]

   for y in range(len(curdata)):
      toadd = str(curdata[y][1]) + ','
      line = addOn(line, toadd)

print >> sys.stdout, line

# Next print out resname

line = " RESNAME='System: Unknown',"

for x in range(len(titrated_residue_nums)):
   toadd = "'Residue: "
   toadd += prmtop_residues[titrated_residue_nums[x]-1] + ' ' + \
            str(titrated_residue_nums[x]) + "',"
   line = addOn(line, toadd)

print >> sys.stdout, line

# Next print out resstate

line = ' RESSTATE='

for x in range(len(titrated_residue_nums)):
   for y in range(len(titrated_residue_names)):
      if prmtop_residues[titrated_residue_nums[x]-1] == titrated_residue_names[y]:
         if residue_states[x] < 0 or residue_states[x] >= len(resdata[y]):
            print >> sys.stderr, 'Error: Residue states must range from 0 to num_states - 1 ' + \
                     'for each residue!'
            sys.exit()
   toadd = str(residue_states[x]) + ','
   line = addOn(line, toadd)

print >> sys.stdout, line

# Next print out stateinf

line = ' '

for x in range(len(titrated_residue_nums)):
   # FIRST_ATOM
   holder = str(prmtop_indices[titrated_residue_nums[x]-1])
   toadd = 'STATEINF(' + str(x) + ')%FIRST_ATOM=' + holder + ', '
   line = addOn(line, toadd)

   # FIRST_CHARGE and FIRST_STATE and NUM_ATOMS and NUM_STATES
   startchrg = 0
   startst = 0
   for y in range(len(titrated_residue_names)):
      if prmtop_residues[titrated_residue_nums[x]-1] == titrated_residue_names[y]:
         natm = len(resdata[y][0]) - 2
         nst = len(resdata[y])
         break
      else:
         startchrg += (len(resdata[y][0]) - 2) * (len(resdata[y]))
         startst += len(resdata[y])
   toadd = 'STATEINF(' + str(x) + ')%FIRST_CHARGE=' + str(startchrg) + ', '
   line = addOn(line, toadd)
   toadd = 'STATEINF(' + str(x) + ')%FIRST_STATE=' + str(startst) + ', '
   line = addOn(line, toadd)
   toadd = 'STATEINF(' + str(x) + ')%NUM_ATOMS=' + str(natm) + ', '
   line = addOn(line, toadd)
   toadd = 'STATEINF(' + str(x) + ')%NUM_STATES=' + str(nst) + ', '
   line = addOn(line, toadd)
   

print >> sys.stdout, line

# Next print out statene

line = ' STATENE='

for x in range(len(resdata)):
   for y in range(len(resdata[x])):
      toadd = str(resdata[x][y][0]) + ','
      line = addOn(line, toadd)

print >> sys.stdout, line

# Next print out TRESCNT

line = ' TRESCNT=' + str(len(titrated_residue_nums)) + ','

# Determine if we need to print out information for hybrid implicit/explicit CpHMD
for i in range(len(prmtop_residues)):
   if prmtop_residues[i] in solvent_ions:
      toadd = 'CPHFIRST_SOL=%s, CPH_IGB=%s, ' % (prmtop_object.parm_data['RESIDUE_POINTER'][i],igb)
      if not ignore_warnings:
         for j in range(i,len(prmtop_residues)):
            if not prmtop_residues[j] in solvent_ions:
               print >> sys.stderr, 'Warning: Residue (%s) found where water or ion was expected!' % prmtop_residues[j]
               sys.exit()
      line = addOn(line, toadd)
      break

print >> sys.stdout, line

# End the namelist

print >> sys.stdout, '/'

print >> sys.stderr, 'CPin created!'
