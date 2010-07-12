#!/usr/bin/env python

###################################################################
#                                                                 #
# This program takes an amber input and uses ptraj to generate    #
# chi1-chi2 plots for specified residues. It formats the output   #
# for use in gnuplot.                                             #
#                                                                 #
# Last update: 10/21/2009   (by jms)                              # 
#  by Jason Swails                                                #
#                                                                 #
################################################################### 

#################   INSTRUCTIONS FOR USE   ########################
#
#  Follow the command line usage message. chi1chi2gen -help gives you
#  the usage message. -gnuplot sets the output for gnuplot (and will
#  also write the corresponding gnuplot scripts to load easily into
#  gnuplot. -origin sets the output for origin. -nobin prevents the
#  script from binning the chi1-chi2 angles.  You will get two sets of
#  data from this script: prefix.residue.dat which you can plot in 2-D
#  to get the individual points on the ramachandran torus projection,
#  and prefix.bin.dat which are the binned data files that plot prob-
#  abilities as a function of CHI1 (x) and CHI2 (y) angles. These are
#  the files plotted by the gnuplot scripts written by this script.
#
###################################################################




import sys, math, os, time  # modules that should exist on system
import utilities            # user-defined modules

tglobalstart = time.time()

chi_defines = {'ALA' : -1, 
               'ARG' : ['N','CA','CB','CG','CD'],
               'ASH' : ['N','CA','CB','CG','OD1'],
               'ASP' : ['N','CA','CB','CG','OD1'],
               'AS4' : ['N','CA','CB','CG','OD1'],
               'ASN' : ['N','CA','CB','CG','OD1'],
               'CYM' : -1,
               'CYS' : ['N','CA','CB','SG','HG'],
               'CYX' : -1,
               'GLH' : ['N','CA','CB','CG','CD'],
               'GLN' : ['N','CA','CB','CG','CD'],
               'GLU' : ['N','CA','CB','CG','CD'],
               'GL4' : ['N','CA','CB','CG','CD'],
               'GLY' : -1,
               'HID' : ['N','CA','CB','CG','ND1'],
               'HIE' : ['N','CA','CB','CG','ND1'],
               'HIP' : ['N','CA','CB','CG','ND1'],
               'ILE' : ['N','CA','CB','CG1','CD1'],
               'LEU' : ['N','CA','CB','CG','CD1'],
               'LYN' : ['N','CA','CB','CG','CD'],
               'LYS' : ['N','CA','CB','CG','CD'],
               'MET' : ['N','CA','CB','CG','SD'],
               'PHE' : ['N','CA','CB','CG','CD1'],
               'PRO' : ['N','CA','CB','CG','CD'],
               'SER' : ['N','CA','CB','OG','HG'],
               'THR' : ['N','CA','CB','OG1','HG1'],
               'TRP' : ['N','CA','CB','CG','CD1'],
               'TYR' : ['N','CA','CB','CG','CD1'],
               'VAL' : ['N','CA','CB','CG1','HG11'],
               'ACE' : -1,
               'NME' : -1,
               'NHE' : -1 }
               
if len(sys.argv) == 1 or '-help' in sys.argv[1] or sys.argv[1] == '-h':
   print 'Usage: chi1chi2gen.py -i inputfile -o output_prefix -p prmtop -y mdcrd1 mdcrd2 ... mdcrdN {-nobin} {-origin || -gnuplot}'
   print 'Default is "gnuplot" format'
   sys.exit()

if '-clea' in sys.argv[1]:
   os.system('rm -f _CHI1CHI2_*')
   sys.exit()

# default values
infile = 'inputfile'
prefix = 'CHI1CHI2_res'
topname = 'prmtop'
bins = [50,50]
mdcrds = []
residues = []
number_frames = 0  # useful so that we can normalize the ramachandran plot to get probabilities
binning = True
gnuformat = True

# get the command line arguments
for x in range(len(sys.argv)):
   if sys.argv[x].lower() == '-i':
      infile = sys.argv[x+1]
   elif sys.argv[x].lower() == '-o':
      prefix = sys.argv[x+1]
   elif sys.argv[x].lower() == '-p':
      topname = sys.argv[x+1]
   elif sys.argv[x].lower() == '-y':
      y = x + 1
      while y < len(sys.argv) and not '-' in sys.argv[y]:
         mdcrds.append(sys.argv[y])
         y = y + 1
   elif sys.argv[x].lower() == '-nobin':
      binning = False 
   elif sys.argv[x].lower() == '-origin':
      gnuformat = False

try:
   inputfile = open(infile,'r')
except IOError:
   print 'Error: Input file "' + infile + '" does not exist!'
   sys.exit()

inputlines = inputfile.readlines()
inputfile.close()

try:
   for x in range(len(inputlines)):
      if inputlines[x].strip().startswith('#'):
         continue
      if '=' in inputlines[x]:
         words = inputlines[x].split('=')
         if len(words) != 2:
            print 'Error: Invalid input file!'
            sys.exit()
         if 'residue' in words[0].lower():
            ress = words[1].split(',')
            for y in range(len(ress)):
               if '-' in ress[y]:
                  lims = ress[y].split('-')
                  if len(lims) != 2 or int(lims[0]) >= int(lims[1]):
                     print 'Error: Invalid input file!'
                     sys.exit()
                  for z in range(int(lims[0].strip()),int(lims[1].strip()) + 1):
                     residues.append(z)
               else:
                  residues.append(int(ress[y].strip()))
         elif 'bin' in words[0].lower():
            bins = words[1].lower().split('x')
            if len(bins) != 2:
               print 'Error: Invalid input file!'
               sys.exit()
            bins[0] = int(bins[0].strip())
            bins[1] = int(bins[1].strip())
except ValueError:
   print 'Error: Invalid input file!'
   sys.exit()

if utilities.fileexists(topname) == -1:
   sys.exit()

residue_names = utilities.getallresinfo(topname,'RESIDUE_LABEL')

if len(mdcrds) == 0:
   mdcrds.append('mdcrd')

residue_number = utilities.resnum(topname)
if residue_number == -1:
   print 'Error: Invalid topology file ' + topname + '!'
   sys.exit()

for x in range(len(residues)):
   if residues[x] < 1:
      print 'Error: You chose a nonsensical residue (residue 0 or less)'
      sys.exit()
   elif residues[x] > residue_number:
      print 'Error: The residue you chose is out of range!'
      sys.exit()

for x in range(len(mdcrds)):
   if utilities.fileexists(mdcrds[x]) == -1:
      sys.exit()

ptraj = utilities.which('ptraj')

if ptraj == 'none':
   print 'Error: ptraj needed for phipsigen.py!'
   sys.exit()
else:
   print 'ptraj Found! Using ' + ptraj

os.system('rm -f _CHI1CHI2_*')
ptrajin = open('_CHI1CHI2_ptraj.in','w')
for x in range(len(mdcrds)):
   ptrajin.write('trajin ' + mdcrds[x] + '\n')

ptrajin.write('\n')

for x in range(len(residues)):

   rsnu = residues[x]
   rsna = residue_names[rsnu-1]

   if chi_defines[rsna] == -1:
      print >> sys.stderr, 'Residue {0} ({1}) omitted. No chi1-chi2 defined for this residue'.format(rsnu,rsna)
   else:
      ptrajin.write('dihedral phires{0} :{0}@{1} :{0}@{2} :{0}@{3} :{0}@{4} out _CHI1CHI2_phires{0}\n'.format(rsnu,chi_defines[rsna][0],
                    chi_defines[rsna][1],chi_defines[rsna][2],chi_defines[rsna][3]))

      ptrajin.write('dihedral psires{0} :{0}@{1} :{0}@{2} :{0}@{3} :{0}@{4} out _CHI1CHI2_psires{0}\n'.format(rsnu,chi_defines[rsna][1],
                     chi_defines[rsna][2],chi_defines[rsna][3],chi_defines[rsna][4]))
ptrajin.close()

tptrajstart = time.time()
print 'Running ptraj to dump dihedral data...'
os.system('ptraj ' + topname + ' _CHI1CHI2_ptraj.in > _CHI1CHI2_ptraj.out 2>&1')
tptrajend = time.time()

print 'Creating combined chi1/chi2 files...'
for x in range(len(residues)):
   phifile = open('_CHI1CHI2_phires' + str(residues[x]),'r')
   psifile = open('_CHI1CHI2_psires' + str(residues[x]),'r')
   combfile = open(prefix + '.' + str(residues[x]) + '.dat', 'w')

   philine = phifile.readline()
   psiline = psifile.readline()

   while philine != '' and psiline != '':
      if x == 0:
         number_frames = number_frames + 1
      words = psiline.split()
      words2 = philine.split()
      if len(words) == 2 and len(words2) == 2:
         combfile.write(words2[1] + ' ' + words[1] + '\n')
      philine = phifile.readline()
      psiline = psifile.readline()

   phifile.close()
   psifile.close()
   combfile.close()

if binning:
   print 'Binning the results...'

   tbinstart = time.time()

   phipsibins = []
   phiinterval = 360.0/bins[0]
   psiinterval = 360.0/bins[1]
   phistart = 180.0/bins[0] - 180.0
   psistart = 180.0/bins[1] - 180.0
   # create a large 1-D array with every bin (phi1psi1 phi1psi2 ... phi1psiN phi2psi1 phi2psi2 ... ... phiNpsiN)
   for x in range(bins[0]):
      for y in range(bins[1]):
         phipsibins.append(0)
   for x in range(len(residues)):
      outputfile = open(prefix + '.' + str(residues[x]) + '.bin.dat','w')
      for y in range(len(phipsibins)):
         phipsibins[y] = 0
      datafile = open(prefix + '.' + str(residues[x]) + '.dat', 'r')
      datalines = datafile.readlines()
      datafile.close()
      # Fill the bins
      for y in range(len(datalines)):
         words = datalines[y].split()
         phival = float(words[0]) + 180.0
         psival = float(words[1]) + 180.0

         binnum = int(math.floor(phival/phiinterval) * bins[1] + math.floor(psival/psiinterval))
         phipsibins[binnum] = phipsibins[binnum] + 1

      phiprintval = phistart
      psiprintval = psistart
      for y in range(len(phipsibins)):
         if y != 0 and y % bins[1] == 0:
            phiprintval = phiprintval + phiinterval
            psiprintval = psistart
            if gnuformat:
               outputfile.write('\n')
         outputfile.write('{0} {1} {2}\n'.format(phiprintval, psiprintval, float(phipsibins[y])/float(number_frames)*360**2/(bins[0]*bins[1])))
         psiprintval = psiprintval + psiinterval

      outputfile.close()

   tbinend = time.time()

   print 'Writing gnuplot scripts...'

   if gnuformat:
      for x in range(len(residues)):
         gnuscript = open(prefix + '.' + str(residues[x]) + '.gnu','w')

         gnuscript.write('unset surface\nset contour base\nset cntrparam levels 20\nset cntrparam bspline\nset cntrparam order 7\n'
                         + 'set view 0,0\nunset ztics\n')
         gnuscript.write('splot \'' + prefix + '.' + str(residues[x]) + '.bin.dat\' w l')
         gnuscript.close()

tglobalend = time.time()

print 'Timings: '
print 'ptraj Dihedral Dump:' + '{0:.3f}'.format( -(tptrajstart-tptrajend)/60,3).rjust(6) + ' min.'
if binning:
   print 'Binning:            ' + '{0:.3f}'.format( -(tbinstart-tbinend)/60,3).rjust(6) + ' min.'
print 'Total time:         ' + '{0:.3f}'.format( -(tglobalstart-tglobalend)/60,3).rjust(6) + ' min.'
