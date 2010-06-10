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
#  Follow the command line usage message. phipsigen -help gives you
#  the usage message. -gnuplot sets the output for gnuplot (and will
#  also write the corresponding gnuplot scripts to load easily into
#  gnuplot. -origin sets the output for origin. -nobin prevents the
#  script from binning the phi-psi angles.  You will get two sets of
#  data from this script: prefix.residue.dat which you can plot in 2-D
#  to get the individual points on the ramachandran torus projection,
#  and prefix.bin.dat which are the binned data files that plot prob-
#  abilities as a function of PHI (x) and PSI (y) angles. These are
#  the files plotted by the gnuplot scripts written by this script.
#
###################################################################




import sys, math, os, time  # modules that should exist on system
import utilities            # user-defined modules

tglobalstart = time.time()
if len(sys.argv) == 1 or '-help' in sys.argv[1]:
   print 'Usage: chi1chi2gen.py -i inputfile -o output_prefix -p prmtop -y mdcrd1 mdcrd2 ... mdcrdN {-nobin} {-origin || -gnuplot}'
   print 'Default is "gnuplot" format'
   sys.exit()

if '-clea' in sys.argv[1]:
   os.system('rm -f _PHIPSI_*')
   sys.exit()

# default values
infile = 'inputfile'
prefix = 'PHIPSI_res'
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

if len(mdcrds) == 0:
   mdcrds.append('mdcrd')

residue_number = utilities.resnum(topname)
if residue_number == -1:
   print 'Error: Invalid topology file ' + topname + '!'
   sys.exit()

for x in range(len(residues)):
   if residues[x] == 2 or residues[x] == residue_number:
      print 'Error: You cannot get phi/psi dihedrals for either terminus!'
      sys.exit()
   elif residues[x] < 1:
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

os.system('rm -f _PHIPSI_*')
ptrajin = open('_PHIPSI_ptraj.in','w')
for x in range(len(mdcrds)):
   ptrajin.write('trajin ' + mdcrds[x] + '\n')

ptrajin.write('\n')

for x in range(len(residues)):

   prevres = str(residues[x] - 1)
   curres  = str(residues[x])
   nextres = str(residues[x] + 1)
   ptrajin.write('dihedral phires' + curres + ' :' + curres + '@N :' + curres + '@CA :' + \
                 curres + '@CB :' + curres + '@CG out _PHIPSI_phires' + curres + '\n')

   ptrajin.write('dihedral psires' + curres + ' :' + curres + '@CA :' + curres + '@CB :' + \
                 curres + '@CG :' + curres + '@OD1 out _PHIPSI_psires' + curres + '\n')
ptrajin.close()

tptrajstart = time.time()
print 'Running ptraj to dump dihedral data...'
os.system('ptraj ' + topname + ' _PHIPSI_ptraj.in > _PHIPSI_ptraj.out 2>&1')
tptrajend = time.time()

print 'Creating combined phi/psi files...'
for x in range(len(residues)):
   phifile = open('_PHIPSI_phires' + str(residues[x]),'r')
   psifile = open('_PHIPSI_psires' + str(residues[x]),'r')
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
         outputfile.write('{0} {1} {2}\n'.format(phiprintval, psiprintval, float(phipsibins[y])/float(number_frames)))
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
