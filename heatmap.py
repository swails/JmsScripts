#!/usr/bin/env python

# A script to take a file containing two columns of data, and, from
# those two columns, construct and plot a 2-D histogram.

# General flow:
# 1. Read in the data file.
# 2. Use numpy's "histogram2d" method to bin the data.
# 3. Visualise the resulting object using matplotlib.

# Written by Benjamin P. Roberts, April 2011

import sys
import locale
from optparse import OptionParser
import math
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import matplotlib.colors as colours

parser = OptionParser(usage="usage: %prog [options] <datafile>")

# Set some default values
# Note: If any of these are changed, remember to also change the corresponding
# text in the help messages (see below).
parser.set_defaults(binsize=0.1, normalise=True, freeEnergy=False,
	unicode=False, raw=False, temperature=None, xub=None, xlb=None,
	xtickspacing=5.0, xlabel="X axis", yub=None, ylb=None,
	ytickspacing=5.0, ylabel="Y axis", resolution=80.0)

parser.add_option("-b", "--binsize", type="float", dest="binsize",
	metavar="NUM", help="Width of each bin ( > 0) [default: 0.1]")
parser.add_option("-e", "--free-energy", action="store_true", dest="freeEnergy",
	help="Generate a free-energy histogram instead of a frequency " + \
	"histogram. A temperature must also be specified with -t.")
parser.add_option("-t", "--temperature", type="float", dest="temperature",
	metavar="NUM", help="Temperature (K) to use when constructing a " + \
	"free energy histogram. Will be ignored unless -e is also specified.")
parser.add_option("--no-normalise", action="store_false", dest="normalise",
	help="Do not normalise values. Will be ignored if -e is specified.")
parser.add_option("--x-tick-spacing", type="float", dest="xtickspacing",
	metavar="NUM", help="Spacing between X axis major ticks (> 0) " + \
	"[default: 5.0]")
parser.add_option("--x-max", type="float", dest="xub", metavar="NUM",
	help="Maximum value to plot on X axis (if unspecified, an " + \
	"appropriate value will be computed from the data and the major " + \
	"tick spacing)")
parser.add_option("--x-min", type="float", dest="xlb", metavar="NUM",
	help="Minimum value to plot on X axis (if unspecified, an " + \
	"appropriate value will be computed from the data and the major " + \
	"tick spacing)")
parser.add_option("--x-label", type="string", dest="xlabel", metavar="STRING",
	help="X-axis label. For information regarding mathematical symbols " + \
	"and so forth, please see " + \
	"http://matplotlib.sourceforge.net/users/mathtext.html " + \
	"[default: the string \"X axis\"]")
parser.add_option("--y-tick-spacing", type="float", dest="ytickspacing",
	metavar="NUM", help="Spacing between Y axis major ticks (> 0) " + \
	"[default: 5.0]")
parser.add_option("--y-max", type="float", dest="yub", metavar="NUM",
	help="Maximum value to plot on Y axis (if unspecified, an " + \
	"appropriate value will be computed from the data and the major " + \
	"tick spacing)")
parser.add_option("--y-min", type="float", dest="ylb", metavar="NUM",
	help="Minimum value to plot on Y axis (if unspecified, an " + \
	"appropriate value will be computed from the data and the major " + \
	"tick spacing)")
parser.add_option("--y-label", type="string", dest="ylabel", metavar="STRING",
	help="Y-axis label. For information regarding mathematical symbols " + \
	"and so forth, please see " + \
	"http://matplotlib.sourceforge.net/users/mathtext.html " + \
	"[default: the string \"Y axis\"]")
parser.add_option("--resolution", type="float", dest="resolution",
	metavar="NUM",
	help="Resolution of final image (dpi) [default: 80.0]")

(options, args) = parser.parse_args()

if (len(args) != 1 \
or (options.freeEnergy == True and options.temperature == None)):
	parser.print_help()
	sys.exit(0)

if (options.binsize <= 0):
	print >> sys.stderr, "Bad value for bin size " + \
	"(should be greater than 0)"
	sys.exit(1)
if (options.xtickspacing <= 0):
	print >> sys.stderr, "Bad value for X tick spacing " + \
	"(should be greater than 0)"
	sys.exit(1)
if (options.ytickspacing <= 0):
	print >> sys.stderr, "Bad value for Y tick spacing " + \
	"(should be greater than 0)"
	sys.exit(1)
if (options.resolution <= 0):
	print >> sys.stderr, "Bad value for resolution " + \
	"(should be greater than 0)"
	sys.exit(1)

# Try to parse X-axis label.
defaultEncoding = sys.getdefaultencoding()
if (defaultEncoding == ""):
	print >> sys.stderr, "Unable to get default character encoding!"
	print >> sys.stderr, "This probably indicates a serious underlying " + \
	"problem."
	sys.exit(1)
fallbackEncodings = { "UTF-8" : "utf-8", "US-ASCII" : "utf-8" }

try:
	uxlabel = options.xlabel.decode(defaultEncoding)
	uylabel = options.ylabel.decode(defaultEncoding)
except UnicodeDecodeError:
	try:
		uxlabel = options.xlabel.decode \
			(fallbackEncodings[locale.nl_langinfo(locale.CODESET)])
		uylabel = options.ylabel.decode \
			(fallbackEncodings[locale.nl_langinfo(locale.CODESET)])
	except KeyError:
		print >> sys.stderr, "Don't know what Python codec to " + \
		"choose for character encoding \"" + \
		locale.nl_langinfo(locale.CODESET) + "\"."
		print >> sys.stderr, "Please change your character " + \
		"encoding to one of:"
		for fe in sorted(fallbackEncodings.keys()):
			print >> sys.stderr, fe
		print >> sys.stderr, "Alternatively, you may add a new " + \
		"mapping, or contact the developers."
		sys.exit(1)
	except UnicodeDecodeError:
		print >> sys.stderr, "Can't encode axis labels using " + \
		"either default or fallback character encodings."
		print >> sys.stderr, "Translation: You're up the creek. " + \
		"How did you enter these labels anyway?"
		sys.exit(1)

try:
	datafile = open(args[0], 'r')
except IOError as exc:
	print >> sys.stderr, exc
	sys.exit(1)

xvals = []
yvals = []
linecount = 0
for line in datafile:
	linecount = linecount + 1
	data = line.split()
	
	if (len(data) > 2):
		print "Warning: Too many entries on line "
		+ str(linecount)
		+ " (expected 2, found "
		+ str(len(data))
		+ "). Skipping."
		continue
	elif (len(data) < 2):
		print "Warning: Too few entries on line "
		+ str(linecount)
		+ " (expected 2, found "
		+ str(len(data))
		+ "). Skipping."
		continue
	else:
		try:
			xval = float(data[0])
		except ValueError as xve:
			print "Warning: On line "
			+ str(linecount)
			+ ", could not convert text \""
			+ data[0]
			+ "\" to a floating-point number. Skipping."
			continue
		
		try:
			yval = float(data[1])
		except ValueError as yve:
			print "Warning: On line "
			+ str(linecount)
			+ ", could not convert text \""
			+ data[1]
			+ "\" to a floating-point number. Skipping."
			continue
		
		xvals.append(xval)
		yvals.append(yval)

datafile.close()

xmax = max(xvals)
xmin = min(xvals)
ymax = max(yvals)
ymin = min(yvals)

if (options.xub == None):
	options.xub = float(int((xmax / options.xtickspacing) + 1.0)
		* int(options.xtickspacing))
if (options.xlb == None):
	options.xlb = float(int(xmin / options.xtickspacing)
		* int(options.xtickspacing))
if (options.yub == None):
	options.yub = float(int((ymax / options.ytickspacing) + 1.0)
		* int(options.ytickspacing))
if (options.ylb == None):
	options.ylb = float(int(ymin / options.ytickspacing)
		* int(options.ytickspacing))

if (options.binsize >= (options.xub-options.xlb) \
or options.binsize >= (options.yub-options.ylb)):
	print >> sys.stderr, "Error: Bins of " + \
	str(options.binsize) + \
	" are too large for this data set!"
	print >> sys.stderr, "That is, all visible data would be within " + \
	"a single bin."
	sys.exit(1)

maxticks = 20

# Set locations of lower and upper X ticks.
# The upper X tick is easy: The X upper bound, rounded down to the nearest
# whole multiple of the X tick spacing.
# The lower X tick is a little harder. If it's already at the X lower bound,
# just use that. Otherwise, if it should be above the X lower bound, round
# up to the nearest whole multiple of the X tick spacing.
# Finally, prepare an array containing the lower and upper X ticks, and
# everything in between.
xut = float(int(options.xub/options.xtickspacing)) * options.xtickspacing
if (options.xlb/options.xtickspacing > int(options.xlb/options.xtickspacing)):
	xlt = float(int((options.xlb/options.xtickspacing) + 1)) \
	* options.xtickspacing
else:
	xlt = options.xlb
xticks = []
xtick = xlt
while (xtick <= xut):
	xticks.append(xtick)
	xtick = xtick + options.xtickspacing
if (len(xticks) > maxticks):
	print >> sys.stderr, "Error: " + str(len(xticks)) + \
	" major ticks requested on X-axis (maximum allowed: "+ str(maxticks) + \
	")"
	print >> sys.stderr, "Try plotting a smaller window " + \
	"or using more widely spaced tick marks."
	sys.exit(1)
	
pyplot.xticks(xticks)

# Do the same for the lower and upper Y ticks.
yut = float(int(options.yub/options.ytickspacing)) * options.ytickspacing
if (options.ylb/options.ytickspacing > int(options.ylb/options.ytickspacing)):
	ylt = float(int((options.ylb/options.ytickspacing) + 1)) \
	* options.ytickspacing
else:
	ylt = options.ylb
yticks = []
ytick = ylt
while (ytick <= yut):
	yticks.append(ytick)
	ytick = ytick + options.ytickspacing
if (len(yticks) > maxticks):
	print >> sys.stderr, "Error: " + str(len(yticks)) + \
	" major ticks requested on Y-axis (maximum allowed: "+ str(maxticks) + \
	")"
	print >> sys.stderr, "Try plotting a smaller window " + \
	"or using more widely spaced tick marks."
	sys.exit(1)

pyplot.yticks(yticks)


# Compute the reciprocal binsize
recbin = 1/options.binsize

# Have to do a bit of legerdemain here. I originally tried modulo
# operations, which, alas, fall foul of floating-point imprecision.
epsilon = 1.0e-10
xbins = int((options.xub-options.xlb)*recbin)
ybins = int((options.yub-options.ylb)*recbin)
if (((options.xub-options.xlb)*recbin) - xbins > epsilon):
	print "Warning: Bins of "
	+ str(options.binsize)
	+ " do not fit evenly into a range of "
	+ str(options.xub-options.xlb)+"."
	print str(xbins)
	+ " bins will be used, for a bin size of "
	+ str((options.xub-options.xlb)/float(xbins))
	+ "."
if (((options.yub-options.ylb)*recbin) - ybins > epsilon):
	print "Warning: Bins of "
	+ str(options.binsize)
	+ " do not fit evenly into a range of "
	+ str(options.yub-options.ylb)+"."
	print str(ybins)
	+ " bins will be used, for a bin size of "
	+ str((options.yub-options.ylb)/float(ybins))
	+ "."

# Recall that, for reasons to do with compatibility, numpy breaks
# convention, putting the independent variable on the vertical axis
# and the dependent on the horizontal by default. So, we'll do a
# judicious swap.
if (options.freeEnergy == True):
	options.normalise = False
h2d, horizedges, vertedges = numpy.histogram2d(yvals, xvals,
	bins=[ybins,xbins], range=[[options.ylb,options.yub],
	[options.xlb,options.xub]], normed=options.normalise)

extent = [vertedges[0], vertedges[-1], horizedges[0], horizedges[-1]]

# Convert to free-energy plot if requested
# Note: The thermochemical calorie (1 cal = 4.184 J) has been used.
# For the gas constant, a value of 8.314472 J K^-1 mol^-1 has been used (NIST
# 2006).
if (options.freeEnergy == True):
	convfact = (8.314472/4.184) * options.temperature * -0.001
	maxfreq = 0
	for i in range(0,ybins):
		for j in range(0,xbins):
			if (h2d[i,j] > maxfreq):
				maxfreq = h2d[i,j]
	for i in range(0,ybins):
		for j in range(0,xbins):
			if (h2d[i,j] > 0):
				h2d[i,j] = math.log(h2d[i,j]/maxfreq) * convfact
			else:
				h2d[i,j] = None
	barLabel = r'$\Delta\mathit{G}$ (kcal mol$^{-1}$)'
	pyplot.imshow(h2d, extent=extent, interpolation='nearest',
		origin='lower')
# If we're not plotting the free-energy (and are just doing frequencies,
# normalised or not), we need to build a colour map that's appropriate for
# that purpose.
else:
	hmdict = { 'red': [	(0.0,1.0,1.0),
				(1.0e-100,0.0,0.0),
				(0.05,0.0,0.0),
				(0.25,1.0,1.0),
				(1.0,1.0,1.0)],
		'green': [	(0.0,1.0,1.0),
				(1.0e-100,1.0,1.0),
				(0.05,1.0,1.0),
				(0.25,1.0,1.0),
				(1.0,0.0,0.0)], 
		'blue': [	(0.0,1.0,1.0),
				(1.0e-100,1.0,1.0),
				(0.05,0.0,0.0),
				(1.0,0.0,0.0)] }
	hm = colours.LinearSegmentedColormap('heatmap', hmdict, 256)
	barLabel = "Frequency"
	pyplot.imshow(h2d, extent=extent, interpolation='nearest',
		origin='lower', cmap=hm)

pyplot.xlabel(uxlabel)
pyplot.ylabel(uylabel)
colourBar = pyplot.colorbar()
colourBar.set_label(barLabel)
defaultres = pyplot.gcf().get_dpi()
if (options.resolution != None):
	pyplot.gcf().set_dpi(options.resolution)

pyplot.gcf().savefig('heatmap.png', dpi=options.resolution, bbox_inches='tight')
