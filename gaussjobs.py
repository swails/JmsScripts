#!/usr/bin/env python

import sys
from pbsjob import PBSjob
from os import path, environ

if len(sys.argv) != 2 and len(sys.argv) != 3:
   print 'gaussjobs.py <gaussian_input> {<job name>}'
   sys.exit()

gaussian_version = "g09"

try:
   l = open(sys.argv[1],'r')
except IOError:
   print >> sys.stderr, "Error: Gaussian input file %s does not exist!" % sys.argv[1]
   sys.exit()

if len(sys.argv) == 3:
   gaussjob.name = sys.argv[2]

l.close()

gaussjob = PBSjob()

gaussjob.workdir = environ["PWD"]

if gaussian_version == "g03":
   gaussjob.addCommand("source $%sroot/em64t/%s/bsd/%s.profile" % (gaussian_version, gaussian_version, gaussian_version))
else:
   gaussjob.addCommand("source $%sroot/em64t/bsd/%s.profile" % (gaussian_version, gaussian_version))
gaussjob.addCommand(gaussian_version + " " + path.split(sys.argv[1])[1])

if len(sys.argv) != 3:
   gaussjob.randomName()

print gaussjob.preview()
gaussjob.submit()
