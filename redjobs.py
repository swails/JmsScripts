#!/usr/bin/env python

import sys
from pbsjob import PBSjob
from os import path, environ

if len(sys.argv) != 1 and len(sys.argv) != 2:
   print 'redjobs.py {job_name}'

red_version = "RED-vIII.4.pl"
gaussian_version = "g09"

redjob = PBSjob()

if len(sys.argv) == 2:
   redjob.name = sys.argv[1]

redjob.workdir = environ["PWD"]

if gaussian_version == "g03":
   redjob.addCommand("source $%sroot/em64t/%s/bsd/%s.profile" % (gaussian_version, gaussian_version, gaussian_version))
else:
   redjob.addCommand("source $%sroot/em64t/bsd/%s.profile" % (gaussian_version, gaussian_version))
redjob.addCommand("mkdir red_gauss_scr")
redjob.addCommand("export GAUSS_SCRDIR=red_gauss_scr")
redjob.addCommand("perl %s/%s > red.log" % (environ["RED"], red_version))
if len(sys.argv) != 2:
   redjob.randomName()

print redjob.preview()
redjob.submit()
