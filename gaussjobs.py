#!/usr/bin/env python

import sys
from pbsjob import PBS_Script
from os import path, getenv
from optparse import OptionParser, OptionGroup

parser = OptionParser(usage="%prog [Options] <gaussian_input>")
parser.add_option('-n', '--name', dest='name', default=None,
                  help='Name given to PBS job')
parser.add_option('-v', '--gaussian-version', dest='version', default='g09',
                  help='Which Gaussian version to use. Default = g09. ' +
                  'Allowed g03 or g09')
parser.add_option('-w', '--walltime', dest='walltime', default='5:00:00',
                  help='Time for job to run in time format. Default is ' +
                  '5 hours ("5:00:00")')
parser.add_option('-t', '--template', dest='template', default=path.join(
                  getenv('HOME'), '.pbsdefaults'), help='PBS template file ' +
                  'to set up default PBS options. Default = ~/.pbsdefaults')
parser.add_option('-p', '--path-to-root', dest='path', default=None,
                  help='Path to Gaussian executable directory. Defaults ' +
                  'to $g09root/g09 or $g03root/g03, depending on version')
group = OptionGroup(parser, 'Job Submission Options',
                    'Details how you would like your job submitted (or just ' +
                    'printed to a job file)')
group.add_option('-f', '--force', dest='force', action='store_true', 
                  default=False, help='Do not ask before submitting job. ' +
                  'Default = False')
group.add_option('-a', '--ask', dest='force', action='store_false',
                help='Ask before submitting job. Overrides previous -f/--force')
group.add_option('-j', '--jobfile', dest='jobfile', metavar='FILE',
                 help='Write submission script to a job file instead of ' +
                 'just submitting it.')
parser.add_option_group(group)
(opt, args) = parser.parse_args()

if len(args) != 1:
   parser.print_help()
   sys.exit(1)

if not opt.version in ('g09', 'g03'):
   print >> sys.stderr, 'Error: Gaussian version must be g03 or g09!'
   sys.exit(1)

if not opt.path:
   opt.path = path.join(getenv('%sroot' % opt.version))

if not path.exists(args[0]):
   print >> sys.stderr, "Error: Gaussian input file %s doesn't exist!" % args[1]
   sys.exit()

# Create the gaussian PBSjob instance
gaussjob = PBS_Script(name=opt.name, template=opt.template)
gaussjob.set_proc_count([1])

# Source the resource file
if not path.exists(path.join(opt.path, 'bsd', '%s.profile' % opt.version)):
   print >> sys.stderr, 'Error: Missing %sprofile script! Looked in %s' % (
            opt.version, path.join(opt.path, 'bsd'))
   sys.exit(1)

gaussjob.add_command('source %s' % (path.join(opt.path, 'bsd', '%s.profile' % 
                     opt.version)))

gaussjob.add_command(opt.version + ' ' + args[0])

# Set a random name if no name had been set
if not opt.name: gaussjob.set_name()

# Do we ask permission before submitting the job or just submit it?
if opt.jobfile: gaussjob.print_submit(opt.jobfile)
elif not opt.force: gaussjob.submit_ask()
else: gaussjob.submit()
