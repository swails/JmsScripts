#!/usr/bin/env python

"""
This program will launch a minimization job to the local PBS server. The
convention followed is to strip the suffix from the topology file and use
that as the prefix for all other output files. This program will append
".min.xxxxx" for each file name, respectively. All trajectories are NetCDF
"""

from pbsjob import PBS_Script
import amber_simulations as am_sim
from optparse import OptionParser
import os, sys

class MinimizationError(Exception): pass

parser = OptionParser(usage='%prog [options] <prmtop> <inpcrd>')

parser.add_option('--maxcyc', dest='maxcyc', default=1000, type='int',
                  help='Number of minimization steps to run')
parser.add_option('--igb', dest='igb', default=5, type='int',
                  help='GB model to run for non-periodic systems. Must be ' +
                  '1, 2, 5, 7, or 8. Default 5')
parser.add_option('--restrain', dest='rst_wt', default=0, type='float',
                  help='Restraint weight to put on restraint mask. 0 means ' +
                  'no restraint. Default 0.')
parser.add_option('--restraint-mask', dest='rst_mask', default='@CA,C,O,N',
                  help='Restraint mask for restrained minimization. Default ' +
                  '"@CA,C,O,N"')
parser.add_option('--nproc', dest='nproc', default='1',
                  help='Number of processors to use. Fields are ,-delimited ' +
                  'and must fit into the processor-count format required by ' +
                  'the local PBS configuration')
parser.add_option('--walltime', dest='walltime', default='30:00',
                  help='How long to ask PBS for resources. Default 30:00')
parser.add_option('--name', dest='job_name', default=None,
                  help='Name to give to PBS job. None will give a random ' +
                  'Final Fantasy character name.')
parser.add_option('--delay', dest='jobid', default=None,
                  help='If present, will delay the present job until after ' +
                  'the give job ID is complete successfully')
parser.add_option('--force', dest='ask', default=True, action='store_false',
                  help='Do not ask confirmation before submitting.')
parser.add_option('--ask', dest='ask', action='store_true',
                  help='Ask confirmation before submitting. Overrides ' +
                  'previous --force. Default behavior.')
parser.add_option('--print-jobfile', dest='jobfile', default=None,
                  help='Print the PBS jobfile instead of submitting it ' +
                  'directly to the queue. Providing no name will bypass this ' +
                  'step.')
parser.add_option('--pbs-nproc', dest='pbs_nproc', default=None,
                  help='In the case that --nproc is not applicable to allowed' +
                  ' resource requests, this will be used for PBS instead.')
parser.add_option('--queue', dest='pbs_queue', default=None,
                  help='To override the default queue in ~/.pbsdefaults')
parser.add_option('--pbs-template', dest='pbs_template',
                  default=os.path.join(os.getenv('HOME'), '.pbsdefaults'),
                  help='PBS template file to use for job. Defaults to ' +
                  '~/.pbsdefaults')
parser.add_option('--no-pbs', dest='pbs', default=True, action='store_false',
                  help='Just run using os.system(), not via PBS')

(opt, args) = parser.parse_args()

if len(args) != 2:
   print >> sys.stderr, 'Bad command-line arguments!'
   parser.print_help()
   sys.exit(1)

amsys = am_sim.AmberSystem(args[0], args[1])

if not amsys.periodic() and not opt.igb in [1,2,5,7,8]:
   raise MinimizationError('Bad igb value (%d)' % opt.igb)

if opt.maxcyc < 0:
   raise MinimizationError('Bad maxcyc value (%d)' % opt.maxcyc)

min_input = am_sim.Minimization(amsys, num_steps=opt.maxcyc, igb=opt.igb,
                                restrained=bool(opt.rst_wt), rst_wt=opt.rst_wt,
                                rst_mask=opt.rst_mask)

# Get the prefix
prefix = os.path.splitext(args[0])[0]

# Print the mdin file
min_input.write_mdin('min.mdin')

# Get the MPI command (mpiexec -n $NPROC, for example)
mpi_cmd = am_sim.get_mpi_cmd()

# Determine if we're doing parallel simulations
nproc = opt.nproc.split(',')
nproc = [int(i) for i in nproc]
parallel = False
for n in nproc:
   if n > 1: parallel = True

# Set up the command string

if parallel: prog_str = '%s pmemd.MPI -O ' % mpi_cmd
else: prog_str = 'pmemd -O '

if opt.pbs:
   # Make the PBS_Script instance
   pbs_job = PBS_Script(template=opt.pbs_template)
   
   # Change the queue if desired
   pbs_job.queue = opt.pbs_queue
   
   # Set the PBS job name
   pbs_job.set_name(opt.job_name)
   
   # Determine processor count
   nproc = opt.nproc.split(',')
   nproc = [int(i) for i in nproc]
   
   if not opt.pbs_nproc: pbs_nproc = nproc[:]
   else:
      pbs_nproc = opt.pbs_nproc.split(',')
      pbs_nproc = [int(i) for i in pbs_nproc]
   
   pbs_job.set_proc_count(pbs_nproc)
   
   # Set walltime
   pbs_job.set_walltime(opt.walltime)
   
   cmd_str = ("%s -i min.mdin -p %s -c %s -r %s.min.rst7 -o %s.min.mdout -inf" +
              " %s.min.mdinfo") % (prog_str, args[0], args[1], prefix, prefix, 
              prefix)

   if opt.rst_wt:
      cmd_str += " -ref %s" % args[1]
   
   pbs_job.add_command(cmd_str)
   
   # If we just want to print the jobfile
   if opt.jobfile:
      pbs_job.print_submit(opt.jobfile)
   
   elif opt.ask:
      pbs_job.submit_ask(after_job=opt.jobid)
   
   else:
      pbs_job.submit(after_job=opt.jobid)
else:
   
   cmd_str = ("%s -i min.mdin -p %s -c %s -r %s.min.rst7 -o %s.min.mdout -inf" +
              " %s.min.mdinfo") % (prog_str, args[0], args[1], prefix, prefix, 
              prefix)

   if opt.rst_wt:
      cmd_str += " -ref %s" % args[1]

   os.system(cmd_str)
