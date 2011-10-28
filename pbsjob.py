#!/usr/bin/env python
from os import getenv, path, linesep
from sys import stderr, stdout, stdin, exit
from optparse import OptionParser
from subprocess import Popen, PIPE
import re

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Make a dictionary of Final Fantasy names
names = {"FF2":["Frionel", "Maria", "Guy", "Minh", "Gordon", "Layla", "Richard",
                "Lionheart"],
         "FF3":["Sara", "Cid", "Desh", "Elia", "Allus", "Dorga", "Unne"],
         "FF4":["Cecil", "Kain", "Rosa", "Cid", "Rydia", "Tellah", "Edward",
                "Yang", "Palom", "Porom", "Edge", "Fusoya", "Golbez"],
         "FF5":["Bartz", "Reina", "Galuf", "Faris", "Krile", "Ex-death",
                "Gilgamesh", "Cid"],
         "FF6":["Terra", "Locke", "Edgar", "Sabin", "Shadow", "Celes", "Cyan",
                "Gau", "Setzer", "Mog", "Umaro", "Gogo", "Leo", "Kefka",
                "Strago", "Relm", "Kefka"],
         "FF7":["Cloud", "Barret", "Tifa", "Aeris", "RedXIII", "Yuffie",
                "CaitSith", "Vincent", "Cid", "Sephiroth"],
         "FF8":["Squall", "Rinoa", "Seifer", "Quistis", "Zell", "Irvine",
                "Selphie"],
         "FF9":["Zidane", "Vivi", "Garnet", "Adelbert", "Freya", "Eiko",
                "Quina", "Amarant", "Beatrix", "Queen", "Kuja", "Tantalus"],
         "FF10":["Tidus", "Yuna", "Auron", "Rikku", "Wakka", "Lulu", "Kimahri",
                 "Maester Seymour"],
         "FF102":["Yuna", "Rikku", "Paine", "Nooj", "Lenne", "ShuyinYuna",
                  "Rikku", "Paine", "Nooj", "Lenne", "Shuyin"],
         "FF12":["Vaan", "Ashe", "Basch", "Penelo", "Balthier", "Fran"],
         "FF13":["Lightning", "Serah", "Hope", "Snow", "Oerba", "Sazh", "Nora"]
        }

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TemplateFile(object):
   """ Template file """
   def __init__(self, my_file):
      """ Make sure that my_file has a write() built-in method """
      if not hasattr(my_file, 'write'):
         raise TypeError('TemplateFile takes a file object!')
      if type(my_file.write).__name__ != 'builtin_function_or_method':
         raise TypeError("file's write() method must be builtin function!")
      self.my_file = my_file
   
   def write(self, my_str):
      """ Write the string to the line with the linesep at the end """
      self.my_file.write(my_str + linesep)
   
   def __getattr__(self, attr):
      return getattr(self.myfile, attr)

   def close(self):
      return self.my_file.close()

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def pbsline(option, value):
   if not value: return ''
   return "#PBS %s %s" % (option, value) + linesep

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PBSClassError(Exception):
   """ Raised if there's a problem in the class """
   pass

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PBSDependError(Exception):
   """ Raised if there's a problem with a job dependency """
   pass

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PBSMissingError(Exception):
   """ Raised if a program is missing """
   pass

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PBSInputError(Exception):
   """ Raised if a user tries to do something not allowed """
   pass

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class QsubError(Exception):
   """ Raised if there's a problem submitting the job """
   pass

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PBS_Script(object):
   " This class contains all of the data for dealing with PBS job scripts "

   # ==============================

   def __init__(self, name=None, commands='', template=None):
      """ Loads the main attributes of the PBS job """
      self.account = ''
      self.date_time = ''
      self.interval = ''
      self.err_path = ''
      self.hold = False
      self.join = ''
      self.keep = ''
      self.walltime = '1:00'
      self.proc_format = 'nodes=%d:ppn=%d'
      self.resources = []
      self.mail_ops = ''
      self.email = ''
      self.name = ''
      self.out_path = ''
      self.priority = ''
      self.queue = ''
      self.shell = '/bin/bash'
      self.workdir = "$PBS_O_WORKDIR"
      self.commands = commands
      # Default template file is ~/.pbsdefaults
      if not template: template = path.join(getenv('HOME'), '.pbsdefaults')
      if not path.exists(template):
         raise PBSClassError('Missing %s file. ' % template +
                             'Run pbsjob.py to create it')
         default = open(template, 'r')
         for line in default:
            parts = line.split('=')
            if len(parts) != 2:
               raise PBSClassError('Corrupt ~/.pbsdefaults file! ' +
                                   'Unrecognized line:\n  > %s' % line.strip())
            parts[0], parts[1] = parts[0].strip(), parts[1].strip()
            if not hasattr(self, parts[0]):
               raise PBSClassError('Unknown attribute %s' % parts[0])
            # The resource list has to be handled separately because it's a list
            if parts[0] == 'resources': self.resources.append(parts[1])
            else: setattr(self, parts[0], parts[1]) # parts[x] are all strings
         default.close()

   # ==============================

   def set_walltime(self, time):
      """ Sets the wall-clock based on a given time string, but check to make
          sure it's a valid time (xx:xx:xx)
      """
      valid_time_hours = re.compile(r'[0-9]{1,2}:[0-9]{2}:[0-9]{2}')
      valid_time_mins = re.compile(r'[0-9]{1,2}:[0-9]{2}')

      if valid_time_hours.match(str(time)):
         if not valid_time_hours.sub('', str(time)):
            # We have a valid time in hours
            self.walltime = str(time)
         else:
            # We have extra characters
            raise PBSInputError('Invalid walltime: %s' % time)
      elif valid_time_mins.match(str(time)):
         if not valid_time_mins.sub('', str(time)):
            # We have a valid time in minutes
            self.walltime = str(time)
         else:
            raise PBSInputError('Invalid walltime: %s' % time)
      else:
         raise PBSInputError('Invalid walltime: %s' % time)

   # ==============================

   def set_proc_count(self, *args):
      """ 
      Any arguments not provided necessary in the self.proc_format will be
      replaced with a 1. If too many are provided, an error is raised
      """
      # First we have to get rid of double-%s, since those are escaped.
      remaining_format = re.sub(r'%{2}', '', self.proc_format)
      # Use REs to find how many %s we can replace
      replace_fields = re.compile('%')
      num_fields = len(replace_fields.findall(remaining_format))
      if num_fields < len(args):
         raise PBSInputError('Too many arguments for processor count!')
      if num_fields > len(args):
         fields = list(args[:])
         for i in range(num_fields - len(args)): fields.append(1)
      else:
         fields = list(args[:])

      self.proc_count = self.proc_format % tuple(fields)

   # ==============================

   def header(self): 
      """ Return header. set_proc_count must have been called """
      if not hasattr(self, 'proc_count'):
         raise PBSInputError('set_proc_count not called!')

      my_header =(pbsline('-A', self.account) + pbsline('-a', self.date_time)
                + pbsline('-c', self.interval) + pbsline('-e', self.err_path)
                + pbsline('-o', self.out_path) + pbsline('-j', self.join)
                + pbsline('-k', self.keep) + pbsline('-m', self.mail_ops)
                + pbsline('-M', self.email) + pbsline('-N', self.name)
                + pbsline('-p', self.priority) + pbsline('-q', self.queue)
                + pbsline('-S', self.shell) + pbsline('-l', self.proc_count)
                + pbsline('-l', self.walltime))
      for item in self.resources:
         my_header += pbsline('-l', item)

      return my_header

   # ==============================

   def _get_sub_script(self):
      """ Returns the string that will be sent as the submission script """
      return self.header() + linesep + 'cd %s' % self.workdir + \
             linesep + self.commands

   # ==============================

   def submit(self, print_result=True):
      """ Submits the job """
      import utilities
      qsub = utilities.which('qsub')

      if not qsub: raise PBSMissingError('Cannot find qsub!')

      sub_script = self._get_sub_script()

      process = Popen([qsub], stdin=PIPE, stdout=PIPE, stderr=PIPE)

      (output, error) = process.communicate(sub_script)

      if process.wait():
         raise QsubError('problem submitting job: %s' % error)

      if print_result: print output

      return output

   # ==============================

   def submit_ask(self, print_result=True, after_job=None):
      """ Submits the job, but prints the jobfile and asks first """
      import utilities
      qsub = utilities.which('qsub')
      qstat = utilities.which('qstat')

      if not qsub: raise PBSMissingError('Cannot find qsub!')

      sub_script = self._get_sub_script()

      stdout.write('Going to submit the following script:' + linesep)
      stdout.write('=============' + linesep + sub_script + linesep
                 + '=============' + linesep + 'OK?  > ')
      response = stdin.readline()

      # Determine if we have to submit this with a dependency
      if after_job:
         # Make sure we have qstat, since that's how we check that the job
         # we're depending on exists in the first place
         if not qstat:
            raise PBSMissingError('Cannot find qstat!')
         process = Popen([qstat, after_job], stdin=PIPE, 
                         stdout=PIPE, stderr=PIPE)
         (output, error) = process.communicate('')
         # If we get a non-zero exit status, that job doesn't exist!
         if process.wait():
            raise PBSDependError('Job %s does not exist. Bad dependency!' % 
                                 after_job)
         # If it does exist, 
      if response.lower() == 'yes':
         process = Popen([qsub], stdin=PIPE, stdout=PIPE, stderr=PIPE)

         (output, error) = process.communicate(sub_script)

         if process.wait():
            raise QsubError('problem submitting job: %s' % error)

         if print_result: print output
         
         return output

      return None

   # ==============================

   def print_submit(self, filename=None):
      " Prints the submission script to a file (or stdout if no file provided) "
      if not filename:
         tofile = stdout
      else:
         if type(filename).__name__ == 'str': tofile = open(filename, 'w')
         elif hasattr(filename, 'write'): tofile = filename

      sub_script = self._get_sub_script()

      tofile.write(sub_script)

      if type(filename).__name__ == 'str': tofile.close()
      else: tofile.flush()

   # ==============================

   def set_name(self, name=None):
      """ 
      Sets the name either from a random Final Fantasy character or 
      the passed name 
      """
      from random import choice
      global names

      if not name:
         key = choice(names.keys())
         self.name = choice(names[key])
      else:
         self.name = name

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == '__main__':
   # If this program is executed, it's to set up a skeleton file to load
   parser = OptionParser()
   parser.add_option('-A', '--account', dest='account', default=None,
                     help='Default account to charge jobs to')
   parser.add_option('-e', '--error', dest='error', default=None,
                     help='Default name for stderr dump')
   parser.add_option('-o', '--output', dest='output', default=None,
                     help='Default name for stdout dump')
   parser.add_option('-j', '--join', dest='join', default=None,
                     help='Join output/error streams. Allowed oe, eo, n')
   parser.add_option('-k', '--keep', dest='keep', default=None,
                     help='Whether error/output streams kept on host. ' +
                          'Allowed e, o, eo, oe, n')
   parser.add_option('-m', '--mail-options', dest='mail_options', default=None,
                     help='Send email notification if job (a)borts, when it ' +
                          '(b)egins, and/or when it (e)nds.')
   parser.add_option('-M', '--emails', dest='addresses', default=None,
                     help='Where email notifications are sent')
   parser.add_option('-S', '--shell', dest='shell', default='/bin/bash',
                     help='Default shell to use to execute script. ' +
                          'Default /bin/bash')
   parser.add_option('-q', '--queue', dest='queue', default=None,
                     help='Default queue to set job to')
   parser.add_option('-l', '--resources', dest='resources', default=None,
                     help='Comma-delimited list of resource lines')
   parser.add_option('-p', '--proc-cnt-fmt', dest='proc_cnt_fmt', 
                     default='nodes=%d:ppn=%d',
                     help='Format of the node/processor count resource ' +
                          'option. Default [nodes=%d:ppn=%d]. %d is a ' +
                          'specifiable number')
   parser.add_option('-t', '--template', dest='template', 
                     default=path.join(getenv('HOME'), '.pbsdefaults'),
                     help='Name of the template file to create. ' +
                          'Default ~/.pbsdefaults')
   (opt, args) = parser.parse_args()

   if len(args) != 0: 
      print >> stderr, 'Error: Unknown command-line argument(s) %s' % args
      parser.print_help()
      exit()

   template = TemplateFile(open(opt.template, 'w'))

   if opt.account: template.write("account = '%s'" % opt.account)
   if opt.error: template.write("err_path = '%s'" % opt.error)
   if opt.output: template.write("out_path = '%s'" % opt.output)
   if opt.join: template.write("join = '%s'" % opt.join)
   if opt.keep: template.write("keep = '%s'" % opt.keep)
   if opt.mail_options: template.write("mail_ops = '%s'" % opt.mail_options)
   if opt.addresses: template.write("email = '%s'" % opt.addresses)
   if opt.shell: template.write("shell = '%s'" % opt.shell)
   if opt.queue: template.write("queue = '%s'" % opt.shell)
   if opt.proc_cnt_fmt: template.write("proc_format = '%s'" % opt.proc_cnt_fmt)
   if opt.resources: 
      template.write("resources = '%s'" % opt.resources.split(','))
   
   template.close()

   template = open(opt.template, 'r')
   stuff = template.read()
   template.close()

   print 'Template file %s created:' % opt.template + linesep
   print stuff
