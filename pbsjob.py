from os import environ
from sys import stderr, stdout

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PBS_usage = """
-a date_time      {When job can be executed}
-A account_string {Where to charge job to}
-c interval       {checkpointing}
-e path           {stderr path}
-h                {hold when submitted}
-j join           {oe, eo, n: mix stderr/out} 
-k keep           {oe, o, e, n: keep stderr/out}
-l resource_list  {walltime, nproc, etc.}
-m mail_options   {a,b,e: when to mail notification}
-M email_addys    {where to mail notifications}
-N name           {name of job}
-o path           {stdout path}
-p priority       {-1024 <= priority <= 1023}
-q destination    {queue}
-S path           {executing shell}
"""

PBS_omissions = """
-C directive_prefix  {prefix to recognize pbs directives}
-I                   {run jobs interactively}
-r y|n               {re-runable}
-u user_list         {who is running job}
-v variable_list     {expands environment of job}
-W more attributes   {dependencies, etc.}
"""

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Add to these when you find systems with new ways of calling
# number of processors
PBS_nproc = ["nproc", "size"]

# Make a dictionary of names
names = {"FF2":["Frionel", "Maria", "Guy", "Minh", "Gordon", "Layla", "Richard", "Lionheart"],
         "FF3":["Sara", "Cid", "Desh", "Elia", "Allus", "Dorga", "Unne"],
         "FF4":["Cecil", "Kain", "Rosa", "Cid", "Rydia", "Tellah", "Edward", "Yang", "Palom", "Porom", "Edge", "Fusoya", "Golbez"],
         "FF5":["Bartz", "Reina", "Galuf", "Faris", "Krile", "Ex-death", "Gilgamesh", "Cid"],
         "FF6":["Terra", "Locke", "Edgar", "Sabin", "Shadow", "Celes", "Cyan", "Gau", "Setzer", 
                "Mog", "Umaro", "Gogo", "Leo", "Kefka", "Strago", "Relm", "Kefka"],
         "FF7":["Cloud", "Barret", "Tifa", "Aeris", "RedXIII", "Yuffie", "CaitSith", "Vincent", "Cid", "Sephiroth"],
         "FF8":["Squall", "Rinoa", "Seifer", "Quistis", "Zell", "Irvine", "Selphie"],
         "FF9":["Zidane", "Vivi", "Garnet", "Adelbert", "Freya", "Eiko", "Quina", "Amarant", "Beatrix", "Queen", "Kuja", "Tantalus"],
         "FF10":["Tidus", "Yuna", "Auron", "Rikku", "Wakka", "Lulu", "Kimahri", "Maester Seymour"],
         "FF102":["Yuna", "Rikku", "Paine", "Nooj", "Lenne", "ShuyinYuna", "Rikku", "Paine", "Nooj", "Lenne", "Shuyin"],
         "FF12":["Vaan", "Ashe", "Basch", "Penelo", "Balthier", "Fran"],
         "FF13":["Lightning", "Serah", "Hope", "Snow", "Oerba", "Sazh", "Nora"] }

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def pbsline(option, value):
   return "#PBS -%s %s \n" % (option, value)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PBSjob:
   # attributes of each PBS job
   date_time = ''
   account = ''
   interval = ''
   err_path = ''
   hold = False
   join = ''
   keep = ''
   walltime = ''
   resources = []
   mail_ops = ''
   email = ''
   name = ''
   out_path = ''
   priority = 1234321
   queue = ''
   shell = '/bin/bash'
   workdir = "$PBS_O_WORKDIR"
   command = []

   valid = False

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def __init__(self, def_script="-1"):
      # first try to load default options
      optionlines = []
      try:
         default_options = open(environ["HOME"] + '/.pbsdefaults','r')
         optionlines = default_options.readlines()
         default_options.close()
      except IOError:
         pass

      for i in range(len(optionlines)):
         optionlines[i] = optionlines[i].strip("#PBS").strip()
         if optionlines[i][0:2] == "-a":
            self.date = optionlines[i][2:].strip()
         elif optionlines[i][0:2] == "-A":
            self.account = optionlines[i][2:].strip()
         elif optionlines[i][0:2] == "-c":
            self.interval = optionlines[i][2:].strip()
         elif optionlines[i][0:2] == "-e":
            self.err_path = optionlines[i][2:].strip()
         elif optionlines[i][0:2] == "-h":
            self.hold = True
         elif optionlines[i][0:2] == "-j":
            self.join = optionlines[i][2:].strip()
         elif optionlines[i][0:2] == "-k":
            self.keep = optionlines[i][2:].strip()
         elif optionlines[i][0:2] == "-l":
            self.resources.append(optionlines[i][2:].strip())
         elif optionlines[i][0:2] == "-m":
            self.mail_ops = optionlines[i][2:].strip()
         elif optionlines[i][0:2] == "-M":
            self.email = optionlines[i][2:].strip()
         elif optionlines[i][0:2] == "-N":
            self.name = optionlines[i][2:].strip()
         elif optionlines[i][0:2] == "-o":
            self.out_path = optionlines[i][2:].strip()
         elif optionlines[i][0:2] == "-p":
            self.priority = optionlines[i][2:].strip()
         elif optionlines[i][0:2] == "-q":
            self.queue = optionlines[i][2:].strip()
         elif optionlines[i][0:2] == "-S":
            self.shell = optionlines[i][2:].strip()
         else:
            print >> stderr, '"%s" unrecognized' % optionlines[i]

      if def_script != "-1":
         self.script(def_script)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def script(self, script_name):
      try:
         script = open(script_name,'r')
      except IOError:
         print >> stderr, 'Error: "%s" script file not found!'
         return -1

      for line in script:
         if line.strip().startswith("#!"):
            self.shell = line.strip().strip("#!")
            continue
         self.command.append(line)
      if len(self.command[len(self.command)-1].strip()) != 0:
         self.command.append('\n')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def addCommand(self, command):
      self.command.append(command + '\n')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def submit(self):
      from os import system
      submit_string = self.preview()
      if self.valid:
         system('qsub << EOF\n%s\nEOF' % submit_string)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def preview(self):
      cmd = ''
      errors = []
      # shell is always specified
      cmd += pbsline('S', self.shell)
      if self.date_time != '':
         cmd += pbsline('a', self.date_time)
      if self.account != '':
         cmd += pbsline('A', self.account)
      if self.interval != '':
         cmd += pbsline('c', self.interval)
      if self.err_path != '':
         cmd += pbsline('e', self.err_path)
      if self.hold:
         cmd += pbsline('h', '')
      if self.join != '':
         cmd += pbsline('j', self.join)
      if self.keep != '':
         cmd += pbsline('k', self.keep)
      if self.mail_ops != '':
         cmd += pbsline('m', self.mail_ops)
      if self.email != '':
         cmd += pbsline('M', self.email)
      if self.name != '':
         cmd += pbsline('N', self.name)
      if self.out_path != '':
         cmd += pbsline('o', self.out_path)
      if self.priority != 1234321:
         if self.priority < -1024 or self.priority > 1023:
            errors.append("Error: -p <priority>; priority must be between -1024 and 1023, inclusive!")
         else:
            cmd += pbsline('p', self.priority)
      if self.queue != '':
         cmd += pbsline('q', self.queue)
      for i in range(len(self.resources)): # print all resources
         cmd += pbsline('l', self.resources[i])

      if len(self.command) == 0:
         errors.append('Error: You must specify command to submit!')

      cmd += 'cd %s\n'%self.workdir

      for i in range(len(self.command)):
         cmd += self.command[i]

      self.valid = True
      for i in range(len(errors)):
         print >> stderr, errors[i]
         self.valid = False

      return cmd

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def makeDefault(self):
      cmd += pbsline('S', self.shell)
      if self.date_time != '':
         cmd += pbsline('a', self.date_time)
      if self.account != '':
         cmd += pbsline('A', self.account)
      if self.interval != '':
         cmd += pbsline('c', self.interval)
      if self.err_path != '':
         cmd += pbsline('e', self.err_path)
      if self.hold:
         cmd += pbsline('h', '')
      if self.join != '':
         cmd += pbsline('j', self.join)
      if self.keep != '':
         cmd += pbsline('k', self.keep)
      if self.mail_ops != '':
         cmd += pbsline('m', self.mail_ops)
      if self.email != '':
         cmd += pbsline('M', self.email)
      if self.name != '':
         cmd += pbsline('N', self.name)
      if self.out_path != '':
         cmd += pbsline('o', self.out_path)
      if self.priority != 1234321:
         if self.priority < -1024 or self.priority > 1023:
            errors.append("Error: -p <priority>; priority must be between -1024 and 1023, inclusive!")
         else:
            cmd += pbsline('p', self.priority)
      if self.queue != '':
         cmd += pbsline('q', self.queue)
      for i in range(len(self.resources)): # print all resources
         cmd += pbsline('l', self.resources[i])

      file = open(environ["HOME"] + '/.pbsdefaults','w')
      file.write(cmd + '\n')
      file.close()

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   def randomName(self, source="any"):
      from random import randint

      ky = names.keys()
      valid_key = True

      if source != "any":
         try:
            test = names[source.upper()][0]
         except KeyError:
            print >> stderr, 'Error: Invalid name key! Choosing a random one.'
            valid_key = False

      if source == "any" or not valid_key:
         kynum = randint(0,len(ky)-1)
         source = ky[kynum]

      nnum = randint(0,len(names[source])-1)

      self.name = names[source][nnum]
      print >> stdout, 'I have chosen the name "%s" for your job...'%self.name

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def walltime(self, time):
      """ Add the default walltime """
      self.resources.append['walltime=%s']

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
