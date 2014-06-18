#!/usr/bin/python

import re, urllib
from optparse import OptionParser

from urllib import urlopen

mailre = re.compile(r'.*\(([SMTWF][a-z]{2}) *([FMJASOND][a-z]{2}) *(\d+) *(\d{4}).*\)')
cudare = re.compile(r'(cuda|gpu|nvidia|k20|gtx|titan|tesla|k10|2050|2070)', re.I)
cudacount = 0
count = 0

parser = OptionParser()
parser.add_option('-n', '--name', default='Jason Swails', dest='name')
parser.add_option('-t', '--tag', default='swails', dest='tag')

opt, arg = parser.parse_args()

allmails = urlopen('http://archive.ambermd.org/all/')

print 'Calculating traffic stats for amber mailing list'

yr_mail_cnt = {}
cuyr_mail_cnt = {}
mnth_mail_cnt = {}
cumnth_mail_cnt = {}
day_mail_cnt = {}
cuday_mail_cnt = {}

for line in allmails:
   rematch = mailre.match(line)
   cuadd = int(bool(cudare.findall(line)))
   if rematch:
      day, month, nday, year = rematch.groups()
      try:
         yr_mail_cnt[int(year)] += 1
         cuyr_mail_cnt[int(year)] += cuadd
      except KeyError:
         yr_mail_cnt[int(year)] = 1
         cuyr_mail_cnt[int(year)] = cuadd
      try:
         mnth_mail_cnt[month] += 1
         cumnth_mail_cnt[month] += cuadd
      except KeyError:
         mnth_mail_cnt[month] = 1
         cumnth_mail_cnt[month] = cuadd
      try:
         day_mail_cnt[day] += 1
         cuday_mail_cnt[day] += cuadd
      except KeyError:
         day_mail_cnt[day] = 1
         cuday_mail_cnt[day] = cuadd

   if opt.tag is not None and opt.tag in line.lower():
      count += 1

print 'Done calculating traffic stats for mailing list:'
print ''
print '     Year |   # Emails'
print '------------------------' + '-'*15
keys = yr_mail_cnt.keys()
keys.sort()
for key in yr_mail_cnt:
   print '%9s | %10d (%7d CUDA)' % (key, yr_mail_cnt[key], cuyr_mail_cnt[key])

print '\n------------------------' + '-'*15 + '\n'
print '    Month |   # Emails'
print '------------------------' + '-'*15
mnths = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep',
         'Oct', 'Nov', 'Dec']
for key in mnths:
   print '%9s | %10d (%7d CUDA)' % (key, mnth_mail_cnt[key], cumnth_mail_cnt[key])

print '\n------------------------' + '-'*15 + '\n'
print '      Day |   # Emails'
print '------------------------' + '-'*15
days = ['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun']
for key in days:
   print '%9s | %10d (%7d CUDA)' % (key, day_mail_cnt[key], cuday_mail_cnt[key])

if opt.name is not None:
   print '\n'
   print '%-20s%d' % (opt.name, count)
