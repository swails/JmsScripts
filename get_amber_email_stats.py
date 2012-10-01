#!/usr/bin/python

import re, urllib

from urllib import urlopen

mailre = re.compile(r'.*\(([SMTWF][a-z]{2}) *([FMJASOND][a-z]{2}) *(\d+) *(\d{4}).*\)')

allmails = urlopen('http://archive.ambermd.org/all/')

print 'Calculating traffic stats for amber mailing list'

yr_mail_cnt = {}
mnth_mail_cnt = {}
day_mail_cnt = {}

for line in allmails:
   rematch = mailre.match(line)
   if rematch:
      day, month, nday, year = rematch.groups()
      try:
         yr_mail_cnt[int(year)] += 1
      except KeyError:
         yr_mail_cnt[int(year)] = 1
      try:
         mnth_mail_cnt[month] += 1
      except KeyError:
         mnth_mail_cnt[month] = 1
      try:
         day_mail_cnt[day] += 1
      except KeyError:
         day_mail_cnt[day] = 1

print 'Done calculating traffic stats for mailing list:'
print ''
print '     Year |   # Emails'
print '------------------------'
keys = yr_mail_cnt.keys()
keys.sort()
for key in yr_mail_cnt:
   print '%9s | %10d' % (key, yr_mail_cnt[key])

print '\n------------------------\n'
print '    Month |   # Emails'
print '------------------------'
for key in mnth_mail_cnt:
   print '%9s | %10d' % (key, mnth_mail_cnt[key])

print '\n------------------------\n'
print '      Day |   # Emails'
print '------------------------'
for key in day_mail_cnt:
   print '%9s | %10d' % (key, day_mail_cnt[key])
