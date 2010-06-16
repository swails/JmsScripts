#!/usr/bin/env python

import sys

if len(sys.argv) != 3:
   print 'Usage: checkCharges.py cpin1 cpin2'
   sys.exit()

cpin1 = open(sys.argv[1],'r')
cpin2 = open(sys.argv[2],'r')

cpin1lines = cpin1.readlines()
cpin2lines = cpin2.readlines()

cpin1.close()
cpin2.close()

cpin1_chrg = []
cpin2_chrg = []

for x in range(len(cpin1lines)):
   if "chrgdat" in cpin1lines[x].lower():
      y = x + 1
      words = cpin1lines[x].strip().split(',')
      words.pop()
      words[0] = words[0][8:]
      for z in range(len(words)):
         cpin1_chrg.append(words[z])

      end = False

      while not end:
         words = cpin1lines[y].strip().split(',')
         words.pop()

         for z in range(len(words)):
            if "PROTCNT" in words[z]:
               end = True
               break
            else:
               cpin1_chrg.append(words[z])

         y+=1

for x in range(len(cpin2lines)):
   if "chrgdat" in cpin2lines[x].lower():
      y = x + 1
      words = cpin2lines[x].strip().split(',')
      words.pop()
      words[0] = words[0][8:]
      for z in range(len(words)):
         cpin2_chrg.append(words[z])

      end = False

      while not end:
         words = cpin2lines[y].strip().split(',')
         words.pop()

         for z in range(len(words)):
            if "PROTCNT" in words[z]:
               end = True
               break
            else:
               cpin2_chrg.append(words[z])

         y+=1

if len(cpin1_chrg) != len(cpin2_chrg):
   print 'Not the same number of charges!'
   sys.exit()

for x in range(len(cpin1_chrg)):
   if float(cpin1_chrg[x]) != float(cpin2_chrg[x]):
      strng = 'FALSE'
   else:
      strng = ''
   print '{0}\t\t{1} \t\t {2}'.format(cpin1_chrg[x],cpin2_chrg[x],strng)
