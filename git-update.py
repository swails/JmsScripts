#!/usr/bin/env python

from optparse import OptionParser, OptionGroup
from subprocess import Popen, PIPE
from os import path
import sys

class GitError(Exception): pass

def get_branch_names():
   """ Gets the branch names from "git branch" """
   parent = '.'
   found_git_dir = False
   while path.isdir(parent):
      if path.isdir(path.join(parent, '.git')):
         found_git_dir = True
         break
      if parent == '.': parent == '..'
      else: parent = path.join(parent, '..')
   
   if not found_git_dir:
      raise GitError('Not in a git repository!')

   process = Popen(['git', 'branch'], stdout=PIPE, stderr=PIPE)

   output, error = process.communicate('')

   if process.wait(): raise GitError('Problem running git branch !!')

   branches = output.split('\n')

   ret_branches = []
   for i, branch in enumerate(branches): 
      branch = branch.replace('*','').strip()
      if branch: ret_branches.append(branch)

   return ret_branches

def pull_only():
   """ Just pulls from tracked branch """
   for branch in get_branch_names():
      print 'Changing to branch', branch
      process = Popen(['git', 'checkout', branch], stdout=PIPE, stderr=PIPE)
      out, err = process.communicate('')
      if process.wait():
         raise GitError('Could not change to branch ' + branch)
      print 'Pulling...'
      process = Popen(['git', 'pull'], stdout=PIPE, stderr=PIPE)
      out, err = process.communicate('')
      if process.wait():
         raise GitError('Problem pulling on branch ' + branch)

def merge(pushto=None):
   """ Merges everything with "master" and pushes to a repo if given """
   for branch in get_branch_names():
      if branch.endswith('-with-patches'): continue
      if branch == 'master': continue
      print 'Changing to branch', branch
      process = Popen(['git', 'checkout', branch], stdout=PIPE, stderr=PIPE)
      out, err = process.communicate('')
      if process.wait():
         raise GitError('Could not change to branch ' + branch)
      print 'Merging with master'
      process = Popen(['git', 'merge', 'master'], stdout=PIPE, stderr=PIPE)
      out, err = process.communicate('')
      if process.wait():
         raise GitError('Problem merging branch %s with master' % branch)
      
      if pushto:
         process = Popen(['git', 'push', pushto], stdout=PIPE, stderr=PIPE)
         out, err = process.communicate('')
         if process.wait():
            raise GitError('Problem pushing to ' + pushto)

def main():

   parser = OptionParser()
   group = OptionGroup(parser, 'Operating Modes',
                       'All operating modes are mutually exclusive')
   group.add_option('--pull-only', dest='pull_only', default=False, 
                    action='store_true', help='Do a "git pull" on every branch')
   group.add_option('--merge-only', dest='merge_only', default=False,
                    action='store_true', help='Merge every branch with master ' +
                    'except those ending in "-with-patches"')
   group.add_option('--merge-push', dest='push_to', default=None,
                    help='Merge every branch with master then push to given ' +
                    'repository')
   parser.add_option_group(group)
   
   opt, arg = parser.parse_args()
   
   if len(arg) != 0:
      print 'Bad command-line flags!'
      parser.print_help()
      sys.exit(1)

   elif opt.pull_only + opt.merge_only + bool(opt.push_to) == 0:
      parser.print_help()
      sys.exit(0)

   elif opt.pull_only + opt.merge_only + bool(opt.push_to) > 1:
      print 'Operating Modes are all mutually exclusive!'
      parser.print_help()
      sys.exit(1)

   if opt.pull_only: pull_only()
   elif opt.merge_only: merge()
   elif opt.push_to: merge(opt.push_to)


if __name__ == '__main__': main()
