#!/usr/bin/env python

#######################################################################
#                                                                     #
# This utility prints the amino acid decomposition for a given prmtop #
#                                                                     #
#    Written by Jason Swails, 04/26/2010                              #
#                                                                     #
#######################################################################


import sys, utilities

if len(sys.argv) != 2 or '-help' in sys.argv[1]:
   print 'AA_decomp.py <prmtop> || <pdb>'
   sys.exit()

if sys.argv[1].lower().endswith("pdb"):
   residue_decomp = utilities.getresdecmp_pdb(sys.argv[1])
else:
   residue_decomp = utilities.getresdecmp_prmtop(sys.argv[1])
# residue indices are : ALA, ARG, ASN, ASP, CYS, GLU, GLN, GLY
#                       HIS, ILE, LEU, LYS, MET, PHE, PRO, SER
#                       THR, TRP, VAL, WAT, UNK

print 'ALA:     ' + str(residue_decomp[0])
print 'ARG:     ' + str(residue_decomp[1])
print 'ASN:     ' + str(residue_decomp[2])
print 'ASP:     ' + str(residue_decomp[3])
print 'CYS:     ' + str(residue_decomp[4])
print 'GLU:     ' + str(residue_decomp[5])
print 'GLN:     ' + str(residue_decomp[6])
print 'GLY:     ' + str(residue_decomp[7])
print 'HIS:     ' + str(residue_decomp[8])
print 'ILE:     ' + str(residue_decomp[9])
print 'LEU:     ' + str(residue_decomp[10])
print 'LYS:     ' + str(residue_decomp[11])
print 'MET:     ' + str(residue_decomp[12])
print 'PHE:     ' + str(residue_decomp[13])
print 'PRO:     ' + str(residue_decomp[14])
print 'SER:     ' + str(residue_decomp[15])
print 'THR:     ' + str(residue_decomp[16])
print 'TRP:     ' + str(residue_decomp[17])
print 'TYR:     ' + str(residue_decomp[18])
print 'VAL:     ' + str(residue_decomp[19])
print 'WAT:     ' + str(residue_decomp[20])
print '\nUNK:     ' + str(residue_decomp[21])       
