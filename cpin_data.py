# This python file is used by cpinutil.py and has all of the data for
# each titratable residue for each value of igb that is desired.

def getData(residue, igb):
   import sys, math

# Factors to adjust TI energies
   KB = 0.00199 # kcal/mol * K
   LN_TO_LOG = math.log(10)
   TEMP = 300

   data = []
   charges = []

   if residue == "AS4":
      pKa = 4.0

      # set the relative energies for the different igb choices
      if igb == 2:
         relene = 26.76101 + KB * LN_TO_LOG * TEMP * pKa
      elif igb == 5:
         relene = 26.43128
      else:
         relene = 0.0


      return   [ [ 0,           # Relative energy  STATE 0 Deprotonated
                   0,           # Relative protonation
                  -0.4157,      # N
                   0.2719,      # H
                   0.0341,      # CA
                   0.0864,      # HA
                  -0.1783,      # CB
                  -0.0122,      # HB2
                  -0.0122,      # HB3
                   0.7994,      # CG
                  -0.8014,      # OD1
                  -0.8014,      # OD2
                   0.0000,      # HD21
                   0.5973,      # C
                  -0.5679,      # O
                   0.0000,      # HD22
                   0.0000,      # HD11
                   0.0000 ],    # HD12

                 [ relene,      # Relative energy  STATE 1 Protonated syn O2
                   1,           # Relative protonation
                  -0.4157,      # N
                   0.2719,      # H
                   0.0341,      # CA
                   0.0864,      # HA
                  -0.0316,      # CB
                   0.0488,      # HB2
                   0.0488,      # HB3
                   0.6462,      # CG
                  -0.5554,      # OD1
                  -0.6376,      # OD2
                   0.4747,      # HD21
                   0.5973,      # C
                  -0.5679,      # O
                   0.0000,      # HD22
                   0.0000,      # HD11
                   0.0000 ],    # HD12

                 [ relene,      # Relative energy  STATE 2 Protonated anti O2
                   1,           # Relative protonation
                  -0.4157,      # N
                   0.2719,      # H
                   0.0341,      # CA
                   0.0864,      # HA
                  -0.0316,      # CB
                   0.0488,      # HB2
                   0.0488,      # HB3
                   0.6462,      # CG
                  -0.5554,      # OD1
                  -0.6376,      # OD2
                   0.0000,      # HD21
                   0.5973,      # C
                  -0.5679,      # O
                   0.4747,      # HD22
                   0.0000,      # HD11
                   0.0000 ],    # HD12

                 [ relene,      # Relative energy  STATE 2 Protonated syn O1
                   1,           # Relative protonation
                  -0.4157,      # N
                   0.2719,      # H
                   0.0341,      # CA
                   0.0864,      # HA
                  -0.0316,      # CB
                   0.0488,      # HB2
                   0.0488,      # HB3
                   0.6462,      # CG
                  -0.6376,      # OD1
                  -0.5554,      # OD2
                   0.0000,      # HD21
                   0.5973,      # C
                  -0.5679,      # O
                   0.0000,      # HD22
                   0.4747,      # HD11
                   0.0000 ],    # HD12

                 [ relene,      # Relative energy  STATE 2 Protonated anti O1
                   1,           # Relative protonation
                  -0.4157,      # N
                   0.2719,      # H
                   0.0341,      # CA
                   0.0864,      # HA
                  -0.0316,      # CB
                   0.0488,      # HB2
                   0.0488,      # HB3
                   0.6462,      # CG
                  -0.6376,      # OD1
                  -0.5554,      # OD2
                   0.0000,      # HD21
                   0.5973,      # C
                  -0.5679,      # O
                   0.0000,      # HD22
                   0.0000,      # HD11
                   0.4747 ] ]   # HD12

   if residue == "GL4":
      pKa = 4.4

      if igb == 2:
         relene = 8.3976585 + KB * TEMP * LN_TO_LOG * pKa
      elif igb == 5:
         relene = 0.0
      else:
         relene = 0.0

      return   [ [ 0,           # Relative energy  STATE 0 Deprotonated
                   0,           # Relative protonation
                  -0.4157,      # N 
                   0.2719,      # H
                   0.0145,      # CA
                   0.0779,      # HA
                  -0.0398,      # CB
                  -0.0173,      # HB2
                  -0.0173,      # HB3
                   0.0136,      # CB
                  -0.0425,      # HG2
                  -0.0425,      # HG3
                   0.8054,      # CD
                  -0.8188,      # OE1
                  -0.8188,      # OE2
                   0.0000,      # HE21
                   0.5973,      # C
                  -0.5679,      # O
                   0.0000,      # HE22
                   0.0000,      # HE11
                   0.0000 ],    # HE12
   
                 [ relene,      # Relative energy  STATE 1 Protonated syn O2
                   1,           # Relative protonation
                  -0.4157,      # N 
                   0.2719,      # H
                   0.0145,      # CA
                   0.0779,      # HA
                  -0.0071,      # CB
                   0.0256,      # HB2
                   0.0256,      # HB3
                  -0.0174,      # CB
                   0.0430,      # HG2
                   0.0430,      # HG3
                   0.6801,      # CD
                  -0.5838,      # OE1
                  -0.6511,      # OE2
                   0.4641,      # HE21
                   0.5973,      # C
                  -0.5679,      # O
                   0.0000,      # HE22
                   0.0000,      # HE11
                   0.0000 ],    # HE12

                 [ relene,      # Relative energy  STATE 2 Protonated anti O2
                   1,           # Relative protonation
                  -0.4157,      # N 
                   0.2719,      # H
                   0.0145,      # CA
                   0.0779,      # HA
                  -0.0071,      # CB
                   0.0256,      # HB2
                   0.0256,      # HB3
                  -0.0174,      # CB
                   0.0430,      # HG2
                   0.0430,      # HG3
                   0.6801,      # CD
                  -0.5838,      # OE1
                  -0.6511,      # OE2
                   0.0000,      # HE21
                   0.5973,      # C
                  -0.5679,      # O
                   0.4641,      # HE22
                   0.0000,      # HE11
                   0.0000 ],    # HE12

                 [ relene,      # Relative energy  STATE 3 Protonated syn O1
                   1,           # Relative protonation
                  -0.4157,      # N 
                   0.2719,      # H
                   0.0145,      # CA
                   0.0779,      # HA
                  -0.0071,      # CB
                   0.0256,      # HB2
                   0.0256,      # HB3
                  -0.0174,      # CB
                   0.0430,      # HG2
                   0.0430,      # HG3
                   0.6801,      # CD
                  -0.6511,      # OE1
                  -0.5838,      # OE2
                   0.0000,      # HE21
                   0.5973,      # C
                  -0.5679,      # O
                   0.0000,      # HE22
                   0.4641,      # HE11
                   0.0000 ],    # HE12

                 [ relene,      # Relative energy  STATE 4 Protonated anti O1
                   1,           # Relative protonation
                  -0.4157,      # N 
                   0.2719,      # H
                   0.0145,      # CA
                   0.0779,      # HA
                  -0.0071,      # CB
                   0.0256,      # HB2
                   0.0256,      # HB3
                  -0.0174,      # CB
                   0.0430,      # HG2
                   0.0430,      # HG3
                   0.6801,      # CD
                  -0.6511,      # OE1
                  -0.5838,      # OE2
                   0.0000,      # HE21
                   0.5973,      # C
                  -0.5679,      # O
                   0.0000,      # HE22
                   0.0000,      # HE11
                   0.4641 ] ]   # HE12

   if residue == "TYR":
      pKa = 9.6
      
      if igb == 2:
         relene = -65.0791 + KB * TEMP * LN_TO_LOG * pKa
      elif igb == 5:
         relene = -78.34438
      else:
         relene = 0.0

      return   [ [ 0,           # Relative energy  STATE 0 Protonated
                   1,           # Relative protonation
                  -0.4157,      # N
                   0.2719,      # H
                  -0.0014,      # CA
                   0.0876,      # HA
                  -0.0152,      # CB
                   0.0295,      # HB2
                   0.0295,      # HB3
                  -0.0011,      # CG
                  -0.1906,      # CD1
                   0.1699,      # HD1
                  -0.2341,      # CE1
                   0.1656,      # HE1
                   0.3226,      # CZ
                  -0.5579,      # OH
                   0.3992,      # HH
                  -0.2341,      # CE2
                   0.1656,      # HE2
                  -0.1906,      # CD2
                   0.1699,      # HD2
                   0.5973,      # C
                  -0.5679 ],    # O 

                 [ relene,      # Relative energy  STATE 1 Deprotonated
                   0,           # Relative protonation
                  -0.4157,      # N
                   0.2719,      # H
                  -0.0014,      # CA
                   0.0876,      # HA
                  -0.0858,      # CB
                   0.0190,      # HB2
                   0.0190,      # HB3
                  -0.2130,      # CG
                  -0.1030,      # CD1
                   0.1320,      # HD1
                  -0.4980,      # CE1
                   0.1320,      # HE1
                   0.7770,      # CZ
                  -0.8140,      # OH
                   0.0000,      # HH
                  -0.4980,      # CE2
                   0.1320,      # HE2
                  -0.1030,      # CD2
                   0.1320,      # HD2
                   0.5973,      # C
                  -0.5679 ] ]   # O 

   if residue == "HIP":
      pKa1 = 6.5
      pKa2 = 7.1

      if igb == 2:
         relene1 = -2.80 - KB * TEMP * LN_TO_LOG * pKa1
         relene2 = -6.46 - KB * TEMP * LN_TO_LOG * pKa2
      elif igb == 5:
         relene1 = -78.34438
         relene2 = -78.34438
      else:
         relene1 = 0.0
         relene2 = 0.0

      return   [ [ 0,           # Relative energy  STATE 0 HIP
                   2,           # Relative Protonation
                  -0.3479,      # N
                   0.2747,      # H
                  -0.1354,      # CA
                   0.1212,      # HA
                  -0.0414,      # CB
                   0.0810,      # HB2
                   0.0810,      # HB3
                  -0.0012,      # CG
                  -0.1513,      # ND1
                   0.3866,      # HD1
                  -0.0170,      # CE1
                   0.2681,      # HE1
                  -0.1718,      # NE2
                   0.3911,      # HE2
                  -0.1141,      # CD2
                   0.2317,      # HD2
                   0.7341,      # C
                  -0.5894 ],    # 0

                 [ relene1,     # Relative energy  STATE 0 HID
                   1,           # Relative Protonation
                  -0.3479,      # N
                   0.2747,      # H
                  -0.1354,      # CA
                   0.1212,      # HA
                  -0.1110,      # CB
                   0.0402,      # HB2
                   0.0402,      # HB3
                  -0.0266,      # CG
                  -0.3811,      # ND1
                   0.3649,      # HD1
                   0.2057,      # CE1
                   0.1392,      # HE1
                  -0.5727,      # NE2
                   0.0000,      # HE2
                   0.1292,      # CD2
                   0.1147,      # HD2
                   0.7341,      # C
                  -0.5894 ],    # 0

                 [ relene2,     # Relative energy  STATE 0 HIE
                   1,           # Relative Protonation
                  -0.3479,      # N
                   0.2747,      # H
                  -0.1354,      # CA
                   0.1212,      # HA
                  -0.1012,      # CB
                   0.0367,      # HB2
                   0.0367,      # HB3
                   0.1868,      # CG
                  -0.5432,      # ND1
                   0.0000,      # HD1
                   0.1635,      # CE1
                   0.1435,      # HE1
                  -0.2795,      # NE2
                   0.3339,      # HE2
                   0.2207,      # CD2
                   0.1862,      # HD2
                   0.7341,      # C
                  -0.5894 ] ]   # 0

   if residue == "LYS":
      pKa = 10.4

      if igb == 2:
         relene = -15.21 + KB * TEMP * LN_TO_LOG * pKa
      elif igb == 5:
         relene = -78.34438
      else:
         relene = 0.0

      return   [ [ relene,      # Relative energy  STATE 0 LYS
                   3,           # Relative protonation
                  -0.3479,      # N
                   0.2747,      # H
                  -0.2400,      # CA
                   0.1426,      # HA
                  -0.0094,      # CB
                   0.0362,      # HB2
                   0.0362,      # HB3
                   0.0187,      # CG
                   0.0103,      # HG2
                   0.0103,      # HG3
                  -0.0479,      # CD
                   0.0621,      # HD2
                   0.0621,      # HD3
                  -0.0143,      # CE
                   0.1135,      # HE2
                   0.1135,      # HE3
                  -0.3854,      # NZ
                   0.3400,      # HZ1
                   0.3400,      # HZ2
                   0.3400,      # HZ3
                   0.7341,      # C
                  -0.5894 ],    # O

                 [ 0,           # Relative energy  STATE 1 LYN
                   2,           # Relative protonation
                  -0.3479,      # N
                   0.2747,      # H
                  -0.2400,      # CA
                   0.1426,      # HA
                  -0.10961,     # CB
                   0.0340,      # HB2
                   0.0340,      # HB3
                   0.06612,     # CG
                   0.01401,     # HG2
                   0.01401,     # HG3
                  -0.03768,     # CD
                   0.01155,     # HD2
                   0.01155,     # HD3
                   0.32604,     # CE
                  -0.03358,     # HE2
                  -0.03358,     # HE3
                  -1.03581,     # NZ
                   0.0000,      # HZ1
                   0.38604,     # HZ2
                   0.38604,     # HZ3
                   0.7341,      # C
                  -0.5894 ] ]   # O

      
   return -1
