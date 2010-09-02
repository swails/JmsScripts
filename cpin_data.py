# If you add titratable residues to getData, also add them here. EXP_PKAS is only
# used to screen residues based on maxpKa and minpKa. flags.
TITRATABLE = "AS4 GL4 HIP TYR LYS  CYS PDA PDC PDG PDG DT  PRA PRC PRG PRG RU"
EXP_PKAS   = "4.0 4.4 6.5 9.6 10.4 8.5 3.9 4.3 2.1 9.2 9.7 3.9 4.3 2.1 9.2 9.3"

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
         relene = 26.77653 + KB * LN_TO_LOG * TEMP * pKa
      elif igb == 5:
         relene = 26.77068 + KB * LN_TO_LOG * TEMP * pKa
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

                 [ relene,      # Relative energy  STATE 3 Protonated syn O1
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

                 [ relene,      # Relative energy  STATE 4 Protonated anti O1
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
         relene = 8.39766 + KB * TEMP * LN_TO_LOG * pKa
      elif igb == 5:
         relene = 8.21030 + KB * TEMP * LN_TO_LOG * pKa
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
         relene = -65.0791  - KB * TEMP * LN_TO_LOG * pKa
      elif igb == 5:
         relene = -64.21455 - KB * TEMP * LN_TO_LOG * pKa
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
         relene1 = -2.84183 - KB * TEMP * LN_TO_LOG * pKa1
         relene2 = -6.35183 - KB * TEMP * LN_TO_LOG * pKa2
      elif igb == 5:
         relene1 = -2.86001 - KB * TEMP * LN_TO_LOG * pKa1
         relene2 = -6.24932 - KB * TEMP * LN_TO_LOG * pKa2
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

                 [ relene1,     # Relative energy  STATE 1 HID
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

                 [ relene2,     # Relative energy  STATE 2 HIE
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
                  -0.2207,      # CD2
                   0.1862,      # HD2
                   0.7341,      # C
                  -0.5894 ] ]   # 0

   if residue == "LYS":
      pKa = 10.4

      if igb == 2:
         relene = -15.26298 + KB * TEMP * LN_TO_LOG * pKa
      elif igb == 5:
         relene = -14.52805 + KB * TEMP * LN_TO_LOG * pKa
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
                   0.01041,     # HG2
                   0.01041,     # HG3
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

   if residue == "CYS":
      pKa = 8.5

      if igb == 2:
         relene = 77.53571 + KB * TEMP * LN_TO_LOG * pKa
      elif igb == 5:
         relene = 76.28363 + KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.0

      return   [ [ relene,      # Relative energy  STATE 0 CYS
                   1,           # Relative protonation
                  -0.4157,      # N
                   0.2719,      # H
                   0.0213,      # CA
                   0.1124,      # HA
                  -0.1231,      # CB
                   0.1112,      # HB2
                   0.1112,      # HB3
                  -0.3119,      # SG
                   0.1933,      # HG
                   0.5973,      # C
                  -0.5679 ],    # O

                 [ 0,           # Relative energy  STATE 1 CYM
                   0,           # Relative protonation
                  -0.4157,      # N
                   0.2719,      # H
                   0.0213,      # CA
                   0.1124,      # HA
                  -0.3593,      # CB
                   0.1122,      # HB2
                   0.1122,      # HB3
                  -0.8844,      # SG
                   0.0000,      # HG
                   0.5973,      # C
                  -0.5679 ] ]   # O
      
   if residue == "PDA":
      pKa = 3.9

      if igb == 2:
         relene = -36.6696 - KB * TEMP * LN_TO_LOG * pKa
      elif igb == 5:
         relene = -36.5912 - KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.0

      return   [ [ relene,      # Relative energy  STATE 0 DA
                   0,           # Relative protonation
                   1.1659,      # P
                  -0.7761,      # O1P
                  -0.7761,      # O2P
                  -0.4954,      # O5'
                  -0.0069,      # C5'
                   0.0754,      # H5'1
                   0.0754,      # H5'2
                   0.1629,      # C4'
                   0.1176,      # H4'
                  -0.3691,      # O4'
                   0.0431,      # C1'
                   0.1838,      # H1'
                  -0.0268,      # N9 
                   0.1607,      # C8 
                   0.1877,      # H8
                  -0.6175,      # N7 
                   0.0725,      # C5 
                   0.6897,      # C6 
                  -0.9123,      # N6 
                   0.4167,      # H61
                   0.4167,      # H62
                  -0.7624,      # N1 
                   0.5716,      # C2 
                   0.0598,      # H2 
                  -0.7417,      # N3
                   0.3800,      # C4
                   0.0713,      # C3'
                   0.0985,      # H3'
                  -0.0854,      # C2'
                   0.0718,      # H2'1
                   0.0718,      # H2'2
                  -0.5232,      # O3'
                   0.0000 ] ,   # H1

                 [ 0,           # Relative energy  STATE 1 PDA
                   1,           # Relative protonation
                   1.1659,      # P
                  -0.7761,      # O1P
                  -0.7761,      # O2P
                  -0.4954,      # O5'
                  -0.0069,      # C5'
                   0.0754,      # H5'1
                   0.0754,      # H5'2
                   0.1629,      # C4'
                   0.1176,      # H4'
                  -0.3691,      # O4'
                   0.2680,      # C1'
                   0.1104,      # H1'
                  -0.1035,      # N9 
                   0.2234,      # C8 
                   0.1948,      # H8
                  -0.5676,      # N7 
                   0.1712,      # C5 
                   0.4302,      # C6 
                  -0.8094,      # N6 
                   0.4418,      # H61
                   0.4418,      # H62
                  -0.3477,      # N1 
                   0.2742,      # C2 
                   0.1717,      # H2 
                  -0.5181,      # N3
                   0.3820,      # C4
                   0.0713,      # C3'
                   0.0985,      # H3'
                  -0.0854,      # C2'
                   0.0718,      # H2'1
                   0.0718,      # H2'2
                  -0.5232,      # O3'
                   0.3584 ] ]   # H1

   if residue == "PDC":
      pKa = 4.3

      if igb == 2:
         relene = -69.3910 - KB * TEMP * LN_TO_LOG * pKa
      elif igb == 5:
         relene = -69.5117 - KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.0

      return   [ [ relene,      # Relative energy  STATE 0 DC
                   0,           # Relative protonation
                   1.1659,      # P 
                  -0.7761,      # O1P
                  -0.7761,      # O2P
                  -0.4954,      # O5'
                  -0.0069,      # C5'
                   0.0754,      # H5'
                   0.0754,      # H5'
                   0.1629,      # C4'
                   0.1176,      # H4'
                  -0.3691,      # O4'
                  -0.0116,      # C1'
                   0.1963,      # H1'
                  -0.0339,      # N1
                  -0.0183,      # C6
                   0.2293,      # H6
                  -0.5222,      # C5
                   0.1863,      # H5
                   0.8439,      # C4
                  -0.9773,      # N4
                   0.4314,      # H41
                   0.4314,      # H42
                  -0.7748,      # N3
                   0.7959,      # C2
                  -0.6548,      # O2
                   0.0713,      # C3'
                   0.0985,      # H3'
                  -0.0854,      # C2'
                   0.0718,      # H2'
                   0.0718,      # H2'
                  -0.5232,      # O3'
                   0.0000 ] ,   # H3

                 [ 0,           # Relative energy  STATE 1 PDC
                   1,           # Relative protonation
                   1.1659,      # P 
                  -0.7761,      # O1P
                  -0.7761,      # O2P
                  -0.4954,      # O5'
                  -0.0069,      # C5'
                   0.0754,      # H5'
                   0.0754,      # H5'
                   0.1629,      # C4'
                   0.1176,      # H4'
                  -0.3691,      # O4'
                   0.2010,      # C1'
                   0.1551,      # H1'
                   0.0920,      # N1
                  -0.0035,      # C6
                   0.2270,      # H6
                  -0.3471,      # C5
                   0.2102,      # H5
                   0.4660,      # C4
                  -0.8412,      # N4
                   0.4561,      # H41
                   0.4561,      # H42
                  -0.2233,      # N3
                   0.4366,      # C2
                  -0.4901,      # O2
                   0.0713,      # C3'
                   0.0985,      # H3'
                  -0.0854,      # C2'
                   0.0718,      # H2'
                   0.0718,      # H2'
                  -0.5232,      # O3'
                   0.3267 ] ]   # H3

   if residue == "PDG":
      pKa1 = 2.1
      pKa2 = 9.2

      if igb == 2:
         relene1 = -115.296 + KB * TEMP * LN_TO_LOG * (pKa1 - pKa2)
         relene2 = -43.3367 + KB * TEMP * LN_TO_LOG * pKa1
         relene3 = -81.6459 - KB * TEMP * LN_TO_LOG * pKa2
      elif igb == 5:
         relene1 = 0.00
         relene2 = 0.00 + KB * TEMP * LN_TO_LOG * pKa1
         relene3 = 0.00 - KB * TEMP * LN_TO_LOG * pKa2
      else:
         relene1 = 0.00
         relene2 = 0.00
         relene3 = 0.00

      return   [ [ 0,           # Relative energy  STATE 0 DG
                   1,           # Relative protonation
                   1.1659,      # P 
                  -0.7761,      # O1P
                  -0.7761,      # O2P
                  -0.4954,      # O5'
                  -0.0069,      # C5'
                   0.0754,      # H5'
                   0.0754,      # H5'
                   0.1629,      # C4'
                   0.1176,      # H4'
                  -0.3691,      # O4'
                   0.0358,      # C1'
                   0.1746,      # H1'
                   0.0577,      # N9
                   0.0736,      # C8
                   0.1997,      # H8
                  -0.5725,      # N7
                   0.1991,      # C5
                   0.4918,      # C6
                  -0.5699,      # 06
                  -0.5053,      # N1
                   0.3520,      # H1
                   0.7432,      # C2
                  -0.9230,      # N2
                   0.4235,      # H21
                   0.4235,      # H22
                  -0.6636,      # N3
                   0.1814,      # C4
                   0.0713,      # C3'
                   0.0985,      # H3'
                  -0.0854,      # C2'
                   0.0718,      # H2'
                   0.0718,      # H2'
                  -0.5232,      # O3'
                   0.0000 ] ,   # H7

                 [ relene1,     # Relative energy  STATE 1 DG'
                   1,           # Relative protonation
                   1.1659,      # P 
                  -0.7761,      # O1P
                  -0.7761,      # O2P
                  -0.4954,      # O5'
                  -0.0069,      # C5'
                   0.0754,      # H5'
                   0.0754,      # H5'
                   0.1629,      # C4'
                   0.1176,      # H4'
                  -0.3691,      # O4'
                   0.1958,      # C1'
                   0.1261,      # H1'
                   0.0493,      # N9
                   0.0220,      # C8
                   0.2235,      # H8
                  -0.3229,      # N7
                  -0.1270,      # C5
                   0.7360,      # C6
                  -0.6505,      # 06
                  -0.8281,      # N1
                   0.0000,      # H1
                   1.0003,      # C2
                  -0.9620,      # N2
                   0.3881,      # H21
                   0.3881,      # H22
                  -0.8188,      # N3
                   0.3122,      # C4
                   0.0713,      # C3'
                   0.0985,      # H3'
                  -0.0854,      # C2'
                   0.0718,      # H2'
                   0.0718,      # H2'
                  -0.5232,      # O3'
                   0.3895 ] ,   # H7

                 [ relene2,     # Relative energy  STATE 2 PDG
                   2,           # Relative protonation
                   1.1659,      # P 
                  -0.7761,      # O1P
                  -0.7761,      # O2P
                  -0.4954,      # O5'
                  -0.0069,      # C5'
                   0.0754,      # H5'
                   0.0754,      # H5'
                   0.1629,      # C4'
                   0.1176,      # H4'
                  -0.3691,      # O4'
                   0.3192,      # C1'
                   0.1045,      # H1'
                   0.0343,      # N9
                   0.0628,      # C8
                   0.2493,      # H8
                  -0.3278,      # N7
                   0.0170,      # C5
                   0.4831,      # C6
                  -0.4870,      # 06
                  -0.4243,      # N1
                   0.3482,      # H1
                   0.7316,      # C2
                  -0.9299,      # N2
                   0.4527,      # H21
                   0.4527,      # H22
                  -0.6279,      # N3
                   0.2529,      # C4
                   0.0713,      # C3'
                   0.0985,      # H3'
                  -0.0854,      # C2'
                   0.0718,      # H2'
                   0.0718,      # H2'
                  -0.5232,      # O3'
                   0.4102 ] ,   # H7

                 [ relene3,     # Relative energy  STATE 3 DDG
                   0,           # Relative protonation
                   1.1659,      # P 
                  -0.7761,      # O1P
                  -0.7761,      # O2P
                  -0.4954,      # O5'
                  -0.0069,      # C5'
                   0.0754,      # H5'
                   0.0754,      # H5'
                   0.1629,      # C4'
                   0.1176,      # H4'
                  -0.3691,      # O4'
                  -0.1125,      # C1'
                   0.2146,      # H1'
                   0.0046,      # N9
                   0.0330,      # C8
                   0.1611,      # H8
                  -0.5631,      # N7
                   0.0216,      # C5
                   0.7451,      # C6
                  -0.6952,      # 06
                  -0.8645,      # N1
                   0.0000,      # H1
                   0.9199,      # C2
                  -0.9469,      # N2
                   0.3545,      # H21
                   0.3545,      # H22
                  -0.8348,      # N3
                   0.3297,      # C4
                   0.0713,      # C3'
                   0.0985,      # H3'
                  -0.0854,      # C2'
                   0.0718,      # H2'
                   0.0718,      # H2'
                  -0.5232,      # O3'
                   0.0000 ] ]   # H7

   if residue == "DT":
      pKa = 9.7

      if igb == 5:
         relene = -9.71931 - KB * TEMP * LN_TO_LOG * pKa
      elif igb == 2:
         relene = -11.0095 - KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.

      return   [ [ 0,           # Relative energy  STATE 0 DT
                   1,           # Relative protonation
                   1.1659,      # P 
                  -0.7761,      # O1P
                  -0.7761,      # O2P
                  -0.4954,      # O5'
                  -0.0069,      # C5'
                   0.0754,      # H5'
                   0.0754,      # H5'
                   0.1629,      # C4'
                   0.1176,      # H4'
                  -0.3691,      # O4'
                   0.0680,      # C1'
                   0.1804,      # H1'
                  -0.0239,      # N1
                  -0.2209,      # C6
                   0.2607,      # H6
                   0.0025,      # C5
                  -0.2269,      # C7
                   0.0770,      # H71
                   0.0770,      # H72
                   0.0770,      # H73
                   0.5194,      # C4
                  -0.5563,      # O4
                  -0.4340,      # N3
                   0.3420,      # H3
                   0.5677,      # C2
                  -0.5881,      # O2
                   0.0713,      # C3'
                   0.0985,      # H3'
                  -0.0854,      # C2'
                   0.0718,      # H2'
                   0.0718,      # H2'
                  -0.5232 ] ,   # O3'

                 [ relene,      # Relative energy  STATE 1 DDT
                   0,           # Relative protonation
                   1.1659,      # P 
                  -0.7761,      # O1P
                  -0.7761,      # O2P
                  -0.4954,      # O5'
                  -0.0069,      # C5'
                   0.0754,      # H5'
                   0.0754,      # H5'
                   0.1629,      # C4'
                   0.1176,      # H4'
                  -0.3691,      # O4'
                  -0.0689,      # C1'
                   0.1595,      # H1'
                  -0.0523,      # N1
                  -0.2190,      # C6
                   0.1737,      # H6
                  -0.0940,      # C5
                  -0.2575,      # C7
                   0.0706,      # H71
                   0.0742,      # H72
                   0.0476,      # H73
                   0.8340,      # C4
                  -0.7116,      # O4
                  -0.9026,      # N3
                   0.0000,      # H3
                   0.7883,      # C2
                  -0.7204,      # O2
                   0.0713,      # C3'
                   0.0985,      # H3'
                  -0.0854,      # C2'
                   0.0718,      # H2'
                   0.0718,      # H2'
                  -0.5232 ] ]   # O3'

   if residue == "PRA":
      pKa = 3.9

      if igb == 5:
         relene = 0.00
      elif igb == 2:
         relene = 55.1918 + KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.

      return   [ [ 0,           # Relative energy STATE 0 RA
                   0,           # Relative protonation
                   1.1662,      # P 
                  -0.7760,      # O1P
                  -0.7760,      # O2P
                  -0.4989,      # O5'
                   0.0558,      # C5'
                   0.0679,      # H5'1
                   0.0679,      # H5'2
                   0.1065,      # C4'
                   0.1174,      # H4'
                  -0.3548,      # O4'
                   0.0394,      # C1' 
                   0.2007,      # H1' 
                  -0.0251,      # N9
                   0.2006,      # C8
                   0.1553,      # H8
                  -0.6073,      # N7
                   0.0515,      # C5
                   0.7009,      # C6
                  -0.9019,      # N6
                   0.4115,      # H61 
                   0.4115,      # H62 
                  -0.7615,      # N1
                   0.5875,      # C2
                   0.0473,      # H2
                  -0.6997,      # N3
                   0.3053,      # C4
                   0.2022,      # C3' 
                   0.0615,      # H3' 
                   0.0670,      # C2' 
                   0.0972,      # H2'1
                  -0.6139,      # O2' 
                   0.4186,      # HO'2
                  -0.5246,      # O3' 
                   0.0000 ] ,   # H1

                 [ relene,      # Relative energy  STATE 1 PRA
                   1,           # Relative protonation
                   1.1662,      # P
                  -0.7760,      # O1P 
                  -0.7760,      # O2P 
                  -0.4989,      # O5' 
                   0.0558,      # C5' 
                   0.0679,      # H5'1
                   0.0679,      # H5'2
                   0.1065,      # C4' 
                   0.1174,      # H4' 
                  -0.3548,      # O4' 
                   0.1348,      # C1' 
                   0.1715,      # H1' 
                  -0.0168,      # N9
                   0.1724,      # C8
                   0.2119,      # H8
                  -0.5546,      # N7
                   0.1609,      # C5
                   0.4241,      # C6
                  -0.8040,      # N6
                   0.4401,      # H61 
                   0.4401,      # H62 
                  -0.3322,      # N1
                   0.2442,      # C2
                   0.1792,      # H2
                  -0.4921,      # N3
                   0.3811,      # C4
                   0.2022,      # C3' 
                   0.0615,      # H3' 
                   0.0670,      # C2' 
                   0.0972,      # H2'1
                  -0.6139,      # O2' 
                   0.4186,      # HO'2
                  -0.5246,      # O3' 
                   0.3554 ] ]   # H1

   if residue == "PRC":
      pKa = 3.9

      if igb == 5:
         relene = 0.00
      elif igb == 2:
         relene = 80.7889 + KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.

      return   [ [ 0,           # Relative energy STATE 0 RA
                   0,           # Relative protonation
                   1.1662,      # P
                  -0.7760,      # O1P 
                  -0.7760,      # O2P 
                  -0.4989,      # O5' 
                   0.0558,      # C5' 
                   0.0679,      # H5'1
                   0.0679,      # H5'2
                   0.1065,      # C4' 
                   0.1174,      # H4' 
                  -0.3548,      # O4' 
                   0.0066,      # C1' 
                   0.2029,      # H1' 
                  -0.0484,      # N1
                   0.0053,      # C6
                   0.1958,      # H6
                  -0.5215,      # C5
                   0.1928,      # H5
                   0.8185,      # C4
                  -0.9530,      # N4
                   0.4234,      # H41 
                   0.4234,      # H42 
                  -0.7584,      # N3
                   0.7538,      # C2
                  -0.6252,      # O2
                   0.2022,      # C3' 
                   0.0615,      # H3' 
                   0.0670,      # C2' 
                   0.0972,      # H2'1
                  -0.6139,      # O2' 
                   0.4186,      # HO'2
                  -0.5246,      # O3' 
                   0.0000 ] ,   # H3

                 [ relene,      # Relative energy STATE 0 PRA
                   1,           # Relative protonation
                   1.1662,      # P
                  -0.7760,      # O1P 
                  -0.7760,      # O2P 
                  -0.4989,      # O5' 
                   0.0558,      # C5' 
                   0.0679,      # H5'1
                   0.0679,      # H5'2
                   0.1065,      # C4' 
                   0.1174,      # H4' 
                  -0.3548,      # O4' 
                   0.2637,      # C1' 
                   0.1075,      # H1' 
                   0.0676,      # N1
                  -0.0924,      # C6
                   0.2894,      # H6
                  -0.3239,      # C5
                   0.2096,      # H5
                   0.4888,      # C4
                  -0.8518,      # N4
                   0.4584,      # H41 
                   0.4584,      # H42 
                  -0.2408,      # N3
                   0.4507,      # C2
                  -0.5003,      # O2
                   0.2022,      # C3' 
                   0.0615,      # H3' 
                   0.0670,      # C2' 
                   0.0972,      # H2'1
                  -0.6139,      # O2' 
                   0.4186,      # HO'2
                  -0.5246,      # O3' 
                   0.3311 ] ]   # H3

   if residue == "PRG":
      pKa1 = 2.1
      pKa2 = 9.2

      if igb == 5:
         relene1 = 0.00 + KB * TEMP * LN_TO_LOG * (pKa1 - pKa2)
         relene2 = 0.00 + KB * TEMP * LN_TO_LOG * pKa1
         relene3 = 0.00 - KB * TEMP * LN_TO_LOG * pKa2
      elif igb == 2:
         relene1 = -53.5572 + KB * TEMP * LN_TO_LOG * (pKa1 - pKa2)
         relene2 = 0.436125 + KB * TEMP * LN_TO_LOG * pKa1
         relene3 = -27.1978 - KB * TEMP * LN_TO_LOG * pKa2
      else:
         relene1 = 0.00
         relene2 = 0.00
         relene3 = 0.00

      return   [ [ 0,           # Relative energy STATE 0 RG
                   1,           # Relative protonation
                   1.1662,      # P
                  -0.7760,      # O1P 
                  -0.7760,      # O2P 
                  -0.4989,      # O5' 
                   0.0558,      # C5' 
                   0.0679,      # H5'1
                   0.0679,      # H5'2
                   0.1065,      # C4' 
                   0.1174,      # H4' 
                  -0.3548,      # O4' 
                   0.0191,      # C1' 
                   0.2006,      # H1' 
                   0.0492,      # N9
                   0.1374,      # C8
                   0.1640,      # H8
                  -0.5709,      # N7
                   0.1744,      # C5
                   0.4770,      # C6
                  -0.5597,      # O6
                  -0.4787,      # N1
                   0.3424,      # H1
                   0.7657,      # C2
                  -0.9672,      # N2
                   0.4364,      # H21 
                   0.4364,      # H22 
                  -0.6323,      # N3
                   0.1222,      # C4
                   0.2022,      # C3' 
                   0.0615,      # H3' 
                   0.0670,      # C2' 
                   0.0972,      # H2'1
                  -0.6139,      # O2' 
                   0.4186,      # HO'2
                  -0.5246,      # O3' 
                   0.0000 ] ,   # H7

                 [ relene1,     # Relative energy STATE 1 RG'
                   1,           # Relative protonation
                   1.1662,      # P
                  -0.7760,      # O1P 
                  -0.7760,      # O2P 
                  -0.4989,      # O5' 
                   0.0558,      # C5' 
                   0.0679,      # H5'1
                   0.0679,      # H5'2
                   0.1065,      # C4' 
                   0.1174,      # H4' 
                  -0.3548,      # O4' 
                   0.0116,      # C1' 
                   0.1784,      # H1' 
                   0.0919,      # N9
                   0.0623,      # C8
                   0.2250,      # H8
                  -0.3772,      # N7
                  -0.0738,      # C5
                   0.7616,      # C6
                  -0.6575,      # O6
                  -0.7848,      # N1
                   0.0000,      # H1
                   0.8684,      # C2
                  -0.9776,      # N2
                   0.4068,      # H21 
                   0.4068,      # H22 
                  -0.6261,      # N3
                   0.1959,      # C4
                   0.2022,      # C3' 
                   0.0615,      # H3' 
                   0.0670,      # C2' 
                   0.0972,      # H2'1
                  -0.6139,      # O2' 
                   0.4186,      # HO'2
                  -0.5246,      # O3' 
                   0.4043 ] ,   # H7

                 [ relene2,     # Relative energy STATE 2 PRG
                   2,           # Relative protonation
                   1.1662,      # P
                  -0.7760,      # O1P
                  -0.7760,      # O2P
                  -0.4989,      # O5'
                   0.0558,      # C5'
                   0.0679,      # H5'1
                   0.0679,      # H5'2
                   0.1065,      # C4'
                   0.1174,      # H4'
                  -0.3548,      # O4'
                   0.2574,      # C1'
                   0.1120,      # H1'
                   0.0524,      # N9
                   0.0169,      # C8
                   0.2812,      # H8
                  -0.3098,      # N7
                   0.0253,      # C5
                   0.4689,      # C6
                  -0.4801,      # O6
                  -0.3683,      # N1
                   0.3389,      # H1
                   0.6529,      # C2
                  -0.9020,      # N2
                   0.4465,      # H21 
                   0.4465,      # H22 
                  -0.5339,      # N3
                   0.2005,      # C4
                   0.2022,      # C3'
                   0.0615,      # H3'
                   0.0670,      # C2'
                   0.0972,      # H2'1
                  -0.6139,      # O2'
                   0.4186,      # HO'2
                  -0.5246,      # O3'
                   0.4107 ] ,   # H7

                 [ relene3,     # Relative energy STATE 3 DRG
                   0,           # Relative protonation
                   1.1662,      # P
                  -0.7760,      # O1P 
                  -0.7760,      # O2P 
                  -0.4989,      # O5' 
                   0.0558,      # C5' 
                   0.0679,      # H5'1
                   0.0679,      # H5'2
                   0.1065,      # C4' 
                   0.1174,      # H4' 
                  -0.3548,      # O4' 
                  -0.2296,      # C1' 
                   0.2269,      # H1' 
                   0.0467,      # N9
                   0.1604,      # C8
                   0.1079,      # H8
                  -0.6073,      # N7
                   0.0834,      # C5
                   0.7356,      # C6
                  -0.6833,      # O6
                  -0.8027,      # N1
                   0.0000,      # H1
                   0.7903,      # C2
                  -0.9414,      # N2
                   0.3628,      # H21 
                   0.3628,      # H22 
                  -0.6429,      # N3
                   0.1464,      # C4
                   0.2022,      # C3' 
                   0.0615,      # H3' 
                   0.0670,      # C2' 
                   0.0972,      # H2'1
                  -0.6139,      # O2' 
                   0.4186,      # HO'2
                  -0.5246,      # O3' 
                   0.0000 ] ]   # H7

   if residue == "RU":
      pKa = 9.3

      if igb == 5:
         relene = 0.00 - KB * TEMP * LN_TO_LOG * pKa
      elif igb == 2:
         relene = -68.4318 - KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.00

      return   [ [ 0,           # Relative energy STATE 0 RU
                   1,           # Relative protonation
                   1.1662,      # P
                  -0.7760,      # O1P 
                  -0.7760,      # O2P 
                  -0.4989,      # O5' 
                   0.0558,      # C5' 
                   0.0679,      # H5'1
                   0.0679,      # H5'2
                   0.1065,      # C4' 
                   0.1174,      # H4' 
                  -0.3548,      # O4' 
                   0.0674,      # C1' 
                   0.1824,      # H1' 
                   0.0418,      # N1
                  -0.1126,      # C6
                   0.2188,      # H6
                  -0.3635,      # C5
                   0.1811,      # H5
                   0.5952,      # C4
                  -0.5761,      # O4
                  -0.3549,      # N3
                   0.3154,      # H3
                   0.4687,      # C2
                  -0.5477,      # O2
                   0.2022,      # C3' 
                   0.0615,      # H3' 
                   0.0670,      # C2' 
                   0.0972,      # H2'1
                  -0.6139,      # O2'
                   0.4186,      # HO'2
                  -0.5246 ] ,   # O3' 

                 [ relene,      # Relative energy STATE 1 DRU
                   0,           # Relative protonation
                   1.1662,      # P
                  -0.7760,      # O1P 
                  -0.7760,      # O2P 
                  -0.4989,      # O5' 
                   0.0558,      # C5' 
                   0.0679,      # H5'1
                   0.0679,      # H5'2
                   0.1065,      # C4' 
                   0.1174,      # H4' 
                  -0.3548,      # O4' 
                  -0.1497,      # C1' 
                   0.1724,      # H1' 
                  -0.0092,      # N1
                  -0.0919,      # C6
                   0.1707,      # H6
                  -0.5518,      # C5
                   0.1711,      # H5
                   0.9793,      # C4
                  -0.7365,      # O4
                  -0.9169,      # N3
                   0.0000,      # H3
                   0.8510,      # C2
                  -0.7725,      # O2
                   0.2022,      # C3' 
                   0.0615,      # H3' 
                   0.0670,      # C2'
                   0.0972,      # H2'1
                  -0.6139,      # O2' 
                   0.4186,      # HO'2
                  -0.5246 ] ]   # O3' 

   return -1
