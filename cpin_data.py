# If you add titratable residues to getData, also add them here. EXP_PKAS is only
# used to screen residues based on maxpKa and minpKa. flags.
TITRATABLE =" AS4 GL4 HIP TYR LYS  CYS DAP DCP DG  DT  AP  CP  G   U   C    A  C5   1UH "
EXP_PKAS   = "4.0 4.4 6.5 9.6 10.4 8.5 3.9 4.3 2.1 9.7 3.9 4.3 2.1 9.3 13.5 13 13.5 0.0"

# This python file is used by cpinutil.py and has all of the data for
# each titratable residue for each value of igb that is desired.

def getData(residue, igb, has_water=False, neighbor_right='none', neighbor_left='none'):
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
         if has_water: relene = 26.6047496 + KB * LN_TO_LOG * TEMP * pKa
         else:         relene = 26.8894581 + KB * LN_TO_LOG * TEMP * pKa
      elif igb == 5:
         if has_water: relene = 26.1881636 + KB * LN_TO_LOG * TEMP * pKa
         else:         relene = 26.5980488 + KB * LN_TO_LOG * TEMP * pKa
      elif igb == 8:
         relene = 26.3448911 + KB * LN_TO_LOG * TEMP * pKa
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
         if has_water: relene = 8.0464682 + KB * TEMP * LN_TO_LOG * pKa
         else:         relene = 8.4057785 + KB * TEMP * LN_TO_LOG * pKa
      elif igb == 5:
         if has_water: relene = 7.6690995 + KB * TEMP * LN_TO_LOG * pKa
         else:         relene = 8.0855764 + KB * TEMP * LN_TO_LOG * pKa
      elif igb == 8:
         relene = 8.3493469 + KB * TEMP * LN_TO_LOG * pKa
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
         if has_water: relene = -65.003415 - KB * TEMP * LN_TO_LOG * pKa
         else:         relene = -65.113428 - KB * TEMP * LN_TO_LOG * pKa
      elif igb == 5:
         if has_water: relene = -64.047229 - KB * TEMP * LN_TO_LOG * pKa
         else:         relene = -64.166385 - KB * TEMP * LN_TO_LOG * pKa
      elif igb == 8:
         relene = -61.3305355 - KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.0

      return   [ [ 0,           # Relative energy  STATE 0 TYR (prot)
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

                 [ relene,      # Relative energy  STATE 1 TYM (deprot)
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
         if has_water:
            relene1 = -2.77641 - KB * TEMP * LN_TO_LOG * pKa1
            relene2 = -6.71562 - KB * TEMP * LN_TO_LOG * pKa2
         else:
            relene1 = -2.84183 - KB * TEMP * LN_TO_LOG * pKa1
            relene2 = -6.58793 - KB * TEMP * LN_TO_LOG * pKa2
      elif igb == 5:
         if has_water:
            relene1 = -2.90517 - KB * TEMP * LN_TO_LOG * pKa1
            relene2 = -6.82684 - KB * TEMP * LN_TO_LOG * pKa2
         else:
            relene1 = -2.86001 - KB * TEMP * LN_TO_LOG * pKa1
            relene2 = -6.70726 - KB * TEMP * LN_TO_LOG * pKa2
      elif igb == 8:
         relene1 = -3.4000 - KB * TEMP * LN_TO_LOG * pKa1
         relene2 = -6.3190 - KB * TEMP * LN_TO_LOG * pKa2
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
         if has_water: relene = -15.1417977 + KB * TEMP * LN_TO_LOG * pKa
         else:         relene = -15.2423959 + KB * TEMP * LN_TO_LOG * pKa
      elif igb == 5:
         if has_water: relene = -14.3162107 + KB * TEMP * LN_TO_LOG * pKa
         else:         relene = -14.5392838 + KB * TEMP * LN_TO_LOG * pKa
      elif igb == 8:
         relene = -18.393654 + KB * TEMP * LN_TO_LOG * pKa
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
         if has_water: relene = 77.5680281 + KB * TEMP * LN_TO_LOG * pKa
         else:         relene = 77.4666763 + KB * TEMP * LN_TO_LOG * pKa
      elif igb == 5:
         if has_water: relene = 76.2827217 + KB * TEMP * LN_TO_LOG * pKa
         else:         relene = 76.2588331 + KB * TEMP * LN_TO_LOG * pKa
      elif igb == 8:
         relene = 71.5804519 + KB * TEMP * LN_TO_LOG * pKa
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
      
   if residue == "DAP":
      pKa = 3.9

      if igb == 2:
         relene = -19.8442 - KB * TEMP * LN_TO_LOG * pKa
      elif igb == 5:
         relene = -19.8442 - KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.0

      return   [ [ relene,      # Relative energy  STATE 0 DA
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

                 [ 0,           # Relative energy  STATE 1 DAP
                   2,           # Relative protonation
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
                   0.0944,      # N9 
                   0.1617,      # C8 
                   0.2281,      # H8
                  -0.5674,      # N7 
                   0.1358,      # C5 
                   0.5711,      # C6 
                  -0.8251,      # N6 
                   0.4456,      # H61
                   0.4456,      # H62
                  -0.5750,      # N1 
                   0.4251,      # C2 
                   0.1437,      # H2 
                  -0.5611,      # N3
                   0.3421,      # C4
                   0.0713,      # C3'
                   0.0985,      # H3'
                  -0.0854,      # C2'
                   0.0718,      # H2'1
                   0.0718,      # H2'2
                  -0.5232,      # O3'
                   0.4301 ] ]   # H1

   if residue == "DCP":
      pKa = 4.3

      if igb == 2:
         relene = -40.526 - KB * TEMP * LN_TO_LOG * pKa
      elif igb == 5:
         relene = -40.526 - KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.0

      return   [ [ relene,      # Relative energy  STATE 0 DC
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

                 [ 0,           # Relative energy  STATE 1 DCP
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
                  -0.0116,      # C1'
                   0.1963,      # H1'
                   0.2167,      # N1
                  -0.0282,      # C6
                   0.2713,      # H6
                  -0.4162,      # C5
                   0.2179,      # H5
                   0.6653,      # C4
                  -0.8590,      # N4
                   0.4598,      # H41
                   0.4598,      # H42
                  -0.4956,      # N3
                   0.5371,      # C2
                  -0.5028,      # O2
                   0.0713,      # C3'
                   0.0985,      # H3'
                  -0.0854,      # C2'
                   0.0718,      # H2'
                   0.0718,      # H2'
                  -0.5232,      # O3'
                   0.4108 ] ]   # H3

   if residue == "DG":
      pKa = 9.2

      if igb == 2:
         relene = -90.0011 - KB * TEMP * LN_TO_LOG * pKa
      elif igb == 5:
         relene = -90.0011 - KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.00

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
                  -0.5232 ] ,   # O3'

                 [ relene,      # Relative energy  STATE 1 DGD
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
                   0.0358,      # C1'
                   0.1746,      # H1'
                  -0.0507,      # N9
                   0.0779,      # C8
                   0.1516,      # H8
                  -0.6122,      # N7
                   0.0806,      # C5
                   0.7105,      # C6
                  -0.7253,      # 06
                  -0.8527,      # N1
                   0.0000,      # H1
                   0.9561,      # C2
                  -0.9903,      # N2
                   0.3837,      # H21
                   0.3837,      # H22
                  -0.8545,      # N3
                   0.2528,      # C4
                   0.0713,      # C3'
                   0.0985,      # H3'
                  -0.0854,      # C2'
                   0.0718,      # H2'
                   0.0718,      # H2'
                  -0.5232 ] ]   # O3'

   if residue == "DT":
      pKa = 9.7

      if igb == 5:
         relene = -56.7729 - KB * TEMP * LN_TO_LOG * pKa
      elif igb == 2:
         relene = -56.7729 - KB * TEMP * LN_TO_LOG * pKa
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

                 [ relene,      # Relative energy  STATE 1 DTD
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
                   0.0680,      # C1'
                   0.1804,      # H1'
                  -0.2861,      # N1
                  -0.1874,      # C6
                   0.2251,      # H6
                  -0.1092,      # C5
                  -0.2602,      # C7
                   0.0589,      # H71
                   0.0589,      # H72
                   0.0589,      # H73
                   0.8263,      # C4
                  -0.7396,      # O4
                  -0.9169,      # N3
                   0.0000,      # H3
                   0.9167,      # C2
                  -0.7722,      # O2
                   0.0713,      # C3'
                   0.0985,      # H3'
                  -0.0854,      # C2'
                   0.0718,      # H2'
                   0.0718,      # H2'
                  -0.5232 ] ]   # O3'

   if residue == "AP":
      pKa = 3.9

      if igb == 5:
         relene = 15.903 + KB * TEMP * LN_TO_LOG * pKa
      elif igb == 2:
         relene = 55.1918 + KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.

      return   [ [ 0,           # Relative energy STATE 0 A
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

                 [ relene,      # Relative energy  STATE 1 AP
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
                   0.0394,      # C1' 
                   0.2007,      # H1' 
                   0.0961,      # N9
                   0.2011,      # C8
                   0.1965,      # H8
                  -0.5569,      # N7
                   0.1136,      # C5
                   0.5845,      # C6
                  -0.8152,      # N6
                   0.4403,      # H61 
                   0.4403,      # H62 
                  -0.5776,      # N1
                   0.4435,      # C2
                   0.1307,      # H2
                  -0.5201,      # N3
                   0.2681,      # C4
                   0.2022,      # C3' 
                   0.0615,      # H3' 
                   0.0670,      # C2' 
                   0.0972,      # H2'1
                  -0.6139,      # O2' 
                   0.4186,      # HO'2
                  -0.5246,      # O3' 
                   0.4310 ] ]   # H1

   if residue == "CP":
      pKa = 4.3

      if igb == 5:
         relene = 40.1407 + KB * TEMP * LN_TO_LOG * pKa
      elif igb == 2:
         relene = 40.1407 + KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.

      return   [ [ 0,           # Relative energy STATE 0 C
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

                 [ relene,      # Relative energy STATE 1 CP
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
                   0.0066,      # C1' 
                   0.2029,      # H1' 
                   0.1954,      # N1
                   0.0028,      # C6
                   0.2366,      # H6
                  -0.4218,      # C5
                   0.2253,      # H5
                   0.6466,      # C4
                  -0.8363,      # N4
                   0.4518,      # H41 
                   0.4518,      # H42 
                  -0.4871,      # N3
                   0.5039,      # C2
                  -0.4753,      # O2
                   0.2022,      # C3' 
                   0.0615,      # H3' 
                   0.0670,      # C2' 
                   0.0972,      # H2'1
                  -0.6139,      # O2' 
                   0.4186,      # HO'2
                  -0.5246,      # O3' 
                   0.4128 ] ]   # H3

   if residue == "G":
      pKa = 9.2

      if igb == 5:
         relene = -96.0454 - KB * TEMP * LN_TO_LOG * pKa
      elif igb == 2:
         relene = -96.0454 - KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.00

      return   [ [ 0,           # Relative energy STATE 0 G
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
                  -0.5246 ] ,   # O3' 

                 [ relene,      # Relative energy STATE 1 GD
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
                   0.0191,      # C1' 
                   0.2006,      # H1' 
                  -0.0623,      # N9
                   0.1479,      # C8
                   0.1137,      # H8
                  -0.6127,      # N7
                   0.0488,      # C5
                   0.7137,      # C6
                  -0.7191,      # O6
                  -0.8557,      # N1
                   0.0000,      # H1
                   0.9976,      # C2
                  -1.0387,      # N2
                   0.3969,      # H21 
                   0.3969,      # H22 
                  -0.8299,      # N3
                   0.1992,      # C4
                   0.2022,      # C3' 
                   0.0615,      # H3' 
                   0.0670,      # C2' 
                   0.0972,      # H2'1
                  -0.6139,      # O2' 
                   0.4186,      # HO'2
                  -0.5246 ] ]   # O3' 

   if residue == "U":
      pKa = 9.3

      if igb == 5:
         relene = -134.883 - KB * TEMP * LN_TO_LOG * pKa
      elif igb == 2:
         relene = -134.883 - KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.00

      return   [ [ 0,           # Relative energy STATE 0 U
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

                 [ relene,      # Relative energy STATE 1 UD
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
                   0.0674,      # C1' 
                   0.1824,      # H1' 
                  -0.2733,      # N1
                   0.0264,      # C6
                   0.1501,      # H6
                  -0.5820,      # C5
                   0.1560,      # H5
                   0.9762,      # C4
                  -0.7808,      # O4
                  -0.9327,      # N3
                   0.0000,      # H3
                   0.8698,      # C2
                  -0.7435,      # O2
                   0.2022,      # C3' 
                   0.0615,      # H3' 
                   0.0670,      # C2'
                   0.0972,      # H2'1
                  -0.6139,      # O2' 
                   0.4186,      # HO'2
                  -0.5246 ] ]   # O3' 

   if residue == "A":
      pKa = 13.3

      if igb == 5:
         relene = 0.00 #- KB * TEMP * LN_TO_LOG * pKa
      elif igb == 2:
         relene = -86.9941 #- KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.00

      return   [ [ 0,           # Relative energy STATE 0 C5 (prot)
                   1,           # Relative protonation
                   1.1662,      # P
                  -0.7760,      # OP1
                  -0.7760,      # OP2
                  -0.4989,      # O5'
                   0.0558,      # C5'
                   0.0679,      # H5'
                   0.0679,      # H5''
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
                   0.0972,      # H2'
                  -0.6139,      # O2'
                   0.4186,      # HO2'
                  -0.5246,  ],  # O3'
                 [ relene,      # Relative energy STATE 0 C5 (prot)
                   0,           # Relative protonation
                   1.1662,      # P
                  -0.7760,      # OP1
                  -0.7760,      # OP2
                  -0.4989,      # O5'
                   0.0558,      # C5'
                   0.0679,      # H5'
                   0.0679,      # H5''
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
                  -0.0311,      # C2'
                   0.0000,      # H2'
                  -1.0000,      # O2'
                   0.0000,      # HO2'
                  -0.5246,  ] ] # O3'
                   
   if residue == "C":
      pKa = 13.5

      if igb == 5:
         relene = 0.00 #- KB * TEMP * LN_TO_LOG * pKa
      elif igb == 2:
         relene = -86.9941 - KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.00

      return   [ [ 0,           # Relative energy STATE 0 C5 (prot)
                   1,           # Relative protonation
                   1.1662,      # P 
                  -0.7760,      # OP1
                  -0.7760,      # OP2
                  -0.4989,      # O5'
                   0.0558,      # C5'
                   0.0679,      # H5'
                   0.0679,      # H5''
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
                   0.0972,      # H2'
                  -0.6139,      # O2'
                   0.4186,      # HO2'
                  -0.5246 ] ,   # O3'

                 [ relene,      # Relative energy STATE 0 C5 (prot)
                   0,           # Relative protonation
                   1.1662,      # P 
                  -0.7760,      # OP1
                  -0.7760,      # OP2
                  -0.4989,      # O5'
                   0.0558,      # C5'
                   0.0679,      # H5'
                   0.0679,      # H5''
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
                  -0.0311,      # C2'
                   0.0000,      # H2'
                  -1.0000,      # O2'
                   0.0000,      # HO2'
                  -0.5246 ] ]   # O3'

   if residue == "C5":
      pKa = 13.5

      if igb == 5:
         relene = 0.00 #- KB * TEMP * LN_TO_LOG * pKa
      elif igb == 2:
         relene = -84.4627 - KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.00

      return   [ [ 0,           # Relative energy STATE 0 C5 (prot)
                   1,           # Relative protonation
                   0.306100,    # HO5' 
                  -0.4989,      # O5'
                   0.0558,      # C5' 
                   0.0679,      # H5' 
                   0.0679,      # H5'' 
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
                   0.0972,      # H2' 
                  -0.6139,      # O2' 
                   0.4186,      # HO2' 
                  -0.5246 ] ,   # O3' 

                 [ relene,      # Relative energy STATE 0 C5 (prot)
                   0,           # Relative protonation
                   0.3061,      # HO5' 
                  -0.4989,      # O5'
                   0.0558,      # C5' 
                   0.0679,      # H5' 
                   0.0679,      # H5'' 
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
                  -0.0311,      # C2' 
                   0.0000,      # H2' 
                  -1.0000,      # O2' 
                   0.0000,      # HO2' 
                  -0.5246 ] ]   # O3' 

   if residue == "1UH":
      pKa = 0

      if igb == 5:
         relene = 0.00 #+ KB * TEMP * LN_TO_LOG * pKa
      elif igb == 2:
         relene = -8.024847 + KB * TEMP * LN_TO_LOG * pKa
      else:
         relene = 0.00

      return   [ [ 0,           # Relative energy STATE 0 1UH (deprot)
                   0,           # Relative protonation
                  -0.0527,      # C1
                  -0.0577,      # C2
                  -0.0704,      # C3
                  -0.0774,      # C4
                  -0.0764,      # C5
                  -0.0794,      # C6
                  -0.7186,      # N7
                   0.1265,      # C8
                  -0.0724,      # C9
                   0.1598,      # C10
                   0.6611,      # C11
                  -0.5699,      # N12
                   0.1244,      # C13
                  -0.1031,      # C14
                  -0.1141,      # C15
                  -0.1071,      # C16
                  -0.6121,      # 017
                   0.1398,      # C18
                   0.1161,      # C19
                   0.0827,      # C20
                  -0.5978,      # 021
                  -0.5549,      # N22
                  -0.0703,      # C23
                   0.6787,      # C24
                  -0.6151,      # O25
                  -0.1106,      # C29
                  -0.1290,      # C30
                  -0.1030,      # C31
                  -0.1950,      # C32
                   0.1111,      # C33
                  -0.0523,      # C34
                  -0.4961,      # 038
                  -0.0598,      # C39
                  -0.1450,      # S74
                  -0.0001,      # C77
                  -0.1260,      # C78
                  -0.1140,      # C79
                  -0.1450,      # C80
                  -0.1110,      # C81
                  -0.1550,      # C82
                   0.3075,      # HNC
                   0.3970,      # HOL
                   0.3175,      # HNM
                   0.4180,      # HO
                   0.0527,      # H
                   0.0387,      # H1
                   0.0357,      # H2
                   0.0407,      # H3
                   0.0407,      # H4
                   0.0377,      # H5
                   0.0377,      # H6
                   0.0417,      # H7
                   0.0627,      # H8
                   0.0567,      # H9          
                   0.0287,      # H10         
                   0.0377,      # H11         
                   0.0747,      # H12         
                   0.0867,      # H13         
                   0.0107,      # H14         
                   0.0667,      # H15         
                   0.0387,      # H16        
                   0.0357,      # H17         
                   0.0687,      # H18         
                   0.0337,      # H19         
                   0.0407,      # H20        
                   0.0507,      # H21        
                   0.0367,      # H22         
                   0.0407,      # H23        
                   0.0247,      # H24         
                   0.0647,      # H25         
                   0.0897,      # H26         
                   0.1137,      # H27         
                   0.0607,      # H28         
                   0.1407,      # H29         
                   0.1500,      # H30         
                   0.1380,      # H31         
                   0.1340,      # H32         
                   0.0607,      # H33        
                   0.0507,      # H34         
                   0.0617,      # H35         
                   0.1430,      # H36         
                   0.1320,      # H37         
                   0.1310,      # H38         
                   0.1310,      # H39         
                   0.1270,      # H40         
                   0.0000 ] ,   # H41         

                 [ relene,      # Relative energy STATE 1 1UH (prot)
                   1,           # Relative protonation
                  -0.0787,      # C1          
                  -0.0657,      # C2          
                  -0.0734,      # C3          
                  -0.0834,      # C4          
                  -0.0794,      # C5          
                  -0.0824,      # C6          
                  -0.6674,      # N7          
                   0.0435,      # C8          
                  -0.0964,      # C9          
                   0.1248,      # C10         
                   0.6471,      # C11         
                  -0.5319,      # N12         
                   0.1184,      # C13         
                  -0.1127,      # C14         
                  -0.1127,      # C15         
                  -0.1127,      # C16        
                  -0.6151,      # O17         
                   0.0658,      # C18         
                   0.1171,      # C19         
                   0.0787,      # C20         
                  -0.5648,      # O21         
                  -0.5559,      # N22         
                  -0.0763,      # C23         
                   0.7087,      # C24         
                  -0.6971,      # O25         
                  -0.1236,      # C29         
                  -0.1340,      # C30         
                  -0.0920,      # C31         
                  -0.1780,      # C32         
                   0.1211,      # C33         
                  -0.0493,      # C34         
                  -0.4911,      # O38        
                  -0.0488,      # C39         
                  -0.1400,      # S74         
                  -0.0281,      # C77         
                  -0.1365,      # C78         
                  -0.1070,      # C79         
                  -0.1240,      # C80         
                  -0.1070,      # C81         
                  -0.1365,      # C82         
                   0.3215,      # HNC         
                   0.4120,      # HOL         
                   0.3445,      # HNM         
                   0.4310,      # HO         
                   0.0917,      # H           
                   0.0527,      # H1          
                   0.0567,      # H2          
                   0.0567,      # H3          
                   0.0527,      # H4          
                   0.0527,      # H5          
                   0.0572,      # H6          
                   0.0572,      # H7          
                   0.0612,      # H8          
                   0.0612,      # H9          
                   0.1027,      # H10         
                   0.0857,      # H11         
                   0.0857,      # H12         
                   0.0977,      # H13         
                   0.0977,      # H14         
                   0.0552,      # H15         
                   0.0552,      # H16        
                   0.0552,      # H17         
                   0.0552,      # H18         
                   0.0552,      # H19         
                   0.0552,      # H20        
                   0.0552,      # H21        
                   0.0552,      # H22         
                   0.0552,      # H23        
                   0.1022,      # H24         
                   0.1022,      # H25         
                   0.1077,      # H26         
                   0.0887,      # H27         
                   0.0922,      # H28         
                   0.0922,      # H29         
                   0.1460,      # H30         
                   0.1500,      # H31         
                   0.1490,      # H32         
                   0.0533,      # H33        
                   0.0533,      # H34         
                   0.0533,      # H35         
                   0.1375,      # H36         
                   0.1440,      # H37         
                   0.1450,      # H38         
                   0.1440,      # H39         
                   0.1375,      # H40         
                   0.5008 ] ]   # H41         

   return -1
