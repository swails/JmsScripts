"""
This module contains objects that deal with system topology and parameters, such
as atoms, residues, bonds, angles, etc.

These are intended for use with the AmberParm class and are used to recompute
the parameter topology file when necessary (i.e., it provides high-level access
to Amber prmtop data)

Designed so it can be used like:
   from chemistry.amber.topologyobjects import *

by Jason Swails
"""

from chemistry.exceptions import BondError, DihedralError, AmberParmError
from chemistry.amber.constants import NATOM

__all__ = ['Atom', 'Bond', 'BondType', 'Angle', 'AngleType', 'Dihedral',
           'DihedralType', 'Residue', 'ResidueList', 'AtomList', 'BondTypeList',
           'AngleTypeList', 'DihedralTypeList', 'TrackedList']

class Atom(object):
   """ 
   An atom. Only use these as elements in AtomList instances, since AtomList
   will keep track of when indexes and other stuff needs to be updated
   """
   #===================================================

   def __init__(self, parm, starting_index):
      self.bond_partners = []
      self.angle_partners = []
      self.dihedral_partners = []
      self.exclusion_partners = [] # For arbitrary exclusions
      self.parm = parm
      self.idx = -1
      self.starting_index = starting_index
      self.load_from_parm()
      self.residue = None
      self.bonds, self.angles, self.dihedrals = [], [], []
      self.marked = 0 # For setting molecules
   
   #===================================================

   def add_data(self):
      """ 
      Writes this atom's data to the AmberParm object. Don't pitch a fit if
      we're missing some useless (unused) flags.
      """
      # Determine how many excluded atoms we have. The only ones that count are
      # those with a smaller index (to avoid double-counting)
      numex = 0
      for atm in self.bond_partners:
         if atm.idx > self.idx: numex += 1
      for atm in self.angle_partners:
         if atm.idx > self.idx: numex += 1
      for atm in self.dihedral_partners:
         if atm.idx > self.idx: numex += 1
      for atm in self.exclusion_partners:
         if atm.idx > self.idx: numex += 1
      # For some reason, existing topology files follow the convention that
      # atoms with no exclusions (because all bonded partners have atom #s lower
      # than theirs) have num_excluded = 1, with a 0 placeholder in
      # EXCLUDED_ATOMS_LIST... Weird.
      if numex == 0: numex = 1
      # Make sure we're indexing from 0
      self.parm.parm_data['ATOM_NAME'][self.idx] = self.atname[:4]
      self.parm.parm_data['CHARGE'][self.idx] = self.charge
      self.parm.parm_data['MASS'][self.idx] = self.mass
      self.parm.parm_data['ATOM_TYPE_INDEX'][self.idx] = self.nb_idx
      self.parm.parm_data['NUMBER_EXCLUDED_ATOMS'][self.idx] = numex
      self.parm.parm_data['AMBER_ATOM_TYPE'][self.idx] = self.attype[:4]
      self.parm.parm_data['JOIN_ARRAY'][self.idx] = 0
      self.parm.parm_data['TREE_CHAIN_CLASSIFICATION'][self.idx] = self.tree[:4]
      self.parm.parm_data['IROTAT'][self.idx] = 0
      self.parm.parm_data['RADII'][self.idx] = self.radii
      self.parm.parm_data['SCREEN'][self.idx] = self.screen

   #===================================================

   def load_from_parm(self):
      """ Load data from the AmberParm class """
      self.atname = self.parm.parm_data['ATOM_NAME'][self.starting_index]
      self.charge = self.parm.parm_data['CHARGE'][self.starting_index]
      self.mass = self.parm.parm_data['MASS'][self.starting_index]
      self.nb_idx = self.parm.parm_data['ATOM_TYPE_INDEX'][self.starting_index]
      self.attype = self.parm.parm_data['AMBER_ATOM_TYPE'][self.starting_index]
      self.tree = self.parm.parm_data['TREE_CHAIN_CLASSIFICATION'] \
                                                          [self.starting_index]
      if 'RADII' in self.parm.parm_data:
         self.radii = self.parm.parm_data['RADII'][self.starting_index]
      else:
         self.radii = 0.0 # dummy number
      if 'SCREEN' in self.parm.parm_data:
         self.screen = self.parm.parm_data['SCREEN'][self.starting_index]
      else:
         self.radii = 0.0 # dummy number
      # Load the positions and velocities if the amberParm object them
      if hasattr(self.parm, 'coords'):
         self.xx = self.parm.coords[self.starting_index*3  ]
         self.xy = self.parm.coords[self.starting_index*3+1]
         self.xz = self.parm.coords[self.starting_index*3+2]
         if self.parm.hasvels:
            self.vx = self.parm.vels[self.starting_index*3  ]
            self.vy = self.parm.vels[self.starting_index*3+1]
            self.vz = self.parm.vels[self.starting_index*3+2]

   #===================================================
      
   def bond_to(self, other):
      """ 
      Log this atom as bonded to another atom. Check if this has already been
      added to the angle list. If so, remove it from there.
      """
      if self == other:
         raise BondError("Cannot bond atom to itself!")
      if other in self.angle_partners:
         del self.angle_partners[self.angle_partners.index(other)]
      if other in self.dihedral_partners:
         del self.dihedral_partners[self.dihedral_partners.index(other)]
      if other in self.bond_partners: return
      self.bond_partners.append(other)

   #===================================================
      
   def angle_to(self, other):
      """
      Log this atom as angled to another atom. Check if this has already been
      added to the bond list. If so, do nothing
      """
      if self == other:
         raise BondError("Cannot angle an atom with itself!")
      if other in self.dihedral_partners:
         del self.dihedral_partners[self.dihedral_partners.index(other)]
      if other in self.bond_partners or other in self.angle_partners:
         return
      self.angle_partners.append(other)
   
   #===================================================

   def dihedral_to(self, other):
      """
      Log this atom as dihedral-ed to another atom. Check if this has already been
      added to the bond or angle list. If so, do nothing
      """
      if self == other:
         raise BondError("Cannot dihedral an atom with itself!")
      if other in self.bond_partners or other in self.angle_partners: return
      if other in self.dihedral_partners: return
      self.dihedral_partners.append(other)
      
   #===================================================

   def exclude(self, other):
      """
      Add one atom to my exclusion list, even if it's not bonded to, angled to, or
      dihedraled to it
      """
      if self == other:
         raise BondError("Cannot exclude an atom from itself")
      if (other in self.bond_partners or other in self.angle_partners or
          other in self.dihedral_partners or other in self.exclusion_partners):
         return
      self.exclusion_partners.append(other)
      # If he is excluded from me, then I am excluded from him
      other.exclude(self)

   #===================================================

   def __eq__(self, other):
      return id(self) == id(other)
      
   def __ne__(self, other):
      return not Atom.__eq__(self, other)

   def __gt__(self, other):
      return self.idx > other.idx

   def __lt__(self, other):
      return self.idx < other.idx

   def __ge__(self, other):
      return not Atom.__lt__(self, other)

   def __le__(self, other):
      return not Atom.__gt__(self, other)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Bond(object):
   """ Bond class. Stores 2 atoms involved and force constant/equil value """

   #===================================================

   def __init__(self, atom1, atom2, bond_type):
      """ Bond constructor """
      # Make sure we're not bonding me to myself
      if atom1 == atom2:
         raise BondError('Cannot bond atom to itself!')
      # Order the atoms so the lowest atom # is first
      self.atom1 = atom1
      self.atom2 = atom2
      # Register each as bonded to the other
      self.atom1.bond_to(atom2)
      self.atom2.bond_to(atom1)
      self.atom1.bonds.append(self)
      self.atom2.bonds.append(self)
      # Load the force constant and equilibrium distance
      self.bond_type = bond_type

   #===================================================

   def write_info(self, parm, key, idx):
      """ Writes the bond info to the topology file. idx starts at 0 """
      parm.parm_data[key][3*idx  ] = 3*(self.atom1.idx)
      parm.parm_data[key][3*idx+1] = 3*(self.atom2.idx)
      parm.parm_data[key][3*idx+2] = self.bond_type.idx + 1
      # Now have this bond type write its info
#     self.bond_type.write_info(parm)

   #===================================================
   
   def __contains__(self, thing):
      """ Quick and easy way to see if an Atom is in this Bond """
      return id(thing) == id(self.atom1) or id(thing) == id(self.atom2)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BondType(object):
   """ A bond type """

   #===================================================

   def __init__(self, k, req, idx):
      """BondType constructor. idx must start from 0!!! """
      self.idx = idx
      self.k = k
      self.req = req

   #===================================================

   def write_info(self, parm):
      """ Writes the bond parameters in the parameter file """
      # If our index is -1 (we're not being used), just return
      if self.idx == -1: return
      parm.parm_data['BOND_FORCE_CONSTANT'][self.idx] = self.k
      parm.parm_data['BOND_EQUIL_VALUE'][self.idx] = self.req

   #===================================================

   def __eq__(self, other):
      return self.k == other.k and self.req == other.req

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Angle(object):
   """ Angle class. Stores 3 atoms involved and force constant/equil value """
      
   #===================================================

   def __init__(self, atom1, atom2, atom3, angle_type):
      """ Angle constructor """
      # Make sure we're not angling me to myself
      if atom1 == atom2 or atom1 == atom3 or atom2 == atom3:
         raise BondError('Cannot angle atom to itself!')
      self.atom1 = atom1
      self.atom2 = atom2
      self.atom3 = atom3
      # Register each as angled to the others
      self.atom1.angle_to(self.atom2)
      self.atom1.angle_to(self.atom3)
      self.atom2.angle_to(self.atom1)
      self.atom2.angle_to(self.atom3)
      self.atom3.angle_to(self.atom1)
      self.atom3.angle_to(self.atom2)
      self.atom1.angles.append(self)
      self.atom2.angles.append(self)
      self.atom3.angles.append(self)
      # Load the force constant and equilibrium angle
      self.angle_type = angle_type

   #===================================================

   def write_info(self, parm, key, idx):
      """ Write the info to the topology file """
      parm.parm_data[key][4*idx  ] = 3*(self.atom1.idx)
      parm.parm_data[key][4*idx+1] = 3*(self.atom2.idx)
      parm.parm_data[key][4*idx+2] = 3*(self.atom3.idx)
      parm.parm_data[key][4*idx+3] = self.angle_type.idx + 1
      # Now have this bond type write its info
#     self.angle_type.write_info(parm)

   #===================================================
   
   def __contains__(self, thing):
      """
      Quick and easy way to see if an AngleType or Atom is in this Angle
      """
      return (id(thing) == id(self.atom1) or id(thing) == id(self.atom2) or
              id(thing) == id(self.atom3))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AngleType(object):
   """ An angle type """
   #===================================================

   def __init__(self, k, theteq, idx):
      """ AngleType constructor. idx must start from 0!!! """
      self.k = k
      self.theteq = theteq
      self.idx = idx

   #===================================================

   def write_info(self, parm):
      """ Writes the bond parameters in the parameter file """
      # If we're not being used (idx == -1) just return
      if self.idx == -1: return
      parm.parm_data['ANGLE_FORCE_CONSTANT'][self.idx] = self.k
      parm.parm_data['ANGLE_EQUIL_VALUE'][self.idx] = self.theteq

   #===================================================

   def __eq__(self, other):
      return self.k == other.k and self.theteq == other.theteq

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Dihedral(object):
   " Dihedral class with 4 atoms involved and force constant/periodicity/phase "
      
   #===================================================

   def __init__(self, atom1, atom2, atom3, atom4, dihed_type, signs):
      """ Dihedral constructor. idx must start from 0!!! """
      # Make sure we're not dihedraling me to myself
      atmlist = [atom1, atom2, atom3, atom4]
      for i in range(len(atmlist)):
         for j in range(i+1, len(atmlist)):
            if atmlist[i] == atmlist[j]:
               raise BondError('Cannot dihedral atom to itself!')
      # Set up instances
      self.atom1 = atom1
      self.atom2 = atom2
      self.atom3 = atom3
      self.atom4 = atom4
      self.atom1.dihedral_to(atom2)
      self.atom1.dihedral_to(atom3)
      self.atom1.dihedral_to(atom4)
      self.atom2.dihedral_to(atom1)
      self.atom2.dihedral_to(atom3)
      self.atom2.dihedral_to(atom4)
      self.atom3.dihedral_to(atom2)
      self.atom3.dihedral_to(atom1)
      self.atom3.dihedral_to(atom4)
      self.atom4.dihedral_to(atom2)
      self.atom4.dihedral_to(atom3)
      self.atom4.dihedral_to(atom1)
      self.atom1.dihedrals.append(self)
      self.atom2.dihedrals.append(self)
      self.atom3.dihedrals.append(self)
      self.atom4.dihedrals.append(self)
      # Load the force constant and equilibrium angle
      self.dihed_type = dihed_type
      self.signs = signs # is our 3rd or 4th term negative?

   #===================================================

   def write_info(self, parm, key, idx):
      """ Write the info to the topology file """
      # In order to preserve multi-term and improper dihedral definitions, we
      # have to make sure that atom index 0 is _never_ in either of the last
      # 2 positions (since you can't have '-0'. Therefore, if atoms 3 or 4 are
      # index 0 and their sign is negative, swap the dihedral
      if (self.atom3.idx == 0 and self.signs[0] == -1) or \
         (self.atom4.idx == 0 and self.signs[1] == -1):
         parm.parm_data[key][5*idx  ] = 3*(self.atom4.idx)
         parm.parm_data[key][5*idx+1] = 3*(self.atom3.idx)
         parm.parm_data[key][5*idx+2] = 3*(self.atom2.idx) * self.signs[0]
         parm.parm_data[key][5*idx+3] = 3*(self.atom1.idx) * self.signs[1]
      else:
         parm.parm_data[key][5*idx  ] = 3*(self.atom1.idx)
         parm.parm_data[key][5*idx+1] = 3*(self.atom2.idx)
         parm.parm_data[key][5*idx+2] = 3*(self.atom3.idx) * self.signs[0]
         parm.parm_data[key][5*idx+3] = 3*(self.atom4.idx) * self.signs[1]
      parm.parm_data[key][5*idx+4] = self.dihed_type.idx + 1

   #===================================================

   def __contains__(self, thing):
      """ 
      Quick and easy way to find out if either a DihedralType or Atom is
      in this Dihedral
      """
      return id(thing) == id(self.atom1) or id(thing) == id(self.atom2) or \
             id(thing) == id(self.atom3) or id(thing) == id(self.atom4)

   #===================================================

   def __eq__(self, thing):
      """
      A dihedral is equivalent if the 4 atoms are the same (or reverse) in order

      Allow comparison with another type of dihedral or with a list of 4 atoms
      (or tuple)
      """
      if isinstance(thing, Dihedral):
         # I'm comparing with another Dihedral here
         return ( (self.atom1 == thing.atom1 and self.atom2 == thing.atom2 and
                   self.atom3 == thing.atom3 and self.atom4 == thing.atom4) or
                  (self.atom1 == thing.atom4 and self.atom2 == thing.atom3 and
                   self.atom4 == thing.atom1) )
      if isinstance(thing, list) or isinstance(thing, tuple):
         # Here, atoms are expected to index from 0 (Python standard) if we
         # are comparing with a list or tuple
         if len(thing) != 4:
            raise DihedralError('comparative %s has %d elements! Expect 4.' % 
                                (type(thing).__name__, len(thing)))
         # Compare starting_index, since we may not have an index right now...
         return ( (self.atom1.starting_index == thing[0] and 
                   self.atom2.starting_index == thing[1] and
                   self.atom3.starting_index == thing[2] and
                   self.atom4.starting_index == thing[3]) or
                  (self.atom1.starting_index == thing[3] and 
                   self.atom2.starting_index == thing[2] and
                   self.atom3.starting_index == thing[1] and
                   self.atom4.starting_index == thing[0]) )

      raise TypeError('Cannot compare Dihedral with %s' % type(thing).__name__)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DihedralType(object):
   """ A type of dihedral """

   #===================================================
   
   def __init__(self, phi_k, per, phase, scee, scnb, idx):
      """ DihedralType constructor """
      self.phi_k = phi_k
      self.per = per
      self.phase = phase
      self.scee = scee
      self.scnb = scnb
      self.idx = idx

   #===================================================

   def write_info(self, parm):
      """ Write out the dihedral parameters """
      # If our idx == -1, we're not being used, so just return
      if self.idx == -1: return
      parm.parm_data['DIHEDRAL_FORCE_CONSTANT'][self.idx] = self.phi_k
      parm.parm_data['DIHEDRAL_PERIODICITY'][self.idx] = self.per
      parm.parm_data['DIHEDRAL_PHASE'][self.idx] = self.phase
      try:
         parm.parm_data['SCEE_SCALE_FACTOR'][self.idx] = self.scee
      except KeyError:
         pass
      try:
         parm.parm_data['SCNB_SCALE_FACTOR'][self.idx] = self.scnb
      except KeyError:
         pass
   
   #===================================================

   def __eq__(self, other):
      return (self.phi_k == other.phi_k and self.per == other.per and 
              self.phase == other.phase and self.scee == other.scee and
              self.scnb == other.scnb)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Residue(object):
   """ Residue class """

   #===================================================

   def __init__(self, resname, idx):
      self.resname = resname
      self.idx = idx

   #===================================================

   def __contains__(self, atom):
      return atom in self.atoms

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ResidueList(list):
   """ Array of Residues. """

   #===================================================

   def __init__(self, parm):
      self.parm = parm
      list.__init__(self, [Residue(self.parm.parm_data['RESIDUE_LABEL'][i], -1)
                    for i in range(self.parm.ptr('nres'))])
   
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AtomList(list):
   """ Array of Atoms """
   #===================================================

   def __init__(self, parm, fill_from=None):
      self.parm = parm
      if fill_from is None:
         list.__init__(self, [Atom(self.parm, i) for i in
                              range(self.parm.ptr('natom'))])
      else:
         list.__init__(self, [0 for i in range(self.parm.ptr('natom'))])
         for i, atm in enumerate(fill_from): self[i] = atm
      self.changed = False

   #===================================================

   def __delitem__(self, idx):
      """ Deletes this atom then re-indexes everybody else """
      self[idx].idx = -1
      list.__delitem__(self, idx)
      self.changed = True

   #===================================================
   
   def unmark(self):
      """ Unmark all atoms in this list """
      for atm in self: atm.marked = 0

   #===================================================
   
   def _index_us(self):
      """ We have deleted an atom, so now we have to re-index everybody """
      for i in range(len(self)): self[i].idx = i

   #===================================================

   def append(self, item):
      """ Don't allow this! """
      raise AmberParmError("Cannot add to an AtomList!")

   #===================================================

   def write_to_parm(self):
      """ Writes all of the atom data to the topology file """
      # Write all of the arrays here
      self.parm.parm_data['POINTERS'][NATOM] = len(self)
      # Array slices are faster than copy() and creating new arrays
      # each time
      zeros = [0 for i in range(self.parm.parm_data['POINTERS'][NATOM])]
      self.parm.parm_data['ATOM_NAME'] = zeros[:]
      self.parm.parm_data['CHARGE'] = zeros[:]
      self.parm.parm_data['MASS'] = zeros[:]
      self.parm.parm_data['ATOM_TYPE_INDEX'] = zeros[:]
      self.parm.parm_data['NUMBER_EXCLUDED_ATOMS'] = zeros[:]
      self.parm.parm_data['AMBER_ATOM_TYPE'] = zeros[:]
      self.parm.parm_data['JOIN_ARRAY'] = zeros[:]
      self.parm.parm_data['TREE_CHAIN_CLASSIFICATION'] = zeros[:]
      self.parm.parm_data['IROTAT'] = zeros[:]
      self.parm.parm_data['RADII'] = zeros[:]
      self.parm.parm_data['SCREEN'] = zeros[:]
      self._index_us()
      self._determine_exclusions()
      for atm in self: 
         atm.add_data()
         atm.starting_index = atm.idx # arrays are updated...

   #===================================================

   def _determine_exclusions(self):
      """
      Figures out the EXCLUDED_ATOMS_LIST. Only do this right before you write
      the topology file, since it's expensive
      """
      self.parm.parm_data['EXCLUDED_ATOMS_LIST'] = []
      # We have to do something different for extra points. See the top of
      # extra_pts.f in the sander src/ directory. Effectively, the EP is
      # excluded from every atom that the atom it's attached to is excluded
      # from. So here we go through and exclude every extra point from every
      # other atom my bonded pair is excluded from:
      for atm in self:
         if not atm.attype[:2] in ['EP', 'LP']: continue
         partner = atm.bond_partners[0]
         # Now add all bond partners
         for patm in partner.bond_partners:
            # Don't add myself
            if patm is atm: continue
            atm.exclude(patm)
         # Now add all angle partners
         for patm in partner.angle_partners: atm.exclude(patm)
         # Now add all dihedral partners
         for patm in partner.dihedral_partners: atm.exclude(patm)
         # Now add all other arbitrary exclusions
         for patm in partner.exclusion_partners:
            if patm is atm: continue
            atm.exclude(patm)

      for atm in self:
         vals_to_add = []
         for member in atm.bond_partners:
            if member.idx > atm.idx: vals_to_add.append(member.idx+1)
         for member in atm.angle_partners:
            if member.idx > atm.idx: vals_to_add.append(member.idx+1)
         for member in atm.dihedral_partners:
            if member.idx > atm.idx: vals_to_add.append(member.idx+1)
         # Enable additional (arbitrary) exclusions
         for member in atm.exclusion_partners:
            if member.idx > atm.idx: vals_to_add.append(member.idx+1)
         vals_to_add.sort()
         # See comment above about numex = 0 --> numex = 1
         if not vals_to_add: vals_to_add = [0]
         self.parm.parm_data['EXCLUDED_ATOMS_LIST'].extend(vals_to_add)

   #===================================================

   def refresh_data(self):
      """
      Re-loads all data in parm.parm_data for each atom in case we changed
      any of it
      """
      for atm in self:
         atm.load_from_parm()

   #===================================================

   def __setitem__(self, idx, thing):
      """
      This means we changed things...
      """
      self.changed = True
      list.__setitem__(self, idx, thing)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _TypeList(list):
   """ Base class for all type lists """

   #===================================================

   def __init__(self, parm):
      """ Constructs a list of bond types from the topology file """
      self.parm = parm
      self._make_array()
      self.changed = False

   #===================================================

   def __delitem__(self, idx):
      """ Removes a bond from the bond type """
      list.__delitem__(self, idx)

   #===================================================
   
   def append(self, item):
      """ Appending changes our TypeList """
      self.changed = True
      list.append(self, item)

   #===================================================

   def extend(self, items):
      """ Extending changes our TypeList """
      self.changed = True
      list.extend(self, items)
   
   #===================================================

   def _make_array(self):
      """ 
      This method fills self with whichever element we need. This MUST be
      overwritten, so I force it here
      """
      raise AmberParmError('Failure to override method')

   #===================================================

   def write_to_parm(self):
      """ Writes the data here to the parm data """
      for item in self: item.write_info(self.parm)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BondTypeList(_TypeList):
   """ Bond type list """

   #===================================================

   def _make_array(self):
      list.__init__(self,
                    [BondType(self.parm.parm_data['BOND_FORCE_CONSTANT'][i],
                               self.parm.parm_data['BOND_EQUIL_VALUE'][i], -1)
                               for i in range(self.parm.ptr('numbnd'))]
                   )
      
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AngleTypeList(_TypeList):
   """ Angle type list """

   #===================================================

   def _make_array(self):
      list.__init__(self,
                    [AngleType(self.parm.parm_data['ANGLE_FORCE_CONSTANT'][i],
                                self.parm.parm_data['ANGLE_EQUIL_VALUE'][i], -1)
                                for i in range(self.parm.ptr('numang')) ])
      
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DihedralTypeList(_TypeList):
   """ Dihedral type list """

   #===================================================

   def _make_array(self):
      if not 'SCEE_SCALE_FACTOR' in self.parm.parm_data.keys() or \
         not 'SCNB_SCALE_FACTOR' in self.parm.parm_data.keys():
         list.__init__(self,
              [DihedralType(self.parm.parm_data['DIHEDRAL_FORCE_CONSTANT'][i],
               self.parm.parm_data['DIHEDRAL_PERIODICITY'][i], 
               self.parm.parm_data['DIHEDRAL_PHASE'][i], 1.2, 2.0, -1)
               for i in range(self.parm.ptr('nptra')) ])
      else:
         list.__init__(self,
              [DihedralType(self.parm.parm_data['DIHEDRAL_FORCE_CONSTANT'][i],
               self.parm.parm_data['DIHEDRAL_PERIODICITY'][i], 
               self.parm.parm_data['DIHEDRAL_PHASE'][i], 
               self.parm.parm_data['SCEE_SCALE_FACTOR'][i],
               self.parm.parm_data['SCNB_SCALE_FACTOR'][i], -1)
               for i in range(self.parm.ptr('nptra')) ])

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TrackedList(list):
   """
   This creates a list type that allows you to see if anything has changed
   """
   def __init__(self, arg=[]):
      self.changed = False
      list.__init__(self, arg)

   def __delitem__(self, idx):
      self.changed = True
      list.__delitem__(self, idx)

   def append(self, stuff):
      self.changed = True
      list.append(self, stuff)

   def extend(self, stuff):
      self.changed = True
      list.extend(self, stuff)

   def __delslice__(self, i, j):
      self.changed = True
      list.__delslice__(self, i, j)

   def __setitem__(self, i, y):
      self.changed = True
      list.__setitem__(self, i, y)

