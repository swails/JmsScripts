#!/usr/bin/env python

#LAST MODIFIED: 12-09-08

#CONVENTION: angles are in DEGREES


from numpy import *
import os, sys, copy, time

#try to load the compiled functions
USELIB = True
try:
  import geometrylib
except ImportError:
  USELIB = False
  

#globals
DegPerRad = 180./pi
RadPerDeg = pi/180.


#======== VECTOR FUNCTIONS ========

def Length(Vec):
  "Returns the length of a vector."
  if USELIB:
    return geometrylib.length(Vec)
  else:
    return sqrt((Vec*Vec).sum())
  
def UnitVec(Vec):
  "Returns a vector of unit length."
  if USELIB:
    return geometrylib.unitvec(Vec)
  else:
    return Vec / Length(Vec)

def RandVec(Dim = 3):
  "Returns a random unit vector on a sphere."
  while True:
    Vec = 2.*random.rand(Dim) - 1.
    if (Vec*Vec).sum() <= 1: break
  return UnitVec(Vec)


#======== ANGLE FUNCTIONS ========

def NormRad(Rad):
  "Normalizes Rad to be within (-pi,pi)."
  if USELIB:
    return geometrylib.normrad(Rad)
  else:
    return mod(Rad + pi, 2*pi) - pi

def NormDeg(Deg):
  "Normalizes Deg to be within (-180.,180.)."
  if USELIB:
    return geometrylib.normdeg(Deg)
  else:
    return mod(Deg + 180., 360.) - 180.

def NearestAngle(Deg, RefDeg):
  "Returns the nearest Deg +- n*360 to RefDeg."
  return mod(Deg - RefDeg + 180., 360.) + RefDeg - 180.

def RadFromTrig(SinVal, CosVal):
  "Determines the angle in radians from sin and cosine values."
  if USELIB:
    return geometrylib.radfromtrig(SinVal, CosVal)
  else:
    Ang = arccos(clip(CosVal, -1, 1))
    if SinVal < 0: Ang = 2*pi - Ang
    return Ang

def DegFromTrig(SinVal, CosVal):
  "Determines the angle in degrees from sin and cosine values."
  if USELIB:
    return geometrylib.degfromtrig(SinVal, CosVal)
  else:
    return RadFromTrig(SinVal, CosVal) * DegPerRad

def Angle(Pos1, Pos2, Pos3):
  "Calculates the angle formed by three positions."
  if USELIB:
    return geometrylib.angle3deg(Pos1, Pos2, Pos3)
  else:
    if all(Pos1 == Pos2) or all(Pos2 == Pos3): return 0.
    Vec21 = UnitVec(Pos1 - Pos2)
    Vec23 = UnitVec(Pos3 - Pos2)
    Phi = clip(dot(Vec21, Vec23), -1., 1.)
    Phi = arccos(Phi)
    #put angle in right bounds
    return NormRad(Phi) * DegPerRad
    
def AngleVec(Vec1, Vec2):
  "Calculates the angle formed by two vectors."
  return Angle(Vec1, zeros_like(Vec1), Vec2)
  
def CosAngleVec(Vec1, Vec2):
    """Calculates the cosine of the angle formed by two vectors."""
    if all(Vec1 == Vec2): return 0.
    Vec1 = UnitVec(Vec1)
    Vec2 = UnitVec(Vec2)
    return clip(dot(Vec1, Vec2), -1., 1.)    

def Dihedral(Pos1, Pos2, Pos3, Pos4):
  "Calculates the dihedral angle formed by four positions."
  if USELIB:
    return geometrylib.dihedral3deg(Pos1, Pos2, Pos3, Pos4)
  else:
    Vec12 = Pos2 - Pos1
    Vec23 = Pos3 - Pos2
    Vec34 = Pos4 - Pos3
    #subtract off component along 23
    Norm12 = cross(Vec12, Vec23)
    Norm34 = cross(Vec23, Vec34)
    #calculate normalization
    Norm = sqrt( sum(Norm12**2) * sum(Norm34**2) )
    #calculate angle
    Phi = clip(dot(Norm12, Norm34) / Norm, -1., 1.)
    Phi = arccos(Phi)
    if dot(Vec12, Norm34) < 0.: Phi = -Phi
    #put angle in right bounds
    return NormRad(Phi) * DegPerRad

def GetVecMapping(Vec1, Vec2):
  """Returns the axis and angle between two vectors.""
Returns Vec, Ang such that Vec1 = dot(Vec2, RotMat(Vec, Ang))"""
  if USELIB:
    return geometrylib.getvecmappingdeg(Vec1, Vec2)
  else:
    Vec1 = UnitVec(Vec1)
    Vec2 = UnitVec(Vec2)
    CosAng = clip(dot(Vec1, Vec2), -1, 1)
    Ang = arccos(CosAng)
    #special case of parallel
    if CosAng == 1:
      return Vec1.copy(), 0.
    #special case of antiparallel; rotate indices
    elif CosAng == -1:
      NewVec = Vec1.copy()
      NewVec[1:] = Vec1[0:-1]
      NewVec[0] = Vec1[-1]
      NewVec = NewVec - dot(NewVec, Vec1)
      return NewVec, Ang * DegPerRad
    else:
      Vec = cross(Vec1, Vec2)
      if (Vec*Vec).sum() == 0:
        return Vec1.copy(), 0.
      else:
        return Vec, Ang * DegPerRad


#======== SPHERE ROUTINES =======

def SpherePoints(NPoints):
  """Returns NPoints points equally spaced on a unit sphere, using the
golden spiral algorithm.  Return array is dimensions (3,NPoints)"""
  if USELIB:
    return geometrylib.spherepoints(NPoints)
  else:
    Points = []
    inc = pi * (3. - sqrt(5.))
    off = 2. / float(NPoints)
    for k in range(NPoints):
      y = k * off - 1. + (off / 2.)
      r = sqrt(max(1. - y*y, 0.))
      phi = k * inc
      Points.append([cos(phi)*r, y, sin(phi)*r])
    return array(Points, float)

def SphereSurfaceAreas(Pos, Radii, NPoints = 1000):
  """Uses the Shrake-Rupley algorithm to compute the surface area
of a set of spheres with center positions in Pos and radii in Radii.
NPoints is the number of points to use on each sphere.
Returns an array of the surface area per sphere."""
  #make a unit sphere
  Points = SpherePoints(NPoints)
  if USELIB:
    return geometrylib.spheresurfaceareas(Pos, Radii, Points)
  else:
    Areas = zeros_like(Radii)
    for (i,iPos) in enumerate(Pos):
      ThisPoints = Points * Radii[i] + iPos
      AreaPerPoint = 4.*pi*Radii[i]**2 / float(NPoints)
      Exposed = ones(NPoints, int)
      for (j,jPos) in enumerate(Pos):
        if i == j: continue
        if all(Exposed == 0): break
        if sum((iPos - jPos)**2) > (Radii[i] + Radii[j])**2: continue
        DistSq = sum((ThisPoints - jPos)**2, axis=1)
        Exposed[DistSq < Radii[j]**2] = 0
      Areas[i] = AreaPerPoint * sum(Exposed)
    return Areas

def SphereVolumes(Pos, Radii, dx = 1.0):
  """Uses a grid-based algorithm to compute the volume of a set of
overlapping spheres with center positions in Pos and radii in Radii.
dx is grid bin width for dividing up space along each axis.
FineScale is used to do a faster calculation by first coarse-
graining space into boxes of dimension dx and then dx/FineScale.
FineScale only works with the library.
Returns an array of the volume per sphere."""
  if USELIB:
    return geometrylib.spherevolumes(Pos, Radii, dx)
  else:
    Vol = zeros_like(Radii)
    RadiiSq = Radii*Radii
    #make grid points for sphere
    GMin = (Pos - Radii[:,newaxis]).min(axis=0)
    GMax = (Pos + Radii[:,newaxis]).max(axis=0)
    dV = dx**3
    for x in arange(GMin[0], GMax[0] + 0.5 * dx, dx):
      for y in arange(GMin[1], GMax[1] + 0.5 * dx, dx):
        for z in arange(GMin[2], GMax[2] + 0.5 * dx, dx):
          GPos = array([x,y,z], float)
          DistSq = sum((Pos - GPos)**2, axis=1)
          Mask = DistSq < RadiiSq
          if not any(Mask): continue
          DistSq[logical_not(Mask)] = 1.e300
          i = DistSq.argmin()
          Vol[i] += dV
    return Vol
      
      


#======== ROTATION ROUTINES ========

def RotMat(Vec, Ang):
  "Returns the rotation matrix about a vector Ang degrees."
  if USELIB:
    return geometrylib.rotmat3deg(Vec, Ang)
  else:
    M = zeros((3,3), float)
    rcos = cos(Ang * RadPerDeg)
    rsin = sin(Ang * RadPerDeg)
    (u, v, w) = UnitVec(Vec)
    M[0,0] =      rcos + u*u*(1-rcos)
    M[1,0] =  w * rsin + v*u*(1-rcos)
    M[2,0] = -v * rsin + w*u*(1-rcos)
    M[0,1] = -w * rsin + u*v*(1-rcos)
    M[1,1] =      rcos + v*v*(1-rcos)
    M[2,1] =  u * rsin + w*v*(1-rcos)
    M[0,2] =  v * rsin + u*w*(1-rcos)
    M[1,2] = -u * rsin + v*w*(1-rcos)
    M[2,2] =      rcos + w*w*(1-rcos)
    return M

def RotMatEuler(Phi, Theta, Psi):
  "Returns the rotation matrix for Euler angles Phi, Theta, and Psi."
  if USELIB:
    return geometrylib.rotmat3eulerdeg(Phi, Theta, Psi)
  else:
    r = zeros((3,3), float)
    Phi, Theta, Psi = Phi*RadPerDeg, Theta*RadPerDeg, Psi*RadPerDeg
    S1, C1 = sin(Phi), cos(Phi)
    S2, C2 = sin(Theta), cos(Theta)
    S3, C3 = sin(Psi), cos(Psi)
    r[0,0] = C1*C3-S1*C2*S3
    r[0,1] = S1*C3+C1*C2*S3
    r[0,2] = S2*S3
    r[1,0] = -C1*S3-S1*C2*C3
    r[1,1] = -S1*S3+C1*C2*C3
    r[1,2] = S2*C3
    r[2,0] = S1*S2
    r[2,1] = -C1*S2
    r[2,2] = C2
    return r

def RotMatQ(q0, q1, q2, q3):
  "Returns the rotation matrix for quarternions q0, q1, q2, an q3."
  if USELIB:
    return geometrylib.rotmat3q(q0, q1, q2, q3)
  else:
    r = zeros((3,3), float)
    r[0,0] = q0*q0 + q1*q1 - q2*q2 - q3*q3
    r[0,1] = 2*(q1*q2 + q0*q3)
    r[0,2] = 2*(q1*q3 - q0*q2)
    r[1,0] = 2*(q1*q2 - q0*q3)
    r[1,1] = q0*q0 - q1*q1 + q2*q2 - q3*q3
    r[1,2] = 2*(q2*q3 + q0*q1)
    r[2,0] = 2*(q1*q3 + q0*q2)
    r[2,1] = 2*(q2*q3 - q0*q1)
    r[2,2] = q0*q0 - q1*q1 - q2*q2 + q3*q3
    return r

def RotateAboutPoint(Pos, RotMatrix, Point):
  "Rotates a position matrix through a point."
  #translate, rotate the positions, and translate back
  if USELIB:
    if len(Pos.shape) == 1:
      return geometrylib.rotatepointaboutpoint(Pos, RotMatrix, Point)
    else:
      return geometrylib.rotatearrayaboutpoint(Pos, RotMatrix, Point)
  else:
    return dot(Pos - Point, RotMatrix) + Point

def EulerFromRotMat(RotMatrix):
  "Gets the Euler angles implied by a rotation matrix."
  if USELIB:
    return geometrylib.eulerfromrotmat3deg(RotMatrix)
  else:
    Theta = arccos(clip(RotMatrix[2,2], -1, 1)) * DegPerRad
    if RotMatrix[2,2] == 1 or (RotMatrix[2,0] == 0 and RotMatrix[2,1] == 0) \
       or (RotMatrix[0,2] == 0 and RotMatrix[1,2] == 0):
      Psi = 0
      Phi = DegFromTrig(RotMatrix[0,1], RotMatrix[0,0])
    else:
      sphi = sqrt(1 - RotMatrix[2,2] * RotMatrix[2,2])
      Phi = DegFromTrig(RotMatrix[2,0] / sphi, -RotMatrix[2,1] / sphi)
      Psi = DegFromTrig(RotMatrix[0,2] / sphi, RotMatrix[1,2] / sphi)
    return Phi, Theta, Psi

def RandRotMat(dAng):
  "Returns a random rotation matrix."
  Vec = RandVec(3)
  Ang =  dAng * (2*random.random() - 1)
  return RotMat(Vec, Ang)	

def RandRotMatEuler(dAng, Phi, Theta, Psi):
  "Returns a random rotation matrix and new Euler angles."
  r = dot(RotMatEuler(Phi, Theta, Psi), RandRotMat(dAng))
  Phi, Theta, Psi = EulerFromRotMat(r)
  return r, Phi, Theta, Psi

def RotMatMapped(Pos1, Pos2):
  """Extracts the rotation matrix from two arrays of three points each.
Returns RotMat such that Pos1 is aligned to dot(Pos2, RotMat)"""
  #check for identity
  if all(Pos1 == Pos2):
    return identity(Pos1.shape[1], Pos1.dtype)
  #align first vector
  Vec, Ang = GetVecMapping(Pos1[1] - Pos1[0], Pos2[1] - Pos2[0])
  rm = RotMat(Vec, Ang)
  #align pos2 to pos1
  Pos2 = dot(Pos2, rm)
  #from the remaining vec, get the perpindicular components to the aligned
  vperp = UnitVec(Pos1[1] - Pos1[0])
  v1 = Pos1[2] - Pos1[1]
  v2 = Pos2[2] - Pos2[1]
  v1 = v1 - dot(v1, vperp) * vperp
  v2 = v2 - dot(v2, vperp) * vperp
  #check to make sure points are not colinear
  if Length(v1) == 0 or Length(v2) == 0:
    return rm
  else:
    Vec, Ang = GetVecMapping(v1, v2)
    rm = dot(rm, RotMat(Vec, Ang))
    return rm

def RotMatRMSD(Pos1, Pos2, Center = True):
  """Extracts the rotation matrix from two arrays of position coordinates,
that would minimize the RMSD when both arrays are aligned.
Returns RotMat such that Pos1 is aligned to dot(Pos2, RotMat)"""
  #check for identity
  if all(Pos1 == Pos2):
    return identity(Pos1.shape[1], Pos1.dtype)
  #make copies
  p1 = Pos1.copy() 
  p2 = Pos2.copy()
  if Center:
    p1 = p1 - p1.mean(axis=0)
    p2 = p2 - p2.mean(axis=0)
  d1, d2 = shape(p1)
  #calculate correlation matrix
  C = dot(p2.transpose(), p1)
  #get singular value decomp
  V, S, Wt = linalg.svd(C)
  #check determinants
  Vdet = linalg.det(V)
  Wtdet = linalg.det(Wt)
  #correct for possible reflection
  Wt[-1,:] = Wt[-1,:] * (Wtdet * Vdet)
  #return rotation matrix
  return dot(V, Wt)


def AlignmentRMSD(Pos1, Pos2, Center = True):
  """Returns the translation vectors, rotation matrix, and sum of
residuals squared (Pos1Vec, Pos2Vec, RotMat, Resid), for aligning
Pos1 to Pos2, such that Pos1 + Pos1Vec is aligned to
dot(Pos2 + Pos2Vec, RotMat)."""
  d1, d2 = Pos1.shape
  #get centers
  if Center:
    Pos1Vec = -average(Pos1, axis=0)
    Pos2Vec = -average(Pos2, axis=0)
    p1 = Pos1 + Pos1Vec
    p2 = Pos2 + Pos2Vec
  else:
    Pos1Vec, Pos2Vec = zeros(d2, Pos1.dtype), zeros(d2, Pos2.dtype)
    p1, p2 = Pos1, Pos2
  #check for identity
  if all(p1 == p2):
    return Pos1Vec, Pos2Vec, identity(d2, Pos1.dtype), 0.
  #calculate E0
  E0 = sum(p1*p1, axis=None) + sum(p2*p2, axis=None)
  #calculate correlation matrix
  C = dot(transpose(p2), p1)
  #get singular value decomp
  V, S, Wt = linalg.svd(C)
  #check determinants
  Vdet = linalg.det(V)
  Wtdet = linalg.det(Wt)
  #if it's a reflection, reflect along lowest eigenvalue
  if Vdet*Wtdet < 0.: S[-1] = -S[-1]
  #compute risiduals
  Residuals = max(E0 - 2. * sum(S), 0.)
  #correct for possible reflection
  Wt[-1,:] = Wt[-1,:] * (Wtdet * Vdet)
  #calculate rotation matrix
  U = dot(V, Wt)
  return Pos1Vec, Pos2Vec, U, Residuals


def RMSD(Pos1, Pos2, Center = True):
  """Returns root mean squared displacement for aligning
Pos1 to Pos2, such that Pos1 + Pos1Vec is aligned to
dot(Pos2 + Pos2Vec, RotMat)."""
  d1, d2 = Pos1.shape
  #get centers
  if Center:
    p1 = Pos1 - average(Pos1, axis=0)
    p2 = Pos2 - average(Pos2, axis=0)
  else:
    p1, p2 = Pos1, Pos2
  #check for identity
  if all(p1 == p2) or d1==0: return 0.
  #calculate E0
  E0 = sum(p1*p1, axis=None) + sum(p2*p2, axis=None)
  #calculate correlation matrix
  C = dot(transpose(p2), p1)
  #get singular value decomp
  V, S, Wt = linalg.svd(C)
  #check determinants
  Vdet = linalg.det(V)
  Wtdet = linalg.det(Wt)
  #if it's a reflection, reflect along lowest eigenvalue
  if Vdet*Wtdet < 0.: S[-1] = -S[-1]
  #compute risiduals
  Residuals = max(E0 - 2. * sum(S), 0.)
  return sqrt(Residuals / d1)

def dRMSD(Pos1, Pos2, Mask = None):
  """Returns the distance-based RMSD."""
  if USELIB:
      if Mask is None:
          return geometrylib.drms(Pos1, Pos2)
      else:
          return geometrylib.drmsmask(Pos1, Pos2, Mask)
  else:
      d1, d2 = Pos1.shape
      Residuals = 0.
      if not Mask is None: Mask = array(Mask, bool)
      for i in range(d1):
          Dists1 = sqrt(sum((Pos1[i:i+1,:] - Pos1)**2, axis=1))
          Dists2 = sqrt(sum((Pos2[i:i+1,:] - Pos2)**2, axis=1))
          if not Mask is None:
              Dists1 = Dists1[Mask[i,:]]
              Dists2 = Dists2[Mask[i,:]]
          Residuals += sum((Dists1-Dists2)**2)
      return sqrt(Residuals / (d1*(d1-1)))
      
      
#======== MISC ROUTINES ========
def GetContactMap(Pos, Radius):
    """Returns a contact map for a position matrix."""
    NInd = len(Pos)
    Map = zeros((NInd, NInd), float)
    for i in range(0, NInd):
      Vecs = Pos[i+1:,:] - Pos[i,:]
      Vecs = (Vecs*Vecs).sum(axis=1)
      Vals = (Vecs < Radius*Radius).astype(float)
      Map[i, i+1:] = Vals
    #update other half of contact map
    Map = Map + Map.transpose()
    return Map
  

def TestLib():
  "Runs comparison tests between the lib and the slower, noncompiled routines."
  NLoop = 1000
  global USELIB
  if not USELIB:
    print "Library currently not loaded."
    return
  v1 = array([1,-3,5], float)
  v2 = array([2,-6,-2], float)
  v3 = array([-1,-1,5], float)
  v4 = array([3,2,1], float)
  M = RotMat(v1, 30.)
  Pos = random.rand(15).reshape((5,3))
  Pos2 = random.rand(15).reshape((5,3))
  Radii = random.rand(5)
  for i in range(0,2):
    if i == 0:
      USELIB = False
      print "Not using library..."
    else:
      USELIB = True
      print "Using library..."
    StartTime = time.time()
    for j in range(NLoop):
      l = []
      l.extend([Length(v1), UnitVec(v1)])
      l.extend([NormRad(1.2*pi), NormRad(-1.2*pi), NormDeg(-185.), NormDeg(185.)])
      l.extend([RadFromTrig(1,0), RadFromTrig(0,-1), RadFromTrig(sin(1.), cos(1.))])
      l.extend([DegFromTrig(1,0), DegFromTrig(0,-1), DegFromTrig(sin(1.), cos(1.))])
      l.extend([Dihedral(v1, v2, v3, v4), Dihedral(v2, v4, v1, v3)])
      l.extend([GetVecMapping(v1, v2)])
      l.extend([abs(RotMat(v1, 51.2)).sum(), abs(RotMatEuler(32., 24., 80.)).sum(), \
            abs(RotMatQ(1., .2, .4, .894)).sum()])
      l.extend([abs(RotateAboutPoint(Pos, M, v2)).sum()])
      l.extend([EulerFromRotMat(M)])
      l.extend([SpherePoints(5)])
    l.extend([SphereSurfaceAreas(Pos, Radii, 1000)])
    l.extend([SphereVolumes(Pos, Radii, 1.)])
    l.extend([SphereVolumes(Pos, Radii, 0.1)])
    l.extend([dRMSD(Pos, Pos2)])
    for x in l: print x
    print "Elapsed time: %.3f sec\n" % (time.time() - StartTime)
  
  
  
