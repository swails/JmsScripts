'''
    Dan Sindhikara sindhikara@gmail.com
    Does PDB related stuff
    Todo: 
1) Integrate better with AMBER modules (e.g. guessatom, atom radii, etc)



'''
import math

class Atom():
    '''
Atom object
Stores atomic information as would be in a pdb line.

Below is PDB format
1-6 7-11 13-16 17 18-20 22 23-26 27 31-38 39-46 47-54 55-60 61-66 77-78 79-80
Record name Integer Atom Character Residue name Character Integer AChar Real(8.3) Real(8.3) Real(8.3) Real(6.2) Real(6.2) LString(2) element LString(2) charge
"ATOM " serial name altLoc resName chainID resSeq iCode
Atom serial number. Atom name. Alternate location indicator. Residue name. Chain identifier. Residue sequence number. Code for insertion of residues. Orthogonal coordinates for X in Angstroms. Orthogonal coordinates for Y in Angstroms. Orthogonal coordinates for Z in Angstroms. Occupancy. Temperature factor. Element symbol, right-justified. Charge on the atom.
    '''
    def __init__(self,serial=1,name="D",altloc=" ",resname="BCD",chainid=1,resseq=1,icode=" ",coord=[0.0,0.0,0.0],occ=1.0,tfac=0.0,element="D",charge=0.0,hetatm=False,ister=False):
        self.name = name[0:4]
        self.serial=serial
        self.altloc=altloc[0]
        self.resname = resname[0:3]
        self.chainid=str(chainid)[0]
        self.resseq=resseq
        self.icode=icode
        self.coord=coord
        self.occ=occ
        self.tfac=tfac
        self.element=element[0]
        self.charge=charge
        self.hetatm=hetatm
        self.ister=ister # this is the last Atom of the unit. (the next line in the pdb will be say "TER"
    def __str__(self):
        '''
    Returns a PDB format string
        '''
        if not self.hetatm :
            pdbstring="ATOM  %5d %4s%1s%3s %1c%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f           %1s%2d\n" % (self.serial,self.name,self.altloc,self.resname,self.chainid,self.resseq,self.icode,self.coord[0],self.coord[1],self.coord[2],self.occ,self.tfac,self.element,self.charge)
        else :
            pdbstring="HETATM%5d %4s%1s%3s %1c%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f           %1s%2d\n" % (self.serial,self.name,self.altloc,self.resname,self.chainid,self.resseq,self.icode,self.coord[0],self.coord[1],self.coord[2],self.occ,self.tfac,self.element,self.charge)
        if self.ister :
            pdbstring="%sTER\n" % pdbstring
        return(pdbstring) 

    def makestring(self):
        '''
    Antiquated function same as __str__
        '''
        print "# Warning, antiquated function, __str__ function added"
        return(str(self))

        return(pdbstring)
        
def centermolecule(listofAtoms):
    '''
Given a list of Atom class objects, returns the same list centered
    '''
    center=[0.0,0.0,0.0]
    for atom in listofAtoms:
        for dim in range(3):
            center[dim]+=atom.coord[dim]
    for dim in range(3):
        center[dim]=center[dim]/float(len(listofAtoms))
    for i in range(len(listofAtoms)):
        for dim in range(3):
            listofAtoms[i].coord[dim]-=center[dim]
    return(listofAtoms)


def readcoordinates(pdbfilename):
    ''' Given a pdbfilename, return list of coordinates (2dlist)
This function is less functional than readpdbintoAtoms()
It only stores coordinates
    '''
    pdbfile=open(pdbfilename)
    coords=[]
    for line in pdbfile:
        ident=line[0:6]
        if ident == "ATOM  " or ident == "HETATM" :
            coords.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
    return(coords)

def readpdbintoAtoms(pdblines,guesselement=False,guesshetatm=False):
    ''' Given a pdb file of arbitrary length parsed into lines,
store the data as a list of "Atom" class objects
    '''

    molecule=[]
    knownelements = ["C","N","O","P","H","S","K","KA"] # put larger ones at the end (so that KA supersedes K)
    knownresidues = ["THR","ILE","ALA","ARG","HIS","HIE","HID","CYS","CYX","LYS","ASP","ASH","GLU","SER","THR","ASN","THR","ASN","GLN","GLY","PRO","LEU","MET","PHE","TRP","TYR","VAL"]
    # The above lists are used to help classify atoms
    if guesselement :
        pass
        #print "Warning! Guesselement only uses single letter atoms"
    terminilist=[]
    for i,line in enumerate(pdblines):
        ident=line[0:6]
        if ident == "ATOM  " or ident == "HETATM" :
            #index in line is one less than column number
            serial=int(line[6:11])
            name=line[12:16]
            altloc=line[16]
            resname=line[17:20]
            chainid=line[21]
            resseq=int(line[22:26])
            icode=line[26]
            coord=[float(line[30:38]),float(line[38:46]),float(line[46:54])]
            occ=float(line[54:60])
            tfac=float(line[60:66])
            if len(line) >= 77 :
                element=line[76:78]
            else:
                element="  "
            if guesselement :
                for knownelement in knownelements:
                    if knownelement in name:
                        element=knownelement
            if len(line) >= 79:
                try :
                    charge=float(line[78:80])
                except :
                    charge=0.0
            else:
                charge=0.0
            hetatm=False
            if ident == "HETATM" :
                hetatm=True
            elif (resname not in knownresidues) and guesshetatm:
                hetatm=True
            molecule.append(Atom(serial,name,altloc,resname,chainid,resseq,icode,coord,occ,tfac,element,charge,hetatm))
        if line[0:3] == "TER" :
            molecule[-1].ister=True
    return(molecule)

def isitinside(coordinate,atoms,buffer=0.0):
    '''
Test if coordinate is inside VDW radius of list of atoms
    '''
    radii={'H' : 0.6, 'O' : 1.6612, 'C' : 1.908, 'S': 2.000, 'N' : 1.8240, 'KA' : 1.7131,  'C0' : 1.7131, 'P' : 2.1000, 'K' : 1.7131}
    inside = False
    for atom in atoms:
        if not atom.ister:
            dist=(coordinate[0]-atom.coord[0])**2+(coordinate[1]-atom.coord[1])**2+(coordinate[2]-atom.coord[2])**2
            dist=math.sqrt(dist)
        
            try:
                 if dist < radii[atom.element[0]]+buffer :
                     inside = True
                     break
            except:
                print "Cannot find radius of atom named: ",atom.element,atom.name
                exit()
    return(inside)

def isitnear(coordinate,atoms,radius):
    '''
Test if coordinate is within specified radius of listed atoms
    '''
    inside = False
    for atom in atoms:
        dist=(coordinate[0]-atom.coord[0])**2+(coordinate[1]-atom.coord[1])**2+(coordinate[2]-atom.coord[2])**2
        dist=math.sqrt(dist)
        if dist < radius:
            inside = True
            break
    return(inside)


def calcmolvolume(atoms,gridspacing,probesize=0.0):
    import modules.djsgrid as grid
    print "This function is under testing and quite slow. This should be replaced by faster, more reliable code"

    volume=0.0
    print "spacing =",gridspacing
    gridvolume=gridspacing**3
    print "grid volume is ",gridvolume
    maxes=[-1000]*3
    mins=[1000]*3
    for atom in atoms:
        if atom.coord[0] < mins[0] :
            mins[0] = atom.coord[0]
        elif atom.coord[0] > maxes[0] :
            maxes[0] = atom.coord[0]
        if atom.coord[1] < mins[1] :
            mins[1] = atom.coord[1]
        elif atom.coord[1] > maxes[1] :
            maxes[1] = atom.coord[1]
        if atom.coord[2] < mins[2] :
            mins[2] = atom.coord[2]
        elif atom.coord[2] > maxes[2] :
            maxes[2] = atom.coord[2]
    # add some buffer
    for i in range(3):
        mins[i]-=5.0
        maxes[i]+=5.0

    spans=[]
    for i in range(3):
        spans.append(int((maxes[i]-mins[i])/gridspacing))

    radmap={'H' : 0.6, 'O' : 1.6612, 'C' : 1.908, 'S': 2.000, 'N' : 1.8240, 'KA' : 1.7131,  'C0' : 1.7131, 'P' : 2.1000, 'K' : 1.7131}
    shelldict={} # construct dictionary of shell indices (these contain all grid points radially from atom center which count as inside the atom
    print " don't need to compute shellindices more than once, do it for a big one and assign shells"
    allshellindices=grid.precomputeshellindices(int((2.5+probesize)/gridspacing))
    for radkey in radmap: 
        myrad=radmap[radkey]+probesize # effective radius
        maxshellindex=int(myrad/gridspacing+0.5) # largest shell index that counts as inside
        shelldict[radkey]=allshellindices[:maxshellindex]

    occupiedflags=[[[False for i in range(spans[0])] for j in range(spans[1])] for k in range(spans[2])]
    testvol=0
    for shellindex in shelldict["C"]:
        for shellpoint in shellindex:
            testvol+=gridvolume
    print "at this resolution, Carbon volume (plus probe) is",testvol
    print "by geometry it's", (4.0/3.0)*3.14159 * (radmap["C"]+probesize) **3

    occupiedpdbfile=open("dumb.pdb","w")
    testedpdbfile=open("tested.pdb","w")
    for atom in atoms:
        myindices=grid.getindicesfromcoord([atom.coord[0],atom.coord[1],atom.coord[2]],mins,[gridspacing]*3)
        #print atom.coord[0],atom.coord[1],atom.coord[2]
        #print grid.getcoordfromindices(myindices,mins,[gridspacing]*3)
        myshells=shelldict[atom.element]
        for shell in myshells: # 0th, 1st, 2nd shell, etc
            for shellindex in shell: # points that lay on a shell
                netindex=[]
                for i in range(3): # indices of each point on shell
                    #print "shellindex,myindex[i]=",shellindex[i],myindices[i]
                    netindex.append(shellindex[i]+myindices[i]) # origin of atom + shell index
                    #print spans
                #print "netindices=",netindex
                #print grid.getcoordfromindices(netindex,mins,[gridspacing]*3)
                coord=grid.getcoordfromindices(netindex,mins,[gridspacing]*3)
                testedpdbfile.write(Atom(coord=coord).makestring())
                if not occupiedflags[netindex[0]][netindex[1]][netindex[2]]:
                    volume+=gridvolume
                    occupiedflags[netindex[0]][netindex[1]][netindex[2]]=True
                    occupiedpdbfile.write(Atom(coord=coord).makestring())
    occupiedpdbfile.close()
    testedpdbfile.close()
    return(volume)




