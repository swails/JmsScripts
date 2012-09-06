#!/usr/bin/env python
'''Created by Daniel Sindhikara, sindhikara@gmail.com
Program converts an AMBER prmtop and PDB to an AMBER mdl file

'''
import sys

def main():
    from chemistry.amber.amberformat import AmberFormat
    from chemistry.amber.readparm import AmberParm
    #need AmberParm because AmberFormat can't load pointers/LJ?
    import chemistry.pdb as pdb 
    minargs=3
    numargs=len(sys.argv)
    if(numargs<minargs) :
        print "# format prmtop_and_pdb_2_mdl.py <prmtop file> <pdbfile>"
    	print "Insufficient arguments, need ",minargs," : ",numargs
    	exit()
    print '''
Warning! This utility gives a mdl file based on the prmtop.
But the results should be inspected afterwards.
Some things that work in MD do not work in RISM.
Consult relevant literature for more information.
    '''

    prmfilename=sys.argv[1]
    pdbfilename=sys.argv[2]
    prmtopdata=AmberParm(prmfilename)
    mdldata=AmberFormat() 
    #TITLE -> TITLE 20a4
    #POINTERS [0:2] -> POINTERS 10i8
    #set(ATOM_TYPE_INDEX) -> ATMTYP 10i8
    #ATOM_NAME->ATMNAME 20a4
    #MASS->MASS 5e16.8
    #CHARGE->CHG 5e16.8
    #unpaired LJ parameters -> LJEPSILON , -> LJSIGMA 5e16.8
    #multiplicity -> MULTI 10I8
    #centered coordinates COORD 5e16.7
    numatoms,numtypes=prmtopdata.parm_data["POINTERS"][0:2]
    mdldata.addFlag("TITLE","20a4",data=prmtopdata.parm_data["TITLE"])
    mdldata.addFlag("POINTERS","10i8",data=prmtopdata.parm_data["POINTERS"][0:2])
    mdldata.addFlag("ATMTYP","10i8",data=set(prmtopdata.parm_data["ATOM_TYPE_INDEX"]))

    atomtypeindices=prmtopdata.parm_data["ATOM_TYPE_INDEX"]
    dupenames=prmtopdata.parm_data["ATOM_NAME"]
    names = [dupenames[atomtypeindices.index(i+1)] for i in range(numtypes)]
    mdldata.addFlag("ATMNAME","20a4",data=names)

    dupemasses=prmtopdata.parm_data["MASS"]
    masses = [dupemasses[atomtypeindices.index(i+1)] for i in range(numtypes)]
    mdldata.addFlag("MASS","5e16.8",data=masses)

    dupecharges=prmtopdata.parm_data["CHARGE"] # need to add back in multiplicative factor
    charges=[dupecharges[atomtypeindices.index(i+1)]*18.2223 for i in range(numtypes)]
    mdldata.addFlag("CHG","5e16.8",data=charges)

    #Nonbonded
    prmtopdata.LoadPointers()
    prmtopdata.fill_LJ()
    mdldata.addFlag("LJEPSILON","5e16.8",data=prmtopdata.LJ_depth)
    mdldata.addFlag("LJSIGMA","5e16.8",data=prmtopdata.LJ_radius)
    multiplicity=[atomtypeindices.count(i+1) for i in range(numtypes)]
    mdldata.addFlag("MULTI","10i8",data=multiplicity)

    #coordinates
    mol=pdb.readpdbintoAtoms(open(pdbfilename).readlines())
    mol=pdb.centermolecule(mol)
    coords=[]
    for atom in mol:
        for dim in range(3):
            coords.append(atom.coord[dim])

    mdldata.addFlag("COORD","5e16.8",data=coords)
    mdlfilename=pdbfilename.split(".pdb")[0]+".mdl" 
    mdldata.writeParm(mdlfilename, overwrite=False)
    print "Finished constructing:",mdlfilename

if __name__ == '__main__' :
    main()
