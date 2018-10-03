
from rdkit import Chem
import sys
import numpy as np
import random
from copy import deepcopy
from rdkithelpers import *
from rdkit.Chem import Descriptors

np.random.seed(2)
random.seed(2)

ringCrossover=False
MAXTRY = 100

maxWeight = 0
MxAtm = 0

debug=False
_ringfilter = lambda x:len(x) in [5,6]


def Init():
    global MxAtm
    if MxAtm==0:
        MxAtm=sys.maxint

###############################################################################
#                            CrossOver Methods                                #
###############################################################################


def Crossover(m1, m2):
    #Kekulize mols:
    for m in (m1, m2):
        Chem.Kekulize(m, True)

    # this is a boolean switch. we can change this later:
    # 1. we can first try crossover and if it did not succeed
    # try a normal crossover
    # 2. we can assign a propability to RingCrossover as with the mutation
    # interface
    if ringCrossover:
        return RingCrossover(m1, m2)
    else:
        return NormalCrossover(m1, m2)


def NormalCrossover(m1, m2):
    #Fragment molecules
    m1fs = GetFragment(m1)
    m2fs = GetFragment(m2)

    #pick fragments for crossover checking:
    # 1. Molecular weight
    Choices = []
    #maxWeight = mprms.maxWeight
    if maxWeight > 0:
        w1 = [Descriptors.MolWt(f) for f in m1fs]
        w2 = [Descriptors.MolWt(f) for f in m2fs]
        for i, j in ((i, j) for i in xrange(2) for j in xrange(2)):
            if w1[i] + w2[j] < maxWeight + 50.0: Choices.append((i, j))
    # 2. Number of Atoms
    #MxAtm=mprms.MxAtm
    if MxAtm > 0:
        a1 = [f.GetNumAtoms() for f in m1fs]
        a2 = [f.GetNumAtoms() for f in m2fs]
        for i, j in ((i, j) for i in xrange(2) for j in xrange(2)):
            if a1[i] + a2[j] < MxAtm + 4: Choices.append((i, j))

    if len(Choices) == 0:
        raise MutateFail()
    else:
        choice = random.choice(Choices)
        mol1 = m1fs[choice[0]]
        mol2 = m2fs[choice[1]]

    # now mol2 has to be connected to mol1
    # but what happens here?
    mol1.SetProp('parent1', m1.GetProp('isosmi'))
    mol2.SetProp('parent2', m2.GetProp('isosmi'))

    #Append molecule 2 to molecule 1.
    newmol = Chem.CombineMols(mol1, mol2)
    newids = Chem.GetMolFrags(newmol)

    # Now, bond the two fragments to create the child molecule.
    # We choose a random pair of possible atoms to do the attachment;
    possibleA1 = [
        atom for atom in GetIAtoms(newids[0], newmol) if EmptyValence(atom) > 0
    ]
    possibleA2 = [
        atom for atom in GetIAtoms(newids[1], newmol) if EmptyValence(atom) > 0
    ]

    if len(possibleA1) == 0 or len(possibleA2) == 0:
        #print "no possible atoms!"
        raise MutateFail()
    newmolRW = Chem.RWMol(newmol)
    atom1 = random.choice(possibleA1)
    atom2 = random.choice(possibleA2)
    newmolRW.AddBond(atom1.GetIdx(), atom2.GetIdx(), Chem.BondType.SINGLE)

    #print "new mol!", Chem.MolToSmiles(newmolRW)
    mol = newmolRW.GetMol()
    return mol


def GetFragment(mol):
    #print "in get fragment with:", Chem.MolToSmiles(mol)
    tried = []
    for itry in xrange(MAXTRY):
        frag = Chem.RWMol(mol)
        #frag.ClearComputedProps()
        ResetProps(frag)
        #Cut an acyclic bond to break molecule in 2
        # get the ring bond ids. this gives a tuple of tuple for every ring
        ringbondids = frag.GetRingInfo().BondRings()
        # flatten this tuple to a 1D set
        ringbondids = list(
            set([bond for subset in ringbondids for bond in subset]))
        # get all the bonds that are NOT ringbonds or already tried
        bonds = [
            bond for bond in mol.GetBonds()
            if bond.GetIdx() not in tried + ringbondids
        ]
        if len(bonds) == 0: break
        # choose a bond
        delbond = random.choice(bonds)
        tried.append(delbond.GetIdx())
        # get atom indices of the bond
        i, j = (delbond.GetBeginAtomIdx(), delbond.GetEndAtomIdx())

        frag.RemoveBond(i, j)
        # get new fragments as two new normal Mol objects
        try:
            frags = Chem.GetMolFrags(frag, asMols=True)
        except ValueError:
            print "frag:", Chem.MolToSmiles(frag)
            print "mol:", Chem.MolToSmiles(mol)
            break
        molfrags = [Chem.MolToSmiles(x, True) for x in frags]
        if len(frags) < 2:
            #print "only one frag"
            continue

        #Check to make sure fragments have at least 2 atoms
        if any(frag.GetNumAtoms() < 2 for frag in frags):
            #print "not enough atoms on frag"
            continue
        #print  molfrags
        # convert to RWMols?
        frags = map(Chem.RWMol, frags)
        return frags
    # If here, we failed to find a good bond to cut
    #print "no nice bond to cut"
    raise MutateFail()






############### RING CROSSOVERS: #############################


totalringcuts = {
        6:(
            ((0,1), (3,4)),
            ((1,2), (4,5)),
            ((2,3), (0,5))
          ),
        5:(
            ((0,1), (3,4)),
            ((0,1), (2,3)),
            ((1,2), (3,4)),
            ((1,2), (0,4)),
            ((2,3), (0,4))
          )
        }
if True:
    cuts24 = (
            ((0,1),(2,3)),
            ((0,1),(4,5)),
            ((2,3),(4,5)),
            ((1,2),(0,5)),
            ((1,2),(3,4)),
            ((0,5),(3,4))
            )
    totalringcuts[6] = totalringcuts[6] + cuts24
if True:
    cuts7 = (
            ((0,1),(3,4)),
            ((0,1),(4,5)),
            ((1,2),(4,5)),
            ((1,2),(5,6)),
            ((2,3),(5,6)),
            ((2,3),(6,0)),
            ((3,4),(6,0))
            )
    totalringcuts[7]=cuts7


def RingCrossover(m1, m2):
    for m in (m1, m2):
        for atom in m.GetAtoms():
            atom.ClearProp('type1')
            atom.ClearProp('type2')

    # 1. cut a ring from mol1 and choose a fragment
    frags1, ftypes = cutring1(Chem.RWMol(m1))

    # 2. try to cut a ring from mol2 with similar bondorders
    frags2 = cutring2(Chem.RWMol(m2), ftypes)

    if debug: print "frags1/frags2:", frags1, frags2
    
    # 3. combined molecule should behave two criteria
    Choices = []
    # 3.1 
    if maxWeight > 0:
        w1 = [Descriptors.MolWt(f) for f in frags1]
        w2 = [Descriptors.MolWt(f) for f in frags2]
        for i, j in ((i, j) for i in xrange(len(w1)) for j in xrange(len(w2))):
            if w1[i] + w2[j] < maxWeight + 50.0: Choices.append((i, j))
    # 3.2 Number of Atoms
    a1 = [f.GetNumAtoms() for f in frags1]
    a2 = [f.GetNumAtoms() for f in frags2]
    for i, j in ((i, j) for i in xrange(len(a1)) for j in xrange(len(a2))):
        if 3 < a1[i] + a2[j] < MxAtm + 4: Choices.append((i, j))

    if len(Choices) == 0:
        raise MutateFail()
    else:
        choice = random.choice(Choices)
        mol1 = frags1[choice[0]]
        mol2 = frags2[choice[1]]

    # 4. Set Parents
    mol1.SetProp('parent1', m1.GetProp('isosmi'))
    mol2.SetProp('parent2', m2.GetProp('isosmi'))

    # 5. combine the fragments to a new molecule
    newmol = combinemols(mol1, mol2, ftypes)
    return newmol

def cutring1(originalmol):
    for atom in originalmol.GetAtoms():
        atom.ClearProp('type1')
    rings = originalmol.GetRingInfo().AtomRings()
    rings = filter(_ringfilter, rings)
    if len(rings)==0:
        raise MutateFail()
    elif len(rings)==1:
        ringtocut = rings[0]
    else:
        try:
            ringtocut = random.choice(rings)
        except ValueError:
            if debug: print rings
            raise
    #print "rings 1:", rings
    # get ways to cut the rings
    ringcuts = totalringcuts[len(ringtocut)]
    nringcuts= len(ringcuts)

    # and choose such a way
    for imanner in np.random.permutation(nringcuts):
        mol = deepcopy(originalmol)
        manner = ringcuts[imanner]
        if debug: print manner

        # and cut mol in these bonds:
        ftypes = []
        for bondtocut in manner:
            i, j = ringtocut[bondtocut[0]], ringtocut[bondtocut[1]]
            ftype= mol.GetBondBetweenAtoms(i, j).GetBondType()
            ftypes.append( ftype )
            for x in (i,j):
                xatom = mol.GetAtomWithIdx(x)
                xatom.SetProp('type1', str(ftype))
            mol.RemoveBond(i, j)
         
        # try to get the two parts
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        if debug: print "molfrags:", [Chem.MolToSmiles(x, True) for x in frags]
        
        if len(frags)==2:
            break
        else:
            if debug: print "not two fragments"
            for frag in frags:
                for atom in frag.GetAtoms():
                    atom.ClearProp('type1')
            for atom in mol.GetAtoms():
                atom.ClearProp('type1')

    if frags is None:
        raise MutateFail()

    return frags, ftypes


def cutring2(originalmol, ftypes):
    frags=None
    for atom in originalmol.GetAtoms():
        atom.ClearProp('type2')
        atom.ClearProp('type1')
    rings = originalmol.GetRingInfo().AtomRings()
    rings = filter(_ringfilter, rings)
    if len(rings)==0:
        raise MutateFail()
    # make rings in random order
    ringindices = np.random.permutation(len(rings))

    # for every ring try to connect it to one of the frags:
    for rindex in ringindices:
        mol = deepcopy(originalmol)
        ring = rings[rindex]
        # for every manner to cut the ring
        ringcuts= totalringcuts[len(ring)]
        nringcuts=len(ringcuts)
        for imanner in np.random.permutation(nringcuts):
            manner = ringcuts[imanner]
            # print manner
            # check if bondtypes match
            # get bondindices
            newftypes = []
            for bondtocut in manner:
                i, j = ring[bondtocut[0]], ring[bondtocut[1]]
                if debug: print "manner:", manner, "bondtocut:", bondtocut, "ring:", ring, "i,j", i,j
                try:
                    btype = mol.GetBondBetweenAtoms(i, j).GetBondType()
                except AttributeError:
                    raise

                newftypes.append( btype )
            if debug: print "newftype:", newftypes, "ftype:", ftypes
            if not set(newftypes)==set(ftypes):
                if debug: print "not similar"
                continue
            else:
                if debug: print "similar!"

            # now if similar cut it and check mol falls in two fragments:
            for bondtocut in manner:
                i, j = ring[bondtocut[0]], ring[bondtocut[1]]
                btype = mol.GetBondBetweenAtoms(i, j).GetBondType()
                for x in (i,j):
                    xatom = mol.GetAtomWithIdx(x)
                    xatom.SetProp('type2', str(btype))
                mol.RemoveBond(i, j)

            # try to get the two parts
            frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
            if debug: print "molfrags 2:", [Chem.MolToSmiles(x, True) for x in frags]
            
            if not len(frags)==2:
                if debug: print "not two fragments"
                # since we cut bonds we go back to original molecule
                mol = deepcopy(originalmol)
                for atom in mol.GetAtoms():
                    atom.ClearProp('type2')
                continue
            else:
                if debug: print "two fragments"
                # now choose one of the two new fragments
                return frags
    
    if frags is None or len(frags)==1:
        raise MutateFail()
    if debug: print "frags 2:", frags
    return


def combinemols(frag1, frag2, ftypes):
    mol = Chem.RWMol(Chem.CombineMols(frag1, frag2))
    #print Chem.MolToSmiles(mol)

    # decide which atoms have to be bonded to what
    props = [ 'type1', 'type2' ]
    bonds = [{'type':None, 'indices':[], 'already':None}, {'type':None, 'indices':[]} ]
    for atom in mol.GetAtoms():
        for prop in props:
            if atom.HasProp(prop):
                atype = atom.GetProp(prop)
                if bonds[0]['type'] is None or bonds[0]['type']==atype and not (
                        len(bonds[0]['indices'])==2 or bonds[0]['already']==prop ):
                    bonds[0]['type']=atype
                    bonds[0]['indices'].append(atom.GetIdx())
                    bonds[0]['already']=prop
                else:
                    bonds[1]['type']=atype
                    bonds[1]['indices'].append(atom.GetIdx())
                break

    for bond in bonds:
        try:
            a1 = mol.GetAtomWithIdx(bond['indices'][0])
            a2 = mol.GetAtomWithIdx(bond['indices'][1])
        except IndexError:
            if debug:
                print bonds
                for atom in mol.GetAtoms():
                    print atom.GetPropsAsDict()
            raise MutateFail
        #print a1, a2
        if bond['type']=='SINGLE':
            bo = Chem.BondType.SINGLE
        else:
            if debug: print bond['type']
            bo = Chem.BondType.DOUBLE
        #print bo
        try:
            mol.AddBond(a1.GetIdx(), a2.GetIdx(), bo)
        except RuntimeError:
            if debug:
                print a1.GetIdx(), a2.GetIdx(), bo, bonds
                for atom in mol.GetAtoms():
                    print atom.GetPropsAsDict()
            raise MutateFail()

    if len(Chem.GetMolFrags(mol))==1:
        if debug: print "succesfull ringcrossover:", Chem.MolToSmiles(mol)
    else:
        raise MutateFail('not a single molecule')

    return mol


