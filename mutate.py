#!/usr/bin/env python
#-*- coding: utf-8 -*-

import random
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkithelpers import *
from helpers import *
from output import stats
import molfails
from molfails import MutateFail
import ACSESS
import mprms


MAXTRY = 100
elements = []
maxWeight = 0
MxAtm = 0

###### mutation probabilities:
mutationtypes = ['FlipAtom', 'FlipBond', 'AddAtom', 'AddBond', 'DelAtom', 'DelBond', 'AddAroRing', 'AddFusionRing', 'CustomMutator']
mutators = dict()

p_FlipAtom = 0.8
p_FlipBond= 0.8
p_AddAtom= 0.5  #Actual probability: (.8*.7)*.5=.28
p_AddBond= 0.1
p_DelAtom = 0.8  #Actual probability: .224
p_DelBond = 0.2  #actual probability=.8*.3=.24
p_AddAroRing = 0.5  #actual probability is rather low sing free double bonds aren't prevalent
p_AddFusionRing = 0.5  #actual probability is rather low sing free double bonds aren't prevalent
p_CustomMutator = 0.0

MutateStereo = False
StereoFlip = 0.2

debug = False

###############################################
# Initialization the module
###############################################


def Init():
    #elements = mprms.elements
    global elements, halogens, mutators
    halogens = [9, 17, 35, 53]

    #do intersection and difference
    halogens = set(elements) & set(halogens)
    elements = set(elements) - set(halogens)
    #change from set to list
    elements = list(elements)
    halogens = list(halogens)
    if halogens:
        print "elements:", elements
        print "halogens:", halogens

    # this is set here. If not given in input 
    # it is deleted in next for-loop!
    try:
        CustomMutator.setnname('nCustom')
        mutators['CustomMutator'] = CustomMutator
    except AttributeError:
        pass

    for mutationtype in mutationtypes:
        p = globals()["p_{}".format(mutationtype)]
        if p==0.0:
            try:
                del mutators[mutationtype]
                print mutationtype, "deleted"
            except KeyError:
                pass
        else:
            mutators[mutationtype].p = p
            print "mutation {:13s} set with probability {:3.2f}".format(mutationtype, p)

    return


###############################################################################
#                            Mutation Methods                                 #
###############################################################################
#     Nota Bene! These mol objects are now based on Chem.RWMol objects!       #
###############################################################################

# Mutation driver
def MakeMutations(candidate):
    ''' The mutations, based on not-aromatic SMILES'''
    # 1. Kekulize:
    try:
        Chem.Kekulize(candidate, True)
    except ValueError:
        print "MakeMutation. Kekulize Error:", Chem.MolToSmiles(candidate)
        raise MutateFail(candidate)
    # 2. Mutate
    candidate = SingleMutate(candidate)
    # 3. SetAromaticity again
    try:
        if aromaticity: Chem.SetAromaticity(candidate)
        candidate = Finalize(candidate, tautomerize=False)
    except ValueError:
        print "MakeMutation. SetAromaticity Error:", Chem.MolToSmiles(
            candidate)
        raise MutateFail(candidate)
    return candidate

def SingleMutate(candidateraw):
    if candidateraw is None: raise ValueError('candidate is none')
    else: candidate = Chem.RWMol(candidateraw)
    global stats

    parent = candidate.GetProp('isosmi')
    ResetProps(candidate)
    change = False
    candidate.SetProp('parent', parent)

    for mutationtype, mutator in mutators.items():
        if random.random() < mutator.p:
            #print "in mutate:",  mutationtype, Chem.MolToSmiles(candidate),
            change = True
            stats[mutator.nname] += 1
            try:
                candidate = mutator(candidate)
            # I don't understand why these errors are not similar?! but this works
            except MutateFail as e:
                stats[mutator.nnameFail] += 1
            except Exception as e:
                if 'MutateFail' in repr(e):
                    #print "exotic MutateFail"
                    stats[mutator.nnameFail] += 1
                else:
                    print "Exception is not MutateFail!"
                    print e, repr(e), type(e)
                    raise

            # should we allow multiple mutation at the same time:
            # if not:
            # break

    if not change:
        stats['nNoMutation'] += 1
        raise MutateFail()
    return candidate


class Mutator(object):
    def __init__(self, mutator, p=0.0, nname=None):
        if not callable(mutator):
            raise TypeError('mutator is not callabe')
        self.mutator = mutator
        self.name = mutator.__name__
        self.setnname(nname)
        self.p = p
        return

    def setnname(self, nname=None):
        if nname is None:
            self.nname = "n{}".format(self.name)
        else:
            self.nname = nname
        self.nnameFail = "{}Fail".format(self.nname)

    def __call__(self, mol):
        Chem.Kekulize(mol, True)
        try:
            mol = self.mutator(mol)
        except MutateFail as e:
            raise
        except RuntimeError as e:
            print "RuntimeError in Mutator:", self.mutator, self.name
            print "repr(e):", repr(e)
            raise MutateFail()
        except Exception as e:
            if 'MutateFail' in repr(e):
                raise MutateFail()
            else:
                raise

        if mol is None:
            raise MutateFail()

        mol = Finalize(mol, tautomerize=False, aromatic=False)

        if type(mol)==Chem.Mol:
            mol = Chem.RWMol(mol)
        return mol

    def __str__(self):
        return "Mutator: {}".format(self.name)

    def __repr__(self):
        return self.__str__()

def FlipBond(mol):
    bonds = list(GetBonds(mol, notprop='group'))
    if len(bonds)==0: raise MutateFail()
    bond = random.choice(bonds)
    oldorder = int(bond.GetBondTypeAsDouble())
    atom1 = bond.GetBeginAtom()
    atom2 = bond.GetEndAtom()
    full = False
    if EmptyValence(atom1) == 0: full = True
    if EmptyValence(atom2) == 0: full = True
    if oldorder == 1 and full:
        raise MutateFail()
    elif oldorder == 1:
        add = 1
    elif oldorder == 3:
        add = -1
    elif oldorder == 2 and full:
        add = -1
    else:
        add = random.choice((-1, 1))
    bond.SetBondType(bondorder[oldorder + add])
    return mol
mutators['FlipBond'] = Mutator(FlipBond)

def FlipAtom(mol):
    atoms = filter(CanChangeAtom, mol.GetAtoms())
    if len(atoms)==0: raise MutateFail()
    atom = random.choice(atoms)

    changed = False
    neighbors=list(atom.GetNeighbors())
    valence = atom.GetExplicitValence() + atom.GetNumRadicalElectrons()

    # Add a halogen, if possible
    # we only add halogens when neighbor and neighbor-neighbor only C
    if valence==1:
        if (len(halogens)>0
        and random.random()>0.6
        and neighbors[0].GetAtomicNum()==6):
            CanAddHalogen=True
            nextneighbors=neighbors[0].GetNeighbors()
            for nb in nextneighbors:
                if (  nb.GetAtomicNum() != 6
                      and nb.GetIdx() != atom.GetIdx()):
                    CanAddHalogen=False
                    break
            if CanAddHalogen:
                atom.SetAtomicNum( random.choice(halogens))
                return mol

    # # # # # # # # #
    # Regular atom switching
    atnum = atom.GetAtomicNum()
    elems = [el for el in elements if el != atnum]

    for itry in range(50):
        cand = random.choice(elems)
        if MaxValence[cand] >= valence:
            atom.SetAtomicNum(cand)
            return mol

    if not changed:
        raise MutateFail()  # if here, mutation failed ...
    return mol
mutators['FlipAtom'] = Mutator(FlipAtom)

def AddBond(mol):
    natoms = mol.GetNumAtoms()
    if natoms < 5: raise MutateFail()  #don't create 4-member rings

    #Don't create non-planar graphs (note that this does not prevent
    #all non-planar graphs, but it does prevent some)
    if natoms >= 3 and mol.GetNumBonds() >= 3 * natoms - 6:
        raise MutateFail()

    # List of atoms that can form new bonds
    atoms = [
        atom for atom in GetAtoms(mol, notprop='protected')
        if EmptyValence(atom) > 0
    ]
    if len(atoms) < 2: raise MutateFail()

    for i in range(MAXTRY):
        atom1 = random.choice(atoms)
        atom2 = random.choice(atoms)
        if atom1 == atom2: continue
        # new rather strict criteria;
        if atom1.IsInRing() and atom2.IsInRing():
            #print "continue because both atoms in ring"
            continue
        # Only make rings of size 5-7
        if not (5 <= len(
                Chem.GetShortestPath(mol, atom1.GetIdx(), atom2.GetIdx())) <=
                7):
            continue
        nfreebonds1 = EmptyValence(atom1)
        nfreebonds2 = EmptyValence(atom2)
        if nfreebonds2 < nfreebonds1: nfreebonds1 = nfreebonds2
        assert nfreebonds1 > 0
        order = random.randrange(1, nfreebonds1 + 1)
        mol.AddBond(atom1.GetIdx(), atom2.GetIdx(), bondorder[order])
        return mol
    # if here, mutation failed
    raise MutateFail()
mutators['AddBond']=Mutator(AddBond)

def DelBond(mol):
    bondringids = mol.GetRingInfo().BondRings()

    # need to flatten bondringids:
    if len(bondringids)>0:
        # fancy manner to flatten a list of tuples with possible duplicates
        #bondringids = set.intersection(*map(set, bondringids))
        bondringids = set.union(*map(set, bondringids))
        bonds = GetIBonds(bondringids, mol, notprop='group')
    else:
        raise MutateFail()
    try:
        bond = random.choice(bonds)
    except IndexError:
        print bondringids
    mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    #print "succesful delbond:", Chem.MolToSmiles(mol)
    return mol
mutators['DelBond'] = Mutator(DelBond)

def AddAtom(mol):
    # 5. add an atom
    atoms = GetAtoms(mol, notprop='protected')
    bonds = GetBonds(mol, notprop='group')
    if len(atoms)<1: raise MutateFail()
    obj = random.choice(atoms + bonds)

    # If passed a bond, insert an atom into it
    if type(obj) == Chem.Bond:
        if obj.HasProp('group'):
            raise MutateFail(mol, 'Grouped bond passed to AddAtom')
        # note that we really nead the atoms here because the Idx may change
        # after adding a new atom to the molecule
        atom1 = obj.GetBeginAtom()
        atom2 = obj.GetEndAtom()
        mol.RemoveBond(atom1.GetIdx(), atom2.GetIdx())
        newatomid = mol.AddAtom(Chem.Atom(6))
        mol.AddBond(atom1.GetIdx(), newatomid, bondorder[1])
        mol.AddBond(atom2.GetIdx(), newatomid, bondorder[1])

    #Otherwise, make a new sidechain (or add a terminal atom)
    elif type(obj) == Chem.Atom and EmptyValence(obj) > 0:
        if obj.HasProp('protected'):
            raise MutateFatal(mol, 'Trying to add atom to protected atom.')

        # choose an atom to add
        atnum = obj.GetAtomicNum()
        elems = [el for el in elements if el != atnum]
        if atnum==6: # only add halogens to carbons
            elems.extend(halogens)
        cand = random.choice(elems)

        # and add it
        newatomid = mol.AddAtom(Chem.Atom(cand))
        mol.AddBond(obj.GetIdx(), newatomid, bondorder[1])
    else:
        raise MutateFail()
    return mol
mutators['AddAtom'] = Mutator(AddAtom)


def DelAtom(mol):
    inismi = Chem.MolToSmiles(mol)
    atoms = filter(CanRemoveAtom, mol.GetAtoms())
    if len(atoms)==1:
        raise MutateFail()

    atom = random.choice(atoms)

    # or GetTotalDegree?
    #Degree = atom.GetDegree() + atom.GetNumRadicalElectrons()
    Degree = atom.GetDegree()

    #If there's only one atom, refuse to remove it
    if len(atom.GetNeighbors()) == 0:
        raise MutateFail(mol)

    #Remove a terminal atom:
    elif Degree == 1:
        for bond in atom.GetBonds():
            mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        mol.RemoveAtom(atom.GetIdx())

    # until here it seems fine

    #Remove an in-chain atom:
    #
    #NOT WORKING WELL
    elif Degree == 2:
        #print Chem.MolToSmiles(mol, True)
        nbr = []
        for neighb in atom.GetNeighbors():
            nbr.append(neighb)
        for bond in atom.GetBonds():
            mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        mol.RemoveAtom(atom.GetIdx())
        # this line give indexerrors:
        mol.AddBond(nbr[0].GetIdx(), nbr[1].GetIdx(), bondorder[1])
    elif True:
        pass

    #Remove a 3-bond atom:
    elif Degree == 3:
        nbr = []
        nbrCarbon = []
        nChoice = 3
        for neighb in atom.GetNeighbors():
            nbr.append(neighb)
            if neighb.GetAtomicNum()==6 and \
                   len(neighb.GetNeighbors())==3:
                nChoice = nChoice + 1
                nbrCarbon.append(neighb)
        choose = random.randrange(0, nChoice, 1)
        for bond in atom.GetBonds():
            mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetBeginAtomIdx())
        mol.RemoveAtom(atom.GetIdx())
        if choose == 0:
            mol.AddBond(nbr[0].GetIdx(), nbr[1].GetIdx(), bondorder[1])
            mol.AddBond(nbr[1].GetIdx(), nbr[2].GetIdx(), bondorder[1])
        elif choose == 1:
            mol.AddBond(nbr[0].GetIdx(), nbr[2].GetIdx(), bondorder[1])
            mol.AddBond(nbr[2].GetIdx(), nbr[1].GetIdx(), bondorder[1])
        elif choose == 2:
            mol.AddBond(nbr[0].GetIdx(), nbr[1].GetIdx(), bondorder[1])
            mol.AddBond(nbr[2].GetIdx(), nbr[0].GetIdx(), bondorder[1])
        else:
            for neighb in nbr:
                if not neighb.GetIdx() == nbrCarbon[choose - 3].GetIdx():
                    mol.AddBond(neighb.GetIdx(),
                                nbrCarbon[choose - 3].GetIdx(), bondorder[1])

    #Remove a 4-bond atom
    elif Degree == 4:
        nbr = []
        for neighb in atom.GetNeighbors():
            nbr.append(neighb)
        for bond in atom.GetBonds():
            mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        mol.RemoveAtom(atom.GetIdx())
        n1 = nbr.pop(random.randrange(0, 4, 1))

        if random.random() < 0.75:  #XA4 -> A-A-A-A
            n2 = nbr.pop(random.randrange(0, 3, 1))
            mol.AddBond(n1.GetIdx(), n2.GetIdx(), bondorder[1])
            n3 = nbr.pop(random.randrange(0, 2, 1))
            mol.AddBond(n2.GetIdx(), n3.GetIdx(), bondorder[1])
            mol.AddBond(n3.GetIdx(), nbr[0].GetIdx(), bondorder[1])
        else:  #XA4 -> A(A3)
            for neighb in nbr:
                mol.AddBond(n1.GetIdx(), neighb.GetIdx(), bondorder[1])

    #Make sure nothing is bonded twice
    oldpairs = []
    todelete = []
    for bond in mol.GetBonds():
        pair = [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]
        pair.sort()
        if pair in oldpairs: todelete.append(bond)
        else: oldpairs.append(pair)
    for bond in todelete:
        mol.RemoveBond(bond.GetIdx())

    # Make sure bonding makes sense
    for atom in mol.GetAtoms():
        try:
            if EmptyValence(atom) < 0:
                for bond in atom.GetBonds():
                    if not bond.HasProp('group'):
                        bond.SetBondType(bondorder[1])
                if EmptyValence(atom) < 0: raise MutateFail(mol)
        except KeyError:
            raise MutateFail(mol)

    if debug: # print a file with the results from these mutations:
        with open('delatoms.smi','a') as f:
            f.write(Chem.MolToSmiles(mol) + '\n')
    return mol
mutators['DelAtom'] = Mutator(DelAtom)

def AddAroRing(mol):
    freesinglebonds = GetFreeBonds(mol, order=1, sides=True)
    freedoublebonds = GetFreeBonds(mol, order=2, notprop='group')
    triplebonds = filter(lambda bond: bond.GetBondType() == bondorder[3],
                         mol.GetBonds())
    correctbonds = freedoublebonds + triplebonds + freesinglebonds
    if len(correctbonds) < 1: raise MutateFail()
    bond = random.choice(correctbonds)

    # if double bond
    if all( atom.GetImplicitValence() for atom in (bond.GetBeginAtom(), bond.GetEndAtom())):
        pass
    elif bond.GetBondType() == bondorder[3]:
        bond.SetBondType(bondorder[2])
    else:
        print "wrong kind of atom given"
        raise MutateFail()

    def AwithLabel(label, idx=True):
        atom = filter(lambda atom: atom.HasProp(label),
                      mol.GetAtoms())[0]
        if idx:
            return atom.GetIdx()
        else:
            return atom

    butadiene = Chem.MolFromSmiles('C=CC=C')
    butadiene.GetAtomWithIdx(0).SetBoolProp('buta1', True)
    butadiene.GetAtomWithIdx(3).SetBoolProp('buta2', True)
    bond.GetBeginAtom().SetBoolProp('ah1', True)
    bond.GetEndAtom().SetBoolProp('ah2', True)
    mol = Chem.RWMol(Chem.CombineMols(mol, butadiene))
    try:
        mol.AddBond(AwithLabel('buta1'), AwithLabel('ah1'), bondorder[1])
        mol.AddBond(AwithLabel('buta2'), AwithLabel('ah2'), bondorder[1])
    except RuntimeError:
        raise MutateFail()
    finally:
        for prop in ['buta1', 'buta2', 'ah1', 'ah2']:
            AwithLabel(prop, idx=False).ClearProp(prop)

    return mol
mutators['AddAroRing'] = Mutator(AddAroRing)

def AddFusionRing(mol):
    try:
        p = Chem.MolFromSmarts('[h]@&=*(@*)@[h]')
        matches = mol.GetSubstructMatches(p)
    except RuntimeError:
        raise MutateFail()

    if len(matches)==0:
        raise MutateFail()

    match = random.choice(matches)

    if mol.GetAtomWithIdx(match[-1]).GetNumRadicalElectrons():
        raise MutateFail()

    def AwithLabel(label, idx=True):
        atom = filter(lambda atom: atom.HasProp(label),
                      mol.GetAtoms())[0]
        if idx:
            return atom.GetIdx()
        else:
            return atom

    propene = Chem.MolFromSmiles('C=CC')
    propene.GetAtomWithIdx(0).SetBoolProp('propane1', True)
    propene.GetAtomWithIdx(2).SetBoolProp('propane2', True)
    mol.GetAtomWithIdx(match[3]).SetBoolProp('ah1', True)
    mol.GetAtomWithIdx(match[0]).SetBoolProp('ah2', True)
    mol = Chem.RWMol(Chem.CombineMols(mol, propene))
    try:
        mol.AddBond(AwithLabel('propane1'), AwithLabel('ah1'), bondorder[1])
        mol.AddBond(AwithLabel('propane2'), AwithLabel('ah2'), bondorder[1])
    except RuntimeError:
        raise MutateFail()
    finally:
        for prop in ['propane1', 'propane2', 'ah1', 'ah2']:
            AwithLabel(prop, idx=False).ClearProp(prop)

    return mol
mutators['AddFusionRing'] = Mutator(AddFusionRing, nname='nAddAroRing')

# this is the interface for a possible user-input mutator
def CustomMutator(mol):
    return mol


