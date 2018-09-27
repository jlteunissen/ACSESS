#!/usr/bin/env python
#-*- coding: utf-8 -*-

import random
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkithelpers import *
from molfails import MutateFail
#import mprms

MAXTRY = 100
elements = []
maxWeight = 0
MxAtm = 0

debug = False

###############################################
# Initialization the module
###############################################


def Init():
    #elements = mprms.elements
    global elements, halogens
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
    return

#########################################
# Returns free valence for an atom
# Obviously, there's a problem if it's negative
## After implicit hydrogens are assigned with the valence model,
## this can be replaced by the implicit hydrogen count
# MaxValence={6:4, 7:3, 8:2, 9:1, 15:3, 17:1, 16:2, 35:1, 50:3, 51:2}
#Jos added hydrogen in this list. I don't know if that is a good idea but for
# openshell a lot of H's will end up in the SMILES strings
MaxValence = {
    1: 1,
    6: 4,
    7: 3,
    8: 2,
    9: 1,
    15: 3,
    17: 1,
    16: 2,
    35: 1,
    50: 3,
    51: 2
}


def EmptyValence(atom):

    # Sulfur can have up to 2 double bonds to oxygen
    # (if it's not aromatic)
    if (atom.GetAtomicNum() == 16 and atom.HasProp('grouprep')) and (
            atom.GetProp('grouprep') == 'sulfone'
            or atom.GetProp('grouprep') == 'sf5'):
        maxv = 6

    elif (atom.GetAtomicNum and atom.HasProp('grouprep')
          and atom.GetProp('grouprep') == 'nitro'):
        maxv = 4

    else:
        try:
            maxv = MaxValence[atom.GetAtomicNum()]
        except KeyError:
            print "Error in EmptyValence"
            print atom.GetAtomicNum()
            raise

    return maxv - atom.GetExplicitValence() - atom.GetNumRadicalElectrons()


###############################################################################
#                            Mutation Methods                                 #
###############################################################################
#     Nota Bene! These mol objects are now based on Chem.RWMol objects!       #
###############################################################################


# 1. Flip a bond order
#@captureMolExceptions
def FlipBond(mol, bond):
    if bond.HasProp('group'):
        raise MutateFail(mol, 'In-group bond passed to FlipBond')
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


#@captureMolExceptions
def SwitchAtom(mol, atom):
    if atom.HasProp('group') and not atom.HasProp('grouprep'):
        raise MutateFail(mol, 'protected or grouped atom passed to switchatom')
    #elif atom.GetNumRadicalElectrons():
    #    raise MutateFail(mol, 'we do not allow radical center to be changed')

    changed = False
    neighbors=list(atom.GetNeighbors())
    valence = atom.GetExplicitValence() + atom.GetNumRadicalElectrons()
    #valence = atom.GetExplicitValence()

    # # # # # # # # # # # # #
    # from Filters import OptSulfone
    # Special rules for sulfurs that can optionally be sulfones
    #if atom.GetAtomicNum()==16 and len(OptSulfone)>0:
    #    SulfCand=set()
    #    for group in OptSulfone:
    #        for match in mol.GetSubstructMatches(group, True):
    #            for mat in GetIAtoms(match, mol):
    #                if mat.target.IsSulfur(): SulfCand.add(mat.target)
    #    if atom in SulfCand and random.random()>0.4:
    #        if atom.HasProp('grouprep'):
    #            RemoveGroup(mol,atom.GetProp('group'))
    #            return
    #        elif not atom.HasProp('group'):
    #            MakeSulfone(mol,atom)
    #            return

    #If it's a group representative, switching it will delete the
    #rest of the group; since they're hard to create, they're also
    #hard to destroy
    #if atom.HasProp('grouprep':
    #    if random.random()>0.3:
    #        RemoveGroup(mol,atom.GetProp('group'))
    #        changed=True
    #    else:
    #        raise MutateFail()

    # # # # # # # # # # # # # #
    #Special cases for Nitro groups and halogens,
    #which can only be on aromatic rings
    #if (  atom.GetExplicitValence()==1
    #      and neighbors[0].IsAromatic() ):

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
                return

    # If not a halogen, maybe make a nitro group
    #    if random.random() < 0.15:
    #        MakeNitro(mol,atom)
    #        return

    # # # # # # # # #
    # Regular atom switching
    atnum = atom.GetAtomicNum()
    elems = [el for el in elements if el != atnum]

    for itry in range(50):
        cand = random.choice(elems)
        if MaxValence[cand] >= valence:
            atom.SetAtomicNum(cand)
            return

    if not changed:
        raise MutateFail()  # if here, mutation failed ...
    return mol


#@captureMolExceptions
def AddBond(mol):
    # NB! This function was overloaded and atoms where possible arguments
    # Not anymore possible to get a cleaner interface

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


########################################################
# Remove a bond. Obviously, the bond must be in a cycle,
# or the molecule would be rent in twain
#@captureMolExceptions
def RemoveBond(mol, bond):
    if bond.HasProp('group'):
        raise MutateFail(mol,
                         'Protected bond (in group) passed to RemoveBond.')
    if not bond.IsInRing():
        raise MutateFail(mol, 'Non-ring bond passed to RemoveBond')
    mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    return mol


#################################################################
# Add an atom, either in a bond, or on a side-chain or terminus
#@captureMolExceptions
def AddAtom(mol, obj):
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
        newatomid = mol.AddAtom(Chem.Atom(6))
        mol.AddBond(obj.GetIdx(), newatomid, bondorder[1])
    else:
        raise MutateFail()
    return mol


#@captureMolExceptions
def RemoveAtom(mol, atom):
    #
    # Deal with special groups
    #
    if atom.HasProp('protected') or atom.HasProp('fixed'):
        raise MutateFatal(mol,
                          'Protected or fixed atom passed to' + " RemoveAtom.")
    Degree = atom.GetDegree() + atom.GetNumRadicalElectrons()
    if debug:
        print "D{:d}".format(Degree),
        print Chem.MolToSmiles(mol), atom.GetNumRadicalElectrons(), atom.GetAtomicNum(),

    # If atom is the representative of a larger group (e.g. N in nitro group)
    # delete the entire group
    if atom.HasProp('grouprep'):
        groupnum = atom.GetProp('group')
        RemoveGroup(mol, groupnum)

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

    try:
        mol = Finalize(mol, aromatic=False)
    except Exception as e:
        print "in Remove Atom with:", Chem.MolToSmiles(mol, True), e

    if debug: # print a file with the results from these mutations:
        with open('delatoms.smi','a') as f:
            f.write(Chem.MolToSmiles(mol) + '\n')
    return mol


def AddArRing(mol, bond):
    # if double bond
    if all( atom.GetImplicitValence() for atom in (bond.GetBeginAtom(), bond.GetEndAtom())):
        pass
    elif bond.GetBondType() == bondorder[3]:
        bond.SetBondType(bondorder[2])
    else:
        print "wrong kind of atom given"
        raise MutateFail

    #print "in AddArRing!", Chem.MolToSmiles(mol)
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
        raise MutateFail
    finally:
        for prop in ['buta1', 'buta2', 'ah1', 'ah2']:
            AwithLabel(prop, idx=False).ClearProp(prop)
    #print "Finished AddFusionRing!", Chem.MolToSmiles(mol)
    return mol


def AddFusionRing(mol, match):

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
        raise MutateFail
    finally:
        for prop in ['propane1', 'propane2', 'ah1', 'ah2']:
            AwithLabel(prop, idx=False).ClearProp(prop)
    #print "Finished AddFusionRing!", Chem.MolToSmiles(mol)
    return mol





