#!/usr/bin/env python
#-*- coding: utf-8 -*-

from filters import NewFilter, NewPatternFilter
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkithelpers import *
from molfails import *
import random
try:
    from GraphPlanar import IsPlanar
except:
    print 'Could not import GraphPlanarity.so'
    print 'Compile by typing "make" in base of source directory.'
    IsPlanar = lambda x: True

DefaultFilters = {}
MxAtm = 0
maxWeight = 0
maxRings = 8
maxRingSize = 11
RingSizeExceptions = 1
UseRuleOf10 = False
SAScore = 0.0

############ FILTERS

DefaultFilters['Too big'] = NewFilter('Too big')

def TooBig(mol):
    global maxWeight
    if MxAtm > 0:
        if mol.GetNumHeavyAtoms() > MxAtm:
            return 'Too Big'
    if maxWeight > 0:
        if Descriptors.MolWt(mol) > maxWeight:
            return True
    return False

DefaultFilters['Too big'].SetFilterRoutine(TooBig)


def CutMoreRings(mol):
    try:
        Sanitize(mol)
        RI = mol.GetRingInfo()
    except Exception as e:
        return True
    if not IsPlanar(mol): return True
    nrings, sa, sb = SSSR(mol, force=True)
    if UseRuleOf10 and RuleOf10(mol) > 10: return True
    if RI.NumRings() > maxRings: return True
    return False


def CutRings(mol):
    import mutate
    MAXTRY = 10
    #non-planar graphs and highly cyclic: keep breaking rings until graph
    #is planar or has an allowed number of rings
    changed = False
    itry = 0
    while CutMoreRings(mol):
        itry += 1
        bondringids = mol.GetRingInfo().BondRings()
        # need to flatten bondringids:
        # fancy manner to flatten a list of tuples with possible duplicates
        bondringids = set.intersection(*map(set, bondringids))
        bonds = GetIBonds(bondringids, mol, notprop='group')
        if len(bonds) == 0: raise MutateFail()
        mol = mutate.DelBond(mol)
        if mol.HasProp('ringcount'): mol.ClearProp('ringcount')
        changed = True
        if itry >= MAXTRY: raise MutateFail()
    return changed


DefaultFilters['Non-planar graph (Euler critereon)'] = NewFilter(
    'Non-planar graph (Euler critereon)')


def NonPlanarEuler(mol):
    # NOTE. Somehow TEMPO crashes on this
    Radical = False
    nNodes = mol.GetNumHeavyAtoms()
    if nNodes >= 3 and mol.GetNumBonds() > 3 * nNodes - 6 and not Radical:
        return 'Non-planar graph (Euler critereon)'


DefaultFilters['Non-planar graph (Euler critereon)'].SetFilterRoutine(
    NonPlanarEuler)
DefaultFilters['Non-planar graph (Euler critereon)'].SetFixRoutine(CutRings)

DefaultFilters['Non-planar graph (Boyes)'] = NewFilter(
    'Non-planar graph (Boyes)')


def NotPlanarBoyes(mol):
    #print "In NotPlanarBoyes",
    return not IsPlanar(mol)


DefaultFilters['Non-planar graph (Boyes)'].SetFilterRoutine(NotPlanarBoyes)
DefaultFilters['Non-planar graph (Boyes)'].SetFixRoutine(CutRings)
DefaultFilters['Too many rings'] = NewFilter('Too many rings')

def TooManyRings(mol):
    # Ring features
    if mol.GetNumAtoms() > maxRings:
        nrings = mol.GetRingInfo().NumRings()
        if nrings > maxRings:
            return 'Too Many Rings: {}'.format(nrings)
    return False

DefaultFilters['Too many rings'].SetFilterRoutine(TooManyRings)
DefaultFilters['Too many rings'].SetFixRoutine(CutRings)

fname = 'SSSR ring bigger than max allowed size'
DefaultFilters[fname] = NewFilter(fname)


def BiggestRing(mol):
    if mol.GetNumAtoms() > maxRingSize:
        nrings, sa, sb = SSSR(mol)
        if sum(nrings[maxRingSize - 2:]) > RingSizeExceptions:
            failname = 'SSSR ring bigger than max allowed size'
            return failname


DefaultFilters[fname].SetFilterRoutine(BiggestRing)


def CutBiggestRings(mol):
    ringatoms, ringbonds = SSSR_GetRings(mol)
    bondids = {bond.GetIdx(): bond for bond in mol.GetBonds()}
    changed = False
    while True:
        nrings = map(int, mol.GetProp('ringcounts').split())
        if not sum(nrings[maxRingSize - 2:]) > RingSizeExceptions:
            break
        toobig = set()
        smallenough = set()
        for ring in ringbonds:
            if len(ring) > maxRingSize: toobig.update(ring)
            else: smallenough.update(ring)
        canremove = toobig - smallenough
        removelist = [
            bondids[id] for id in canremove
            if not bondids[id].HasProp('mygroup')
        ]
        if len(canremove) == 0: raise MutateFail()
        else:
            delbond = random.choice(removelist)
            i, j = (delbond.GetBeginAtomIdx(), delbond.GetEndAtomIdx())
            mol.RemoveBond(i, j)
        changed = True
        ringatoms, ringbonds = SSSR_GetRings(mol, True)
    return changed


DefaultFilters[fname].SetFixRoutine(CutBiggestRings)

LookUpFilter = NewFilter('Compound not in lookup table')


def lu(mol):
    smi = Chem.MolToSmiles(mol, True)
    return (smi not in lookUpTable)
LookUpFilter.SetFilterRoutine(lu)

LipinskiFilter = NewFilter('Lipinski violation')
def lp_routine(mol):
    return LipinskiRuleOf5(mol) > 1

def LipinskiRuleOf5(mol):
    PropCalc(mol)
    ofs.flush()
    return mol.GetProp('Lipinski violations')
LipinskiFilter.SetFilterRoutine(lp_routine)

RuleOf10Filter = NewFilter('Rule of 10')
def r10f_routine(mol):
    return RuleOf10(mol) > 10

def RuleOf10(mol):
    OEPerceiveChiral(mol)
    ringcount, na, nb = SSSR(mol)
    nStereos = OECount(mol, OEIsChiralAtom()) + OECount(mol, OEIsChiralBond())
    return sum(ringcount) + nStereos
RuleOf10Filter.SetFilterRoutine(r10f_routine)

SAScoreFilter = NewFilter("SA-Score synthetic accessibility")
def sascore_filt(mol):
    import os
    import SAS as sa
    score = sa.CalcSAScore(mol)
    if score > SAScore:
        return 'SAScore: ' + str(score)
    else:
        return False
SAScoreFilter.SetFilterRoutine(sascore_filt)
