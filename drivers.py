#!usr/bin/env python
#-*- coding: utf-8 -*-
'''
 This forms an middle layer between the main algorithm steps and the lowlevel functions
'''
from molfails import *
from rdkithelpers import *
import mutate, crossover
import filters
import mprms
import random, sys
random.seed(2)
from output import logtime, StartTimer, EndTimer, stats
import output

###### number of mutations
nMut = 0
nCross = 0


##### workflow switches:
startFilter = 2
startTautomer = None
startGenStruc = 0
KeepNoGeomPool = True

##### other variables
maxPool = 1000
EdgeLen = 10
EdgeRatio = 0.1
resonance = False

_tautomerEnumerator = None


debug = False

def Init():
    global _tautomerEnumerator
    if hasattr(mprms, 'canonicalTautomer') and mprms.canonicalTautomer:
        from molvs.tautomer import TautomerEnumerator
        _tautomerEnumerator = TautomerEnumerator()
        print "_tautomerEnumerator:", _tautomerEnumerator
    return


def SetIterationWorkflow(gen):
    Tautomerizing = ( (not startTautomer is None) and gen >= startTautomer)
    Filter = (gen >= startFilter)
    GenStruc = (gen >= startGenStruc)
    return Tautomerizing, Filter, GenStruc


@logtime()
def DriveMutations(lib):
    global stats
    nDups, nExcp, nCand = (0, 0, 0)

    # 1. CROSSOVERS:
    print "crossovers...",
    sys.stdout.flush()
    newmols = []
    for i in xrange(min(nCross, len(lib) - 1)):
        try:
            if i < nCross * EdgeRatio and EdgeLen > 0:
                mol1 = random.choice(lib[:EdgeLen])
                mol2 = random.choice(lib + newmols)
            else:
                mol1 = random.choice(lib + newmols)
                mol2 = random.choice(lib + newmols)
            #candidate = mutate.Crossover(mol1, mol2)
            candidate = crossover.Crossover(mol1, mol2)
            if type(candidate) == Chem.RWMol: candidate = candidate.GetMol()
            candidate = Finalize(candidate, aromatic=False)
            newmols.append(candidate)
            stats['nCand'] += 1
        except MutateFail:
            stats['nExcp'] += 1
    newmols = filter(Sane, newmols)

    nbefore = len(newmols)
    newmols = RemoveDuplicates(newmols)
    stats['nDups'] += nbefore - len(newmols)

    # 2. MUTATIONS
    print "mutating...",
    sys.stdout.flush()
    for i in xrange(nMut):
        # choose random mol:
        if i < nMut * EdgeRatio and EdgeLen > 0:
            candidate = Chem.Mol(random.choice(lib[:EdgeLen]))
        else:
            candidate = Chem.Mol(random.choice(lib + newmols))

        # optionally select a random tautomer
        if _tautomerEnumerator:
            candidate = SelectTautomer(candidate)

        # optionally select a random resonance_structure
        if resonance:
            candidate = SelectResonanceStructure(candidate)

        # try to find a correct mutation
        try:
            candidate = mutate.MakeMutations(candidate)
            if type(candidate) == Chem.RWMol:
                candidate = candidate.GetMol()
            try:
                candidate = Finalize(candidate, aromatic=False)
            except Exception as e:
                raise MutateFail(candidate)
            newmols.append(candidate)
            stats['nCand'] += 1
        except MutateFail:
            stats['nExcp'] += 1
    newmols = filter(Sane, newmols)
    nbefore = len(newmols)
    newmols = RemoveDuplicates(newmols)
    stats['nDups'] = nbefore - len(newmols)
    if debug:
        with open('mutatelib', 'w') as f:
            for mol in newmols:
                f.write(Chem.MolToSmiles(mol) + '\n')

    return newmols


###############################################
# Drive filtering
def DriveFilters(lib, Filtering, GenStrucs):
    global stats

    def torwmol(mol):
        if not type(mol) == Chem.RWMol: return Chem.RWMol(mol)
        else: return mol

    # every mol is converted to a RWMol otherwise the 
    # fixers will have no effect
    lib = map(torwmol, lib)
    nbefore = len(lib)

    if not (Filtering or GenStrucs):
        return lib

    ####### Parallel filtering ######
    if Filtering and mprms.mpi:
        import parallel as pl
        if GenStrucs:
            pl.ScatterFixFilterStruc(lib)
        else:
            pl.ScatterFixFilter(lib)

    #################### Serial filtering #####################
    elif Filtering:
        print "filtering...",
        sys.stdout.flush()
        StartTimer('Filters')
        for mol in lib:
            if not mol.HasProp('filtered'):
                filters.FixAndFilter(mol)
                assert mol.HasProp('failedfilter')
                if debug:
                    if not mol.HasProp('failedfilter'):
                        print "Jos Error",
                    else:
                        if mol.GetProp('failedfilter'):
                            print "FaFi:", mol.GetProp('failedfilter'),
            else:
                assert mol.GetBoolProp('filtered') == True
                mol.SetProp('failedfilter', '')
        nbefore = len(lib)
        lib = RemoveDuplicates(lib)
        stats['nDups'] += nbefore - len(lib)
        EndTimer('Filters')

        #### Generate 3d structures if requested ###
        if GenStrucs:
            StartTimer('GenStrucs')
            print "genstrucs..."
            sys.stdout.flush()
            for i, mol in enumerate(lib):
                if sys.stdout.isatty():
                    sys.stdout.write("{}%\r".format(int(100 * i / len(lib))))
                    sys.stdout.flush()
                if (mol.HasProp('hasstructure')
                        or mol.GetProp('failedfilter')):
                    continue

                from helpers import xyzfromrdmol
                try:
                    xyz = xyzfromrdmol(mol)
                    mol.SetProp('hasstructure', xyz)
                    if filters.GeomFilter(mol) == 'SAV':
                        mol.SetProp('failedfilter', 'SAV')
                except NoGeom:  #if structure can't be generated, filter it
                    mol.SetProp('failedfilter', 'failed to generate geometry')
                except MutateFatal:
                    from pprint import pprint as pp
                    print "Mutate Fatal Error with:", pp(mol.__dict__)
                    print "smiles:", Chem.MolToSmiles(mol), "end smi"
                    mol.SetProp('failedfilter', 'failed to generate geometry')
                    raise
            EndTimer('GenStrucs')

    ## Remove filtered compounds from library ##
    if Filtering:
        for mol in lib:
            try:
                failed = mol.GetProp('failedfilter')
            except KeyError:
                print "no failed filter:", Chem.MolToSmiles(mol)
            if failed:
                stats['nFilt'] += 1
                try:
                    smi = Chem.MolToSmiles(mol)
                except NotImplementedError:
                    smi = mol._smiles
                output.filterFile.write(smi + '  ' + failed + '\n')
        lib = [mol for mol in lib if not mol.GetProp('failedfilter')]
        output.filterFile.flush()
    return lib


def DrivePoolFilters(pool, Filtering, GenStruc, Tautomerizing, gen):
    pool = filter(bool, pool)
    if Tautomerizing and gen == startTautomer:
        print 'pool - tautomerizing ...',
        sys.stdout.flush()
        pool = map(Tautomerize, pool)
    if GenStruc and gen == startGenStruc and not KeepNoGeomPool:
        print 'restarting and filtering pool...',
        sys.stdout.flush()
        pool = [mol for mol in pool if mol.GetBoolProp('hasstructure')]
        DriveFilters(pool, Filtering, GenStruc)
    if (Filtering and gen == startFilter) or (GenStruc and gen == startGenStruc
                                              and KeepNoGeomPool):
        print 'pool -',
        pool = DriveFilters(pool, Filtering, GenStruc)
    return pool


@logtime()
def ExtendPool(pool, lib, newmol):
    def key(x):
        if not x.HasProp('selected'): x.SetIntProp('selected', 0)
        return x.GetIntProp('selected')

    from random import shuffle
    if len(pool) + len(newmol) <= maxPool:
        pool += newmol
        return RemoveDuplicates(pool)
    else:
        newpool = RemoveDuplicates(
            newmol + lib)  #this is the minumum pool content
        if len(newpool) < maxPool:
            shuffle(pool)
            oldpool = sorted(set(pool) - set(lib), key=key, reverse=True)
            newpool = newpool + oldpool[:maxPool - len(newpool)]
        return RemoveDuplicates(newpool)


@logtime()
def DriveSelection(pool, subsetSize):
    print "selecting...",
    sys.stdout.flush()
    #1. select maximin algorithm.
    if mprms._similarity:
        from similarity import FPMaximin
        lib = FPMaximin(pool, mprms.subsetSize)
    else:
        from distance import Maximin
        lib = Maximin(pool, mprms.subsetSize)
    return lib


############ EXTRA FUNCTIONS: #############


##############################################################
# Simple utility to remove duplicate structures
# Tests for identical molecules by generating smiles strings
# Take the molecule with more information generated
def GetScore(mol):
    score = 0
    for prop in ['filtered', 'hasstructure', 'Objective', 'tautomerized']:
        try:
            score += bool(mol.GetProp(prop))
        except KeyError:
            pass
    if not mol.HasProp('selected'): mol.SetIntProp('selected', 0)
    return score


def RemoveDuplicates(lib):
    def verify_isosmi(mol):
        if not mol.HasProp('isosmi'):
            mol.SetProp('isosmi', Chem.MolToSmiles(mol))
    map(verify_isosmi, lib)
    lib.sort(key=lambda x: x.GetProp('isosmi'))

    i = 1
    while i < len(lib):
        #print len(lib),
        if lib[i].GetProp('isosmi') == lib[i - 1].GetProp('isosmi'):
            iscore = GetScore(lib[i])
            i1score = GetScore(lib[i - 1])

            if iscore >= i1score:
                lib[i].SetIntProp('selected', lib[i].GetIntProp('selected') +
                                  lib[i - 1].GetIntProp('selected'))
                lib.pop(i - 1)
            else:
                lib[i - 1].SetIntProp('selected',
                                      lib[i].GetIntProp('selected') +
                                      lib[i - 1].GetIntProp('selected'))
                lib.pop(i)

        else:
            i += 1
    return lib

def SelectTautomer(candidate):
    smi1 = Chem.MolToSmiles(candidate)
    tautomers = _tautomerEnumerator(candidate)
    candidate = random.choice(tautomers)
    Chem.Kekulize(candidate, True)
    smi2 = Chem.MolToSmiles(candidate)
    if not smi1==smi2 and len(tautomers)>1:
        #print "from canonical {} using random tautomer {} of total {}".format(smi1, smi2, len(tautomers))
        pass
    return candidate
