#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os
import sys
import cPickle as pickle
import gzip
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from molfails import NoGeom
'''
this is the helper module that contains helping functions that relates to
both chemical and non-chemical perspectives
'''
########################################
#       PRE-DEFINED VARIABLES
########################################

atomNum = {
    'Ac': 89,
    'Ag': 47,
    'Al': 13,
    'Am': 95,
    'Ar': 18,
    'As': 33,
    'At': 85,
    'Au': 79,
    'B': 5,
    'Ba': 56,
    'Be': 4,
    'Bh': 107,
    'Bi': 83,
    'Bk': 97,
    'Br': 35,
    'C': 6,
    'Ca': 20,
    'Cd': 48,
    'Ce': 58,
    'Cf': 98,
    'Cl': 17,
    'Cm': 96,
    'Cn': 112,
    'Co': 27,
    'Cr': 24,
    'Cs': 55,
    'Cu': 29,
    'Db': 105,
    'Ds': 110,
    'Dy': 66,
    'Er': 68,
    'Es': 99,
    'Eu': 63,
    'F': 9,
    'Fe': 26,
    'Fm': 100,
    'Fr': 87,
    'Ga': 31,
    'Gd': 64,
    'Ge': 32,
    'H': 1,
    'He': 2,
    'Hf': 72,
    'Hg': 80,
    'Ho': 67,
    'Hs': 108,
    'I': 53,
    'In': 49,
    'Ir': 77,
    'K': 19,
    'Kr': 36,
    'La': 57,
    'Li': 3,
    'Lr': 103,
    'Lu': 71,
    'Md': 101,
    'Mg': 12,
    'Mn': 25,
    'Mo': 42,
    'Mt': 109,
    'N': 7,
    'Na': 11,
    'Nb': 41,
    'Nd': 60,
    'Ne': 10,
    'Ni': 28,
    'No': 102,
    'Np': 93,
    'O': 8,
    'Os': 76,
    'P': 15,
    'Pa': 91,
    'Pb': 82,
    'Pd': 46,
    'Pm': 61,
    'Po': 84,
    'Pr': 59,
    'Pt': 78,
    'Pu': 94,
    'Ra': 88,
    'Rb': 37,
    'Re': 75,
    'Rf': 104,
    'Rg': 111,
    'Rh': 45,
    'Rn': 86,
    'Ru': 44,
    'S': 16,
    'Sb': 51,
    'Sc': 21,
    'Se': 34,
    'Sg': 106,
    'Si': 14,
    'Sm': 62,
    'Sn': 50,
    'Sr': 38,
    'Ta': 73,
    'Tb': 65,
    'Tc': 43,
    'Te': 52,
    'Th': 90,
    'Ti': 22,
    'Tl': 81,
    'Tm': 69,
    'U': 92,
    'Uuh': 116,
    'Uuo': 118,
    'Uup': 115,
    'Uuq': 114,
    'Uus': 117,
    'Uut': 113,
    'V': 23,
    'W': 74,
    'Xe': 54,
    'Y': 39,
    'Yb': 70,
    'Zn': 30,
    'Zr': 40
}
elements = {}
for element, num in atomNum.iteritems():
    elements[num] = element

########################################
#       NON-CHEMICAL RELATED
########################################


# read a zip file directly into memory, then iterate through it.
# much faster than using gzip.open, which is broken
def FastGZ(fname):
    import cStringIO
    io_method = cStringIO.StringIO
    import subprocess
    p = subprocess.Popen(["zcat", fname], stdout=subprocess.PIPE)
    fh = io_method(p.communicate()[0])
    assert p.returncode == 0

    return fh


def Depickle(fname, gz=False):
    if fname.split('.')[-1] == 'gz':
        gz = True
    if gz:
        myfile = FastGZ(fname)
    else:
        myfile = open(fname, 'r')
    tmp = pickle.load(myfile)
    myfile.close()

    return tmp


def Enpickle(obj, fname, protocal=2, gz=False):
    if fname.split('.')[-1] == 'gz':
        gz = True
    if gz:
        myfile = gzip.open(fname, 'wb')
    else:
        myfile = open(fname, 'wb')
    pickle.dump(obj, myfile, protocal)
    myfile.close()


# get SMILES name of specific atom
def GetSmilesName(atom):
    global elements
    symbol = elements[atom.GetAtomicNum()]
    if atom.GetIsAromatic():
        symbol = symbol.lower()
    else:
        symbol = symbol.upper()

    return symbol


# extract number in a string
def ExtractNum(string):
    outstr = ''
    for char in string:
        if char.isdigit():
            outstr += charg

    return int(outstr)


# get the file format
def GetFileFormat(filename):
    fields = filename.split('.')
    if fields[-1].lower() == 'gz':
        return fields[-2].lower()
    else:
        return fields[-1].lower()


def GetBaseName(filename, stripstr=True):
    if stripstr == True:
        return filename.split('/')[-1].split('.')[0]
    elif stripstr:
        return filename.split('/')[-1].rstrip(stripstr)
    else:
        return name.split('/')[-1]


# set a unique scratch directory (specifically for ET-MEI cluster use)
_scratch_loc = '/scr/' + os.environ['USER'] + '/'
_scratch_dir = ''


def SetScratchDir(location=None):
    global _scratch_dir, _scratch_dir

    if location is not None:
        _scratch_dir = location
        return _scratch_dir

    #scratch already set
    if _scratch_dir != '':
        return _scratch_dir

    #parallel case
    try:
        from mpi4py import MPI
        myRank = MPI.COMM_WORLD.Get_rank()
    except ImportError:
        myRank = 0

    try:
        jobID = os.environ['JOB_ID']
    except KeyError:
        jobID = 'off_queue_job'

    _scratch_dir = _scratch_loc + str(jobID) + '_' + str(myRank) + '/'
    if not os.path.isdir(_scratch_dir):
        os.makedirs(_scratch_dir)

    print 'Scratch directory set to: ' + _scratch_dir
    return _scratch_dir


_warnScratch = False


def GetScratchDir():
    global _warnScratch

    if _scratch_dir == '' and not _warnScratch:
        print 'WARNING: scratch dir not set, current directory will be used.'
        _warnScratch = True

    return _scratch_dir


def Normalize(vec):
    n = np.sqrt(np.dot(vec, vec))
    return vec / n


def DumpMols(lib, gen=None, filename=None):
    if filename: # force filename to have .smi extension
        if not filename[-4:]=='.smi':
            filename = '{}.smi'.format(filename)
    else:
        if gen is None:
            filename = 'pool.smi'
        else:
            filename = 'it{}.smi'.format(gen)
    w = Chem.SmilesWriter(filename)
    for mol in lib:
        w.write(mol)
    return


def FinishSelection(lib):
    for i, mol in enumerate(lib):
        try:
            mol.SetIntProp('selected', mol.GetIntProp('selected') + 1)
        except KeyError:
            # if mol doesn't have 'selected' set to 1
            mol.SetIntProp('selected', 1)


########################################
#       CHEMICAL RELATED
########################################

from contextlib import contextmanager


@contextmanager
def custom_redirection(fileobj):
    old = sys.stdout
    sys.stdout = fileobj
    try:
        yield fileobj
    finally:
        sys.stdout = old


def xyzfromrdmol(rdmol, string=False, **kwargs):
    atomN = {
        1: 'H',
        5: 'B',
        6: 'C',
        7: 'N',
        8: 'O',
        9: 'F',
        14: 'Si',
        15: 'P',
        16: 'S',
        17: 'Cl',
        35: 'Br',
        53: 'I'
    }
    with open('omega.mol', 'a') as out:
        with custom_redirection(out):
            molcoords = Compute3DCoords(rdmol, **kwargs)
            if string:
                xyz = molcoords
            else:
                reorder = lambda xyz: "{3} {0} {1} {2}".format(*xyz)
                cartesian = map(reorder, molcoords)
                xyz = "\n".join(cartesian)
                xyz += '\n\n'
    return xyz


# calculate 3D coordinate of a molecule using RDKit build-in functions
# note that output of this function is often redirected to omega.out
def Compute3DCoords(mol, ff='ETKDG', RDGenConfs=False, **kwargs):
    ''' different types of 3D coordinate acquisition are available
    here the ETKDG / UFF / MMFF are provided.
    '''

    #mol: rdkit RWMol or Mol
    if RDGenConfs:
        try:
            from rdkithelpers import ConformerGenerator
            generator = ConformerGenerator(max_conformers=1, **kwargs)
            molAddH = generator(mol)
        except RuntimeError as e:
            raise NoGeom
    else:
        try:
            molAddH = Chem.AddHs(mol)
            if ff == 'ETKDG':
                #calculate mol 3D coordinate
                AllChem.EmbedMolecule(molAddH, AllChem.ETKDG())
            elif ff == 'UFF':
                AllChem.EmbedMolecule(molAddH, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
                AllChem.UFFOptimizeMolecule(molAddH)
            elif ff == 'MMFF':
                AllChem.EmbedMolecule(molAddH, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
                AllChem.MMFFOptimizeMolecule(molAddH)
            else:
                raise TypeError('force field not recognized')
        except ValueError:
            raise NoGeom

    molStr = Chem.MolToMolBlock(molAddH)

    #parse mol string file into list
    molCoords = []

    molStr2List = molStr.split('\n')

    #get number of atoms
    num = int(molStr2List[3].split()[0])

    for i in xrange(4, 4 + num):
        # x, y, z, atomname format
        coords = molStr2List[i].split()[0:4]
        molCoords.append(coords)

    return molCoords


# get the Murckoscaffold of a molecule
def GetMurckoScaffold(mol):
    #mol: rdkit RWMol or Mol
    from rdkit.Chem.Scaffolds import MurckoScaffold

    scaffold = MurckoScaffold.MakeScaffoldGeneric(mol)

    #return scaffold rdkit.mol object
    return scaffold

class FakeModule(object):
    def __init__(self):
        try:
            self.molfails = __import__('molfails')
        except ImportError as e:
            pass
        return

ACSESS = FakeModule()

if __name__=="__main__":
    import sys
    smiles = sys.argv[1]
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    print xyzfromrdmol(mol)

