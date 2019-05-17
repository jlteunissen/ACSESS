#!/usr/bin/env python
#-*- coding: utf-8 -*-
''' This module acts as an interface between ACSESS and CINDES

    The result is that all Chem.Mol objects have an objective attribute

    by J.L. Teunissen
'''
import os, sys
import numpy as np
import molfails

from rdkit import Chem

from helpers import Compute3DCoords, xyzfromrdmol
#s3d.UseBalloon=True
#s3d.S3DInit()

from CINDES.utils.molecule import SmiMolecule
from CINDES.utils.writings import log_io
table = None
run = None
minimize = False  #Minimizing or maximizing the objective function?
RDGenConfs = False
pool_multiplier=10

def Init():
    global run, table
    from CINDES.inputreader import read_input
    from CINDES.utils.table import set_table

    # 1. read input
    run = read_input('INPUT')

    # 2. set optimum from mprms file so they are automatically similar
    if minimize:
        run.optimum = 'minimum'
    else:
        run.optimum = 'maximum'

    # 3. set table
    table = set_table(run)

    # 4. print input
    print run

    return run

def calculate(rdmols, QH2=False, gen=0):

    # -2 prepare optionally for logging to file
    global table, run
    print "len(table):{}".format(len(table)),
    print "n rdmols:{}".format(len(rdmols))
    #if len(table) < 10: print table

    # -1 imports
    from CINDES.evaluation.calculator import evaluate_mols

    # 1. RDKit.Chem.Mol(rdmol) objects have to have a molecular specification
    # CINDES.utils.molecule.SmiMolecule objects
    mols = [SmiMolecule(Chem.MolToSmiles(rdmol, True)) for rdmol in rdmols]

    # 2. set xyz coordinates of molecules:
    for mol, rdmol in zip(mols, rdmols):
        mol.rdmol=rdmol
        try:
            mol.xyz = xyzfromrdmol(rdmol, RDGenConfs=RDGenConfs, pool_multiplier=pool_multiplier)
        except (ValueError, molfails.NoGeom) as e:
            print "no geom constructed for now ignored molecule:", Chem.MolToSmiles(rdmol)
            print e
            mol.discard()
            continue

    # 3. Do the actual calculation:
    #mols_all, nnewcalcs, made_pred = evaluate_mols(run, mols, table, count=gen, nsite=0)
    evaluate_mols(run, mols, table, count=gen, nsite=0)

    # 4. Do the administration & set 'Objective' Property for RWMol objects
    table = loggings(mols, table)

    return

@log_io()
def loggings(mols, table):
    # Optional. log results to screen
    from CINDES.algorithms.loggings import log_screen, log_table

    log_screen(mols)

    # set the results as attributes from the rdmols
    for mol in mols:
        try:
            mol.rdmol.SetDoubleProp('Objective', float(mol.Pvalue))
        except AttributeError as e:
            print "mol has no rdmol?:", mol
            print "mol.ignored?", mol.IsDiscarded
            print e

    print "logging table..."
    table = log_table(mols, table)
    return table


def xyzfromstring(string):
    # split to get the first line which has the number of atoms
    splitted = string.split('\n', 2)
    nxyz = int(splitted[0])
    # take only as many atoms as are in one configuration:
    xyz = '\n'.join(splitted[2].split('\n', nxyz)[:-1])
    # add an empty line
    xyz += '\n\n'
    return xyz


if __name__ == "__main__":
    import sys

    class Unbuffered(object):
        def __init__(self, stream):
            self.stream = stream

        def write(self, data):
            self.stream.write(data)
            self.stream.flush()

        def __getattr__(self, attr):
            return getattr(self.stream, attr)

    sys.stdout = Unbuffered(sys.stdout)

    import Filters as fl
    fl.FLInit()

    if sys.argv[1][-4:] == '.smi':
        print "testset SMILES file detected\n\tperforming test calculation"
        run = getrun()
        with open(sys.argv[1]) as f:
            smiles = [line.strip('\n') for line in f.readlines()]
            print smiles
        rdmols = []
        for smi in smiles:
            mol = Chem.MolFromSmiles(smi, True)
            rdmols.append(mol)
        calculate(rdmols, run, log=True)
    print "done"
