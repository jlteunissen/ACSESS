#!/usr/bin/env python
#-*- coding: utf-8 -*-
debug = False

from rdkit import Chem
from rdkit.Chem import AllChem
import traceback
from molfails import MutateFail
import numpy as np
from copy import deepcopy

########### Some global variables
bondorder = {
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    3: Chem.BondType.TRIPLE
}
aromaticity = False
verbose  = True
canonicalTautomer = False

############ FROM CANONICAL #############################
def Finalize(mol, tautomerize=None, aromatic=aromaticity):
    ''' This function makes sure that a new mutated/crossover/fixfilter molecule:
        - corrects the implicit valence
        - removes data from parent molecules.
        - optionally converts the molecule to its canonical tautomer
        - checks if it is a valid molecule, i.e. Sanitize step. (without aromaticity)
    '''
    # 1.
    try:
        mol.UpdatePropertyCache()
    except ValueError:
        raise MutateFail(mol)

    # 2.
    RW = type(mol) == Chem.RWMol

    # 3.
    ResetProps(mol)

    # 4.
    if canonicalTautomer and not tautomerize is False:
        mol = Tautomerize(mol, aromatic)

    # 5.
    try:
        Sanitize(mol, aromatic)
    except Exception as e:
        print "Error in Finalize with", Chem.MolToSmiles(mol, False),
        if verbose: print e,
        if debug:
            for item in traceback.extract_stack():
                print item

    # 6.
    if aromatic: Chem.SetAromaticity(mol)

    # 7.
    if RW and not type(mol) == Chem.RWMol:
        return Chem.RWMol(mol)
    else:
        return mol


def ResetProps(mol):
    ''' makes a fresh molecule without properties inherited from there parents
        by deleting all listed properties and resetting the SMILES string. '''
    isosmi = Chem.MolToSmiles(mol, True)
    for prop in [
            'filtered', 'hasstructure', 'tautomerized', 'minimized',
            'selected', 'failed', 'failedfilter', 'Objective'
    #        'selected', 'failed', 'Objective'
    ]:
        mol.ClearProp(prop)
    mol.SetProp('isosmi', isosmi)
    return


def Sane(mol, *args, **kwargs):
    ''' A SanitizeCheck function; function returns False if molecule is not sane '''
    try:
        Sanitize(mol, *args, **kwargs)
        return True
    except:
        return False


#def ToRWMol(mol):
#    if type(mol)==Chem.RWMol: return mol

############ NEW RDKIT HELPERS: #########################
flatten = lambda X: tuple(set(i for x in X for i in x))


def MolAtIsInGroup(mol, groupid):
    requirement = lambda atom: atom.HasProp('group') and atom.GetProp('group') == groupid
    return filter(requirement, mol.GetAtoms())


def MolBondIsInGroup(mol, groupid):
    requirement = lambda bond: bond.HasProp('group') and bond.GetProp('group') == groupid
    return filter(requirement, mol.GetBonds())


CanRemoveAtom = lambda atom: not atom.HasProp('protected') and not atom.HasProp('fixed')
CanChangeAtom = lambda atom: (not atom.HasProp('group')) or atom.HasProp('grouprep') or atom.GetNumRadicalElectrons()


def GetSmallestRingSize(atom):
    return min([i for i in range(1, 20) if atom.IsInRingSize(i)])


# Three Indices based getters:
def GetIAtoms(indices, mol, notprop=None):
    if notprop:
        requirement = lambda atom: atom.GetIdx() in indices and not atom.HasProp(notprop)
    else:
        requirement = lambda atom: atom.GetIdx() in indices
    return filter(requirement, mol.GetAtoms())


def GetIBonds(indices, mol, notprop=None):
    if notprop:
        requirement = lambda bond: bond.GetIdx() in indices and not bond.HasProp(notprop)
    else:
        requirement = lambda bond: bond.GetIdx() in indices
    return filter(requirement, mol.GetBonds())


def GetAtomIBonds(atomids, mol):
    ''' gets bonds based on the atom indices '''
    bonds = []
    for bond in mol.GetBonds():
        if all(
                index in atomids for index in (bond.GetBeginAtomIdx(),
                                               bond.GetEndAtomIdx())):
            bonds.append(bond)
    return bonds


def GetAtoms(mol, notprop=None):
    if notprop:
        requirement = lambda atom: not atom.HasProp(notprop)
    else:
        requirement = lambda atom: True
    return filter(requirement, mol.GetAtoms())


def GetBonds(mol, notprop=None):
    if notprop:
        requirement = lambda bond: not bond.HasProp(notprop)
    else:
        requirement = lambda bond: True
    return filter(requirement, mol.GetBonds())


def GetXAtoms(mol, num=6, notprop=None):
    '''gets the number of atoms in mol with atomic num. default carbon, 6 and not has prop notprop
    for every atom it is checked if it has property indicated with the string notprop
    '''
    if notprop:
        requirement = lambda atom: atom.GetAtomicNum() == num and not atom.HasProp(notprop)
    else:
        requirement = lambda atom: atom.GetAtomicNum() == num
    return filter(requirement, mol.GetAtoms())


def GetFreeBonds(mol, order=None, notprop=None, sides=False):
    ''' This function is used by aromatic ring addition. It returns the bonds with specified order
    and on both bond.atoms at least one hydrogen. '''
    # 1. select bonds with required bondorder
    if order:
        IsOrder = lambda bond: bond.GetBondType() == bondorder[order]
    else:
        IsOrder = lambda x: True
    ordbonds = filter(IsOrder, mol.GetBonds())
    # 2. select bonds which have at least on H at each atom
    HasHs= lambda bond:all(
            atom.GetImplicitValence() - atom.GetNumRadicalElectrons()
            for atom in (bond.GetBeginAtom(), bond.GetEndAtom())
            )
    withHbonds = filter(HasHs, ordbonds)
    # 3. test if not has prop notprop
    if notprop:
        bonds = filter(lambda bond: not bond.HasProp(notprop), withHbonds)
    else:
        bonds = withHbonds
    # 4. test if single bond has two double bonds connected:
    if order==1 and sides:
        hasdoubleside    = lambda atom:2.0 in map(lambda x:x.GetBondTypeAsDouble(), atom.GetBonds())
        hastwodoublesides= lambda bond: all( map( hasdoubleside, ( bond.GetBeginAtom(), bond.GetEndAtom() ) ) )
        bonds = filter( hastwodoublesides, bonds)
    return bonds

#########################################
# Returns free valence for an atom
# Obviously, there's a problem if it's negative
## After implicit hydrogens are assigned with the valence model,
## this can be replaced by the implicit hydrogen count
# MaxValence={6:4, 7:3, 8:2, 9:1, 15:3, 17:1, 16:2, 35:1, 50:3, 51:2}
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

########## Set List properties


def SetListProp(mol, name, iterable):
    ''' Since rdkit molecules can only store single values, list properties
    are stored as strings. Values are separated by a space
    NB: float precision is hardcoded to 20 decimals'''
    try:
        string = ' '.join(['{:.20f}'.format(x) for x in iterable])
    except ValueError:
        print iterable
        print x, type(x)
        for x in iterable:
            try:
                float(x)
            except ValueError:
                print x, type(x)
                raise
        raise
    mol.SetProp(name, string)


def GetListProp(mol, name):
    string = mol.GetProp(name)
    return map(float, string.split())


########## function potentially useful to avoid aromaticity


def Sanitize(mol, aromatic=False):
    '''The rdkit sanitize step with the option to switch off aromaticity'''
    if aromatic or aromaticity:
        Chem.SanitizeMol(mol)
    else:
        Chem.SanitizeMol(
            mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_SETAROMATICITY)
        #sanitizeOps=Chem.SANITIZE_ALL^Chem.SANITIZE_KEKULIZE^\
        #Chem.SANITIZE_SETAROMATICITY^Chem.SANITIZE_CLEANUP^\
        #Chem.SANITIZE_CLEANUPCHIRALITY)
    return




################ Resonance Structures for radical and cations:

def Resonate(mol, matches):
    # make single displacements of [+1]-*=[+0] to [+0]=*-[+1]
    newmols = []
    for match in matches:
        newmol = deepcopy(mol)
        Chem.Kekulize(newmol, True)
        a0 = newmol.GetAtomWithIdx(match[0])
        a0.SetFormalCharge(0)
        for i, ai in enumerate(match[1:],1):
            bond = newmol.GetBondBetweenAtoms(match[i-1], match[i])
            bondtype = bond.GetBondType()
            if bondtype==Chem.BondType.SINGLE:
                bond.SetBondType(Chem.BondType.DOUBLE)
            else:
                assert bondtype==Chem.BondType.DOUBLE
                bond.SetBondType(Chem.BondType.SINGLE)
        alast = newmol.GetAtomWithIdx(match[-1])
        #nH = alast.GetNumImplicitHs() + alast.GetNumExplicitHs()
        alast.SetFormalCharge(1)
        #alast.SetNumExplititHs(nH)
        newmols.append(newmol)
    return newmols

def Resonate_p2(mol, matches):
    newmols = []
    for match in matches:
        newmol = deepcopy(mol)
        Chem.Kekulize(newmol, True)

        # switch positive charge
        a0 = newmol.GetAtomWithIdx(match[0])
        a0.SetFormalCharge(0)
        a1 = newmol.GetAtomWithIdx(match[1])
        a1.SetFormalCharge(1)

        # make bond a double bond
        bond = newmol.GetBondBetweenAtoms(match[0], match[1])
        assert bond.GetBondType()==Chem.BondType.SINGLE
        bond.SetBondType(Chem.BondType.DOUBLE)

        newmols.append(newmol)
    return newmols

def FullyResonate(mol, idx):
    # find a conjugated cation.
    cationpattern1 = Chem.MolFromSmarts('[+1]~*=*')
    cationpattern2 = Chem.MolFromSmarts('[#7&+1]-[#16&X2]')

    radpositions = set()
    radpositions.add(idx)

    Chem.Kekulize(mol, True)

    # make a selection of resonant structures including our initial molecule:
    #mols = [mol]
    molsD= {idx:mol}
    while True:
        l1 = len(radpositions)
        for mol in molsD.values():

            # resonate pattern 1
            matches = mol.GetSubstructMatches(cationpattern1)
            # restrict matches if atoms were already conjugated:
            for match in matches:
                if not match[-1] in radpositions:
                    newstruct = Resonate(mol, [match])[0]
                    radpositions.add(match[-1])
                    molsD[match[-1]]=newstruct

            # resonate pattern 2
            matchesp2 = mol.GetSubstructMatches(cationpattern2)
            for matchp2 in matchesp2:
                if not matchp2[-1] in radpositions:
                    newstruct = Resonate_p2(mol, [matchp2])[0]
                    radpositions.add(matchp2[-1])
                    molsD[matchp2[-1]]=newstruct

        l2 = len(radpositions)
        # if no new conjugated atoms where found:
        if l1==l2:
            break
    #print "molsD:"
    #for i, imol in molsD.iteritems():
    #    print i, Chem.MolToSmiles(imol)
    return molsD

def SelectResonanceStructure(mol):
    smi1 = Chem.MolToSmiles(mol)
    doublet = False
    idx = None

    # convert radical to cation
    for a in mol.GetAtoms():
        if a.GetNumRadicalElectrons()==1:
            doublet = True
            idx = a.GetIdx()
            a.SetNumRadicalElectrons(0)
            a.SetFormalCharge(1)
    if idx is None:
        print "no idx!", Chem.MolToSmiles(mol)

    # here do resonance enumeration
    if False: # this should work in RDKit but it doesn't
        resSuppl = list(Chem.ResonanceMolSupplier(mol, flags=Chem.KEKULE_ALL))
        print "n resonance:", len(resSuppl), "for", Chem.MolToSmiles(mol)
        newmol = np.random.choice(resSuppl)
    else:
        newmols = FullyResonate(mol, idx).values()
        #print "n resonance:", len(newmols), "for", Chem.MolToSmiles(mol)
        newmol = np.random.choice(newmols)

    # reconvert to radical
    if doublet:
        for a in newmol.GetAtoms():
            if not a.GetFormalCharge()==0:
                a.SetNumRadicalElectrons(1)
                a.SetFormalCharge(0)
    Chem.Kekulize(newmol, True)

    smi2 = Chem.MolToSmiles(newmol)
    if not smi1==smi2:
        pass
        #print "from {} took resonant {}".format(smi1, smi2)

    try:
        Finalize(newmol)
    except (ValueError, MutateFail) as e:
        return mol
    return newmol


########## TAUTOMERIZING:


def Tautomerize(mol, aromatic=aromaticity):
    try:
        if mol.GetBoolProp('tautomerized'): return mol
    except KeyError:
        pass


    Chem.SanitizeMol(mol)
    if not (aromatic or aromaticity):
        Chem.Kekulize(mol, True)

    smi1 = Chem.MolToSmiles(mol)
    from molvs import Standardizer
    s = Standardizer()
    try:
        molnew = s.standardize(mol)
    except ValueError as e:
        raise MutateFail(mol)

    if not aromatic:
        Chem.Kekulize(molnew, True)
    smi2 = Chem.MolToSmiles(molnew)

    if smi1 == smi2:
        # we return mol because it contains some properties
        # tautomerized mols need to get the props again
        mol.SetBoolProp('tautomerized', True)
        return mol
    else:
        if mol.HasProp('failedfilter'):
            ff = mol.GetProp('failedfilter')
            molnew.SetProp('failedfilter', ff)
        #print "tautomerized:", smi1, 'to:', smi2
        with open('tautomerized.smi', 'a') as f:
            f.write("{} {}\n".format(smi1, smi2))
        molnew.SetBoolProp('tautomerized', True)
        return molnew


####### Complex Ring Functions
ringsearch = {}


def SSSR(mol, force=False):

    if (not force) and mol.HasProp('ringcounts'):
        ringcounts = map(int, mol.GetProp('ringcounts').split())
        sharedatoms = mol.GetIntProp('sharedatoms')
        sharedbonds = mol.GetIntProp('sharedbonds')
        return (ringcounts, ) + (sharedatoms, sharedbonds)
    #(map(int, mol.GetProp('ringcounts').split()),)+mol.GetIntProp('sharedcounts'

    ## SSSR ring analysis
    RI = mol.GetRingInfo()

    AssignedAtoms = set()
    AssignedBonds = set()
    SharedAtoms = set()
    SharedBonds = set()

    #nringatom= oe.OECount(mol, AtInRing )
    #nringbond= oe.OECount(mol, BondInRing )
    nringatom = len(flatten(RI.AtomRings()))
    nringbond = len(flatten(RI.BondRings()))

    nRings = [0] * max(nringatom, 8)

    #Loop over all possible ring sizes
    for i in xrange(3, nringatom + 1):
        if (len(AssignedAtoms) == nringatom
                and len(AssignedBonds) == nringbond):
            break

        if not ringsearch.has_key(i):
            ringsearch[i] = Chem.MolFromSmarts('*~1' + '~*' * (i - 1) + '1')

        #Find all instances of ring size i
        #matches=ringsearch[i].Match(mol,True)
        matches = mol.GetSubstructMatches(ringsearch[i])
        for match in matches:
            atomids = set(match)
            bondids = set(
                [bond.GetIdx() for bond in GetAtomIBonds(match, mol)])
            #bondids=set( [bond.target.GetIdx()
            #              for bond in match.GetBonds() ] )

            #Count this ring only if some of its atoms or bonds
            #have not already been assigned to a smaller ring
            if not (atomids.issubset(AssignedAtoms)
                    and bondids.issubset(AssignedBonds)):
                nRings[i - 3] += 1
                SharedAtoms.update(atomids.intersection(AssignedAtoms))
                SharedBonds.update(bondids.intersection(AssignedBonds))
                AssignedAtoms.update(atomids)
                AssignedBonds.update(bondids)

    mol.SetProp('ringcounts', " ".join(map(str, nRings)))
    mol.SetIntProp('sharedatoms', len(SharedAtoms))
    mol.SetIntProp('sharedbonds', len(SharedBonds))
    return nRings, len(SharedAtoms), len(SharedBonds)


#Smallest set of smallest rings analysis
#Just returns counts of various ring sizes. While SSSR in general
#is not invariant, these counts are.
ringsearch = {}


def SSSR_GetRings(mol, force=False):

    if not force and mol.HasProp('numSSSRrings'):
        ringatoms = []
        ringbonds = []
        numring = mol.GetProp('numSSSRrings')
        for i in xrange(numring):
            ringatoms.append(
                set(map(int,
                        mol.GetProp('AtomSSSR_' + str(i)).split())))
            ringbonds.append(
                set(map(int,
                        mol.GetProp('BondSSSR_' + str(i)).split())))
        return ringatoms, ringbonds

    ## SSSR ring analysis
    RI = mol.GetRingInfo()
    #oe.OEFindRingAtomsAndBonds(mol)
    AssignedAtoms = set()
    AssignedBonds = set()
    SharedAtoms = set()
    SharedBonds = set()

    nringatom = len(flatten(RI.AtomRings()))
    nringbond = len(flatten(RI.BondRings()))
    nRings = [0] * max(nringatom, 8)

    #Loop over all possible ring sizes
    ringatoms = []
    ringbonds = []
    for i in xrange(3, nringatom + 1):
        if (len(AssignedAtoms) == nringatom
                and len(AssignedBonds) == nringbond):
            break
        if not ringsearch.has_key(i):
            ringsearch[i] = Chem.MolFromSmarts('*~1' + '~*' * (i - 1) + '1')

        #Find all instances of ring size i
        matches = mol.GetSubstructMatches(ringsearch[i])
        for match in matches:
            if len(AssignedAtoms)==nringatom and \
               len(AssignedBonds)==nringbond:
                break
            atomids = set(match)
            try:
                bondids = set(
                    [bond.GetIdx() for bond in GetAtomIBonds(match, mol)])
            except AttributeError:
                print "mol:", mol
                print "matches:", matches
                raise

            #Count this ring only if some of its atoms or bonds
            #have not already been assigned to a smaller ring
            if not (atomids.issubset(AssignedAtoms)
                    and bondids.issubset(AssignedBonds)):
                ringatoms.append(atomids)
                ringbonds.append(bondids)
                nRings[i - 3] += 1
                SharedAtoms.update(atomids.intersection(AssignedAtoms))
                SharedBonds.update(bondids.intersection(AssignedBonds))
                AssignedAtoms.update(atomids)
                AssignedBonds.update(bondids)

    mol.SetProp('ringcounts', " ".join(map(str, nRings)))
    mol.SetIntProp('sharedatoms', len(SharedAtoms))
    mol.SetIntProp('sharedbonds', len(SharedBonds))
    mol.SetIntProp('numSSSRrings', len(ringatoms))
    for i in xrange(len(ringatoms)):
        mol.SetProp('AtomSSSR_' + str(i), " ".join(
            map(str, list(ringatoms[i]))))
        mol.SetProp('BondSSSR_' + str(i), " ".join(
            map(str, list(ringbonds[i]))))
    return ringatoms, ringbonds


"""
Conformer generation.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "3-clause BSD"


class ConformerGenerator(object):
    """
    Generate molecule conformers.
    Procedure
    ---------
    1. Generate a pool of conformers.
    2. Minimize conformers.
    3. Prune conformers using an RMSD threshold.
    Note that pruning is done _after_ minimization, which differs from the
    protocol described in the references.
    References
    ----------
    * http://rdkit.org/docs/GettingStartedInPython.html
      #working-with-3d-molecules
    * http://pubs.acs.org/doi/full/10.1021/ci2004658
    Parameters
    ----------
    max_conformers : int, optional (default 1)
        Maximum number of conformers to generate (after pruning).
    rmsd_threshold : float, optional (default 0.5)
        RMSD threshold for pruning conformers. If None or negative, no
        pruning is performed.
    force_field : str, optional (default 'uff')
        Force field to use for conformer energy calculation and
        minimization. Options are 'uff', 'mmff94', and 'mmff94s'.
    pool_multiplier : int, optional (default 10)
        Factor to multiply by max_conformers to generate the initial
        conformer pool. Since conformers are pruned after energy
        minimization, increasing the size of the pool increases the chance
        of identifying max_conformers unique conformers.
    """
    def __init__(self, max_conformers=1, rmsd_threshold=0.5, force_field='mmff94s',
                 pool_multiplier=10):
        self.max_conformers = max_conformers
        if rmsd_threshold is None or rmsd_threshold < 0:
            rmsd_threshold = -1.
        self.rmsd_threshold = rmsd_threshold
        self.force_field = force_field
        print "forcefield:", force_field
        self.pool_multiplier = pool_multiplier

    def __call__(self, mol):
        """
        Generate conformers for a molecule.
        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        return self.generate_conformers(mol)

    def generate_conformers(self, mol):
        """
        Generate conformers for a molecule.
        This function returns a copy of the original molecule with embedded
        conformers.
        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """

        # initial embedding
        print "embedding...",
        mol = self.embed_molecule(mol)
        if not mol.GetNumConformers():
            msg = 'No conformers generated for molecule'
            if mol.HasProp('_Name'):
                name = mol.GetProp('_Name')
                msg += ' "{}".'.format(name)
            else:
                smi = Chem.MolToSmiles(mol)
                msg += '. ({smi})'.format(smi=smi)
            raise RuntimeError(msg)

        # minimization and pruning
        print "minimizing...",
        self.minimize_conformers(mol)
        print "pruning...",
        mol = self.prune_conformers(mol)

        return mol

    def embed_molecule(self, mol):
        """
        Generate conformers, possibly with pruning.
        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        mol = Chem.AddHs(mol)  # add hydrogens
        n_confs = self.max_conformers * self.pool_multiplier
        AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, pruneRmsThresh=-1.)
        return mol

    def get_molecule_force_field(self, mol, conf_id=None, **kwargs):
        """
        Get a force field for a molecule.
        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        conf_id : int, optional
            ID of the conformer to associate with the force field.
        kwargs : dict, optional
            Keyword arguments for force field constructor.
        """
        if self.force_field == 'uff':
            ff = AllChem.UFFGetMoleculeForceField(
                mol, confId=conf_id, **kwargs)
        elif self.force_field.startswith('mmff'):
            AllChem.MMFFSanitizeMolecule(mol)
            mmff_props = AllChem.MMFFGetMoleculeProperties(
                mol, mmffVariant=self.force_field)
            ff = AllChem.MMFFGetMoleculeForceField(
                mol, mmff_props, confId=conf_id, **kwargs)
        else:
            raise ValueError("Invalid force_field " +
                             "'{}'.".format(self.force_field))
        return ff

    def minimize_conformers(self, mol):
        """
        Minimize molecule conformers.
        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        for conf in mol.GetConformers():
            ff = self.get_molecule_force_field(mol, conf_id=conf.GetId())
            ff.Minimize()

    def get_conformer_energies(self, mol):
        """
        Calculate conformer energies.
        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        Returns
        -------
        energies : array_like
            Minimized conformer energies.
        """
        energies = []
        for conf in mol.GetConformers():
            ff = self.get_molecule_force_field(mol, conf_id=conf.GetId())
            energy = ff.CalcEnergy()
            energies.append(energy)
        energies = np.asarray(energies, dtype=float)
        return energies

    def prune_conformers(self, mol):
        """
        Prune conformers from a molecule using an RMSD threshold, starting
        with the lowest energy conformer.
        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        Returns
        -------
        A new RDKit Mol containing the chosen conformers, sorted by
        increasing energy.
        """
        if self.rmsd_threshold < 0 or mol.GetNumConformers() <= 1:
            return mol
        energies = self.get_conformer_energies(mol)
        rmsd = self.get_conformer_rmsd(mol)

        sort = np.argsort(energies)  # sort by increasing energy
        keep = []  # always keep lowest-energy conformer
        discard = []
        for i in sort:

            # always keep lowest-energy conformer
            if len(keep) == 0:
                keep.append(i)
                continue

            # discard conformers after max_conformers is reached
            if len(keep) >= self.max_conformers:
                discard.append(i)
                continue

            # get RMSD to selected conformers
            this_rmsd = rmsd[i][np.asarray(keep, dtype=int)]

            # discard conformers within the RMSD threshold
            if np.all(this_rmsd >= self.rmsd_threshold):
                keep.append(i)
            else:
                discard.append(i)

        # create a new molecule to hold the chosen conformers
        # this ensures proper conformer IDs and energy-based ordering
        new = Chem.Mol(mol)
        new.RemoveAllConformers()
        conf_ids = [conf.GetId() for conf in mol.GetConformers()]
        for i in keep:
            conf = mol.GetConformer(conf_ids[i])
            new.AddConformer(conf, assignId=True)
        return new

    @staticmethod
    def get_conformer_rmsd(mol):
        """
        Calculate conformer-conformer RMSD.
        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        rmsd = np.zeros((mol.GetNumConformers(), mol.GetNumConformers()),
                        dtype=float)
        for i, ref_conf in enumerate(mol.GetConformers()):
            for j, fit_conf in enumerate(mol.GetConformers()):
                if i >= j:
                    continue
                rmsd[i, j] = AllChem.GetBestRMS(mol, mol, ref_conf.GetId(),
                                                fit_conf.GetId())
                rmsd[j, i] = rmsd[i, j]
        return rmsd
