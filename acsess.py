#!/usr/bin/env python
#-*- coding: utf-8 -*-
'''
This is the main function of the ACSESS
'''
debug = False
##############################
# Import statements
##############################
import sys, random
from rdkit import Chem
sys.path.append('.')
import mprms
import init
import os
import subprocess
import drivers as dr
import output
from output import stats
import objective
from helpers import DumpMols, FinishSelection
from distance import AveNNDistance
from similarity import NNSimilarity

# set global variables:
_iterhead = "\n-------------------- Iteration {0} ----------------\n"
_gridAssign = None
writeInterval=1

def initiate():
    global startiter, lib, pool, _gridAssign, cd, writeInterval
    ##################################
    # 1. read input and initialize
    ##################################
    init.Initialize()
    ##################################
    # 2. Get starting library
    ##################################
    startiter, lib, pool = init.StartLibAndPool(mprms.restart)

    if mprms.cellDiversity:
        import celldiversity as cd
        siml,mylib,_gridAssign = cd.GridDiversity( [], lib+pool )
        #wrotepool=set( m.GetData('isosmi') for m in mylib+pool )

    # this file is load before initialization so we have to change it here:
    if hasattr(mprms, 'writeInterval'):
        writeInterval = mprms.writeInterval
    return


def evolve():
    global startiter, lib, pool, _iterhead, _gridAssign

    ###################################################
    ##########                              ###########
    ##########           MAIN LOOP          ###########
    ##########                              ###########
    ###################################################
    for gen in xrange(startiter, mprms.nGen):
        # 0. Decide workflow for current iteration:
        Tautomerizing, Filtering, GenStrucs = dr.SetIterationWorkflow(gen)

        # 1. PRELOGGING
        print _iterhead.format(gen)
        stats.update({'gen': gen, 'nPool': len(pool), 'nLib': len(lib)})

        # 2.MUTATIONS AND CROSSOVERS
        newlib = dr.DriveMutations(lib)

        # 3. FILTERS
        newlib = dr.DriveFilters(newlib, Filtering, GenStrucs)

        # 3b. Once Filtering/GenStrucs is switched on Pool has to be filtered too:
        pool = dr.DrivePoolFilters(pool, Filtering, GenStrucs, Tautomerizing,
                                   gen)

        # 4. OBJECTIVE EVALUATION
        if mprms.optimize:
            print "acsess.py len(newlib, pool):{} {}".format(
                len(newlib), len(pool))
            pool = objective.EvaluateObjective(newlib + pool, gen)
        else:
            pool = dr.ExtendPool(pool, lib, newlib)

        # 5. SELECTION
        siml = None
        if mprms.cellDiversity:
            if mprms.optimize:
                siml, lib, _gridAssign = cd.GridDiversity(lib, newlib, molgrid=_gridAssign)
            else:
                siml, lib, _gridAssign = cd.GridDiversity_JITFilter(lib, newlib, molgrid=_gridAssign)
        else:
            if mprms.optimize:
                oldN = len(pool)
                lib, pool = objective.SelectFittest(pool, mprms.subsetSize, gen)
                stats['nUnFit'] = oldN - len(pool)
            elif len(pool) > mprms.subsetSize:
                lib = dr.DriveSelection(pool, mprms.subsetSize)
            else:
                lib = [mol for mol in pool]
        FinishSelection(lib)
        if len(lib) == 0:
            raise RuntimeError('no molecules left')

        # 6. DIVERSITY IS:
        if siml is None:
            if mprms._similarity:
                siml = NNSimilarity(lib)
            else:
                # library diversity should not be assessed by normalization of only
                # the lib. so either no normalization or a normalization based on 
                # the whole pool.
                siml = AveNNDistance(lib, norm=False)
        print '\nLIBRARY DIVERSITY: ', siml

        # 7. POSTLOGGING
        with open('mylib.smi', 'w') as f:
            for i, mol in enumerate(lib):
                f.write(Chem.MolToSmiles(mol) + ' {:d}\n'.format(i))
        if gen % writeInterval == 0 or gen == mprms.nGen - 1:
            print "writeInterval:", writeInterval
            DumpMols(lib, gen)
        DumpMols(pool)
        stats['diversity'] = siml
        output.PrintTimings()
        output.PrintStat()
        print "processing time:",
        print [ p for p in subprocess.check_output(['ps','-fu','vsc10010']).split('\n')[:-1] if str(os.getpid()) in p.split()[1] ][0].split()[6]

    output.PrintTotalTimings()
    print "DONE"
    return


if __name__ == "__main__":
    import sys
    import signal
    import inspect

    def signal_term_handler(signal, frame):
        if signal or frame:
            print 'got kill signal!'
        output.PrintTimings()
        output.PrintStat()
        output.PrintTotalTimings()
        try:
            f_locals = inspect.trace()[-1][0].f_locals
            DumpMols(f_locals['lib'], f_locals['gen'])
            DumpMols(f_locals['pool'])
        except Exception as e:
            print "didn't manage to dump pool and mylib after kill signal", e
            print "Error Termination"
        sys.exit(0)

    signal.signal(signal.SIGTERM, signal_term_handler)

    class Unbuffered(object):
        def __init__(self, stream):
            self.stream = stream

        def writelines(self, datas):
            self.stream.writelines(datas)
            self.stream.flush()

        def write(self, data):
            self.stream.write(data)
            self.stream.flush()

        def __getattr__(self, attr):
            return getattr(self.stream, attr)

    sys.stdout = Unbuffered(sys.stdout)

    class RunACSESS(object):
        def __init__(self):
            initiate()

        def evolve(self):
            try:
                evolve()
            except KeyboardInterrupt:
                print "\n\n\t\t##############################\n"+\
                          "\t\t# Catched KeyboardInterrupt! #\n"+\
                          "\t\t##############################"
                signal_term_handler(None, None)
            else:
                print "Normal Termination"

    run = RunACSESS()
    run.evolve()
