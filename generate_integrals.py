import numpy as np
import scipy
import sys

import pyscf
from pyscf import gto, scf, ao2mo, fci, cc
from pyscf_helper import *

r0 = 1.0977
molecule = '''
N      0.00       0.00       0.00
N      0.00       0.00       {}'''.format(r0)

charge = 0
spin  = 0
basis_set = 'ccpvdz'
basis_set = 'sto-3g'

orb_basis = 'scf'
cas = True
cas_nstart = 2
cas_nstop = 28
cas_nel = 10
loc_nstart = 2
loc_nstop = 10

pmol = PyscfHelper()
#pmol.init(molecule,charge,spin,basis_set,orb_basis,cas,cas_nstart,cas_nstop,cas_nel,loc_nstart,loc_nstop)
pmol.init(molecule,charge,spin,basis_set,orb_basis)

np.save('ints_0b.npy',pmol.ecore)
np.save('ints_1b.npy',pmol.h)
np.save('ints_2b.npy',pmol.g)

test = True
if test:
    pyscf.lib.num_threads(1)
    mol = gto.M(verbose=3)
    mol.nelectron = 6 
    mol.incore_anyway = True
    cisolver = fci.direct_spin1.FCI(mol)
    #e, ci = cisolver.kernel(h1, eri, h1.shape[1], 2, ecore=mol.energy_nuc())
    e, ci = cisolver.kernel(pmol.h, pmol.g, pmol.h.shape[1], mol.nelectron, ecore=pmol.ecore)
    print(" FCI:        %12.8f"%e)
