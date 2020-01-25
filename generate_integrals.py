import numpy as np
import scipy
import sys

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

