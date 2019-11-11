import numpy as np
from pyscf import gto, scf, mcscf, ao2mo
from openfermion import FermionOperator
from itertools import product
from collections import namedtuple

MoleInput = namedtuple('MoleInput', ['mol', 'norb', 'nelec'])

def unfold_compact_h2mat(h2mat, norb):
    retval = np.empty([norb]*4, float)
    r=0
    for r1 in range(norb):
        for r2 in range(r1+1):
            c=0
            for c1 in range(norb):
                for c2 in range(c1+1):
                    retval[r1,r2,c1,c2] = h2mat[r,c]
                    retval[r2,r1,c1,c2] = h2mat[r,c]
                    retval[r1,r2,c2,c1] = h2mat[r,c]
                    retval[r2,r1,c2,c1] = h2mat[r,c]
                    c+=1
            r+=1
    return retval


def calculate_moint_and_energy(moleInput):
    mf = scf.RHF(moleInput.mol)
    mf.max_cycle = 2000
    mf.kernel()
    
    norb = moleInput.norb
    nelec = moleInput.nelec
    mc = mcscf.CASCI(mf, norb, nelec)
    h1, hcore = mc.get_h1eff()
    eri = mc.get_h2eff()
    mc.kernel()

    #mo_coeff = mean_field.mo_coeff
    #hint_spacial = mo_coeff.T.dot(mean_field.get_hcore()).dot(mo_coeff)
    #n_orbs = hint_spacial.shape[0]
    #gint_spacial = ao2mo.full(mol, mo_coeff, compact=False).reshape([n_orbs]*4)
    n_site=norb*2

    hint_spatial = h1
    gint_spatial = unfold_compact_h2mat(eri, norb)

    hint = np.zeros([n_site]*2, float)
    gint = np.zeros([n_site]*4, float)

    hint[0::2, 0::2] = hint_spatial
    hint[1::2, 1::2] = hint_spatial
    gint[0::2, 0::2, 0::2, 0::2] = gint_spatial
    gint[1::2, 1::2, 0::2, 0::2] = gint_spatial
    gint[0::2, 0::2, 1::2, 1::2] = gint_spatial
    gint[1::2, 1::2, 1::2, 1::2] = gint_spatial

    return hcore, hint, gint


def construct_ham_1body(hint, orb_qubit_map=None):
    if orb_qubit_map is None: orb_qubit_map = list(range(hint.shape[0]))
    
    ham_1body = FermionOperator()
    for (p, p_qubit), (q, q_qubit) in product(enumerate(orb_qubit_map), repeat=2):
        ham_1body += FermionOperator(((p_qubit, 1), (q_qubit, 0)), hint[p,q])

    return ham_1body


def construct_ham_2body(gint, orb_qubit_map=None):
    if orb_qubit_map is None: orb_qubit_map = list(range(gint.shape[0]))
    
    ham_2body = FermionOperator()
    for (p,p_qubit), (q,q_qubit), (r,r_qubit), (s,s_qubit) in product(enumerate(orb_qubit_map), repeat=4):
        if p!=r and q!=s:
            ham_2body += FermionOperator(
                ((p_qubit, 1),(r_qubit, 1),(s_qubit, 0),(q_qubit, 0)), gint[p,q,r,s]) * 0.5

    return ham_2body


def construct_exact_ham(moleInput, orb_qubit_map=None):
    hcore, hint, gint = calculate_moint_and_energy(moleInput)
    ham = FermionOperator.identity() * hcore
    ham += construct_ham_1body(hint, orb_qubit_map)
    ham += construct_ham_2body(gint, orb_qubit_map)
    return ham


moleInput_example = {}

mol = gto.Mole()
mol.basis = "sto-6g"
mol.verbose = 0
mol.spin = 0
mol.atom = [["Li", (0.0, 0.0, 0.0)],["H", (0.0, 0.0, 1.5)]]
mol.build()
moleInput_example['LiH'] = MoleInput(mol, norb=5, nelec=2)


if __name__=="__main__":
    r=1.5
    
    mol = gto.Mole()
    mol.basis = "sto-6g"
    mol.verbose = 0
    mol.spin = 0
    mol.atom = [["Li", (0.0, 0.0, 0.0)],["H", (0.0, 0.0, r)]]
    mol.build()

    norb  = 5
    nelec = 2

    moleInput = MoleInput(mol, norb, nelec)
    hcore, hint, gint = calculate_moint_and_energy(moleInput)
    ham = construct_exact_ham(moleInput)

    import openfermion
    print('n_qubit:', openfermion.count_qubits(ham))



