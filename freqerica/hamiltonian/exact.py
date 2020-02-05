from openfermion import FermionOperator, jordan_wigner, count_qubits
from itertools import product

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


def construct_exact_ham(orbprop, orb_qubit_map=None):
    ham = FermionOperator.identity() * orbprop.hcore
    ham += construct_ham_1body(orbprop.hint, orb_qubit_map)
    ham += construct_ham_2body(orbprop.gint, orb_qubit_map)
    return ham


class ExactHam:
    def __init__(self, orbprop, orb_qubit_map=None):
        self.ham_fop = construct_exact_ham(orbprop)
        self.ham_qop = jordan_wigner(self.ham_fop)
        self.n_site  = count_qubits(self.ham_qop)
