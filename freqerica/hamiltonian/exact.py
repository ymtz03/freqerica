from openfermion import FermionOperator
import openfermion
from itertools import product

from ..circuit.trotter import TrotterStep
from ..op import symm


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
    def __init__(self, orbprop, orb_qubit_map=None, symm_eigval_map=None):
        ham_fop = construct_exact_ham(orbprop)
        ham_qop = openfermion.jordan_wigner(ham_fop)
        n_site  = openfermion.count_qubits(ham_qop)

        if symm_eigval_map is not None:
            symm_paulistr = symm.symmetry_pauli_string(orbprop, operation_list=list(symm_eigval_map))
            for symmop, qop in symm_paulistr.items():
                print('[ham_qop, {:5s}] = [ham_qop, {}] = {}'.format(symmop, str(qop), openfermion.commutator(ham_qop, qop)))

            self.remover = symm.SymmRemover(n_site, list(symm_paulistr.values()))
            self.remover.set_eigvals([symm_eigval_map[symm] for symm in symm_paulistr])
            ham_qop_tapered = self.remover.remove_symm_qubits(ham_qop)
            print(self.remover)
            print('len(ham         ) = ', len(ham_qop.terms))
            print('len(ham[tapered]) = ', len(ham_qop_tapered.terms))

            ham_qop = ham_qop_tapered

        nterm = [0]*(n_site+1) # number of paulistrs which length is n (n=0,1,..,n_site).
        for k in ham_qop.terms: nterm[len(k)] += 1

        self.ham_fop = ham_fop
        self.ham_qop = ham_qop
        self.n_site = n_site
        self.nterm = nterm

    def set_state(self, state):
        if hasattr(self, 'remover'):
        
            from ..circuit.symm import SymmRemoveClifford
            symm_remove_circuits = SymmRemoveClifford(self.n_site, self.remover)
    
            from ..op.util import paulistr
            state_symm_removed = state.copy()
            for symm_remove_circuit in symm_remove_circuits.circuit_list:
                symm_remove_circuit.update_quantum_state(state_symm_removed)
            #for sd, coeff in enumerate(state_symm_removed.get_vector()):
            #    if abs(coeff)<1e-10: continue
            #    print('{:010b} : {:+.3e}'.format(sd, coeff))

            statevec = state_symm_removed.get_vector().reshape([2]*(self.n_site))    
            n_qubit_tapered = self.n_site
            slices = [slice(None)]*self.n_site
            for qop_tgtpauli in self.remover.targetpauli_qop_list:
                index_tgtpauli = paulistr(qop_tgtpauli)[0][0]
                slices[-(index_tgtpauli+1)] = 0
                n_qubit_tapered -= 1
        
            slices_debug = [slice(None)]*self.n_site
            for qop_tgtpauli in self.remover.targetpauli_qop_list:
                index_tgtpauli = paulistr(qop_tgtpauli)[0][0]
                slices_debug[-(index_tgtpauli+1)] = 1
            print(slices)
            #print(slices_debug)
            #for sd, coeff in enumerate(statevec[slices].reshape(-1)):
            #    if abs(coeff)<1e-10: continue
            #    print('{:010b} : {:+.3e}'.format(sd, coeff))

            import numpy as np
            print(np.dot(statevec[slices].reshape(-1), statevec[slices_debug].reshape(-1)))
            #print(statevec[slices]-statevec[slices_debug])
            #assert np.max(abs(statevec[slices]-statevec[slices_debug])) < 1e-10

            from ..circuit.universal import prepstate
            statevec_tapered = statevec[slices].reshape(-1)
            statevec_tapered /= np.linalg.norm(statevec_tapered)
            state_tapered = prepstate(n_qubit_tapered, statevec_tapered)

            state = state_tapered

        self.state = state.copy()
        

    def init_propagate(self, dt):
        self.trotterstep = TrotterStep(self.n_site, -1j * dt * self.ham_qop)

    def propagate(self, state):
        self.trotterstep._circuit.update_quantum_state(state)
