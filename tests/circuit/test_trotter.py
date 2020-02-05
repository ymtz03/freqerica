# -*- coding: utf-8 -*-

import unittest
import freqerica.circuit.trotter


class BasicTestSuite(unittest.TestCase):
    """Basic test cases."""

    def test_trotterstep_circuit(self):
        from freqerica.op.symbol import WrappedExpr as Symbol
        from openfermion import FermionOperator, jordan_wigner, get_sparse_operator
        from sympy import Array
        import numpy as np
        np.random.seed(100)
        import scipy
        import qulacs

        n_orb = 2

        const = Symbol('const')
        T1 = [[None for _ in range(n_orb)] for _ in range(n_orb)]
        for p in range(n_orb):
            for q in range(p, n_orb):
                t = Symbol('t{}{}'.format(p,q))
                T1[p][q] = T1[q][p] = t
        T1 = Array(T1)
        print(T1)

        const_value = np.random.rand()
        T1_value = np.random.rand(n_orb,n_orb)*0.01
        T1_value += T1_value.T
        print(const_value)
        print(T1_value)

        def op1e(const, Tmat):
            fop = FermionOperator('', const)
            for p in range(n_orb):
                for q in range(n_orb):
                    fop += FermionOperator( ((2*p  , 1),(2*q  , 0)), Tmat[p,q] )
                    fop += FermionOperator( ((2*p+1, 1),(2*q+1, 0)), Tmat[p,q] )
            return fop

        fop_symbol = op1e(const, T1)
        qop_symbol = jordan_wigner(fop_symbol)
        print(fop_symbol)
        print(qop_symbol)

        n_qubit = n_orb*2

        # c1 : TrotterStep with symbolic qop
        c1 = freqerica.circuit.trotter.TrotterStep(n_qubit, (-1j)*qop_symbol)
        print(c1)
        symbol_number_pairs = [(const, const_value), (T1, T1_value)]
        c1.subs(symbol_number_pairs)
        print(c1)

        # c2 : TrotterStep with numerical qop
        fop_number = op1e(const_value, T1_value)
        qop_number = jordan_wigner(fop_number)
        c2 = freqerica.circuit.trotter.TrotterStep(n_qubit, (-1j)*qop_number)
        print(c2)

        # c3 : oldtype with numerical qop
        c3 = qulacs.QuantumCircuit(n_qubit)
        freqerica.circuit.trotter.trotter_step_2nd_order(c3, (-1j)*qop_number)



        s0 = qulacs.QuantumState(n_qubit)
        s0.set_Haar_random_state()
        s1 = s0.copy()
        s2 = s0.copy()
        s3 = s0.copy()
        from freqerica.util.qulacsnize import convert_state_vector
        sv = convert_state_vector( n_qubit, s0.get_vector() )

        corr1 = []
        corr2 = []
        corr3 = []
        corrv = []
        for t in range(100):
            corr1.append( qulacs.state.inner_product(s0, s1) )
            corr2.append( qulacs.state.inner_product(s0, s2) )
            corr3.append( qulacs.state.inner_product(s0, s3)*np.exp(-1j*qop_number.terms[()]*t) )
            corrv.append( np.dot(np.conjugate(convert_state_vector(n_qubit, s0.get_vector())), sv) )
            c1._circuit.update_quantum_state(s1)
            c2._circuit.update_quantum_state(s2)
            c3.update_quantum_state(s3)
            sv = scipy.sparse.linalg.expm_multiply(-1j * get_sparse_operator(qop_number), sv)

        #import matplotlib.pyplot as plt
        #plt.plot(np.array(corr1).real, label='1')
        #plt.plot(np.array(corr2).real, label='2')
        #plt.plot(np.array(corr3).real, label='3')
        #plt.plot(np.array(corrv).real, label='4')
        #plt.legend()
        #plt.show()

        # print('s1', s1.get_vector())
        # print('s2', s2.get_vector())
        # print('s3', s3.get_vector())

if __name__ == '__main__':
    unittest.main()
