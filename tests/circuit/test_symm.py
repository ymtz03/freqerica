# -*- coding: utf-8 -*-

import unittest
import freqerica.circuit.symm


class BasicTestSuite(unittest.TestCase):
    """Basic test cases."""

    def test_SymmRemoveClifford(self):
        n_qubit = 4
        dim = 2**n_qubit

        from qulacs import QuantumState
        from qulacs.state import inner_product

        from freqerica.circuit.universal import prepstate
        civec = {0b0011:+0.5,
                 0b0101:+0.5j,
                 0b1001:-0.5j,
                 0b0110:-0.5}
        s0 = prepstate(n_qubit, civec)
        print('s0:', s0)

        from openfermion import QubitOperator
        ZZZZ = QubitOperator('Z0 Z1 Z2 Z3')

        from freqerica.operator.symm import SymmRemover
        remover = SymmRemover(n_qubit, [ZZZZ])
        print(remover)
        
        from freqerica.circuit.symm import SymmRemoveClifford
        cliffords = SymmRemoveClifford(n_qubit, remover)
        s1 = s0.copy()
        cliffords.circuit_list[0].update_quantum_state(s1)
        print('s1:', s1)

        import numpy as np
        s1vec = s1.get_vector()
        s1vec = s1vec.reshape(*(2,2,2,2))[:,:,:,0].reshape(-1)
        s1vec /= np.linalg.norm(s1vec)
        s2 = prepstate(n_qubit-1, s1vec)
        print('s2:', s2)
        
        assert True


if __name__ == '__main__':
    unittest.main()
