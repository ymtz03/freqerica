# -*- coding: utf-8 -*-

import unittest
import freqerica.circuit.universal


class BasicTestSuite(unittest.TestCase):
    """Basic test cases."""

    def test_prepstate(self):
        n_qubit = 4
        dim = 2**n_qubit

        from qulacs import QuantumState
        from qulacs.state import inner_product
        s0 = QuantumState(n_qubit)
        s0.set_Haar_random_state()

        s1 = freqerica.circuit.universal.prepstate(n_qubit, s0.get_vector())
        print(inner_product(s0, s1))        


        civec = {0b0011:.5, 0b0101:+.5j, 0b1001:-.5j, 0b0110:-.5}
        s2 = freqerica.circuit.universal.prepstate(n_qubit, civec)
        print(s2)
        
        assert True


if __name__ == '__main__':
    unittest.main()
