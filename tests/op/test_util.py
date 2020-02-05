# -*- coding: utf-8 -*-

import unittest
import freqerica.op.util


class BasicTestSuite(unittest.TestCase):
    """Basic test cases."""

    def test_paulistr(self):
        from openfermion import QubitOperator, commutator
        n_qubit = 3
        XXX = QubitOperator('X0 X1 X2')
        YYY = QubitOperator('Y0 Y1 Y2')
        YZX = QubitOperator('Y0 Z1 X2')
        ZZI = QubitOperator('Z0 Z1   ')
        IZZ = QubitOperator('   Z1 Z2')
        IXX = QubitOperator('   X1 X2')

        ZZZZ = QubitOperator('Z0 Z1 Z2 Z3')
        ZZII = QubitOperator('Z0 Z1      ')
        IIZZ = QubitOperator('      Z2 Z3')
        XXII = QubitOperator('X0 X1      ')
        XYII = QubitOperator('X0 Y1      ')
        YXII = QubitOperator('Y0 X1      ')
        YYII = QubitOperator('Y0 Y1      ')
        IIXX = QubitOperator('      X2 X3')
        IIXY = QubitOperator('      X2 Y3')
        IIYX = QubitOperator('      Y2 X3')
        IIYY = QubitOperator('      Y2 Y3')
        XXXX = QubitOperator('X0 X1 X2 X3')
        YYYY = QubitOperator('Y0 Y1 Y2 Y3')

        freqerica.op.util.paulistr(ZZZZ)
        freqerica.op.util.paulistr(IXX)
        
        assert True


    def test_listupCSFs(self):
        from openfermion import QubitOperator, commutator
        ZZZZ  = QubitOperator('Z0 Z1 Z2 Z3')
        ZZZZ2 = QubitOperator('Z0 Z1 Z4 Z6')
        
        from freqerica.op.symm import SymmRemover
        norb = 4
        n_qubit = norb*2
        remover = SymmRemover(n_qubit, [ZZZZ, ZZZZ2])
        print(remover)
        remover.set_eigvals([-1, +1])
        print(remover)
        
        wfs = freqerica.op.util.listupCSFs(norb, mult=1, na=2, nb=2, remover=remover)
        #print(wfs)
        
        assert True

    def test_remove(self):
        wf = {0b0110: 1.0} # (1b)(2a)
        #wf = {0b0101: 1.0} # (1a)(2a)
        #wf = {0b1010: 1.0} # (1b)(2b)
        wf_1let = freqerica.op.util.remove(wf, norb=4, mult_tgt=1, mult_rmv=3) # remove 3let
        wf_3let = freqerica.op.util.remove(wf, norb=4, mult_tgt=3, mult_rmv=1) # remove 1let

        f = freqerica.op.util.printwf
            
        print('1let', f(wf_1let), sep='\n')
        print('3let', f(wf_3let), sep='\n')

        wf = {0b10100101: 1.0} # (1a)(2a)(3b)(4b)
        wf_1let = freqerica.op.util.remove(wf     , norb=4, mult_tgt=1, mult_rmv=3) # remove 3let
        wf_1let = freqerica.op.util.remove(wf_1let, norb=4, mult_tgt=1, mult_rmv=5) # remove 5tet
        wf_3let = freqerica.op.util.remove(wf     , norb=4, mult_tgt=3, mult_rmv=1) # remove 1let
        wf_3let = freqerica.op.util.remove(wf_3let, norb=4, mult_tgt=3, mult_rmv=5) # remove 5tet
        wf_5tet = freqerica.op.util.remove(wf     , norb=4, mult_tgt=5, mult_rmv=1) # remove 1let
        wf_5tet = freqerica.op.util.remove(wf_5tet, norb=4, mult_tgt=5, mult_rmv=3) # remove 3tet

        print('1let')
        print(f(wf_1let))

        print('3let')
        print(f(wf_3let))

        print('5let')
        print(f(wf_5tet))
        
        assert True


if __name__ == '__main__':
    unittest.main()
