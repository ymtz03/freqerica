# -*- coding: utf-8 -*-

import unittest
import freqerica.operator.symm


class BasicTestSuite(unittest.TestCase):
    """Basic test cases."""

    def test_reduce_qubits(self):
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

        print()
        n_qubit=4
        ham  = XXII + XYII*2 + YXII*3 + YYII*4 + IIXX*5 + IIXY*6 + IIYX*7 + IIYY*8 + XXXX*9 + YYYY*10
        symmetry_list = [ZZZZ, ZZII, IIZZ]
        remover = freqerica.operator.symm.SymmRemover(n_qubit, symmetry_list)
        print(remover)
        ham_c = remover.operate_clifford(ham)
        ZZZZ_c = remover.operate_clifford(ZZZZ)
        ZZII_c = remover.operate_clifford(ZZII)
        IIZZ_c = remover.operate_clifford(IIZZ)

        print('[H, ZZZZ] = ', commutator(ham, ZZZZ))
        print('[H, ZZII] = ', commutator(ham, ZZII))
        print('[H, IIZZ] = ', commutator(ham, IIZZ))
        print('U(ZZZZ)U = ', remover.operate_clifford(ZZZZ))
        print('U(ZZII)U = ', remover.operate_clifford(ZZII))
        print('U(IIZZ)U = ', remover.operate_clifford(IIZZ))
        print('H = \n', ham)
        print('UHU = \n', ham_c)
        print('[UHU, U(ZZZZ)U] = ', commutator(ham_c, ZZZZ_c))
        print('[UHU, U(ZZII)U] = ', commutator(ham_c, ZZII_c))
        print('[UHU, U(IIZZ)U] = ', commutator(ham_c, IIZZ_c))

        remover.set_eigvals([+1,-1,+1])
        ham_cr = remover.replace_pauli_with_eigenvalue(ham_c)
        print('UHU(replaced) = \n', ham_cr)

        ham_crt = remover.taper_qubits_off(ham_cr)
        print('UHU(tapered) = \n', ham_crt)

        ham_remove = remover.remove_symm_qubits(ham)
        print('H_removed = \n', ham_crt)
        
        
        assert True


if __name__ == '__main__':
    unittest.main()
