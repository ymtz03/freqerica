# -*- coding: utf-8 -*-

import unittest
import freqerica.circuit.rotorb


def get_random_antihermite_mat(N):
    import numpy as np
    m  = np.random.rand(N*N).reshape(N,N)*1j
    m += np.random.rand(N*N).reshape(N,N)
    for i in range(N):
        m[i,i] = m[i,i].imag * 1j
        for j in range(N):
            m[i,j] = -np.conjugate(m[j,i])
    return m


class BasicTestSuite(unittest.TestCase):
    """Basic test cases."""

    def test_decomp_unitary(self):
        import numpy as np
        import scipy.linalg
        from pprint import pprint
        np.set_printoptions(linewidth=300)
        np.random.seed(100)

        def single_test():
            N=30
            #umat = get_random_umat(N)
            k = get_random_antihermite_mat(N)
            umat = scipy.linalg.expm(k)
            steps, phases_diag = freqerica.circuit.rotorb.decomp_unitary(umat)
            #pprint(steps)

            operations = []
            for step in steps: operations += list(step)
            mats = freqerica.circuit.rotorb.get_rotation_matrices(N, operations)

            #print('kappa\n', k, sep='')
            #print('umat\n', umat, sep='')
            #print('dot(U,U.T)\n', umat.dot(np.conjugate(umat.T)),sep='')
            for m in mats:
                umat = np.dot(m, umat)
            #print('r...ru\n', umat, sep='')
            #print('exp(phase)\n', np.exp(1j*phases_diag), sep='')
            diff = np.exp(1j*phases_diag)-np.diag(umat)
            #print('diff(phase)\n', diff, sep='')
            
            print('[trial{:02d}] maxdiff : {:.1e}'.format(i,max(np.abs(diff))))
            self.assertLess(max(np.abs(diff)), 1e-10)

        print()
        for i in range(10):
            with self.subTest(i=i):
                single_test()

    def test_rotorb_circuit(self):
        import numpy as np
        import scipy.linalg
        np.set_printoptions(linewidth=300)

        def single_test(i):
            N = 8

            #k = np.random.rand(N,N)
            #k -= (k+k.T)/2  # random real antisymmetric matrix
            k = get_random_antihermite_mat(N)
            umat = scipy.linalg.expm(k)
            #print(umat)

            c1 = freqerica.circuit.rotorb.OrbitalRotation(umat)

            from openfermion import FermionOperator, jordan_wigner, get_sparse_operator
            fop = FermionOperator()
            for p in range(N):
                for q in range(N):
                    fop += FermionOperator( ((p,1),(q,0)), k[p,q] )
            qop = jordan_wigner(fop)
            #print(qop)

            #from freqerica.circuit.trotter import TrotterStep
            from freqerica.util.qulacsnize import convert_state_vector
            import qulacs
            #M = 100
            #c2 = TrotterStep(N, qop/M)

            s1 = qulacs.QuantumState(N)
            s1.set_Haar_random_state()
            #s2 = s1.copy()
            s3 = convert_state_vector(N, s1.get_vector())

            c1._circuit.update_quantum_state(s1)
            #for i in range(M): c2._circuit.update_quantum_state(s2)

            s3 = scipy.sparse.linalg.expm_multiply(get_sparse_operator(qop), s3)

            #ip12 = qulacs.state.inner_product(s1, s2)
            ip13 = np.conjugate(convert_state_vector(N, s1.get_vector())).dot(s3)
            #ip23 = np.conjugate(convert_state_vector(N, s2.get_vector())).dot(s3)
            #print('dot(s1,s2)', abs(ip12), ip12)
            #print('dot(s1,s3)', abs(ip13), ip13)
            #print('dot(s2,s3)', abs(ip23), ip23)

            print('[trial{:02d}] fidelity-1 : {:+.1e}   dot=({:+.5f})'.format(i, abs(ip13)-1, ip13))

        print()
        for i in range(10):
            with self.subTest(i=i):
                single_test(i)

            
                
if __name__ == '__main__':
    unittest.main()
