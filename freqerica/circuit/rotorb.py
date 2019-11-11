import qulacs
import numpy as np
from freqerica.context import rotation_factor

def decomp_unitary(umat):
    """ decompose unitary matrix U to 2x2 rotation matrices {R_k} """
    N = umat.shape[0]
    umat = umat + 0j  # deepcopy umat and change dtype to complex
    steps = []
    
    for istep in range(2*N-3):
        rotations = []
        
        icol_tgt = max(0, istep-N+2)
        irow_tgt = N - istep + icol_tgt*2 - 1

        while irow_tgt < N:
            irow_eliminator = irow_tgt-1
            
            if abs(umat[irow_eliminator, icol_tgt])!=0:
                delta = 0.5*( np.angle(umat[irow_eliminator, icol_tgt]) - np.angle(umat[irow_tgt, icol_tgt]) )
                angle = -np.arctan(abs(umat[irow_tgt, icol_tgt])/abs(umat[irow_eliminator, icol_tgt]))
                #angle = -np.arctan(umat[irow_tgt, icol_tgt].real/umat[irow_eliminator, icol_tgt].real)
            else:
                delta = 0
                angle = np.pi/2
                
            rotations.append((irow_tgt, delta, angle))

            row_eliminator_new = umat[irow_eliminator] * np.exp(-1j*delta) * np.cos(angle) - umat[irow_tgt] * np.exp(+1j*delta) * np.sin(angle)
            row_tgt_new        = umat[irow_eliminator] * np.exp(-1j*delta) * np.sin(angle) + umat[irow_tgt] * np.exp(+1j*delta) * np.cos(angle)
            umat[irow_eliminator] = row_eliminator_new
            umat[irow_tgt]        = row_tgt_new

            icol_tgt+=1
            irow_tgt+=2


        steps.append(rotations)

    phases_diag = np.angle(np.diag(umat))
        
    return steps, phases_diag



class OrbitalRotation:
    """
    Examples:
        .. code-block:: python

            or = OrbitalRotation(umat)
            or._circuit.update_quantum_state(state)
    """

    
    def __init__(self, umat):
        self._n_qubit = umat.shape[0]
        self._circuit = qulacs.ParametricQuantumCircuit(self._n_qubit)
        
        steps, phases_diag = decomp_unitary(umat)

        # 1-qubit phase rotation layer
        for i in range(self._n_qubit):
            self._circuit.add_RZ_gate(i, -0.5*phases_diag[i]*rotation_factor)

        # 2-qubit operations
        for step in reversed(steps):
            for irow_tgt, delta, angle in reversed(step):
                target = [irow_tgt-1, irow_tgt]
                self._circuit.add_parametric_multi_Pauli_rotation_gate(target, [1, 2], +0.5*angle*rotation_factor) # X[p-1] Y[p]
                self._circuit.add_parametric_multi_Pauli_rotation_gate(target, [2, 1], -0.5*angle*rotation_factor) # Y[p-1] X[p]
                self._circuit.add_RZ_gate(irow_tgt-1, -0.5*delta*rotation_factor)
                self._circuit.add_RZ_gate(irow_tgt  , +0.5*delta*rotation_factor)

    def __str__(self):
        return str(self._circuit)

    

######## codes below are for debugging

def get_rotation_matrices(N, operations):
    mats = []
    for ir, delta, angle in operations:
        ms = np.eye(N, dtype=complex)
        ms[ir  ,ir  ] = np.exp(+1j*delta)
        ms[ir-1,ir-1] = np.exp(-1j*delta)
        mats.append(ms)
        
        mr = np.eye(N)
        mr[ir,ir] = mr[ir-1,ir-1] = np.cos(angle)
        mr[ir  ,ir-1] = +np.sin(angle)
        mr[ir-1,ir  ] = -np.sin(angle)
        mats.append(mr)

    return mats

def get_random_umat(N):
    m = np.random.rand(N*N).reshape(N,N)
    for i in range(N):
        m[i] /= np.dot(m[i], m[i])**0.5
        for j in range(i+1, N):
            m[j] -= np.dot(m[i], m[j]) * m[i]

    return m

def get_random_antihermite_mat(N):
    m  = np.random.rand(N*N).reshape(N,N)*1j
    m += np.random.rand(N*N).reshape(N,N)
    for i in range(N):
        m[i,i] = m[i,i].imag * 1j
        for j in range(N):
            m[i,j] = -np.conjugate(m[j,i])
    return m


def test_decomp_unitary():
    import scipy.linalg
    from pprint import pprint
    np.set_printoptions(linewidth=300)

    N=6
    #umat = get_random_umat(N)
    k = get_random_antihermite_mat(N)
    umat = scipy.linalg.expm(k)
    steps, phases_diag = decomp_unitary(umat)
    pprint(steps)

    operations = []
    for step in steps: operations += list(step)
    mats = get_rotation_matrices(N, operations)

    print('kappa\n', k, sep='')
    print('umat\n', umat, sep='')
    print('dot(U,U.T)\n', umat.dot(np.conjugate(umat.T)),sep='')
    for m in mats:
        umat = np.dot(m, umat)
    print('r...ru\n', umat, sep='')
    print('exp(phase)\n', np.exp(1j*phases_diag), sep='')
    print('diff(phase)\n', np.exp(1j*phases_diag)-np.diag(umat), sep='')

    
def test_rotorb_circuit():
    import scipy.linalg
    np.set_printoptions(linewidth=300)

    N = 8

    #k = np.random.rand(N,N)
    #k -= (k+k.T)/2  # random real antisymmetric matrix
    k = get_random_antihermite_mat(N)
    umat = scipy.linalg.expm(k)
    print(umat)

    c1 = OrbitalRotation(umat)

    from openfermion import FermionOperator, jordan_wigner, get_sparse_operator
    fop = FermionOperator()
    for p in range(N):
        for q in range(N):
            fop += FermionOperator( ((p,1),(q,0)), k[p,q] )
    qop = jordan_wigner(fop)
    #print(qop)

    from trotter_qulacs import TrotterStep
    from util_qulacs import convert_state_vector
    M = 100
    c2 = TrotterStep(N, qop/M)

    s1 = qulacs.QuantumState(N)
    s1.set_Haar_random_state()
    s2 = s1.copy()
    s3 = convert_state_vector(N, s1.get_vector())

    c1._circuit.update_quantum_state(s1)
    for i in range(M): c2._circuit.update_quantum_state(s2)

    s3 = scipy.sparse.linalg.expm_multiply(get_sparse_operator(qop), s3)

    ip12 = qulacs.state.inner_product(s1, s2)
    ip13 = np.conjugate(convert_state_vector(N, s1.get_vector())).dot(s3)
    ip23 = np.conjugate(convert_state_vector(N, s2.get_vector())).dot(s3)
    print('dot(s1,s2)', abs(ip12), ip12)
    print('dot(s1,s3)', abs(ip13), ip13)
    print('dot(s2,s3)', abs(ip23), ip23)
    
if __name__=="__main__":
    #test_decomp_unitary()
    test_rotorb_circuit()
