import numpy as np
from global_qulacs import rotation_factor

def decomp_unitary(umat):
    N = umat.shape[0]
    P = np.arange(N)
    P = P^(P>>1)

    # eliminate (P[N-1-j], P[i]) element by rotate P[N-1-j]'th row and P[N-2-j]'th row.
    operations = []
    for i in range(N-1):
        for j in range(N-1-i):
            a = umat[P[N-2-j]]
            c = umat[P[N-1-j]]

            angle = np.pi if a[P[i]]==0 else np.arctan(-c[P[i]]/a[P[i]])
            operations.append((P[N-2-j], P[N-1-j], angle))

            a_new = a * np.cos(angle) - c * np.sin(angle)
            c_new = a * np.sin(angle) + c * np.cos(angle)

            umat[P[N-2-j]] = a_new
            umat[P[N-1-j]] = c_new

    return operations, np.diag(umat)

def make_umat_by_operations(N, operations):
    mats = []
    for r1, r2, angle in reversed(operations):
        m = np.eye(N)
        m[r1,r1] = m[r2,r2] = np.cos(angle)
        m[r1,r2] =  np.sin(angle)
        m[r2,r1] = -np.sin(angle)
        mats.append(m)

    return mats

def random_umat(N):
    m = np.random.rand(N*N).reshape(N,N)
    for i in range(N):
        m[i] /= np.dot(m[i], m[i])**0.5
        for j in range(i+1, N):
            m[j] -= np.dot(m[i], m[j]) * m[i]

    return m

def toffoli(circuit, c1, c2, tgt):
    circuit.add_H_gate(tgt)
    circuit.add_CNOT_gate(c2, tgt)
    circuit.add_Tdag_gate(tgt)
    circuit.add_CNOT_gate(c1, tgt)
    circuit.add_T_gate(tgt)
    circuit.add_CNOT_gate(c2, tgt)
    circuit.add_Tdag_gate(tgt)
    circuit.add_CNOT_gate(c1, tgt)
    circuit.add_T_gate(c2)
    circuit.add_T_gate(tgt)
    circuit.add_CNOT_gate(c1, c2)
    circuit.add_H_gate(tgt)
    circuit.add_T_gate(c1)
    circuit.add_Tdag_gate(c2)
    circuit.add_CNOT_gate(c1, c2)
    
def control_m_gate(circuit, n_qubits, ctrs, tgt):
    assert len(ctrs)>=1 and (len(ctrs)<=2 or len(ctrs)<=n_qubits-2)
    assert tgt not in ctrs
    
    m = len(ctrs)

    if m==1:
        circuit.add_CNOT_gate(ctrs[0], tgt)
        #print('CNOT({},{})'.format(ctrs[0], tgt))
        return

    if m==2:
        toffoli(circuit, ctrs[0], ctrs[1], tgt)
        #print('Toffoli({},{},{})'.format(ctrs[0], ctrs[1], tgt))
        return

    if m<=n_qubits//2+n_qubits%2:
        ancillae=[]
        for i in range(n_qubits):
            if i not in ctrs and i != tgt:
                ancillae.append(i)
                if len(ancillae)==m-2: break

        toffoli(circuit, ctrs[-1], ancillae[-1], tgt)
        for i in range(m-3):
            toffoli(circuit, ctrs[-2-i], ancillae[-2-i], ancillae[-1-i])
        toffoli(circuit, ctrs[0], ctrs[1], ancillae[0])
        for i in reversed(range(m-3)):
            toffoli(circuit, ctrs[-2-i], ancillae[-2-i], ancillae[-1-i])
        toffoli(circuit, ctrs[-1], ancillae[-1], tgt)

        for i in range(m-3):
            toffoli(circuit, ctrs[-2-i], ancillae[-2-i], ancillae[-1-i])
        toffoli(circuit, ctrs[0], ctrs[1], ancillae[0])
        for i in reversed(range(m-3)):
            toffoli(circuit, ctrs[-2-i], ancillae[-2-i], ancillae[-1-i])

        return

    else:
        for i in range(n_qubits):
            if i not in ctrs and i != tgt:
                ancilla = i
                break

        m1 = m-m//2
        control_m_gate(circuit, n_qubits, ctrs[:m1], ancilla)
        control_m_gate(circuit, n_qubits, ctrs[m1:]+[ancilla], tgt)
        control_m_gate(circuit, n_qubits, ctrs[:m1], ancilla)
        control_m_gate(circuit, n_qubits, ctrs[m1:]+[ancilla], tgt)
        return
    
def control_n_RY_gate(circuit, n_qubits, ctrs, tgt, angle):
    assert tgt not in ctrs
    
    if len(ctrs)==0:
        circuit.add_RY_gate(tgt,  angle * rotation_factor)

    elif len(ctrs)<=2 or len(ctrs) < n_qubits-1:
        circuit.add_RY_gate(tgt,  angle/2 * rotation_factor)
        control_m_gate(circuit, n_qubits, ctrs, tgt)
        circuit.add_RY_gate(tgt, -angle/2 * rotation_factor)
        control_m_gate(circuit, n_qubits, ctrs, tgt)
    
    else:
        control_n_RY_gate(circuit, n_qubits, ctrs[:1], tgt,  angle/2)
        control_m_gate(circuit, n_qubits, ctrs[1:], tgt)
        control_n_RY_gate(circuit, n_qubits, ctrs[:1], tgt, -angle/2)
        control_m_gate(circuit, n_qubits, ctrs[1:], tgt)
        
def arbitrary_unitary_gate(circuit, n_qubits, umat):
    operations, diag = decomp_unitary(umat.copy())
    qubits = list(range(n_qubits))

    for r1, r2, angle in reversed(operations):
        tgt_bin = r1^r2
        if r1 & tgt_bin : angle *= -1 ##

        for i in range(n_qubits):
            if (tgt_bin >>i) & 1 : tgt=i
            elif not ((r1 >> i) & 1):
                circuit.add_X_gate(i)

        ctrs = qubits.copy()
        ctrs.remove(tgt)

        control_n_RY_gate(circuit, n_qubits, ctrs, tgt, angle)

        for i in ctrs:
            if not ((r1 >> i) & 1):
                circuit.add_X_gate(i)


#######################################################################

def highest_order_pos(n):
    assert n>0
    retval=1
    while(retval <=n): retval = retval<<1
    return retval>>1

def get_ancestors(civec):
    ancestors = []
    child=civec
    while True:
        parent={}
        for k,v in child.items():
            key_parent = k-highest_order_pos(k)
            parent.setdefault(key_parent, 0)
            parent[key_parent] += v**2
        
        for k,v in parent.items():
            parent[k] = v**0.5
        
        ancestors.append(parent)
        if len(parent)==1 and 0 in parent: break
        child=parent
    return ancestors

def prepare_civec_circuit(circuit, n_qubits, civec, orb_qubit_map=None):
    if orb_qubit_map==None:
        orb_qubit_map = list(range(n_qubits))
    
    all_components = civec.copy()
    for a in get_ancestors(civec):
        all_components.update(a)
        if len(a)==1: key_root = list(a.keys())[0]
    
    for k_child in sorted(all_components.keys()):
        if k_child==key_root: continue
            
        tgt_bin = highest_order_pos(k_child)
        k_parent = k_child - tgt_bin
        c_parent = all_components[k_parent]
        c_child = all_components[k_child]
        
        angle = np.arcsin(min(max(c_child/c_parent,-1.0),1.0))
        #print('{:2d} {:2d} {} {} {} {}'.format(k_child, k_parent, c_child, c_parent, c_child/c_parent, angle))

        for i in range(n_qubits):
            if tgt_bin>>i & 1 :
                tgt=i
                break
        ctrs=list(range(tgt))
            
        for i in ctrs:
            if k_parent>>i & 1==0: circuit.add_X_gate(orb_qubit_map[i])
            
        #print(ctrs, tgt)
        control_n_RY_gate(circuit, n_qubits, [orb_qubit_map[i] for i in ctrs], orb_qubit_map[tgt], angle)
            
        for i in ctrs:
            if k_parent>>i & 1==0: circuit.add_X_gate(orb_qubit_map[i])
            
        all_components[k_parent] = max((c_parent**2 - c_child**2),0)**0.5
        
    #print(all_components)
