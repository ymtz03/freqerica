from openfermion import QubitOperator

# symmtable['irrep']['Symmetry Operation'] = sign
NA=None
symmtable = {
    'A1' : {'Rz2':+1, 'sxz':+1, 'syz':+1 }, 
    'A2' : {'Rz2':+1, 'sxz':-1, 'syz':-1 },
    'B1' : {'Rz2':-1, 'sxz':+1, 'syz':+1 },
    'B2' : {'Rz2':-1, 'sxz':-1, 'syz':-1 },
    'E1x': {'Rz2':-1, 'sxz':+1, 'syz':-1 },
    'E1y': {'Rz2':-1, 'sxz':-1, 'syz':+1 },
    'E2x': {'Rz2':+1, 'sxz':NA, 'syz':NA },
    'E2y': {'Rz2':+1, 'sxz':NA, 'syz':NA },
    'E3x': {'Rz2':-1, 'sxz':NA, 'syz':NA },
    'E3y': {'Rz2':-1, 'sxz':NA, 'syz':NA },
    'E4x': {'Rz2':+1, 'sxz':NA, 'syz':NA },
    'E4y': {'Rz2':+1, 'sxz':NA, 'syz':NA },
    'E5x': {'Rz2':-1, 'sxz':NA, 'syz':NA },
    'E5y': {'Rz2':-1, 'sxz':NA, 'syz':NA },
}

def symmetry_pauli_string(orbprop, operation_list):
    operation_paulistr_table = {}
    
    na = QubitOperator(())
    nb = QubitOperator(())
    for i in range(orbprop.ncas):
        na *= QubitOperator((i*2  , 'Z'))
        nb *= QubitOperator((i*2+1, 'Z'))
    if 'na' in operation_list: operation_paulistr_table['na'] = na
    if 'nb' in operation_list: operation_paulistr_table['nb'] = nb

    ncore, ncas = orbprop.ncore, orbprop.ncas
    for op in operation_list:
        if op=='na' or op=='nb': continue
        
        qop_paulistr = QubitOperator(())
        for i, iorb in enumerate(range(ncore, ncore+ncas)):
            irrep = orbprop.irreps[iorb]
            char = symmtable[irrep][op]
            if char==-1: qop_paulistr *= QubitOperator([(i*2, 'Z'),(i*2+1, 'Z')])
            if char==NA:
                qop_paulistr=None
                break
        if qop_paulistr is not None:
            operation_paulistr_table[op] = qop_paulistr

    return operation_paulistr_table

    
"""
A1	0	A1	0
A2	1	A2	1
E1x	2	B1	2
E1y	3	B2	3
E2x	10	A1	0
E2y	11	A2	1
E3x	12	B1	2
E3y	13	B2	3
E4x	20	A1	0
E4y	21	A2	1
E5x	22	B1	2
E5y	23	B2	3

A1g	0	Ag	0
A2g	1	B1g	1
A1u	5	B1u	5
A2u	4	Au	4
E1gx	2	B2g	2
E1gy	3	B3g	3
E1uy	6	B2u	6
E1ux	7	B3u	7
E2gx	10	Ag	0
E2gy	11	B1g	1
E2ux	15	B1u	5
E2uy	14	Au	4
E3gx	12	B2g	2
E3gy	13	B3g	3
E3uy	16	B2u	6
E3ux	17	B3u	7
E4gx	20	Ag	0
E4gy	21	B1g	1
E4ux	25	B1u	5
E4uy	24	Au	4
E5gx	22	B2g	2
E5gy	23	B3g	3
E5uy	26	B2u	6
E5ux	27	B3u	7

C2h	XOR	ID	C2v	XOR	ID	D2	XOR	ID
Ag	00	0	A1	00	0	A1	00	0
Bg	01	1	A2	01	1	B1	01	1
Au	10	2	B1	10	2	B2	10	2
Bu	11	3	B2	11	3	B3	11	3

Cs	XOR	ID	Ci	XOR	ID	C2	XOR	ID
A’	0	0	Ag	0	0	A	0	0
B”	1	1	Au	1	1	B	1	1

D2h	XOR	ID
A1g	000	0
B1g	001	1
B2g	010	2
B3g	011	3
A1u	100	4
B1u	101	5
B2u	110	6
B3u	111	7
"""

import numpy as np
from . import util
class SymmRemover:
    def __init__(self, n_qubit, symmetry_qop_list):
        # qops in symmetry_qop_list are assumed as elements of a stabilizer group.

        self.n_qubit = n_qubit
        
        self.argsym_qop_list = symmetry_qop_list
        self.n_argsym = len(self.argsym_qop_list)
        self.argsym_binrep_matrix = np.array([SymmRemover.binary_string_rep(n_qubit, util.paulistr(qop)) for qop in self.argsym_qop_list])
        
        matrix = self.argsym_binrep_matrix.copy()
        self.targetpauli_qop_list = []

        # gensym[i] = prod_i argsym[i]^{expansion_matrix[j,i]}
        # binrep(gensym[j]) = sum_i expansion_matrix[j,i] * binrep(argsym[i])
        self.expansion_matrix = np.eye(self.n_argsym, dtype=int) 

        # 行基本変形を繰り返して独立なsymをみつける
        r = 0
        for i in range(n_qubit*2):
            for j in range(r, self.n_argsym):
                if matrix[j,i]==0: continue

                qop_targetpauli = QubitOperator((i, 'Z')) if i<n_qubit else QubitOperator((i-n_qubit, 'X'))
                self.targetpauli_qop_list.append(qop_targetpauli)

                row_j = matrix[j].copy()
                row_j_emat = self.expansion_matrix[j].copy()
                
                for k in range(self.n_argsym):
                    if k!=j and matrix[k,i]&1: # k行目のi列目が1ならば
                        matrix[k]^=row_j       # k行目からj行目を引く
                        self.expansion_matrix[k]^=row_j_emat # expansion_matrix にも同じ操作をする
                matrix[j] = matrix[r] # r行目とスワップ
                matrix[r] = row_j
                self.expansion_matrix[j] = self.expansion_matrix[r]
                self.expansion_matrix[r] = row_j_emat
                r += 1
                break

        self.rank = r
        self.gensym_qop_list = []
        for j in range(self.rank):
            qop_gensym = QubitOperator(())
            for i in range(self.n_argsym):
                if self.expansion_matrix[j,i]==1: qop_gensym *= self.argsym_qop_list[i]
            self.gensym_qop_list.append(qop_gensym)
        
        self.clifford_operation_list = []
        rsqrt2 = 0.5**0.5
        for i in range(self.rank):
            qop_gensym = self.gensym_qop_list[i]
            qop_targetpauli = self.targetpauli_qop_list[i]
            qop_clifford_operation = rsqrt2*(qop_gensym + qop_targetpauli)
            self.clifford_operation_list.append(qop_clifford_operation)

            
    def set_eigvals(self, eigvals_argsym):
        assert len(eigvals_argsym) == self.n_argsym
        self.eigvals_argsym = eigvals_argsym

        # calculate eigvals_gensym
        self.eigvals_gensym = [1]*self.rank
        for j in range(self.rank):
            for i in range(self.n_argsym):
                if self.expansion_matrix[j,i]==1: self.eigvals_gensym[j] *= self.eigvals_argsym[i]

                
    def remove_symm_qubits(self, qop):
        qop = self.operate_clifford(qop)
        qop = self.replace_pauli_with_eigenvalue(qop)
        qop = self.taper_qubits_off(qop)
        return qop
        
                
    def operate_clifford(self, qop):
        all_clifford = QubitOperator(())
        for c in self.clifford_operation_list:
            all_clifford = util.cleanup(all_clifford * c)
        return util.cleanup(all_clifford * qop * all_clifford)

    
    def replace_pauli_with_eigenvalue(self, qop):
        targetpauli_eigval_table = {util.paulistr(self.targetpauli_qop_list[j])[0] : self.eigvals_gensym[j]  for j in range(self.rank)}
                
        qop_replaced = QubitOperator()
        for pauli_string, coeff in qop.terms.items():
            pauli_string_replaced = []
            for pauli in pauli_string:
                if pauli in targetpauli_eigval_table:
                    coeff *= targetpauli_eigval_table[pauli]
                else:
                    pauli_string_replaced.append(pauli)
            qop_replaced += QubitOperator(pauli_string_replaced, coeff)

        return qop_replaced


    def taper_qubits_off(self, qop):
        indices_targetpauli = {util.paulistr(qop_tgtpauli)[0][0] for qop_tgtpauli in self.targetpauli_qop_list}
        index_pack_map = []
        index_new = 0
        for index_old in range(self.n_qubit):
            if index_old in indices_targetpauli:
                index_pack_map.append(None)
            else:
                index_pack_map.append(index_new)
                index_new+=1

        qop_tapered = QubitOperator()
        for pauli_string, coeff in qop.terms.items():
            pauli_string_new = []
            for index_old, axis in pauli_string:
                index_new = index_pack_map[index_old]
                pauli_string_new.append((index_new, axis))
            qop_tapered += QubitOperator(pauli_string_new, coeff)

        return qop_tapered

                

    def __str__(self):
        def f(qop):
            ps = util.paulistr(qop)
            sign = '+' if qop.terms[ps]>0 else '-'
            s = str(qop)
            return sign + s[s.find('['):]
        
        s = ''
        s += '='*40 + '\n'
        s += 'Symmetries(Argument [nargsym={:2d}]):\n'.format(self.n_argsym)
        for i in range(self.n_argsym):
            eigval = '{:+2d}'.format(self.eigvals_argsym[i]) if hasattr(self, 'eigvals_argsym') else 'NA'
            s += 'argsym({:2d}):  {}  {}  {}\n'.format(i, eigval, f(self.argsym_qop_list[i]), str(self.argsym_binrep_matrix[i]))
        s += '\n'
            
        s += 'Symmetries(Internal [rank   ={:2d}]):\n'.format(self.rank)
        for i in range(self.rank):
            eigval = '{:+2d}'.format(self.eigvals_gensym[i]) if hasattr(self, 'eigvals_gensym') else 'NA'
            s += 'gensym({:2d}):  {}  {}  '.format(i, eigval, f(self.gensym_qop_list[i]))
            for i, d in enumerate(self.expansion_matrix[i]):                
                if d==1: s+='{} '.format(f(self.argsym_qop_list[i]))
            s+='\n'

        s += '='*40 + '\n'
        return s
        
        

    @staticmethod
    def binary_string_rep(n_qubit, pauli_string):
        a = np.zeros(n_qubit*2, int)
        for i, pauli in pauli_string:
            if pauli in 'XY' : a[i        ] = 1
            if pauli in 'ZY' : a[i+n_qubit] = 1
        return a

    @staticmethod
    def binary_string_rep_inv(binary_string):
        n_qubit = len(binary_string)//2
        qop = QubitOperator(())
        for i, k in enumerate(binary_string):
            if k==0: continue
            if i < n_qubit:
                qop *= QubitOperator((i        , 'X'))
            else:
                qop *= QubitOperator((i-n_qubit, 'Z'))
        return qop
