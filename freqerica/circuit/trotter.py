import qulacs
from qulacs.gate import to_matrix_gate, RZ
from util import break_operators_into_subsets_dummy
from math import pi
from global_qulacs import rotation_factor
import numpy as np

from time import time
class ElpTime:
    init_tot = 0.0
    init_0 = 0.0
    init_1 = 0.0
    init_2 = 0.0
    init_3 = 0.0
    init_4 = 0.0
    symbol = 0.0
    coef_list_0 = 0.0
    coef_list_1 = 0.0
    coef_list_2 = 0.0
    coef_list_3 = 0.0
    circuit = 0.0

def pauli_product_exponentiate_circuit(circuit, qoperator, control_bit=None):
    assert len(qoperator.terms) == 1, 'len(qoperator.terms) == {}'.format(len(qoperator.terms))

    pauli_product = list(qoperator.terms)[0]
    # e.x. pauli_product = ((0,'X'), (1, 'Z'))
    if len(pauli_product)==0: return

    coeff = qoperator.terms[pauli_product]
#kura>    assert coeff.real == 0, 'coeff.real == {}'.format(coeff.real)
#    print('pauli_product',pauli_product,'coeff',coeff)
    assert abs(coeff.real) < 1e-10, 'coeff.real == {0:10.3e}'.format(coeff.real)
    theta = qoperator.terms[pauli_product].imag

    # 1.
    relevant_qubits = []

    gates = []
    for index_qubit, pauli_axis in pauli_product:
        relevant_qubits += [index_qubit]
        if pauli_axis=='X':
            circuit.add_H_gate(index_qubit)
        if pauli_axis=='Y':
            circuit.add_RX_gate(index_qubit, -pi/4 * rotation_factor)

    # 2.
    pairs_cnot = [(relevant_qubits[i], relevant_qubits[i+1]) for i in range(len(relevant_qubits)-1)]

    for pair_cnot in pairs_cnot:
        circuit.add_CNOT_gate(*pair_cnot)

    rz_gate = RZ(relevant_qubits[-1], theta * rotation_factor)
    if control_bit==None:
        circuit.add_gate(rz_gate)
    else:
        rz_mat_gate = to_matrix_gate(rz_gate)
        rz_mat_gate.add_control_qubit(control_bit, 1)
        circuit.add_gate(rz_mat_gate)

    for pair_cnot in reversed(pairs_cnot):
        circuit.add_CNOT_gate(*pair_cnot)

    # 3.
    gates = []
    for index_qubit, pauli_axis in pauli_product:
        if pauli_axis=='X':
            circuit.add_H_gate(index_qubit)
        if pauli_axis=='Y':
            circuit.add_RX_gate(index_qubit, pi/4 * rotation_factor)


def trotter_step_1st_order(circuit, qubit_operator, control_bit=None, function_break_operator=break_operators_into_subsets_dummy, *, reverse_order=False):
    subsets_qubit_operator = function_break_operator(qubit_operator)

    for subset in subsets_qubit_operator:
        for term in subset:
            pauli_product_exponentiate_circuit(circuit, term, control_bit)

def trotter_step_2nd_order(circuit, qubit_operator, control_bit=None, function_break_operator=break_operators_into_subsets_dummy):
    subsets_qubit_operator = function_break_operator(qubit_operator)

    for subset in subsets_qubit_operator[:-1]:
        for term in subset:
            pauli_product_exponentiate_circuit(circuit, .5*term, control_bit)

    for term in subsets_qubit_operator[-1]:
        pauli_product_exponentiate_circuit(circuit, term, control_bit)

    for subset in reversed(subsets_qubit_operator[:-1]):
        for term in subset:
            pauli_product_exponentiate_circuit(circuit, .5*term, control_bit)


from symbol_for_openfermion import WrappedExpr
class TrotterStep:
    """TrotterStepは、qubit operatorの指数関数をTrotter展開によりシミュレートする量子回路を扱うクラスです。
    qubit operatorの係数にはWrappedExprクラスのインスタンスである数式を含むことができます。

    Args:
        n_qubit (int): 量子ビット数
        qubit_operator (openfermion.QubitOperator or list[tuple[pauli_string, coefficient]]): qubit operator
        order (int, optional): Trotter order (1 or 2)

    Examples:
        .. code-block:: python

            n_site = 4
            fop = FermionOperator((), const)
            for p in range(n_site):
                for q in range(n_site):
                    fop += FermionOperator('p^ q', WrappedExpr('t{}{}'.format(p,q)))
            qop = jordan_wigner(fop)

            ts = TrotterStep(n_site, qop)

    Note:
        hoge

    Attributes:
        circuit (qulacs.ParametricQuantumCircuit): パラメータを含む量子回路

    """

    def __init__(self, n_qubit, qubit_operator, order=2, reverse_exp_sequence=False):

        self._n_qubit = n_qubit

        ElpTime.init_tot -= time()
        ElpTime.init_0 -= time()
        import openfermion
        if isinstance(qubit_operator, openfermion.QubitOperator):
            pauli_string_and_amplitude_list = list(qubit_operator.terms.items())
        elif isinstance(qubit_operator, dict):
            pauli_string_and_amplitude_list = list(qubit_operator.items())

        self._n_term = len(pauli_string_and_amplitude_list)
        self._pauli_string_list   = []
        self._amplitude_expr_list = []
        for pauli_string, amplitude_expr in pauli_string_and_amplitude_list[::-1 if reverse_exp_sequence else 1]:
            self._pauli_string_list.append(pauli_string)
            self._amplitude_expr_list.append(amplitude_expr)
        ElpTime.init_0 += time()

        # initialize parametrized 1st- or 2nd-order trotter circuit
        ElpTime.init_1 -= time()
        self._circuit = qulacs.ParametricQuantumCircuit(self._n_qubit)

        assert order in (1,2), 'order({}) not in (1,2)'.format(order)
        if order==1:
            self._iterm_list = list(range(self._n_term))
            self._rotation_coeff_list = [1.0]*self._n_term
        elif order==2:
            self._iterm_list = list(range(self._n_term-1)) + [self._n_term-1] + list(range(self._n_term-2, -1, -1))
            self._rotation_coeff_list = [0.5]*(self._n_term-1) + [1.0] + [0.5]*(self._n_term-1)

        for iterm in self._iterm_list:
            pauli_string = self._pauli_string_list[iterm]
            target = [index for index, axis in pauli_string]
            pauli_ids = [{'X':1, 'Y':2, 'Z':3}[axis] for index, axis in pauli_string]
            self._circuit.add_parametric_multi_Pauli_rotation_gate(target, pauli_ids, 0.0)
        ElpTime.init_1 += time()

        ElpTime.init_2 -= time()
        self._symbols = set()
        for amplitude_expr in self._amplitude_expr_list:
            if isinstance(amplitude_expr, WrappedExpr):
                self._symbols.update(amplitude_expr.expr.free_symbols)

        self._symbol_number_pairs = {}
        self._amplitude_value_list = [0]*self._n_term
        self._angle_list = [0]*self._n_term
        ElpTime.init_2 += time()

        ElpTime.init_3 -= time()
        import sympy #from sympy import expand, Matrix
        _coef_symbols = []
        for iterm, amp_expr in enumerate(self._amplitude_expr_list):
            if isinstance(amp_expr, WrappedExpr):
                #print(sympy.expand(amp_expr.expr))
                symbols_sorted = sorted(list(self._symbols),key=str)
                _coef_symbols.append([complex(sympy.expand(amp_expr.expr).coeff(symbol)) for symbol in symbols_sorted])
        self._coef_symbols = np.array(_coef_symbols)

        self.subs([] if len(self._symbols)==0 else [(symbol, 0.0) for symbol in self._symbols])
        ElpTime.init_3 += time()

        ElpTime.init_4 -= time()
#only if you need expantion coefficient for T(\theta)>        from sympy import solve, Eq
#only if you need expantion coefficient for T(\theta)>        equation_list = []
#only if you need expantion coefficient for T(\theta)>        for iterm, amp_expr in enumerate(self._amplitude_expr_list):
#only if you need expantion coefficient for T(\theta)>            amp_expr = amp_expr.expr if isinstance(amp_expr, WrappedExpr) else amp_expr
#only if you need expantion coefficient for T(\theta)>            p = WrappedExpr('p{}'.format(iterm)).expr
#only if you need expantion coefficient for T(\theta)>            equation_list.append(Eq(amp_expr, p))
#only if you need expantion coefficient for T(\theta)>        symbol_list = list(self._symbols)
#only if you need expantion coefficient for T(\theta)>        self._fop_coeff_expr = solve(equation_list, symbol_list)
        ElpTime.init_4 += time()
        ElpTime.init_tot += time()

#PRINT        print("TrotterStep:{0:5.2f}".format(ElpTime.init_tot),end="")
#PRINT        print(" ( init_0:{0:5.2f}".format(ElpTime.init_0),end="")
#PRINT        print(" | init_1:{0:5.2f}".format(ElpTime.init_1),end="")
#PRINT        print(" | init_2:{0:5.2f}".format(ElpTime.init_2),end="")
#PRINT        print(" | init_3:{0:5.2f}".format(ElpTime.init_3),end="")
#PRINT        print(" | init_4:{0:5.2f}".format(ElpTime.init_4),end="")
#PRINT        print(")")

    def subs(self, symbol_number_pairs):
        """引数として渡されたテーブルに基づいて量子回路のパラメータを更新します。
        """

        from sympy import Array
        from sympy.utilities.iterables import flatten

        ElpTime.symbol -= time()
        # update self._symbol_number_pairs
        if isinstance(symbol_number_pairs, dict): symbol_number_pairs = list(symbol_number_pairs.items())

        symbol_list, number_list = [], []
        for old, new in symbol_number_pairs:
            if isinstance(old, Array):
                assert old.shape == new.shape
                symbol_list += flatten(old)
                number_list += list(new.flatten()) # assume new as numpy.ndarray
            else:
                symbol_list += [old]
                number_list += [new]
#PRINT        print('slist', symbol_list)
#PRINT        print('nlist', number_list)

        update_pairs = {}
        for symbol, number in zip(symbol_list, number_list):
            symbol = symbol.expr if isinstance(symbol, WrappedExpr) else symbol
            if symbol not in self._symbols: continue
            assert (symbol not in update_pairs) or (update_pairs[symbol]==number) ,'{} {} {}'.format(symbol, update_pairs[symbol], number)  # assert that no conflicts occur
            update_pairs[symbol] = number
#PRINT        print('update_pairs', update_pairs)
        self._symbol_number_pairs.update(update_pairs)
        ElpTime.symbol += time()

        # update self._substituted_coeff_list & self._angle_list
#SLOW>
#SLOW>        for iterm in range(self._n_term):
#SLOW>            ElpTime.coef_list_0 -= time()
#SLOW>            amplitude_expr = self._amplitude_expr_list[iterm]
#SLOW>            ElpTime.coef_list_0 += time()
#SLOW>            ElpTime.coef_list_1 -= time()
#SLOW>            if isinstance(amplitude_expr, WrappedExpr):
#SLOW>                amp_value = complex(amplitude_expr.subs(self._symbol_number_pairs))
#SLOW>            else:
#SLOW>                amp_value = complex(amplitude_expr)
#SLOW>            self._amplitude_value_list[iterm] = amp_value
#SLOW>            ElpTime.coef_list_1 += time()
#SLOW>            ElpTime.coef_list_2 -= time()
#SLOW>            self._angle_list[iterm] = float(amp_value.imag)
#SLOW>            ElpTime.coef_list_2 += time()
#BGN:FAST
        if isinstance(self._amplitude_expr_list[0], WrappedExpr):
            ElpTime.coef_list_3 -= time()
            symbols_sorted = sorted(list(self._symbols),key=str)
            subs_number_vec = np.array([self._symbol_number_pairs.get(symbol) for symbol in symbols_sorted])
            _amplitude_value_vec = np.einsum('ij,j->i', self._coef_symbols, subs_number_vec)
            np.allclose(np.array(self._amplitude_value_list),_amplitude_value_vec)
            self._amplitude_value_list = _amplitude_value_vec.tolist()
            self._angle_list = _amplitude_value_vec.imag.tolist()
            ElpTime.coef_list_3 += time()
        else:
            for iterm in range(self._n_term):
                amplitude_value = complex(self._amplitude_expr_list[iterm])
                self._angle_list[iterm] = float(amplitude_value.imag)
#END:FAST
        ElpTime.circuit -= time()
        # update params of self._circuit
        for iparam, iterm in enumerate(self._iterm_list):
            angle = self._angle_list[iterm]
            rotation_coeff = self._rotation_coeff_list[iterm]
            param = rotation_coeff * angle * rotation_factor
            self._circuit.set_parameter(iparam, param)
        ElpTime.circuit += time()

#PRINT        print("symbol:{0:5.2f}".format(ElpTime.symbol),end="")
#PRINT        print(" | coef_list_0:{0:5.2f}".format(ElpTime.coef_list_0),end="")
#PRINT        print(" | coef_list_1:{0:5.2f}".format(ElpTime.coef_list_1),end="")
#PRINT        print(" | coef_list_2:{0:5.2f}".format(ElpTime.coef_list_2),end="")
#PRINT        print(" | coef_list_3:{0:5.2f}".format(ElpTime.coef_list_3),end="")
#PRINT        print(" | circuit:{0:5.2f}".format(ElpTime.circuit),end="")
#PRINT        print("")

    def __repr__(self):
        return str(self)

    def __str__(self):
        len_pauli_string_print = sum([len(str(isite))+2 for isite in range(self._n_qubit)])

        s = '='*(len_pauli_string_print+101) + '\n'
        s += str(self._circuit)
        s += 'n_qubits     ({:3d})\n'.format(self._n_qubit)
        s += 'free_symbols ({:3d}) : {}\n'.format(len(self._symbols), str(self._symbols))
        for symbol, number in sorted(self._symbol_number_pairs.items(), key=lambda x:str(x[0])):
            s += '{:6} : {:.8e}\n'.format(str(symbol), number)

        s += '-'*(len_pauli_string_print+101) + '\n'

        s += ' '*7+'PauliString' + ' '*(len_pauli_string_print-8) +  'Amplitude[Expr]'
        s += ' '*37 + 'Amplitude[Value]' +'\n'# + ' '*13 + 'Angle(deg)\n'
        for iterm in range(self._n_term):
            pauli_string    = self._pauli_string_list[iterm]
            amplitude_expr  = self._amplitude_expr_list[iterm]
            amplitude_value = self._amplitude_value_list[iterm]
            angle           = self._angle_list[iterm]

            s += '{:3d}:   '.format(iterm)

            ipauli = 0
            for index in range(self._n_qubit):
                if ipauli < len(pauli_string) and index == pauli_string[ipauli][0]:
                    index, axis = pauli_string[ipauli]
                    s += axis + str(index) + ' '
                    ipauli += 1
                else:
                    s += ' ' * (len(str(index))+2)

            str_amp_expr = str(amplitude_expr)
            if len(str_amp_expr)>50: str_amp_expr = str_amp_expr[:47] + '...'
            s += '   {:50}  '.format(str_amp_expr)
            s += '{:+.10e} '.format(amplitude_value.imag)
            s += ' '*11 if amplitude_value.real==0 else '({:+1.0e} ~ real part)   '.format(amplitude_value.real)
            #s += '{:+.7f}\n'.format(angle)
            s += '\n'

        s += '='*(len_pauli_string_print+101)
        return s


def _test():
    from symbol_for_openfermion import WrappedExpr as Symbol
    from openfermion import FermionOperator, jordan_wigner
    x = Symbol('x')
    y = Symbol('y')
    fop  = FermionOperator('3^ 0 3', x)
    fop += FermionOperator('3^ 1 2', y**2)
    fop *= FermionOperator('2^', x + y*4 -2)
    c1 = TrotterStep(n_qubit=4, qubit_operator=jordan_wigner(fop))
    print(c1)
    c1.subs([(x, 1.), (y, 5.)])
    print(c1)

    x, y = 1., 5.
    fop  = FermionOperator('3^ 0 3', x)
    fop += FermionOperator('3^ 1 2', y**2)
    fop *= FermionOperator('2^', x + y*4 -2)
    c1 = TrotterStep(n_qubit=4, qubit_operator=jordan_wigner(fop))
    print(c1)


def _test2():
    from symbol_for_openfermion import WrappedExpr as Symbol
    from openfermion import FermionOperator, jordan_wigner, get_sparse_operator
    from sympy import Array
    np.random.seed(100)
    import scipy

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
    c1 = TrotterStep(n_qubit, (-1j)*qop_symbol)
    print(c1)
    symbol_number_pairs = [(const, const_value), (T1, T1_value)]
    c1.subs(symbol_number_pairs)
    print(c1)

    # c2 : TrotterStep with numerical qop
    fop_number = op1e(const_value, T1_value)
    qop_number = jordan_wigner(fop_number)
    c2 = TrotterStep(n_qubit, (-1j)*qop_number)
    print(c2)

    # c3 : oldtype with numerical qop
    c3 = qulacs.QuantumCircuit(n_qubit)
    trotter_step_2nd_order(c3, (-1j)*qop_number)



    s0 = qulacs.QuantumState(n_qubit)
    s0.set_Haar_random_state()
    s1 = s0.copy()
    s2 = s0.copy()
    s3 = s0.copy()
    from util_qulacs import convert_state_vector
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

    import matplotlib.pyplot as plt
    plt.plot(np.array(corr1).real, label='1')
    plt.plot(np.array(corr2).real, label='2')
    plt.plot(np.array(corr3).real, label='3')
    plt.plot(np.array(corrv).real, label='4')
    plt.legend()
    plt.show()

    # print('s1', s1.get_vector())
    # print('s2', s2.get_vector())
    # print('s3', s3.get_vector())


def _test3():
    n_qubit = 1
    c1 = qulacs.ParametricQuantumCircuit(n_qubit)
    pauli_string = ((0,'Z'),)
    target = [index for index, axis in pauli_string]
    pauli_ids = [{'X':1, 'Y':2, 'Z':3}[axis] for index, axis in pauli_string]
    angle = 1.0
    c1.add_parametric_multi_Pauli_rotation_gate(target, pauli_ids, angle)

    s1 = qulacs.QuantumState(n_qubit)
    c1.update_quantum_state(s1)
    print(s1)


def _example():
    from symbol_for_openfermion import WrappedExpr as Symbol
    from openfermion import FermionOperator, jordan_wigner, get_sparse_operator
    import sympy
    np.random.seed(100)
    import scipy

    n_orb = 2

    const = Symbol('const')
    T1 = [[None for _ in range(n_orb)] for _ in range(n_orb)]
    for p in range(n_orb):
        for q in range(p, n_orb):
            t = Symbol('t{}{}'.format(p,q))
            T1[p][q] = T1[q][p] = t
    T1 = sympy.Array(T1)
    print('T1:')
    print(T1)

    const_value = np.random.rand()
    T1_value = np.random.rand(n_orb,n_orb)*0.01
    T1_value += T1_value.T
    print('const_value = ', const_value)
    print('T1_value = ')
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
    c1 = TrotterStep(n_qubit, (-1j)*qop_symbol)
    print(c1)
    symbol_number_pairs = [(const, const_value), (T1, T1_value)]
    c1.subs(symbol_number_pairs)
    print(c1)

    # c2 : TrotterStep with numerical qop
    fop_number = op1e(const_value, T1_value)
    qop_number = jordan_wigner(fop_number)
    c2 = TrotterStep(n_qubit, (-1j)*qop_number)
    print(c2)

    s0 = qulacs.QuantumState(n_qubit)
    s0.set_Haar_random_state()
    s1 = s0.copy()
    s2 = s0.copy()

    corr1 = []
    corr2 = []
    for t in range(100):
        corr1.append( qulacs.state.inner_product(s0, s1) )
        corr2.append( qulacs.state.inner_product(s0, s2) )
        c1._circuit.update_quantum_state(s1)
        c2._circuit.update_quantum_state(s2)

    import matplotlib.pyplot as plt
    plt.plot(np.array(corr1).real, '--', label='c1')
    plt.plot(np.array(corr2).real, '-.', label='c2')
    plt.legend()
    plt.show()

if __name__=='__main__':
    _example()
