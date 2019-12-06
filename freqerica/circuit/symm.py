import qulacs
from ..operator.util import paulistr

class SymmRemoveClifford:
    def __init__(self, n_qubit, remover):
        self.n_qubit = []

        self.circuit_list = []
        for isym in range(remover.rank):
            targetpauli = paulistr(remover.targetpauli_qop_list[isym])
            gensym      = paulistr(remover.gensym_qop_list[isym]     )

            index_tgt, axis_tgt = targetpauli[0]
            indices_gensym = [index for index, axis in gensym]
            pauli_gensym = [axis for index, axis in gensym]

            assert axis_tgt=='X'
            assert index_tgt in indices_gensym
            assert 'X' not in pauli_gensym
            assert 'Y' not in pauli_gensym

            indices_gensym_wo_tgt = [index for index in indices_gensym if index!=index_tgt][::-1]
            ctrl_tgt_pair = [pair for pair in zip(indices_gensym_wo_tgt[:-1], indices_gensym_wo_tgt[1:])]
            
            circuit = qulacs.QuantumCircuit(n_qubit)
            for ctrl, tgt in ctrl_tgt_pair:
                circuit.add_CNOT_gate(ctrl, tgt)
            circuit.add_CNOT_gate(indices_gensym_wo_tgt[-1], index_tgt)
            circuit.add_CZ_gate(indices_gensym_wo_tgt[-1], index_tgt)
            for ctrl, tgt in reversed(ctrl_tgt_pair):
                circuit.add_CNOT_gate(ctrl, tgt)
            circuit.add_H_gate(index_tgt)
                
            self.circuit_list.append(circuit)
