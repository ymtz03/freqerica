from qulacs import Observable
import numpy as np
from itertools import product

def convert_qoperator_to_observable(n_qubits, qop):
    observable = Observable(n_qubits)

    for term in qop.terms:
        operators = ["{} {}".format(axis, index_qubit) for index_qubit, axis in term]
        #print(operators)
        observable.add_operator(qop.terms[term].real, ' '.join(operators))

    return observable

def convert_state_vector(n_qubits, state_vector):
    return state_vector.reshape([2]*n_qubits).transpose(tuple(range(n_qubits-1,-1,-1))).reshape(-1)
