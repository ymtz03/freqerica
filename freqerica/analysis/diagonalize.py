import openfermion
import numpy as np
from collections import namedtuple

Result = namedtuple('Result', ['list_energy', 'list_number', 'list_2S', 'eigvec'])

def _get_op_number(n_site):
    op_number = openfermion.FermionOperator()
    for i in range(n_site):
        op_number += openfermion.FermionOperator(((i,1),(i,0)))
    return op_number

def _get_op_Sz(n_site):
    assert n_site%2==0
    
    op_Sz = openfermion.FermionOperator()
    for i in range(n_site//2):
        op_Sz += openfermion.FermionOperator(((2*i  , 1), (2*i  , 0)), +0.5) # alpha
        op_Sz += openfermion.FermionOperator(((2*i+1, 1), (2*i+1, 0)), -0.5) # beta
    return op_Sz

def _get_op_Sp(n_site):
    assert n_site%2==0

    op_Sp = openfermion.FermionOperator()
    for i in range(n_site//2):
        op_Sp += openfermion.FermionOperator(((2*i  , 1), (2*i+1, 0)))
    return op_Sp

def _get_op_Sm(n_site):
    assert n_site%2==0

    op_Sm = openfermion.FermionOperator()
    for i in range(5):
        op_Sm += openfermion.FermionOperator(((2*i+1, 1), (2*i  , 0)))
    return op_Sm

def _get_op_S2(n_site):
    op_Sp = _get_op_Sp(n_site)
    op_Sm = _get_op_Sm(n_site)
    op_Sz = _get_op_Sz(n_site)
    
    op_S2 = op_Sp * op_Sm + op_Sz * (op_Sz  - openfermion.FermionOperator((),1) )
    return op_S2


def exact_diagonalize(ham_qop, n_site, jobname=None):
    from os.path import isfile
    if jobname and (isfile(jobname+'_list_energy.npy') and
                    isfile(jobname+'_list_number.npy') and
                    isfile(jobname+'_list_2S.npy') and
                    isfile(jobname+'_eigvec.npy')):
        print("Load results of exact diagonalize (jobname='{}')".format(jobname))
        list_energy = np.load(jobname+'_list_energy.npy')
        list_number = np.load(jobname+'_list_number.npy')
        list_2S     = np.load(jobname+'_list_2S.npy')
        eigvec      = np.load(jobname+'_eigvec.npy')
        return Result(list_energy, list_number, list_2S, eigvec)
    
    op_number_tensor = openfermion.get_sparse_operator(_get_op_number(n_site))
    op_S2_tensor = openfermion.get_sparse_operator(_get_op_S2(n_site))

    ham_tensor = openfermion.get_sparse_operator(ham_qop)
    list_energy, eigvec = np.linalg.eigh(ham_tensor.todense())

    nstate = len(list_energy)
    list_number = np.empty(nstate, int)
    list_2S = np.empty(nstate, float)
    
    for istate in range(nstate):
        list_number[istate] = round(openfermion.expectation(op_number_tensor, eigvec[:,istate]).real)
        list_2S[istate] = round((openfermion.expectation(op_S2_tensor, eigvec[:,istate]).real + 0.25)**0.5 *2 -1)

    if jobname:
        print("Save results of exact diagonalize (jobname='{}')".format(jobname))
        np.save(jobname+'_list_energy.npy', list_energy)
        np.save(jobname+'_list_number.npy', list_number)
        np.save(jobname+'_list_2S.npy'    , list_2S)
        np.save(jobname+'_eigvec.npy'     , np.asarray(eigvec))
            
    return Result(list_energy, list_number, list_2S, np.asarray(eigvec))

