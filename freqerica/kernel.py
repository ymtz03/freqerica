import numpy as np
import scipy
from time import time
from collections import namedtuple

import qulacs
import openfermion

from .circuit import trotter
from .util.qulacsnize import convert_state_vector
from .circuit.universal import prepare_civec_circuit, prepstate
from .operator import ham, util
from .output import freqgraph

EstimateCorrelationResult = namedtuple('EstimateCorrelationResult', ['corr_exact', 'corr_trotter',
                                                                     'statevec_exact', 'statevec_trotter',
                                                                     'time_energy_fidelity'])

def estimate_correlation(ham_qop, state, dt, max_trotter_step, savefilename=None):
    const_term = ham_qop.terms[()]
    n_site = openfermion.count_qubits(ham_qop)
    print('n_site : ', n_site)

    circuit_time_evo = qulacs.QuantumCircuit(n_site)
    trotter.trotter_step_2nd_order(circuit_time_evo, -1j * dt * ham_qop)

    ham_tensor = openfermion.get_sparse_operator(ham_qop)
    
    state_vec_exact = convert_state_vector(n_site, state.get_vector())

    steps = np.arange(max_trotter_step+1)

    save_state_vec_exact   = np.empty((len(steps), 2**n_site), np.complex128)
    save_state_vec_trotter = np.empty((len(steps), 2**n_site), np.complex128)
    save_time_energy_fidelity = np.empty((len(steps), 4), float)
    save_time_energy_fidelity[:, 0] = steps * dt

    for istep, n_trotter_step in enumerate(steps):
        state_vec_trotter = convert_state_vector(n_site, state.get_vector())

        # calculate energy and fidelity
        energy_exact   = openfermion.expectation(ham_tensor, state_vec_exact)
        energy_trotter = openfermion.expectation(ham_tensor, state_vec_trotter)
        fidelity = np.abs(np.dot(np.conjugate(state_vec_exact), state_vec_trotter))**2

        # save
        save_state_vec_exact[istep]   = state_vec_exact
        save_state_vec_trotter[istep] = state_vec_trotter
        save_time_energy_fidelity[istep, 1] = energy_exact.real
        save_time_energy_fidelity[istep, 2] = energy_trotter.real
        save_time_energy_fidelity[istep, 3] = fidelity

        # time propergation
        circuit_time_evo.update_quantum_state(state)
        state_vec_exact = scipy.sparse.linalg.expm_multiply(-1j * dt * ham_tensor, state_vec_exact)

    corr_exact   = np.dot(save_state_vec_exact[0]  , save_state_vec_exact.T)
    corr_trotter = np.dot(save_state_vec_trotter[0], save_state_vec_trotter.T) * np.exp(-1j * dt * steps * const_term)

    if savefilename:
        np.save(savefilename+'_ham_tensor', ham_tensor.todense())
        np.save(savefilename+'_corr_exact.npy', corr_exact)
        np.save(savefilename+'_corr_trotter.npy', corr_trotter)
        np.save(savefilename+'_state_vec_exact.npy', save_state_vec_exact)
        np.save(savefilename+'_state_vec_trotter.npy', save_state_vec_trotter)
        np.save(savefilename+'_time_energy_fidelity.npy', save_time_energy_fidelity)

    return EstimateCorrelationResult(corr_exact, corr_trotter, save_state_vec_exact, save_state_vec_trotter,
                                     save_time_energy_fidelity)

def kernel(mol, norb=None, nelec=None,
           dt=1.0, max_trotter_step=100,
           jobname='noname',
           civec=None, # dict or 'HF', 'S', 'SD', 'SDT'
           mult=None,
           max_excitation=None,
           symm_eigval_map={},
           use_symm_remover=False,
           ft_energy_range=np.arange(-.5, 2.5, 0.002),
):

    out = freqgraph.Outputter()

    print('kernel invoked')
    print('mol.atom :', mol.atom)
    print('mol.symmetry :', mol.symmetry)
    print('mol.topgroup :', mol.topgroup)
    print('mol.groupname :', mol.groupname)

    ### 1 : Construct Hamiltonian Operator
    orbprop, energy_casci = ham.calculate_moint_and_energy(ham.MoleInput(mol, norb, nelec))
    ham_fop = ham.construct_exact_ham(orbprop)
    ham_qop = openfermion.jordan_wigner(ham_fop)
    n_site  = openfermion.count_qubits(ham_qop)

    for iorb in range(len(orbprop.irreps)):
        orbtype = 'core' if iorb < orbprop.ncore else 'active' if iorb < orbprop.ncore+orbprop.ncas else 'virtual'
        print('{:03d}  {:4s}  {:7s}  {:+.6e}'.format(iorb, orbprop.irreps[iorb], orbtype, orbprop.mo_energy[iorb]))

    ### 2 : Tapering of qubits using Molecular Symmetry
    from .operator import symm
    #symm_paulistr = symm.symmetry_pauli_string(orbprop, operation_list=['Rz2', 'sxz', 'syz'])
    symm_paulistr = symm.symmetry_pauli_string(orbprop, operation_list=list(symm_eigval_map))
    for symmop, qop in symm_paulistr.items():
        print('[ham_qop, {:5s}] = [ham_qop, {}] = {}'.format(symmop, str(qop), openfermion.commutator(ham_qop, qop)))

    remover = symm.SymmRemover(n_site, list(symm_paulistr.values()) )
    remover.set_eigvals([symm_eigval_map[symm] for symm in symm_paulistr])
    ham_qop_tapered = remover.remove_symm_qubits(ham_qop)
    print(remover)
    print('len(ham         ) = ', len(ham_qop.terms))
    print('len(ham[tapered]) = ', len(ham_qop_tapered.terms))

    nterm         = [0]*(norb*2+1)
    nterm_tapered = [0]*(norb*2+1)
    for k in ham_qop.terms        : nterm        [len(k)] += 1
    for k in ham_qop_tapered.terms: nterm_tapered[len(k)] += 1
    out.nterm         = nterm
    out.nterm_tapered = nterm_tapered
    
    # ### 2 : Prepare Prepare-CIVec-Circuit
    # circuit_state_prep = qulacs.QuantumCircuit(n_site)
    # if civec==None:
    #     civec={}
    #     #-> coeff = (2/(n_site*(n_site-1)))**0.5
    #     #-> for e1 in range(n_site):
    #     #->     for e2 in range(e1):
    #     #->         civec[1<<e1 | 1<<e2] = coeff
    #     n_orb = n_site//2
    #     coeff = ((n_orb)*(n_orb-1)//2 + n_orb)**(-.5) # Comb(norb, 2)**(-1/2)
    #     for e1 in range(n_site//2):
    #         for e2 in range(e1, n_site//2):
    #             if e1==e2:
    #                 civec[1<<(e1*2) | 2<<(e2*2)] = +coeff                    
    #             else:
    #                 civec[1<<(e1*2) | 2<<(e2*2)] = +coeff*(0.5**0.5)
    #                 civec[1<<(e2*2) | 2<<(e1*2)] = -coeff*(0.5**0.5)
    # 
    # det_hf = sum([1<<i for i in range(nelec)])
    # civec = {det_hf:1.0}
    # prepare_civec_circuit(circuit_state_prep, n_site, civec)
    # 
    # 
    # ### 3 : Initialize Quantum State
    # state = qulacs.QuantumState(n_site)
    # circuit_state_prep.update_quantum_state(state)

    ### 3 : Initialize Quantum State
    if mult is None: mult = 1 if nelec%2==0 else 2
    na = (nelec+mult-1)//2
    nb = nelec - na
    wfs = util.listupCSFs(norb, mult, na, nb, remover, max_excitation=max_excitation)
    print(wfs)
    wf_initial = {}
    coeff_csf = len(wfs)**(-0.5)
    for csf in wfs:
        for sd, coeff_sd in csf.items():
            if sd not in wf_initial: wf_initial[sd]=0
            wf_initial[sd] += coeff_sd*coeff_csf

    wf_initial = util.mixCSFs(wfs)
    print(wf_initial)
    state = prepstate(norb*2, wf_initial)
    print('Target Multiplicity : ', mult)
    print('norm(state) : ', state.get_norm())
        

    ### 3' : Operate clifford operation
    from .circuit.symm import SymmRemoveClifford
    symm_remove_circuits = SymmRemoveClifford(norb*2, remover)
    
    from .operator.util import paulistr
    state_symm_removed = state.copy()
    for symm_remove_circuit in symm_remove_circuits.circuit_list:
        symm_remove_circuit.update_quantum_state(state_symm_removed)
    #for sd, coeff in enumerate(state_symm_removed.get_vector()):
    #    if abs(coeff)<1e-10: continue
    #    print('{:010b} : {:+.3e}'.format(sd, coeff))

    statevec = state_symm_removed.get_vector().reshape([2]*(norb*2))    
    n_qubit_tapered = norb*2
    slices = [slice(None)]*(norb*2)
    for qop_tgtpauli in remover.targetpauli_qop_list:
        index_tgtpauli = paulistr(qop_tgtpauli)[0][0]
        slices[-(index_tgtpauli+1)] = 0
        n_qubit_tapered -= 1
        
    slices_debug = [slice(None)]*(norb*2)
    for qop_tgtpauli in remover.targetpauli_qop_list:
        index_tgtpauli = paulistr(qop_tgtpauli)[0][0]
        slices_debug[-(index_tgtpauli+1)] = 1
    print(slices)
    #print(slices_debug)
    #for sd, coeff in enumerate(statevec[slices].reshape(-1)):
    #    if abs(coeff)<1e-10: continue
    #    print('{:010b} : {:+.3e}'.format(sd, coeff))
    
    print(np.dot(statevec[slices].reshape(-1), statevec[slices_debug].reshape(-1)))
    #print(statevec[slices]-statevec[slices_debug])
    #assert np.max(abs(statevec[slices]-statevec[slices_debug])) < 1e-10

    statevec_tapered = statevec[slices].reshape(-1)
    statevec_tapered /= np.linalg.norm(statevec_tapered)
    state_tapered = prepstate(n_qubit_tapered, statevec_tapered)

    ### 4 : Estimate correlation function
    if use_symm_remover:
        result_corr = estimate_correlation(ham_qop_tapered, state_tapered, dt, max_trotter_step, savefilename=None)
    else:
        result_corr = estimate_correlation(ham_qop, state, dt, max_trotter_step, savefilename=None)


    ### 5 : Calculate eigenenergies
    from .analysis.diagonalize import exact_diagonalize
    list_ene, list_num, list_2S, eigvec = exact_diagonalize(ham_qop, n_site, jobname)
    energy_1let = list_ene[(list_num==nelec)*(list_2S==0)]

    from .analysis import prony_like
    corr_exact_extend = prony_like.calc_g_list(result_corr.corr_exact)
    corr_trott_extend = prony_like.calc_g_list(result_corr.corr_trotter)
    phase_exact  , Avec_exact   = prony_like.main(corr_exact_extend)
    phase_trotter, Avec_trotter = prony_like.main(corr_trott_extend)
    energy_and_amp_exact = prony_like.estimate(corr_exact_extend, dt, hint_energy=energy_casci)
    energy_and_amp_trott = prony_like.estimate(corr_trott_extend, dt, hint_energy=energy_casci)

    from .analysis.ft import ft
    #energy_range = np.arange(-8.1, -6.7, 0.0001)
    energy_range = energy_casci + ft_energy_range
    #spectrum = ft(dt*np.arange(-max_trotter_step, max_trotter_step+1), corr_trott_extend, energy_range)
    #spectrum = ft(dt*np.arange(-max_trotter_step, max_trotter_step+1), corr_trott_extend, energy_range, use_window=False)
    spectrum = ft(dt*np.arange(-max_trotter_step, max_trotter_step+1), corr_exact_extend, energy_range)

    #freqgraph.draw(energy_4elec_1let, phase_exact, Avec_exact, phase_trotter, Avec_trotter, dt, energy_range, spectrum)
    out.orbprop       = orbprop
    out.energy        = energy_1let
    out.phase_exact   = phase_exact
    out.Avec_exact    = Avec_exact
    out.phase_trotter = phase_trotter
    out.Avec_trotter  = Avec_trotter
    out.prony_exact   = energy_and_amp_exact
    out.prony_trott   = energy_and_amp_trott
    out.dt            = dt
    out.energy_range  = energy_range
    out.spectrum      = spectrum
    out.corr_exact    = result_corr.corr_exact
    out.corr_trotter  = result_corr.corr_trotter
    out.save()

    print("kernel finished")

    
if __name__=='__main__':
    moleInput = ham.moleInput_example['LiH']
    main(moleInput, dt=2., max_trotter_step=100)
