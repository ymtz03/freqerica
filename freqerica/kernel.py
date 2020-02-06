import numpy as np
from time import time
from collections import namedtuple

import qulacs
import openfermion

from .hamiltonian import exact
from .circuit import trotter
from .util.qulacsnize import convert_state_vector
from .circuit.universal import prepare_civec_circuit, prepstate
from .op import orbital, util
from .output import freqgraph

# EstimateCorrelationResult = namedtuple('EstimateCorrelationResult', ['corr_exact', 'corr_trotter',
#                                                                      'statevec_exact', 'statevec_trotter',
#                                                                      'time_energy_fidelity'])

EstimateCorrelationResult = namedtuple('EstimateCorrelationResult',
                                       ['corr', 'state_vec', 'energy'])

def estimate_correlation(ham_qop, state, dt, max_trotter_step, outputter, savefilename=None):
    const_term = ham_qop.terms[()]
    n_site = openfermion.count_qubits(ham_qop)
    outputter.n_qubit = n_site
    outputter.n_trott_step = max_trotter_step

    #circuit_time_evo = qulacs.QuantumCircuit(n_site)
    #trotter.trotter_step_2nd_order(circuit_time_evo, -1j * dt * ham_qop)
    trotterstep = trotter.TrotterStep(n_site, -1j * dt * ham_qop)
    circuit_time_evo = trotterstep._circuit
    outputter.ngate = trotterstep.count_gates()
    print(trotterstep)
    print(circuit_time_evo)

    ham_tensor = openfermion.get_sparse_operator(ham_qop)
    
    state_vec_exact = convert_state_vector(n_site, state.get_vector())

    steps = np.arange(max_trotter_step+1)

    save_state_vec_exact   = np.empty((len(steps), 2**n_site), np.complex128)
    save_state_vec_trotter = np.empty((len(steps), 2**n_site), np.complex128)
    save_time_energy_fidelity = np.empty((len(steps), 4), float)
    save_time_energy_fidelity[:, 0] = steps * dt

    time_sim = 0
    
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
        time_bgn = time()
        circuit_time_evo.update_quantum_state(state)
        time_sim += time()-time_bgn
        state_vec_exact = scipy.sparse.linalg.expm_multiply(-1j * dt * ham_tensor, state_vec_exact)

    corr_exact   = np.dot(save_state_vec_exact[0]  , save_state_vec_exact.T)
    corr_trotter = np.dot(save_state_vec_trotter[0], save_state_vec_trotter.T)# * np.exp(-1j * dt * steps * const_term)

    outputter.time_sim = time_sim
    
    if savefilename:
        np.save(savefilename+'_ham_tensor', ham_tensor.todense())
        np.save(savefilename+'_corr_exact.npy', corr_exact)
        np.save(savefilename+'_corr_trotter.npy', corr_trotter)
        np.save(savefilename+'_state_vec_exact.npy', save_state_vec_exact)
        np.save(savefilename+'_state_vec_trotter.npy', save_state_vec_trotter)
        np.save(savefilename+'_time_energy_fidelity.npy', save_time_energy_fidelity)

    return EstimateCorrelationResult(corr_exact, corr_trotter, save_state_vec_exact, save_state_vec_trotter,
                                     save_time_energy_fidelity)

def estimate_corr_new__(ham, dt, max_trotter_step, outputter, savefilename=None):    
    n_site = openfermion.count_qubits(ham.ham_qop)
    #state = ham.state
    #state_vec_exact = convert_state_vector(n_site, state.get_vector())

    steps = np.arange(max_trotter_step+1)

    save_state_vec_exact   = np.empty((len(steps), 2**n_site), np.complex128)
    save_state_vec_trotter = np.empty((len(steps), 2**n_site), np.complex128)
    save_time_energy_fidelity = np.empty((len(steps), 4), float)
    save_time_energy_fidelity[:, 0] = steps * dt

    time_sim = 0

    ham.init_propagate(dt)
    ham_tensor = openfermion.get_sparse_operator(ham.ham_qop)    
    
    for istep, n_trotter_step in enumerate(steps):
        state_vec_trotter = convert_state_vector(n_site, ham.state.get_vector())

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
        time_bgn = time()
        ham.propagate()
        time_sim += time()-time_bgn
        state_vec_exact = scipy.sparse.linalg.expm_multiply(-1j * dt * ham_tensor, state_vec_exact)

    corr_exact   = np.dot(save_state_vec_exact[0]  , save_state_vec_exact.T)
    corr_trotter = np.dot(save_state_vec_trotter[0], save_state_vec_trotter.T)

    outputter.n_qubit = n_site
    outputter.n_trott_step = max_trotter_step
    outputter.ngate = ham.trotterstep.count_gates()
    outputter.time_sim = time_sim
    
    if savefilename:
        np.save(savefilename+'_ham_tensor', ham_tensor.todense())
        np.save(savefilename+'_corr_exact.npy', corr_exact)
        np.save(savefilename+'_corr_trotter.npy', corr_trotter)
        np.save(savefilename+'_state_vec_exact.npy', save_state_vec_exact)
        np.save(savefilename+'_state_vec_trotter.npy', save_state_vec_trotter)
        np.save(savefilename+'_time_energy_fidelity.npy', save_time_energy_fidelity)

    return EstimateCorrelationResult(corr_exact, corr_trotter, save_state_vec_exact, save_state_vec_trotter,
                                     save_time_energy_fidelity)

def estimate_corr_new(ham, dt, max_trotter_step, out, savefilename=None):    
    steps = np.arange(max_trotter_step+1)

    save_state_vec = np.empty((len(steps), 2**ham.n_qubit), np.complex128)
    save_time = steps*dt
    save_energy = np.empty(len(steps), float)

    time_sim = 0

    ham.init_propagate(dt)
    
    for istep, n_trotter_step in enumerate(steps):
        state_vec = ham.get_state_vec()

        # calculate energy
        energy = ham.get_energy()

        # save
        save_state_vec[istep] = state_vec
        save_energy[istep] = energy

        # time propergation
        time_bgn = time()
        ham.propagate()
        time_sim += time()-time_bgn

    corr = np.dot(save_state_vec[0]  , save_state_vec.T)

    out.n_qubit = ham.n_qubit
    out.n_trott_step = max_trotter_step
    out.ngate = ham.ngate
    out.time_sim = time_sim
    
    if savefilename:
        np.save(savefilename+'_ham_tensor', ham_tensor.todense())
        np.save(savefilename+'_corr_exact.npy', corr_exact)
        np.save(savefilename+'_corr_trotter.npy', corr_trotter)
        np.save(savefilename+'_state_vec_exact.npy', save_state_vec_exact)
        np.save(savefilename+'_state_vec_trotter.npy', save_state_vec_trotter)
        np.save(savefilename+'_time_energy_fidelity.npy', save_time_energy_fidelity)

    return EstimateCorrelationResult(corr, save_state_vec, save_energy)


def kernel(mol, norb=None, nelec=None,
           dt=1.0, max_trotter_step=100,
           jobname='noname',
           civec=None, # dict or 'HF', 'S', 'SD', 'SDT'
           mult=None,
           max_excitation=None,
           symm_eigval_map={},
           use_symm_remover=False,
           ft_energy_range=np.arange(-.5, 1.5, 0.002),
):

    out = freqgraph.Outputter()

    print('kernel invoked')

    ### 0 : Prepare molecular orbitals and integrals
    orbprop, energy_casci = orbital.calculate_moint_and_energy(orbital.MoleInput(mol, norb, nelec))
    for iorb in range(len(orbprop.irreps)):
        orbtype = 'core' if iorb < orbprop.ncore else 'active' if iorb < orbprop.ncore+orbprop.ncas else 'virtual'
        print('{:03d}  {:4s}  {:7s}  {:+.6e}'.format(iorb, orbprop.irreps[iorb], orbtype, orbprop.mo_energy[iorb]))
    
    ### 1 : Construct Hamiltonian Operator
    ham = ham_exact = exact.ExactHam(orbprop)
    if use_symm_remover:
        ham = exact.ExactHam(orbprop, symm_eigval_map=symm_eigval_map)
    #TODO-> ham = hamiltonian.Rankone(orbprop)
    #TODO-> ham = hamiltonian.Jastrow(orbprop)
    ham_ref = exact.ReferenceHam(orbprop)
    
    out.nterm         = ham_exact.nterm
    out.nterm_tapered = ham.nterm
    
    ### 3 : Initialize Quantum State
    if mult is None: mult = 1 if nelec%2==0 else 2
    na = (nelec+mult-1)//2
    nb = nelec - na
    wfs = util.listupCSFs(norb, mult, na, nb, ham.remover if hasattr(ham, 'remover') else None, max_excitation=max_excitation)
    #print(wfs)
    wf_initial = {}
    coeff_csf = len(wfs)**(-0.5)
    for csf in wfs:
        for sd, coeff_sd in csf.items():
            if sd not in wf_initial: wf_initial[sd]=0
            wf_initial[sd] += coeff_sd*coeff_csf

    wf_initial = util.mixCSFs(wfs)
    #print(wf_initial)
    state = prepstate(norb*2, wf_initial)
    print('Target Multiplicity : ', mult)
    print('norm(state) : ', state.get_norm())

    ham_ref.set_state(state)
    ham.set_state(state)
        
    ### 4 : Estimate correlation function
    result_corr_ref = estimate_corr_new(ham_ref, dt, max_trotter_step, out, savefilename=None)
    result_corr = estimate_corr_new(ham, dt, max_trotter_step, out, savefilename=None)

    ### 5 : Calculate eigenenergies
    from .analysis.diagonalize import exact_diagonalize
    list_ene, list_num, list_2S, eigvec = exact_diagonalize(ham.ham_qop, ham.n_site, jobname)
    energy_1let = list_ene[(list_num==nelec)*(list_2S==0)]

    from .analysis import prony_like
    corr_exact_extend = prony_like.calc_g_list(result_corr_ref.corr)
    corr_trott_extend = prony_like.calc_g_list(result_corr.corr)
    phase_exact  , Avec_exact   = prony_like.main(corr_exact_extend)
    phase_trotter, Avec_trotter = prony_like.main(corr_trott_extend)
    energy_and_amp_exact = prony_like.estimate(corr_exact_extend, dt, hint_energy=energy_casci)
    energy_and_amp_trott = prony_like.estimate(corr_trott_extend, dt, hint_energy=energy_casci)

    from .analysis.ft import ft
    energy_range = energy_casci + ft_energy_range
    time_range = dt*np.arange(-max_trotter_step, max_trotter_step+1)
    spectrum_exact = ft(time_range, corr_exact_extend, energy_range)
    spectrum_trott = ft(time_range, corr_trott_extend, energy_range)

    out.orbprop        = orbprop
    out.energy         = energy_1let
    out.phase_exact    = phase_exact
    out.Avec_exact     = Avec_exact
    out.phase_trotter  = phase_trotter
    out.Avec_trotter   = Avec_trotter
    out.prony_exact    = energy_and_amp_exact
    out.prony_trott    = energy_and_amp_trott
    out.dt             = dt
    out.energy_range   = energy_range
    out.spectrum_exact = spectrum_exact
    out.spectrum_trott = spectrum_trott
    out.corr_exact     = result_corr_ref.corr
    out.corr_trotter   = result_corr.corr
    out.save()

    print("kernel finished")
