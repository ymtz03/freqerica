import numpy as np
import scipy
from time import time
from collections import namedtuple

import qulacs
import openfermion

from .circuit import trotter
from .util.qulacsnize import convert_state_vector
from .circuit.universal import prepare_civec_circuit
from .operator import ham


EstimateCorrelationResult = namedtuple('EstimateCorrelationResult', ['corr_exact', 'corr_trotter',
                                                                     'statevec_exact', 'statevec_trotter',
                                                                     'time_energy_fidelity'])

def estimate_correlation(ham_qop, state, dt, max_trotter_step, savefilename=None):
    const_term = ham_qop.terms[()]
    n_site = openfermion.count_qubits(ham_qop)

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
           dt=1.0, max_trotter_step=200,
           jobname='noname',
           civec=None # dict or 'HF', 'S', 'SD', 'SDT'
):
    print('kernel invoked')
    print(mol.atom)

    ### 1 : Construct Hamiltonian Operator
    ham_fop = ham.construct_exact_ham(ham.MoleInput(mol, norb, nelec))
    ham_qop = openfermion.jordan_wigner(ham_fop)
    n_site  = openfermion.count_qubits(ham_qop)

    ### 2 : Prepare Prepare-CIVec-Circuit
    circuit_state_prep = qulacs.QuantumCircuit(n_site)
    if civec==None:
        civec={}
        #-> coeff = (2/(n_site*(n_site-1)))**0.5
        #-> for e1 in range(n_site):
        #->     for e2 in range(e1):
        #->         civec[1<<e1 | 1<<e2] = coeff
        n_orb = n_site//2
        coeff = ((n_orb)*(n_orb-1)//2 + n_orb)**(-.5) # Comb(norb, 2)**(-1/2)
        for e1 in range(n_site//2):
            for e2 in range(e1, n_site//2):
                if e1==e2:
                    civec[1<<(e1*2) | 2<<(e2*2)] = +coeff                    
                else:
                    civec[1<<(e1*2) | 2<<(e2*2)] = +coeff*(0.5**0.5)
                    civec[1<<(e2*2) | 2<<(e1*2)] = -coeff*(0.5**0.5)
    
    prepare_civec_circuit(circuit_state_prep, n_site, civec)

    ### 3 : Initialize Quantum State
    state = qulacs.QuantumState(n_site)
    circuit_state_prep.update_quantum_state(state)

    result_corr = estimate_correlation(ham_qop, state, dt, max_trotter_step, savefilename=None)

    ###################################################
    state = qulacs.QuantumState(n_site)
    circuit_state_prep.update_quantum_state(state)
    state_vec_exact = convert_state_vector(n_site, state.get_vector())
    
    import matplotlib.pyplot as plt
    plt.figure(figsize=(15,6))
    
    from .analysis.diagonalize import exact_diagonalize
    list_ene, list_num, list_2S, eigvec = exact_diagonalize(ham_qop, n_site, jobname)
    energy_4elec_1let = list_ene[(list_num==2)*(list_2S==0)]
    #energy_4elec_3let = list_ene[(list_num==2)*(list_2S==2)]
    plt.vlines(energy_4elec_1let, 0, .5, lw=1)
    #plt.vlines(energy_4elec_3let, 0, 1, 'r', lw=1, label='3let')
    weight_4elec_1let = np.abs(np.dot(state_vec_exact, eigvec[:, (list_num==2)*(list_2S==0)]))**2
    print(np.sum(weight_4elec_1let))
    #-> plt.plot(energy_4elec_1let, weight_4elec_1let, 'rx', label='exact')
    #-> np.save(jobname + '_energy_4elec_1let.npy', energy_4elec_1let)

    from .analysis import prony_like
    corr_exact_extend = prony_like.calc_g_list(result_corr.corr_exact)
    #phase, Avec = prony_like.main(corr_exact_extend, l=32)
    phase, Avec = prony_like.main(corr_exact_extend)
    plt.plot(-(phase+2*np.pi)/dt, Avec, 'b+', ms=12)
    plt.plot(-(phase+4*np.pi)/dt, Avec, 'b+', ms=12, label='prony (exact g(k))')
    plt.plot(-(phase+6*np.pi)/dt, Avec, 'b+', ms=12)
    l = len(phase)
    m = 5
    prony_save = np.empty((2,l*m), float)
    for i in range(m):
        prony_save[0, i*l:i*l+l] = -(phase+m*np.pi)/dt
        prony_save[1, i*l:i*l+l] = Avec
    np.save(jobname + '_prony_save_exact.npy', prony_save)
    

    corr_trotter_extend = prony_like.calc_g_list(result_corr.corr_trotter)
    phase, Avec = prony_like.main(corr_trotter_extend)
    plt.plot(-(phase+2*np.pi)/dt, Avec, 'g1', ms=12)
    plt.plot(-(phase+4*np.pi)/dt, Avec, 'g1', ms=12, label='prony (trotter g(k))')
    plt.plot(-(phase+6*np.pi)/dt, Avec, 'g1', ms=12)
    prony_save = np.empty((2,l*m), float)
    for i in range(m):
        prony_save[0, i*l:i*l+l] = -(phase+m*np.pi)/dt
        prony_save[1, i*l:i*l+l] = Avec
    np.save(jobname + '_prony_save_trotter.npy', prony_save)

    
    from .analysis.ft import ft
    energy_range = np.arange(-8.1, -6.7, 0.0001)
    spectrum = ft(dt*np.arange(-max_trotter_step, max_trotter_step+1), corr_trotter_extend, energy_range)
    plt.plot(energy_range, spectrum.real/max(spectrum.real)/2, label='FT')

    plt.xlim(energy_range[0], energy_range[-1])
    plt.ylim(-0.02, 0.52)
    plt.xticks(energy_4elec_1let, rotation=70)
    plt.yticks(np.arange(0, 0.55, 0.05))
    plt.xlabel(r'$energy\ (a.u.)$')
    plt.ylabel(r'$A_j$')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(jobname+'_figure_energy.png')
    #plt.show()

    #time_plot = np.arange(-max_trotter_step*dt, (max_trotter_step+1)*dt, dt)
    #plt.plot(time_plot, corr_trotter_extend.real, 'r-' , label='trotter.real')
    #plt.plot(time_plot, corr_trotter_extend.imag, 'r--', label='trotter.imag')
    #plt.plot(time_plot, corr_exact_extend  .real, 'b-' , label='exact.real')
    #plt.plot(time_plot, corr_exact_extend  .imag, 'b--', label='exact.imag')
    #-> time_plot = np.arange(0, (max_trotter_step+1)*dt, dt)
    #-> plt.plot(time_plot, result_corr.corr_trotter.real, 'r-' , lw=1, label='trotter.real')
    #-> plt.plot(time_plot, result_corr.corr_trotter.imag, 'r--', lw=1, label='trotter.imag')
    #-> plt.plot(time_plot, result_corr.corr_exact  .real, 'b-' , lw=1, label='exact.real')
    #-> plt.plot(time_plot, result_corr.corr_exact  .imag, 'b--', lw=1, label='exact.imag')
    #-> plt.xlabel('k')
    #-> plt.ylabel('g(k)')
    #-> plt.grid()
    #-> plt.legend()
    #-> plt.savefig(jobname+'_figure_corr.png')
    #-> plt.show()
    
    
if __name__=='__main__':
    moleInput = ham.moleInput_example['LiH']
    main(moleInput, dt=2., max_trotter_step=100)
