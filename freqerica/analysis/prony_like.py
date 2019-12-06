import numpy as np
from scipy.optimize import lsq_linear

def calc_g_list(corr, K=None):
    if K==None: K=len(corr)-1
    g_list =  np.empty(2*K+1, complex)
    g_list[:K] = np.conjugate(corr[K:0:-1])
    g_list[K:] = corr[:K+1]
    return g_list

def Gmat(l, K, a, g_list):
    assert l>0 and K>0 and a in (0,1)
    
    G = np.empty((l, 2*K-l), complex)
    for i in range(l):
        G[i] = g_list[i+a:2*K-l+i+a]
    return G

def main(g_list, l=None, verbose=True, complexcoeff=False):
    """g_list = [g(-K), g(-K+1), ... g(K-1), g(K)]"""

    K = (len(g_list)-1)//2
    l = l if l!=None else K

    G0 = Gmat(l, K, 0, g_list)
    G1 = Gmat(l, K, 1, g_list)

    # Tmat = np.linalg.pinv(G0.T).dot(G1.T).T
    Tmat = np.empty((l,l), complex)
    for i in range(l):
        result = lsq_linear(G0.T, G1[i])
        Tmat[i] = result.x
    #print('maxabs: ',np.max(abs(Tmat_-Tmat)))
    #print(Tmat)
    #print(Tmat_)
        
    
    eigval, eigvec = np.linalg.eig(Tmat)
    phase = np.angle(eigval)%(2*np.pi)

    Bmat = np.exp(1j*np.arange(l).reshape(-1,1)*phase)
    if complexcoeff:
        #Avec = np.linalg.pinv(Bmat).dot(g_list[K:K+l])
        result = lsq_linear(Bmat[:,:-1]-Bmat[:,-1:], g_list[K:K+l]-Bmat[:,-1], bounds=(0,1))
        Avec = np.concatenate((result.x, np.array([1-np.sum(result.x)])))
        
    else:
        Bmat_extend = np.concatenate((Bmat.real, Bmat.imag), axis=0)
        g_list_extend = np.concatenate((g_list[K:K+l].real, g_list[K:K+l].imag))
        #Avec = np.linalg.pinv(Bmat_extend).dot(g_list_extend)
        result = lsq_linear(Bmat_extend[:,:-1]-Bmat_extend[:,-1:], g_list_extend - Bmat_extend[:,-1], bounds=(0,1))
        Avec = np.concatenate((result.x, np.array([1-np.sum(result.x)])))

    #print('K : {:4d}'.format(K))
    #print('l : {:4d}'.format(l))
    
    return phase, Avec

def estimate(corr, dt, hint_energy, l=None, verbose=True, complexcoeff=False):
    phase, Avec = main(corr, l=l, verbose=verbose, complexcoeff=complexcoeff)

    # e^{-iEt} = e^{i(phase + 2pi*k)} (0<= phase <2pi, k is integer)
    # phase + 2pi*k = -Et
    # E = -(phase + 2pi*k)/t

    #        k=k0+1      k=k0        k=k0-1      k=k0-2
    # ---|-----------|-----------|-----------|-----------|---> E
    #     2pi<------0 2pi<------0 2pi<------0 2pi<------0      phase

    k0 = int(round(-hint_energy*dt/2/np.pi))
    k_range = range(k0-2, k0+2)
    l = len(phase)
    energy_and_amp = np.zeros((2, l*len(k_range)), float)
    for i,k in enumerate(range(k0-2, k0+2)):
        energy_and_amp[0, i*l:(i+1)*l] = -(phase + 2*np.pi*k)/dt
        energy_and_amp[1, i*l:(i+1)*l] = Avec

    return energy_and_amp
    

if __name__=='__main__':
    A = np.array([0.1, 0.2, 0.3, 0.4])
    p = np.array([0.3, 0.4, 0.7, 0.9])*2*np.pi

    K = 5
    k_list = np.arange(-K,K+1)
    g_list = np.sum(A.reshape(-1,1)*np.exp(1j*k_list*p.reshape(-1,1)), axis=0)

    phase, Avec = main(g_list, l=5, complexcoeff=True)
    print('-----------------------------')
    print('    Phase  Phase/2pi  Coeff')
    print('-----------------------------')
    for i, (p, c) in enumerate(zip(phase, Avec)):
        print('{:2d}  {:5.3f}  {:9.7f}  {:5.3f}'.format(i, p, p/(2*np.pi), c))
    print('-----------------------------')
    
    print('Avec  :', Avec, '(sum={})'.format(sum(Avec)))
