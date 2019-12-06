import numpy as np

def window_func(time):
    return np.exp(-(time/time[-1]*2)**2)

def ft(time, corr, energy_range, use_window=True):
    if use_window:
        return np.sum(np.exp(1j*time*energy_range.reshape(-1,1))*corr*window_func(time), axis=1)
    return np.sum(np.exp(1j*time*energy_range.reshape(-1,1))*corr, axis=1)
