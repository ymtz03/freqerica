import os.path
import shutil

class Outputter:
    def __init__(self):
        with open(os.path.dirname(__file__)+'/template/data.js') as f:
            self.template = f.read()
            
    def save(self):
        shutil.copyfile(os.path.dirname(__file__)+'/template/result.html', 'result.html')
        shutil.copyfile(os.path.dirname(__file__)+'/template/style.css', 'style.css')
        shutil.copyfile(os.path.dirname(__file__)+'/template/chart.js', 'chart.js')

        def shortstr(x):
            return '{: .5e}'.format(x)
    
        with open('data.js','w') as f:
            data_save = self.template[:]
            data_save = data_save.replace('/*dt*/'           , shortstr(self.dt))
            data_save = data_save.replace('/*energy*/'       , ','.join(map(shortstr, self.energy           )))
            data_save = data_save.replace('/*energy_range*/' , ','.join(map(shortstr, self.energy_range     )))
            data_save = data_save.replace('/*spectrum*/'     , ','.join(map(shortstr, self.spectrum.real    )))
            data_save = data_save.replace('/*prony_p_exact*/', ','.join(map(shortstr, self.phase_exact      )))
            data_save = data_save.replace('/*prony_A_exact*/', ','.join(map(shortstr, self.Avec_exact       )))
            data_save = data_save.replace('/*prony_p_trott*/', ','.join(map(shortstr, self.phase_trotter    )))
            data_save = data_save.replace('/*prony_A_trott*/', ','.join(map(shortstr, self.Avec_trotter     )))
            data_save = data_save.replace('/*corr_re_exact*/', ','.join(map(shortstr, self.corr_exact.real  )))
            data_save = data_save.replace('/*corr_im_exact*/', ','.join(map(shortstr, self.corr_exact.imag  )))
            data_save = data_save.replace('/*corr_re_trott*/', ','.join(map(shortstr, self.corr_trotter.real)))
            data_save = data_save.replace('/*corr_im_trott*/', ','.join(map(shortstr, self.corr_trotter.imag)))
            f.write(data_save)
            

def draw(energy, phase_exact, Avec_exact, phase_trotter, Avec_trotter, dt, energy_range, spectrum):    
    import matplotlib.pyplot as plt
    import numpy as np
    plt.figure(figsize=(15,6))
    
    plt.vlines(energy, 0, .5, lw=1)
    
    plt.plot(-(phase_exact+2*np.pi)/dt, Avec_exact, 'b+', ms=12)
    plt.plot(-(phase_exact+4*np.pi)/dt, Avec_exact, 'b+', ms=12, label='prony (exact g(k))')
    plt.plot(-(phase_exact+6*np.pi)/dt, Avec_exact, 'b+', ms=12)

    plt.plot(-(phase_trotter+2*np.pi)/dt, Avec_trotter, 'g1', ms=12)
    plt.plot(-(phase_trotter+4*np.pi)/dt, Avec_trotter, 'g1', ms=12, label='prony (trotter g(k))')
    plt.plot(-(phase_trotter+6*np.pi)/dt, Avec_trotter, 'g1', ms=12)
    
    plt.plot(energy_range, spectrum.real/max(spectrum.real)/2, label='FT')

    plt.xlim(energy_range[0], energy_range[-1])
    plt.ylim(-0.02, 0.52)
    plt.xticks(energy, rotation=70)
    plt.yticks(np.arange(0, 0.55, 0.05))
    plt.xlabel(r'$energy\ (a.u.)$')
    plt.ylabel(r'$A_j$')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    #plt.savefig(jobname+'_figure_energy.png')
    #plt.show()
