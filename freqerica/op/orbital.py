import numpy as np
from pyscf import gto, scf, mcscf, ao2mo, symm
from collections import namedtuple

MoleInput = namedtuple('MoleInput', ['mol', 'norb', 'nelec'])
OrbitalProperty = namedtuple('OrbitalProperty', ['hcore', 'hint', 'gint', 'irreps', 'mo_energy', 'mo_occ', 'mo_coeff', 'ncore', 'ncas', 'nelecas'])

def unfold_compact_h2mat(h2mat, norb):
    retval = np.empty([norb]*4, float)
    r=0
    for r1 in range(norb):
        for r2 in range(r1+1):
            c=0
            for c1 in range(norb):
                for c2 in range(c1+1):
                    retval[r1,r2,c1,c2] = h2mat[r,c]
                    retval[r2,r1,c1,c2] = h2mat[r,c]
                    retval[r1,r2,c2,c1] = h2mat[r,c]
                    retval[r2,r1,c2,c1] = h2mat[r,c]
                    c+=1
            r+=1
    return retval


def calculate_moint_and_energy(moleInput):
    mol = moleInput.mol
    
    mf = scf.RHF(moleInput.mol)
    mf.max_cycle = 2000
    hf_energy = mf.kernel()
    print('hf_energy : ', hf_energy)
    
    print('mf.mo_coeff : \n', mf.mo_coeff)
    print('mf.mo_energy : \n', mf.mo_energy)
    print('mf.mo_occ : \n', mf.mo_occ)
    irreps = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mf.mo_coeff) if mol.symmetry else ['A']*moleInput.norb
    print('irreps :', irreps)
    
    norb = moleInput.norb
    nelec = moleInput.nelec
    mc = mcscf.CASCI(mf, norb, nelec)
    h1, hcore = mc.get_h1eff()
    eri = mc.get_h2eff()
    result_casci = mc.kernel()
    energy_casci = result_casci[0]
    print('energy_casci :', energy_casci)
    print(mc.ncore)
    print('hcore', hcore)

    #mo_coeff = mean_field.mo_coeff
    #hint_spacial = mo_coeff.T.dot(mean_field.get_hcore()).dot(mo_coeff)
    #n_orbs = hint_spacial.shape[0]
    #gint_spacial = ao2mo.full(mol, mo_coeff, compact=False).reshape([n_orbs]*4)
    n_site=norb*2

    hint_spatial = h1
    gint_spatial = unfold_compact_h2mat(eri, norb)

    hint = np.zeros([n_site]*2, float)
    hint[0::2, 0::2] = hint_spatial
    hint[1::2, 1::2] = hint_spatial

    gint = np.zeros([n_site]*4, float)
    gint[0::2, 0::2, 0::2, 0::2] = gint_spatial
    gint[1::2, 1::2, 0::2, 0::2] = gint_spatial
    gint[0::2, 0::2, 1::2, 1::2] = gint_spatial
    gint[1::2, 1::2, 1::2, 1::2] = gint_spatial

    return OrbitalProperty(hcore, hint, gint, irreps, mf.mo_energy, mf.mo_occ, mf.mo_coeff, mc.ncore, mc.ncas, mc.nelecas), energy_casci
