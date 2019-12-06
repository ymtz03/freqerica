import openfermion

def break_operators_into_subsets(qubit_operator):
    """Break a qubit operator into subsets, where the terms in each subset all mutually commute."""

    qop_zero = openfermion.QubitOperator('',0)
    subsets = []
    for term in qubit_operator:
        #print('term ',term)
        new_subset = True
        for subset in subsets:
            #print('subset ',subset)
            is_commutable_with_all_term_in_subset = True
            for term_in_subset in subset:
                #print('term_in_subset ',term_in_subset)
                if openfermion.utils.commutator(term, term_in_subset) != qop_zero:
                    is_commutable_with_all_term_in_subset = False
                    #print(term,'and',term_in_subset,'dont commute.')
                    break
            
            if is_commutable_with_all_term_in_subset:
                subset += term
                new_subset = False
                break
                
        if new_subset:
            subsets.append(term)
            
    return subsets


def break_operators_into_subsets_dummy(qubit_operator):
    """Break a qubit operator into subsets, where the terms in each subset all mutually commute."""

    qop_zero = openfermion.QubitOperator('',0)
    subsets = []
    for term in qubit_operator:
        subsets.append(term)
            
    return subsets


def paulistr(qop):
    assert len(qop.terms)==1
    return list(qop.terms)[0]

def cleanup(qop, thresh=1e-12):
    retval = openfermion.QubitOperator()
    terms = retval.terms
    for k, v in qop.terms.items():
        if abs(v)>thresh: terms[k] = v 
    return retval

def listupCSFs(norb, mult, na, nb, remover=None, max_excitation=None):
    """mult, na, nb, removerが指定する対称性 をもつCSFを全列挙する
    
    max_excitationが指定されている場合、
    HF配置（対称性を全て無視して下の軌道からna+nb個の電子を入れた配置）からの
    励起レベルがmax_excitation以下のCSFだけを残す
    """

    from .symm import SymmRemover
    if remover is None: remover = SymmRemover(norb*2, [])

    if max_excitation is None: max_excitation = norb*2
        
    from itertools import combinations, product
    o = [1<<index for index in range(norb*2)]
    occ_a_list = combinations(o[0::2], na)
    occ_b_list = combinations(o[1::2], nb)
    seed_a_list = {sum(oa) for oa in occ_a_list}
    seed_b_list = {sum(ob) for ob in occ_b_list}
    
    seed_list = []
    for oa, ob in product(seed_a_list, seed_b_list):
        seed = oa+ob

        excitation_level = na+nb
        for i in range(na+nb):
            if (seed>>i)&1: excitation_level-=1
        if excitation_level>max_excitation: continue

        # check if the seed is satisfies the symmetry condition
        have_appropriate_symm = True
        for isym in range(remover.rank):
            gensym      = paulistr(remover.gensym_qop_list[isym])
            
            indices_gensym = [index for index, axis in gensym]
            pauli_gensym = [axis for index, axis in gensym] # for assertion only

            assert 'X' not in pauli_gensym
            assert 'Y' not in pauli_gensym

            parity = 1
            for i in indices_gensym:
                if seed&(1<<i): parity *= -1

            eigval_gensym = remover.eigvals_gensym[isym]
            if eigval_gensym!=parity:
                #print(('!!-> {:0' + str(norb*2) + 'b} is NOT satisfy symm condition').format(seed))
                have_appropriate_symm = False
                break

        if not have_appropriate_symm: continue
        seed_list.append(oa+ob)
        #print(('{:0' + str(norb*2) + 'b}').format(seed))

    #s2_tgt = mult-1 # 1let->0, 2let->1, 3let->2, ...
    mult_max = 1 + min(na+nb, norb*2-(na+nb))
    remove_mult_list = [mult_rmv for mult_rmv in range(mult_max, -1, -2) if mult_rmv!=mult]
    wfs = []
    
    for seed in seed_list:
        wf = {seed: 1.0}
        for mult_rmv in remove_mult_list:
            wf = remove(wf, norb, mult, mult_rmv)

        r = sum([c*c for c in wf.values()])**(-0.5)
        for sd in wf: wf[sd]*=r
        
        wfs.append(wf)

    return wfs

def remove(wf, norb, mult_tgt, mult_rmv):
    # wf := {determinant: coeff}
    # op_projection := (S^2 - Sr(Sr+1))/(St(St+1)-Sr(Sr+1)) = (S-S+ + Sz(Sz+1) - Sr(Sr+1))/(St(St+1)-Sr(Sr+1))

    wf_arg = dict(wf)
    
    # S+ operation
    wf_new = {}
    for sd, coeff in wf.items():
        for orb in range(norb):
            # move beta elec to alpha
            a = 1<<(orb*2  )
            b = 1<<(orb*2+1)
            if (not sd & a) and sd & b:
                sd_new = sd + a - b
                if sd_new not in wf_new: wf_new[sd_new] = 0
                wf_new[sd_new] += coeff
    wf = wf_new

    # S- operation
    wf_new = {}
    for sd, coeff in wf.items():
        for orb in range(norb):
            # move alpha elec to beta
            a = 1<<(orb*2  )
            b = 1<<(orb*2+1)
            if (not sd & b) and sd & a:
                sd_new = sd + b - a
                if sd_new not in wf_new: wf_new[sd_new] = 0
                wf_new[sd_new] += coeff
    wf = wf_new

    ssq_rmv = (mult_rmv*mult_rmv-1)/4
    ssq_tgt = (mult_tgt*mult_tgt-1)/4

    for sd, coeff in wf_arg.items():
        sz = 0
        for orb in range(norb):
            if (sd>>(orb*2  ))&1: sz+=0.5
            if (sd>>(orb*2+1))&1: sz-=0.5
        if sd not in wf: wf[sd] = 0
        wf[sd] += coeff*(sz*(sz+1) - ssq_rmv)

    thresh = 1e-10
    for sd_0 in [sd for sd in wf if abs(wf[sd])<thresh]:
        del wf[sd_0]

    denom = ssq_tgt - ssq_rmv
    for sd in wf:
        wf[sd]/=denom
    
    return wf_new


def mixCSFs(csf_list):
    def inner_product(wf1, wf2):
        retval = 0
        for sd in wf1:
            if sd in wf2: retval+=wf1[sd]*wf2[sd]
        return retval

    wf = {}
    for csf in csf_list:
        if abs(inner_product(wf, csf)) > 0.1: continue
        
        for sd, coeff in csf.items():
            if sd not in wf: wf[sd] = 0
            wf[sd] += coeff

    r = sum([c*c for c in wf.values()])**(-0.5)
    for sd in wf: wf[sd]*=r
    return wf

def printwf(wf):
    return '\n'.join(['{:08b} : {:+.3e}'.format(sd, coeff) for sd, coeff in wf.items()])
