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
