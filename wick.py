from itertools import combinations
from itertools import product
from copy import copy
from copy import deepcopy
import sys
#from paldus_classes import ugg
OCC_CRE = 1
OCC_ANI = 2
VIR_CRE = 3
VIR_ANI = 4
ALPHA = 1
BETA = 2
virtual  = ["a", "b", "c", "d", "e", "f", "g", \
                "a>", "b>", "c>", "d>", "e>", "f>", \
                "a>>", "b>>", "c>>", "d>>", "e>>", "f>>"]

occupied = ["i", "j", "k", "l", "m", "n", \
                "i>", "j>", "k>", "l>", "m>", "n>",\
                "i>>", "j>>", "k>>", "l>>", "m>>", "n>>"]

def unique_pairs(n):
    """Compute all partitions of a set of numbers {1, 2, 3, ..., n} into unique pairs.
    Example (n=4):
    [[(1,2), (3,4)],
    [(1,3), (2,4)],
    [(1,4), (2,3)]]
    """
    if n < 2 or n % 2 > 0:
        print("UNIQUE_PAIRS accepts only positive, even integers.")
        sys.exit(1)
    elements = list(range(0, n))
    return list(two_element_subsets(elements))


def two_element_subsets(elements):
    """Compute all the set of all possible partitions of a set (containing an even number of elements)
    into two-element subsets. The resulting subsets have no common elements."""

    if len(elements) > 2:
        a = elements[0]
        for i in range(1, len(elements)):
            b = elements[i]
            ab = (a, b)
            new_elements = elements[1:i] + elements[i+1:]
            for x in two_element_subsets(new_elements):
                yield [ab] + x
    else:
        yield [(elements[0], elements[1])]
        return



def spin_split(e, e_type):
    """ Translate list representing ugg string
    E_ab E_kj T_ck ...
    of type [['a', 'b'], ['k', 'j'], ['c', 'k']]
    from direct definition
    E_ai = (a*_alpha i_alpha + a*_beta i_beta)
    T_ai = (a*_alpha i_alpha - a*_beta i_beta)
    to list of list of dictionaries:
    e_with_spin = [[{a*_alpha b_alpha},{a*_beta b_beta}], 
    [{k*_beta j_beta},{k*_beta j_beta}], 
    [{c*_beta k_beta},{c*_beta k_beta}]]
    
    Return product of e_with_spin
    """

    e_with_spin = []

    for i in range(0, len(e)):
        mini1 = {}
        mini2 = {}
        mini1['idx1'] = e[i][0]
        mini1['idx2'] = e[i][1]
        mini2['idx1'] = e[i][0]
        mini2['idx2'] = e[i][1]
        if e_type[i] == 's':
            mini1['sign'] = 1
            mini2['sign'] = 1
        elif e_type[i] == 't0':
            mini1['sign'] = 1
            mini2['sign'] = -1
        else:
            print('WICK:spin_split: unsupported operator type')
            sys.exit(1)
        if e[i][0] in virtual:
            mini1['idx1_type'] = VIR_CRE
            mini2['idx1_type'] = VIR_CRE
        elif e[i][0] in occupied:
            mini1['idx1_type'] = OCC_CRE
            mini2['idx1_type'] = OCC_CRE
        if e[i][1] in virtual:
            mini1['idx2_type'] = VIR_ANI
            mini2['idx2_type'] = VIR_ANI
        elif e[i][1] in occupied:
            mini1['idx2_type'] = OCC_ANI
            mini2['idx2_type'] = OCC_ANI

        mini1['idx1_spin'] = ALPHA
        mini1['idx2_spin'] = ALPHA
        mini2['idx1_spin'] = BETA
        mini2['idx2_spin'] = BETA
        e_with_spin.append([mini1, mini2])
    c = list(product(*e_with_spin))
    return c

def print_ca_string_single(c):

    if c['sign'] == 1:
        s = "+"
    elif c['sign'] == -1:
        s = "-"
    for j in range(0, len(c['string'])):
        if c['string'][j]['spin'] == ALPHA:
            s += c['string'][j]['idx'].upper()
        elif c['string'][j]['spin'] == BETA:
            s += c['string'][j]['idx'].lower()
        if c['string'][j]['type'] == OCC_CRE or c['string'][j]['type'] == VIR_CRE:
            s += "* "
        else:
            s += " "

    print(s)


def print_ca_string(c):

    for i in range(0, len(c)):
        if c[i]['sign'] == 1:
            s = "+"
        elif c[i]['sign'] == -1:
            s = "-"
        for j in range(0, len(c[i]['string'])):
            if c[i]['string'][j]['spin'] == ALPHA:
                s += c[i]['string'][j]['idx'].upper()
            elif c[i]['string'][j]['spin'] == BETA:
                s += c[i]['string'][j]['idx'].lower()
            if c[i]['string'][j]['type'] == OCC_CRE or c[i]['string'][j]['type'] == VIR_CRE:
                s += "* "
            else:
                s += " "

        print(s)

def generate_ca_string(c):
    
    ca_string = []
    for i in range(0, len(c)):
        ca_mini = {}
        ca_mini['sign'] = 1
        ca_mini['string'] = []
        for j in range(0, len(c[i])):
            ca_mini['sign'] *= c[i][j]['sign']
            minidict = {}
            minidict['idx'] = c[i][j]['idx1']
            minidict['type'] = c[i][j]['idx1_type']
            minidict['spin'] = c[i][j]['idx1_spin']
            ca_mini['string'].append(deepcopy(minidict))
            minidict = {}
            minidict['idx'] = c[i][j]['idx2']
            minidict['type'] = c[i][j]['idx2_type']
            minidict['spin'] = c[i][j]['idx2_spin']
            ca_mini['string'].append(deepcopy(minidict))
        ca_string.append(deepcopy(ca_mini))

#    print_ca_string(ca_string)
    return(ca_string)

def exec_contr(d1, d2):
    "Returns True if contraction is not zero"

    if d1['type'] == OCC_CRE and d2['type'] == OCC_ANI:
        if d1['spin'] != d2['spin']:
            return False
        else:
            return True
    elif d1['type'] == VIR_ANI and d2['type'] == VIR_CRE:
        if d1['spin'] != d2['spin']:
            return False
        else:
            return True
    else:
        return False

def get_sign(sign1, contractions):

    new_order = []
    for i in contractions:
        for j in i:
            new_order.append(j)

    n = len(new_order)
    old_order = []
    for i in range(0, n):
        old_order.append(i)

    n_transp = 0
    for i in range(0, n):
        if old_order[i] != new_order[i]:
            pos = old_order.index(new_order[i])
            if pos == n -1:
                temp_list = old_order[i:pos]
            else:
                temp_list = old_order[i:pos] + old_order[pos+1:n]

            n_transp += pos - i
            old_order = new_order[0:i+1] + temp_list

    if n_transp % 2 > 0:
        sign2 = -1
    else:
        sign2 = 1
    sign = sign1 * sign2

    return sign
        

def remove_null_contractions(ca_string, contractions):

    res = []
    for i in range(0, len(ca_string)):
#        print_ca_string_single(ca_string[i])
        for j in range(0, len(contractions)):
            delta_set = set()
            for k in range(0, len(contractions[j])):
                pos1 = contractions[j][k][0]
                pos2 = contractions[j][k][1]
                is_nonzero  = exec_contr(ca_string[i]['string'][pos1], ca_string[i]['string'][pos2])
                if is_nonzero == False:
                    delta_set = set()
                    break
                idx1 = ca_string[i]['string'][pos1]['idx']
                idx2 = ca_string[i]['string'][pos2]['idx']
                delta_set.add(tuple(sorted((idx1, idx2))))

            if len(delta_set) != 0:
                minires = {}
                minires['delta'] = delta_set
                minires['n_factor'] = get_sign(ca_string[i]['sign'], contractions[j])
                added = False
                for m in range(0, len(res)):
                    if minires['delta'] == res[m]['delta']:
                        res[m]['n_factor'] += minires['n_factor']
                        added = True
                        break
                if added == False:
                    res.append(minires)
#                print(minires)


    result = []
    for x in range(0, len(res)):
        if res[x]['n_factor'] != 0:
            delta_list = []
            for k in res[x]['delta']:
                delta_list.append(list(k))
            minires = {}
            minires['delta'] = delta_list
            minires['n_factor'] = res[x]['n_factor']
            result.append(minires)
            


    # print('result')
    # for x in result:
    #     print(x)

    return result
    
    
                
            
def generate_and_execute_contractions(ca_string):
    
    n = len(ca_string[0]['string'])
    """Generate all possible contractions"""

    contractions = unique_pairs(n)

    """Remove contractions equal to zero"""
    nonzero_contractions = remove_null_contractions(ca_string, contractions)

    return nonzero_contractions

    
def integrate_wick_basic(e, e_type):

    c = spin_split(e, e_type)
    ca_string = generate_ca_string(c)

    nonzero_contractions = generate_and_execute_contractions(ca_string)
    return nonzero_contractions

