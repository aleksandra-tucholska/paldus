from params import *
import paldus_classes
from paldus_basic import *
from paldus_classes import ugg
from paldus_commutators import *
from copy import deepcopy
from eomccjac import jacobian_loop
from paldus_classes import arithmetic_string
from paldus_classes import disambiguate
from collections import deque
from paldus_classes import pair_permutations
import math
from itertools import product
from itertools import permutations
from paldus_classes import integrate
from paldus_classes import virtual, occupied, general
from fortran_code import *
from templates import  *
from joblib import Parallel, delayed
import sys
import time
import pickle
from joblib import Parallel, delayed



def slater_rules_driver(op, diff_det):

    if op == 21:
        r1 = ebig(['p', 'r'], ['q', 's'], op)
        print('Two electron operator')
        for x in r1:
            print(x)
    if op == 22:
        r1 = ebig(['p', 'r'], ['q', 's'], op)
    elif op == 3:
        r1 = ebig(['p', 'r', 't'], ['q', 's', 'u'], op)
        print('Three electron operator')
        for x in r1:
            print(x)
    elif op == 4:
        r1 = ebig(['p', 'r', 't', 'p>'], ['q', 's', 'u', 'q>'], op)
        print('Four electron operator')

    elif op == 5:
        r1 = ebig(['p', 'r', 't', 'p>', 'r>'], ['q', 's', 'u', 'q>', 's'], op)               
        print('Five electron operator')
    elif op == 6:
        r1 = ebig(['p', 'r', 't', 'p>', 'r>', 't>'], ['q', 's', 'u', 'q>', 's', 'u>'], op)            
        print('Six electron operator')

    # print('')
    # for x in r1:
    #     print(x)
    # print('')
    
    if diff_det == 1:
        eop = ugg()
        eop.operator_idx.append(['a', 'i'])
        eop.operator_type.append('s')
        r1 = r1.fromright(eop)
    elif diff_det == 2:
        eop = ugg()
        eop.operator_idx.append(['a', 'i'])
        eop.operator_idx.append(['b', 'j'])
        eop.operator_type.append('s')
        eop.operator_type.append('s')
        r1 = r1.fromright(eop)
    elif diff_det == 3:
        eop = ugg()
        eop.operator_idx.append(['a', 'i'])
        eop.operator_idx.append(['b', 'j'])
        eop.operator_idx.append(['c', 'k'])
        eop.operator_type.append('s')
        eop.operator_type.append('s')
        eop.operator_type.append('s')
        r1 = r1.fromright(eop)
    elif diff_det == 4:
        eop = ugg()
        eop.operator_idx.append(['a', 'i'])
        eop.operator_idx.append(['b', 'j'])
        eop.operator_idx.append(['c', 'k'])
        eop.operator_idx.append(['d', 'l'])
        eop.operator_type.append('s')
        eop.operator_type.append('s')
        eop.operator_type.append('s')
        eop.operator_type.append('s')
        r1 = r1.fromright(eop)
    elif diff_det == 5:
        eop = ugg()
        eop.operator_idx.append(['a', 'i'])
        eop.operator_idx.append(['b', 'j'])
        eop.operator_idx.append(['c', 'k'])
        eop.operator_idx.append(['d', 'l'])
        eop.operator_idx.append(['a>', 'i>'])
        eop.operator_type.append('s')
        eop.operator_type.append('s')
        eop.operator_type.append('s')
        eop.operator_type.append('s')
        eop.operator_type.append('s')
        r1 = r1.fromright(eop)


    rint = Parallel(n_jobs=30,verbose=50)(delayed(integrate)(r1[i]) for i in range(0, len(r1)))

    rsimp = arithmetic_string()
    for i in range(0, len(rint)):
        rint[i].exec_delta()
        rsimp_small = simplify(rint[i])
        for y in rsimp_small:
            rsimp.append(y)

    rsimp = simplify(rsimp)

    print('')
    print('wynik')
    if len(rsimp) == 0:
        print(0)
    else:
        for x in range(0, len(rsimp)):
            if (x%6 == 0 and x!= 0 and x!= len(rsimp)-1):
                print(rsimp[x], '**')
            else:
                print(rsimp[x])
        




def etwo(ulst, llst):

    res1 = ugg()
    res1.operator_idx.append([ulst[0], llst[0]])
    res1.operator_idx.append([ulst[1], llst[1]])
    res1.operator_type.append('s')
    res1.operator_type.append('s')

    res2 = ugg()
    res2.operator_idx.append([ulst[0], llst[1]])
    res2.delta.append([ulst[1], llst[0]])
    res2.operator_type.append('s')
    res2.num_factor = -1.0


    res = arithmetic_string(res1, res2)

    return res

def ethree(ulst, llst):

    res = arithmetic_string()
    
    res1 = etwo(ulst[0:2], llst[0:2])
    for i in range(0, len(res1)):
        res1[i].operator_idx.append([ulst[2], llst[2]])
        res1[i].operator_type.append('s')
        res.append(res1[i])
    print('res1', res1)
    
    res2 = etwo(ulst[0:2], [llst[2],llst[1]])
    print('res2', res2, len(res2))
    for i in range(0, len(res2)):
        res2[i].delta.append([ulst[2], llst[0]])
        if (i==0):
            res2[i].num_factor = -1.0
        else:
            res2[i].num_factor = 1.0
        res.append(res2[i])
    print('res22', res2)
    res3 = etwo(ulst[0:2], [llst[0],llst[2]])
    for i in range(0, len(res3)):
        res3[i].delta.append([ulst[2], llst[1]])
        if (i==0):
            res3[i].num_factor = -1.0
        else:
            res3[i].num_factor = 1.0
        res.append(res3[i])


    return res


def efour(ulst, llst):

    res = arithmetic_string()
    ost = len(ulst) - 1
    a = ulst[0]
    b = ulst[1]
    c = ulst[2]
    d = ulst[3]
    m = llst[0]
    n = llst[1]
    p = llst[2]
    q = llst[3]


    res1 = ethree(ulst[0:ost], llst[0:ost])
    for i in range(0, len(res1)):
        res1[i].operator_idx.append([ulst[ost], llst[ost]])
        res1[i].operator_type.append('s')
        res.append(res1[i])

    res2 = ethree(ulst[0:ost], [q, n, p])
    for i in range(0, len(res2)):
        res2[i].delta.append([d, m])
        res2[i].num_factor = -1.0 * res2[i].num_factor
        res.append(res2[i])

    res3 = ethree(ulst[0:ost], [m, q, p])
    for i in range(0, len(res3)):
        res3[i].delta.append([d, n])
        res3[i].num_factor = -1.0 * res3[i].num_factor
        res.append(res3[i])

    res4 = ethree(ulst[0:ost], [m, n, q])
    for i in range(0, len(res4)):
        res4[i].delta.append([d, p])
        res4[i].num_factor = -1.0 * res4[i].num_factor
        res.append(res4[i])

    return res

def efive(ulst, llst):

    res = arithmetic_string()
    ost = len(ulst) -1

    a = ulst[0]
    b = ulst[1]
    c = ulst[2]
    d = ulst[3]
    e = ulst[4]
    m = llst[0]
    n = llst[1]
    p = llst[2]
    q = llst[3]
    r = llst[4]
    

    res1 = efour(ulst[0:ost], llst[0:ost])
    for i in range(0, len(res1)):
        res1[i].operator_idx.append([ulst[ost], llst[ost]])
        res1[i].operator_type.append('s')
        res.append(res1[i])

    res2 = efour(ulst[0:ost], [r, n, p, q])
    for i in range(0, len(res2)):
        res2[i].delta.append([e, m])
        res2[i].num_factor = -1.0
        res.append(res2[i])

    res3 = efour(ulst[0:ost], [m, r, p, q])
    for i in range(0, len(res3)):
        res3[i].delta.append([e, n])
        res3[i].num_factor = -1.0
        res.append(res3[i])

    res4 = efour(ulst[0:ost], [m, n, r, q])
    for i in range(0, len(res4)):
        res4[i].delta.append([e, p])
        res4[i].num_factor = -1.0
        res.append(res4[i])

    res5 = efour(ulst[0:ost], [m, n, p, r])
    for i in range(0, len(res5)):
        res5[i].delta.append([e, q])
        res5[i].num_factor = -1.0
        res.append(res5[i])

    return res

def esix(ulst, llst):

    res = arithmetic_string()
    ost = len(ulst) - 1
    f = ulst[ost]
    s = llst[ost]
    
    res1 = efive(ulst[0:ost], llst[0:ost])

    for i in range(0, len(res1)):
        res1[i].operator_idx.append([f, s])
        res1[i].operator_type.append('s')
        res.append(res1[i])

    return res

def ebig(ulst, llst, op):

    n = len(ulst)

    if n == 2:
        r = etwo(ulst, llst)
        for x in r:
            aa = []
            x.summation = ulst + llst
            x.num_factor = x.num_factor / 2.0
            if op == 21:
                x.coefficient = [SLAT2_SYM1]
            elif op == 22:
                x.coefficient = [SLAT2_SYM2]
            for y in range(0, len(ulst)):
                aa.append(ulst[y])
                aa.append(llst[y])
            x.coefficient_idx.append(aa)
    if n == 3:
        r = ethree(ulst, llst)
        for x in r:
            aa = []
            x.summation = ulst + llst
            x.num_factor = x.num_factor/ (2.0*3.0)
            x.coefficient = [SLAT3]
            for y in range(0, len(ulst)):
                aa.append(ulst[y])
                aa.append(llst[y])
            x.coefficient_idx.append(aa)


    if n == 4:
        r = efour(ulst, llst)
        for x in r:
            aa = []
            x.summation = ulst + llst
            x.num_factor = x.num_factor/ (2.0 * 3.0 * 4.0)
            x.coefficient = [SLAT4]
            for y in range(0, len(ulst)):
                aa.append(ulst[y])
                aa.append(llst[y])
            x.coefficient_idx.append(aa)

    if n == 5:
        r = efive(ulst, llst)
        for x in r:
            aa = []
            x.summation = ulst + llst
            x.num_factor = x.num_factor/ (2.0 * 3.0 * 4.0*5.0)
            x.coefficient = [ANY_OP]
            for y in range(0, len(ulst)):
                aa.append(ulst[y])
                aa.append(llst[y])
            x.coefficient_idx.append(aa)


    if n == 6:
        r = esix(ulst, llst)
        for x in r:
            x.summation = ulst + llst
            x.num_factor = x.num_factor/ (2.0 * 3.0 * 4.0 * 5.0 * 6.0)
            x.coefficient = [ANY_OP]
            for y in range(0, len(ulst)):
                aa.append(ulst[y])
                aa.append(llst[y])
            x.coefficient_idx.append(aa)


    return r
