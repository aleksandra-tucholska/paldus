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
from paldus_commutators import *
from paldus_density import *
from joblib import Parallel, delayed
import sys
import time
import pickle



def check_conditions_og_b(ik, nt, n, ml, mr):

    swr = n + nt + mr
    if (swr == ml):
        app = True
    else:
        app = False

    return app, swr


# omega = (e^-T X e^T + e^-T W Og e^T)
#          --og_a---      ---og_b----

def generate_omega_operators(theory, og_ord):

    og_a = generate_og_a_comb(theory, og_ord)
    print('')
    og_b = generate_og_b_comb(theory, og_ord)

    # print("A")
    # for x in og_a:
    #     print(x)
    # print("B")
    # for x in og_b:
    #     print(x)
    # sys.exit(0)

    og_a_comm = commutators_og_a(og_a, og_ord)
    og_b_comm = commutators_og_b(og_b, og_ord)

    # for x in og_a_comm:
    #     for y in x:
    #         print(y)
    # for x in og_b_comm:
    #     for y in x:
    #         print(y)
    # sys.exit(0)

    print('integrate og a')
    rsimp_a = integrate_og(og_a_comm)
    print('integrate og b')
    rsimp_b = integrate_og(og_b_comm)

    print('wynik_a')
    for x in rsimp_a:
        nv = 0
        no = 0
        for i in x.summation:
            if i in virtual:
                nv += 1
            elif i in occupied:
                no += 1
        nv += og_ord
        no += og_ord
        print(x, nv, no)
    print('wynik_b')
    for x in rsimp_b:
        nv = 0
        no = 0
        for i in x.summation:
            if i in virtual:
                nv += 1
            elif i in occupied:
                no += 1
        nv += og_ord
        no += og_ord
        print(x, nv, no)

    int_a, noint_a = generate_best_intermediates(rsimp_a, 'wmom', 0, True)
    print('')
    print('')
    print('noint a')
    for x in noint_a:
        print(x)
    sys.exit(0)
    print('int a')
    for x in int_a:
        print(x)

    int_b, noint_b = generate_best_intermediates(rsimp_b, 'wmom', 0, True)
    print('noint b')
    for x in noint_b:
        print(x)
    print('int b')
    for x in int_b:
        print(x)
    sys.exit(0)

def generate_og_a_comb(theory, og_ord):
    # omega = (e^-T X e^T + e^-T W e^T)


    if (theory == 'ccsd'):
        maxexc = 2
        maxmu = 2
    elif (theory == 'cc3'):
        maxexc = 3
        maxmu = 3

    og_a = []

    for i in range(0, 4):
        list_t = comb(i, maxexc + 1)
        lent  = len(list_t)
#        print('ln', lent, list_t)

        for j in range(0, lent):
            nt = particle(list_t[j])

            for n in range(-1, 2):
                app, sw = check_conditions_Wlr(i, nt, n, og_ord)

                if (sw >= 1 and sw <= og_ord):
                    if app == True:
                        minidict = {}
                        minidict['T_list'] = list_t[j]
                        minidict['sw'] = sw
                        minidict['og_ord'] = og_ord
                        og_a.append(minidict)
                        
                        print('mini1', minidict)


    return og_a


def generate_og_b_comb(theory, og_ord):

    # omega = (e^-T X e^T + e^-T W e^T)                                                                                                                  
    if (theory == 'ccsd'):
        maxexc = 2
        maxmu = 2
    elif (theory == 'cc3'):
        maxexc = 3
        maxmu = 3

    og_b = []

    for i in range(0, 4):
        list_t = comb(i, maxexc + 1)
        lent  = len(list_t)
#        print('ln', lent, list_t)


        for j in range(0, lent):
            if (1 in list_t[j]):
                continue
            nt = particle(list_t[j])
            for n in range(-2, 3):
                for m in range(1,4):
                    app, sw = check_conditions_og_b(i, nt, n, og_ord, m)
                    if (app == True):
                        minidict = {}
                        minidict['T_list'] = list_t[j]
                        minidict['Omega'] = m
                        minidict['sw'] = sw
                        minidict['og_ord'] = og_ord
                        og_b.append(minidict)
                        print('mini2', minidict)

    return og_b

def commutators_og_b(og_b, og_ord):
    
    og_b_comm = []

    if (og_ord) == 1:
            mu = deepcopy(mu1)
    elif (og_ord) == 2:
            mu =  deepcopy(mu2)
    elif (og_ord) == 3:
            mu =  deepcopy(mu3)
    mu.transpose()

    for x in range(0, len(og_b)):

        ogt = og_b[x]['T_list']
        lt = len(ogt)
        if og_b[x]['Omega'] == 1:
            om = omega1
        elif og_b[x]['Omega'] == 2:
            om = omega2
        elif og_b[x]['Omega'] == 3:
            om = omega3


        was = arithmetic_string(flukt_potentialt)
        
        for w in was:
            if lt !=0:            
                for i in range(0, lt):
                    w1e = evaluate(w, stoo(ogt[i], "t", False))
            else:
                w1e = w
            for j in range(0, len(w1e)):
                mucp = deepcopy(mu)
                disambiguate(mucp, w1e[j])
                w1e[j] = w1e[j].fromleft(mucp)
                disambiguate(om, w1e[j])
                w1e[j] = w1e[j].fromright(om)
                             
            og_b_comm.append(w1e)


    for x in og_b_comm:
        for y in x:
            print('----------------', y)

    return og_b_comm



def commutators_og_a(og_a, og_ord):
    
    og_a_comm = []
    mu_list = []

    if (og_ord) == 1:
            mu =  deepcopy(mu1)
    elif (og_ord) == 2:
            mu =  deepcopy(mu2)
    elif (og_ord) == 3:
            mu =  deepcopy(mu3)
    mu.transpose()


    for x in range(0, len(og_a)):
        print('la', x, og_a[x])
        
        ogt = og_a[x]['T_list']
        lt = len(ogt)
        
        w = arithmetic_string(obs)
        if lt !=0:
            for i in range(0, lt):
                w = evaluate(w, stoo(ogt[i], "t", False))
        for j in range(0, len(w)):
            print('przed', w[j])
            mucp = deepcopy(mu)
            disambiguate(mucp, w[j])
            w[j] = w[j].fromleft(mucp)
            print(x, j, 'po', 'w[j]', w[j])


        og_a_comm.append(w)

    for x in og_a_comm:
        for y in x:
            print('----------------', y)

    return og_a_comm

        

def integrate_og(og_comm):
    
    og_int_list = []
    for x in og_comm:
        for y in x:
            og_int_list.append(y)
            print(y)

    # Z = []
    # for i in range(0, len(og_int_list)):
    #     print('calk', og_int_list[i])
    #     rint = og_int_list[i].integrate()
    #     rint.exec_delta()
    #     if (len(rint) > 0):
    #         print(rint)
    # sys.exit(0)
        

    start = time.time()
    z = Parallel(n_jobs=30,verbose=100)(delayed(integrate)(og_int_list[i]) for i in range(0, len(og_int_list)))
    end = time.time()
    time_new = abs(end-start)
    print('czas caly', time_new)
    print('len', len(z))

    rsimp = arithmetic_string()
    for i in range(0, len(z)):
        og_rint = deepcopy(z[i])
        og_rint.exec_delta()
        og_rsimp = simplify(og_rint)
        if len(og_rsimp) > 0 :
            for x in og_rsimp:
                rsimp.append(x)
        

    print('wynik')
    print(len(rsimp))
    rsimp = simplify(rsimp)
    print(len(rsimp))
    for x in rsimp:
        print(x)

    return rsimp


