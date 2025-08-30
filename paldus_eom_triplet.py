from params import *
import paldus_classes
from paldus_basic import *
from paldus_classes import ugg
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
from templates import *
import sys
import time
import pickle


def triplet_bra_part(r, a, i, b, j, c, k, no_fx=None):
    
    rrint1 = r.integrate(bra = [a, i, b, j, c, k], braspin = ["s", "t0", "s"])
    rrint2 = r.integrate(bra = [a, i, c, k, b, j], braspin = ["s", "t0", "s"])
    rrint  = rrint1.scale(1.0/8.0) + rrint2.scale(1.0/8.0)
    rrsimp = simplify(rrint, no_fx)
#    rrsimp = simplify(rrint)

    return rrsimp

def bra_ort9(r, a, i, b, j, c, k):

    print('rint1')
    rint1 = triplet_bra_part(r, a, i, b, j, c, k)
    print('rint2')
    rint2 = triplet_bra_part(r, a, k, b, j, c, i)
    print('rint3')
    rint3 = triplet_bra_part(r, a, j, b, i, c, k)
    
    print('rint4')
    rint4 = triplet_bra_part(r, b, i, a, j, c, k)
    print('rint5')
    rint5 = triplet_bra_part(r, b, k, a, j, c, i)
    print('rint6')
    rint6 = triplet_bra_part(r, b, j, a, i, c, k)
    print('rint7')
    rint7 = triplet_bra_part(r, c, i, b, j, a, k)
    print('rint8')
    rint8 = triplet_bra_part(r, c, k, b, j, a, i)
    print('rint9')
    rint9 = triplet_bra_part(r, c, j, b, i, a, k)
    
    rint0 = rint1.scale(9.0) + rint2 + rint3 +  \
        rint4 + rint5.scale(-1.0) + rint6.scale(-1.0) + \
        rint7 + rint8.scale(-1.0) + rint9.scale(-1.0) 
    
    rint = rint0.scale(1.0/10.0)
#    rint.clear_fixed()
#    rint.establish_fixed()
    ars9 = preprep_for_fortran(rint, True)
#    print(fixed)
#    ars9 = rint

    return ars9

def bra_ort3_virt(r, a, i, b, j, c, k):

    rint1 = triplet_bra_part(r, a, i, b, j, c, k)
    rint2 = triplet_bra_part(r, a, k, b, j, c, i)
    rint3 = triplet_bra_part(r, a, j, b, i, c, k)
    
    rint = rint1 + rint2 + rint3
    ars3_virt = preprep_for_fortran(rint)
    
    return ars3_virt

def bra_aac(KET, r, a, i, b, j, c, k):

    rint1 = triplet_bra_part(r, a, i, b, j, c, k, True)
    rint2 = triplet_bra_part(r, c, i, b, j, a, k, True)
    rint3 = triplet_bra_part(r, c, k, b, j, a, i, True)
    rint4 = triplet_bra_part(r, c, j, b, i, a, k, True)
    print(fixed)
    rint = rint1.scale(10.0) + rint2 + rint3.scale(-1.0) + rint4.scale(-1.0)
    rint = rint.scale(1.0/10.0)
    for x in range(0, len(rint)):
        rint[x].constraints['a'] = set('a')
        rint[x].constraints['c'] = set('c')
        rint[x].constraints['i'] = set('i')
        rint[x].constraints['j'] = set('j')
        rint[x].constraints['k'] = set('k')
        if KET == 2:
            rint[x].constraints['d'] = set('d')
            rint[x].constraints['e'] = set('e')
            rint[x].constraints['l'] = set('l')
            rint[x].constraints['m'] = set('m')
        elif KET == 1:
            rint[x].constraints['d'] = set('d')
            rint[x].constraints['l'] = set('l')
    print('fxfxfx', fixed)
    ars3_virt = preprep_for_fortran(rint, True)
    print('fixedpo')

    return ars3_virt

def bra_aba(r, a, i, b, j, c, k):

    rint1 = triplet_bra_part(r, a, i, b, j, c, k, True)
    rint2 = triplet_bra_part(r, b, i, a, j, c, k, True)
    rint3 = triplet_bra_part(r, b, k, a, j, c, i, True)
    rint4 = triplet_bra_part(r, b, j, a, i, c, k, True)
    print('fixed', fixed)
    rint = rint1.scale(10.0) + rint2 + rint3.scale(-1.0) + rint4.scale(-1.0)
    rint = rint.scale(1.0/10.0)
    ars3_virt = preprep_for_fortran(rint, True)

    return ars3_virt

def bra_iik(KET, r, a, i, b, j, c, k):

    rint1 = triplet_bra_part(r, a, i, b, j, c, k, True)
    rint2 = triplet_bra_part(r, a, k, b, j, c, i, True)
    rint3 = triplet_bra_part(r, b, k, a, j, c, i, True)
    rint4 = triplet_bra_part(r, c, k, b, j, a, i, True)

    rint = rint1.scale(10.0) + rint2 + rint3.scale(-1.0) + rint4.scale(-1.0)
    rint = rint.scale(1.0/10.0)

    for x in range(0, len(rint)):
        rint[x].constraints['a'] = set('a')
        rint[x].constraints['b'] = set('b')
        rint[x].constraints['c'] = set('c')
        rint[x].constraints['i'] = set('i')
        rint[x].constraints['k'] = set('k')
        if KET == 2:
            rint[x].constraints['d'] = set('d')
            rint[x].constraints['e'] = set('e')
            rint[x].constraints['l'] = set('l')
            rint[x].constraints['m'] = set('m')
        elif KET == 1:
            rint[x].constraints['d'] = set('d')
            rint[x].constraints['l'] = set('l')
    ars3_virt = preprep_for_fortran(rint, True)

    return ars3_virt

def bra_iji(r, a, i, b, j, c, k):

    rint1 = triplet_bra_part(r, a, i, b, j, c, k, True)
    rint2 = triplet_bra_part(r, a, j, b, i, c, k, True)
    rint3 = triplet_bra_part(r, b, j, a, i, c, k, True)
    rint4 = triplet_bra_part(r, c, j, b, i, a, k, True)

    rint = rint1.scale(10.0) + rint2 + rint3.scale(-1.0) + rint4.scale(-1.0)
    rint = rint.scale(1.0/10.0)
    ars3_virt = preprep_for_fortran(rint, True)

    return ars3_virt

def bra_aac_iik(KET, r, a, i, b, j, c, k):

    rint1 = triplet_bra_part(r, a, i, b, j, c, k, True)
    rint2 = triplet_bra_part(r, c, k, b, j, a, i, True)

    rint = rint1.scale(10.0) + rint2.scale(-1.0) 
    rint = rint.scale(1.0/10.0)
    for x in range(0, len(rint)):
        rint[x].constraints['a'] = set('a')
        rint[x].constraints['c'] = set('c')
        rint[x].constraints['i'] = set('i')
        rint[x].constraints['k'] = set('k')
        if KET == 2:
            rint[x].constraints['d'] = set('d')
            rint[x].constraints['e'] = set('e')
            rint[x].constraints['l'] = set('l')
            rint[x].constraints['m'] = set('m')
        elif KET == 1:
            rint[x].constraints['d'] = set('d')
            rint[x].constraints['l'] = set('l')
    ars3_virt = preprep_for_fortran(rint, True)

    return ars3_virt

def bra_aba_iik(KET, r, a, i, b, j, c, k):

    rint1 = triplet_bra_part(r, a, i, b, j, c, k, True)
    rint2 = triplet_bra_part(r, b, k, a, j, c, i, True)

    rint = rint1.scale(10.0) + rint2.scale(-1.0)
    rint = rint.scale(1.0/10.0)
    for x in range(0, len(rint)):
        rint[x].constraints['a'] = set('a')
        rint[x].constraints['b'] = set('b')
        rint[x].constraints['i'] = set('i')
        rint[x].constraints['k'] = set('k')
        if KET == 2:
            rint[x].constraints['d'] = set('d')
            rint[x].constraints['e'] = set('e')
            rint[x].constraints['l'] = set('l')
            rint[x].constraints['m'] = set('m')
        elif KET == 1:
            rint[x].constraints['d'] = set('d')
            rint[x].constraints['l'] = set('l')
    ars3_virt = preprep_for_fortran(rint, True)

    return ars3_virt

def bra_aac_iji(KET, r, a, i, b, j, c, k):

    rint1 = triplet_bra_part(r, a, i, b, j, c, k, True)
    rint2 = triplet_bra_part(r, b, j, a, i, c, k, True)

    rint = rint1.scale(10.0) + rint2.scale(-1.0)
    rint = rint.scale(1.0/10.0)
    for x in range(0, len(rint)):
        rint[x].constraints['a'] = set('a')
        rint[x].constraints['c'] = set('c')
        rint[x].constraints['i'] = set('i')
        rint[x].constraints['j'] = set('j')
        if KET == 2:
            rint[x].constraints['d'] = set('d')
            rint[x].constraints['e'] = set('e')
            rint[x].constraints['l'] = set('l')
            rint[x].constraints['m'] = set('m')
        elif KET == 1:
            rint[x].constraints['d'] = set('d')
            rint[x].constraints['l'] = set('l')
    ars3_virt = preprep_for_fortran(rint, True)

    return ars3_virt


def bra_aba_iji(r, a, i, b, j, c, k):

    rint1 = triplet_bra_part(r, a, i, b, j, c, k, True)
    rint2 = triplet_bra_part(r, c, j, b, i, a, k, True)

    rint = rint1.scale(10.0) + rint2.scale(-1.0)
    rint = rint.scale(1.0/10.0)
    ars3_virt = preprep_for_fortran(rint, True)

    return ars3_virt


def bra_ort3_occ(r, a, i, b, j, c, k):

    rint1 = triplet_bra_part(r, a, i, b, j, c, k)
    rint4 = triplet_bra_part(r, b, i, a, j, c, k)
    rint7 = triplet_bra_part(r, c, i, b, j, a, k)
    
    rint = rint1 + rint4 + rint7

    ars3_occ = preprep_for_fortran(rint)
    
    return ars3_occ



def task_eom_triplet(BRA, KET, theory, trans, pick, lmfold = "", rmfold = ""):


    if BRA == 3:
        task_eom_triple_triplet(BRA, KET, 'cc3', 'trans', pick, \
                                    lmfold = lmfold, rmfold = rmfold)
        return

    open_flags = ""
    if pick == "dump":
        open_flags = "wb"        
    elif pick == "load":
        open_flags = "rb"


    if BRA == 1:
        if KET == 1:
            out = open('./pickle/triplet/ars11.pkl', open_flags)
        elif KET == 2:
            if rmfold == "p":
                out = open('./pickle/triplet/ars12p.pkl', open_flags)
            elif rmfold == "m":
                out = open('./pickle/triplet/ars12m.pkl', open_flags)
        elif KET == 3:
            out = open('./pickle/triplet/ars13.pkl', open_flags)
    elif BRA == 2:
        print('z')
        if KET == 1:
            print('a1')
            if lmfold == "p":
                print( 'pwero')
                out = open('./pickle/triplet/ars2p1.pkl', open_flags)
            elif lmfold == "m":
                out = open('./pickle/triplet/ars2m1.pkl', open_flags)
        elif KET == 2:
            if rmfold == "p":
                if lmfold == "p":
                    out = open('./pickle/triplet/ars2p2p.pkl', open_flags)
                elif lmfold == "m":
                    out = open('./pickle/triplet/ars2m2p.pkl', open_flags)
            elif rmfold == "m":
                if lmfold == "p":
                    out = open('./pickle/triplet/ars2p2m.pkl', open_flags)
                elif lmfold == "m":
                    out = open('./pickle/triplet/ars2m2m.pkl', open_flags)
        elif KET == 3:
            if lmfold == "p":
                out = open('./pickle/triplet/ars2p3.pkl', open_flags)
            elif lmfold == "m":
                out = open('./pickle/triplet/ars2m3.pkl', open_flags)
    # elif BRA == 3:
    #     if KET == 1:
    #         out = open('./pickle/triplet/ars31.pkl', open_flags)
    #     if KET == 2:
    #         if rmfold == "p":
    #             out = open('./pickle/triplet/ars32p.pkl', open_flags)
    #         if rmfold == "m":
    #             out = open('./pickle/triplet/ars32m.pkl', open_flags)

    if pick == "dump":
        if BRA == 1:
            if KET == 1:
                nu = operat1(["b", "j"], "t0")
                r = evaluate(hamiltoniant, nu) + evaluate(hamiltoniant, t2, nu)
                
            elif KET == 2:
                if rmfold == "p":
                    
                    nu1 = operat2([["b", "j"], ["c", "k"]], ["t0", "s"])
                    nu2 = operat2([["c", "k"], ["b", "j"]], ["t0", "s"])
                    r = evaluate(hamiltoniant, nu1) +  evaluate(hamiltoniant, nu2)

                    for x in r:
                        print(x)
                    print('')

                elif rmfold == "m":
                    
                    nu1 = operat2([["b", "j"], ["c", "k"]], ["t0", "s"])
                    nu2 = operat2([["c", "k"], ["b", "j"]], ["t0", "s"])
                    nu2.scale(-1.0)
                    r = evaluate(hamiltoniant, nu1) +  evaluate(hamiltoniant, nu2)                  

            elif KET == 3:
                
                nu1 = operat3([["b", "j"], ["c", "k"], ["d", "l"]], ["s", "t0", "s"])
                nu2 = operat3([["b", "j"], ["d", "l"], ["c", "k"]], ["s", "t0", "s"])

                r = evaluate(hamiltonian, nu1) + evaluate(hamiltonian, nu2)

            print('wynik r')
            for x in r:
                print(x)

            rint = r.integrate(bra = ['a', 'i'], braspin = ["t0"]).scale(0.5)
            for x in rint:
                print(x)


            ars = preprep_for_fortran(rint)
            print('WYNIK')
            for x in ars:
                print(x)
            print('')
            
            pickle.dump(ars, out)
            sys.exit(0)

        elif BRA == 2:

            if KET == 1:
                nu = operat1(["c", "k"], "t0")
                if theory == 'ccsd':
                    r = evaluate(hamiltoniant, nu) + evaluate(hamiltoniant, t2, nu) 
                elif theory == 'cc3':
                    r = evaluate(hamiltoniant, nu) + evaluate(hamiltoniant, t2, nu)\
                        +  evaluate(hamiltoniant, t3, nu)
                    #                    r = evaluate(hamiltoniant, t3, pluszR1)

            elif KET == 2:
                if rmfold == "p":
                    nu1 = operat2([["c", "k"], ["d", "l"]], ["t0", "s"])
                    nu2 = operat2([["d", "l"], ["c", "k"]], ["t0", "s"])
                    r = evaluate(hamiltoniant, nu1) + evaluate(hamiltoniant, t2, nu1) + \
                        evaluate(hamiltoniant, nu2) + evaluate(hamiltoniant, t2, nu2)
                elif rmfold == "m":
                    nu1 = operat2([["c", "k"], ["d", "l"]], ["t0", "s"])
                    nu2 = operat2([["d", "l"], ["c", "k"]], ["t0", "s"])
                    nu2.scale(-1.0)
                    r = evaluate(hamiltoniant, nu1) + evaluate(hamiltoniant, t2, nu1) + \
                        evaluate(hamiltoniant, nu2) + evaluate(hamiltoniant, t2, nu2)

            elif KET == 3:
                 nu1 = operat3([["c", "k"], ["d", "l"], ["e", "m"]], ["s", "t0", "s"])
                 nu2 = operat3([["c", "k"], ["e", "m"], ["d", "l"]], ["s", "t0", "s"])
                 r = evaluate(hamiltoniant, nu1) + evaluate(hamiltoniant, nu2)

                # nu1 = operat2([["c", "k"], ["d", "l"]], ["t0", "s"])
                # nu2 = operat2([["d", "l"], ["c", "k"]], ["t0", "s"])
                # nu2.scale(-1.0)
                # r = arithmetic_string(nu1, nu2) 
                # print(g_1t)
                # print(pluszR3)
                # print('')
                # r = evaluate(ht, pluszR3)
                # print('po evaluate')
                # for x in r:
                #     print(x)
                # print('')

            if lmfold == "p":
                rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ["t0", "s"]).scale(1.0/8.0)
                rint2 = r.integrate(bra = ['b', 'j', 'a', 'i'], braspin = ["t0", "s"]).scale(1.0/8.0)
                rint  = rint1 + rint2
            elif lmfold == "m":
                rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ["t0", "s"]).scale(1.0/8.0)
                rint2 = r.integrate(bra = ['b', 'j', 'a', 'i'], braspin = ["t0", "s"]).scale(-1.0/8.0)
                rint  = rint1 + rint2


            print('wynik przed upraszczaniem')
            for x in rint:
                print(x)
            print('po exec delta')
            for x in rint:
                x.exec_delta()

            rint.cleanup()
            print('po cleanup')

            print('przed symplifikacja')

            print('po uproszczeniu')
            rsimp = simplify(rint)
            print('po simplify')
            print('wynik')
            for x in rsimp:
                print(x)
            rint = deepcopy(rsimp)        

            ars = preprep_for_fortran(rint)

            pickle.dump(ars, out)
            sys.exit(0)

        # elif BRA == 3:

        #     if KET == 1:
        #         nu = operat1(["d", "l"], "t0")
        #         r = evaluate(hamiltoniant, t2, nu)
        #     elif KET == 2:
        #         nu1 = operat2([["d", "l"], ["e", "m"]], ["t0", "s"])
        #         nu2 = operat2([["e", "m"], ["d", "l"]], ["t0", "s"])
        #         if rmfold == "m":
        #             nu2.scale(-1.0)
        #         r = evaluate(hamiltoniant, nu1) + evaluate(hamiltoniant, nu2) 

        #     rint1 = triplet_bra_part(r, 'a', 'i', 'b', 'j', 'c', 'k')
        #     rint2 = triplet_bra_part(r, 'a', 'k', 'b', 'j', 'c', 'i')
        #     rint3 = triplet_bra_part(r, 'a', 'j', 'b', 'i', 'c', 'k')
            
        #     rint4 = triplet_bra_part(r, 'b', 'i', 'a', 'j', 'c', 'k')
        #     rint5 = triplet_bra_part(r, 'b', 'k', 'a', 'j', 'c', 'i')
        #     rint6 = triplet_bra_part(r, 'b', 'j', 'a', 'i', 'c', 'k')
            
        #     rint7 = triplet_bra_part(r, 'c', 'i', 'b', 'j', 'a', 'k')
        #     rint8 = triplet_bra_part(r, 'c', 'k', 'b', 'j', 'a', 'i')
        #     rint9 = triplet_bra_part(r, 'c', 'j', 'b', 'i', 'a', 'k')
            
            
        #     rint0 = rint1.scale(9.0) + rint2 + rint3 +  \
        #         rint4 + rint5.scale(-1.0) + rint6.scale(-1.0) + \
        #         rint7 + rint8.scale(-1.0) + rint9.scale(-1.0) 

        #     rint = rint0.scale(1.0/10.0)

        #     ars = preprep_for_fortran(rint)
        #     pickle.dump(ars, out)

    elif pick == "load":
        ars = pickle.load(out)


    ars.cleanup()

    print('wynikf')
    for x in ars:
        print(x)


    # eom_func_easy(ars, BRA, KET, theory, trans,
    #               triplet = "triplet", lmfold=lmfold, rmfold=rmfold)
#    sys.exit(0)
    delta_list, name1, name_list, arg_list = eom_func(ars, BRA, KET, theory, trans, 
                                                                      triplet = "_triplet", lmfold=lmfold, rmfold=rmfold)

    # print('arg_list')
    # print(name1)
    # for x in name_list:
    #     print(x)

    jacobian_loop(theory, arg_list, 'eom-cc', (BRA,KET), delta_list, name1, name_list,
                  triplet = "_triplet", lmfold=lmfold, rmfold=rmfold)



def task_eom_triple_triplet(BRA, KET, theory, trans, pick, \
                                lmfold = "", rmfold = ""):

    open_flags = ""
    if pick == "dump":
        open_flags = "wb"        
    elif pick == "load":
        open_flags = "rb"

    # nu1 = operat3([["d", "l"], ["e", "m"], ["f", "n"]], ["s", "t0", "s"])
    # nu2 = operat3([["d", "l"], ["f", "n"], ["e", "m"]], ["s", "t0", "s"])

    # r = arithmetic_string(nu1) + arithmetic_string(nu2)
    # print('r')
    # for x in r:
    #     print(x)
    # print('startb9')
    # ars9 = bra_ort9(r, 'a', 'i', 'b', 'j', 'c', 'k')
    # print('ars9')
    # for x in ars9:
    #     print(x)
    # print('')
    # sys.exit(0)


        
    if KET == 1:
        o9 = open('./pickle/triplet/ars9_31.pkl', open_flags)
        o_aac = open('./pickle/triplet/ars_aac_31.pkl', open_flags)
        o_aba =open('./pickle/triplet/ars_aba_31.pkl', open_flags)
        o_iik =open('./pickle/triplet/ars_iik_31.pkl', open_flags)
        o_iji =open('./pickle/triplet/ars_iji_31.pkl', open_flags)
        o_aac_iik =open('./pickle/triplet/ars_aac_iik_31.pkl', open_flags)
        o_aac_iji =open('./pickle/triplet/ars_aac_iji_31.pkl', open_flags)
        o_aba_iik =open('./pickle/triplet/ars_aba_iik_31.pkl', open_flags)
        o_aba_iji =open('./pickle/triplet/ars_aba__iji_31.pkl', open_flags)

        # o3_virt = open('./pickle/triplet/ars3_virt_31.pkl', open_flags)
        # o3_occ = open('./pickle/triplet/ars3_occ_31.pkl', open_flags)
        # o1 = open('./pickle/triplet/ars1_31.pkl', open_flags)
    elif KET == 2:
        if rmfold == "p":
            o9 = open('./pickle/triplet/ars9_32p.pkl', open_flags)
            o_aac =open('./pickle/triplet/ars_aac_32p.pkl', open_flags)
            o_aba =open('./pickle/triplet/ars_aba_32p.pkl', open_flags)
            o_iik =open('./pickle/triplet/ars_iik_32p.pkl', open_flags)
            o_iji =open('./pickle/triplet/ars_iji_32p.pkl', open_flags)
            o_aac_iik =open('./pickle/triplet/ars_aac_iik_32p.pkl', open_flags)
            o_aac_iji =open('./pickle/triplet/ars_aac_iji_32p.pkl', open_flags)
            o_aba_iik =open('./pickle/triplet/ars_aba_iik_32p.pkl', open_flags)
            o_aba_iji =open('./pickle/triplet/ars_aba__iji_32p.pkl', open_flags)


            # o3_virt = open('./pickle/triplet/ars3_virt_32p.pkl', open_flags)
            # o3_occ = open('./pickle/triplet/ars3_occ_32p.pkl', open_flags)
            # o1 = open('./pickle/triplet/ars1_32p.pkl', open_flags)
        if rmfold == "m":
            o9 = open('./pickle/triplet/ars9_32m.pkl', open_flags)
            o_aac =open('./pickle/triplet/ars_aac_32m.pkl', open_flags)
            o_aba =open('./pickle/triplet/ars_aba_32m.pkl', open_flags)
            o_iik =open('./pickle/triplet/ars_iik_32m.pkl', open_flags)
            o_iji =open('./pickle/triplet/ars_iji_32m.pkl', open_flags)
            o_aac_iik =open('./pickle/triplet/ars_aac_iik_32m.pkl', open_flags)
            o_aac_iji =open('./pickle/triplet/ars_aac_iji_32m.pkl', open_flags)
            o_aba_iik =open('./pickle/triplet/ars_aba_iik_32m.pkl', open_flags)
            o_aba_iji =open('./pickle/triplet/ars_aba__iji_32m.pkl', open_flags)

            # o3_virt = open('./pickle/triplet/ars3_virt_32m.pkl', open_flags)
            # o3_occ = open('./pickle/triplet/ars3_occ_32m.pkl', open_flags)
            # o1 = open('./pickle/triplet/ars1_32m.pkl', open_flags)
            
    if pick == 'dump':
        
        
        if KET == 1:
            nu = operat1(["d", "l"], "t0")
            r = evaluate(hamiltoniant, t2, nu)
            #r = evaluate(hamiltoniant, t2, pluszR1)
        elif KET == 2:
            nu1 = operat2([["d", "l"], ["e", "m"]], ["t0", "s"])
            nu2 = operat2([["e", "m"], ["d", "l"]], ["t0", "s"])
            if rmfold == "m":
                nu2.scale(-1.0)
            r = evaluate(hamiltoniant, nu1) + evaluate(hamiltoniant, nu2) 
            #r = evaluate(hamiltoniant, pluszR2m)

        # print('r')
        # for x in r:
        #     print(x)
        # print('')
        # sys.exit(0)
        ars9 = bra_ort9(r, 'a', 'i', 'b', 'j', 'c', 'k')
        print('ars9')
        for x in ars9:
            print(x)
        print('')
        pickle.dump(ars9, o9)
#        sys.exit(0)

        ars_aac = bra_aac(KET, r, 'a', 'i', 'a', 'j', 'c', 'k')
        print('aac')
        for x in ars_aac:
            print(x)
        print('')
        pickle.dump(ars_aac, o_aac)

        ars_aba = bra_aba(r, 'a', 'i', 'b', 'j', 'a', 'k')
        print('ars_aba')
        for x in ars_aba:
            print(x)
        print('')        
        pickle.dump(ars_aba, o_aba)


        ars_iik = bra_iik(KET, r, 'a', 'i', 'b', 'i', 'c', 'k')
        print('ars_iik')
        for x in ars_iik:
            print(x)
        print('')
        pickle.dump(ars_iik, o_iik)

        ars_iji = bra_iji(r, 'a', 'i', 'b', 'j', 'c', 'i')
        print('ars_iji')
        for x in ars_iji:
            print(x)
        print('')
        pickle.dump(ars_iji, o_iji)

        ars_aac_iik = bra_aac_iik(KET, r, 'a', 'i', 'a', 'i', 'c', 'k')
        print('ars_aac_iik')
        for x in ars_aac_iik:
            print(x)
        print('')
        pickle.dump(ars_aac_iik, o_aac_iik)

        ars_aba_iik = bra_aba_iik(KET, r, 'a', 'i', 'b', 'i', 'a', 'k')
        print('ars_aba_iik')
        for x in ars_aba_iik:
            print(x)
        print('')
        pickle.dump(ars_aba_iik, o_aba_iik)

        ars_aac_iji = bra_aac_iji(KET, r, 'a', 'i', 'a', 'j', 'c', 'i')
        print('ars_aac_iji')
        for x in ars_aac_iji:
            print(x)
        print('')
        pickle.dump(ars_aac_iji, o_aac_iji)

        ars_aba_iji = bra_aba_iji(r, 'a', 'i', 'b', 'j', 'a', 'i')
        print('ars_aba_iji')
        for x in ars_aba_iji:
            print(x)
        print('')
        pickle.dump(ars_aba_iji, o_aba_iji)

        sys.exit(0)


#         # r9 integration
#         ars9 = bra_ort9(r, 'a', 'i', 'b', 'j', 'c', 'k')
#         pickle.dump(ars9, o9)
#         # r3_virt integration
#         ars3_virt = bra_ort3_virt(r, 'a', 'i', 'b', 'j', 'c', 'k')
# #        ars3_virt = bra_ort9(r, 'a', 'i', 'b', 'j', 'c', 'k')
#         pickle.dump(ars3_virt, o3_virt)
#         # r3_occ integration
#         ars3_occ = bra_ort3_occ(r, 'a', 'i', 'b', 'j', 'c', 'k')
# #        ars3_occ = bra_ort9(r, 'a', 'i', 'b', 'j', 'c', 'k')
#         pickle.dump(ars3_occ, o3_occ)
#         # r1 integration
#         ars10 = triplet_bra_part(r, 'a', 'i', 'b', 'j', 'c', 'k')
# #        ars1 = bra_ort9(r, 'a', 'i', 'b', 'j', 'c', 'k')
#         ars1 = preprep_for_fortran(ars10)
#         pickle.dump(ars1, o1)
        
    elif pick == 'load':
        # ars9 = pickle.load(o9)
        # ars3_virt = pickle.load(o3_virt)
        # ars3_occ = pickle.load(o3_occ)
        # ars1 = pickle.load(o1)

        #        ars = pickle.load(o9)

        # for x in ars:
        #     print(x)

        ars9 = pickle.load(o9)
        ars_aac = pickle.load(o_aac)
        ars_aba = pickle.load(o_aba)
        ars_iik = pickle.load(o_iik)
        ars_iji = pickle.load(o_iji)
        ars_aac_iik = pickle.load(o_aac_iik)
        ars_aac_iji = pickle.load(o_aac_iji)
        ars_aba_iik = pickle.load(o_aba_iik)
        ars_aba_iji = pickle.load(o_aba_iji)

    # for x in range(0, len(ars9)):
    #     ars9[x].constraints['a'] = set('a')
    #     ars9[x].constraints['b'] = set('b')
    #     ars9[x].constraints['c'] = set('c')
    #     ars9[x].constraints['i'] = set('i')
    #     ars9[x].constraints['j'] = set('j')
    #     ars9[x].constraints['k'] = set('k')
    #     if KET == 2:
    #         ars9[x].constraints['d'] = set('d')
    #         ars9[x].constraints['e'] = set('e')
    #         ars9[x].constraints['l'] = set('l')
    #         ars9[x].constraints['m'] = set('m')
    #     elif KET == 1:
    #         ars9[x].constraints['d'] = set('d')
    #         ars9[x].constraints['l'] = set('l')

    # for x in range(0, len(ars_aba)):
    #     ars_aba[x].constraints['a'] = set('a')
    #     ars_aba[x].constraints['b'] = set('b')
    #     ars_aba[x].constraints['i'] = set('i')
    #     ars_aba[x].constraints['j'] = set('j')
    #     ars_aba[x].constraints['k'] = set('k')
    #     if KET == 2:
    #         ars_aba[x].constraints['d'] = set('d')
    #         ars_aba[x].constraints['e'] = set('e')
    #         ars_aba[x].constraints['l'] = set('l')
    #         ars_aba[x].constraints['m'] = set('m')
    #     elif KET == 1:
    #         ars_aba[x].constraints['d'] = set('d')
    #         ars_aba[x].constraints['l'] = set('l')

    # for x in range(0, len(ars_iji)):
    #     ars_iji[x].constraints['a'] = set('a')
    #     ars_iji[x].constraints['b'] = set('b')
    #     ars_iji[x].constraints['c'] = set('c')
    #     ars_iji[x].constraints['i'] = set('i')
    #     ars_iji[x].constraints['j'] = set('j')
    #     if KET == 2:
    #         ars_iji[x].constraints['d'] = set('d')
    #         ars_iji[x].constraints['e'] = set('e')
    #         ars_iji[x].constraints['l'] = set('l')
    #         ars_iji[x].constraints['m'] = set('m')
    #     elif KET == 1:
    #         ars_iji[x].constraints['d'] = set('d')
    #         ars_iji[x].constraints['l'] = set('l')

    # for x in range(0, len(ars_aba_iji)):
    #     ars_aba_iji[x].constraints['a'] = set('a')
    #     ars_aba_iji[x].constraints['b'] = set('b')
    #     ars_aba_iji[x].constraints['i'] = set('i')
    #     ars_aba_iji[x].constraints['j'] = set('j')
    #     if KET == 2:
    #         ars_aba_iji[x].constraints['d'] = set('d')
    #         ars_aba_iji[x].constraints['e'] = set('e')
    #         ars_aba_iji[x].constraints['l'] = set('l')
    #         ars_aba_iji[x].constraints['m'] = set('m')
    #     elif KET == 1:
    #         ars_aba_iji[x].constraints['d'] = set('d')
    #         ars_aba_iji[x].constraints['l'] = set('l')


    #             #v9    v1         v2      v3     v4       v5            v6          v7          v8
    #             #1     2           3     4       5        6             7            8           9
    ars_list = [ars9, ars_aac, ars_aba, ars_iik, ars_iji, ars_aac_iik, ars_aac_iji, ars_aba_iik, ars_aba_iji]

    # k = 0
    # for x in ars_list:
    #     k += 1
    #     print(k)
    #     z = 0
    #     for y in x:            
    #         print(z, y)
    #         z += 1

    # sys.exit(0)
    #    ars_list = [ars9, ars1, ars3_occ, ars3_virt]
    # eom_func_triple_triplet_easy(ars_list, BRA, KET, theory, trans,
    #               triplet = "triplet", lmfold=lmfold, rmfold=rmfold)

    # eom_func_easy(ars, BRA, KET, theory, trans,                                                                                                                      
    #               triplet = "triplet", lmfold=lmfold, rmfold=rmfold)  
    # sys.exit(0)

 
    delta_list, name1, name_list, vanish_case1, arg_list = eom_func_triple_triplet(ars_list, BRA, KET, theory, trans, 
                                                                           triplet = "_triplet", lmfold=lmfold, rmfold=rmfold)

    jacobian_loop(theory, arg_list, 'eom-cc', (BRA,KET), delta_list, name1, name_list, vanish_case1, 
                  triplet = "_triplet", lmfold=lmfold, rmfold=rmfold)

