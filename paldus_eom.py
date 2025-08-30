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

def task_eom(BRA, KET,theory, trans, pick):

    if BRA == 3:
        task_eom_triple(BRA, KET, 'cc3', 'trans', pick)
        return

    print(pick)

    open_flags = ""
    if pick == "dump":
        open_flags = "wb"        
    elif pick == "load":
        open_flags = "rb"


    if BRA == 1:
        if KET == 1:
            out = open('./pickle/ars11.pkl', open_flags)
        elif KET == 2:
            out = open('./pickle/ars12.pkl', open_flags)
        elif KET == 3:
            out = open('./pickle/ars13.pkl', open_flags)
    elif BRA == 2:
        if KET == 1:
            if theory == 'ccsd':
                out = open('./pickle/ars21.pkl', open_flags)
            elif theory == 'cc3':
#                out = open('./pickle/ars21.pkl', open_flags)
                out = open('./pickle/ars21_cc3.pkl', open_flags)
        elif KET == 2:
            out = open('./pickle/ars22.pkl', open_flags)
        elif KET == 3:
            out = open('./pickle/ars23.pkl', open_flags)

    if pick == "dump":
        if BRA == 1:

            if KET == 1:
                r = evaluate(hamiltoniant, nubj) + evaluate(hamiltoniant, t2, nubj)
            elif KET == 2:
                r = evaluate(hamiltoniant, nubjck)
            elif KET == 3:
                r = evaluate(hamiltonian, nubjckdl)

            rint = r.integrate(bra = ['a', 'i'], braspin = ['s']).scale(0.5)
            for x in rint:
                print(x)
            ars = preprep_for_fortran(rint)
            print('ars')
            for x in ars:
                print(x)

            pickle.dump(ars, out)

        elif BRA == 2:

            if KET == 1:
                if theory == 'ccsd': 
                    r = evaluate(hamiltoniant, nuck) + evaluate(hamiltoniant, t2, nuck) 
                elif theory == 'cc3':
                    r = evaluate(hamiltoniant, nuck) + evaluate(hamiltoniant, t2, nuck) + evaluate(hamiltoniant, t3, nuck)
            elif KET == 2:
                r = evaluate(hamiltoniant, nuckdl) + evaluate(hamiltoniant, t2, nuckdl)
            elif KET == 3:
                r = evaluate(hamiltoniant, nuckdlem)


            rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s', 's'])
            rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], braspin = ['s', 's'])
            rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
    
            rsimp = simplify(rint)
            rint = deepcopy(rsimp)

            r2 = deepcopy(rint)
            for x in r2:
                x.new_delta("a", "b")
                x.new_delta("i", "j")
            rint += r2.scale(-1./2.)

            ars = preprep_for_fortran(rint)
            for x in ars:
                print(x)

            pickle.dump(ars, out)
            
    elif pick == "load":
        ars = pickle.load(out)

    # i = 0
    # for x in ars:
    #     i+=1
    #     print(i, x)

    print('potrojne')
    # Usun wszystkie wyrazy z amplitudami potrojnymi,
    # zostana one doliczone w specjalnym module.
    # if BRA == 2:
    #     if KET == 1:
    #         if theory == 'cc3':
    #             for x in ars:
    #                 for i in range(0, len(x.coefficient)):
    #                     if x.coefficient[i] == CC_AMPLITUDE:
    #                         ln, at = x.amplitude_type(i)
    #                         if ln != '3':
    #                             print(x)
    #                             x.num_factor = 0

    # print('')
    # i =0
    # for x in ars:
    #     i+=1
    #     print(i, x)
    # sys.exit(0)
    # print(len(ars))
    # sys.exit(0)
    ars.cleanup()

    # print(len(ars))
    # sys.exit(0)


    print('wynik')
    kk = 1
    for x in ars:
        print(kk, x)
        kk += 1

    print('')

    delta_list, name1, name_list, arg_list = eom_func(ars, BRA, KET ,theory, trans)

    jacobian_loop(theory, arg_list, 'eom-cc', (BRA,KET), delta_list, name1, name_list)

def task_eom_triple(BRA, KET, theory, trans, pick):

    open_flags = ""
    if pick == "dump":
        open_flags = "wb"        
    elif pick == "load":
        open_flags = "rb"


    if KET == 1:
        o0 = open('./pickle/ars0.pkl', open_flags)
        o6 = open('./pickle/ars6.pkl', open_flags)
        o1 = open('./pickle/ars1.pkl', open_flags)
        o2 = open('./pickle/ars2.pkl', open_flags)
        o3 = open('./pickle/ars3.pkl', open_flags)
        o4 = open('./pickle/ars4.pkl', open_flags)
        o5 = open('./pickle/ars5.pkl', open_flags)
    elif KET == 2:
        o0 = open('./pickle/ars0_32.pkl', open_flags)
        o6 = open('./pickle/ars6_32.pkl', open_flags)
        o1 = open('./pickle/ars1_32.pkl', open_flags)
        o2 = open('./pickle/ars2_32.pkl', open_flags)
        o3 = open('./pickle/ars3_32.pkl', open_flags)
        o4 = open('./pickle/ars4_32.pkl', open_flags)
        o5 = open('./pickle/ars5_32.pkl', open_flags)

    if pick == 'dump':
        
        if KET == 1:
            r = evaluate(hamiltoniant, t2, nudl)
            #r = evaluate(hamiltoniant, t2, polnR1)
        elif KET == 2:
            r = evaluate(hamiltoniant, nudlem)
            
        print('r0 integration')
        r01 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./3.)
        r02 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./6.)
        r0 = r01 + r02
        ars0 = preprep_for_fortran(r0)
        print('ars0')
        for x in ars0:
            print(x)
        pickle.dump(ars0, o0)

        print('r6 integration')
        r61 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./3.)
        r62 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./6.)
        r6 = r61 + r62
        print('pluszak przed preprep')
        for x in r6:
            print(x)
        
        ars6 = preprep_for_fortran(r6)
        print('ars6')
        for x in ars6:
            print(x)
        print('')
        pickle.dump(ars6, o6)
        print('')

        print('r1 integration')
        rint1 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./4.)
        rint2 = r.integrate(bra = ['a', 'k', 'b', 'i', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./12.)
        rint3 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./6.)
        rint4 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./6)
        rint5 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./12.)
        r1 = rint1 + rint2 + rint3 + rint4 + rint5
        ars1 = preprep_for_fortran(r1)
        print('ars1')
        for x in ars6:
            print(x)
            print('')

        pickle.dump(ars1, o1)

        print('r2 integration')
        rint1 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./12.)
        rint2 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./4.)
        rint3 = r.integrate(bra = ['a', 'i', 'b', 'k', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./6.)
        rint4 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./6.)
        rint5 = r.integrate(bra = ['a', 'k', 'b', 'i', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./12.)
        r2 = rint1 + rint2 + rint3 + rint4 + rint5
        ars2 = preprep_for_fortran(r2)
        print('ars2')
        for x in ars6:
            print(x)
            print('')

        pickle.dump(ars2, o2)

        print('r3 integration')
        rint1 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./6.)
        rint2 = r.integrate(bra = ['a', 'i', 'b', 'k', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./6.)
        rint3 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./3.)
        rint4 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./6.)
        rint5 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./6.)
        r3 = rint1 + rint2 + rint3 + rint4 + rint5
        ars3 = preprep_for_fortran(r3)
        print('ars3')
        for x in ars6:
            print(x)
            print('')

        pickle.dump(ars3, o3)

        print('r4 integration')
        rint1 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./6.)
        rint2 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./6.)
        rint3 = r.integrate(bra = ['a', 'k', 'b', 'i', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./6.)
        rint4 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./3.)
        rint5 = r.integrate(bra = ['a', 'i', 'b', 'k', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./6.)
        r4 = rint1 + rint2 + rint3 + rint4 + rint5
        ars4 = preprep_for_fortran(r4)
        print('ars4')
        for x in ars6:
            print(x)
            print('')
            
        pickle.dump(ars4, o4)

        print('r5 integration')
        rint1 = r.integrate(bra = ['a', 'k', 'b', 'i', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./12.)
        rint2 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./12.)
        rint3 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./6.)
        rint4 = r.integrate(bra = ['a', 'i', 'b', 'k', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./6.)
        rint5 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./4.)
        r5 = rint1 + rint2 + rint3 + rint4 + rint5
        ars5 = preprep_for_fortran(r5)
        print('ars5')
        for x in ars6:
            print(x)
            print('')

        pickle.dump(ars5, o5)

    elif pick == 'load':
        ars0 = pickle.load(o0)
        ars6 = pickle.load(o6)        
        ars1 = pickle.load(o1)
        ars2 = pickle.load(o2)
        ars3 = pickle.load(o3)
        ars4 = pickle.load(o4)
        ars5 = pickle.load(o5)

        print('pluszuszusz')
        for x in ars6:
            print(x)

    ars_list = [ars0, ars1, ars2, ars3, ars4, ars5, ars6]

    for x in range(0, len(ars_list)):
        for y in range(0, len(ars_list[x])):
            # print(ars_list[x][y])
             if len(ars_list[x][y].summation) > 0:
                 print('tak')
             

    delta_list, name1, name_list, vanish_case1, arg_list = eom_func_triple(ars_list, BRA, KET, theory, trans)
    print('to tutaj')
    print(len(delta_list))
    print(name1)
    print(len(name_list))
    print(len(vanish_case1))
    print(len(arg_list))

    for x in range(0, len(delta_list)):
        print(delta_list[x], name_list[x], arg_list[x])


    print('sr')
    print('vanish', vanish_case1)
    for x in range(0, len(arg_list)):
        print(x, arg_list[x], name_list[x])

    jacobian_loop(theory, arg_list, 'eom-cc', (BRA,KET), delta_list, name1, name_list, vanish_case1)

   
def integrate_triple_exc_bra(a1, a2, a3, a4, a5, r0):
    
    r = deepcopy(r0)
    
    rint1 = r.integrate(bra = ['a', 'i', 'b', 'k', 'c', 'j'], braspin = ['s', 's', 's']).scale(a1)
    rint2 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k'], braspin = ['s', 's', 's']).scale(a2)
    rint3 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i'], braspin = ['s', 's', 's']).scale(a3)
    rint4 = r.integrate(bra = ['a', 'k', 'b', 'i', 'c', 'j'], braspin = ['s', 's', 's']).scale(a4)
    rint5 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's']).scale(a5)

    rint = rint1 + rint2 + rint3 + rint4 + rint5

    return rint
