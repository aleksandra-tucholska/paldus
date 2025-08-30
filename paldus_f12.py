from params import *
import paldus_classes
from paldus_classes import ugg
from copy import deepcopy
from paldus_classes import arithmetic_string
from paldus_classes import disambiguate
from paldus_classes import pair_permutations
import math
from paldus_classes import integrate
from paldus_classes import virtual, occupied, general
from fortran_code import *
from fortran_code_f12 import *
from templates import *
from random import shuffle
from factor import *
from paldus_basic import *
import sys
import time
import pickle
from joblib import Parallel, delayed
from paldus_latex import *
from datetime import datetime
from paldus_disconnected import *
import os, sys

import cProfile
MAX_TERMS_PER_BATCH = 150
MAX_TERMS_PER_BATCH_OUTER = 800
MEM_DISK_THRESH = 5
MEM_LITTLE_THRESH = 0.5

##
# @file paldus_f12.py
#
# @brief Generates various expressions in f12 approach
#
# @section description_paldus_f12 Description
# generates ccsd-f12 according to Shiozaki paper DOI: 10.1039/b803704n
# 
#
# @section libraries_sensors Libraries/Modules
# - random standard library (https://docs.python.org/3/library/random.html)
#   - Access to randint function.
#
#
# @section author_sensors Author(s)
# - Created by Aleksandra Tucholska
#
#
#


def execute_ccsd_f12(theory, pick, arg, canonical):
    ## dupa
    ##
    """! 
    @param theory - 'ccsd' or 'cc3'\\
    @param param pick - 'dump' - compute integrals and dump to the disk.
       'load' - load integrals from the disk and maniputate, simplity, rearrange.
    @param  arg - 't1', 't2', or 'tf' - this name will be used as prefix in all file names.
    @param canonical - if True, canonical orbitals are used, f_a^i = 0.
    """
#    r1 = evaluate(hamiltonian_comp_fock, t2fa)
#    r1 = evaluate(fock_exp1c, t2fa)
#    r1 = evaluate(hamiltonian_comp, t1, t1).scale(0.5)

    open_flags = ""
    if pick == "dump":
        open_flags = "wb"
    elif pick == "load":
        open_flags = "rb"


    if pick == "dump":
        # r1 = hamiltonian_comp + evaluate(hamiltonian_comp, t1)+ evaluate(hamiltonian_comp, t2)+ evaluate(hamiltonian_comp, t2fa)\
        #     + evaluate(hamiltonian_comp, t1, t2) + evaluate(hamiltonian_comp, t1, t2fa)+ evaluate(hamiltonian_comp, t1, t1).scale(0.5)+\
        #     evaluate(hamiltonian_comp, t1, t1, t1).scale(0.16666666666666) \
        #     + evaluate(hamiltonian_comp, t1, t1, t2).scale(0.5) + evaluate(hamiltonian_comp, t1, t1, t2fa).scale(0.5)\
        #     + evaluate(hamiltonian_comp, t2, t2fa).scale(0.5)\
        #     + evaluate(hamiltonian_comp, t2, t2).scale(0.5) + evaluate(hamiltonian_comp, t2fa, t2fa).scale(0.5) \
        #     + evaluate(hamiltonian_comp, t1, t1, t1, t1).scale(0.0416666666)

        # r2 = hamiltonian_comp + evaluate(hamiltonian_comp, t1)+ evaluate(hamiltonian_comp, t2)+ evaluate(hamiltonian_comp, t2fa)\
        #     + evaluate(hamiltonian_comp, t1, t2) + evaluate(hamiltonian_comp, t1, t2fa)+ evaluate(hamiltonian_comp, t1, t1).scale(0.5)+\
        #     evaluate(hamiltonian_comp, t1, t1, t1).scale(0.16666666666666) \
        #     + evaluate(hamiltonian_comp, t1, t1, t2).scale(0.5) + evaluate(hamiltonian_comp, t1, t1, t2fa).scale(0.5)\
        #     + evaluate(hamiltonian_comp, t2, t2fa).scale(0.5)\
        #     + evaluate(hamiltonian_comp, t2, t2).scale(0.5) + evaluate(hamiltonian_comp, t2fa, t2fa).scale(0.5) \
        #     + evaluate(hamiltonian_comp, t1, t1, t1, t1).scale(0.0416666666)

        # rgemi = hamiltonian_comp + evaluate(hamiltonian_comp, t1)+ evaluate(hamiltonian_comp, t2)+ evaluate(hamiltonian_comp, t2fa)\
        #     + evaluate(hamiltonian_comp, t1, t2) + evaluate(hamiltonian_comp, t1, t2fa)+ evaluate(hamiltonian_comp, t1, t1).scale(0.5)+\
        #     evaluate(hamiltonian_comp, t1, t1, t1).scale(0.16666666666666) \
        #     + evaluate(hamiltonian_comp, t1, t1, t2).scale(0.5) + evaluate(hamiltonian_comp, t1, t1, t2fa).scale(0.5)\
        #     + evaluate(hamiltonian_comp, t2, t2fa).scale(0.5)\
        #     + evaluate(hamiltonian_comp, t2, t2).scale(0.5) + evaluate(hamiltonian_comp, t2fa, t2fa).scale(0.5) \
        #     + evaluate(hamiltonian_comp, t1, t1, t1, t1).scale(0.0416666666)


        r0 = hamiltoniant_comp + evaluate(hamiltoniant_comp, t2)+ evaluate(hamiltoniant_comp, t2fa)
        
        r1 = hamiltoniant_comp + evaluate(hamiltoniant_comp, t2)+ evaluate(hamiltoniant_comp, t2fa)
        
        r2 = hamiltoniant_comp + evaluate(hamiltoniant_comp, t2)+ evaluate(hamiltoniant_comp, t2fa)\
            + evaluate(hamiltoniant_comp, t2, t2fa).scale(0.5)\
            + evaluate(hamiltoniant_comp, t2fa, t2).scale(0.5)\
            + evaluate(hamiltoniant_comp, t2, t2).scale(0.5) + evaluate(hamiltoniant_comp, t2fa, t2fa).scale(0.5)

        
        # r2 = hamiltoniantc + evaluate(hamiltoniantc, t2)+ evaluate(hamiltoniantc, t2fa)\
        #     + evaluate(hamiltoniantc, t2, t2fa).scale(0.5)\
        #     + evaluate(hamiltoniantc, t2fa, t2).scale(0.5)\
        #     + evaluate(hamiltoniantc, t2, t2).scale(0.5) + evaluate(hamiltoniantc, t2fa, t2fa).scale(0.5)

        

        rf =  hamiltoniant_comp + evaluate(hamiltoniant_comp, t2)+ evaluate(hamiltoniant_comp, t2fa)\
            + evaluate(hamiltoniant_comp, t2, t2fa)\
            + evaluate(hamiltoniant_comp, t2fa, t2fa).scale(0.5)

        # to jest bez operatora focka
        # ale nie moge tak zrobic bo calki z libinta sa z operatorem focka
        # rf =  hamiltoniantc + evaluate(hamiltoniantc, t2)+ evaluate(hamiltoniantc, t2fa)\
        #     + evaluate(hamiltoniantc, t2, t2fa)\
        #     + evaluate(hamiltoniantc, t2fa, t2fa).scale(0.5)

        
        # r1 = hamiltoniant + evaluate(hamiltoniant, t2)
        # r2 = hamiltoniant + evaluate(hamiltoniant, t2) +evaluate(hamiltoniant, t2, t2).scale(0.5)
        
        #                                                                                                                                                                 
        # 'Molecular electronic structure theory' vol 2,  Trygve Helgaker                                                                                                 
        # p. 692 (170 evince) eq 13.7.58 and 13.7.59                                                                                                                      
        # Omega vector function in biorthogonal basis.                                                                                                                    
        #                                                                                                                                                                
        # -------------------------------------------------                                                                                                              
        # print('to jest r1', len(r1))
        # for x in r1:
        #     print(x)
        # print('')


        pick_f12 = open('./pickle/amp_{arg}_f12.pkl'.format(arg=arg), open_flags)
        
        if arg == "energy":
            rint1 = r0.integrate()
            pickle.dump(rint1, pick_f12)
        elif arg == "t1":
            rint1 = r1.integrate(bra = ['a', 'i'], braspin = ['s']).scale(0.5)
            print('wynik')
            for x in rint1:
                print(x)
            pickle.dump(rint1, pick_f12)
        elif arg == "t2":
            print('r21')
            start = time.time()
            rint21  = r2.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s', 's']).scale(1.0/3.0)
            end = time.time()
            print('czas byl', end-start)
            print('r22')
            rint22 = r2.integrate(bra = ['a', 'j', 'b', 'i'], braspin = ['s', 's']).scale(1.0/6.0)
            rint2 = (rint21 + rint22)
            pickle.dump(rint2, pick_f12)
            # pickle.dump(rint21, pick_f12)
        elif arg == "tf":
            k =1                                                                                                                   
            r3 = arithmetic_string()                                                                                       
            for x in rf:                                                                                                            
                print(k, x)                                                                                                                 
                disambiguate(x, gemi)                                                                                                         
                y = x.fromleft(gemi)                                        
                r3 = r3 + arithmetic_string(y)                                                             
                print(k,  y)                                                                                          
                print('')                                                                                                             
                k += 1                                                                                                                               
            print(len(r3))   

            print('oto gemi')                                                                                   
            for x in r3:                                                                     
                print(x)                                                                                          
            print('')
            print('to jest po calkowaniu')

            rint3 = r3.integrate()
        
            pickle.dump(rint3, pick_f12)

        sys.exit(0)
    elif pick == "load":

        load_f12 = open('./pickle/amp_{arg}_f12.pkl'.format(arg=arg), open_flags)
        rint1 = pickle.load(load_f12)


        if arg == "energy":
            prefix = "erg"
        elif arg == "t1":
            prefix = "t1"
        elif arg == "t2":
            prefix = "t2"
        elif arg == "tf":
            prefix = "tf"

    print("I'm removing all instances of ff{akbl}, where a and b are both wirtual indices")
    print("According to definition of F in Shiozaki article, these instances are equal to zero")
    for x in rint1:
        for j in range(0, len(x.coefficient)):
            y = x.coefficient[j]
            if y == F12_TWOEL or y ==F12_TWOEL_COMBO:
                n_of_virt = 0
                for z in x.coefficient_idx[j]:
                    if z in virtual:
                        n_of_virt +=1
                if n_of_virt >=2:
                    x.num_factor =0
#                    print('usuwam', x)
    print('oto rint')
    k = 0
    for x in rint1:
        print(k, x)
        k += 1

    # print('wydrukuj wybrane rint')

    # rra = arithmetic_string()
    # for x in rint1:
    #     if len(x.summation) == 0 or len(x.summation)==len(x.delta):
    #         rra = rra + arithmetic_string(x)
    #         print(x)
    # rs = simplify(rra)
    # print('rsrs')
    # for x in rs:
    #     print(x)
    # sys.exit(0)
        
    print('koniec rint1', len(rint1))

    rsimp  =simplify(rint1)
    print('oto rsimp')
    k = 0
    for x in rsimp:
        print(k, x)
        k += 1


    rint1 = rsimp
    print('regroup', len(rint1))     
    k = 0
    for x in rint1:
        print(k, x)
        k += 1
    
    print('zaczynam sprawdzanie L i F')
    rgroup = find_L(rint1, 0)                                                                                                       
    print(len(rgroup), 'RGROUPPPPPPP', '1111')
    rgroup2 = find_F(rgroup, 0)
    print('dupor')
    print(len(rgroup2), '2222')
    rgroup3 = find_L(rgroup2, 1)
    print(len(rgroup3), '3333')
    rgroup4 = find_F(rgroup3, 1)
    print(len(rgroup4), '4444')
    rgroup5 = find_L(rgroup4, 0)
    print(len(rgroup5), '5555')
    rgroup6 = find_F(rgroup5, 0)

    rgroup = deepcopy(rgroup6)
    print('po find L i F', len(rgroup))                                                                                                               

    # rgroup = deepcopy(rint1)
    k = 0
    for x in rgroup:                                                                                                        
        print(k, x)                                                                                                                
        k += 1
    print('')                                                                                                 

    rsimp2 = simplify(rgroup)                                                                                                         
    print('po rsimp ost', len(rsimp2))                                                                                                
    for x in rsimp2:
        print(x)                                                                                                                                  
    print('')                     

    rint1 = rsimp2
    for x in rint1:
        print(x)
    print('')
    rint1.exec_delta()
#    rint2.exec_delta()

    rsimp1 = simplify(rint1)
 #   rsimp2 = simplify(rint2)


    print('po rsimp przed cleanup', len(rsimp1))
    z = 0
    for x in rsimp1:
        print(z, x)
        z += 1


    # for x in rsimp1:
    #     x.optimize_f12()
    # for x in rsimp2:
    #     x.optimize_f12()
    #
    # Fixme: Add tunable elimination of non-diago
    # #                                                                                                                                          
    # # DELETE ALL TERMS WITH NONDIAGONAL FOCK MATRIX ELEMENTS
    # #
    # for x in rsimp:   
    #     for y in range(0, len(x.coefficient)):  
    #         if x.coefficient[y] == FOCK_MATRIX:        
    #                 x.num_factor = 0.0                                                                                                                             
    rsimp1.cleanup()
    # rsimp2.cleanup()

    # function_template_ccsd(rsimp1, 1)
    # function_template_ccsd(rsimp2, 2)
    
    # if arg == "t2" or arg == "tf":
    #     res = find_permutations(rsimp1)
    #     print('po find permutations')
    #     k = 0
    #     for p in res:
    #         print(k, p)
    #         k += 1
    #     print('')
    # else:
    res = rsimp1
    
    print('po rsimp po cleanup', len(rsimp1))
    z = 0
    for x in rsimp1:
        print(z, x)
        z += 1

    """! Identify intermediates
    """
    ##! Identify intermediates

    # res = rsimp1
    print(' a teraz identify interm')
    res2 = identify_interm_V_f12(res)
    print('po find V')
    k = 0
    for p in res2:
        print(k, p)
        k += 1
    print('')
    res3 = identify_interm_X_f12(res2)
    print('po find X')
    k = 0
    for p in res3:
        print(k, p)
        k +=  1
    print('')
    res4 = identify_interm_B_f12(res3)
    print('po find B')
    k = 0
    for p in res4:
        print(k, p)
        k += 1
    print('')
    res5 = identify_interm_P_f12(res4)
    print('po find P')
    k = 0
    for p in res5:
        print(k, p)
        k += 1
    print('')
    s = arithmetic_string()
    for x in res5:
        s = s + arithmetic_string(x)

    print('')
    ss = s.cabstransform()

    print('po cabs transform')
    z = 0
    for x in ss:        
        print(z, x)
        z += 1

        
    #-----Delete all nondiagonal elements of Fock matrix if canonical MO are used (Brillouin Condition) and also if GBC and EBC
    if canonical == True:
        for x in range(0, len(ss)):
            for i in range(0, len(ss[x].coefficient)):
                coef = ss[x].coefficient[i]
                if coef == FOCK_MATRIX or coef == FOCK_MATRIX_TRANS:                    
                    if ss[x].coefficient_idx[i][0] in virtual and ss[x].coefficient_idx[i][1] in virtual or \
                       ss[x].coefficient_idx[i][0] in CABS and ss[x].coefficient_idx[i][1] in CABS or \
                       ss[x].coefficient_idx[i][0] in occupied and ss[x].coefficient_idx[i][1] in occupied:
                        ss[x].delta.append([ss[x].coefficient_idx[i][0], ss[x].coefficient_idx[i][1]])
                        stay = True
                    else:
                        stay = False
                    if stay == False:
                        ss[x].num_factor = 0

    ss.exec_delta()
    ss.cleanup()

    print('po usunieciu fock')
    z = 0
    for x in ss:
        print(z, x)
        z+=1
                       
                    

    res6 = identify_interm_Ft_f12(ss)
    print('i to jest ostateczny wyniczor przed simp')
    res_ost = arithmetic_string()
    for x in res6:
        res_ost = res_ost + arithmetic_string(x)

    rsimp_ost = simplify(res_ost)

    res_ost = rsimp_ost

    print('i to jest ostateczny wyniczor po symplifikacji')    
    k = 0
    for x in res_ost:
        print(k, "   ", x)
        k += 1

    print('')

    print('koszty bez intermediates')

    cost_old = 0
    z = 0
    idx = 0
    costidxm = [0,0,0]
    for x in res_ost:
        memidx, costidx = memcost(x)
        cost_all = 888**(costidx[0]) * 200**(costidx[1]) * 30**(costidx[2])
        if cost_all > cost_old:
            cost_old = deepcopy(cost_all)
            idx = deepcopy(z)
            costidxm = deepcopy(costidx)
        print(z, costidx[0], ',', costidx[1], ',', costidx[2], costidx[0]+costidx[1]+costidx[2], cost_all)
        z += 1
    vspace(0)
    cstr = f"c^({costidxm[0]})v^({costidxm[1]})o^({costidxm[2]})"
    print('najdrozsza petla bez interm to', idx, cost_old/(10**9), cstr)

    vspace(0)

    print('')


    k = 0
    for x in res_ost:
        oc = 0
        vt = 0
        cbs = 0
        excl = []
        for y in x.summation:
            for idx in y:
                if idx not in excl:
                    if idx in occupied:
                        oc += 1
                    elif idx in virtual:
                        vt += 1
                    elif idx in CABS:
                        cbs +=1
                    excl.append(idx)
        s = "o^{oc}v^{vt}v'^{cbs}".format(oc=oc, vt=vt, cbs=cbs)
        s1 = "o^{oc}v^{vt}v'^{cbs}".format(oc=oc+2, vt=vt+2, cbs=cbs)
        #        print(k, "   ", s, s1, "     ", oc+vt+cbs, "    ", oc+vt+cbs+4)
        print("&", x, "\quad      ", s1, "\quad    ", oc+vt+cbs+4, "\\\\")
        k += 1


    if arg == "energy":
        excl = []
    elif arg == "t1":
        excl = ['a', 'i']
    elif arg == "t2":
        excl = ['a', 'i', 'b', 'j']
    elif arg == "tf":
        excl = ['i', 'j', 'k', 'l']


    print('')
    print('res ost')
    for x in res_ost:
        print(x)
    
        
    basket_outer, basket_nointerm, list_of_int, n_max, all_hash, mem_dict, interm_fx, xfx_dict, list_of_int_xfx  = factorize_driver(res_ost, excl, COST=6)
    print('basket_outer', basket_outer)
    print('basket_nointerm', basket_nointerm)
    intermediates_to_fortran(theory, prefix, basket_outer, basket_nointerm, list_of_int, n_max, all_hash, mem_dict, interm_fx, xfx_dict)


    # basket_outer, list_of_int, n_max, all_hash, mem_dict, interm_fx, xfx_dict  = factorize_driver(res_ost, excl=['a', 'i', 'b', 'j'], COST=6)
    # intermediates_to_fortran(list_of_int, basket_outer, n_max, all_hash, mem_dict, interm_fx, xfx_dict)
    

    sys.exit(0)



def evaluate_s_f12_operator_separate_functions(pick):

    now = datetime.now()
    dt_string = now.strftime("%B-%d-%Y-%H-%M-%S")
    date = "date and time = {dt}".format(dt=dt_string)
    print(dt_string)

    os.system(f"touch latex/S_op/S_operators_f12{dt_string}.tex")
    f = open(f"latex/S_op/S_operators_f12{dt_string}.tex", 'w')
    s_preamble = tex_preamble
    s_beginning = """
    \\begin{document} \n"""
    f.write(s_preamble)
    f.write(s_beginning)

    open_flags = ""
    if pick == "dump":
        open_flags = "wb"
    elif pick == "load":
        open_flags = "rb"


    # ev_s13 = evaluate(t1c, t2) + evaluate(t1c, t2fa)

    # print('evvvv')
    # for x in ev_s13:
    #     print(x)
    # print('po evvv')
    # s13 = int_and_simp(ev_s13, 1)

    # print('S_1(2)')
    # for x in s13:
    #     print(x)
    # print('')
    # sys.exit(0)

    if pick == "dump":

        ev_s12 = arithmetic_string(t1)
        ev_s13 = evaluate(t1c, t2) + evaluate(t1c, t2fa) 

#        ev_s13 =  evaluate(t1c, t2fa) 
        
        ev_s22 = arithmetic_string(t2) + arithmetic_string(t2fa)

        t2fac = deepcopy(t2fa)
        t2fac.transpose()
        
        ev_s23 = evaluate(t2c, t2, t2).scale(0.5) + evaluate(t2c, t2, t2fa).scale(0.5) + evaluate(t2c, t2fa, t2).scale(0.5) \
            + evaluate(t2c, t2fa, t2fa).scale(0.5) + evaluate(t2fac, t2fa, t2).scale(0.5) + evaluate(t2fac, t2, t2fa).scale(0.5)\
            + evaluate(t2fac, t2fa, t2fa).scale(0.5) + evaluate(t2fac, t2, t2).scale(0.5)

        # ev_s23 = evaluate(t2c, t2, t2fa).scale(0.5)

        # ev_s34 = evaluate(t1c, t2, t2).scale(0.5) + evaluate(t1c, t2, t2fa).scale(0.5) + evaluate(t1c, t2fa, t2).scale(0.5) \
        #     + evaluate(t1c, t2fa, t2fa).scale(0.5)



        print('ev_s12 przed calkowaniem')
        k = 1
        for x in ev_s12:
            print(k, x)
            k += 1
        print('')

        print('ev_s13 przed calkowaniem')
        k = 1
        for x in ev_s13:
            print(k, x)
            k += 1
        print('')

        # print('ev_s22 przed calkowaniem')
        # k = 1
        # for x in ev_s22:
        #     print(k, x)
        #     k += 1
        # print('')
        
        print('ev_s23 przed calkowaniem')
        k = 1
        for x in ev_s23:
            print(k, x)
            k += 1
        print('')

        # print('ev_s34 przed calkowaniem')
        # k = 1
        # for x in ev_s34:
        #     print(k, x)
        #     k += 1
        # print('')

        # s12 = int_and_simp(ev_s12, 1)
        # ev_s13a = arithmetic_string(ev_s13[0]) < - to bylo do debug
        s13 = int_and_simp(ev_s13, 1, complete=True)
        s22 = int_and_simp(ev_s22, 2, complete=True)
        s23 = int_and_simp(ev_s23, 2, complete=True)
        # s34 = int_and_simp(ev_s34, 3, complete=True)

        print('S_1(2)')
        for x in s22:
            print(x)

        print('')
        print('S_1(3)')
        for x in s13:
            print(x)
        print('')

        s_v = arithmetic_string()
        s_c = arithmetic_string()

        for x in s13:
            x.clear_fixed()
            x.operator_idx.append(["α", "i"])
            x.operator_type.append("s")
            x.summation.append("α")
            x.summation.append("i")
            x.exec_delta()
            cc = 0
            for i in range(0, len(x.operator_idx)):
                for j in x.operator_idx[i]:
                    if j in completev:
                        cc += 1
            if cc == 1:
                s_c = s_c + arithmetic_string(x)
            elif cc == 0:
                s_v = s_v + arithmetic_string(x)

        print('s_c')
        z = 0
        for x in s_c:
            print(z, x)
            z+= 1
        print('s_v')
        z = 0
        for x in s_v:
            print(z, x)
            z+= 1


        # print('S_2(2)')
        # for x in s22:
        #     print(x)
        # print('')
        print('S_2(3)')
        for x in s23:
            print(x)
        print('')

        print('S_2(3) po')
        s_vv = arithmetic_string()
        s_vc = arithmetic_string()
        s_cc = arithmetic_string()
        for x in s23:
#            print(fixed)
            x.clear_fixed()
            x.operator_idx.append(["α", "i"])
            x.operator_idx.append(["β", "j"])
            x.operator_type.append("s")
            x.operator_type.append("s")
            x.summation.append("α")
            x.summation.append("β")
            x.summation.append("i")
            x.summation.append("j")
            x.exec_delta()
            cc = 0
            for i in range(0, len(x.operator_idx)):
                for j in x.operator_idx[i]:
                    if j in completev:
                        cc += 1
            if cc == 2:
                s_cc = s_cc + arithmetic_string(x)
            elif cc == 1:
                s_vc = s_vc + arithmetic_string(x)
            elif cc == 0:
                s_vv = s_vv + arithmetic_string(x)
        print('s_cc')
        z = 0
        for x in s_cc:
            print(z, x)
            z+= 1
        print('s_vc')
        z = 0
        for x in s_vc:
            print(z, x)
            z+= 1
        print('s_vv')
        z = 0
        for x in s_vv:
            print(z, x)
            z+= 1

        # print('S_3(4)')
        # for x in s34:
        #     print(x)
        # print('')

        s1 = s13# s12 + s13

        #s2 = s23#s22 + s23

        # s3 = s34
        pick_s1_f12_c = open('./pickle/amp_s1_f12_c.pkl', open_flags)
        pick_s1_f12_v = open('./pickle/amp_s1_f12_v.pkl', open_flags)
                
        pick_s2_f12_vv = open('./pickle/amp_s2_f12_vv.pkl', open_flags)
        pick_s2_f12_vc = open('./pickle/amp_s2_f12_vc.pkl', open_flags)
        pick_s2_f12_cc = open('./pickle/amp_s2_f12_cc.pkl', open_flags)

        # pick_s2_f12 = open('./pickle/amp_s2_f12.pkl', open_flags)
        # pick_s3_f12 = open('./pickle/amp_s3_f12.pkl', open_flags)
        

        pickle.dump(s_c, pick_s1_f12_c)
        pickle.dump(s_v, pick_s1_f12_v)
        
        pickle.dump(s_vv, pick_s2_f12_vv)
        pickle.dump(s_vc, pick_s2_f12_vc)
        pickle.dump(s_cc, pick_s2_f12_cc)
        # pickle.dump(s2, pick_s2_f12)
        # pickle.dump(s3, pick_s3_f12)
        sys.exit(0)        
    elif pick == "load":
        
        load_s1_f12 = open('./pickle/amp_s1_f12.pkl', open_flags)
        s1 = pickle.load(load_s1_f12)

        load_s2_f12 = open('./pickle/amp_s2_f12.pkl', open_flags)
        s2 = pickle.load(load_s2_f12)

        load_s3_f12 = open('./pickle/amp_s3_f12.pkl', open_flags)
        s3 = pickle.load(load_s3_f12)

    print("I'm removing all instances of ff{akbl}, where a and b are both wirtual indices")
    print("According to definition of F in Shiozaki article, these instances are equal to zero")

    # print('')
    # for x in s3:
    #     print(x)
    # sys.exit(0)
              
    
    print('po s1 load')
    for x in s1:
        print(x)
        
        for j in range(0, len(x.coefficient)):
            y = x.coefficient[j]
            if y == F12_TWOEL or y ==F12_TWOEL_COMBO:
                n_of_virt = 0
                for z in x.coefficient_idx[j]:
                    if z in virtual:
                        n_of_virt +=1
                if n_of_virt >=2:
                    print('ten usune', x)
                    x.num_factor =0

    for x in s2:
        for j in range(0, len(x.coefficient)):
            y = x.coefficient[j]
            if y == F12_TWOEL or y ==F12_TWOEL_COMBO:
                n_of_virt = 0
                for z in x.coefficient_idx[j]:
                    if z in virtual:
                        n_of_virt +=1
                if n_of_virt >=2:
#                    print('ten usune', x)
                    x.num_factor =0
                # else:
                #     print('ten zostaje', x)

    print('---------------------------------------')
    print('')
    for x in s3:
        for j in range(0, len(x.coefficient)):
            y = x.coefficient[j]
            if y == F12_TWOEL or y ==F12_TWOEL_COMBO:
                n_of_virt = 0
                for z in x.coefficient_idx[j]:
                    if z in virtual:
                        n_of_virt +=1
                if n_of_virt >=2:
#                    print('ten usune', x)
                    x.num_factor =0
                # else:
                #     print('ten zostaje', x)




    rsimp_s1 = simplify(s1)
    rsimp_s2 = simplify(s2)
    rsimp_s3 = simplify(s3)

    rsimp_s1.exec_delta()
    rsimp_s2.exec_delta()
    rsimp_s3.exec_delta()

    rsimp_s1.cleanup()
    rsimp_s2.cleanup()
    rsimp_s3.cleanup()



    print('amplitudy s1')
    for h in rsimp_s1:
        print(h)
    print('')


    sfile = "S2_temp.tex"
    os.system(f"touch latex/{sfile}")
    f2 = open(f"latex/{sfile}", 'w')
    f2.write(tex_preamble)
    f2.write("\\begin{document}\n\equa{\n")

    print('amplitudy s2')
    k = 1
    for h in rsimp_s2:
        print(k, h)
        k+= 1

    sys.exit(0)
    print('amplitudy s2')
    for h in rsimp_s2:
        f2.write(str(h)+"\\\\&\n")
        print(h)
    print('')
    f2.write("}\n\end{document}")
    f2.close()
    
    print('amplitudy s3')
    for h in rsimp_s3:
        print(h)
    print('')


    res_s1x = identify_interm_X_f12(rsimp_s1)
    res_s2x = identify_interm_X_f12(rsimp_s2)
    res_s3x = identify_interm_X_f12(rsimp_s3)

    # print('s przed czymkolwiek, tylko calkowanie')
    # for x in res_s2x:
    #     print(x)
    # sys.exit(0)

    
    print('po find X s1')
    k = 0
    for p in res_s1x:
        print(k, p)
        k +=  1
    print('')
    print('po find X s2')
    k = 0

    sfile = "S2_tempX.tex"
    os.system(f"touch latex/{sfile}")
    f2 = open(f"latex/{sfile}", 'w')
    f2.write(tex_preamble)
    f2.write("\\begin{document}\n\equa{\n")

    print('')
    print('po find X s3')

    for p in res_s2x:
        f2.write(str(p)+"\\\\&\n")
        print(k, p)
        k +=  1
    f2.write("}\n\end{document}")
    f2.close()

    print('')
    print('po find X s3')
    k = 0
    for p in res_s3x:
        print(k, p)
        k +=  1
    print('')

    print('ppp')
    s1 = arithmetic_string()
    for x in res_s1x:
        print(x)
        s1 = s1 + arithmetic_string(x)

    s2 = arithmetic_string()
    for x in res_s2x:
        s2 = s2 + arithmetic_string(x)

    s3 = arithmetic_string()
    for x in res_s3x:
        s3 = s3 + arithmetic_string(x)

    print('przed cabs transform s1')
    print('')
    cabs_s1 = s1.cabstransform()
    print('po cabs transform s1 s1')
    z=0
    for x in cabs_s1:        
        print('aa', z, x)
        z += 1

    cabs_s2 = s2.cabstransform()

    cabs_s3 = s3.cabstransform()
    cabs_s1.cleanup()
    cabs_s2.cleanup()

    print('po cabs transform s1')
    z = 0
    for x in s1:
        print(x)
    print('')
    for x in cabs_s1:        
        print(z, x)
        z += 1
 
    print('po cabs transform s2', len(cabs_s2))
    z = 0
    for x in cabs_s2:        
        print(z, x)
        z += 1

    # print('po cabs transform s3')
    # z = 0
    # for x in cabs_s3:        
    #     print(z, x)
    #     z += 1
    
    res_s1_ft = identify_interm_Ft_f12(cabs_s1)
    print('ident ft t2')
    for x in res_s1_ft:
        print(x)
        
    res_s2_ft = identify_interm_Ft_f12(cabs_s2)
    res_s3_ft = identify_interm_Ft_f12(cabs_s3)

    print('to po find Ft s1')
    
    res_ost_s1 = arithmetic_string()
    for x in res_s1_ft:
        res_ost_s1 = res_ost_s1 + arithmetic_string(x)

    print('to po find Ft s2', len(res_s2_ft))

    res_ost_s2 = arithmetic_string()
    for x in res_s2_ft:
        res_ost_s2 = res_ost_s2 + arithmetic_string(x)

#    print('to po find Ft s3')

    res_ost_s3 = arithmetic_string()
    for x in res_s3_ft:
        res_ost_s3 = res_ost_s3 + arithmetic_string(x)


    rsimp_ost_s1 = simplify(res_ost_s1)
    s1_preamble = """
    \\begin{center}
    Wzory orbitalne do 3 rzędu MBPT na $S_1(f12) = T_1 + \hat{\mathcal{P}}_1\left ([T_1^{\dagger}, T_2+T_2'] \\right )  $\\
    \\end{center}
    """
    #f.write(s1_preamble)
    #latex_with_cost(f, 1, rsimp_ost_s1, name=False, f12=True)    
    print('SYMPLI-SPRAWDZAM')
    rsimp_ost_s2 = simplify(res_ost_s2)
    s2_preamble = """
    \\newpage\n
    \\begin{center}
    Wzory orbitalne do 3 rzędu na $S_2(f12) = T_2 + T_2'+ \\frac12\hat{\mathcal{P}}_2\left ([[(T_2+T_2')^{\dagger}, T_2+T_2'], T_2+T_2'] \\right ) $ zapisane w formacie: wyraz, $V^iv^jo^k$ koszt obliczeniowy podzielony na $V$ - indeksy CABS, $v$ - indeksy wirtualne, $o$ - indeksy occupied, całkowity koszt obliczeniowy $N^{m}$.\\
    \\end{center}
    """
    #f.write(s2_preamble)
#    latex_with_cost(f, 1, rsimp_ost_s2, name=False, f12=True)
    #just_latex(f, 1, rsimp_ost_s2, name=False, f12=True)
    

    rsimp_ost_s3 = simplify(res_ost_s3)
    s3_preamble = """
    \\newpage\n
    \\begin{center}
    Wzory orbitalne do 4 rzędu na $S_3(f12) = T_2 + T_2'+ \\frac12\hat{\mathcal{P}}_2\left ([[T_1^{\dagger}, T_2+T_2'], T_2+T_2'] \\right ) $ zapisane w formacie: wyraz, $V^iv^jo^k$ koszt obliczeniowy podzielony na $V$ - indeksy CABS, $v$ - indeksy wirtualne, $o$ - indeksy occupied, całkowity koszt obliczeniowy $N^{m}$.\\
    \\end{center}
    """
    #f.write(s3_preamble)
    #latex_with_cost(f, 1, rsimp_ost_s3, name=False, f12=True)
    #just_latex(f, 1, rsimp_ost_s3, name=False, f12=True)
    

    print('i to jest ostateczny wyniczor po symplifikacji s1')    
    k = 0

    for x in rsimp_ost_s1:
        print(k, "   ", x)
        k += 1

    print('i to jest ostateczny wyniczor po symplifikacji s2', len(rsimp_ost_s2))    
    k = 0
    for x in rsimp_ost_s2:
        fx_ret = x.establish_fixed_and_return()
        # print(k, "   ", x, fx_ret)
        print(x+"\\\\")
        k += 1

    # print('i to jest ostateczny wyniczor po symplifikacji s3')    
    # k = 0
    # for x in rsimp_ost_s3:
    #     print(k, "   ", x)
    #     k += 1

    # print('')

    print('koszty bez intermediates s1')
    for x in rsimp_ost_s1:
        memidx, costidx = memcost(x)
        print(costidx[0], ',', costidx[1], ',', costidx[2], costidx[0] + costidx[1]+costidx[2])

    print('')

    print('koszty bez intermediates s2')
    for x in rsimp_ost_s2:
        memidx, costidx = memcost(x)
        print(x, costidx[0], ',', costidx[1], ',', costidx[2], costidx[0] + costidx[1]+costidx[2])

    print('')

    # print('koszty bez intermediates s3')
    # for x in rsimp_ost_s3:
    #     memidx, costidx = memcost(x)
        
    #     print(costidx[0], ',', costidx[1], ',', costidx[2], costidx[0] + costidx[1]+costidx[2])

    print('')


    excl_s1 = ['a', 'i']
    excl_s2 = ['a', 'i', 'b', 'j']
    excl_s3 = ['a', 'i', 'b', 'j', 'c', 'k']


    # teraz trzeba podzielić s na rózne s w zależnosci od indeksow CABS
    # w klasie s2_vv są te s2, które mają dwa indeksy virt
    # w klasie s2_vC są te s2, które mają dwa indeksy cabs
    # w klasie s2_CC są te s2, które mają jeden indeks virt i jeden cabs
    
    s2_vv, s2_Cv, s2_vC, s2_CC = divide_s_for_f12_classes(rsimp_ost_s2)
    rsimp_ost_s2 = s2_vv + s2_vC + s2_CC
    vspace(0)
    print('część S2, gdzie S2 mają tylko indeksy vv ma dlugosc:', len(s2_vv))
    print('i oto ona')
    k = 1
    for x in s2_vv:
        fx_ret = x.establish_fixed_and_return()
#        print(k, x, fx_ret)
        print(x+"\\\\")
        k += 1
    vspace(1)
    print('część S2, gdzie S2 mają tylko indeksy Cv ma dlugosc:', len(s2_Cv))
    print('i oto ona')
    k = 1
    for x in s2_Cv:
        fx_ret = x.establish_fixed_and_return()
#        print(k, x, fx_ret)
        print(x+"\\\\")
        k += 1
    vspace(1)

    print('część S2, gdzie S2 mają tylko indeksy vC ma dlugosc:', len(s2_vC))
    print('i oto ona')
    k = 1
    for x in s2_vC:
        fx_ret = x.establish_fixed_and_return()
#        print(k, x, fx_ret)
        print(x+"\\\\")
        k += 1
    vspace(1)
    print('część S2, gdzie S2 mają tylko indeksy CC ma dlugosc:', len(s2_CC))
    print('i oto ona')
    k = 1
    for x in s2_CC:
        fx_ret = x.establish_fixed_and_return()
 #       print(k, x, fx_ret)
        print(x+"\\\\")
        k += 1
    vspace(1)
    print('to koniec ss')


    prefix = 's2'
    
    print('bede generowac intermediates dla wyrazow s2, których jest', len(s2))
    #    rsimp_ost_s2 = rsimp_ost_s2[0:40]
#    rsimp_ost_s2 = s2_vC
    vspace(0)
    for x in range(0, len(rsimp_ost_s2)):

        rsimp_ost_s2[x].clear_fixed()
        rsimp_ost_s2[x].establish_fixed()
        
    vspace(1)
    # basket_outer_s2, basket_nointerm_s2, list_of_int_s2, n_max_s2, all_hash_s2, mem_dict_s2, \
    #     interm_fx_s2, xfx_dict_s2, list_of_int_s2_xfx  = factorize_driver(rsimp_ost_s2, excl_s2, COST=6, filel=f)

    basket_outer_s2, basket_nointerm_s2, list_of_int_s2, n_max_s2, all_hash_s2, mem_dict_s2, \
        interm_fx_s2, xfx_dict_s2, list_of_int_s2_xfx  = factorize_driver(rsimp_ost_s2, excl_s2, COST=6, filel=f, s_ampli=True)

    print('baskett')
    for x in basket_outer_s2:
        print(x)

    print('interm length')
    for x in list_of_int_s2:
        print(len(x))

    print('a', len(rsimp_ost_s2), len(basket_outer_s2), len(list_of_int_s2_xfx))
    s2_vv, s2_Cv, s2_vC, s2_CC = divide_s_for_f12_classes(basket_outer_s2)
    # print(list_of_int_s2_xfx)
    vspace(0)
    print('część S2, gdzie S2 mają tylko indeksy vv ma dlugosc:', len(s2_vv))
    print('i oto ona')
    k = 1
    for x in s2_vv:
        fx_ret = x.establish_fixed_and_return()
        #        print(k, x, fx_ret)
        print(x+"\\\\")
        k += 1
    vspace(1)
    print('część S2, gdzie S2 mają tylko indeksy Cv ma dlugosc:', len(s2_Cv))
    print('i oto ona')
    k = 1
    for x in s2_Cv:
        fx_ret = x.establish_fixed_and_return()
        #        print(k, x, fx_ret)
        print(x+"\\\\")
        k += 1
    vspace(1)

    print('część S2, gdzie S2 mają tylko indeksy vC ma dlugosc:', len(s2_vC))
    print('i oto ona')
    k = 1
    for x in s2_vC:
        fx_ret = x.establish_fixed_and_return()
        #        print(k, x, fx_ret)
        print(x+"\\\\")
        k += 1
    vspace(1)
    print('część S2, gdzie S2 mają tylko indeksy CC ma dlugosc:', len(s2_CC))
    print('i oto ona')
    k = 1
    for x in s2_CC:
        fx_ret = x.establish_fixed_and_return()
        #print(k, x, fx_ret)
        print(x+"\\\\")
        k += 1
    vspace(1)
    print('to koniec ss')
    sys.exit(0)

    

    print('bbasket_outer', len(basket_outer_s2))
    for x in range(0, len(basket_outer_s2)):
          print(basket_outer_s2[x])
    print('')
    print('bbasket noninterm s2', len(basket_nointerm_s2))
    for x in basket_nointerm_s2:
        print(x)
    print('llist_of_int', len(list_of_int_s2))
    for x in range(0, len(list_of_int_s2)):
        print(len(list_of_int_s2[x]))
        for y in list_of_int_s2[x]:
            print(y)

    print('hehelo')
    sys.exit(0)
    # prefix = 's3'
    # basket_outer_s3, basket_nointerm_s3, list_of_int_s3, n_max_s3, all_hash_s3, mem_dict_s3, \
    #     interm_fx_s3, xfx_dict_s3, list_of_int_s3_xfx  = factorize_driver(rsimp_ost_s3, excl_s3, COST=6, filel=f)

    f.write("\end{document}")
    f.close()
    print('finito-end')
    sys.exit(0)
    
    intermediates_to_fortran('ccsd', prefix, basket_outer_s2, basket_nointerm_s2, list_of_int_s2, n_max_s2, \
                             all_hash_s2, mem_dict_s2, interm_fx_s2, xfx_dict_s2)




    for x in range(0, len(s21)):
        for y in s21[x].coefficient:
            at_least_one = True
            if y == F12_TWOEL:
                at_least_one = False
                for i in s21[x].coefficient_idx:
                    if i in completev:
                        at_least_one = True
                        break
            if at_least_one == False:
                s21[x].num_factor == 0

    s21.cleanup()


    print('no to zaczynam s22')
    print('len', len(s22))
    for x in range(0, len(s22)):
        for j in range(0, len(s22[x].coefficient)):
            y = s22[x].coefficient[j]
            at_least_one = True
            if y == F12_TWOEL:
                print('tak jest f12_twoel')
                print(s22[x])
                at_least_one = False
                print('at_least_one', at_least_one)
                for i in s22[x].coefficient_idx[j]:
                    print(i)
                    if i in completev:
                        at_least_one = True
                        print('at_least_one', at_least_one)
                        break
            if at_least_one == False:
                print('tak, falsz')
                s22[x].num_factor = 0

    s22.cleanup()
    print('len2', len(s22))


    print('s22 wynik')
    for x in s22:
        print(x)


    factorize_driver(s22, ['a', 'i', 'b', 'j'], 6)




def divide_s_for_f12_classes(rsimp_ost_s2):

    s2_vv = arithmetic_string()
    s2_vC = arithmetic_string()
    s2_Cv = arithmetic_string()
    s2_CC = arithmetic_string()
    
    for x in rsimp_ost_s2:
        # esthablish fixed indices and return them as the list fx_ret
        fx_ret = x.establish_fixed_and_return()
        print('fxret', fx_ret)
        # establish whether S2 has two virtual, two cabs or 1 virt and 1 cabs indices.
        vc = virt_class(fx_ret)
        if vc == 1:
            s2_vv.append(x)
        elif vc == 2:
            print('takdwa')
            if 'A' in fx_ret:
                print('takA')
                s2_Cv.append(x)
            else:
                print('takB')
                s2_vC.append(x)
        elif vc == 3:
            s2_CC.append(x)
        else:
            print(fx_ret)
            print('cos bardzo dziwnego')
            sys.exit(0)
    print('alla')
    print(len(s2_vv), len(s2_Cv), len(s2_vC), len(s2_CC))

    return s2_vv, s2_Cv, s2_vC, s2_CC

def virt_class(fx_ret):

    v = 0
    c = 0
    for x in fx_ret:
        if x in virtual:
            v += 1
        if x in CABS:
            c += 1
    if v == 2 and c == 0:
        return 1
    elif v == 1 and c == 1:
        return 2
    elif c == 2 and v == 0:
        return 3

def factorize_driver(ars, excl=[], COST=[], is_x_excluded = False, filel=[], idx_str='', s_ampli=False, disconnected = False):

    start_time = time.time()
    
    print('FACTORIZE_DRIVER START')
    print('indices excluded from factorization:', excl)
    k = 1
    kkk = 0
    biglevel_super = []
    bigparameters_super = []
    basket_nointerm = []
    idx_interm = []
    idx_noninterm = []
    idxidx = -1
    l = 0
    ars_copy = arithmetic_string()
    big_x_copy = deepcopy(ars)
    toll=0
    allii = 0
    allss = 0
    lw = 0
    li = 0
    for ii, x in enumerate(ars):
        excl_this = x.establish_fixed_and_return()

        #subtract Pp and Pm  - one nie tworza intermediate, zostana dodane po faktoryzacji do czesci outer.


        x_copy = big_x_copy[ii] #deepcopy(x)
        x_copy.coefficient = []
        x_copy.coefficient_idx = []
        x_temp = []
        for i in range(0, len(x.coefficient)):
            if (x.coefficient[i] != "Pp" and x.coefficient[i] != "Pm"):
                x_copy.coefficient.append(x.coefficient[i])
                x_copy.coefficient_idx.append(x.coefficient_idx[i])
            else:
                x_temp.append([x.coefficient[i], x.coefficient_idx[i]])

        l +=1

        t1 = time.time()

        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0
        n5=0
        n6=0
        n7=0
        n8=0
        n9=0
        n10=0
        n11 =0
        n12=0

        nn = [0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12]
        biglevel_best, bigparameters_best, haveinterm = factorize(nn, x_copy, excl_this, COST, is_x_excluded, disconnected=disconnected)
        print(nn)
        t2 = time.time()
        if (abs(t2-t1) > 0.1):
            print('pluszyna', ii, len(ars), t2-t1)
            if ((t2-t1) > lw):
                lw = deepcopy(t2-t1)
                li = ii
            allii +=1
        elif (abs(t2-t1) < 0.1 and abs(t2-t1)>0.01):
            print('snieznina', ii, len(ars), t2-t1)
            allss +=1
            
        toll += (t2-t1)

        if len(x_temp)>0:            
            for j in x_temp:
                x_copy.coefficient.append(j[0])
                x_copy.coefficient_idx.append(j[1])

        ars_copy.append(x_copy)

        idxidx += 1
        if haveinterm == True:
            idx_interm.append(idxidx)
            biglevel_super.append(biglevel_best)
            bigparameters_super.append(bigparameters_best)

            k += 1
            s = 1
            if len(biglevel_best) == 0:
                kkk += 1
            for f in range(0, len(bigparameters_best)):
                z = bigparameters_best[f]
                zz = biglevel_best[f]
            s+= 1
        else:
            basket_nointerm.append(x)
            idx_noninterm.append(idxidx)
            
    
    ars = deepcopy(ars_copy)
    print('sfaktoryzowane')
    print('allii', allii)
    print('lw', lw, li)
    print('allss', allss)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time1: {elapsed_time} seconds")
    print('toll', toll)
#    sys.exit(0)
    start_time = time.time()


    # print('oto sa piramidy')
    # for i in range(0, len(idx_interm)):
    #     print(i, 'interm i jego piramida', ars[idx_interm[i]], len(biglevel_super[i]))
    #     for j in biglevel_super[i]:
    #         print(j)
    #     print('')

    print('')
    print('TERAZ CZYTAM PIRAMIDY')
    k = 0

    basket_super = []
    basket_super_disconnected = []
    interm_dict = {}
    interm_hash = {}
    basket_outer_all = []
    basket_idx_fixed = []
    basket_original_only_with_interm = []

    all_hash = {}
    xfx_dict = {}
    interm_fx = {}
    ni = 1
    lenall = 0

    print('noninterm')
    for x in basket_nointerm:
        print(x)
    print('')

    print('rangeragne', len(idx_interm))
    # print(idx_interm)

    # ten slownik zawiera ogolna calkowita ilosc wystapien danego intermediate w kazdym
    # z zestawow w baskecie    
    interm_dict_all_disc = {}
    # ten slownik zawiera informacje o tym, ze dany intermediate w ogole wystapil w danym
    # baskecie, ale nie wiadomo w ktorym zestawie
    interm_dict_once_per_basket_disc = {}
    
    disconnected_super = []
    #  w tablicy idx_interm sa indeksy tych wyrazow ktore maja intermediate
    #  w tablicy idx_interm_disc sa indeksy tych wyrazow ktore zawieraja disconnected.

    idx_interm_disc = []
    for i in range(0, len(idx_interm)):
        # vspace(0)
        basket_original_only_with_interm.append(ars[idx_interm[i]])
        # print(i, 'teraz taki interm', ars[idx_interm[i]], idx_interm[i], 'nnii', ni)
        # print('mam dla niego następującą ilość piramid:', len(biglevel_super[i]), biglevel_super[i])
        # print()
        x = ars[idx_interm[i]]
        excl_this = x.establish_fixed_and_return()

        worm = []
        outer_worm = []
        idxfx_worm = []

        jjj = 0
        disconnected_mini = []

        interm_dict_mini = {}
        # print('w ramach tego basketa', idx_interm[i])
        for minilevel in biglevel_super[i]:
            disconnected_this = [0]*len(minilevel)
#            print('jjj', jjj, minilevel, len(biglevel_super[i]), [0]*len(minilevel))
            ni, disconn_out = read_pyramid(minilevel, x, worm, outer_worm, idxfx_worm, interm_dict,
                                           all_hash, xfx_dict, interm_fx, \
                                           interm_hash, ni, disconnected_this, excl_this, 6, disconnected)
 #           print('po przeczytaniu piramidy dostaje numer intermediate', ni)
  #          print('Alevel nr', jjj, 'czy byly disconnected?', disconn_out, disconnected_this)
            # if disconn_out:
            #     interm_idx.append(idx_interm[i])
            disconnected_mini.append(disconnected_this)
            jjj += 1
#        print('len(worm)', len(worm), 'powinien byc równy len(minilevel)',len(biglevel_super[i]), 'chyba ze nie bylo' )
        if worm == []:
            print('worm byl zero')
            print('idx interm i', i, idx_interm[i], x)
            sys.exit(0)
            # print('')
        # print('i oto on worm dla wyrazu', x)
        # for hh in range(0, len(worm)):
        #     print(worm[hh], idxfx_worm[hh])
        # print('a oto outer worm')
        # for hh in outer_worm:
        #     print(hh)
    
        basket_super.append(worm)
        disconnected_super.append(disconnected_mini)
        basket_outer_all.append(outer_worm)
        # print('idxfxworm', idxfx_worm)
        # print('worm', worm)
        basket_idx_fixed.append(idxfx_worm)
        # print(i, 'wormworm ', x)
        # for gg in worm:
        #     print('worm', gg)
        #     print('')
        # print('outer worm')
        # for gg in outer_worm:
        #     print('outer_worm', gg)
        #     print('')

    all_len = 0
    for x in basket_super:
        all_len += len(x)
    print('PIRAMIDY PRZECZYTANE', 'len wszystkiego:', all_len,',', len(basket_super), len(basket_outer_all))



    
    k = 0
    print(len(interm_dict))
    print(len(all_hash))
    
    for x in interm_dict:
        print(all_hash[x], interm_dict[x], interm_hash[x])
        k+=1

    sorted_interms = sorted(interm_dict.items(), key=lambda x: x[1], reverse=True)



    # Print the sorted dictionary items
    print('sorteed')
    for key, value in sorted_interms:
        print(value, all_hash[key], interm_hash[key], '                     ', isdisc(interm_hash[key]))
    print('')


    
#-------------------------------to na chwile wykomentowanae
# ale tutaj pozniej nalezy dodawac punkty dla not disconnected

    k = 1
    basket_interm = []
    basket_outer = []
    basket_xfx = []

    # n_of_disconnected = 0
    rest = 0

    print('all baskets')
    for x in range(0, len(basket_super)):
        print('basket nr', x, 'o dlugosci', len(basket_super[x]))
        for y in range(0, len(basket_super[x])):
            for z in range(0, len(basket_super[x][y])):
                print(all_hash[basket_super[x][y][z].binary_hash])
            print('')
                
    basket_points_super = []
    for x in range(0, len(basket_super)):
        print('basket nr', x, 'o dlugosci', len(basket_super[x]))
        # print('disconnected dlan', disconnected_super[x])
        # ts =0
        # for z in disconnected_super[x]:
        #     ts += sum(z)
        # if ts>0:
        #     n_of_disconnected += 1
        # else:
        #     rest += 1
        print('original', basket_original_only_with_interm[x])
        basket_points = []
        for y in range(0, len(basket_super[x])):
            print('outer', basket_outer_all[x][y])
            print('')
            points = 0
            for z in range(0, len(basket_super[x][y])):
                points += interm_dict[basket_super[x][y][z].binary_hash]
                print(x, y, z, basket_super[x][y][z], all_hash[basket_super[x][y][z].binary_hash], interm_dict[basket_super[x][y][z].binary_hash])
            print(points)
            basket_points.append(points)
            print('')
        print(basket_points)
        basket_points_super.append(basket_points)
            #    sys.exit(0)

    print('wszystkich wyrazow bylo' ,len(ars))
    print('tych co nie maja interm', len(basket_nointerm))
    print('tych co maja interm', len(basket_super))
#    print('disconnected', n_of_disconnected)
    print('a roznych interm jest', len(interm_dict))
    print('rest', rest)

    #------------------------ CALCULATE INTERMEDIATE LEVEL ----------------------------------------

    interm_whichlevel_dict = {} # 0 - base interm, 1 - interm that has at least one interm0, 2 - interm that has at least one interm1...
    maxl = 10

    for interm in interm_dict:
        notzero = False
        for coef in interm_hash[interm].coefficient:
            if 'interm' in coef:
                notzero = True
                break
        if not notzero:
            interm_whichlevel_dict[all_hash[interm]] = 0

    for lev in range(1, maxl):
        for interm in interm_dict:
            addthis = False
            maxlev = 0
            for coef in interm_hash[interm].coefficient:
                thislev = 0
                if 'interm' in coef:
                    if coef in interm_whichlevel_dict:
                        addthis = True
                        thislev +=  interm_whichlevel_dict[coef] + 1
                        if thislev>maxlev:
                            maxlev = deepcopy(thislev)
            if addthis:
                interm_whichlevel_dict[all_hash[interm]] = maxlev


    for key in interm_whichlevel_dict:
        print('interm', key, 'is at level', interm_whichlevel_dict[key])

    print(len(basket_original_only_with_interm))
    print(len(basket_super))

          

    # ------------------------------------ wybieranie disconnected -------------------------------------------
    if disconnected:

        print('pprzed disc1')
        for key, value in sorted_interms:
            print(value, all_hash[key], interm_hash[key], '                     ', isdisc(interm_hash[key]))
        print('')


        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"zzElapsed time1: {elapsed_time} seconds")
        start_time = time.time()
        sorted_disc_interms, basket_super_disc, basket_outer_disc, basket_original_disc, \
            sorted_nodisc_interms, basket_super_nodisc, basket_outer_nodisc, basket_original_nodisc, disc_super_new, basket_points_disc, \
            basket_points_nodisc = \
                identify_disconnected(all_hash, basket_super, basket_outer_all, sorted_interms, interm_dict, interm_hash,\
                                      disconnected_super, basket_original_only_with_interm)

        print('len(basket_super)', len(basket_super), len(basket_super_disc), len(basket_super_nodisc))

        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed disconeccted identify time2: {elapsed_time} seconds")

        start_time = time.time()
        old_new_dict = {}
        n_interm = 0

        print('pprzed disc2')
        for key, value in sorted_interms:
            print(value, all_hash[key], interm_hash[key], '                     ', isdisc(interm_hash[key]))
        print('')


        old_new_dict, n_interm = pick_best_basket(sorted_disc_interms, basket_super_disc, basket_outer_disc, \
                                        basket_original_disc, basket_points_disc, interm_dict,  \
                                                  interm_hash, disc_super_new, \
                                                  interm_whichlevel_dict, all_hash, old_new_dict, n_interm, True)

        end_time = time.time()
        elapsed_time = end_time - start_time

        print(f"Elapsed pick best disconnected time303: {elapsed_time} seconds")
#        sys.exit(0)


        # ------------------------------------ wybieranie connected -------------------------------------------

        print('pprzed disc3')
        for key, value in sorted_interms:
            print(value, all_hash[key], interm_hash[key], '                     ', isdisc(interm_hash[key]))
        print('')
        start_time = time.time()
        old_new_dict, n_interm = pick_best_basket(sorted_nodisc_interms, basket_super_nodisc, basket_outer_nodisc, \
                         basket_original_nodisc, basket_points_nodisc, interm_dict,  \
                                                  interm_hash, disc_super_new, interm_whichlevel_dict, \
                                                  all_hash, old_new_dict, n_interm, False)

        end_time = time.time()
        elapsed_time = end_time - start_time

        print(f"Elapsed pick best connected time304: {elapsed_time} seconds")


        sys.exit(0)



    
    #----------------------------------------------------------------------------------------------


    k = 1
    basket_interm = []
    basket_outer = []
    basket_xfx = []

    print('TERAZ WYBIERANIE BASKETA', len(basket_super))
    for x in range(0, len(basket_super)):
        
        basket_points = []
        for y in range(0, len(basket_super[x])):
            points = 0
            for z in range(0, len(basket_super[x][y])):
                points += interm_dict[basket_super[x][y][z].binary_hash]
            basket_points.append(points)

        print('basket nr', x, 'o dlugosci', len(basket_super[x]))
        for y in basket_super[x]:
            for z in y:
                print('w tym baskecie jest zestaw', z, len(y))
            print('albo')
        print('')
        print(basket_points)
        print('--------------------------------')
        print('')

    
    print('TERAZ WYBIERANIE BASKETA', len(basket_super))
    for x in range(0, len(basket_super)):
        print('basket nr', x, 'o dlugosci', len(basket_super[x]))
        for y in basket_super[x]:
            for z in y:
                print('w tym baskecie jest zestaw', z, len(y))
            print('albo')
        print('')
        print('--------------------------------')
        print('')
        basket_points = []
        k += 1
        for y in range(0, len(basket_super[x])):
            points = 0
            print('przyznaje punkty dla zestawu', y)
            for z in range(0, len(basket_super[x][y])):
                print('od tego dodam', basket_super[x][y][z], interm_dict[basket_super[x][y][z].binary_hash])
                points += interm_dict[basket_super[x][y][z].binary_hash]
                print('SPR POINT', points, basket_super[x][y][z], basket_idx_fixed[x][y])
            print('dodaje points', points)
            basket_points.append(points)

 #           print('SPR POINT', points)
#        print('')

        if basket_points != []:
            idx = basket_points.index(max(basket_points))
            print(basket_points)
            print(idx)
            #print(basket_super[x][idx])
            basket_interm.append(basket_super[x][idx])
            basket_xfx.append(basket_idx_fixed[x][idx])
            print('basket_xfx', basket_xfx, basket_interm)
            print(x, 'dodaje do interm-sss')
            www = 0
            for ww in basket_interm:
                www += 1
                print(www, ww[0])
            print('')
            
            for gg in basket_super[x][idx]:
                print(x, gg)
            print('')
            for jj in basket_super[x]:
                for ii in jj:
                    print('to mam do wyboru', ii, 'a wybieram', idx)
                print('')

            if len(basket_outer_all[x][idx].coefficient)==0:
                sys.exit(0)
            basket_outer.append(basket_outer_all[x][idx])
        else:
            print('ars', ars[x])
            basket_interm.append([])
            basket_xfx.append([])                        
            basket_outer.append([])
            print('nie ma')
            sys.exit(0)
    print("TOTOTO", len(basket_interm))


    print('basket outer zaraz po stworzeniugo')
    for x in basket_outer:
        print(x)

    n_max = 0

    for i in range(0, len(basket_interm)):

        if (len(basket_interm[i])) > n_max:

            n_max = deepcopy(len(basket_interm[i]))


    print('OSTATECZNE')
    list_of_int= []
    list_of_int_xfx= []
    list_of_names = []
    list_of_int_mem = []
    list_of_int_disk = []
    list_of_int_mem_names = []
    list_of_int_disk_names = []
    for n in range(0, n_max):
        list_of_names.append([])
        list_of_int.append([])
        list_of_int_xfx.append([])
        list_of_int_mem.append([])
        list_of_int_disk.append([])
        list_of_int_mem_names.append([])
        list_of_int_disk_names.append([])


    mem_dict = {}
    cost_dict = {}
    lenall2 = 0
    la = 0
    lb = 0
    # w tym momencie basket_interm to jest
    # [[interm1, interm2], [interm3], [interm4, interm5, interm6].....] <-- tylko numery nie sa po kolei, tylko losowe
    print('TWORZE LISTE INTERMEDIATES')
    # for i in range(0, len(basket_interm)):
    #     for j in range(0, len(basket_interm[i])):
    #         print('teraz sprawdzam wyraz', basket_interm[i][j])
    #     print('')
    # print('')
    memall = 0
    whichlevelmax = 0 
    for i in range(0, len(basket_interm)):
        whichlevel = 0
        for j in range(0, len(basket_interm[i])):
            lenall2 += 1
            print('teraz sprawdzam wyraz', basket_interm[i][j])
            memidx, costidx = memcost(basket_interm[i][j])
            memreal = (888**(memidx[0]) * 232**(memidx[1]) * 20**memidx[2])*8/(10**9)
            # print('memreal', memreal, memidx)
            memall += memreal
  #          print('COSTIDX', costidx, memreal)
            # Sprawdzam jaki jest koszt obliczeniowy i pamieciowy danego intermediate)
            if memreal > 10.0:
                print(memidx[0], memidx[1])
                print(basket_interm[i][j], memreal)
                print('DUZE GOWNO')
            # Dodaje do slownika koszty danego intermedaite. Bdzie dopasowane do hasha.
            mem_dict[basket_interm[i][j].binary_hash] = memreal
            cost_dict[basket_interm[i][j].binary_hash] = costidx
            
            # teraz sprawdzam ktorego rzedu jest dany intermediate                                                                                            
            for coef in basket_interm[i][j].coefficient:
                if 'interm' in coef:
                    whichlevel += 1
                    break
            if whichlevel > whichlevelmax:
                whichlevelmax = whichlevel
            if basket_interm[i][j] not in list_of_int[whichlevel]:
                print(i, 'dodaje do', whichlevel, len(list_of_int[whichlevel]))
                list_of_int[whichlevel].append(basket_interm[i][j])
                list_of_int_xfx[whichlevel].append(basket_xfx[i][j])
                list_of_names[whichlevel].append(all_hash[basket_interm[i][j].binary_hash])
                print('dany interm i jego imie', basket_xfx[i][j], basket_interm[i][j], all_hash[basket_interm[i][j].binary_hash])
                print(all_hash[basket_interm[i][j].binary_hash], ' = ', basket_interm[i][j])
#                print('dodaje', lenall2)
                la += 1
            else:
                lb += 1

    # n_max = whichlevelmax + 1
    print('memall', memall)

    print(lenall2, la, lb)

    print('koszt noninterm')
    for x in basket_nointerm:
        memidx, costidx = memcost(x)
        memreal = (888**(memidx[0]) * 232**(memidx[1]) * 20**memidx[2])*8/(10**9)
        costoverall = costidx[0]+costidx[1]+costidx[2]
        print(costoverall, costidx, x)

    print('basket outer', len(basket_outer))
    for x in basket_outer:
        print(x)
        
    print('intermediates len')
    for x in range(0, n_max):
        print(len(list_of_int[x]))
    print('')



    print('basket interm')
    for x in range(0, n_max):
        for y in range(0, len(list_of_int[x])):
            print(x, list_of_int[x][y], list_of_int_xfx[x][y], all_hash[list_of_int[x][y].binary_hash])
        print('')
    print('---------------------------------------------------')

    #-----------------------------SORTOWANIE SUM PO KOLEI--------------------------------------------#

    for x in range(0, n_max):
        all_summation_list = []
        interm_mini = []
        interm_mini_xfx = []
        for y in list_of_int[x]:
            y.summation.sort()
            if y.summation not in all_summation_list:
                all_summation_list.append(y.summation)
        all_summation_list = sorted(all_summation_list, key=lambda i: len(i))
        for s in all_summation_list:
            print('ssss', s)
            for y in range(0, len(list_of_int[x])):
                print('no i co to za gowno', y)
                if list_of_int[x][y].summation == s:
                    print('ddoaje')
                    interm_mini.append(list_of_int[x][y])
                    interm_mini_xfx.append(list_of_int_xfx[x][y])
        list_of_int[x] = interm_mini
        list_of_int_xfx[x] = interm_mini_xfx

    print('a po sortowaniu')
    for x in range(0, n_max):
        for y in range(0, len(list_of_int[x])):
            print(list_of_int[x][y], list_of_int_xfx[x][y], all_hash[list_of_int[x][y].binary_hash])
        print('')
    print('---------------------------------------------------')

    #------------------------ZMIANA NAZW INTERMEDIATE NA KOLEJNE-------------------------------------#

    print('TERAZ ZMIENIE NAZWY')
    new_names_dict = {}
    bn = "l"
    unique_dict = []
    # for x in range(0, n_max):
    #     for y in range(0, len(list_of_int[x])):
    #         print(list_of_int[x][y],  all_hash[list_of_int[x][y].binary_hash])
    #     print('')
    # sys.exit(0)

    print('Intermediates przed zmianami nazwy')
    for x in range(0, n_max):
        print('level', x)
        k = 0
        for y in range(0, len(list_of_int[x])):
            print(x, y,list_of_int[x][y], all_hash[list_of_int[x][y].binary_hash]) 

    fixed_in_outer = deepcopy(list_of_int_xfx)

    n_of_real_occurences = {}
    print('plusz')
    bh = 'l'
    for x in range(0, n_max):
        print('level', x)
        k = 0
        unique_dict_mini = identify_unique_interm(list_of_int[x], fixed_in_outer[x], all_hash, 1, new_names_dict, disconnected)
        distinct_values_mini = set(unique_dict_mini.values())
        vspace(0)
        print('REAL number of distinct intermediates')
        vspace(0)
        unique_dict.append(unique_dict_mini)
        print('uniq mini', unique_dict_mini)
        for y in range(0, len(list_of_int[x])):
            print('kalara', x, y, unique_dict[x][y], list_of_int[x][y].coefficient[0], list_of_int[x][y], all_hash[list_of_int[x][y].binary_hash])
            if (unique_dict[x][y] == y):
                new_names_dict[all_hash[list_of_int[x][y].binary_hash]] = bn+str(k)
                print('nowa nazwa to bedzie',   bn+str(k), 'a stara to ', all_hash[list_of_int[x][y].binary_hash])
                k += 1
            else:
                new_names_dict[all_hash[list_of_int[x][y].binary_hash]] =  new_names_dict[all_hash[list_of_int[x][unique_dict[x][y]].binary_hash]]
                print('ten juz byl wiec bedzie sie nazywal', all_hash[list_of_int[x][y].binary_hash], \
                      new_names_dict[all_hash[list_of_int[x][y].binary_hash]], \
                      new_names_dict[all_hash[list_of_int[x][unique_dict[x][y]].binary_hash]])

        bn = bn + 'l'

    for key in new_names_dict:
        print(key, new_names_dict[key])

    print('')
    print('-----------------')
    print('')
    all_hash_new = {}
    for x in range(0, n_max):
        k = 0
        for y in range(0, len(list_of_int[x])):
            #------- zmiana nazwy w all_hash--------
            print('zmiana nazwy w all hash z lewelu', x, 'numer', y)
            print('sssyk', list_of_int[x][y])
            print('all hash', all_hash[list_of_int[x][y].binary_hash])
            print('names', new_names_dict[all_hash[list_of_int[x][y].binary_hash]])
            all_hash[list_of_int[x][y].binary_hash] = new_names_dict[all_hash[list_of_int[x][y].binary_hash]]
            if x > 0:
                print('zmieniam nazwy w lewelu', x)
                #------- dodatkowo zmiana nazwy uzywanych intermediate w ugg
                for j in range(0, len(list_of_int[x][y].coefficient)):

                    if list_of_int[x][y].coefficient[j] in list(new_names_dict.keys()):
                        print('tak', list_of_int[x][y].coefficient[j], 'jest w slowniku',\
                              list_of_int[x][y], 'nowa nazwa to', new_names_dict[list_of_int[x][y].coefficient[j]])
                        list_of_int[x][y].coefficient[j] = new_names_dict[list_of_int[x][y].coefficient[j]]
                        print('i wyszlo', list_of_int[x][y])


    print('Intermediates po zmianach nazw')
    for x in range(0, n_max):
        print('level', x)
        k = 0
        for y in range(0, len(list_of_int[x])):
            print(x, y,list_of_int[x][y], all_hash[list_of_int[x][y].binary_hash]) 

#    sys.exit(0)
    print('OUTER')
    print('new_names_basket_outer', new_names_dict)
    for x in range(0, len(basket_outer)):
        print('przed', basket_outer[x])
        for j in range(0, len(basket_outer[x].coefficient)):
            print('sprawdzam coefficient', basket_outer[x].coefficient[j])
            if basket_outer[x].coefficient[j] in list(new_names_dict.keys()):
                print('trala')
                print('tak', basket_outer[x].coefficient[j], 'jest w slowniku', \
                      basket_outer[x], new_names_dict[basket_outer[x].coefficient[j]])
                basket_outer[x].coefficient[j] = new_names_dict[basket_outer[x].coefficient[j]]
                # print('zzzz', new_names_dict[basket_outer[x].coefficient[j]])
        print('po', basket_outer[x])

    print('outer po zmianach', len(basket_nointerm))
    for x in range(0, len(basket_outer)):
        print(x, basket_outer[x])

    print('i do tego Intermediates po zmianach nazw')
    list_of_int_unique = []
    for x in range(0, n_max):
        print('level', x, len(list_of_int[x]))
        # list_of_int_unique.append()
        k = 0
        for y in range(0, len(list_of_int[x])):
            print(x, y,list_of_int[x][y], all_hash[list_of_int[x][y].binary_hash]) 
            

    print('NOINTERM')
    print(new_names_dict)

    for x in range(0, len(basket_nointerm)):
        print('przed', basket_nointerm[x])
        for j in range(0, len(basket_nointerm[x].coefficient)):                              
            if basket_nointerm[x].coefficient[j] in list(new_names_dict.keys()):
                print('tak', basket_nointerm[x].coefficient[j], 'jest w slowniku', \
                    basket_nointerm[x], new_names_dict[basket_nointerm[x].coefficient[j]])
                basket_nointerm[x].coefficient[j] = new_names_dict[basket_nointerm[x].coefficient[j]]
        print('po', basket_nointerm[x])


            
    for x in new_names_dict:
        print(x, new_names_dict[x])
    print('KONIEC ZMIAN NAZW')

    for x in range(0, n_max):
        for y in range(0, len(list_of_int[x])):
            print(list_of_int[x][y], list_of_int_xfx[x][y], all_hash[list_of_int[x][y].binary_hash])
        print('')
    print('---------------------------------------------------')

    #-----------------------------------------------------------------------------------------------#

    print( 'zmiana slownikow')

    vspace(0)

    # ten fragment kodu jest chyba niepotrzebny. ------------------------
    # indeksy fixed mozna odczytac z intermediatu
    # ------------------------------------------
    # print('stare int')
    # k = 0
    # for x in list_of_int:
    #     for y in x:
    #         print(k, y, all_hash[y.binary_hash])
    #         k+=1
    # print('aaa')
    # list_of_int_xfx = []
    # for x in range(0, n_max):
    #     print('level', x)
    #     k = 0
    #     list_of_int_xfx.append([])
    #     new_list_of_int, new_list_of_fx = identify_unique_interm(list_of_int[x], fixed_in_outer[x], all_hash, 2, new_names_dict)
    #     print('nowe')
    #     k = 0
    #     for y in new_list_of_int:
    #         print(k, y, all_hash[y.binary_hash])            
    #         k+=1
    #     print('')
    #     print('newnew')
    #     for y in range(0, len(list_of_int[x])):
    #         mem_dict[new_list_of_int[y].binary_hash] =  mem_dict[list_of_int[x][y].binary_hash]
    #         cost_dict[new_list_of_int[y].binary_hash] =  cost_dict[list_of_int[x][y].binary_hash]
    #         interm_fx[new_list_of_int[y].binary_hash] =  new_list_of_fx[y]#interm_fx[list_of_int[x][y].binary_hash]
    #         xfx_dict[new_list_of_int[y].binary_hash] = xfx_dict[list_of_int[x][y].binary_hash]
    #         print('pluszek',list_of_int[x][y], new_list_of_int[y], interm_fx[new_list_of_int[y].binary_hash], new_list_of_fx[y])
            
    #     list_of_int[x] = new_list_of_int
    #     list_of_int_xfx[x] = new_list_of_fx
            
        
    print( 'zmiana slownikow2')



    print('TERAZ ZOSTAWIAM JEDYNIE UNIQUE')
    list_of_int_copy = deepcopy(list_of_int)
#    list_of_int_xfx_copy = deepcopy(list_of_int_xfx)
    list_of_int = []
 #   list_of_int_xfx= []
    used_list = []
    for x in range(0, n_max):

        list_of_int.append([])
  #      list_of_int_xfx.append([])

        for y in range(0, len(list_of_int_copy[x])):

            if all_hash[list_of_int_copy[x][y].binary_hash] not in used_list:

                list_of_int[x].append(list_of_int_copy[x][y])
   #             list_of_int_xfx[x].append(list_of_int_xfx_copy[x][y])
                used_list.append(all_hash[list_of_int_copy[x][y].binary_hash])
#                print(list_of_int_copy[x][y], list_of_int_xfx_copy[x][y], all_hash[list_of_int_copy[x][y].binary_hash])
                print(list_of_int_copy[x][y], all_hash[list_of_int_copy[x][y].binary_hash])
            else:
                print(y, 'ten hash juz byl uzyty', list_of_int_copy[x][y], all_hash[list_of_int_copy[x][y].binary_hash], used_list)

    print('UNIQUE')
    for x in range(0, n_max):
        for y in list_of_int[x]:
            print(y, all_hash[y.binary_hash])
        print('')
    print('koniec unique')

    
    

    print('TERAZ STWORZE BACZE')
    memall = 0
    for x in range(0, n_max):
        memmini = 0
        print('parampam', x, len(list_of_int[x]))
        k = 0
        for y in range(0, len(list_of_int[x])):
            memidx, costidx = memcost(list_of_int[x][y])
            # memreal = 888**memidx[0] * 232**memidx[1]* 20**memidx[2]*8/(10**9)
            memreal = 888**memidx[0] * 232**memidx[1]* 20**memidx[2]*8/(10**9)
            memall  += memreal
            memmini += memreal
            print(k, list_of_int[x][y], all_hash[list_of_int[x][y].binary_hash], list_of_int_xfx[x][y], costidx, sum(costidx), memidx, '     ', memreal)
            print(memreal, memmini, memall)
            k+= 1
        print('memmini', memmini)

    print('memall', memall)

    print('basket outer')
    for x in basket_outer:
        print(x)

    for x in range(0, n_max):
        print('parampam', x, len(list_of_int[x]))

    print('-------------------------------------------------')
    print('teraz policze koszty')
    print('-------------------------------------------------')

    for x in range(0, n_max):
        print('level nr ', x)
        cost_old = 0
        z = 0
        costidxm = [0,0,0]
        for y in list_of_int[x]:
            memidx, costidx = memcost(y)
            cost_all = 888**(costidx[0]) * 200**(costidx[1]) * 30**(costidx[2])
#            print(888**(costidx[0]), 200**(costidx[1]), 30**(costidx[2]))
            if cost_all > cost_old:
                cost_old = deepcopy(cost_all)
                idx = deepcopy(z)
                costidxm = deepcopy(costidx)
            print(z, costidx[0], ',', costidx[1], ',', costidx[2], costidx[0]+costidx[1]+costidx[2], cost_all)
            z += 1
        vspace(0)
        cstr = f"c^({costidxm[0]})v^({costidxm[1]})o^({costidxm[2]})"
        print('najdrozsza petla to', idx, cost_old/(10**9), cstr)
        vspace(0)
    print('cost for outer')

    cost_old = 0
    z = 0
    idx = 0
    costidxm = [0,0,0]
    for x in basket_outer:
        memidx, costidx = memcost(x)
        cost_all = 888**(costidx[0]) * 200**(costidx[1]) * 30**(costidx[2])
        if cost_all > cost_old:
            cost_old = deepcopy(cost_all)
            idx = deepcopy(z)
            costidxm = deepcopy(costidx)
        print(z, costidx[0], ',', costidx[1], ',', costidx[2], costidx[0]+costidx[1]+costidx[2], cost_all)
        z += 1
    vspace(0)
    cstr = f"c^({costidxm[0]})v^({costidxm[1]})o^({costidxm[2]})"
    print('najdrozsza petla outer to', idx, cost_old/(10**9), cstr)
    vspace(0)




    basket_all = basket_nointerm + basket_outer
    interm_all = []
    interm_all_xfx = []
    for x in range(0, n_max):
        all_summation_list = []
        for y in range(0, len(list_of_int[x])):
            interm_all.append(list_of_int[x][y])
            interm_all_xfx.append(list_of_int_xfx[x][y])
            print('grrrrrrrrrrr', list_of_int[x][y], list_of_int_xfx[x][y])

        #     y.summation.sort()
        #     if y.summation not in all_summation_list:
        #         all_summation_list.append(y.summation)
        # all_summation_list = sorted(all_summation_list, key=lambda i: len(i))
        # for s in all_summation_list:
        #     for y in list_of_int[x]:
        #         if y.summation == s:
                    # interm_all.append(y)

        # for y in range(0, len(list_of_int[x])):
        #     interm_all.append(list_of_int[x][y])
    print('OUTER-dupa')
    if filel != []:
        if len(basket_all) != 0:
#            filel.write('\\newpage\n')
            filel.write('\\begin{center}\n')
            #zz = "Wzory orbitalne na macierz gestości $\gamma_{{{idx_str}}}$\n".format(idx_str=idx_str)
            zz = "Wzory orbitalne na S $\gamma_{{{idx_str}}}$\n".format(idx_str=idx_str)                          
            filel.write(zz)
            filel.write('\\end{center}\n')
    for x in basket_all:
        print(x)
    print('***********************************************************************')
    #latex_with_cost(filel, 1, basket_all, name=False, f12=True)
    # s = """ wzory na s
    # \\newpage
    # """
    # filel.write(s)

    if s_ampli:
        s2_vv, s2_Cv, s2_vC, s2_CC = divide_s_for_f12_classes(basket_all)
        zz1 = " wzory na $s^{ab}_{ij}$"
        filel.write(zz1)
        just_latex(filel, 1, s2_vv, name=False, f12=True)
        zz2 = " wzory na $s^{Aa}_{ij}$"
        filel.write(zz2)
        just_latex(filel, 1, s2_Cv, name=False, f12=True)
        zz3 = " wzory na $s^{aB}_{ij}$"
        filel.write(zz3)
        just_latex(filel, 1, s2_vC, name=False, f12=True)
        zz4 = " wzory na $s^{AB}_{ij}$"
        filel.write(zz4)
        just_latex(filel, 1, s2_CC, name=False, f12=True)


    just_latex(filel, 1, basket_all, name=False, f12=True, idx_str=idx_str)
    if filel != []:
        if len(interm_all) != 0:
 #           filel.write('\\newpage\n')
            filel.write('\\begin{center}\n')
            filel.write('Konieczne do obliczenia tego rzędu intermediates')
            filel.write('\\end{center}\n')
    print('INTERMEDIATES')
    for x in interm_all:
        print(x)
    for x in interm_all_xfx:
        print(x)
    print('***********************************************************************')
    for ww in interm_all:
        print(ww)

#    latex_with_cost(filel, 1, interm_all, name=True, all_hash=all_hash, f12=True)
    just_latex(filel, 1, interm_all, name=True, all_hash=all_hash, f12=True, interm_all_xfx=interm_all_xfx)
    if filel != []:
        filel.write('\\end{document}\n')
    # return basket_outer, basket_nointerm, list_of_int, n_max, all_hash, mem_dict, interm_fx, xfx_dict
    return basket_outer, basket_nointerm, list_of_int, n_max, all_hash, mem_dict, interm_fx, xfx_dict, list_of_int_xfx

    
def isdisc(res):

    for x in res.coefficient_idx:
        for y in x:
            if y not in res.summation:
                return False

    return True

def identify_unique_interm(list_of_int, list_of_int_xfx, all_hash, version, new_names_dict, disconnected=False):

    unique_dict = {}
    print('dla tego lewelu mamy')    

    new_list = []
    max_summation = 0
    for j in range(0, len(list_of_int)):
        y = list_of_int[j]
        x = deepcopy(y)
        print('xissimo', x)
        if len(x.summation) > max_summation:
            max_summation = len(x.summation)
        x.clear_fixed()
        x.establish_fixed()
        x.summation += fixed        
        print('jaaa', x, fixed, list_of_int_xfx[j])
        d = simplify(arithmetic_string(x), drag_fixed=deepcopy(list_of_int_xfx[j]))#fixed))

        go = True
        if disconnected:
            if len(d[0].df) > 0:
                go  = True
            else:
                go = False

        if go:
            for i in d[0].df[1]:
                if i in d[0].summation:
                    d[0].summation.remove(i)
            print('ddissimo', d[0], d[0].df)
            new_list.append([d[0], d[0].df[1]])
        else:
            new_list.append([d[0], []])


    used_interm = []
    unique_list = []
    for i in range(0, len(new_list)):
        if i not in used_interm:
            print('dodaje do unique i used',  all_hash[list_of_int[i].binary_hash])
            unique_list.append(i)
            used_interm.append(i)
            unique_dict[i] = deepcopy(i)
            for j in range(i+1, len(new_list)):
                print('sprawdzam',  all_hash[list_of_int[j].binary_hash])
                if j not in used_interm:
                    if new_list[i][0].summation == new_list[j][0].summation:
                        if new_list[i][1] == new_list[j][1]:
                            if new_list[i][0] == new_list[j][0]:
                                print('te dwa sa takie same', i, j, all_hash[list_of_int[i].binary_hash], all_hash[list_of_int[j].binary_hash],\
                                      list_of_int[i], list_of_int[j])
                                # print('te dwa sa takie same',  all_hash[list_of_int[i].binary_hash], all_hash[list_of_int[j].binary_hash], new_names_dict[all_hash[list_of_int[i].binary_hash]], \
                                #       new_names_dict[all_hash[list_of_int[J].binary_hash]])
                                # te dwa intermediate sa takie same
                                used_interm.append(j)
                                unique_dict[j] = deepcopy(i)

    new_list_of_int = []
    new_list_of_fx = []
    print('unique version', version)
    for key in unique_dict:
        print(key, new_list[key][0], new_list[key][1], list_of_int[key], unique_dict[key])
        if version == 2:
            all_hash[new_list[key][0].binary_hash] = all_hash[list_of_int[key].binary_hash]
            #list_of_int[key] = new_list[key][0]
            new_list_of_int.append(new_list[key][0])
            new_list_of_fx.append(new_list[key][1])


                        
    if version == 2:
        return new_list_of_int, new_list_of_fx
    else:
        return unique_dict


def extract_int_no(name):
    ns = '1234567890'
    no =""
    for s in name:
        if s in ns:
                no += s

    return no

def intermediates_to_fortran(theory, prefix, basket_outer, basket_nointerm, list_of_int, n_max, all_hash, mem_dict, interm_fx, xfx_dict):

    '''! this function translates intermediates to Fortran language'''
    """! 
    param theory - 'ccsd' or 'cc3'\\
    """
    
    list_of_names = []
    list_of_int_mem = []
    list_of_int_disk = []
    list_of_int_mem_names = []
    list_of_int_disk_names = []
    for n in range(0, n_max):
        list_of_names.append([])
        list_of_int_mem.append([])
        list_of_int_disk.append([])
        list_of_int_mem_names.append([])
        list_of_int_disk_names.append([])



    print(len(list_of_int), n_max)
    all_interm_dict = {}
    mem_mem= []
    mem_disk = []
    for x in range(0, n_max):
        mem_mem.append([])
        mem_disk.append([])

    print('teraz bedzie podzial pamieci')
    mem_ram_total = 0
    mem_disk_total = 0

    
    for x in range(0, n_max):
        used_dict_temp = []

        print('len(list_of_int[n])', x, len(list_of_int[x]))
        mem_memt = 0
        mem_diskt = 0
        minidict = {}
        for l in range(0, len(list_of_int[x])):
            minidict = {}
            minidict['lv'] = x # <- level of intermediate. Level 0 needs to be computed bf level 1
            k = list_of_int[x][l]
            name_int = all_hash[k.binary_hash]
            mem_int = mem_dict[k.binary_hash]
            
            minidict['name'] = name_int
            minidict['if'] = interm_fx[k.binary_hash]
            minidict['ugg'] = k
            minidict['no'] = extract_int_no(name_int)
            minidict['filename'] = 'file_'+name_int
            minidict['xfx'] = xfx_dict[k.binary_hash]
            
            if k.binary_hash not in used_dict_temp:
                print('name_int', name_int, k)
                if (mem_int) < MEM_LITTLE_THRESH:  #<-- these intermediates are going to be stored in memory
                    minidict['mem'] = 'ram'
                    list_of_int_mem[x].append(k)
                    used_dict_temp.append(k.binary_hash)
                    list_of_int_mem_names[x].append(name_int)
                    mem_memt += mem_int
                    mem_ram_total += mem_int
           #         print('RAM', mem_int, MEM_LITTLE_THRESH)
                else:
                    minidict['mem'] = 'disk'
                    used_dict_temp.append(k.binary_hash)
                    list_of_int_disk[x].append(k)
                    list_of_int_disk_names[x].append(name_int)
                    mem_diskt += mem_int
                    mem_disk_total += mem_int
            #        print('DISK', mem_int, MEM_LITTLE_THRESH)
            all_interm_dict[name_int] = deepcopy(minidict)
            mem_mem[x] = mem_memt
            mem_disk[x] = mem_diskt
            if mem_disk[x] < MEM_DISK_THRESH:
                for l in range(0, len(list_of_int_disk[x])):
                    k = list_of_int_disk[x][l]
                    name_int = all_hash[k.binary_hash]
                    all_interm_dict[name_int]['mem'] = 'ram'
                    k = list_of_int_disk[x][l]
                    list_of_int_mem[x].append(k)
                    list_of_int_mem_names[x].append(list_of_int_disk_names[x][l])
                mem_mem[x] += mem_disk[x]
                list_of_int_disk[x] = []
                list_of_int_disk_names[x] = []
                mem_disk[x] = 0

    print('calkowita pemiec potrzebna RAM')
    print(mem_ram_total)
    print('calkowita pemiec potrzebna DISK')
    print(mem_disk_total)

    # print('i intermedaites wszystkie z pamieciami')
    # k = 0
    # for x in range(0, n_max):
    #     kk = 0        
    #     for l in range(0, len(list_of_int[x])):
    #         # print(k, kk, list_of_int_mem_names[x][l], list_of_int_mem[x][l], list_of_int[x][l])
    #         print(list_of_int[x][l])
    #         kk += 1
    #     k += 1
    # print('')
    # sys.exit(0)

    big_batch_list = []
    for n in range(0, n_max):
        big_batch_list.append([])
        blen = len(list_of_int[n])

        nstep = blen // MAX_TERMS_PER_BATCH
        rest = blen%MAX_TERMS_PER_BATCH
        step = MAX_TERMS_PER_BATCH

        for j in range(0, nstep):
            big_batch_list[n].append(list_of_int[n][j*step:(j+1)*step])

        if rest != 0:
            big_batch_list[n].append(list_of_int[n][nstep*step:])


    print('INTERMEDIATES LEVELS')
    for n in range(0, n_max):
        for y in big_batch_list[n]:
            for x in y:
                memidx, costidx = memcost(x)
                print(n, x, 'KOSZT z Duzego batcha', costidx)
            print('')

    print('BASKET NOINTERM')
    for x in basket_nointerm:
        print(x)
        basket_outer.append(x)

    print('BASKET OUTER')
    for x in basket_outer:
        print(x)

    # latex_S_amplitudes(n_max, list_of_int, all_hash, basket_outer)
    # print('la')
    # print('BASKET NOINTERM')
    # for x in basket_nointerm:
    #     print(x)
    #     basket_outer.append(x)


    batch_outer = []
    blen = len(basket_outer)
    print('BLEN', blen)
    nstep = blen // MAX_TERMS_PER_BATCH_OUTER
    rest = blen%MAX_TERMS_PER_BATCH_OUTER
    step = MAX_TERMS_PER_BATCH_OUTER
    for j in range(0, nstep):
        batch_outer.append(basket_outer[j*step:(j+1)*step])
    if rest != 0:
        batch_outer.append(basket_outer[nstep*step:])

    for x in batch_outer:
        print(len(x))
        for y in x:
            print(y)

    unit_decl0 = []
    print('to bedzie unit decl')
    int_app = "{prefix}_ccsd_f12".format(prefix=prefix)
    """! dodajemy deklaracje  unitow dla 
    wszystkich intermediates, ktore beda odczytywane z dysku
    """
    for n in range(0, n_max):
        unit_decl_temp = ""
        for l in range(0, len(list_of_int_disk[n])):
            k = list_of_int_disk[n][l]
            name_int = all_hash[k.binary_hash]
            if (all_interm_dict[name_int]['mem'] == 'disk'):
                unit_decl_temp += "integer :: u{name_int}_{int_app}\n".format(name_int=name_int, int_app=int_app)
        unit_decl0.append(unit_decl_temp)
        # print(unit_decl0)
        # print( '')

    for n in range(0, n_max):
        print('uunit', n )
        print(unit_decl0[n])
        print('')
        
    s_decl = ""
    for n in range(0, n_max):
        declname = prefix+"_lev"+str(n)
        # if theory == 'ccsd':
        s_decl += """use decl_{declname}_interm_ccsd_f12                                                                                             
                          """.format(declname=declname)
        # elif theory == 'cc3':
        #     s_decl += """use decl_{declname}_interm_cc3_f12}
        #                   """.format(declname=declname)


    print('n_max', n_max)
    dname = ''
    for n in range(0, n_max):
        new_decl = True
        declname = prefix+"_lev"+str(n)
        for x in range(0, len(big_batch_list[n])):
            # if n == 0:
            unit_decl = unit_decl0[n]
            # else:
            #     unit_decl = ""
            partname = prefix+"lev"+str(n)+"_batch"+str(x)
            new_file = True
            last_file = True
            last_decl = False
            if x ==  len(big_batch_list[n]) - 1:
                last_decl = True
            dname = function_template_batch_intermediates_f12(big_batch_list[n][x], x, prefix, partname, declname, all_hash, all_interm_dict,\
                                                              new_file, last_file, \
                                                              new_decl, last_decl, s_decl, unit_decl)
            new_decl = False

    for x in range(0, len(batch_outer)):
        
        partname = "batch"+str(x)
        new_file = True
        last_file = True
        print('sraka', s_decl)
        # block_oo, block_ov, block_vo, block_vv = fill_dm_blocks_wm(batch_outer[x])

        batch_outer_ars = arithmetic_string()
        for y in batch_outer[x]:
            batch_outer_ars += arithmetic_string(y)
            print('bb', y)
        function_template_batch_f12(batch_outer_ars, prefix, x, \
                                    partname, all_hash, theory, new_file, last_file, dname, s_decl, all_interm_dict)


    print('KONCZE')
    sys.exit(0)







    

