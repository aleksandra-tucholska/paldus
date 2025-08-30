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
from fortran_code_factor import *
from templates import  *
from factor import *
from joblib import Parallel, delayed
from paldus_f12 import factorize_driver
import sys
import time
import pickle

MAX_TERMS_PER_BATCH = 150
MAX_TERMS_PER_BATCH_OUTER = 800
MEM_DISK_THRESH = 5
MEM_LITTLE_THRESH = 0.015

# Wm =  <e(-S)e(T*) X e(-T*) e(S) P(e(-S)e(T*) \mu_n e(-T*) e(S))| P(e(S*)\mu_l e(-S*))>                                                         
        #--------Wm1----------#  #-------------Wm2-------------#  #-------Wm3-------#

# Gm =  <e(S*)e(-T) E_pq e(T) e(-S*) P(e(S*)e(-T) E_rs e(T) e(-S*))>                                                         
        #--------Gm1----------#  #-------------Gm2-------------#  

def compute_cost_with_intermediates(Wm_intermediates_list):
    cost_virt = []
    print(len(Wm_intermediates_list))
    for x in Wm_intermediates_list:

        s_idx_int1_v = []
        s_idx_int1_o = []
        s_idx_int2_v = []
        s_idx_int2_o = []
        s_idx_nonint_v = []
        s_idx_nonint_o = []
        s_idx_org_v = []
        s_idx_org_o = []
        fx_idx_int1_v = []
        fx_idx_int1_o = []
        fx_idx_int2_v = []
        fx_idx_int2_o = []
        fx_idx_nonint_v = []
        fx_idx_nonint_o = []
        fx_idx_org_v = []
        fx_idx_org_o = []

        for y in x['interm1'].summation:
            if y in virtual:
                s_idx_int1_v.append(y)
            else:
                s_idx_int1_o.append(y)

        for y in x['interm2'].summation:
            if x in virtual:
                s_idx_int2_v.append(y)
            else:
                s_idx_int2_o.append(y)
        for y in x['noninterm'].summation:
            if y in virtual:
                s_idx_nonint_v.append(y)
            else:
                s_idx_nonint_o.append(y)
        for y in x['original'].summation:
            if y in virtual:
                s_idx_org_v.append(y)
            else:
                s_idx_org_o.append(y)

        for y in x['interm1'].coefficient_idx:
            for z in y:
                if z not in x['interm1'].summation:
                    if z in virtual  and z not in fx_idx_int1_v:
                        fx_idx_int1_v.append(z)
                    elif z in occupied  and z not in fx_idx_int1_o:
                        fx_idx_int2_o.append(z)
        for y in x['interm2'].coefficient_idx:
            for z in y:
                if z not in x['interm2'].summation:
                    if z in virtual and z not in fx_idx_int2_v:
                        fx_idx_int2_v.append(z)
                    elif z in occupied and z not in fx_idx_int2_o:
                        fx_idx_int2_o.append(z)
        for y in x['noninterm'].coefficient_idx:
            for z in y:
                if z not in x['noninterm'].summation:
                    if z in virtual and z not in fx_idx_nonint_v:
                        fx_idx_nonint_v.append(z)
                    elif z in occupied and z not in fx_idx_nonint_o:
                        fx_idx_nonint_o.append(z)
        for y in x['original'].coefficient_idx:
            for z in y:
                if z not in x['original'].summation:
                    if z in virtual and z not in fx_idx_org_v:
                        fx_idx_org_v.append(z)
                    elif z in occupied and z not in fx_idx_org_o:
                        fx_idx_org_o.append(z)
        s = """"original cost: v^{a}o^{b}, nonint cost: v^{c}o^{d}, interm1 cost: v^{e}o^{f}, interm2 cost: v^{g}o^{h} """.format(a=len(s_idx_org_v + fx_idx_org_v), b = len(s_idx_org_o + fx_idx_org_o), c = len(s_idx_nonint_v + fx_idx_nonint_v), d = len(s_idx_nonint_o + fx_idx_nonint_o), e = len(s_idx_int1_v + fx_idx_int1_v), f = len(s_idx_int1_o + fx_idx_int1_o), g = len(s_idx_int2_v + fx_idx_int2_v), h = len(s_idx_int2_o + fx_idx_int2_o))
        print('ssss')
        print(s)
        cost_virt.append(len(s_idx_int1_v + fx_idx_int1_v))
        cost_virt.append(len(s_idx_int2_v + fx_idx_int2_v))

    if (cost_virt != []):
        print('maxcost_virtual', max(cost_virt))
    else:
        print()


def mem_of_all_intermediates(Wm_with_int):

    print('COMPUTING MEMORY COST OF ALL INTERMEDIATES')
    print('------------------------------------------')
    print('Number of intermediates', len(Wm_with_int))
    
    list_all = {}

    for x in range(0, len(Wm_with_int)):
        int1 = Wm_with_int[x]
        addlist = []
        if int1['interm'].num_factor != 0:
            print('int1', int1['interm'])
            addlist = int_mem(int1['interm'])
            print(addlist)
        keys = list(list_all.keys())
        print('KEYS', keys)
        if addlist != []:
            stradd = str(addlist)
            if stradd not in keys:
                print('nie ma')
                list_all[stradd] = 1
            else:
                print('jest')
                list_all[stradd] += 1

    print(list_all)
    print('Estimated memory for 10 occupied and 200 virtual')

    cost_all = 0
    for key in list_all:
        oo = int(key[4])
        vv = int(key[1])
        x = list_all[key]
        cost = (10**oo * 200**vv * 64 / 10**9) * x
        s = """ there is {x} terms of form v^{v}o^{o}   ----   {cost}""".format(x = list_all[key], v = key[1], o = key[4], cost=cost)
        print(s)
        cost_all += cost
    print('COST MEM ALL', cost_all)
        
    print('------------------------------------------')
    
def int_mem(interm):

    addedo = 0
    addedv = 0
    addedol = []
    addedvl = []

    for j in interm.coefficient_idx:
        for i in j:
            if i not in interm.summation:
                if i in virtual:
                    if i not in addedvl:
                        addedv += 1
                        addedvl.append(i)
                elif i in occupied:
                    if i not in addedol:
                        addedo += 1
                        addedol.append(i)

    return [addedv, addedo]

def cost_after_intermediates(xx):
    s3 = []
    for x in range(0, len(xx.coefficient) - 1):
        for y in range(x+1, len(xx.coefficient)):
            s = deepcopy(xx.coefficient_idx[x])
            for yy in xx.coefficient_idx[y]:
                if yy not in s:
                    s.append(yy)
            t = []
            for i in range(0, len(xx.coefficient)):
                if i != x and i != y:
                    for j in xx.coefficient_idx[i]:
                        if j not in t:
                            t.append(j)
            s2 = []
            for z in s:
                if z in xx.summation:
                    if z not in t:
                        s2.append(z)
            if len(s2) > len(s3):
                s3 = s2
            elif len(s2) == len(s3):
                t1 = []
                for z in s2:
                    if z in virtual:
                        t1.append(z)
                t2 = []
                for z in s3:
                    if z in virtual:
                        t2.append(z)
                if len(t1) > len(t2):
                    s3 = s2
    svi = []
    svo = []
    for z in s3:
        if z in virtual:
            svi.append(z)
        elif z in occupied:
            svo.append(z)
            
    return svi, svo

def temper():

    S1capol  = arithmetic_string(s1c).scale(0.5)
    S2ca  = arithmetic_string(s2c)
    YT2T2 = evaluate(nobs, t2, t2)
    r = []
    r.append(S1capol * S2ca * YT2T2)

    for i in range(0, len(r)):
        rint = r[i].integrate()
        rint.exec_delta()
        rsimp = simplify(rint)

        block_oo_diag, block_vv_diag, rsimp_without_diag = fill_dm_blocks_diag(rsimp)

        block_oo, block_ov, block_vo, block_vv = fill_dm_blocks(rsimp_without_diag)

        for x in block_ov:
            print(x)
            # s = []
            # for y in x.summation:
            #     s.append(y)
            # for y in x.coefficient_idx:
            #     for z in y:
            #         if z not in s:
            #             s.append(z)
            # sv = []
            # so = []
            # for y in s:
            #     if y in occupied:
            #         so.append(y)
            #     elif y in virtual:
            #         sv.append(y)
            # if len(x.coefficient) > 2:
            #     svi, svo = cost_after_intermediates(deepcopy(x))
            #     print(x)
            # else:
            #     print(x)

            
def generate_cumulant_ground_f12_driver(pick, method, only_quadratic=False):

    f = open("Cumulant_f12.tex", 'w')
    s_preamble = tex_preamble
    s_beginning = """
    \\begin{document}"""
    f.write(s_preamble)
    f.write(s_beginning)

    for order in range(2, 6):
        print('teraz bedzie rzad', order)
        s2 = "\n"
        if order > 0:
            s2 += "\\newpage \n"
        s2 += f"Kumulanta w rzędzie {order}\\\\"
        s2 += "\n"
        f.write(s2)
        generate_cumulant_ground(pick, order, method, f, withCABS=True, only_quadratic=only_quadratic)
        f.write("\n")
        
        f.write("\\hrule")
        f.write("\n")

       
    f.write("\end{document}")
    f.close()

def generate_cumulant_ground_f12(pick, method, mbpt, f):

    # < 0 | e(-S*)e(-T) X e(T) e(-S*)P(e(-S*)e(-T) X e(T) e(-S*)) | 0>
    # < 0 | e(-S*)e(-T) X e(T) e(-S*) | 0>
    # < 0 | e(-S*)e(-T) X e(T) e(-S*) | 0>


    # Gm2_ground  P(e(S*)e(-T) X e(T) e(-S*))
    # Gm1_ground  e(S*)e(-T) X e(T) e(-S*)

    
    open_flags = ""

    if pick == "dump":
        open_flags = "wb"
    elif pick == "load":
        open_flags = "rb"

    print('GENERATING Gm2')
    Gm2t = generate_Wm_ground_f12(mbpt, 'ccsd', Gm1 = False, Gm2=True)
    print('GENERATING Gm1')
    Gm1t = generate_Wm_ground_f12(mbpt, 'ccsd', Gm1=True, Gm2 = False)

    print('gm2')
    for x in Gm2t:
        print(x)
    print('')
    print('gm1')
    for x in Gm1t:
        print(x)

    Gm_middle = generate_Gm(Gm1t, Gm2t, mbpt, False, "ccsd")
    print('dla tego rzedu mam', mbpt, len(Gm_middle))

    print('gmmuahaha')
    for x in Gm_middle:
        print(x)

    print('gmmuahaha-latex')

    latex_Gmiddle(Gm_middle, mbpt, f)

    #sys.exit(0)

    # Gm1 = []
    # Gm2 = []

    # for x in Gm2t:
    #     print(x['mbpt'], x)
    #     if x['mbpt'] <= mbpt:
    #         Gm2.append(x)
    # for x in Gm1t:
    #     print(x['mbpt'], x)
    #     if x['mbpt'] <=mbpt:
    #         Gm1.append(x)

    # print('Gm2')
    # for x in Gm2:
    #     print(x)
    # print('Gm1')
    # for x in Gm1:
    #     print(x)
                    
                   
    # G2_middle, G2_commutators = commutators_density_ground_f12(obsc, Gm2, 'T_list', 'S_list')
    # G1_middle, G1_commutators = commutators_density_ground_f12(obscy, Gm1, 'T_list', 'S_list')

    # print(len(G2_middle))
    # print(len(G2_commutators))
    # latex_Gm(G1_middle, G2_middle,  mbpt, f)


        
    # print('czesc Gm2')
    # G2_f12 = arithmetic_string()
    # for x in G2_commutators:
    #     for y in x:
    #         print(y)
    #         G2_f12.append(y)
    #     print('')

    # print('czesc Gm1')                
    # G1_f12 = arithmetic_string()
    # for x in G1_commutators:
    #     for y in x:
    #         print(y)
    #         G1_f12.append(y)
    #     print('')


    # print('gl_middle')
    # for x in Gl_middle:
    #     print(x)


    # print('gr_middle')
    # for x in Gr_middle:
    #     print(x)


    if pick == "dump":
        # rint = integrate_Gm(G1_commutators, G2_commutators, G1_middle, G2_middle, method, mbpt)

        # # rint = W_f12.integrate()
        
        # rint.exec_delta()

        # print('result przed simplify ground ground')
        # kk = 1
        # for yy in rint:
        #     print(kk, yy)
        #     kk += 1
        # print('')
        # pick_cumulant_ccsd_f12 = open('./pickle/cumulant_ccsd_f12_{mbpt}.pkl'.format(mbpt=mbpt), open_flags)
        
        # pickle.dump(rint, pick_cumulant_ccsd_f12)

        a=1
    elif pick == "load":
        
        load_cumulant_ccsd_f12 = open('./pickle/cumulant_ccsd_f12_{mbpt}.pkl'.format(mbpt=mbpt), open_flags)
        rint = pickle.load(load_cumulant_ccsd_f12)


        rsimp = simplify(rint)
        print('result pooo simplify ground ground')
        kk = 1
        for yy in rsimp:
            print(kk, yy)
            kk += 1
        print('')

        print('teraz bede zmieniac x * y na x z czterema indeksami')

        for i in range(0, len(rsimp)):
            x_j = 100
            y_j = 100
            print('przed', rsimp[i])
            coef_list = []
            for j in range(0, len(rsimp[i].coefficient)):
                if rsimp[i].coefficient[j] == OBSERVABLE_X:
                    x_j = deepcopy(j)
                    for k in rsimp[i].coefficient_idx[j]:
                        coef_list.append(k)
                    if x_j > y_j:
                        x_j = x_j - 1
                if rsimp[i].coefficient[j] == OBSERVABLE_Y:
                    y_j = deepcopy(j)
                    for k in rsimp[i].coefficient_idx[j]:
                        coef_list.append(k)
                    if y_j > x_j:
                        y_j = y_j - 1
            del rsimp[i].coefficient[x_j]
            del rsimp[i].coefficient[y_j]
            del rsimp[i].coefficient_idx[x_j]
            del rsimp[i].coefficient_idx[y_j]

            rsimp[i].coefficient.append(OBSERVABLE_X)
            rsimp[i].coefficient_idx.append(coef_list)
            print('po', rsimp[i])
           
        print('')
                    


        print("I'm removing all instances of ff{akbl}, where a and b are both wirtual indices")
        print("According to definition of F in Shiozaki article, these instances are equal to zero")
        for x in rsimp:
            for j in range(0, len(x.coefficient)):
                y = x.coefficient[j]
                if y == F12_TWOEL or y ==F12_TWOEL_COMBO:
                    n_of_virt = 0
                    for z in x.coefficient_idx[j]:
                        if z in virtual:
                            n_of_virt +=1
                    if n_of_virt >=2:
                        x.num_factor =0

        rsimp.exec_delta()
        rsimp.cleanup()
        res_x = identify_interm_X_f12(rsimp)
        print('po find X')
        k = 0
        for p in res_x:
            print(k, p)
            k +=  1
        print('')
        res_x = identify_interm_Z_f12(res_x, OBSERVABLE_X)
        print('po find Z')
        k = 0
        for p in res_x:
            print(k, p)
            k +=  1
        print('')

        res = arithmetic_string()
        for x in res_x:
            res = res + arithmetic_string(x)
            print(x)
        print('')
        print('przed cabs transform')
        cabs_res = res.cabstransform()
        print('po cabs transform s1')
        z = 0
        for x in cabs_res:        
            print(z, x)
            z += 1
        res_ft = identify_interm_Ft_f12(cabs_res)
        print('to po find Ft')

        res_ost = arithmetic_string()
        for x in res_ft:
            res_ost.append(x)

        rsimp_ost = simplify(res_ost)

        print('i to jest ostateczny wyniczor po symplifikacji')    
        k = 0
        for x in rsimp_ost:
            print(k, "   ", x)
            k += 1

        print('koszty bez intermediates cumulant')
        for x in rsimp_ost:
            memidx, costidx = memcost(x)
            print(costidx[0], ',', costidx[1], ',', costidx[2], costidx[0] + costidx[1]+costidx[2])

        f.write('\\newpage')
        print('wynik ostateczny latex')
        latex_with_cost(f,1, rsimp_ost, name=False)
        print('')

        prefix = 'cumulant_{mbpt}'.format(mbpt=mbpt)
        excl = []

        basket_outer, basket_nointerm, list_of_int, n_max, all_hash, mem_dict, \
            interm_fx, xfx_dict, list_of_int_xfx  = factorize_driver(rsimp_ost, excl, COST=6, is_x_excluded=True, filel=f)
        
        # intermediates_to_fortran('ccsd', prefix, basket_outer, basket_nointerm, list_of_int, n_max, \
        #                          all_hash, mem_dict, interm_fx, xfx_dict)

        # basket_outer_s2, basket_nointerm_s2, list_of_int_s2, n_max_s2, all_hash_s2, mem_dict_s2, \
        #     interm_fx_s2, xfx_dict_s2  = factorize_driver(rsimp_ost_s2, [], COST=6)


def generate_cumulant_ground_driver(pick, method):

    
    f = open(f"Cumulant_{method}.tex", 'w')
    s_preamble = tex_preamble
    s_beginning = """
    \\begin{document}"""
    f.write(s_preamble)
    f.write(s_beginning)

    for order in range(1, 2):
        print('teraz bedzie rzad', order)
        s2 = "\n"
        if order > 0:
            s2 += "\\newpage \n"
        s2 += f"Kumulanta w rzędzie {order}\\\\"
        s2 += "\n"
        f.write(s2)
        if method == 'ccd':
            generate_cumulant_ground(pick, order, 'ccd', f)
        elif method == 'cc3':
            generate_cumulant_ground(pick, order, 'cc3', f)            
        elif method == 'f12':
            generate_cumulant_ground(pick, order, method, f, withCABS=True, only_quadratic=True)
            
        f.write("\n")
        
        f.write("\\hrule")
        f.write("\n")

        
    f.write("\end{document}")
    f.close()

            


        
def generate_cumulant_ground(pick, mbpt, method, f, withCABS=False, only_quadratic=False):

    # < 0 | e(-S*)e(-T) X e(T) e(-S*)P(e(-S*)e(-T) X e(T) e(-S*)) | 0>
    # < 0 | e(-S*)e(-T) X e(T) e(-S*) | 0>
    # < 0 | e(-S*)e(-T) X e(T) e(-S*) | 0>


    # Gm2_ground  P(e(S*)e(-T) X e(T) e(-S*))
    # Gm1_ground  e(S*)e(-T) X e(T) e(-S*)

    
    open_flags = ""

    if pick == "dump":
        open_flags = "wb"
    elif pick == "load":
        open_flags = "rb"

    if method == 'ccd':
        print('GENERATING Gm2')
        Gm2t = generate_Wm_ground_ccd(mbpt, 'ccsd', Gm1 = False, Gm2=True)
        print('GENERATING Gm1')
        Gm1t = generate_Wm_ground_ccd(mbpt, 'ccsd', Gm1=True, Gm2 = False)
    elif method == 'cc3':
        print('GENERATING Gm2')
        Gm2t = generate_Wm_ground_ccd(mbpt, 'cc3', Gm1 = False, Gm2=True)
        print('GENERATING Gm1')
        Gm1t = generate_Wm_ground_ccd(mbpt, 'cc3', Gm1=True, Gm2 = False)
    elif method == 'f12':
        print('GENERATING Gm2')
        Gm2t = generate_Wm_ground_f12(mbpt, 'ccsd', Gm1 = False, Gm2=True)
        print('GENERATING Gm1')
        Gm1t = generate_Wm_ground_f12(mbpt, 'ccsd', Gm1=True, Gm2 = False)

    print('------------------------------------------------------')
    print('gm2ratatata', method)
    for x in Gm2t:
        print(x)
    print('')
    print('gm1')
    for x in Gm1t:
        print(x)
    print('------------------------------------------------------')
    
    Gm_middle = generate_Gm(Gm1t, Gm2t, mbpt, False, method)
    print('dla rzedu ', mbpt, 'mam', len(Gm_middle), 'wyrazow')

    print('')
    print('a tak wyglada Gm_middle')
    print('')
    Gm_middle_only_quadratic = []
    for x in Gm_middle:
        s_order = len(x['Gm1_S_list']) + len(x['Gm1_T_list']) + len(x['Gm2_S_list']) + len(x['Gm2_T_list'])
        if s_order <= 2:
            Gm_middle_only_quadratic.append(x)
        print(x)

    if only_quadratic:
        print('TYLKO WYRAZY KWADRATOWE')
        Gm_middle =  Gm_middle_only_quadratic
    for x in Gm_middle:
        print(x)

    print('len(Gm_middle)', len(Gm_middle))

#    Gm_middle = [Gm_middle[2]]
    for x in Gm_middle:
        print(x)

    latex_Gmiddle(Gm_middle, mbpt, f, withCABS)

    Commutators_dict = {}
    # in order to create a two particle matrix, we need to find all combinations of pqrs
    # ie abcd, abci, abic, abji....ijkl
    # For each of the four indices _ _ _ _ we need to pick a name from occupied or virtual part
    if method == 'ccd' or method == 'ccsd' or method == 'cc3':
        choosing_pool = [virtual, occupied]
    elif method == 'f12':
#        choosing_pool = [CABS, virtual, occupied]
         choosing_pool = [completev, occupied]

    banned = []
    Oneel_exch_dict = {}
    Oneel_delta_dict = {}
    for x1 in choosing_pool:
        a1 = free_idx(x1, [])
        banned.append(a1)
        for x2 in choosing_pool:
            a2 = free_idx(x2, [a1])
            banned.append(a2)
            for x3 in choosing_pool:
                a3 = free_idx(x3, [a1,a2])
                banned.append(a3)                    
                for x4 in choosing_pool:
                    a4 = free_idx(x4, [a1,a2,a3])
                    banned.append(a4)
                    op1 = operat1([a1, a2], 's')
                    op2 = operat1([a3, a4], 's')
                    op3 = operat1([a1, a4], 's')
                    op4 = operat1([a3, a2], 's')
                    # if op1.operator_idx == [['a', 'i']] and op2.operator_idx == [['b', 'j']] or\
                    #    op1.operator_idx == [['a', 'i']] and op2.operator_idx == [['j', 'b']] :
                    # if op1.operator_idx == [['i', 'a']] and op2.operator_idx == [['j', 'b']] or\
                    #    op1.operator_idx == [['i', 'a']] and op2.operator_idx == [['b', 'j']] :
                    go = True
                    if go:
                        print('')
                        print('glk', [a1, a2], [a3, a4])
                        idx_str = a1+a2+a3+a4
                        print('TERAZ LICZE KOMUTATORY', idx_str)
                        print('')
                        Gm_commutators = commutators_cumulant(op1, op2, method, Gm_middle)
                        vspace(0)
                        print('KOMUTATORY POLICZONE dla ', idx_str)
                        print('length of the result is', len(Gm_commutators))
                        vspace(0)
                        # for ff in Gm_commutators:
                        #     print('lfff', ff, len(ff))
                        #     if idx_str == 'abcd':
                        #         for x in ff:
                        #             print(x)
                        #     # print('lfff', ff, len(ff))
                        # print('')

                        # Evaluate 1/2 <0|E_ps|0><0|E_qr|0> and -delta_{qr}<0|Eps|0>
                        # rsimp_oneel_exch = 1/2 <0|E_ps|0><0|E_qr|0>
                        # rsimp_oneel_delta = -delta_{qr}<0|Eps|0>
                        rsimp_oneel_exch = arithmetic_string()
                        rsimp_oneel_delta = arithmetic_string()

                        fdens = open(f"dens_do_cumul_{method}.tex", 'w')
                        vspace(0)
                        print('a to z macierzy gestosci')
                        vspace(0)
                        for mbpt1 in range(0, mbpt+1):
                            
                            print('przed mbpt1 dla', idx_str,'to', mbpt1, mbpt-mbpt1, a1, a2, a3, a4)
                            # print('op3', op3, op4)
                            if method == 'ccd' or method == 'cc3':
                                rsimp_eps, W_middle1 = \
                                    generate_density_Wm_ground_mbpt(pick, method, mbpt1, \
                                                                    fdens, onlyrsimp = True, operator=op3, cumulative = False)
                                rsimp_eqr, W_middle2 = \
                                    generate_density_Wm_ground_mbpt(pick, method, mbpt-mbpt1, \
                                                                    fdens, onlyrsimp = True, operator=op4, cumulative = False)
                                
                            elif method == 'f12':
                                rsimp_eps, W_middle1, ml = \
                                    generate_density_matrix_ground_f12(pick, method, mbpt1, fdens, onlyrsimp = True, operator=op3)
                                rsimp_eqr, W_middle2, ml = \
                                    generate_density_matrix_ground_f12(pick, method, mbpt-mbpt1, fdens, onlyrsimp = True, operator=op4)

                            print('finitissimo dla minirzedu', mbpt1, mbpt-mbpt1, mbpt)

                            if len(rsimp_eps)==0:
                                print('rsimp_eps zero')
                            else:
                                print('rsimp_eps NOT zero')
                                for x in rsimp_eps:
                                    print(x)
                                print('')
                                
                            if len(rsimp_eqr)==0:
                                print('rsimp_eqr zero')
                            else:
                                print('rsimp_eqr NOT zero')
                                for x in rsimp_eqr:
                                    print(x)
                                print('')



                            kkk = 1
                            for elem1 in rsimp_eps:
                                for elem2 in rsimp_eqr:
                                    print(kkk, 'srelele1', elem1, elem2)
                                    disambiguate(elem1, elem2)
                                    print(kkk, 'srelele2', elem1, elem2)
                                    kkk += 1
                                    r0 = elem2.fromleft(elem1)
                                    r0.num_factor = r0.num_factor * 0.5                                
                                    rsimp_oneel_exch.append(r0)

                            # if len(rsimp_oneel_exch) > 0:
                            #     str = '\n {idx_str} iloczyn jednoelektronowych macierzy gestosci  \n'.format(idx_str=idx_str)
                            #     f.write(str)
                            #     latex_WmWm(W_middle1, W_middle2, [a1, a4], [a3, a2], mbpt, f)


                        # print('mbpt', mbpt1, a1, a2, a3, a4)
                        print('op3-przed-delta', op3, op4, idx_str)
                        if method == 'ccd' or method == 'cc3':
                            print('generuje delta delta dla rzedu', mbpt)
                            rsimp_eps, W_middle = generate_density_Wm_ground_mbpt(pick, method, mbpt, fdens, onlyrsimp = True, operator=op3, cumulative = False)
                            print('i wygenerowane')
                            for x in rsimp_eps:
                                print(x)
                            print('koniec printa')
                            
                        elif method == 'f12':
                            rsimp_eps, W_middle, ml = generate_density_matrix_ground_f12(pick, method, mbpt, fdens, onlyrsimp = True, operator=op3)
                            
                        elem2 = ugg()
                        elem2.delta = [[a2, a3]]
                        # rsimp_eqr = arithmetic_string(rsimp_eqr)
                        # print('rssss',idx_str, rsimp_eps)
                        print('wm', W_middle)
                        kkk = 1
                        for elem1 in rsimp_eps:
                            print(kkk, 'frelele1', elem1, elem2)
                            disambiguate(elem1, elem2)
                            print(kkk , 'frelele2', elem1, elem2)
                            kkk += 1
                            r0 = elem2.fromleft(elem1)
                            r0.num_factor = r0.num_factor * (-1.0)
                            # print('r00001', mbpt, idx_str,  r0)
                            rsimp_oneel_delta.append(r0)

                        # for elem in rsimp_oneel_delta:
                        #     print('r0000', mbpt1, mbpt-mbpt1, idx_str, elem)

                        # if len(rsimp_oneel_delta) > 0:
                        #     str = '\n {idx_str} iloczyn z delta \n'.format(idx_str=idx_str)
                        #     f.write(str)
                        #     latex_WmWm_delta(W_middle, [a3, a2], [a1, a4], mbpt, f)
                            # for t1 in rsimp_oneel_delta:
                            #     f.write('\n')
                            #     se = '{t1}'.format(t1=t1)
                            #     f.write(se)
                            #     f.write('\n')



                        # print('')

                        idx_str = a1+a2+a3+a4
                        print('a witam idx_str', idx_str)
                        print('END KLASY')
                        vspace(1)
                        print('wyraz rsimp_oneel_exch z ', idx_str)
                        vspace(0)
                        for g1 in  rsimp_oneel_exch:
                            print (g1)
                        vspace(1)
                        print('wyraz rsimp_oneel_delt z', idx_str)
                        vspace(0)
                        for g1 in  rsimp_oneel_delta:
                            print (g1)
                        vspace(1)
                        print('endend' , '- dodaje do slownika!')
                        Commutators_dict[idx_str] = Gm_commutators
                        Oneel_exch_dict[idx_str] = rsimp_oneel_exch
                        Oneel_delta_dict[idx_str] = rsimp_oneel_delta
                        print('')
                        print('WYSZLY NASTEPUJACE KOMUTATORY z wyrazu z P', idx_str)
                        print('')
                        for x in range(0, len(Gm_middle)):
                            print('G_middle', Gm_middle[x])
                            print('komutatory Gm_middle1')
                            for y in Gm_commutators[x][0]:
                                print(y)
                            print('komutatory Gm_middle2')
                            for y in Gm_commutators[x][1]:
                                print(y)

                        # print('')
                        # print('WYSZLY NASTEPUJACE KOMUTATORY z wyrazu exch-ttio', idx_str)
                        # for t1 in rsimp_oneel_exch:
                        #     print('too', t1)
                        # print('WYSZLY NASTEPUJACE KOMUTATORY z wyrazu delta-ssio', idx_str)
                        # for t1 in rsimp_oneel_delta:
                        #     print('sioo', t1)
                        # print('')
                        # if (idx_str == 'abcd'):
                        #     sys.exit(0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Gm_commutators = commutators_cumulant_ccd(obsy, obs, Gm_middle)
    # print('')
    # print('WYSZLY NASTEPUJACE KOMUTATORY')
    # print('')
    # for x in range(0, len(Gm_middle)):
    #     print('G_middle', Gm_middle[x])
    #     print('komutatory Gm_middle2')
    #     for y in Gm_commutators[x][0]:
    #         print(y)
    #     print('komutatory Gm_middle1')
    #     for y in Gm_commutators[x][1]:
    #         print(y)

    #     print('')

    # G2_middle, G2_commutators = commutators_density_ground_ccd(obsc, Gm2, 'T_list', 'S_list')
    # G1_middle, G1_commutators = commutators_density_ground_ccd(obscy, Gm1, 'T_list', 'S_list')

    # print('czesc Gm2')
    # G2_f12 = arithmetic_string()
    # for x in G2_commutators:
    #     for y in x:
    #         print(y)
    #         G2_f12.append(y)
    #     print('')

    # print('czesc Gm1')                
    # G1_f12 = arithmetic_string()
    # for x in G1_commutators:
    #     for y in x:
    #         print(y)
    #         G1_f12.append(y)
    #     print('')

    
    if pick == "dump":
        print('brabra')
        rint_big = integrate_Gm_one(Commutators_dict, Oneel_exch_dict, Oneel_delta_dict, Gm_middle, "ccsd", mbpt)

        for key in rint_big:
            rint = rint_big[key]
            rint.exec_delta()
            print('to bede wpisywac key', key)
            for x in rint:
                print(x)
            print('')
            if len(rint) > 0:

                latex_int = key+"\n"
            
                latex_int +="\equa{\n"

                for x in rint:
                    latex_int += f"& {x} \\\\"
                    latex_int += "\n"                                    
                latex_int += "\n}"
                f.write(latex_int)

        
        print('')
        if method == 'ccd':
            pick_cumulant = open('./pickle/cumulant_ccd_{mbpt}.pkl'.format(mbpt=mbpt), open_flags)
        elif method == 'cc3':
            pick_cumulant = open('./pickle/cumulant_cc3_{mbpt}.pkl'.format(mbpt=mbpt), open_flags)
        elif method == 'f12':
            pick_cumulant = open('./pickle/cumulant_f12_{mbpt}.pkl'.format(mbpt=mbpt), open_flags)            
            
        
        pickle.dump(rint_big, pick_cumulant)
        

    elif pick == "load":
        if method == 'ccd':
            load_cumulant = open('./pickle/cumulant_ccd_{mbpt}.pkl'.format(mbpt=mbpt), open_flags)
        elif method == 'cc3':
            load_cumulant = open('./pickle/cumulant_cc3_{mbpt}.pkl'.format(mbpt=mbpt), open_flags)
        elif method == 'f12':
            load_cumulant = open('./pickle/cumulant_f12_{mbpt}.pkl'.format(mbpt=mbpt), open_flags)
            
        rint_big = pickle.load(load_cumulant)

        print('LOADLOADLOAD')
        for key in rint_big:
            rint = rint_big[key]
            print('key', key)
            for y in rint:
                print(y)
            print('')
            #    print('pick', pick)

    vspace(0)
    print('i to jest nowa czesc')

    for key in rint_big:
        rint = rint_big[key]
        print('klucz', key)
        fx_list = [char for char in key]
        res_x = identify_interm_X_f12(rint)
        res_x = identify_interm_Ft_f12(res_x)
        print('po ft')
        for i in range(0, len(res_x)):
            res_x[i].summation = []
            # rint[i].coefficient.insert(0, OBSERVABLE_AX)
            # rint[i].coefficient_idx.insert(0, fx_list)
            # rint[i].summation += fx_list
            # rint[i].exec_delta()
            print(res_x[i])
        print('')
    sys.exit(0)

    

    sys.exit(0)
        # rsimp = simplify(rint)
        # print('result pooo simplify ground ground')
        # kk = 1
        # for yy in rsimp:
        #     print(kk, yy)
        #     kk += 1
        # print('')

        # print('teraz bede zmieniac x * y na x z czterema indeksami')

        # for i in range(0, len(rsimp)):
        #     x_j = 100
        #     y_j = 100
        #     print('przed', rsimp[i])
        #     coef_list = []
        #     for j in range(0, len(rsimp[i].coefficient)):
        #         if rsimp[i].coefficient[j] == OBSERVABLE_X:
        #             print('tak, mam observable')
        #             x_j = deepcopy(j)
        #             print('xj', x_j)
        #             for k in rsimp[i].coefficient_idx[j]:
        #                 coef_list.append(k)
        #             print('coef list', coef_list)
        #             if x_j > y_j:
        #                 x_j = x_j - 1
        #                 print('xj2', x_j)
        #         if rsimp[i].coefficient[j] == OBSERVABLE_Y:
        #             y_j = deepcopy(j)
        #             for k in rsimp[i].coefficient_idx[j]:
        #                 coef_list.append(k)
        #             if y_j > x_j:
        #                 y_j = y_j - 1
        #     print('po dodaniu', rsimp[i])
        #     del rsimp[i].coefficient[x_j]
        #     del rsimp[i].coefficient[y_j]
        #     del rsimp[i].coefficient_idx[x_j]
        #     del rsimp[i].coefficient_idx[y_j]

        #     rsimp[i].coefficient.append(OBSERVABLE_X)
        #     rsimp[i].coefficient_idx.append(coef_list)
        #     print('po', rsimp[i])
           
        # print('')
                    

        # # rsimp.exec_delta()
        # # rsimp.cleanup()

        # # rsimp_ost = simplify(rsimp)
        # rsimp_ost = rsimp
        # print('i to jest ostateczny wyniczor po symplifikacji')    
        # k = 0
        # for x in rsimp_ost:
        #     print(k, "   ", x)
        #     k += 1

        # print('koszty bez intermediates cumulant')
        # for x in rsimp_ost:
        #     memidx, costidx = memcost(x)
        #     print(costidx[0], ',', costidx[1], ',', costidx[2], costidx[0] + costidx[1]+costidx[2])

        # f.write('\\newpage')
        # print('wynik ostateczny latex')
        # latex_with_cost(f,1, rsimp_ost, name=False)
        # print('')

        # prefix = 'cumulant_{mbpt}'.format(mbpt=mbpt)
        # excl = []

        # basket_outer, basket_nointerm, list_of_int, n_max, all_hash, mem_dict, \
        #     interm_fx, xfx_dict  = factorize_driver(rsimp_ost, excl, COST=6, is_x_excluded=True, filel=f)


        # intermediates_to_fortran('ccsd', prefix, basket_outer, basket_nointerm, list_of_int, n_max, \
        #                          all_hash, mem_dict, interm_fx, xfx_dict)

        # basket_outer_s2, basket_nointerm_s2, list_of_int_s2, n_max_s2, all_hash_s2, mem_dict_s2, \
        #     interm_fx_s2, xfx_dict_s2  = factorize_driver(rsimp_ost_s2, [], COST=6)


        

def generate_density_matrix_ground_f12_driver(pick, method):

    f = open("Density_matrix_f12.tex", 'w')
    s_preamble = tex_preamble
    s_beginning = """
    \\begin{document}"""
    f.write(s_preamble)
    f.write(s_beginning)

    for order in range(2, 3):
        print('teraz bedzie rzad', order)
        s2 = "\n"
        if order > 0:
            s2 += "\\newpage \n"
        s2 += f"Macierz gęstości w rzędzie {order}\\\\"
        s2 += "\n"
        f.write(s2)

        list_of_int = []
#        choosing_pool = [CABS, virtual, occupied]
        choosing_pool = [completev, occupied]
        print(choosing_pool)

        banned = []
        start = 0
        for x1 in choosing_pool:
            a1 = free_idx(x1, [])
            banned.append(a1)
            for x2 in choosing_pool:
                a2 = free_idx(x2, [a1])
                
#                op1 = operat1([a1, a2], 's')
#                op1 = operat1x([a1, a2], 's')
                op1 = operat1xs([a1, a2], 's')
                print('opo', op1)

                idx_str = a1+a2
                print('TERAZ LICZE KOMUTATORY', idx_str)
                # f.write("\n")
                # zz = "{idx_str}".format(idx_str=idx_str)
                # f.write(zz)
                # f.write("\n")
                rsimp, W_middle, minilist_of_int = generate_density_matrix_ground_f12(pick, \
                            method, order, f, onlyrsimp = False, operator=op1, start=start, idx_str=idx_str)
                print('minilisti', minilist_of_int)
                if minilist_of_int != []:
                    list_of_int.append(minilist_of_int)
                start +=1
                print('gammaNO I ZROBIONE', idx_str)
                print('')

        for x in list_of_int:
            for y in list_of_int:
                print(y)

                
        f.write("\n")
        
        f.write("\\hrule")
        f.write("\n")

        
    f.write("\end{document}")
    f.close()

def generate_density_matrix_ground_f12(pick, method, mbpt, f, onlyrsimp = False, operator=obs, start=0, idx_str=''):

    print('LICZE DENSITY DLA OPERATORA', operator)
    open_flags = ""
    if pick == "dump":
        open_flags = "wb"
    elif pick == "load":
        open_flags = "rb"

    Wm = generate_Wm_ground_f12(mbpt, 'ccsd')
    print('')
    W_middle_f12_all_mbpt = generate_W_middle_for_Wm1(Wm)

    W_middle_f12 = []
    for x in W_middle_f12_all_mbpt:
        if x['mbpt'] == mbpt: 
            print(x, x['mbpt'])
            W_middle_f12.append(x)

    vspace(0)
    print('oto jest latex')
    vspace(0)
    latex_W_middle_Wm1(W_middle_f12, mbpt, False, f, start)

    #
    # !!!!!!!! WYBIERAM TYLKO DANY WYRAZ Z RZEDU
    #
#    W_middle_f12 = [W_middle_f12[3]]
#sys.exit(0)
    #
    #
    #

    print('oto jest latex SKROCONY')
    vspace(0)
    latex_W_middle_Wm1(W_middle_f12, mbpt, False, f, start)


    
    print('koniec latexa')
    W_middle, W_commutators, rsimp = commutators_density_ground_f12(operator, W_middle_f12, 'Wm1_T_list', 'Wm1_S_list')
    print('po wyjsciu z commutators')
    if onlyrsimp == True:
        print('Return Only Rsimp, do not group, do not change s operators')
        Wm1_commutators = []
        minilist = []
        return rsimp, W_middle, minilist
        

    print("I'm removing all instances of ff{akbl}, where a and b are both wirtual indices")
    print("According to definition of F in Shiozaki article, these instances are equal to zero")
    for x in rsimp:
        for j in range(0, len(x.coefficient)):
            y = x.coefficient[j]
            if y == F12_TWOEL or y ==F12_TWOEL_COMBO:
                n_of_virt = 0
                for z in x.coefficient_idx[j]:
                    if z in virtual:
                        n_of_virt +=1
                if n_of_virt >=2:
                    x.num_factor =0

    rsimp.exec_delta()
    rsimp.cleanup()
    
    res_x = identify_interm_X_f12(rsimp)
    print('po find X')
    k = 0
    for p in res_x:
        print(k, p)
        k +=  1
    print('')

#    res_x = identify_interm_Z_f12(res_x)
    print('po find Z')
    k = 0
    for p in res_x:
        print(k, p)
        k +=  1
    print('')


    res = arithmetic_string()
    for x in res_x:
        res = res + arithmetic_string(x)
    print('')
    cabs_res = res.cabstransform()
    print('po cabs transform s1')
    z = 0
    for x in cabs_res:        
        print(z, x)
        z += 1
    res_ft = identify_interm_Ft_f12(cabs_res)
    print('to po find Ft')
    z = 0
    for x in res_ft:        
        print(z, x)
        z += 1
    
    res_ost = arithmetic_string()
    for x in res_ft:
        res_ost.append(x)
        
    rsimp_ost = simplify(res_ost)
    
    print('i to jest ostateczny wyniczor po symplifikacji')    
    k = 0
    for x in rsimp_ost:
        print(k, "   ", x)
        k += 1

    print('koszty bez intermediates density')
    for x in rsimp_ost:
        memidx, costidx = memcost(x)
        print(costidx[0], ',', costidx[1], ',', costidx[2], costidx[0] + costidx[1]+costidx[2])

    # f.write('\\newpage')
    # print('wynik ostateczny latex')
    # latex_with_cost(f,1, rsimp_ost, name=False)
    # print('')

    prefix = 'dens_{mbpt}'.format(mbpt=mbpt)
    excl = []
    if idx_str != "":
        for x in idx_str:
            excl.append(x)
    print(excl)
    basket_outer, basket_nointerm, list_of_int, n_max, all_hash, mem_dict, interm_fx, xfx_dict, list_of_int_xfx  = factorize_driver(rsimp_ost, excl, COST=6, is_x_excluded=True, filel=f, idx_str=idx_str)

    print('all leny')
    print(len(basket_outer))
    print(len(basket_nointerm))
    print(len(list_of_int))

    print('finfin outer')
    for x in basket_outer:
        print(x)
    print('')
    print('finfin noint')
    for x in basket_nointerm:
        print(x)
    print('')

    print('interm_dict')
    interm_dict_list = []
    for x in range(0, len(list_of_int)):
        interm_dict = {}
        for y in range(0, len(list_of_int[x])):
            interm_dict['name'] = all_hash[list_of_int[x][y].binary_hash]
            interm_dict['xfx'] = list_of_int_xfx[x][y]
            interm_dict['ugg'] = list_of_int[x][y]
        print(interm_dict)
        interm_dict_list.append(interm_dict)
    # sys.exit(0)

    # intermediates_to_fortran('ccsd', prefix, basket_outer, basket_nointerm, list_of_int, n_max, \
    #                          all_hash, mem_dict, interm_fx, xfx_dict)
    
    # basket_outer_s2, basket_nointerm_s2, list_of_int_s2, n_max_s2, all_hash_s2, mem_dict_s2, \
    #     interm_fx_s2, xfx_dict_s2  = factorize_driver(rsimp_ost_s2, [], COST=6)

    return rsimp_ost, W_middle_f12, interm_dict_list

def execute_density_matrix_ground_ground(method, MBPT_order):

    # Density matrix in CCSD and CC3 correct through the third 
    # order MBPT
    #
    # CCSD
    # Y + [S1*, Y] + [Y,T1] 
    # + [S1*,[Y,T2]] + [S2*,[Y,T2]]
    # CC3                                                                                       
    # Y + [S1*, Y] + [Y,T1]                                                                
    # + [S1*,[Y,T2]] + [S2*,[Y,T2]] 
    # + [S2*,[Y,T3] 

    # Density matrix in CCSD and CC3 terms up to the fourth
    # order MBPT, but not correct through the fourth order
    #
    # CCSD
    # Y + [S1*, Y] + [Y,T1] 
    # + [S1*,[Y,T2]] + [S2*,[Y,T2]]+ [S1*,[Y,T1]]
    # + 0.5[S2*,[[Y,T1], T2]]
    # CC3                                                                                       
    # Y + [S1*, Y] + [Y,T1]                                                                
    # + [S1*,[Y,T2]] + [S2*,[Y,T2]] + [S1*,[Y,T1]]
    # + 0.5[S2*,[[Y,T1], T2]]
    # + [S2*,[Y,T3] + [S3*,[Y,T3]] + 0.5[S3*,[[Y,T2],T2]] 


    print('zaczynam1')
    if(MBPT_order == 0):
        Y     = arithmetic_string(nobs)
        S1ca  = arithmetic_string(s1c)
        S1capol  = arithmetic_string(s1c).scale(0.5)
        S1casz  = arithmetic_string(s1c).scale(1.0/6.0)
        S2ca  = arithmetic_string(s2c)
        S3ca  = arithmetic_string(t3c)
        S3capol  = arithmetic_string(t3c).scale(0.5)
        YT1   = evaluate(nobs, t1)
        YT2   = evaluate(nobs, t2)
        YT3   = evaluate(nobs, t3)
        YT1T1 = evaluate(nobs, t1, t1)
        YT1T2 = evaluate(nobs, t1, t2)
        YT2T2 = evaluate(nobs, t2, t2)
        if(method =='ccsd'):
            print('zaczynam2')
            r = []
            r.append(Y)                                                # 0                                                   
            r.append(S1ca * Y + YT1 + S2ca * YT2)                      # 2                                                     
            r.append(S1ca * YT2)                                       # 3                                                          
            r.append(S1ca * YT1 + S2ca * YT1T2)                        # 4                                                 
            r.append(S1capol * S1ca * YT2 + S1capol * S2ca * YT2T2)    # 5                                                 
            r.append(S1capol * YT1T1)                                  # 6                                                 
            r.append(S1capol * S1ca * YT1T2)                           # 7                                                   
            r.append(S1casz * S1capol * S1ca *YT2T2)                      # 8          
            

        elif(method == 'cc3'):
            r = []
            r.append(Y)                                                # 0
            r.append(S1ca * Y + YT1 + S2ca * YT2)                      # 2
            r.append(S1ca * YT2 + S2ca * YT3)                          # 3
            r.append(S1ca * YT1 + S2ca * YT1T2 +  S3ca * YT3 + \
                         S3capol * YT2T2)                             # 4
            r.append(S1capol * S1ca * YT2 + S1capol * S2ca * YT2T2 + \
                        S1ca * S2ca * YT3)                            # 5
            r.append(S1capol * YT1T1 + S1capol * S1ca * YT3)           # 6
            r.append(S1capol * S1ca * YT1T2)                           # 7
            r.append(S1casz * S1capol * S1ca * YT2T2 + \
                         S1casz * S1ca * S1ca * YT3)                   # 8

        # print('zaczynam3')
        # for x in r:
        #     for y in x:
        #         print('xyxy', y)
        
        for i in range(0, len(r)):
#            print('czaz', i, len(r))
            if(i==0):
                j = 0
            else:
                j = i + 1
            print('-----------------------------------------------------------')
            print(i, j)
            print('-----------------------------------------------------------')
            # if (j==3):
            #     sys.exit(0)
            print('zaczynam5')
            rint = r[i].integrate()
            print('zaczynam6')
            rint.exec_delta()
            print('zaczynam7')
            for x in rint:
                print('xxx', x)
            rsimp = simplify(rint)
            print('zaczynam8')
            

            block_oo_diag, block_vv_diag, rsimp_without_diag = fill_dm_blocks_diag(rsimp)
            
            block_oo, block_ov, block_vo, block_vv = fill_dm_blocks(rsimp_without_diag)
            for x in block_oo:
                s = []
                for y in x.summation:
                    s.append(y)
                for y in x.coefficient_idx:
                    for z in y:
                        if z not in s:
                            s.append(z)
                sv = []
                so = []
                for y in s:
                    if y in occupied:
                        so.append(y)
                    elif y in virtual:
                        sv.append(y)
                print(x)
                if len(x.coefficient) > 2:
                    svi, svo = cost_after_intermediates(deepcopy(x))
                    print(j, 'sniez-oo', x, '         ', sv, so, len(sv), len(so),svi, svo, len(svi), len(svo))
                else:
                    print(j, 'sniez-oo', x, '         ', sv, so, len(sv), len(so))
            print('block_ov')
            for x in block_ov:
                s = []
                for y in x.summation:
                    s.append(y)
                for y in x.coefficient_idx:
                    for z in y:
                        if z not in s:
                            s.append(z)
                sv = []
                so = []
                for y in s:
                    if y in occupied:
                        so.append(y)
                    elif y in virtual:
                        sv.append(y)
                if len(x.coefficient) > 2:
                    svi, svo = cost_after_intermediates(deepcopy(x))
                    print(j, 'sniez-ov', x, '         ', sv, so, len(sv), len(so),svi, svo, len(svi), len(svo))
                else:
                    print(j, 'sniez-ov', x, '         ', sv, so, len(sv), len(so))

            print('block_vo')
            for x in block_vo:
                s = []
                for y in x.summation:
                    s.append(y)
                for y in x.coefficient_idx:
                    for z in y:
                        if z not in s:
                            s.append(z)
                sv = []
                so = []
                for y in s:
                    if y in occupied:
                        so.append(y)
                    elif y in virtual:
                        sv.append(y)
                if len(x.coefficient) > 2:
                    svi, svo = cost_after_intermediates(deepcopy(x))
                    print(j, 'sniez-vo', x, '         ', sv, so, len(sv), len(so),svi, svo, len(svi), len(svo))
                else:
                    print(j, 'sniez-vo', x, '         ', sv, so, len(sv), len(so))
            print('block_vv')
            for x in block_vv:
                s = []
                for y in x.summation:
                    s.append(y)
                for y in x.coefficient_idx:
                    for z in y:
                        if z not in s:
                            s.append(z)
                sv = []
                so = []
                for y in s:
                    if y in occupied:
                        so.append(y)
                    elif y in virtual:
                        sv.append(y)
                if len(x.coefficient) > 2:
                    svi, svo = cost_after_intermediates(deepcopy(x))
                    print(j, 'sniez-vv', x, '         ', sv, so, len(sv), len(so),svi, svo, len(svi), len(svo))
                else:
                    print(j, 'sniez-vv', x, '         ', sv, so, len(sv), len(so))
        #        sys.exit(0)

            function_template_density_matrix_ground(j, block_oo, block_ov, block_vo, block_vv, \
                                                        block_oo_diag, block_vv_diag, method)

    else:


        if(MBPT_order == 3):
            if(method == 'ccsd'):
                r = (arithmetic_string(nobs)
                     + arithmetic_string(s1c) * arithmetic_string(nobs)
                     + evaluate(nobs, t1)
                     + arithmetic_string(s1c) * evaluate(nobs, t2)
                     + arithmetic_string(s2c) * evaluate(nobs, t2))

            elif(method == 'cc3'):
                r = (arithmetic_string(nobs)
                     + arithmetic_string(s1c) * arithmetic_string(nobs)
                     + evaluate(nobs, t1)
                     + arithmetic_string(s1c) * evaluate(nobs, t2)
                     + arithmetic_string(s2c) * evaluate(nobs, t2)
                     + arithmetic_string(s2c) * evaluate(nobs, t3))

            else:
                print('NO SUCH METHOD')
                sys.exit(0)

        elif(MBPT_order == 4):
            Y     = arithmetic_string(nobs)
            T1a    = arithmetic_string(t1)
            S1ca  = arithmetic_string(s1c)
            S1capol  = arithmetic_string(s1c).scale(0.5)
            S1casz  = arithmetic_string(s1c).scale(1.0/6.0)
            S2ca  = arithmetic_string(s2c)
            S3ca  = arithmetic_string(t3c)
            S3capol  = arithmetic_string(t3c).scale(0.5)
            YT1   = evaluate(nobs, t1)
            YT2   = evaluate(nobs, t2)
            YT3   = evaluate(nobs, t3)
            YT1T1 = evaluate(nobs, t1, t1)
            YT1T2 = evaluate(nobs, t1, t2)
            YT2T2 = evaluate(nobs, t2, t2)

            if(method =='ccsd'):
                r = (arithmetic_string(Y)                             # 0
                     + (S1ca * Y + YT1 + S2ca * YT2)                      # 2
                     + (S1ca * YT2)                                       # 3
                     + (S1ca * YT1 + S2ca * YT1T2))                        # 4

            elif(method == 'cc3'):
                r = (arithmetic_string(Y)                             # 0
                     + (S1ca * Y + YT1 + S2ca * YT2)                      # 2
                     + (S1ca * YT2 + S2ca * YT3)                          # 3
                     + (S1ca * YT1 + S2ca * YT1T2 +  S3ca * YT3 + \
                            S3capol * YT2T2))                              # 4
            else:
                print('NO SUCH METHOD')
                sys.exit(0)

        print('PIERWSZE')
        for x in r:
            print(x)
        print('')

        rint = r.integrate()
        rint.exec_delta()
        print('przed simplify')
        for x in rint:
            print(x)
        print('')
        rsimp = simplify(rint)
        print('')
        print('przed usunieciem diag')
        print('')

        for x in rsimp:
            print(x)
        print('')

        block_oo_diag, block_vv_diag, rsimp_without_diag = fill_dm_blocks_diag(rsimp)
        print('Blok oo')
        print('')
        for x in block_oo_diag:
            print(x)
        print('')
        print('Blok vv')
        for x in block_vv_diag:
            print(x)
        print('rsimp po usunieciu diag')
        for x in rsimp_without_diag:
            print(x)
        print('')

        block_oo, block_ov, block_vo, block_vv = fill_dm_blocks(rsimp_without_diag)
        # function_template_density_matrix_ground(MBPT_order, block_oo, block_ov, block_vo, block_vv, \
        #                                             block_oo_diag, block_vv_diag, method)
        print('')
        print('Blok oo')
        print('')
        for x in block_oo:
            print(x)
        print('')
        print('Blok ov')
        print('')
        for x in block_ov:
            print(x)
        print('')
        print('Blok vo')
        print('')
        for x in block_vo:
            print(x)
        print('')
        print('Blok vv')
        print('')
        for x in block_vv:
            print(x)


def find_largest_cost(oo, ov, vo, vv, interm):
    
    cost_list_oo = compute_cost(oo, 'oo', interm)
    cost_list_ov = compute_cost(ov, 'ov', interm)
    cost_list_vo = compute_cost(vo, 'vo', interm)
    cost_list_vv = compute_cost(vv, 'vv', interm)


    print('cost oo')
    for x in range(0, len(oo)):
        print(oo[x], cost_list_oo[x])
    print('cost ov')
    for x in range(0, len(ov)):
        print(ov[x], cost_list_ov[x])

    print('cost vo')
    for x in range(0, len(vo)):
        print(vo[x], cost_list_vo[x])

    print('cost vv')
    for x in range(0, len(vv)):
        print(vv[x], cost_list_vv[x])


    print('cost_list_oo', cost_list_oo)
    print('cost_list_ov', cost_list_ov)
    print('cost_list_vo', cost_list_vo)
    print('cost_list_vv', cost_list_vv)
    cost_list = cost_list_oo + cost_list_ov + cost_list_vo + cost_list_vv
    
    cost_old = 0
    
    for x in cost_list:
        cost_new = x[0] * 20 + x[1] * 0.1
        if cost_new > cost_old:
            cost_old = cost_new
            cost = x

    return cost_old, cost


def compute_cost(ars, block, interm):

    cost_list = []
    if interm == False:

        for x in ars:    
            if block == 'oo':
                oro = 2
                orv = 0
            elif block == 'ov' or block == 'vo':
                oro = 1
                orv = 1
            elif block == 'vv':
                oro = 0
                orv = 2
            for q in x.summation:
                if q in virtual:
                    orv += 1
                elif q in occupied:
                    oro += 1
            cost_list.append([orv, oro])
    if interm == True:
        for x in ars:    
            oro = 0
            orv = 0
        
            for q in x['idx_fx']:
                if q in virtual:
                    orv += 1
                elif q in occupied:
                    oro += 1

            for q in x['idx_sum']:
                if q in virtual:
                    orv += 1
                elif q in occupied:
                    oro += 1
                cost_list.append([orv, oro])

    return cost_list

def compute_cost_interm(ars):

    cost_list = []
    for x in ars:
        orv = 0
        oro = 0
        for q in x['interm'].summation:
            if q in virtual:
                        orv += 1
            elif q in occupied:
                        oro += 1
        cost_list.append([orv, oro])

    return cost_list


def compute_cost_direct_ugg(e):

    orv = 0
    oro = 0
    for q in e.summation:
        if q in virtual:
            orv += 1
        elif q in occupied:
            oro += 1
    cost_list = [orv, oro]

    return cost_list


def compute_cost_direct(ars):

    cost_list = []
    for x in ars:
        orv = 0
        oro = 0
        for q in x.summation:
            if q in virtual:
                orv += 1
            elif q in occupied:
                oro += 1
        cost_list.append([orv, oro])
        
    return cost_list

def compute_cost_ov(oo, ov, vo, vv):

    cost_list = []

    cost_oo = compute_cost_direct(oo)
    cost_ov = compute_cost_direct(ov)
    cost_vo = compute_cost_direct(vo)
    cost_vv = compute_cost_direct(vv)
    for x in range(0, len(cost_oo)):
        cost_oo[x] = [cost_oo[x][0],cost_oo[x][1]+2]
    for x in range(0, len(cost_ov)):
        cost_ov[x] = [cost_ov[x][0]+1,cost_ov[x][1]+1]
    for x in range(0, len(cost_vo)):
        cost_vo[x] = [cost_vo[x][0]+1,cost_vo[x][1]+1]
    for x in range(0, len(cost_vv)):
        cost_vv[x] = [cost_vv[x][0]+2,cost_vv[x][1]]

    return cost_oo, cost_ov, cost_vo, cost_vv

# Functions that creates density matrices
def execute_transition_gamma_dm(method, pick, mbpt, maxcluster, multiplicity, cumulative):

    #snorlax

    W_gamma = generate_commutators_for_gamma(mbpt, maxcluster, method, cumulative)

    print('SAME KOMUTATORY')
    for x in W_gamma:
        print('la', x)

    sys.exit(0)
    mbpt = 3

    open_flags = ""
    if pick == "dump":
        open_flags = "wb"        
    elif pick == "load":
        open_flags = "rb"

    if method == 'ccsd':
        out = open('./pickle/gr_exc/gamma_dm_ccsd.pkl', open_flags)
    elif method == 'cc3':
        out = open('./pickle/gr_exc/gamma_dm_cc3.pkl', open_flags)

            
    if pick == 'dump':
        
        r1, r21, r22, r31a, r31b = compute_gamma_dm(method)
        
        
        rsimp = int_and_simp_transition_gamma_dm(r1, r21, r22, r31a, r31b, method, multiplicity)
        
        pickle.dump(rsimp, out)
    elif pick == 'load':
        rsimp = pickle.load(out)
        
#    block_oo, block_ov, block_vo, block_vv = fill_dm_blocks_wm(intermediates)         
#    block_oo, block_ov, block_vo, block_vv = fill_dm_blocks_wm(rsimp)        

    block_oo, block_ov, block_vo, block_vv = fill_dm_blocks(rsimp)  
    
    cost_without_int, cost_lst = find_largest_cost(block_oo, block_ov, block_vo, block_vv, False)


    intermediates_dict_oo, block_oo, sum_list, k  = generate_intermediates_gaxi(block_oo, 'gamma', method, 0, 'oo')
    intermediates_dict_ov, block_ov, sum_list, k  = generate_intermediates_gaxi(block_ov, 'gamma', method, k, 'ov')
    intermediates_dict_vo, block_vo, sum_list, k  = generate_intermediates_gaxi(block_vo, 'gamma', method, k, 'vo')
    intermediates_dict_vv, block_vv, sum_list ,k  = generate_intermediates_gaxi(block_vv, 'gamma', method, k, 'vv')

    cost_with_int, cost_lst_with_int = find_largest_cost(block_oo, block_ov, block_vo, block_vv, False)
    cost_int, cost_lst_int = find_largest_cost(intermediates_dict_oo, intermediates_dict_ov, intermediates_dict_vo, intermediates_dict_vv, True)

    print(cost_without_int, cost_with_int, cost_int)
    print(cost_lst, cost_lst_with_int, cost_lst_int)

    intermediates_dict =  intermediates_dict_oo +  intermediates_dict_ov +  intermediates_dict_vo + intermediates_dict_vv

    for x in intermediates_dict:
        print(x['interm'], len(x['noninterm'].coefficient))


    print('block_oo')
    for x in block_oo:
        print(x)
    print('block_ov')
    for x in block_ov:
        print(x)
    print('block_vo')
    for x in block_vo:
        print(x)
    print('block_vv')
    for x in block_vv:
        print(x)

    print(len(block_oo))
    print(len(block_ov))
    print(len(block_vo))
    print(len(block_vv))
    
    function_template_wm_intermediates(intermediates_dict, method, 'gamma', multiplicity, mbpt)
    function_template_transition_gamma_dm_with_intermediates(block_oo, block_ov, block_vo, block_vv, method)

def execute_transition_xi_dm(method, pick):

    open_flags = ""
    if pick == "dump":
        open_flags = "wb"        
    elif pick == "load":
        open_flags = "rb"

    if method == 'ccsd':
        out = open('./pickle/gr_exc/xi_dm_ccsd.pkl', open_flags)
    elif method == 'cc3':
        out = open('./pickle/gr_exc/xi_dm_cc3.pkl', open_flags)

            
    if pick == 'dump':
        
        r1, r2, r3 = compute_xi_dm(method)

        rsimp = int_and_simp_transition_xi_dm(r1, r2, r3, method)
        
        pickle.dump(rsimp, out)

    elif pick == 'load':
        rsimp = pickle.load(out)


    block_oo, block_ov, block_vo, block_vv = fill_dm_blocks(rsimp) 



    cost_without_int, cost_lst = find_largest_cost(block_oo, block_ov, block_vo, block_vv, False)

    intermediates_dict_oo, block_oo, sum_list, k  = generate_intermediates_gaxi(block_oo, 'xi', method, 0, 'oo')
    intermediates_dict_ov, block_ov, sum_list, k  = generate_intermediates_gaxi(block_ov, 'xi', method, k, 'ov')
    intermediates_dict_vo, block_vo, sum_list, k  = generate_intermediates_gaxi(block_vo, 'xi', method, k, 'vo')
    intermediates_dict_vv, block_vv, sum_list, k  = generate_intermediates_gaxi(block_vv, 'xi', method, k, 'vv')


    cost_with_int, cost_lst_with_int = find_largest_cost(block_oo, block_ov, block_vo, block_vv, False)
    cost_int, cost_lst_int = find_largest_cost(intermediates_dict_oo, intermediates_dict_ov, intermediates_dict_vo, intermediates_dict_vv, True)

    print(cost_without_int, cost_with_int, cost_int)
    print(cost_lst, cost_lst_with_int, cost_lst_int)




    intermediates_dict =  intermediates_dict_oo +  intermediates_dict_ov +  intermediates_dict_vo + intermediates_dict_vv

    for x in intermediates_dict:
        print(x['interm'], len(x['noninterm'].coefficient))
    sys.exit(0)


    # intermediates_dict, intermediates, sum_list  = generate_intermediates_gaxi(rsimp, 'xi', method)

    # block_oo, block_ov, block_vo, block_vv = fill_dm_blocks_wm(intermediates)

    # print('block_oo')
    # for x in block_oo:
    #     print(x)
    # print('block_ov')
    # for x in block_ov:
    #     print(x)
    # print('block_vo')
    # for x in block_vo:
    #     print(x)
    # print('block_vv')
    # for x in block_vv:
    #     print(x)

    function_template_wm_intermediates(intermediates_dict, method, 'xi')
    function_template_transition_xi_dm_with_intermediates(block_oo, block_ov, block_vo, block_vv, method)

def int_and_simp_transition_gamma_dm(r1, r21, r22, r31a, r31b, method, multiplicity):

#-----------------------------------------------
    if (multiplicity == 1) :
        e1 = eomr1
    elif (multiplicity == 3) :
        e1 = eomrr1_triplet
        
    r1j = r1 * arithmetic_string(e1)

    rint1 = r1j.integrate()
    rint1.exec_delta()

    rsimp1 = simplify(rint1)

#-----------------------------------------------
    rint2 = r22.integrate(ket = ['a', 'i', 'b', 'j'], ketspin = ['s', 's'])
    
    rint2.exec_delta()
    rsimp2 = simplify(rint2)

    tauaibj = ugg()
    tauaibj.operator_idx = [["a", "i"], ["b", "j"]]
    tauaibj.operator_type = ['s', 's']
    ketcm = evaluate(s1c, tauaibj) 
    ketcms = simplify(ketcm)
    ketcmr = r21 * ketcms
    rint21 = ketcmr.integrate() 
    rint2 += rint21

    rint2.exec_delta()

    rsimp2 = simplify(rint2)


    for x in range(0, len(rsimp2)):
        rsimp2[x].coefficient.append(EOM_CC_AMPLITUDE_R)
        rsimp2[x].coefficient_idx.append(['a','i', 'b', 'j'])
        rsimp2[x].summation.append('a')
        rsimp2[x].summation.append('b')
        rsimp2[x].summation.append('i')
        rsimp2[x].summation.append('j')
        rsimp2[x].num_factor *= 1./2.
        
    rsimp2.exec_delta()

    rsimp22 = simplify(rsimp2)

#-----------------------------------------------
    if method == 'ccsd':
        rsimp = rsimp1 + rsimp22
 
    if method == 'cc3':
        tauaibjck = ugg()
        tauaibjck.operator_idx = [["a", "i"], ["b", "j"], ['c', 'k']]
        tauaibjck.operator_type = ['s', 's', 's']
        ketcm2a = evaluate(s2c, tauaibjck)
        ketcm2b = evaluate(s1c, tauaibjck)

        ketcm2as = simplify(ketcm2a)
        ketcm2bs = simplify(ketcm2b)
        ketcm2ar = r31a * ketcm2as
        ketcm2br = r31b * ketcm2bs

        rint3a = ketcm2ar.integrate()
        rint3a.exec_delta()
        rsimp3a = simplify(rint3a)

        rint3b = ketcm2br.integrate()
        rint3b.exec_delta()
        rsimp3b = simplify(rint3b)

        rsimp3 = rsimp3a + rsimp3b

        for x in range(0, len(rsimp3)):
            rsimp3[x].coefficient.append(EOM_CC_AMPLITUDE_R)
            rsimp3[x].coefficient_idx.append(['a','i', 'b', 'j', 'c', 'k'])
            rsimp3[x].summation.append('a')
            rsimp3[x].summation.append('b')
            rsimp3[x].summation.append('c')
            rsimp3[x].summation.append('i')
            rsimp3[x].summation.append('j')
            rsimp3[x].summation.append('k')
            rsimp3[x].num_factor *= 1./6.
            
        rsimp3.exec_delta()
        rsimp3.clear_fixed()
        rsimp33 = simplify(rsimp3)


        rsimp = rsimp1 + rsimp2 + rsimp33

    return rsimp


def int_and_simp_transition_xi_dm(r1, r2, r3, method):

    e = deepcopy(eoml1)
    e.transpose()
    r1j = arithmetic_string(e) * r1
    rint1 = r1j.integrate().scale(0.5)

    rint1.exec_delta()

    rsimp1 = simplify(rint1)

    rint2a = r2.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s', 's'])
    rint2b = r2.integrate(bra = ['a', 'j', 'b', 'i'], braspin = ['s', 's'])

    rint2a.exec_delta()
    rint2b.exec_delta()

    rint2  = rint2a.scale(1./3.) + rint2b.scale(1./6.)

    temp = deepcopy(rint2)
    for x in temp:
        x.new_delta("a", "b")
        x.new_delta("i", "j")
    rint2 += temp.scale(-1./2.)

    rint2.exec_delta()

    rsimp2 = simplify(rint2)

    
    print('rsimp2')
    for x in rsimp2:
        print(x)

    for x in range(0, len(rsimp2)):
        rsimp2[x].coefficient.append(EOM_CC_AMPLITUDE_L)
        rsimp2[x].coefficient_idx.append(['a','i', 'b', 'j'])
        rsimp2[x].summation.append('a')
        rsimp2[x].summation.append('b')
        rsimp2[x].summation.append('i')
        rsimp2[x].summation.append('j')
        rsimp2[x].num_factor *= 1./2.
        

    rsimp2.exec_delta()

    rsimp22 = simplify(rsimp2)

    if method == 'ccsd':
        rsimp = rsimp1 + rsimp22

    elif method == 'cc3':

        rint3 = r3.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's'])
        rint3.exec_delta()
        rsimp3 = simplify(rint3)

        for x in range(0, len(rsimp3)):
            rsimp3[x].coefficient.append(EOM_CC_AMPLITUDE_L)
            rsimp3[x].coefficient_idx.append(['a','i', 'b', 'j', 'c', 'k'])
            rsimp3[x].summation.append('a')
            rsimp3[x].summation.append('b')
            rsimp3[x].summation.append('c')
            rsimp3[x].summation.append('i')
            rsimp3[x].summation.append('j')
            rsimp3[x].summation.append('k')
            rsimp3[x].num_factor *= 1./6.

        rsimp3.exec_delta()
        rsimp3.clear_fixed()
        rsimp33 = simplify(rsimp3)

        rsimp = rsimp1 + rsimp22 + rsimp33
        
    return rsimp


# Auxilliary functions

def fill_dm_blocks_diag(rsimp):

    block_oo_diag = arithmetic_string()
    block_vv_diag = arithmetic_string()
    block_oo_diag_simp = arithmetic_string()
    block_vv_diag_simp = arithmetic_string()

    rsimp_without_diag = arithmetic_string()

    for x in range(0,len(rsimp)):
        print('x', x, rsimp[x])
        for i in range(0, len(rsimp[x].coefficient)):
            if(rsimp[x].coefficient[i]==OBSERVABLE_X or rsimp[x].coefficient[i]==NONSYM_MTRX or rsimp[x].coefficient[i]==OBSERVABLE_X_ASYM):
                jj = deepcopy(i)
                obs_idx = rsimp[x].coefficient_idx[i]
                block = find_block(obs_idx)
                print('s2', rsimp[x])
                if(block != 'oo_diag' and block != 'vv_diag'): 
                    print('s3', rsimp[x])
                    rsimp_without_diag.append(deepcopy(rsimp[x]))
                sm = []
                print('s4', rsimp[x])
                for j in rsimp[x].summation:
                    if j not in rsimp[x].coefficient_idx[i]:
                        sm.append(j)
                rsimp[x].summation = sm
        print('s5', rsimp[x])
        rsimp[x].remove_coefficient(OBSERVABLE_X)
        rsimp[x].remove_coefficient(NONSYM_MTRX)
        rsimp[x].remove_coefficient(OBSERVABLE_X_ASYM)
        print('iiii', jj)
        rsimp[x].remove_coefficient_idx(jj)
        print('s6', rsimp[x])
        if block == 'oo_diag':
            print('ooo')
            rsimp[x].recreate_fixed(obs_idx, 'oo')
            block_oo_diag.append(rsimp[x])
            rsimp[x].clear_fixed()
            rsimp[x].establish_fixed()
            block_oo_diag_simp = simplify(block_oo_diag)
        elif block == 'vv_diag':
            print('vooo', rsimp[x])
            rsimp[x].recreate_fixed(obs_idx, 'vv')
            print('1a', rsimp[x])
            block_vv_diag.append(rsimp[x])
            print('1b', rsimp[x])
            rsimp[x].clear_fixed()
            print('1c', rsimp[x])
            rsimp[x].establish_fixed()
            block_vv_diag_simp = simplify(block_vv_diag)

    return block_oo_diag, block_vv_diag, rsimp_without_diag

def fill_dm_blocks_dens_ground(rsimp):

    print(len(rsimp))


    test = arithmetic_string()

    block_oo = arithmetic_string()
    block_ov = arithmetic_string()
    block_vo = arithmetic_string()
    block_vv = arithmetic_string()
    block = ''

    oo = 0
    vv = 0
    ov = 0
    vo = 0

    for x in range(0,len(rsimp)):

        laciaty = deepcopy(rsimp[x])
        for i in range(0, len(rsimp[x].coefficient)):
            if(rsimp[x].coefficient[i]==OBSERVABLE_X or rsimp[x].coefficient[i]==NONSYM_MTRX or rsimp[x].coefficient[i]==OBSERVABLE_X_ASYM):
                obs_idx = deepcopy(rsimp[x].coefficient_idx[i])
                block = find_block(obs_idx)

        if block == 'oo':
            block_oo.append(rsimp[x])
        elif block =='ov':
            block_ov.append(rsimp[x])
        elif block =='vo':
            block_vo.append(rsimp[x])
        elif block =='vv':
            block_vv.append(rsimp[x])

    return block_oo, block_ov, block_vo, block_vv



def fill_dm_blocks(rsimp):

    print(len(rsimp))


    test = arithmetic_string()

    block_oo = arithmetic_string()
    block_ov = arithmetic_string()
    block_vo = arithmetic_string()
    block_vv = arithmetic_string()
    block_oo_simp = arithmetic_string()
    block_ov_simp = arithmetic_string()
    block_vo_simp = arithmetic_string()
    block_vv_simp = arithmetic_string()
    block = ''

    oo = 0
    vv = 0
    ov = 0
    vo = 0

    for x in range(0,len(rsimp)):

        laciaty = deepcopy(rsimp[x])
        for i in range(0, len(rsimp[x].coefficient)):
            if(rsimp[x].coefficient[i]==OBSERVABLE_X or rsimp[x].coefficient[i]==NONSYM_MTRX or rsimp[x].coefficient[i]==OBSERVABLE_X_ASYM):
                obs_idx = deepcopy(rsimp[x].coefficient_idx[i])
                block = find_block(obs_idx)

                sm = []
                for j in rsimp[x].summation:
                    if j not in rsimp[x].coefficient_idx[i]:
                        sm.append(j)
                rsimp[x].summation = sm

                rsimp[x].remove_coefficient_idx(i)
        rsimp[x].remove_coefficient(OBSERVABLE_X)
        rsimp[x].remove_coefficient(NONSYM_MTRX)
        rsimp[x].remove_coefficient(OBSERVABLE_X_ASYM)

        if block == 'oo':

            oo += 1
            rsimp[x].recreate_fixed(obs_idx, 'oo')
            block_oo.append(rsimp[x])
            rsimp[x].clear_fixed()
            rsimp[x].establish_fixed()
            block_oo_simp = simplify(block_oo)

        elif block == 'ov':
            rsimp[x].recreate_fixed(obs_idx, 'ov')
            block_ov.append(rsimp[x])
            rsimp[x].clear_fixed()
            rsimp[x].establish_fixed()
            block_ov_simp = simplify(block_ov)
        elif block == 'vo':
            rsimp[x].recreate_fixed(obs_idx, 'vo')
            block_vo.append(rsimp[x])
            rsimp[x].clear_fixed()
            rsimp[x].establish_fixed()
            block_vo_simp = simplify(block_vo)
        elif block == 'vv':
            rsimp[x].recreate_fixed(obs_idx, 'vv')
            block_vv.append(rsimp[x])
            rsimp[x].clear_fixed()
            rsimp[x].establish_fixed()
            block_vv_simp = simplify(block_vv)
        elif block == 'po':
            rsimp = rsimp[x].general_split(obs_idx, 'po')
            rsimp[0].recreate_fixed(obs_idx, 'oo')
            rsimp[1].recreate_fixed(obs_idx, 'vo')
            block_oo.append(rsimp[0])
            rsimp[0].clear_fixed()
            rsimp[0].establish_fixed()
            block_oo_simp = simplify(block_oo)
            block_vo.append(rsimp[1])
            rsimp[1].clear_fixed()
            rsimp[1].establish_fixed()
            block_vo_simp = simplify(block_vo)
        elif block == 'op':
            rsimp = rsimp[x].general_split(obd_idx, 'op')
            rsimp[0].recreate_fixed(obs_idx, 'oo')
            rsimp[1].recreate_fixed(obs_idx, 'ov')
            block_oo.append(rsimp[0])
            rsimp[0].clear_fixed()
            rsimp[0].establish_fixed()
            block_oo_simp = simplify(block_oo)
            rsimp[1].clear_fixed()
            rsimp[1].establish_fixed()
            block_ov.append(rsimp[1])

            block_ov_simp = simplify(block_ov)
        elif block == 'pp':
            rsimp = rsimp[x].general_split(obd, idx, 'pp')
            rsimp[0].recreate_fixed(obs_idx, 'oo')
            rsimp[1].recreate_fixed(obs_idx, 'ov')
            rsimp[2].recreate_fixed(obs_idx, 'vo')
            rsimp[3].recreate_fixed(obs_idx, 'vv')

            block_oo.append(rsimp[0])
            rsimp[0].clear_fixed()
            rsimp[0].establish_fixed()
            block_oo_simp = simplify(block_oo)

            block_ov.append(rsimp[1])
            rsimp[1].clear_fixed()
            rsimp[1].establish_fixed()
            block_ov_simp = simplify(block_ov)

            block_vo.append(rsimp[2])
            rsimp[2].clear_fixed()
            rsimp[2].establish_fixed()
            block_vo_simp = simplify(block_vo)

            block_vv.append(rsimp[3])
            rsimp[3].clear_fixed()
            rsimp[3].establish_fixed()
            block_vv_simp = simplify(block_vv)

        test.append(rsimp[x])

    print('oo', len(block_oo_simp))
    print('ov', len(block_ov_simp))
    print('vo', len(block_vo_simp))
    print('vv', len(block_vv_simp))

    return block_oo_simp, block_ov_simp, block_vo_simp, block_vv_simp


def fill_dm_blocks_wm(rsimp):
    print(len(rsimp))

    test = arithmetic_string()

    # block_oo_diag = arithmetic_string()
    # block_vv_diag = arithmetic_string()
    block_oo = arithmetic_string()
    block_ov = arithmetic_string()
    block_vo = arithmetic_string()
    block_vv = arithmetic_string()
    block = ''

    # \sum_{pq}x_{pq}D_{pq} 
    # D_{pq} is divied into four blocks oo, ov, vo and vv. x coefficient is removed 
    # from all the terms in rsimp. Indices of x are replaced by pq in all cases. 

    k = 0
    pqlist = ['p', 'q']
#    print(len(rsimp))
    for x in range(0,len(rsimp)):
        rsimp_temp = deepcopy(rsimp[x])
        print('---RSIMP TEMP-----')
        print(rsimp_temp)
        for i in range(0, len(rsimp[x].coefficient)):
            if(rsimp[x].coefficient[i]==OBSERVABLE_X or rsimp[x].coefficient[i]==NONSYM_MTRX or rsimp[x].coefficient[i]==OBSERVABLE_X_ASYM):
                obs_idx = deepcopy(rsimp[x].coefficient_idx[i])
                block = find_block_wm(obs_idx)
                # print(block, 'block')
                sm = []
                for j in rsimp[x].summation:
                    if j not in rsimp[x].coefficient_idx[i]:
                        sm.append(j)
                rsimp[x].summation = sm                
                rsimp[x].remove_coefficient_idx(i)

        if block == 'vv':
            if len(rsimp_temp.summation) >=8:
                print('pzed remove', rsimp_temp)
        rsimp[x].remove_coefficient(OBSERVABLE_X)
        rsimp[x].remove_coefficient(NONSYM_MTRX)
        rsimp[x].remove_coefficient(OBSERVABLE_X_ASYM)
        print('POOOO', rsimp[x])
        if block == 'oo':
            rsimp[x].multisubst(obs_idx, pqlist)
            k += 1
            block_oo.append(rsimp[x])
        elif block == 'ov':
            rsimp[x].multisubst(obs_idx, pqlist)
            k += 1
            block_ov.append(rsimp[x])
        elif block == 'vo':
            rsimp[x].multisubst(obs_idx, pqlist)
            k += 1
            block_vo.append(rsimp[x])
        elif block == 'vv':
            rsimp[x].multisubst(obs_idx, pqlist)
            k += 1
            block_vv.append(rsimp[x])
        else:
            print(block)
            print('ten wyraz nie nalezy do zadnego bloku')
            print(rsimp[x])
            sys.exit(0)
    # print('simplify oo')
    # for i in block_oo:
    #     print(i)
    # print('')
    block_oo = simplify(block_oo)
    # print('')
    # for i in block_oo:
    #     print(i)


#    print('simplify ov')
    block_ov = simplify(block_ov)
 #   print('simplify vo')
    block_vo = simplify(block_vo)
   # print('simplify vv')
  #  print('blok vv przed simplify')
   # for x in block_vv:
    #    if len(x.summation) >= 6:
     #       print(x)
    block_vv = simplify(block_vv)
   # print('blok vv po simplify')
    #for x in block_vv:
#        if len(x.summation) >= 6:
     #   print(x)

    # print('blockoo')
    # for x in range(0, len(block_oo)):
    #     print(x, block_oo[x])

    print('oo', len(block_oo))
    print('ov', len(block_ov))
    print('vo', len(block_vo))
    print('vv', len(block_vv))

    # return block_oo_diag, block_vv_diag, block_oo, block_ov, block_vo, block_vv
    return block_oo, block_ov, block_vo, block_vv

def find_block(lst):

    if(lst[0] in occupied and lst[1] in occupied):
        if(lst[0] == lst[1]):
            block = 'oo_diag'
        else:
            block = 'oo'
    elif(lst[0] in occupied and lst[1] in virtual):
        block = 'ov'
    elif(lst[0] in virtual and lst[1] in occupied):
        block = 'vo'
    elif(lst[0] in virtual and lst[1] in virtual):
        if(lst[0] == lst[1]):
            block = 'vv_diag'
        else:
            block = 'vv'
    elif(lst[0] in general):
        if(lst[1] in occupied):
            block = 'po'
        elif(lst[1] in virtual):
            block = 'pv'
        elif(lst[1] in general):
            block = 'pp'
    elif(lst[1] in general):
        if(lst[0] in occupied):
            block = 'op'
        elif(lst[0] in virtual):
            block = 'vp'

    return block

def find_block_wm(lst):

    if(lst[0] in occupied and lst[1] in occupied):
        block = 'oo'
    elif(lst[0] in occupied and lst[1] in virtual):
        block = 'ov'
    elif(lst[0] in virtual and lst[1] in occupied):
        block = 'vo'
    elif(lst[0] in virtual and lst[1] in virtual):
        block = 'vv'
    elif(lst[0] in general):
        if(lst[1] in occupied):
            block = 'po'
        elif(lst[1] in virtual):
            block = 'pv'
        elif(lst[1] in general):
            block = 'pp'
    elif(lst[1] in general):
        if(lst[0] in occupied):
            block = 'op'
        elif(lst[0] in virtual):
            block = 'vp'

    return block


def compute_gamma_dm(method):

    if method == 'ccsd':
        r1 = (arithmetic_string(obs)      
              + evaluate(s1c, obs)        
              + evaluate(s2c, obs)        
              + evaluate(obs, t1, s2c).scale(-1.0) 
              + evaluate(obs, t2, s2c).scale(-1.0)     
              # + evaluate(obs, t2, s2c, s2c).scale(-1.0)     
              )
    elif method == 'cc3':
        r1 = (arithmetic_string(obs)      
              + evaluate(s1c, obs)        
              +  evaluate(s2c, obs)
              + evaluate(obs, t1, s2c).scale(-1.0) 
              + evaluate(obs, t2, s2c).scale(-1.0)     
              + evaluate(obs, t2, t3c).scale(-1.0)     
              # + evaluate(obs, t2, s2c, s2c).scale(-1.0)     
              )

    if method == 'ccsd':
        r22 = ( evaluate(s2c, obs)
                + evaluate(obs, s1c, s2c).scale(-1.0)
                + evaluate(obs, t1, s2c).scale(-1.0)
                # + evaluate(obs, t2, s2c, s2c).scale(-0.5)
                )
    elif method == 'cc3':
        r22 = ( evaluate(s2c, obs)
                + evaluate(t3c, obs)
                + evaluate(obs, s1c, s2c).scale(-1.0)
                + evaluate(obs, t1, s2c).scale(-1.0)
                + evaluate(obs, t2, t3c).scale(-1.0)
                # + evaluate(obs, t2, s2c, s2c).scale(-0.5)
               )

    r21 = (arithmetic_string(obs)
           + evaluate(s2c, obs)
           )
    
    r32 = (evaluate(t3c, obs)
           + evaluate(obs, s1c, s2c).scale(-1.0)
           + evaluate(obs, s2c, s2c).scale(-0.5)
           # + evaluate(obs, t2, s2c, s2c).scale(-0.5)
           )

    r31a = (arithmetic_string(obs)
           + evaluate(s1c, obs)
           + evaluate(s2c, obs)
           # + evaluate(obs, t2, s2c).scale(-1.0)
           )

    r31b = evaluate(s2c, obs)
    
    return r1, r21, r22, r31a, r31b

def compute_xi_dm(method):

    r1 = arithmetic_string(obs) + evaluate(obs, t1) + evaluate(obs, t2)

    if method == 'ccsd':
        r2 = evaluate(obs, t2) + evaluate(obs, t1, t2)
        r3 =  ugg()
        r3.num_factor = 0.0
    elif method == 'cc3':
        r2 = evaluate(obs, t2) + evaluate(obs, t3) + evaluate(obs, t1, t2)
        r3 = evaluate(obs, t3) + evaluate(obs, t2, t2).scale(0.5) + evaluate(obs, t1, t2)
    else:
        print('No such method')
        sys.exit(0)
        

    return r1, r2, r3
   
# ----------------------------------Transition  moments between excited states-------------------------------------

def stoo(nr, st, c, idx = []):

    # translate number nr to the operator
    # st can be s or t
    # c = True if operator should be complex conjugate

    if st == "tf":
        return t2fa

    if st == "sf":
        return s2fa

    if st == "t":
        if nr == 1:
            if c == True:
                return t1c
            else:
                return t1
        if nr == 2:
            if c == True:
                return t2c
            else:
                return t2
        if nr == 3:
            if c == True:
                return t3c
            else:
                return t3
    if st == "s":
        if nr == 1:
            if c == True:
                return s1c
            else:
                return s1
        if nr == 2:
            if c == True:
                return s2c
            else:
                return s2
        if nr == 3:
            if c == True:
                return s3c
            else:
                return s3
    if st == "sf12":
        print('nr', nr, c)
        if nr == 1:
            if c == True:
            #    return s1Ac
                return s1fac
            else:
                return s1fa
#                return s1A
        if nr == 2:
            if c == True:
                return t2fac
            else:
                return t2fa

            
            # if idx == 1:
            #     if c == True:
            #         return s2Aac
            #     else:
            #         return s2Aa
            # elif idx == 2:
            #     if c == True:
            #         return s2bBc
            #     else:
            #         return s2bB
            # elif idx == 3:
            #     if c == True:
            #         return s2ABc
            #     else:
            #         return s2AB

    return False


def stoo_f12(nr, st, c, idx = []):

    # translate number nr to the operator
    # st can be s or t
    # c = True if operator should be complex conjugate


    if st == "s_v":
        return s_v1, s_v2 

    if st == "s_c":
        return s_c1, s_c2

    if st == "s_A":
        return s_A1, s_A2

    if st == "s":
        return t2c        

    if st == "s_cc":
        return t2fac

    if st == "s_vv":
        return t2c

    return False

def generate_commutators_overlap(maxmbpt, theory, cumulative):
    
    Wm2 = generate_Wm2(maxmbpt, theory)
    Wm3 = generate_Wm3(maxmbpt, theory)

    print('Wm2')
    for x in Wm2:
        print(x)

    print('Wm3')
    for x in Wm3:
        print(x)


    Wm_overlap = generate_overlap(Wm2, Wm3, maxmbpt, cumulative)
    
    return Wm_overlap


def generate_commutators_for_gamma(maxmbpt, maxcluster, theory, cumulative):

    print('GENERATING WM3')
    Wm3t = generate_Wm3(maxmbpt, theory)
    print('GENERATING WM1')
    Wm1t = generate_Wm1(maxmbpt, theory)

    Wm3 = []
    Wm1 = []

    for x in Wm3t:
        print(x['mbpt'], x)
        if x['mbpt'] <=maxmbpt:
            Wm3.append(x)
    for x in Wm1t:
        print(x['mbpt'], x)
        if x['mbpt'] <=maxmbpt:
            Wm1.append(x)

    print('Wm3')
    for x in Wm3:
        print(x)
    print('Wm1')
    for x in Wm1:
        print(x)

    print('GENERATING WM')
    Wg = generate_Wg(Wm1, Wm3, maxmbpt, maxcluster, cumulative)


    for x in Wg:
        print('aa', x)
    print('GENERATING W_MIDDLE')
    W_gamma = generate_W_gamma(Wg)
    
    
    return W_gamma

def generate_commutators_for_cumulant(maxmbpt, maxcluster, theory, cumulative):

    # Gm =  <e(S*)e(-T) E_pq e(T) e(-S*) P(e(S*)e(-T) E_rs e(T) e(-S*))>                                                         
        #--------Gm1----------#  #-------------Gm2-------------#  

    print('GENERATING GM2')
    Gm2t = generate_Gm2(maxmbpt, theory)
    print('GENERATING GM1')
    Gm3t = generate_Gm1(maxmbpt, theory)



    Gm2 = []
    Gm1 = []

    for x in Wm2t:
        print(x['mbpt'], x)
        if x['mbpt'] <= maxmbpt:
            Wm2.append(x)
    for x in Wm3t:
        print(x['mbpt'], x)
        if x['mbpt'] <=maxmbpt:
            Wm3.append(x)
    for x in Wm1t:
        print(x['mbpt'], x)
        if x['mbpt'] <=maxmbpt:
            Wm1.append(x)

    print('Wm2')
    for x in Wm2:
        print(x)
    print('Wm3')
    for x in Wm3:
        print(x)
    print('Wm1')
    for x in Wm1:
        print(x)
        
    print('GENERATING WM')
    Wm, Wm_nd = generate_Wm(Wm1, Wm2, Wm3, maxmbpt, maxcluster, cumulative)
    Um = generate_Um(Wm1, Wm2, Wm3, maxmbpt)
    Wl, Wr = generate_WlWr(maxmbpt, theory)


    for x in Wm:
        print(x)



    print('GENERATING W_MIDDLE')
    W_middle = generate_W_middle(Wm)
    
#    W_left, W_middle, W_right = generate_W_left_W_middle_W_right(Wl, Wm, Wr, maxmbpt)
 #   U_middle = generate_U_middle(Wl, Um, Wr, maxmbpt)
    W_left = []
    W_right = []
    U_middle = []
    
    return W_left, W_middle, W_right, U_middle



def generate_commutators_for_trans_exc(maxmbpt, maxcluster, theory, cumulative):

    print('GENERATING WM2')
    Wm2t = generate_Wm2(maxmbpt, theory)
    print('GENERATING WM3')
    Wm3t = generate_Wm3(maxmbpt, theory)
    print('GENERATING WM1')
    Wm1t = generate_Wm1(maxmbpt, theory)



    Wm2 = []
    Wm3 = []
    Wm1 = []

    for x in Wm2t:
        print(x['mbpt'], x)
        if x['mbpt'] <= maxmbpt:
            Wm2.append(x)
    for x in Wm3t:
        print(x['mbpt'], x)
        if x['mbpt'] <=maxmbpt:
            Wm3.append(x)
    for x in Wm1t:
        print(x['mbpt'], x)
        if x['mbpt'] <=maxmbpt:
            Wm1.append(x)

    print('Wm2')
    for x in Wm2:
        print(x)
    print('Wm3')
    for x in Wm3:
        print(x)
    print('Wm1')
    for x in Wm1:
        print(x)
        
    print('GENERATING WM')
    Wm, Wm_nd = generate_Wm(Wm1, Wm2, Wm3, maxmbpt, maxcluster, cumulative)
    Um = generate_Um(Wm1, Wm2, Wm3, maxmbpt)
    Wl, Wr = generate_WlWr(maxmbpt, theory)


    for x in Wm:
        print(x)



    print('GENERATING W_MIDDLE')
    W_middle = generate_W_middle(Wm)
    
#    W_left, W_middle, W_right = generate_W_left_W_middle_W_right(Wl, Wm, Wr, maxmbpt)
 #   U_middle = generate_U_middle(Wl, Um, Wr, maxmbpt)
    W_left = []
    W_right = []
    U_middle = []
    
    return W_left, W_middle, W_right, U_middle

def commutators_quadra_Wlr(W_left, W_right):

    mu_list = []
    Wl_commutators = []
    Wr_commutators = []

    for x in range(0, len(W_left)):
        
        wt = W_left[x]['T_list']

        ls = len(wt)

        if (W_left[x]['n']) == 1:
            mu_list.append(1)
        elif (W_left[x]['n']) == 2:
            mu_list.append(2)
        elif (W_left[x]['n']) == 3:
            mu_list.append(3)

       # Wleft
        w = arithmetic_string(obs)
        if ls !=0:
            for i in range(0, ls):
                w = evaluate(w, stoo(wt[i], "t", False))


        Wl_commutators.append(w)
        Wr_commutators.append(w)


    return Wl_commutators, Wr_commutators, mu_list

def commutators_quadra_overlap(Wm_overlap):
    Wm_overlap_commutators = []
    Wm2_commutators = []
    Wm3_commutators = []
    Wm2_mu_list = []
    Wm3_mu_list = []

    for x in range(0, len(Wm_overlap)):

        wm2t = Wm_overlap[x]['Wm2_T_list']
        wm2s = Wm_overlap[x]['Wm2_S_list']

        wm3s = Wm_overlap[x]['Wm3_S_list']

        n_wm2t = len(set(permutations(wm2t, len(wm2t))))
        n_wm2s = len(set(permutations(wm2s, len(wm2s))))
        n_wm3s = len(set(permutations(wm3s, len(wm3s))))


        lt2 = len(wm2t)
        ls2 = len(wm2s)

        ls3 = len(wm3s)

        wm2_mu = ugg()
        if (Wm_overlap[x]['Wm2_n']) == 1:
#            wm2_mu = mu1
            Wm2_mu_list.append(1)
            eoml = eomrl1
        elif (Wm_overlap[x]['Wm2_n']) == 2:
#            wm2_mu = mu2
            Wm2_mu_list.append(2)
            eoml = eomrl2
        elif (Wm_overlap[x]['Wm2_n']) == 3:
#            wm2_mu = mu3
            Wm2_mu_list.append(3)
            eoml = eomrl3

        wm3_mu = ugg()
        eomr = ugg()
        if (Wm_overlap[x]['Wm3_n']) == 1:
            # wm3_mu = mu1p
            Wm3_mu_list.append(1)
            eomr = eomrr1
        elif (Wm_overlap[x]['Wm3_n']) == 2:
            # wm3_mu = mu2p
            Wm3_mu_list.append(2)
            eomr = eomrr2
        elif (Wm_overlap[x]['Wm3_n']) == 3:
            eomr = eomrr3
            Wm3_mu_list.append(3)

        # Wm2
        if lt2 == 0:
#            wm2_inner = arithmetic_string(wm2_mu)
            wm2_inner = arithmetic_string(eoml)
        else:
#            wm2_inner = arithmetic_string(wm2_mu)
            wm2_inner = arithmetic_string(eoml)
            for i in range(0, lt2):
                wm2_inner = evaluate(wm2_inner, stoo(wm2t[i], "t", True))

        wm2_outer = wm2_inner
        if ls2 !=0:
            for i in range(0, ls2):
                wm2_outer = evaluate(wm2_outer, stoo(wm2s[i], "s", False))

            wm2_outer = wm2_outer.scale(-1.0)
            wm2_outer = wm2_outer.scale(n_wm2t)
            wm2_outer =wm2_outer.scale(n_wm2s)
        Wm2_commutators.append(wm2_outer)

        # Wm3
        wm3_outer = arithmetic_string(eomr)
        # wm3_outer = arithmetic_string(wm3_mu)
        if ls3 !=0:
            for i in range(0, ls3):
                wm3_outer = evaluate(wm3_outer, stoo(wm3s[i], "s", True))

            wm3_outer = wm3_outer.scale(-1.0)
            wm3_outer =wm3_outer.scale(n_wm3s)
        Wm3_commutators.append(wm3_outer)



    return Wm2_commutators,  Wm3_commutators, Wm2_mu_list, Wm3_mu_list


def commutators_quadra_overlap_triplet(Wm_overlap):

    Wm_overlap_commutators = []
    Wm2_commutators = []
    Wm3_commutators = []
    Wm2_mu_list = []
    Wm3_mu_list = []

    for x in range(0, len(Wm_overlap)):

        wm2t = Wm_overlap[x]['Wm2_T_list']
        wm2s = Wm_overlap[x]['Wm2_S_list']

        wm3s = Wm_overlap[x]['Wm3_S_list']


        n_wm2t = len(set(permutations(wm2t, len(wm2t))))
        n_wm2s = len(set(permutations(wm2s, len(wm2s))))

        n_wm3s = len(set(permutations(wm3s, len(wm3s))))

    

        lt2 = len(wm2t)
        ls2 = len(wm2s)

        ls3 = len(wm3s)

        if (Wm_overlap[x]['Wm2_n']) == 1 and (Wm_overlap[x]['Wm3_n']) == 1:
            wm2_mu = eomrl1_triplet
            eomr = eomrr1_triplet
            Wm2_mu_list.append(1)
            Wm3_mu_list.append(1)

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (Wm_overlap[x]['Wm2_n']) == 1 and (Wm_overlap[x]['Wm3_n']) == 2:
            wm2_mu = eomrl1_triplet
            eomr_plus = eomrr2_triplet_plus
            eomr_minus = eomrr2_triplet_minus
            Wm2_mu_list.append(1)
            Wm2_mu_list.append(1)
            Wm3_mu_list.append('2p')
            Wm3_mu_list.append('2m')


            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr_plus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            wm3_outer = wm3_inner_outer(eomr_minus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (Wm_overlap[x]['Wm2_n']) == 2 and (Wm_overlap[x]['Wm3_n']) == 1:
            wm2_mu_plus = eomrl2_triplet_plus
            wm2_mu_minus = eomrl2_triplet_minus
            eomr = eomrr1_triplet
            Wm2_mu_list.append('2p')
            Wm2_mu_list.append('2m')
            Wm3_mu_list.append(1)
            Wm3_mu_list.append(1)

            wm2_outer = wm2_inner_outer(wm2_mu_plus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            wm2_outer = wm2_inner_outer(wm2_mu_minus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (Wm_overlap[x]['Wm2_n']) == 2 and (Wm_overlap[x]['Wm3_n']) == 2:
            wm2_mu_plus = eomrl2_triplet_plus
            wm2_mu_minus = eomrl2_triplet_minus
            eomr_plus = eomrr2_triplet_plus
            eomr_minus = eomrr2_triplet_minus
            Wm2_mu_list.append('2p')
            Wm2_mu_list.append('2p')
            Wm2_mu_list.append('2m')
            Wm2_mu_list.append('2m')

            Wm3_mu_list.append('2p')
            Wm3_mu_list.append('2p')
            Wm3_mu_list.append('2m')
            Wm3_mu_list.append('2m')

            wm2_outer = wm2_inner_outer(wm2_mu_plus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            wm2_outer = wm2_inner_outer(wm2_mu_plus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            wm2_outer = wm2_inner_outer(wm2_mu_minus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            wm2_outer = wm2_inner_outer(wm2_mu_minus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr_plus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            wm3_outer = wm3_inner_outer(eomr_minus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            wm3_outer = wm3_inner_outer(eomr_plus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            wm3_outer = wm3_inner_outer(eomr_minus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))


    return Wm2_commutators,  Wm3_commutators, Wm2_mu_list, Wm3_mu_list


def look_for_t(wm1t, wm1s, wm2t, wm2s, wm3s, t):

    add = True

    for y in wm1t:
        if y == t:
            add = False
    for y in wm1s:
        if y== t:
            add = False
    for y in wm2t:
        if y== t:
            add = False
    for y in wm2s:
        if y== t:
            add = False
    for y in wm3s:
        if y== t:
            add = False

    if t == 1:
        add = not add
    
    return add

def cut_commutators_Wm(W_middle, maxpt):

    W_middle_out = []
    len(W_middle)

    for x in range(0, len(W_middle)):
        print(x, W_middle[x])
        wm1t = W_middle[x]['Wm1_T_list']
        wm1s = W_middle[x]['Wm1_S_list']

        wm2t = W_middle[x]['Wm2_T_list']
        wm2s = W_middle[x]['Wm2_S_list']

        wm3s = W_middle[x]['Wm3_S_list']
        pt = W_middle[x]['mbpt']

        n_wm1t = len(set(permutations(wm1t, len(wm1t))))
        n_wm1s = len(set(permutations(wm1s, len(wm1s))))

        n_wm2t = len(set(permutations(wm2t, len(wm2t))))
        n_wm2s = len(set(permutations(wm2s, len(wm2s))))

        n_wm3s = len(set(permutations(wm3s, len(wm3s))))

        add = True

        if (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 1:
            if pt == 4:
                add = look_for_t(wm1t, wm1s, wm2t, wm2s, wm3s, 1)
            else:
                add = True
        elif (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 2:
            if pt == 4:
                add = look_for_t(wm1t, wm1s, wm2t, wm2s, wm3s, 1)
            else:
                add = True
        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 1:
            if pt == 4:
                add = look_for_t(wm1t, wm1s, wm2t, wm2s, wm3s, 1)
            else:
                add = True
        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 2:
            if pt == 4:
                add = look_for_t(wm1t, wm1s, wm2t, wm2s, wm3s, 1)
            else:
                add = True

        elif (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 3:
            if (pt == 3 or pt == 2 or pt == 4):
                add = False
            else:
                add = True
        elif (W_middle[x]['Wm2_n'] + W_middle[x]['Wm3_n']) == 5:
            if (pt == 4):
                add = False
                #if (pt == 3 or pt == 2 or pt ==4):
            # if (pt == 2):
            #     add = True
                
            # if (pt == 3):
            #     add = True

            if (pt == 3) or (pt == 2):
                add1 = look_for_t(wm1t, wm1s, wm2t, wm2s, wm3s, 1)
                add2 = look_for_t(wm1t, wm1s, wm2t, wm2s, wm3s, 3)
                if add1 == True and add2 == True:
                    add = True
                else:
                    add = False
            #     exc = sum(wm1t) + sum(wm1s) + sum(wm2t) + sum(wm2s) + sum(wm3s)
                # if exc == 6:
                #     # print('excccc 6')
                #     # sys.exit(0)
                #     add = False

        elif (W_middle[x]['Wm2_n']==1 and W_middle[x]['Wm3_n'] == 3) or (W_middle[x]['Wm2_n']==3 and W_middle[x]['Wm3_n'] == 1):
            if (pt == 4):
                add = False
            if (pt == 2):
                add = True
            if (pt == 3):
                add1 = look_for_t(wm1t, wm1s, wm2t, wm2s, wm3s, 1)
                add2 = look_for_t(wm1t, wm1s, wm2t, wm2s, wm3s, 3)
                if add1 == True and add2 == True:
                    add = True
                else:
                    add = False
                # exc = sum(wm1t) + sum(wm1s) + sum(wm2t) + sum(wm2s) + sum(wm3s)
                # if exc == 6:
                #     print('excccc 6')
                #     sys.exit(0)
                #     add = False

           

        if add == True:
            print(x, 'dodaje')
            W_middle_out.append(W_middle[x])
        else:
            print(x, 'nie dodaje')
    print(len(W_middle_out))


    return W_middle_out



def commutators_Gm(Gm):

    Gm_commutators = []
    Gm1_commutators = []
    Gm2_commutators = []

    for x in range(0, len(Gm)):

        wm1t = Gm[x]['Gm1_T_list']
        wm1s = Gm[x]['Gm1_S_list']

        wm2t = Gm[x]['Gm2_T_list']
        wm2s = Gm[x]['Gm2_S_list']

        wm2p = Gm[x]['Gm2_P_list']

        n_wm1t = len(set(permutations(wm1t, len(wm1t))))
        n_wm1s = len(set(permutations(wm1s, len(wm1s))))

        n_wm2t = len(set(permutations(wm2t, len(wm2t))))
        n_wm2s = len(set(permutations(wm2s, len(wm2s))))

        lt1 = len(wm1t)
        ls1 = len(wm1s)

        lt2 = len(wm2t)
        ls2 = len(wm2s)


        if (Gm[x]['Gm2_n']) == 1 and (Gm[x]['Gm3_n']) == 1:
            Gm_out.append(Gm[x])
            #wm2_mu = mu1
            wm2_mu = eomrl1
            eomr = eomrr1
            Gm2_mu_list.append(1)
            Gm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s) 
            Gm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Gm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Gm3_commutators.append(deepcopy(wm3_outer))

        elif (Gm[x]['Gm2_n']) == 3 and (Gm[x]['Gm3_n']) == 3:

            Gm_out.append(Gm[x])
            #wm2_mu = mu3
            wm2_mu = eomrl3
            eomr = eomrr3
            Gm2_mu_list.append(3)
            Gm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Gm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Gm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Gm3_commutators.append(deepcopy(wm3_outer))

        elif (Gm[x]['Gm2_n']) == 1 and (Gm[x]['Gm3_n']) == 3:
            Gm_out.append(Gm[x])
            #wm2_mu = mu1
            wm2_mu = eomrl1
            eomr = eomrr3
            Gm2_mu_list.append(1)
            Gm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Gm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Gm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Gm3_commutators.append(deepcopy(wm3_outer))

        elif (Gm[x]['Gm2_n']) == 3 and (Gm[x]['Gm3_n']) == 1:
            Gm_out.append(Gm[x])
            #wm2_mu = mu3
            wm2_mu = eomrl3
            eomr = eomrr1
            Gm2_mu_list.append(3)
            Gm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Gm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Gm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Gm3_commutators.append(deepcopy(wm3_outer))

        elif (Gm[x]['Gm2_n']) == 1 and (Gm[x]['Gm3_n']) == 2:
            Gm_out.append(Gm[x])
            #wm2_mu = mu1
            wm2_mu = eomrl1
            eomr = eomrr2
            Gm2_mu_list.append(1)
            Gm3_mu_list.append(2)


            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Gm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Gm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Gm3_commutators.append(deepcopy(wm3_outer))

        elif (Gm[x]['Gm2_n']) == 3 and (Gm[x]['Gm3_n']) == 2:

            Gm_out.append(Gm[x])

            #wm2_mu = mu3
            wm2_mu = eomrl3
            eomr = eomrr2
            Gm2_mu_list.append(3)
            Gm3_mu_list.append(2)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Gm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Gm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Gm3_commutators.append(deepcopy(wm3_outer))


            for er in Gm3_commutators:
                print(er)
        elif (Gm[x]['Gm2_n']) == 2 and (Gm[x]['Gm3_n']) == 1:
            Gm_out.append(Gm[x])
            #wm2_mu = mu2
            wm2_mu = eomrl2
            eomr = eomrr1
            Gm2_mu_list.append(2)
            Gm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Gm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Gm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Gm3_commutators.append(deepcopy(wm3_outer))

        elif (Gm[x]['Gm2_n']) == 2 and (Gm[x]['Gm3_n']) == 3:
            Gm_out.append(Gm[x])
            # disambiguate(eomrr2, eomrr3)
            #wm2_mu = mu2
            wm2_mu = eomrl2
            eomr = eomrr3
            Gm2_mu_list.append(2)
            Gm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Gm1_commutators.append(deepcopy(wm1_outer))
            # for x in Gm1_commutators:
            #     print('wm1', x)
                
            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Gm2_commutators.append(deepcopy(wm2_outer))
            # for x in Gm2_commutators:
            #     print('wm2', x)

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Gm3_commutators.append(deepcopy(wm3_outer))
            # for x in Gm3_commutators:
            #     for y in x:
            #         print('wm3', y)
            #     print('--')

            
        elif (Gm[x]['Gm2_n']) == 2 and (Gm[x]['Gm3_n']) == 2:

            Gm_out.append(Gm[x])
            #wm2_mu = mu2
            wm2_mu = eomrl2
            eomr = eomrr2
            Gm2_mu_list.append(2)
            Gm3_mu_list.append(2)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Gm1_commutators.append(deepcopy(wm1_outer))

            
            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Gm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Gm3_commutators.append(deepcopy(wm3_outer))


    print(len( Gm1_commutators), len(Gm2_commutators), len(Gm3_commutators), len(Gm2_mu_list),  len(Gm3_mu_list))


    return Gm_out, Gm1_commutators,  Gm2_commutators,  Gm3_commutators, Gm2_mu_list, Gm3_mu_list



def commutators_quadra_Wm(W_middle):

    W_middle_commutators = []
    Wm1_commutators = []
    Wm2_commutators = []
    Wm3_commutators = []
    Wm2_mu_list = []
    Wm3_mu_list = []
    W_middle_out = []

    for x in range(0, len(W_middle)):

        wm1t = W_middle[x]['Wm1_T_list']
        wm1s = W_middle[x]['Wm1_S_list']

        wm2t = W_middle[x]['Wm2_T_list']
        wm2s = W_middle[x]['Wm2_S_list']

        wm3s = W_middle[x]['Wm3_S_list']

        n_wm1t = len(set(permutations(wm1t, len(wm1t))))
        n_wm1s = len(set(permutations(wm1s, len(wm1s))))

        n_wm2t = len(set(permutations(wm2t, len(wm2t))))
        n_wm2s = len(set(permutations(wm2s, len(wm2s))))

        n_wm3s = len(set(permutations(wm3s, len(wm3s))))

    
        lt1 = len(wm1t)
        ls1 = len(wm1s)

        lt2 = len(wm2t)
        ls2 = len(wm2s)

        ls3 = len(wm3s)



        if (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 1:
            W_middle_out.append(W_middle[x])
            #wm2_mu = mu1
            wm2_mu = eomrl1
            eomr = eomrr1
            Wm2_mu_list.append(1)
            Wm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s) 
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 3:

            W_middle_out.append(W_middle[x])
            #wm2_mu = mu3
            wm2_mu = eomrl3
            eomr = eomrr3
            Wm2_mu_list.append(3)
            Wm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 3:
            W_middle_out.append(W_middle[x])
            #wm2_mu = mu1
            wm2_mu = eomrl1
            eomr = eomrr3
            Wm2_mu_list.append(1)
            Wm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 1:
            W_middle_out.append(W_middle[x])
            #wm2_mu = mu3
            wm2_mu = eomrl3
            eomr = eomrr1
            Wm2_mu_list.append(3)
            Wm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 2:
            W_middle_out.append(W_middle[x])
            #wm2_mu = mu1
            wm2_mu = eomrl1
            eomr = eomrr2
            Wm2_mu_list.append(1)
            Wm3_mu_list.append(2)


            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 2:

            W_middle_out.append(W_middle[x])

            #wm2_mu = mu3
            wm2_mu = eomrl3
            eomr = eomrr2
            Wm2_mu_list.append(3)
            Wm3_mu_list.append(2)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))


            for er in Wm3_commutators:
                print(er)
        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 1:
            W_middle_out.append(W_middle[x])
            #wm2_mu = mu2
            wm2_mu = eomrl2
            eomr = eomrr1
            Wm2_mu_list.append(2)
            Wm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 3:
            W_middle_out.append(W_middle[x])
            # disambiguate(eomrr2, eomrr3)
            #wm2_mu = mu2
            wm2_mu = eomrl2
            eomr = eomrr3
            Wm2_mu_list.append(2)
            Wm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))
            # for x in Wm1_commutators:
            #     print('wm1', x)
                
            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            # for x in Wm2_commutators:
            #     print('wm2', x)

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            # for x in Wm3_commutators:
            #     for y in x:
            #         print('wm3', y)
            #     print('--')

            
        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 2:

            W_middle_out.append(W_middle[x])
            #wm2_mu = mu2
            wm2_mu = eomrl2
            eomr = eomrr2
            Wm2_mu_list.append(2)
            Wm3_mu_list.append(2)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            
            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))


    print(len( Wm1_commutators), len(Wm2_commutators), len(Wm3_commutators), len(Wm2_mu_list),  len(Wm3_mu_list))


    return W_middle_out, Wm1_commutators,  Wm2_commutators,  Wm3_commutators, Wm2_mu_list, Wm3_mu_list


def generate_t1_in_2_mbpt():

    r = evaluate(flukt_potential, t2)

    rint = r.integrate(bra = ['a', 'i'], braspin =['s']).scale(0.5)
    rint.exec_delta()
    
    rsimp = simplify(rint)
    rsimp.cleanup()
    for x in rsimp:
        print(x)

    for x in rsimp:
        x.operator_idx.append(['a', 'i'])
        x.operator_type.append('s')
        print('xxx', x.operator_idx)



    print('')
    print('to wynik')
    print('')
    for x in rsimp:
        print(x)
    sys.exit(0)

def generate_t1_in_3_mbpt():

    r = evaluate(flukt_potential, t2)  +  evaluate(flukt_potential, t1, t2) \
        +  evaluate(flukt_potential, t2, t2).scale(0.5)

    rint = r.integrate(bra = ['a', 'i'], braspin =['s']).scale(0.5)
    rint.exec_delta()
    
    rsimp = simplify(rint)
    rsimp.cleanup()
    for x in rsimp:
        print(x)

    for x in rsimp:
        x.operator_idx.append(['a', 'i'])
        x.operator_type.append('s')
        print('xxx', x.operator_idx)



    print('')
    print('to wynik')
    print('')
    for x in rsimp:
        print(x)
    sys.exit(0)


def generate_t2_in_2_mbpt():

    r = evaluate(flukt_potential, t2)

    # rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin =['s', 's'])
    # rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], braspin =['s', 's'])
    # rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
    rint = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin =['s', 's'])
    rint.exec_delta()
    
    rsimp = simplify(rint)
    rsimp.cleanup()
    for x in rsimp:
        print(x)

    for x in rsimp:
        x.num_factor *= 1./2.
        x.operator_idx.append(['a', 'i'])
        x.operator_idx.append(['b', 'j'])
        x.operator_type.append('s')
        x.operator_type.append('s')
        print('xxx', x.operator_idx)



    print('')
    print('to wynik')
    print('')
    for x in rsimp:
        print(x)
    sys.exit(0)


def generate_t2_in_3_mbpt():

    r = evaluate(flukt_potential, t2) +  evaluate(flukt_potential, t1, t1).scale(0.5)+  evaluate(flukt_potential, t1, t2) \
        +  evaluate(flukt_potential, t2, t2).scale(0.5)
        
    
    # rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin =['s', 's'])
    # rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], braspin =['s', 's'])
    # rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
    rint = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin =['s', 's'])
    rint.exec_delta()
    
    rsimp = simplify(rint)
    rsimp.cleanup()
    for x in rsimp:
        print(x)

    for x in rsimp:
        x.num_factor *= 1./2.
        x.operator_idx.append(['a', 'i'])
        x.operator_idx.append(['b', 'j'])
        x.operator_type.append('s')
        x.operator_type.append('s')
        print('xxx', x.operator_idx)



    print('')
    print('to wynik')
    print('')
    for x in rsimp:
        print(x)
    sys.exit(0)


def generate_t2_in_1_mbpt():

    r = arithmetic_string(flukt_potential)

    rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin =['s', 's'])
    rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], braspin =['s', 's'])
    rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
    rint.exec_delta()
    
    rsimp = simplify(rint)
    rsimp.cleanup()
    
    for x in rsimp:
        x.num_factor *= 1./2.
        x.operator_idx.append(['a', 'i'])
        x.operator_idx.append(['b', 'j'])
        x.operator_type.append('s')
        x.operator_type.append('s')
        print('xxx', x.operator_idx)



    print('')
    print('to wynik')
    print('')
    for x in rsimp:
        print(x)
    sys.exit(0)


def generate_density_Wm_ground_driver(pick, method):

    f = open("Density_matrix.tex", 'w')
    s_preamble = tex_preamble
    s_beginning = """
    \\begin{document}"""
    f.write(s_preamble)
    f.write(s_beginning)

    for order in range(2, 3):
        print('teraz bedzie rzad', order)
        s2 = "\n"
        if order > 0:
            s2 += "\\newpage \n"
        s2 += f"Macierz gęstości w rzędzie {order}\\\\"
        s2 += "\n"
        f.write(s2)
        rsimp, W_middle = generate_density_Wm_ground_mbpt(pick, method, order, f, onlyrsimp = True, operator=obs)
        f.write("\n")
        
        f.write("\\hrule")
        f.write("\n")

        
    f.write("\end{document}")
    f.close()

    
def generate_density_Wm_ground_mbpt(pick, method, mbpt, f, onlyrsimp, operator=obs, cumulative = True):

    print('jestem wewnatrz density wmwmwm', cumulative)
    # Tworze liste mozliwych T1 i S1
    if method == 'ccd' or method == 'ccsd':
        Wm1 = generate_Wm_ground_ccd(mbpt, 'ccsd', cumulative=cumulative)
    elif method == 'cc3':
        # print('cumulative-aa', cumulative)
        Wm1 = generate_Wm_ground_ccd(mbpt, 'cc3', cumulative=cumulative)
        # print('FF wm1', Wm1, mbpt)
    W_middle_Wm1 = generate_W_middle_for_Wm1(Wm1)
    for x in W_middle_Wm1:
        print('a', x)
    # W_middle = W_middle_Wm1
    W_middle = []
    #----------------------WERSJA Z T2
    if method == 'ccd':
        for x in W_middle_Wm1:
            if 1 not in x['Wm1_S_list']:
                if 1 not in x['Wm1_T_list']:
                    W_middle.append(x)

    #----------------------WERSJA Z T1 i T2
    if method == 'ccsd' or method == 'cc3':
        for x in W_middle_Wm1:
            W_middle.append(x)

    print('przed latex')
    # licze komutatory i nastpenie to calkuje
    latex_W_middle_Wm1(W_middle, mbpt, False, f)
    # # sys.exit(0)
    # W_middle = [W_middle[1]]

    # UWAGA!!!!!!!!!!!!!!!1 - tu sie podaje jako drugi parametr, jakie amplitudy beda zamieniane, 1=s1 czy 2=s2
    
    W_middle, W_commutators, rsimp = commutators_density_ground_mbpt(W_middle, 2, f, onlyrsimp, operator)
    print('')
    print('rsimp ost w rzad ',mbpt)
    for x in rsimp:
        print(x)
    print('')

    return rsimp, W_middle

def compute_mbpt_for_density_Wm_ground_mbpt(x):

    mbpt = 0
    for i in range(0, len(x.coefficient)):
        if x.coefficient[i] == S_AMPLITUDE or x.coefficient[i] == CC_AMPLITUDE:
            ln, at = x.amplitude_type(i)
            if ln == '1':
                mbpt += 1
            elif ln == '2':
                mbpt += 1

    return mbpt

def commutators_density_ground_mbpt(W_middle, s_op, f, onlyrsimp, operator):

    W_middle_commutators = []
    Wm1_commutators = []
    W_middle_out = []
    print('mbpt')
    W_middle2 = W_middle
    # W_middle2 = []
    # # usuwam wszystkie T1 i S1, bo tak chciał Polnon
    # for x in W_middle:
    #     if not(1 in x['Wm1_T_list'] or 1 in x['Wm1_S_list']):
    #         W_middle2.append(x)
            
    print('wmiddle wmiddle')
    for x in W_middle2:
        print(x)
    print('')

    W_middle = deepcopy(W_middle2)
    # licze komutatory wewnetrzne i zewnetrzne
    for x in range(0, len(W_middle)):
        print('a teraz rozwazam takie cos', x, len(W_middle))
        wm1t = W_middle[x]['Wm1_T_list']
        wm1s = W_middle[x]['Wm1_S_list']

        n_wm1t = len(set(permutations(wm1t, len(wm1t))))
        n_wm1s = len(set(permutations(wm1s, len(wm1s))))

        lt1 = len(wm1t)
        ls1 = len(wm1s)

        W_middle_out.append(W_middle[x])
        
        wm1_outer = density_mbpt_inner_outer(operator, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
        
        Wm1_commutators.append(deepcopy(wm1_outer))

    print('W_middle_out')
    for x in W_middle_out:
        print(x)

        # przekształcam listw w arithmetic string
    density_int = arithmetic_string()
    for i in range(0, len(Wm1_commutators)):
        for j in range(0, len(Wm1_commutators[i])):
            print(Wm1_commutators[i][j])
            density_int.append(Wm1_commutators[i][j])

    # calkowanie
    print('integrate')
    rint = density_int.integrate()
    print('po integrate')
    print('rint')
    for x in rint:
        print(x)

    # simplifikacja
    rsimp = simplify(rint)
    rsimp_for_out = deepcopy(rsimp)

    if onlyrsimp == True:
        print('Return Only Rsimp, do not group, do not change s operators')
        Wm1_commutators = []
        return W_middle_out, Wm1_commutators, rsimp_for_out
    else:
        

        print('po simplify rra')

        for x in rsimp:
            print(x)
            #        print(x.summation)
        print('')
    #    sys.exit(0)
        # generowanie s21 i s23 bo tylko te chcial polnon
        s21, s23, s12, s13 = generate_W_middle_with_explicit_S()

        print('s21')
        for x in s21:
            print(x)

        print('s23')
        for x in s23:
            print(x)

        print('s12')
        for x in s12:
            print(x)

        print('s13')
        for x in s13:
            print(x)

        print('')


        # # Wersja z T1 i T2
        # # tworze kopie, zeby w pierwszej kopii za s wstawic s21  a w drugiej s23
        # rsimp1 = deepcopy(rsimp)

        # # chamska zamiama t2 na s2 i t1 na s1
        # for x in range(0, len(rsimp1)):
        #     for j in range(0, len(rsimp1[x].coefficient)):
        #         coef = rsimp1[x].coefficient[j]
        #         if coef == S_AMPLITUDE:
        #             rsimp1[x].coefficient[j] = CC_AMPLITUDE

        # print('teraz jedzie zamiana S na s13')
        # restr_aibj = ['a', 'i']
        # rsimp2 = arithmetic_string()
        # for x in range(0, len(rsimp)):
        #     print('biore comut', rsimp[x])
        #     for j in range(0, len(rsimp[x].coefficient)):
        #         coef = rsimp[x].coefficient[j]
        #         coef_idx = rsimp[x].coefficient_idx[j]
        #         coef_sum = rsimp[x].summation
        #         restr2 = coef_idx
        #         restr3 = coef_sum

        #         if coef == S_AMPLITUDE:
        #             for k in range(0, len(s13)):
        #                 s23_copy = deepcopy(s13)
        #                 #1 zamien indeksy sumowania tak, żeby nie wolno bylo skorzystac z inkdesow restr2:'cbkl' i abij
        #                 print('1   ', s23_copy[k], restr_aibj, restr2)
        #                 s23_copy[k].rename_summation_with_restriction(restr_aibj, restr2)
        #                 print('2   ',s23_copy[k])
        #                 s23_copy[k].multisubst(restr_aibj, restr2)
        #                 print('3   ',s23_copy[k], restr2, restr3)
        #                 s23_copy[k].rename_summation_with_restriction(restr2, restr3)
        #                 print(i, j, k, rsimp[x], s23_copy[k])
        #                 # tworze wyraz bedacy suma s i rsimp[x]
        #                 ars_temp = deepcopy(rsimp[x])
        #                 print('ars_temp', ars_temp)
        #                 ars_temp.summation += s23_copy[k].summation
        #                 print('ars_temp sum', ars_temp)
        #                 ars_temp.coefficient += s23_copy[k].coefficient
        #                 # print('ars_tempcf', ars_temp)
        #                 ars_temp.coefficient_idx += s23_copy[k].coefficient_idx
        #                 print('ars_tempidx', ars_temp)
        #                 ars_temp.num_factor *= s23[k].num_factor
        #                 print('ars_temp nm', ars_temp)
        #                 # usuwam operator s
        #                 for cc in range(0, len(ars_temp.coefficient)):
        #                     if ars_temp.coefficient[cc] == S_AMPLITUDE:
        #                         ars_temp.coefficient.pop(cc)
        #                         ars_temp.coefficient_idx.pop(cc)
        #                         break
        #                 print('ars temps', ars_temp) 
        #                 rsimp2.append(ars_temp)

        # sys.exit(0)


        # tworze kopie, zeby w pierwszej kopii za s wstawic s21  a w drugiej s23
        rsimp1 = deepcopy(rsimp)

        print('do tego bede wstawiac')
        print('rsimp1')
        for x in rsimp1:
            print(x)
        print('')
        print('')

        # chamska zamiama t2 na s2
        for x in range(0, len(rsimp1)):
            for j in range(0, len(rsimp1[x].coefficient)):
                coef = rsimp1[x].coefficient[j]
                if coef == S_AMPLITUDE:
                    rsimp1[x].coefficient[j] = CC_AMPLITUDE

        rsimp2 = substitute_evaluated_s_for_dummy_s(s23, rsimp, s_op)

        rsimp_all = rsimp1
        for x in rsimp2:
            rsimp_all.append(x)

        rsimp_all_s = simplify(rsimp_all)
        rsimp_all_s2 = simplify(rsimp_all_s)
        rsimp_all_s = simplify(rsimp_all_s2)

        # print('laciaty')    
        # rs_temp  =arithmetic_string(rsimp_all_s[71]) + arithmetic_string(rsimp_all_s[103])

        # print(rs_temp)
        # rs = simplify(rs_temp)
        # print(rs_temp)

        # sys.exit(0)


        print('to jest wynik', len(rsimp_all_s))

        k = 1
        for x in rsimp_all_s:
            print('&', x, '\\\\')
            k += 1
        # sys.exit(0)
        # block_oo, block_ov, block_vo, block_vv = fill_dm_blocks_dens_ground(rsimp_all_s)

        # print('block_oo')
        # for x in block_oo:
        #     print(x)
        # print('')
        # print('block_ov')
        # for x in block_ov:
        #     print(x)
        # print('')
        # print('block_vo')
        # for x in block_vo:
        #     print(x)
        # print('')
        # print('block_vv')
        # for x in block_vv:
        #     print(x)
        # print('')

        group(rsimp_all_s)
    #    sys.exit(0)


        # procedura ktora znajduje wymienne wklady ktore chcial Polnon
        print('WITAMMM')

        rsimp_tau = find_tau(rsimp_all_s)

        # find_exch_t(rsimp_all)



        latex_with_cost(f, 1, rsimp_tau, name=False)                    
        #f.write(rsimp_all_s)

        # sys.exit(0)
        # restr_abij = ['a', 'i', 'b', 'j']
        # Wm1_commutators_copy = []
        # for i in range(0, len(Wm1_commutators)):
        #     for j in range(0, len(Wm1_commutators[i])):
        #         print('biore komutator Wm1', Wm1_commutators[i][j])
        #         mbpt_order = compute_mbpt_for_density_Wm_ground_mbpt(Wm1_commutators[i][j])
        #         if mbpt_order == 2:
        #             # wstawiam S2 = T2 + P2([T2*,[T2,T2]])
        #             # trzeba zamienic indeksy
        #             Wm1_commutators_copy = deepcopy(Wm1_commutators[i][j])
        #             for k in range(0, len(Wm1_commutators[i][j].coefficient)):
        #                 print('sprawdzam współczynnik', Wm1_commutators[i][j].coefficient[k])
        #                 # Dla każdego współczynnika coef ktory jest s-em, musze podmienic kazdy s z listy s21 i s23
        #                 coef = Wm1_commutators[i][j].coefficient[k]
        #                 coef_idx = Wm1_commutators[i][j].coefficient_idx[k]
        #                 coef_sum = Wm1_commutators[i][j].summation
        #                 if coef == S_AMPLITUDE:
        #                     print('mam wspolczynnik s', k, coef, coef_idx, coef_sum)
        #                     restr2 = coef_idx
        #                     restr3 = coef_sum
        #                     s21_copy = deepcopy(s21)
        #                     for l in range(0, len(s21_copy)):
        #                         #1 zamien indeksy sumowania tak, żeby nie wolno bylo skorzystac z inkdesow restr2:'cbkl' i abij
        #                         print('1   ', s21_copy[l], restr_abij, restr2)
        #                         s21_copy[l].rename_summation_with_restriction(restr_abij, restr2)
        #                         print('2   ',s21_copy[l])
        #                         s21_copy[l].multisubst(restr_abij, restr2)
        #                         print('3   ',s21_copy[l])
        #                         s21_copy[l].rename_summation_with_restriction(restr2, restr3)
        #                         print(i, j, k, l, Wm1_commutators[i][j], s21_copy[l])


        #         elif mbpt_order == 4:
        #             # wstawiam S2 = T2
        #             a = 1

        #         #print(y,  compute_mbpt_for_density_Wm_ground_mbpt(y))

        # sys.exit(0)
        return W_middle_out, Wm1_commutators, rsimp_for_out


def substitute_evaluated_s_for_dummy_s(s23, rsimp, s_op):

    # petle w ktorych zmieniane jest s2 na s23 policzone wczesniej
    #1. zamien indeksy sumowania tak, żeby nie wolno bylo skorzystac z inkdesow cbkl i abij
    #   \sum_{edmn}t_{jn}^{bd}t_{mi}^{ae}t_{nm}^{ed}E_{jb}E_{ia}
    #2. zamien indeksy abij -> cbkl
    #\sum_{edmn}t_{ln}^{bd}t_{mi}^{ce}t_{nm}^{ed}E_{lb}E_{kc}
    #3. zamien indeksy symowania tak, zeby nie wolno było korzystać z indeksów paij


    if s_op == 2:
        restr_aibj = ['a', 'i', 'b', 'j']
        s_all = deepcopy(s23)
    elif s_op == 1:
        restr_aibj = ['a', 'i']
        s_all = deepcopy(s13)
    rsimp2 = arithmetic_string()
    for x in range(0, len(rsimp)):
        print('biore comut', rsimp[x])
        for j in range(0, len(rsimp[x].coefficient)):
            coef = rsimp[x].coefficient[j]
            coef_idx = rsimp[x].coefficient_idx[j]
            coef_sum = rsimp[x].summation
            restr2 = coef_idx
            restr3 = coef_sum
            if coef == S_AMPLITUDE:
                print('mam wspolczynnik s')
                for k in range(0, len(s_all)):
                    s23_copy = deepcopy(s_all)
                    #1 zamien indeksy sumowania tak, żeby nie wolno bylo skorzystac z inkdesow restr2:'cbkl' i abij
                    print('1   ', s23_copy[k], restr_aibj, restr2)
                    s23_copy[k].rename_summation_with_restriction(restr_aibj, restr2)
                    print('2   ',s23_copy[k])
                    s23_copy[k].multisubst(restr_aibj, restr2)
                    print('3   ',s23_copy[k], restr2, restr3)
                    s23_copy[k].rename_summation_with_restriction(restr2, restr3)
                    print(x, j, k, rsimp[x], s23_copy[k])
                    # tworze wyraz bedacy suma s i rsimp[x]
                    ars_temp = deepcopy(rsimp[x])
                    ars_temp.summation += s23_copy[k].summation
                    ars_temp.coefficient += s23_copy[k].coefficient
                    ars_temp.coefficient_idx += s23_copy[k].coefficient_idx
                    ars_temp.num_factor *= s23[k].num_factor
                    print('ars temp', ars_temp)
                    # usuwam operator s
                    for cc in range(0, len(ars_temp.coefficient)):
                        if ars_temp.coefficient[cc] == S_AMPLITUDE:
                            ars_temp.coefficient.pop(cc)
                            ars_temp.coefficient_idx.pop(cc)
                            break
                    print('ars temp po usunieciu s', ars_temp)

                    rsimp2.append(ars_temp)
            print('')
            print('')

    return rsimp2

def group(block):

    idx_list_big = []
    for i in range(0, len(block)):
        term_i = block[i]
        virt_list = []
        occ_list = []
        idx_lst = []
        count_dict = {}
        for j in range(0, len(term_i.coefficient)):
            if term_i.coefficient[j] == CC_AMPLITUDE:
                occ_idx_list = [term_i.coefficient_idx[j][1],term_i.coefficient_idx[j][3]]
                virt_idx_list = [term_i.coefficient_idx[j][0],term_i.coefficient_idx[j][2]]
                # virt_list = virt_list + virt_idx_list
                # occ_list = occ_list + occ_idx_list
                occ_idx_list.sort()
                virt_idx_list.sort()
                idx_lst.append(virt_idx_list)
                idx_lst.append(occ_idx_list)
        
                # idx_lst = idx_lst + occ_idx_list + virt_idx_list
                # for k in term_i.summation:
                #     count_dict[k] = idx_lst.count(k)
        print(term_i)
        print(idx_lst)
        print('')
        idx_list_big.append(idx_lst)

    banned = []
    ars_group = []
    idx_big = []
    for i in range(0, len(block)):
        if i not in banned:
            list_i = idx_list_big[i]
            ars_group.append([block[i]])
            banned.append(i)
            idx_big.append([i+1])
            for j in range(i+1, len(block)):
                if j not in banned:
                    list_j = idx_list_big[j]
                    if list_j == list_i:
                        ars_group[-1].append(block[j])
                        banned.append(j)
                        idx_big[-1].append(j+1)

    for i in range(0, len(ars_group)):
        for j in ars_group[i]:
            print('&', j, '\qquad', idx_big[i], '\qquad\\\\')
        print('\\\\')
                

#    sys.exit(0)
                

def find_tau(r):

    # Find instances of 2*t^ab_ij t^bc_ji in one term and t^ab_ij t^bc_ij in other term
    # replace this with tau = 2t^ab_ij(2 t^bc_ji - t^bc_ij)

    origin_list = []

    for i in range(0, len(r)):
        origin_list.append({i+1})
    
    print('------------------------------------------')
    for xx in range(0, 3):
        big_banned= []
        for i in range(0, len(r)):
            if i not in big_banned:
                print('riii', r[i])
                num_i = r[i].num_factor
                for j in range(i+1, len(r)):
                    if j not in big_banned and i not in big_banned:
                        num_j = r[j].num_factor
                        print('rjj', r[j])
                        # Find can
                        cn1 = compare_num(2.0*num_i, num_j, 'minus')
                        cn2 = compare_num(num_i, 2.0*num_j, 'minus')
                        cn3 = compare_num(num_i, 0.0, 'minus')
                        if (cn1 or cn2) and not cn3:                
                            term_i = r[i]
                            term_j = r[j]
                            banned_i = []
                            banned_j = []
                            done = False
                            print('sprawdzam wyrazy ktore maja rowne wspolczynniki',i, j)
                            print(r[i])
                            print(r[j])
                            print('')
                            for i1 in range(0, len(term_i.coefficient)):
                                if (i1 not in banned_i):
                                    if term_i.coefficient[i1] == CC_AMPLITUDE:
                                        for j1 in range(0, len(term_j.coefficient)):
                                            if (j1 not in banned_j):
                                                if term_j.coefficient[i1] == CC_AMPLITUDE:
                                                    if term_i.coefficient_idx[i1] == term_j.coefficient_idx[j1]:
                                                        banned_i.append(i1)
                                                        banned_j.append(j1)
                                                        done = True
                                                        break
                                if done:
                                    break

                            if done:
                                print('tak, maja rowne wspolczynniki i  jedna amp')
                                done2 = False
                                for i1 in range(0, len(term_i.coefficient)):
                                    if (i1 not in banned_i):
                                        if term_i.coefficient[i1] == CC_AMPLITUDE:
                                            occ_idx_list_i = [term_i.coefficient_idx[i1][1],term_i.coefficient_idx[i1][3]]
                                            virt_idx_list_i = [term_i.coefficient_idx[i1][0],term_i.coefficient_idx[i1][2]]
                                            occ_idx_list_i_sort = deepcopy(occ_idx_list_i)
                                            occ_idx_list_i_sort.sort()
                                            for j1 in range(0, len(term_j.coefficient)):
                                                if (j1 not in banned_j):
                                                    if term_j.coefficient[j1] == CC_AMPLITUDE:
                                                        occ_idx_list_j = [term_j.coefficient_idx[j1][1],term_j.coefficient_idx[j1][3]]
                                                        virt_idx_list_j = [term_j.coefficient_idx[j1][0],term_j.coefficient_idx[j1][2]]
                                                        occ_idx_list_j_sort = deepcopy(occ_idx_list_j)
                                                        occ_idx_list_j_sort.sort()
                                                        c1 = (virt_idx_list_i == virt_idx_list_j)
                                                        c2 = (occ_idx_list_i_sort==occ_idx_list_j_sort)
                                                        c3 = (occ_idx_list_i==occ_idx_list_j)
                                                        if c1 and c2 and not c3:
                                                            banned_i.append(i1)
                                                            banned_j.append(j1)
                                                            print('te dwa sa rowne12 i maja tau')
                                                            print(r[i])
                                                            print(r[j], i, j, banned_i, banned_j)
                                                            print('')
                                                            big_banned.append(i)
                                                            big_banned.append(j)
                                                            print(big_banned)
                                                            if cn1:
                                                                r[j].num_factor = 0.0
                                                                r[i].coefficient[banned_i[1]] = CC_TAU
                                                                origin_list[i] = origin_list[i]|{j+1}
                                                                origin_list[i] = origin_list[i]|origin_list[j]
                                                                origin_list[j] = {999}

                                                            if cn2:
                                                                r[i].num_factor = 0.0
                                                                r[j].coefficient[banned_j[1]] = CC_TAU
                                                                origin_list[j] = origin_list[j]|{i+1}
                                                                origin_list[j] = origin_list[j]|origin_list[j]
                                                                origin_list[i] = {999}

                                                                



                                                            done2 = True
                                                            print('przerywam1')
                                                            break
                                    if done2:
                                        print('przerywam2')
                                        break
                        
    r.cleanup()
    origin_list_cleanup = []
    for x in origin_list:
        if x != {999}:
            origin_list_cleanup.append(x)
    

    print('po znaleznienu tau', len(r), len(origin_list_cleanup))
    k = 0
    for x in r:
        print('&', x, '\mbox{powstalem z wyrazow}', origin_list_cleanup[k], '\\')
        k += 1
                            
    return r



def find_exch_t(r):
    print('zaczynam find_exch')

    t_exch_pair_big_list = [] 
    
    for i in range(0, len(r)):
        term = r[i]
        print(term)
        banned = []
        t_exch_pair_list = []
        for j in range(0, len(term.coefficient)):
            if (j not in banned):
                if term.coefficient[j] == CC_AMPLITUDE:
                    tidx_j = term.coefficient_idx[j]
                    occ_idx_list_j = [term.coefficient_idx[j][1],term.coefficient_idx[j][3]]
                    virt_idx_list_j = [term.coefficient_idx[j][0],term.coefficient_idx[j][2]]
                    occ_idx_list_j_sort = deepcopy(occ_idx_list_j)
                    occ_idx_list_j_sort.sort()

                    for k in range(j+1, len(term.coefficient)):
                        if k not in banned:
                            if term.coefficient[k] == CC_AMPLITUDE:
                                tidx_k = term.coefficient_idx[k]
                                occ_idx_list_k = [term.coefficient_idx[k][1],term.coefficient_idx[k][3]]
                                virt_idx_list_k = [term.coefficient_idx[k][0],term.coefficient_idx[k][2]]
                                occ_idx_list_k_sort = deepcopy(occ_idx_list_k)
                                occ_idx_list_k_sort.sort()
                                if occ_idx_list_j_sort == occ_idx_list_k_sort:
                                    if occ_idx_list_j != occ_idx_list_k:
                                        set_virt_k = set(virt_idx_list_k)
                                        common_elements = set_virt_k.intersection(virt_idx_list_j)
                                        if len(common_elements) == 1:
                                            print(j, k, tidx_j, tidx_k)
                                            print('occ_idx_list_j', occ_idx_list_j, occ_idx_list_j_sort)
                                            print('occ_idx_list_k', occ_idx_list_k, occ_idx_list_k_sort)
                                            print('virt_idx_list_j', virt_idx_list_j)
                                            print('virt_idx_list_k', virt_idx_list_k)
                                            print('common', common_elements)
                                            t_exch_pair_list.append([j+1, k+1])
                                            banned.append(j)
                                            banned.append(k)
                                            print('dodaje', j, k)
                                            
        t_exch_pair_big_list.append(t_exch_pair_list)
        print(t_exch_pair_list)                                    

        print('')

    r_sorted = []
    t_sorted = []
    for j in range(0, 5):
        print('ILOSC T EXCH', j)
        r_sorted.append([])
        t_sorted.append([])
        for i in range(0, len(r)):
            if len(t_exch_pair_big_list[i]) == j:
                r_sorted[j].append(r[i])
                t_sorted[j].append(t_exch_pair_big_list[i])
                print(r[i], t_exch_pair_big_list[i], i)
        print('')
                


def commutators_density_ground_f12(obs, W_middle, t_list, s_list):

    W_middle_commutators = []
    Wm1_commutators = []
    W_middle_out = []
    
    # licze komutatory wewnetrzne i zewnetrzne

    for x in range(0, len(W_middle)):
        for y in range(0, len(W_middle[x][t_list])):
            if W_middle[x][t_list][y] == 3:
                W_middle[x][t_list][y] = 'tf'

    for x in range(0, len(W_middle)):
        wm1t = W_middle[x][t_list]
        wm1s = W_middle[x][s_list]
        # print(wm1t)
        n_wm1t = len(set(permutations(wm1t, len(wm1t))))
        n_wm1s = len(set(permutations(wm1s, len(wm1s))))
        n = W_middle[x]['n']


        lt1 = len(wm1t)
        ls1 = len(wm1s)

        W_middle_out.append(W_middle[x])
        

        # print(W_middle[x])
        # print('')
        print('z wyrazu numer', x, W_middle[x])
        wm1_outer = density_f12_inner_outer2(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s, n)

        Wm1_commutators.append(deepcopy(wm1_outer))


    print('')
    print('dla tego operatora mam')

    density_int = arithmetic_string()
    for i in range(0, len(Wm1_commutators)):
        for j in range(0, len(Wm1_commutators[i])):
            print(Wm1_commutators[i][j])
            density_int.append(Wm1_commutators[i][j])
    print('')
    # calkowanie
    print('integrate')
    for x in density_int:
        print(x)

    vspace(0)
    rint = density_int.integrate()
    print('po integratea')
    
    for x in rint:
        print(x)
    print('')
    # simplifikacja
    rsimp = simplify(rint)
    print('rsimp')
    for x in rsimp:
        print(x)
    print('len(rsimp)', len(rsimp))
    print('')
    return W_middle_out, Wm1_commutators, rsimp


def commutators_cumulant(obsg1, obsg2, method, W_middle):

    W_middle_commutators = []
    Wm1_commutators = []
    W_middle_out = []
    
    # licze komutatory wewnetrzne i zewnetrzne
    

    G_commutators = []
    print('len(W_middle) is', len(W_middle))
    print('')
    print('and the W_middle itself is:')
    print('')    
    for x in W_middle:
        print(x)
        print('end of W_middle-------------------')
    vspace(0)
    print('The matrix element is')
    print('')
    print('Gm =  <e(S*)e(-T) E_pq e(T) e(-S*) P(e(S*)e(-T) E_rs e(T) e(-S*))>')
    print('           #--------Gm1----------#  #-------------Gm2-------------# ' )
    vspace(1)
    
    for x in range(0, len(W_middle)):
        # print(W_middle[x])
        G_commutators.append([])
        # First evaluate the Gm2 part
        print('EVALUATE GM2 part nr', x)
        print(W_middle[x])
        
        vspace(0)
        wm2t = W_middle[x]['Gm2_T_list']
        wm2s = W_middle[x]['Gm2_S_list']
        n = W_middle[x]['Gm2_n']
        
        # print(wm2t)
        n_wm2t = len(set(permutations(wm2t, len(wm2t))))
        n_wm2s = len(set(permutations(wm2s, len(wm2s))))

        lt2 = len(wm2t)
        ls2 = len(wm2s)

        if method == 'ccd' or method == 'cc3':
            wm2_outer = density_ccd_inner_outer(obsg2, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s, n)
        elif method == 'f12':
            wm2_outer = density_f12_inner_outer2(obsg2, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s, n)


        # Check if excitation level of Gm2 part is equal to the conditions
        # that were derived for this combination of S and T amplitudes

        # print('wm2_outer-przed', wm2_outer)

        
        wm2_outer = check_excitation_level(W_middle[x], wm2_outer)

        # Czynniki numeryczne wynikajace z komutatorow
        nf2t = compute_num_factor(W_middle[x]['Gm2_T_list'])
        nf2s = compute_num_factor(W_middle[x]['Gm2_S_list'])

        # Czynnik numeryczny z operatora P_3
        # nfP = W_middle[x]['Gm2_P_list']
        # print('nfP', nfP)
        nfP = 1.0
        num = 1/(nf2t*nf2s*nfP)

        if len(wm2_outer) > 0:
            print('Gm2 part', W_middle[x]['Gm2_T_list'], W_middle[x]['Gm2_S_list'],nf2t, nf2s, num)
            for xx in range(0, len(wm2_outer)):
                print('przed num', wm2_outer[xx])        
                wm2_outer[xx].num_factor *= num
                print(wm2_outer[xx])
            print('')
        else:
            print('Gm2 part is', 0)
        print('------------------------------------------------------------')
            
        # Next evaluate the Gm1 part
        vspace(1)
        print(x, 'EVALUATE GM1 part')
        print('operator is', obsg1)
        print(W_middle[x])
        wm1t = W_middle[x]['Gm1_T_list']
        wm1s = W_middle[x]['Gm1_S_list']

        n = W_middle[x]['Gm1_n']
        # print(wm1t)
        # print(wm1s)
        # print(wm1t)
        n_wm1t = len(set(permutations(wm1t, len(wm1t))))
        n_wm1s = len(set(permutations(wm1s, len(wm1s))))

        lt1 = len(wm1t)
        ls1 = len(wm1s)
        print('METHOD', method)
        if method == 'ccd' or method == 'cc3':
            wm1_outer = density_ccd_inner_outer(obsg1, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s, n)
        elif method == 'f12':
            wm1_outer = density_f12_inner_outer2(obsg1, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s, n)


        # Czynniki numeryczne wynikajace z komutatorow
        nf1t = compute_num_factor(W_middle[x]['Gm1_T_list'])
        nf1s = compute_num_factor(W_middle[x]['Gm1_S_list'])
        num1 = 1/(nf1t*nf1s)

        if len(wm1_outer) > 0:
            print('Gm1 part is', W_middle[x]['Gm1_S_list'], W_middle[x]['Gm1_T_list'])
            for xx in range(0, len(wm1_outer)):
                wm1_outer[xx].num_factor *= num1
                print(wm1_outer[xx])

            print('')
            print('&&&&')
            print('')
        else:
            print('Gm1 part is', 0)

        G_commutators[x].append(wm1_outer)
        G_commutators[x].append(wm2_outer)

    return G_commutators

def check_excitation_level(W_middle_one, wm2_outer):

    for x in range(0, len(wm2_outer)):
        exc_lev = compute_excitation_level(wm2_outer[x])
        if exc_lev != W_middle_one['Gm2_P_list']:
            print('excexcexc', exc_lev ,W_middle_one['Gm2_P_list'])
            wm2_outer[x].num_factor = 0.0

    return wm2_outer
        


def compute_excitation_level(elem):

    exc_lev = 0
    for i in range(0, len(elem.operator_idx)):
        if elem.operator_idx[i][0] in virtualall and elem.operator_idx[i][1] in occupied:
            exc_lev += 1 
        elif elem.operator_idx[i][0] in occupied and elem.operator_idx[i][1] in virtualall:
            exc_lev += -1

    return exc_lev
# def commutators_cumulant_f12(obs, obsy, G_middle):

#     G_middle_commutators = []
#     Wm1_commutators = []
    
#     # licze komutatory wewnetrzne i zewnetrzne

#     for x in range(0, len(G_middle)):
#         for y in range(0, len(G_middle[x]['Gm1_T_list'])):
#             if G_middle[x][t_list][y] == 3:
#                 G_middle[x][t_list][y] = 'tf'
#         for y in range(0, len(G_middle[x]['Gm2_T_list'])):
#             if G_middle[x][t_list][y] == 3:
#                 G_middle[x][t_list][y] = 'tf'

#     for x in range(0, len(G_middle)):
#         wm1t = G_middle[x]['Gm1_T_list']
#         wm1s = G_middle[x]['Gm2_T_list']

#         n_wm1t = len(set(permutations(wm1t, len(wm1t))))
#         n_wm1s = len(set(permutations(wm1s, len(wm1s))))

#         lt1 = len(wm1t)
#         ls1 = len(wm1s)

#         wm1_outer = density_f12_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)

#         Wm1_commutators.append(deepcopy(wm1_outer))

#     for x in range(0, len(G_middle)):
#         wm1t = G_middle[x]['Gm2_T_list']
#         wm1s = G_middle[x]['Gm2_T_list']

#         n_wm1t = len(set(permutations(wm1t, len(wm1t))))
#         n_wm1s = len(set(permutations(wm1s, len(wm1s))))

#         lt1 = len(wm1t)
#         ls1 = len(wm1s)

#         wm1_outer = density_f12_inner_outer(obsy, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)

#         Wm1_commutators.append(deepcopy(wm1_outer))

#     return Gm1_commutators, Gm2_commutators





# def commutators_quadra_Wm(W_middle):

#     #-----------------------------------------------------------------------------------------------------------------------------
#     # Commutators_quadra_Wm evaluates all the commutators in W_middle. This subroutine gives back
#     # three lists of arithmetic string. The mu_lists contains the excitation of mu for all the elements
#     # in the Wm_c list. 
#     #
#     # len(Wm2_c) = len(Wm2_mu_list)
#     # len(Wm3_c) = len(Wm3_mu_list)
#     #
#     # Commutators_quadra_Wlr does the same for W_left and W_right. There is only one mu_list as Wl_c and Wr_c share the same
#     # mu_list.
#     #
#     #-----------------------------------------------------------------------------------------------------------------------------

#     W_middle_commutators = []
#     Wm1_commutators = []
#     Wm2_commutators = []
#     Wm3_commutators = []
#     Wm2_mu_list = []
#     Wm3_mu_list = []

#     for x in range(0, len(W_middle)):
# #        print('evaluating', x, 'of', len(W_middle))
#         start = time.time()
# #        print(x, W_middle[x])

#         wm1t = W_middle[x]['Wm1_T_list']
#         wm1s = W_middle[x]['Wm1_S_list']

#         wm2t = W_middle[x]['Wm2_T_list']
#         wm2s = W_middle[x]['Wm2_S_list']

#         wm3s = W_middle[x]['Wm3_S_list']

        
#         n_wm1t = len(set(permutations(wm1t, len(wm1t))))
#         n_wm1s = len(set(permutations(wm1s, len(wm1s))))

#         n_wm2t = len(set(permutations(wm2t, len(wm2t))))
#         n_wm2s = len(set(permutations(wm2s, len(wm2s))))
        
#         n_wm3s = len(set(permutations(wm3s, len(wm3s))))


    
#         lt1 = len(wm1t)
#         ls1 = len(wm1s)

#         lt2 = len(wm2t)
#         ls2 = len(wm2s)

#         ls3 = len(wm3s)

#         wm2_mu = ugg()
#         if (W_middle[x]['Wm2_n']) == 1:
#             wm2_mu = mu1
#             Wm2_mu_list.append(1)
#         elif (W_middle[x]['Wm2_n']) == 2:
#             wm2_mu = mu2
#             Wm2_mu_list.append(2)
#         elif (W_middle[x]['Wm2_n']) == 3:
#             if (W_middle[x]['Wm3_n']) == 3:
#                 print('excluded 33')
#             else:
#                 wm2_mu = mu3
#                 Wm2_mu_list.append(3)

#         wm3_mu = ugg()
#         eomr = ugg()
#         if (W_middle[x]['Wm3_n']) == 1:
#             # wm3_mu = mu1p
#             Wm3_mu_list.append(1)
#             eomr = eomrr1
#         elif (W_middle[x]['Wm3_n']) == 2:
#             # wm3_mu = mu2p
#             Wm3_mu_list.append(2)
#             eomr = eomrr2
#         elif (W_middle[x]['Wm3_n']) == 3:
#             if (W_middle[x]['Wm2_n']) == 3:
#                 print('excluded 33')
#             else:
#                 eomr = eomrr3
#                 Wm3_mu_list.append(3)

#         if lt1 == 0:
#             wm1_inner = arithmetic_string(obs)
#         else:
#             wm1_inner = arithmetic_string(obs)
#             for i in range(0, lt1):
#                 wm1_inner = evaluate(wm1_inner, stoo(wm1t[i], "t", True))
                
#         wm1_outer = wm1_inner
#         if ls1 !=0:
#             for i in range(0, ls1):
#                 wm1_outer = evaluate(wm1_outer, stoo(wm1s[i], "s", False))

#             wm1_outer = wm1_outer.scale(-1.0)
#             wm1_outer = wm1_outer.scale(n_wm1t)
#             wm1_outer = wm1_outer.scale(n_wm1s)
#         if len(wm1_outer) == 0:
#             print('wm1zero')
#             sys.exit(0)
#         Wm1_commutators.append(wm1_outer)
#         end = time.time()
#         start = time.time()
#         # Wm2
#         if lt2 == 0:
#             wm2_inner = arithmetic_string(wm2_mu)
#         else:
#             wm2_inner = arithmetic_string(wm2_mu)
#             for i in range(0, lt2):
#                 wm2_inner = evaluate(wm2_inner, stoo(wm2t[i], "t", True))

#         wm2_outer = wm2_inner
#         if ls2 !=0:
#             for i in range(0, ls2):
#                 wm2_outer = evaluate(wm2_outer, stoo(wm2s[i], "s", False))


#             wm2_outer = wm2_outer.scale(-1.0)
#             wm2_outer = wm2_outer.scale(n_wm2t)
#             wm2_outer = wm2_outer.scale(n_wm2s)
#         if len(wm2_outer) == 0:
#             print('wm2zero')
#             sys.exit(0)
#         Wm2_commutators.append(wm2_outer)
#         end = time.time()
#         start = time.time()

#         # Wm3
#         wm3_outer = arithmetic_string(eomr)
#         if ls3 !=0:
#             for i in range(0, ls3):
#                 wm3_outer = evaluate(wm3_outer, stoo(wm3s[i], "s", True))

#             wm3_outer = wm3_outer.scale(-1.0)
#             wm3_outer = wm3_outer.scale(n_wm3s)
#         if len(wm3_outer) == 0:
#             print('wm3zero')
#             sys.exit(0)

#         Wm3_commutators.append(wm3_outer)
#         end = time.time()
    
#     print(len( Wm1_commutators), len(Wm2_commutators), len(Wm3_commutators), len(Wm2_mu_list),  len(Wm3_mu_list))


#     return Wm1_commutators,  Wm2_commutators,  Wm3_commutators, Wm2_mu_list, Wm3_mu_list


def commutators_quadra_Wm_triplet(W_middle):

    W_middle_commutators = []
    Wm1_commutators = []
    Wm2_commutators = []
    Wm3_commutators = []
    Wm2_mu_list = []
    Wm3_mu_list = []
    W_middle_out = []

    for x in range(0, len(W_middle)):
        wm1t = W_middle[x]['Wm1_T_list']
        wm1s = W_middle[x]['Wm1_S_list']

        wm2t = W_middle[x]['Wm2_T_list']
        wm2s = W_middle[x]['Wm2_S_list']

        wm3s = W_middle[x]['Wm3_S_list']

        n_wm1t = len(set(permutations(wm1t, len(wm1t))))
        n_wm1s = len(set(permutations(wm1s, len(wm1s))))

        n_wm2t = len(set(permutations(wm2t, len(wm2t))))
        n_wm2s = len(set(permutations(wm2s, len(wm2s))))

        n_wm3s = len(set(permutations(wm3s, len(wm3s))))

    
        lt1 = len(wm1t)
        ls1 = len(wm1s)

        lt2 = len(wm2t)
        ls2 = len(wm2s)

        ls3 = len(wm3s)



        if (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 1:
            W_middle_out.append(W_middle[x])
            wm2_mu = mu1_triplet
            eomr = eomrr1_triplet
            Wm2_mu_list.append(1)
            Wm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s) 
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 3:

            W_middle_out.append(W_middle[x])
            wm2_mu = mu3_triplet
            eomr = eomrr3_triplet
            Wm2_mu_list.append(3)
            Wm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 3:
            W_middle_out.append(W_middle[x])
            wm2_mu = mu1_triplet
            eomr = eomrr3_triplet
            Wm2_mu_list.append(1)
            Wm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 1:
            W_middle_out.append(W_middle[x])
            wm2_mu = mu3_triplet
            eomr = eomrr1_triplet
            Wm2_mu_list.append(3)
            Wm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 2:
            W_middle_out.append(W_middle[x])
            wm2_mu = mu1_triplet
            eomr_plus = eomrr2_triplet_plus
            eomr_minus = eomrr2_triplet_minus
            Wm2_mu_list.append(1)
            Wm2_mu_list.append(1)
            Wm3_mu_list.append('2p')
            Wm3_mu_list.append('2m')

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer)) # The double adding of wm1_outer etc is intentional
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr_plus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            wm3_outer = wm3_inner_outer(eomr_minus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
#            print(wm3_outer)

        elif (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 2:
            W_middle_out.append(W_middle[x])
            wm2_mu = mu3_triplet
            eomr_plus = eomrr2_triplet_plus
            eomr_minus = eomrr2_triplet_minus
            Wm2_mu_list.append(3)
            Wm2_mu_list.append(3)
            Wm3_mu_list.append('2p')
            Wm3_mu_list.append('2m')

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer)) # The double adding of wm1_outer etc is intentional                                                                             
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr_plus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            wm3_outer = wm3_inner_outer(eomr_minus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 1:
            W_middle_out.append(W_middle[x])
            wm2_mu_plus = mu2_triplet
            wm2_mu_minus = mu2_triplet
            eomr = eomrr1_triplet
            Wm2_mu_list.append('2p')
            Wm2_mu_list.append('2m')
            Wm3_mu_list.append(1)
            Wm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu_plus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            wm2_outer = wm2_inner_outer(wm2_mu_minus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 3:
            W_middle_out.append(W_middle[x])
            wm2_mu_plus = mu2_triplet
            wm2_mu_minus = mu2_triplet
            eomr = eomrr3_triplet
            Wm2_mu_list.append('2p')
            Wm2_mu_list.append('2m')
            Wm3_mu_list.append(3)
            Wm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu_plus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            wm2_outer = wm2_inner_outer(wm2_mu_minus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            Wm3_commutators.append(deepcopy(wm3_outer))


        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 2:
            W_middle_out.append(W_middle[x])
            wm2_mu_plus = mu2_triplet
            wm2_mu_minus = mu2_triplet
            eomr_plus = eomrr2_triplet_plus
            eomr_minus = eomrr2_triplet_minus
            Wm2_mu_list.append('2p')
            Wm2_mu_list.append('2p')
            Wm2_mu_list.append('2m')
            Wm2_mu_list.append('2m')

            Wm3_mu_list.append('2p')
            Wm3_mu_list.append('2m')
            Wm3_mu_list.append('2p')
            Wm3_mu_list.append('2m')

            wm1_outer = wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))
            Wm1_commutators.append(deepcopy(wm1_outer))
            Wm1_commutators.append(deepcopy(wm1_outer))
            Wm1_commutators.append(deepcopy(wm1_outer))

            
            wm2_outer = wm2_inner_outer(wm2_mu_plus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            wm2_outer = wm2_inner_outer(wm2_mu_plus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            wm2_outer = wm2_inner_outer(wm2_mu_minus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            wm2_outer = wm2_inner_outer(wm2_mu_minus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr_plus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            wm3_outer = wm3_inner_outer(eomr_minus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            wm3_outer = wm3_inner_outer(eomr_plus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            wm3_outer = wm3_inner_outer(eomr_minus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))


    print(len( Wm1_commutators), len(Wm2_commutators), len(Wm3_commutators), len(Wm2_mu_list),  len(Wm3_mu_list))

    return W_middle_out, Wm1_commutators,  Wm2_commutators,  Wm3_commutators, Wm2_mu_list, Wm3_mu_list




def commutators_quadra_Wm_singlet_triplet(W_middle):

    W_middle_commutators = []
    Wm1_commutators = []
    Wm2_commutators = []
    Wm3_commutators = []
    Wm2_mu_list = []
    Wm3_mu_list = []
    W_middle_out = []
    
    # for x in range(0, len(W_middle)):
    #     if (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 3:
    #         print('')
    #     else:
    #         W_middle_out.append(W_middle[x])
    W_middle_out = W_middle

    for x in range(0, len(W_middle)):
        wm1t = W_middle[x]['Wm1_T_list']
        wm1s = W_middle[x]['Wm1_S_list']

        wm2t = W_middle[x]['Wm2_T_list']
        wm2s = W_middle[x]['Wm2_S_list']

        wm3s = W_middle[x]['Wm3_S_list']

        n_wm1t = len(set(permutations(wm1t, len(wm1t))))
        n_wm1s = len(set(permutations(wm1s, len(wm1s))))

        n_wm2t = len(set(permutations(wm2t, len(wm2t))))
        n_wm2s = len(set(permutations(wm2s, len(wm2s))))

        n_wm3s = len(set(permutations(wm3s, len(wm3s))))

    
        lt1 = len(wm1t)
        ls1 = len(wm1s)

        lt2 = len(wm2t)
        ls2 = len(wm2s)

        ls3 = len(wm3s)


        # Wm =  <e(-S)e(T*) X e(-T*) e(S) P(e(-S)e(T*) \mu_n e(-T*) e(S))| P(e(S*)\mu_l e(-S*))>                     
        #--------Wm1----------#  #-------------Wm2-------------#  #-------Wm3-------#    
          
        #  Wm2 = \mu_n (singletowe)
        #  Wm3 = \mu_l (trypletowe)

        if (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 1:
#            wm2_mu = mu1
            wm2_mu = eomrl1
            eomr = eomrr1_triplet
            Wm2_mu_list.append(1)
            Wm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s) 
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 3:
            #wm2_mu = mu3
            wm2_mu = eomrl3
            eomr = eomrr3_triplet
            Wm2_mu_list.append(3)
            Wm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 3:
            print('----13 jeden')
            #wm2_mu = mu1
            em2_mu = eomrl1
            eomr = eomrr3_triplet
            Wm2_mu_list.append(1)
            Wm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 1:
            print('----31 jeden')
            #wm2_mu = mu3
            wm2_mu = eomrl3
            eomr = eomrr1_triplet
            Wm2_mu_list.append(3)
            Wm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 2:
            print('----12')
            #wm2_mu = mu1
            wm2_mu = eomrl1
            eomr_plus = eomrr2_triplet_plus
            eomr_minus = eomrr2_triplet_minus
            Wm2_mu_list.append(1)
            Wm2_mu_list.append(1)
            Wm3_mu_list.append('2p')
            Wm3_mu_list.append('2m')

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer)) # The double adding of wm1_outer etc is intentional
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr_plus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            wm3_outer = wm3_inner_outer(eomr_minus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 2:
            print('----32 dwa')

            #wm2_mu = mu3
            wm2_mu = eomrl3
            eomr_plus = eomrr2_triplet_plus
            eomr_minus = eomrr2_triplet_minus
            Wm2_mu_list.append(3)
            Wm2_mu_list.append(3)
            Wm3_mu_list.append('2p')
            Wm3_mu_list.append('2m')

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer)) # The double adding of wm1_outer etc is intentional                             
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr_plus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            wm3_outer = wm3_inner_outer(eomr_minus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 1:
            print('----21')
            #wm2_mu = mu2
            wm2_mu = eomrl2
            eomr = eomrr1_triplet
            Wm2_mu_list.append(2)
            Wm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            
            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 3:
            print('----23 jeden')
            #wm2_mu = mu2
            wm2_mu = eomrl2
            eomr = eomrr3_triplet
            Wm2_mu_list.append(2)
            Wm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))


        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 2:
            print('----22')
            #wm2_mu = mu2
            wm2_mu = eomrl2
            eomr_plus = eomrr2_triplet_plus
            eomr_minus = eomrr2_triplet_minus
            Wm2_mu_list.append(2)
            Wm2_mu_list.append(2)

            Wm3_mu_list.append('2p')
            Wm3_mu_list.append('2m')

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))
            Wm1_commutators.append(deepcopy(wm1_outer))

            
            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr_plus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            wm3_outer = wm3_inner_outer(eomr_minus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))


    print(len( Wm1_commutators), len(Wm2_commutators), len(Wm3_commutators), len(Wm2_mu_list),  len(Wm3_mu_list))
    # print('lenn')
    # print('wm1')

    # for x in Wm1_commutators:
    #     print(x)
    # print('wm2')
    # for x in Wm2_commutators:
    #     print(x)
    # print('wm3')
    # for x in Wm3_commutators:
    #     print(x)





    return W_middle_out, Wm1_commutators,  Wm2_commutators,  Wm3_commutators, Wm2_mu_list, Wm3_mu_list


def commutators_quadra_Wm_singlet_triplet_full(W_middle):

    W_middle_commutators = []
    Wm1_commutators = []
    Wm2_commutators = []
    Wm3_commutators = []
    Wm2_mu_list = []
    Wm3_mu_list = []

    for x in range(0, len(W_middle)):
        wm1t = W_middle[x]['Wm1_T_list']
        wm1s = W_middle[x]['Wm1_S_list']

        wm2t = W_middle[x]['Wm2_T_list']
        wm2s = W_middle[x]['Wm2_S_list']

        wm3s = W_middle[x]['Wm3_S_list']

        n_wm1t = len(set(permutations(wm1t, len(wm1t))))
        n_wm1s = len(set(permutations(wm1s, len(wm1s))))

        n_wm2t = len(set(permutations(wm2t, len(wm2t))))
        n_wm2s = len(set(permutations(wm2s, len(wm2s))))

        n_wm3s = len(set(permutations(wm3s, len(wm3s))))

    
        lt1 = len(wm1t)
        ls1 = len(wm1s)

        lt2 = len(wm2t)
        ls2 = len(wm2s)

        ls3 = len(wm3s)


        # Wm =  <e(-S)e(T*) X e(-T*) e(S) P(e(-S)e(T*) \mu_n e(-T*) e(S))| P(e(S*)\mu_l e(-S*))>                     
        #--------Wm1----------#  #-------------Wm2-------------#  #-------Wm3-------#    
          
        #  Wm2 = \mu_n (singletowe)
        #  Wm3 = \mu_l (trypletowe)

        if (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 1:
            wm2_mu = eomrl1
            eomr = eomrr1_triplet
            Wm2_mu_list.append(1)
            Wm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s) 
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 3:
            wm2_mu = eomrl3
            eomr = eomrr3_triplet
            Wm2_mu_list.append(3)
            Wm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 3:
            wm2_mu = eomrl1
            eomr = eomrr3_triplet
            Wm2_mu_list.append(1)
            Wm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 1:
            wm2_mu = eomrr3
            eomr = eomrr1_triplet
            Wm2_mu_list.append(3)
            Wm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 2:
            wm2_mu = eomrl1
            eomr_plus = eomrr2_triplet_plus
            eomr_minus = eomrr2_triplet_minus
            Wm2_mu_list.append(1)
            Wm2_mu_list.append(1)
            Wm3_mu_list.append('2p')
            Wm3_mu_list.append('2m')

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer)) # The double adding of wm1_outer etc is intentional
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr_plus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            wm3_outer = wm3_inner_outer(eomr_minus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 2:

            wm2_mu = eomrl3
            eomr_plus = eomrr2_triplet_plus
            eomr_minus = eomrr2_triplet_minus
            Wm2_mu_list.append(3)
            Wm2_mu_list.append(3)
            Wm3_mu_list.append('2p')
            Wm3_mu_list.append('2m')

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer)) # The double adding of wm1_outer etc is intentional                             
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr_plus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            wm3_outer = wm3_inner_outer(eomr_minus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 1:
            wm2_mu = eomrl2
            eomr = eomrr1_triplet
            Wm2_mu_list.append(2)
            Wm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            
            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 3:
            wm2_mu = eomrl2
            eomr = eomrr3_triplet
            Wm2_mu_list.append(2)
            Wm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))


        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 2:
            wm2_mu = eomrl2
            eomr_plus = eomrr2_triplet_plus
            eomr_minus = eomrr2_triplet_minus
            Wm2_mu_list.append(2)
            Wm2_mu_list.append(2)

            Wm3_mu_list.append('2p')
            Wm3_mu_list.append('2m')

            wm1_outer = wm1_inner_outer(aobst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))
            Wm1_commutators.append(deepcopy(wm1_outer))

            
            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr_plus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            wm3_outer = wm3_inner_outer(eomr_minus, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))


    print(len( Wm1_commutators), len(Wm2_commutators), len(Wm3_commutators), len(Wm2_mu_list),  len(Wm3_mu_list))
    print('lenn')
    print('wm1')
    for x in Wm1_commutators:
        print(x)
    print('wm2')
    for x in Wm2_commutators:
        print(x)
    print('wm3')
    for x in Wm3_commutators:
        print(x)



    return Wm1_commutators,  Wm2_commutators,  Wm3_commutators, Wm2_mu_list, Wm3_mu_list



def commutators_quadra_Wm_triplet_singlet(W_middle):

    W_middle_commutators = []
    Wm1_commutators = []
    Wm2_commutators = []
    Wm3_commutators = []
    Wm2_mu_list = []
    Wm3_mu_list = []

    for x in range(0, len(W_middle)):
        wm1t = W_middle[x]['Wm1_T_list']
        wm1s = W_middle[x]['Wm1_S_list']

        wm2t = W_middle[x]['Wm2_T_list']
        wm2s = W_middle[x]['Wm2_S_list']

        wm3s = W_middle[x]['Wm3_S_list']

        n_wm1t = len(set(permutations(wm1t, len(wm1t))))
        n_wm1s = len(set(permutations(wm1s, len(wm1s))))

        n_wm2t = len(set(permutations(wm2t, len(wm2t))))
        n_wm2s = len(set(permutations(wm2s, len(wm2s))))

        n_wm3s = len(set(permutations(wm3s, len(wm3s))))

    
        lt1 = len(wm1t)
        ls1 = len(wm1s)

        lt2 = len(wm2t)
        ls2 = len(wm2s)

        ls3 = len(wm3s)


        # Wm =  <e(-S)e(T*) X e(-T*) e(S) P(e(-S)e(T*) \mu_n e(-T*) e(S))| P(e(S*)\mu_l e(-S*))>                     
        #--------Wm1----------#  #-------------Wm2-------------#  #-------Wm3-------#    
          
        #  Wm2 = \mu_n (singletowe)
        #  Wm3 = \mu_l (trypletowe)

        if (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 1:
            wm2_mu = mu1_triplet
            eomr = eomrr1
            Wm2_mu_list.append(1)
            Wm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(obst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s) 
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 3:
            wm2_mu = mu3_triplet
            eomr = eomrr3
            Wm2_mu_list.append(3)
            Wm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(obst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 3:
            wm2_mu = mu1_triplet
            eomr = eomrr3
            Wm2_mu_list.append(1)
            Wm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(obst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 1:
            wm2_mu = mu3_triplet
            eomr = eomrr1
            Wm2_mu_list.append(3)
            Wm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(obst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 1 and (W_middle[x]['Wm3_n']) == 2:
            wm2_mu = mu1_triplet
            eomr = eomrr2
            Wm2_mu_list.append(1)
            Wm3_mu_list.append(2)

            wm1_outer = wm1_inner_outer(obst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))

        elif (W_middle[x]['Wm2_n']) == 3 and (W_middle[x]['Wm3_n']) == 2:

            wm2_mu = mu3_triplet
            eomr = eomrr2
            Wm2_mu_list.append(3)
            Wm3_mu_list.append(2)

            wm1_outer = wm1_inner_outer(obst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
 
        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 1:
            wm2_mu_plus = mu2_triplet
            wm2_mu_minus = mu2_triplet
            eomr = eomrr1
            Wm2_mu_list.append('2p')
            Wm2_mu_list.append('2m')
            Wm3_mu_list.append(1)
            Wm3_mu_list.append(1)

            wm1_outer = wm1_inner_outer(obst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu_plus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            wm2_outer = wm2_inner_outer(wm2_mu_minus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            Wm3_commutators.append(deepcopy(wm3_outer))


        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 3:

            wm2_mu_plus = mu2_triplet
            wm2_mu_minus = mu2_triplet
            eomr = eomrr3
            Wm2_mu_list.append('2p')
            Wm2_mu_list.append('2m')
            Wm3_mu_list.append(3)
            Wm3_mu_list.append(3)

            wm1_outer = wm1_inner_outer(obst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))
            Wm1_commutators.append(deepcopy(wm1_outer))

            wm2_outer = wm2_inner_outer(wm2_mu_plus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            wm2_outer = wm2_inner_outer(wm2_mu_minus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            Wm3_commutators.append(deepcopy(wm3_outer))


        elif (W_middle[x]['Wm2_n']) == 2 and (W_middle[x]['Wm3_n']) == 2:
            wm2_mu_plus = mu2_triplet
            wm2_mu_minus = mu2_triplet
            eomr = eomrr2

            Wm2_mu_list.append('2p')
            Wm2_mu_list.append('2m')
            Wm3_mu_list.append(2)
            Wm3_mu_list.append(2)

            wm1_outer = wm1_inner_outer(obst, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
            Wm1_commutators.append(deepcopy(wm1_outer))
            Wm1_commutators.append(deepcopy(wm1_outer))

            
            wm2_outer = wm2_inner_outer(wm2_mu_plus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))
            wm2_outer = wm2_inner_outer(wm2_mu_minus, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s)
            Wm2_commutators.append(deepcopy(wm2_outer))

            wm3_outer = wm3_inner_outer(eomr, wm3s, ls3, n_wm3s)
            Wm3_commutators.append(deepcopy(wm3_outer))
            Wm3_commutators.append(deepcopy(wm3_outer))


    print(len( Wm1_commutators), len(Wm2_commutators), len(Wm3_commutators), len(Wm2_mu_list),  len(Wm3_mu_list))

    return Wm1_commutators,  Wm2_commutators,  Wm3_commutators, Wm2_mu_list, Wm3_mu_list


def wm3_inner_outer(eomr, wm3s, ls3, n_wm3s):

    wm3_outer = arithmetic_string(eomr)

    # for i in range(0, ls3):
    #     sss = stoo(wm3s[i], "s", True)
    #     print('przed', sss)
    #     for x in wm3_outer:
    #         # disambiguate(x, sss)
    #         disambiguate(sss, x)
    #         print(x, sss)

    #         zzz = evaluate(x, sss)
    #         print('')
    #         for r in zzz:
    #             print('zzz', r)
    #         print('')
    #     print('iii', i, wm3_outer, sss)

    # sys.exit(0)

    if ls3 !=0:
        for i in range(0, ls3):
            wm3_outer = evaluate(wm3_outer, stoo(wm3s[i], "s", True))
            # sss = stoo(wm3s[i], "s", True)
            # wm3_outer = evaluate(wm3_outer, sss)
        wm3_outer = wm3_outer.scale(-1.0)
        wm3_outer = wm3_outer.scale(n_wm3s)


    return wm3_outer

def wm1_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s):

    if lt1 == 0:
        wm1_inner = arithmetic_string(obs)
    else:
        wm1_inner = arithmetic_string(obs)
        for i in range(0, lt1):
            wm1_inner = evaluate(wm1_inner, stoo(wm1t[i], "t", True))

    wm1_outer = wm1_inner
    if ls1 !=0:
        for i in range(0, ls1):
            wm1_outer = evaluate(stoo(wm1s[i], "s", False), wm1_outer)
            # sss = stoo(wm1s[i], "s", False)
            # # for x in wm1_outer:
            # #     disambiguate(sss, x)
            # wm1_outer = evaluate(wm1_outer, sss)

        # wm1_outer = wm1_outer.scale(-1.0)
        wm1_outer = wm1_outer.scale(n_wm1t)
        wm1_outer = wm1_outer.scale(n_wm1s)
        
    print('wm1')
    for x in wm1_outer:
        print(x)

    return wm1_outer


def gm_inner_outer(obs, gm1t, gm1s, lt1, ls1, n_gm1t, n_gm1s):

    if lt1 == 0:
        gm1_inner = arithmetic_string(obs)
    else:
        gm1_inner = arithmetic_string(obs)
        for i in range(0, lt1):
            gm1_inner = evaluate(gm1_inner, stoo(gm1t[i], "t", False))

    gm1_outer = gm1_inner
    if ls1 !=0:
        for i in range(0, ls1):
            gm1_outer = evaluate(gm1_outer, stoo(gm1s[i], "s", True))

        gm1_outer = gm1_outer.scale(-1.0)
        gm1_outer = gm1_outer.scale(n_gm1t)
        gm1_outer = gm1_outer.scale(n_gm1s)
        
    print('gm1')
    for x in gm1_outer:
        print(x)

    return gm1_outer

def density_f12_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s, n):

    vspace(0)
    print('the operator is: ', obs)
    vspace(0)
    
    obs_exc = 0
    if obs.operator_idx[0][0]  in virtualall and obs.operator_idx[0][1] in occupied:
        obs_exc = 1
    elif obs.operator_idx[0][0]  in virtualall and obs.operator_idx[0][1] in virtualall:
        obs_exc = 0
    elif obs.operator_idx[0][0]  in occupied and obs.operator_idx[0][1] in occupied:
        obs_exc  = 0
    elif obs.operator_idx[0][0]  in occupied and obs.operator_idx[0][1] in virtualall:
        obs_exc = -1

    print('the excitation of the operator is', obs_exc, 'and should be', n)
    
    if obs_exc != n:
        wm1_outer = arithmetic_string()
        # print('nie bedzie', obs_exc, n)
        print('exit')
        return wm1_outer


    if lt1 == 0:
        print('wewnetrznego stringa e(-T) X e(T) nie ma, bo Gm1_T_list = []')
        wm1_inner = arithmetic_string(obs)
        # if wm1_inner is None:
        #     print('wm1_none')
        # else:
        #     print('wm1_inner', wm1_inner)
    else:
        vspace(0)
        print('wewnetrzny string e(-T) X e(T) jest, lt1=', lt1)
        vspace(0)
        wm1_inner = arithmetic_string(obs)
        print(wm1_inner)
        # if wm1_inner is None:
        #     print('wm1_none')
        # else:            
        #     print('wm1_inner1', wm1_inner, wm1t, lt1)
        
        for i in range(0, lt1):
            if wm1t[i] == 'tf':
                vspace(0)
                print('teraz amplituda T\'')
                vspace(0)
                wm1_inner = evaluate(wm1_inner, stoo(2, "tf", False))
            else:
                vspace(0)
                print('teraz amplituda zwykle T')
                vspace(0)
                print('rr', wm1_inner)
                print(stoo(wm1t[i], "t", False))
                wm1_inner = evaluate(wm1_inner, stoo(wm1t[i], "t", False))
#                print('rr', wm1_inner)

        vspace(0)
        # print('typ', type(wm1_inner))
        # if wm1_inner is None:
        #     print('wm1_none')
        # else:
        #     print('wm1_inner', wm1_inner)
        #     print('la')

    wm1_outer = wm1_inner
    
    print('wm1_inner', wm1_inner)
    if ls1 !=0:
        print('Gm2_S_list != []')
        for i in range(0, ls1):
            if wm1s[i] == 1:
                # print('wr0', wm1_outer)
                w1 = deepcopy(wm1_outer)
                g1 = stoo(wm1s[i], "s", True)
                print('wuuu0', w1, g1)
                for x in wm1_outer:
                    disambiguate(g1, x)
                wm1_outer1 = evaluate(g1, wm1_outer)
                g2 = stoo(wm1s[i], "sf12", True)
                print('wuuu', w1, g2)
                for x in w1:
                    disambiguate(g2, x)

                # print('wuuu', w1, g2)
                wm1_outer2 = evaluate(g2, w1)
                
                wm1_outer = wm1_outer1 + wm1_outer2
                # print('wr', wm1_outer1)
                # print('wr', wm1_outer2)
            if wm1s[i] == 2:
                print('jestem tu 1')
                # for f in wm1_outer:
                #     print(f)
                
                w1 = deepcopy(wm1_outer)
                print('w1', w1)
                w2 = deepcopy(wm1_outer)
                w3 = deepcopy(wm1_outer)

                g1 = stoo(wm1s[i], "s", True)
                print('samo g1', g1)
                for x in wm1_outer:
                    print('x przed dis', x)
                    disambiguate(g1, x)
                    print('x', g1, x)
                
                g2 = stoo(wm1s[i], "sf12", True, idx=1)
                print('samo g2', g2)
                # print('g2', g2)
                for x in w1:
                    disambiguate(g2, x)
                    print('g2', g2, x)
                g3 = stoo(wm1s[i], "sf12", True, idx=2)
                # print('g3', g3)
                for x in w2:
                    disambiguate(g3, x)
                    # print(g3, x)
                g4 = stoo(wm1s[i], "sf12", True, idx=3)
                # print('g4', g4)
                for x in w3:
                    disambiguate(g4, x)
                    # print(g4, x)

                print('ggg1', g1)
                vspace(0)
                print('----------------------------')
                vspace(0)
                print('wm1_outer')
                for x in wm1_outer:
                    print(x)
                vspace(0)
                print('ggg1', g1)
                vspace(0)

                wm1_outer1 = evaluate(g1, wm1_outer)
                if len(wm1_outer) > 0:
                    for x in wm1_outer1:
                        print('out1', x)

                vspace(0)
                # print('eval2')
                wm1_outer2 = evaluate(g2, w1)
                vspace(0)
                if len(wm1_outer2) > 0:
                    for x in wm1_outer2:                                        
                        print('out2', x)
                wm1_outer3 = evaluate(g3, w2)
                # if len(wm1_outer3) > 0:
                #     print('out3', wm1_outer3)
                wm1_outer4 = evaluate(g4, w3)
                # if len(wm1_outer4) > 0:
                #     print('out4', wm1_outer4)
#                wm1_outer = wm1_outer1 + wm1_outer2 + wm1_outer3 + wm1_outer4

                wm1_outer = wm1_outer1 + wm1_outer2
                
        if wm1_outer is not None:
            print('wm1out', wm1_outer)
            # if i%2 != 0:
            #     wm1_outer = wm1_outer.scale(-1.0)


        wm1_outer = wm1_outer.scale(n_wm1t)
        wm1_outer = wm1_outer.scale(n_wm1s)

    vspace(0)
    print('caly wynik to ')
    for x in wm1_outer:
        print(x)
    vspace(1)

    # print('wm1-outer-haaa')
    # for x in wm1_outer:
    #     print(x)
    # print('')
    # print('nie')
    return wm1_outer


def density_f12_inner_outer2(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s, n, CCAABBSS = False):

    vspace(0)
    print('the operator is: ', obs)
    vspace(0)
    
    obs_exc = 0
    if obs.operator_idx[0][0]  in virtualall and obs.operator_idx[0][1] in occupied:
        obs_exc = 1
    elif obs.operator_idx[0][0]  in virtualall and obs.operator_idx[0][1] in virtualall:
        obs_exc = 0
    elif obs.operator_idx[0][0]  in occupied and obs.operator_idx[0][1] in occupied:
        obs_exc  = 0
    elif obs.operator_idx[0][0]  in occupied and obs.operator_idx[0][1] in virtualall:
        obs_exc = -1

    print('the excitation of the operator is', obs_exc, 'and should be', n)
    
    if obs_exc != n:
        wm1_outer = arithmetic_string()
        # print('nie bedzie', obs_exc, n)
        print('exit')
        return wm1_outer


    if lt1 == 0:
        print('wewnetrznego stringa e(-T) X e(T) nie ma, bo Gm1_T_list = []')
        wm1_inner = arithmetic_string(obs)
        # if wm1_inner is None:
        #     print('wm1_none')
        # else:
        #     print('wm1_inner', wm1_inner)
    else:
        vspace(0)
        print('wewnetrzny string e(-T) X e(T) jest, lt1=', lt1)
        vspace(0)
        wm1_inner = arithmetic_string(obs)
        print(wm1_inner)
        # if wm1_inner is None:
        #     print('wm1_none')
        # else:            
        #     print('wm1_inner1', wm1_inner, wm1t, lt1)
        
        for i in range(0, lt1):
            if wm1t[i] == 'tf':
                vspace(0)
                print('teraz amplituda T\'')
                vspace(0)
                wm1_inner = evaluate(wm1_inner, stoo(2, "tf", False))
            else:
                vspace(0)
                print('teraz amplituda zwykle T')
                vspace(0)
                print('rr', wm1_inner)
                print(stoo(wm1t[i], "t", False))
                wm1_inner = evaluate(wm1_inner, stoo(wm1t[i], "t", False))
#                print('rr', wm1_inner)

        vspace(0)
        # print('typ', type(wm1_inner))
        # if wm1_inner is None:
        #     print('wm1_none')
        # else:
        #     print('wm1_inner', wm1_inner)
        #     print('la')

    wm1_outer = wm1_inner
    

    if ls1 !=0:
        print('Gm2_S_list != []')

        for i in range(0, ls1):
            if wm1s[i] == 1:
                
                w_0 = deepcopy(wm1_outer)

                w_v1 = deepcopy(wm1_outer)
                w_v2 = deepcopy(wm1_outer)
                w_c1 = deepcopy(wm1_outer)
                w_c2 = deepcopy(wm1_outer)

                g_0 = stoo(wm1s[i], "t", True) # <-tu specjalnie jest t
                g_v1, g_v2 = stoo_f12(wm1s[i], "s_v", True)
                g_c1, g_c2 = stoo_f12(wm1s[i], "s_c", True)
                

                for x in w_0:
                    disambiguate(g_0, x)

                for x in w_v1:
                    disambiguate(g_v1, x)
                for x in w_v2:
                    disambiguate(g_v2, x)

                for x in w_c1:
                    disambiguate(g_c1, x)
                for x in w_c2:
                    disambiguate(g_c2, x)

                z_0 = evaluate(g_0, w_0)

                z_v1 = evaluate(g_v1, w_v1)
                z_v2 = evaluate(g_v2, w_v2)
                
                z_c1 = evaluate(g_c1, w_c1)
                z_c2 = evaluate(g_c2, w_c2)
                
                wm1_outer = z_0 + z_v1 + z_v2 + z_c1 + z_c2

            if wm1s[i] == 2:
                print('jestem tu 1')
                
                w_vv = deepcopy(wm1_outer)
                w_cc = deepcopy(wm1_outer)

                g_vv = stoo_f12(wm1s[i], "s_vv", True)
                g_cc = stoo_f12(wm1s[i], "s_cc", True)
                
                print('samo g_vv', g_vv)
                for x in w_vv:
                    print('x przed dis', x)
                    disambiguate(g_vv, x)
                    print('x', g_vv, x)

                print('samo g_cc', g_cc)
                for x in w_cc:
                    print('x przed dis', x)
                    disambiguate(g_cc, x)
                    print('x', g_cc, x)
                
                z_vv = evaluate(g_vv, w_vv)
                z_cc = evaluate(g_cc, w_cc)

                print('z_vv')
                if len(z_vv) > 0:
                    for x in z_vv:
                        print('out1', x)
                vspace(0)


                print('z_cc')
                if len(z_cc) > 0:
                    for x in z_cc:
                        print('out1', x)
                vspace(0)

                wm1_outer = z_vv + z_cc
                
        if wm1_outer is not None:
            for x in wm1_outer:
                print('wm1out', x)


        wm1_outer = wm1_outer.scale(n_wm1t)
        wm1_outer = wm1_outer.scale(n_wm1s)

    vspace(0)
    print('caly wynik to ')
    for x in wm1_outer:
        print(x)
    vspace(1)

    return wm1_outer


def density_ccd_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s, n):

    print('LICZE TERAZ DOKLADNIE', obs, wm1t, wm1s, lt1, ls1)
    obs_exc = 0
    if obs.operator_idx[0][0]  in virtual and obs.operator_idx[0][1] in occupied:
        obs_exc = 1
    elif obs.operator_idx[0][0]  in virtual and obs.operator_idx[0][1] in virtual:
        obs_exc = 0
    elif obs.operator_idx[0][0]  in occupied and obs.operator_idx[0][1] in occupied:
        obs_exc  = 0
    elif obs.operator_idx[0][0]  in occupied and obs.operator_idx[0][1] in virtual:
        obs_exc = -1

    
    print('OBS_EXC', wm1t, wm1s, obs_exc, n)
    # if obs_exc != n:
    #     wm1_outer = arithmetic_string()

    #     return wm1_outer
        
    print('obs_exc', obs_exc, n)
    wm1_outer = arithmetic_string()

    if lt1 == 0:
        print('la0')
        wm1_inner = arithmetic_string(obs)
    else:
        wm1_inner = arithmetic_string(obs)
        print('la1', wm1_inner)
        for i in range(0, lt1):
            wm1_inner = evaluate(wm1_inner, stoo(wm1t[i], "t", False))


    # print('hhhhhh', len(wm1_inner), wm1_inner)
    wm1_outer = wm1_inner
    # print( 'wm1_outererererere', wm1_outer)
    # for x in wm1_inner:
    #     print('wm1_inner', x)
    if ls1 !=0:
        for i in range(0, ls1):
            wm1_outer = evaluate(stoo(wm1s[i], "s", True), wm1_outer)
            # if i%2 != 0:
            #     wm1_outer = wm1_outer.scale(-1.0)


        wm1_outer = wm1_outer.scale(n_wm1t)
        wm1_outer = wm1_outer.scale(n_wm1s)

    print('wm1')
    for x in wm1_outer:
        print(x)
    print('-----------')
    return wm1_outer


def density_mbpt_inner_outer(obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s):

    if lt1 == 0:
        wm1_inner = arithmetic_string(obs)
    else:
        wm1_inner = arithmetic_string(obs)
        for i in range(0, lt1):
            wm1_inner = evaluate(wm1_inner, stoo(wm1t[i], "t", False))

    print('witam', obs, wm1t, wm1s, lt1, ls1, n_wm1t, n_wm1s)
    
    wm1_outer = wm1_inner
    if ls1 !=0:
        for i in range(0, ls1):
            wm1_outer = evaluate(stoo(wm1s[i], "s", True), wm1_outer)
            # print('laaaaaaaaaaakaaaaa', i, i%2)
            # if i%2 != 0:
            #     wm1_outer = wm1_outer.scale(-1.0)
            #     print('tak skaluje')
            # sss = stoo(wm1s[i], "s", False)
            # # for x in wm1_outer:
            # #     disambiguate(sss, x)
            # wm1_outer = evaluate(wm1_outer, sss)

        # nie trzeba zmieniac znaku, poniewaz jest zmieniona kolejnosc w drugim komutatorze
        # zgodnie z
        # e^(-T) X e^T = X + [X,T] + ...
        # e^(T) X e^(-T) = X - [X,T] + ... + = X + [T, X]+..
        wm1_outer = wm1_outer.scale(n_wm1t)
        wm1_outer = wm1_outer.scale(n_wm1s)
        
    # print('wm1')
    # if wm1t == [2, 2] and wm1s == [2, 2]:
    #     for x in wm1_outer:
    #         print(x)
    #     sys.exit(0)
    # print('outer outer', wm1_outer)
    return wm1_outer


def wm2_inner_outer(wm2_mu, wm2t, wm2s, lt2, ls2, n_wm2t, n_wm2s):

    if lt2 == 0:
        wm2_inner = arithmetic_string(wm2_mu)
    else:
        wm2_inner = arithmetic_string(wm2_mu)
        for i in range(0, lt2):
            wm2_inner = evaluate(wm2_inner, stoo(wm2t[i], "t", True))
            # sss = stoo(wm2t[i], "t", True)
            # print('sss przed', sss)
            # for x in wm2_inner:
            #     print('x przed', x)
            #     disambiguate(sss, x)
            #     print('x po', x)                
            # wm2_inner = evaluate(wm2_inner, sss)

    wm2_outer = wm2_inner
    if ls2 !=0:
        for i in range(0, ls2):
            wm2_outer = evaluate(stoo(wm2s[i], "s", False), wm2_outer)
            # sss = stoo(wm2s[i], "s", False)
            # # for x in wm2_outer:
            # #     disambiguate(sss, x)
            # wm2_outer = evaluate(wm2_outer, sss)

        # wm2_outer = wm2_outer.scale(-1.0)
        wm2_outer = wm2_outer.scale(n_wm2t)
        wm2_outer = wm2_outer.scale(n_wm2s)

    return wm2_outer


def commutators_quadra_Um(U_middle):

    U_middle_commutators = []
    Um1_commutators = []
    Um2_commutators = []
    Um3_commutators = []
    Um2_mu_list = []
    Um3_mu_list = []

    for x in range(0, len(U_middle)):

        um1t = U_middle[x]['Um1_T_list']
        um1s = U_middle[x]['Um1_S_list']

        um2t = U_middle[x]['Um2_T_list']
        um2s = U_middle[x]['Um2_S_list']
        um3s = U_middle[x]['Um3_S_list']
    
        lt1 = len(um1t)
        ls1 = len(um1s)

        lt2 = len(um2t)
        ls2 = len(um2s)

        ls3 = len(um3s)

        um2_mu = ugg()
        if (U_middle[x]['Um2_n']) == 1:
            um2_mu = mu1
            Um2_mu_list.append(1)
        elif (U_middle[x]['Um2_n']) == 2:
            um2_mu = mu2
            Um2_mu_list.append(2)
        elif (U_middle[x]['Um2_n']) == 3:
            um2_mu = mu3
            Um2_mu_list.append(3)

        um3_mu = ugg()
        eomr = ugg()
        if (U_middle[x]['Um3_n']) == 1:
            # um3_mu = mu1p
            Um3_mu_list.append(1)
            eomr = eomrr1
        elif (U_middle[x]['Um3_n']) == 2:
            # um3_mu = mu2p
            Um3_mu_list.append(2)
            eomr = eomrr2
        elif (U_middle[x]['Um3_n']) == 3:
            eomr = eomrr3
            Um3_mu_list.append(3)

    

        # Um1
        if lt1 == 0:
            um1_inner = arithmetic_string(obs)
        else:
            um1_inner = arithmetic_string(obs)
            for i in range(0, lt1):
                um1_inner = evaluate(um1_inner, stoo(um1t[i], "t", True))



        um1_outer = um1_inner
        if ls1 !=0:
            for i in range(0, ls1):
                um1_outer = evaluate(um1_outer, stoo(um1s[i], "s", False))

            um1_outer = um1_outer.scale(-1.0)
        Um1_commutators.append(um1_outer)


        # Um2
        if lt2 == 0:
            um2_inner = arithmetic_string(um2_mu)
        else:
            um2_inner = arithmetic_string(um2_mu)
            for i in range(0, lt2):
                um2_inner = evaluate(um2_inner, stoo(um2t[i], "t", True))


        um2_outer = um2_inner
        if ls2 !=0:
            for i in range(0, ls2):
                um2_outer = evaluate(um2_outer, stoo(um2s[i], "s", False))

            um2_outer = um2_outer.scale(-1.0)
        Um2_commutators.append(um2_outer)


        # Um3
        um3_outer = arithmetic_string(eomr)
        # um3_outer = arithmetic_string(um3_mu)
        if ls3 !=0:

            for i in range(0, ls3):
                um3_outer = evaluate(um3_outer, stoo(um3s[i], "s", True))

            um3_outer = um3_outer.scale(-1.0)
        Um3_commutators.append(um3_outer)


    return Um1_commutators,  Um2_commutators,  Um3_commutators, Um2_mu_list, Um3_mu_list


def simplify_W_middle(W_middle_int):

    # To simplify, first divide W_middle_int in parts:
    # W_middle_int1 with fixed ['a', 'i']
    # W_middle_int2 with fixed ['a', 'i', 'b', 'j']
    # W_middle_int3 with fixed ['a', 'i', 'b', 'j', 'c', 'k']


    # W_middle_int1 = arithmetic_string()
    # W_middle_int2 = arithmetic_string()
    # W_middle_int3 = arithmetic_string()


    # for i in range(0, len(W_middle_int)):
    #     W_middle_int[i].clear_fixed()
    #     W_middle_int[i].establish_fixed()
    #     if len(fixed) == 2:
    #         W_middle_int1.append(W_middle_int[i])
    #     elif len(fixed) == 4:
    #         W_middle_int2.append(W_middle_int[i])
    #     elif len(fixed) == 6:
    #         W_middle_int3.append(W_middle_int[i])

    # W_middle_int1simp = simplify(W_middle_int1)
    # W_middle_int2simp = simplify(W_middle_int2)
    # W_middle_int3simp = simplify(W_middle_int3)

    # W_middle_int_simp = W_middle_int1simp

    # for x in W_middle_int2simp:
    #     W_middle_int_simp.append(x)
    # for x in W_middle_int3simp:
    #     W_middle_int_simp.append(x)

    W_middle_int_simp = simplify(W_middle_int)
    return W_middle_int_simp

def integrate_quadra_Wlr(Wl_c, Wr_c, mu_list):



    ln = len(Wl_c)

    for i in range(0, ln):
        Wl_c[i].transpose()

    Wl_int = arithmetic_string()
    Wr_int = arithmetic_string()


    for i in range(0, ln):
        ln2 = len(Wl_c[i])

        for i2 in range(0, ln2):
        
            if mu_list[i] == 1:
                eoml = deepcopy(eoml1)
                print(Wl_c[i][i2])
                print(eoml)
                disambiguate(Wl_c[i][i2], eoml)

                r1l = Wl_c[i][i2].fromright(eoml)
                r1r = Wr_c[i][i2].fromleft(eoml)

                rintl = r1l.integrate().scale(0.5)
                rintr = r1r.integrate().scale(0.5)

                rintl.exec_delta()
                rintr.exec_delta()

                rsimpl = simplify(rintl)
                rsimpr = simplify(rintr)


            elif mu_list[i] == 2:


                rintla = Wl_c[i][i2].integrate(ket = ['a', 'i', 'b', 'j'], ketspin = ['s', 's'])
                rintlb = Wl_c[i][i2].integrate(ket = ['a', 'j', 'b', 'i'], ketspin = ['s', 's'])

                rintla.exec_delta()
                rintlb.exec_delta()

                rintl  = rintla.scale(1./3.) + rintlb.scale(1./6.)

                temp = deepcopy(rintl)
                for x in temp:
                    x.new_delta("a", "b")
                    x.new_delta("i", "j")
                rintl += temp.scale(-1./2.)

                rintl.exec_delta()

                rsimpl = simplify(rintl)

                rintra = Wr_c[i][i2].integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s', 's'])
                rintrb = Wr_c[i][i2].integrate(bra = ['a', 'j', 'b', 'i'], braspin = ['s', 's'])

                rintra.exec_delta()
                rintrb.exec_delta()

                rintr  = rintra.scale(1./3.) + rintrb.scale(1./6.)

                temp = deepcopy(rintr)
                for x in temp:
                    x.new_delta("a", "b")
                    x.new_delta("i", "j")
                rintr += temp.scale(-1./2.)

                rintr.exec_delta()

                rsimpr = simplify(rintr)

                for x in range(0, len(rsimpl)):
                    rsimpl[x].coefficient.append(EOM_CC_AMPLITUDE_L)
                    rsimpl[x].coefficient_idx.append(['a','i', 'b', 'j'])
                    rsimpl[x].summation.append('a')
                    rsimpl[x].summation.append('b')
                    rsimpl[x].summation.append('i')
                    rsimpl[x].summation.append('j')
                    rsimpl[x].num_factor *= 1./2.


                for x in range(0, len(rsimpr)):
                    rsimpr[x].coefficient.append(EOM_CC_AMPLITUDE_L)
                    rsimpr[x].coefficient_idx.append(['a','i', 'b', 'j'])
                    rsimpr[x].summation.append('a')
                    rsimpr[x].summation.append('b')
                    rsimpr[x].summation.append('i')
                    rsimpr[x].summation.append('j')
                    rsimpr[x].num_factor *= 1./2.
    

                rsimpl.exec_delta()
                rsimpl = simplify(rsimpl)
                rsimpr.exec_delta()
                rsimpr = simplify(rsimpr)

            elif mu_list[i] == 3:
                rintl = Wl_c[i][i2].integrate(ket = ['a', 'i', 'b', 'j', 'c', 'k'], ketspin = ['s', 's', 's'])
                rintl.exec_delta()
                rsimpl = simplify(rintl)

                rintr = Wr_c[i][i2].integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's'])
                rintr.exec_delta()
                rsimpr = simplify(rintr)

    
                for x in range(0, len(rsimpl)):
                    rsimpl[x].coefficient.append(EOM_CC_AMPLITUDE_L)
                    rsimpl[x].coefficient_idx.append(['a','i', 'b', 'j', 'c', 'k'])
                    rsimpl[x].summation.append('a')
                    rsimpl[x].summation.append('b')
                    rsimpl[x].summation.append('c')
                    rsimpl[x].summation.append('i')
                    rsimpl[x].summation.append('j')
                    rsimpl[x].summation.append('k')
                    rsimpl[x].num_factor *= 1./6.

                for x in range(0, len(rsimpr)):
                    rsimpr[x].coefficient.append(EOM_CC_AMPLITUDE_L)
                    rsimpr[x].coefficient_idx.append(['a','i', 'b', 'j', 'c', 'k'])
                    rsimpr[x].summation.append('a')
                    rsimpr[x].summation.append('b')
                    rsimpr[x].summation.append('c')
                    rsimpr[x].summation.append('i')
                    rsimpr[x].summation.append('j')
                    rsimpr[x].summation.append('k')
                    rsimpr[x].num_factor *= 1./6.

                rsimpl.exec_delta()
                rsimpl.clear_fixed()
                rsimpl = simplify(rsimpl)

                rsimpr.exec_delta()
                rsimpr.clear_fixed()
                rsimpr = simplify(rsimpr)

            for x in rsimpl:
                Wl_int.append(x)
            for x in rsimpl:
                Wr_int.append(x)

    Wl_int = simplify(Wl_int)
    Wr_int = simplify(Wr_int)

    return Wl_int, Wr_int

def integrate_quadra_Wm_overlap(Wm2_c, Wm3_c, Wm2_mu_list, Wm3_mu_list):
    
    
    ln = len(Wm2_c) # len(Wm1_c)=len(Wm2_c)=len(Wm3_c)

    wywalam = 0
    for i in range(0, ln):
        Wm2_c[i].transpose()
    sniez = 0
    W_middle_int = arithmetic_string()
    for i in range(0, ln):

        ln2 = len(Wm2_c[i])
        ln3 = len(Wm3_c[i])
        kk = ln2 * ln3
        
#        print('z wyrazu numer:', i, 'na', ln, 'bedzie do policzenia', kk, 'calek')
        k = 0
        for i2 in range(0, ln2):
            for i3 in range(0, ln3):
                k = k + 1
                # print('wm2', Wm2_c[i][i2])
                # print('wm3', Wm3_c[i][i3])
                # print('')
                Wm2copy = deepcopy(Wm2_c[i][i2])
                Wm3copy = deepcopy(Wm3_c[i][i3])
                disambiguate(Wm2copy, Wm3copy)

                # print('')
                r1 = Wm3copy
                r  = r1.fromleft(Wm2copy)

                if 1==1:
                    sniez += 1

                    print('wyraz' , i, 'na', ln, 'licze calke nr:', k, 'z', kk, '', len(r.operator_idx), '- operatorow')
                    # print('----------------------------------------------------------------------------')
                    # print(r)
                    # print('----------------------------------------------------------------------------')

                    start = time.time()
                    rint = r.integrate()
                    end = time.time()
                    # print('time elapsed=', end-start)

                    rint.exec_delta()
                    rsimp = simplify(rint)
#                    for x in rsimp:
#                    print('ZAPISUJE')
                    for x in rsimp:
 #                       print(x)
                        W_middle_int.append(x)
                else:
                    wywalam += 1
                    print('wywalam', r)
    
    print('wywalam', wywalam)
  #  print('przed simplify')
    # for x in W_middle_int:
    #     if len(x.summation) == 2:
    #         print(x)

    W_middle_int_simp = simplify(W_middle_int)

    print('simplify')
    for x in W_middle_int_simp:
        print(x)


    return W_middle_int_simp

def integrate_Gm(Gm1_c, Gm2_c, Gm1_middle, Gm2_middle, theory, mbpt):

    print('na wstepie gm')
    for x in Gm1_middle:
        print(x)
    print('')
    for x in Gm2_middle:
        print(x)
    print('++++++++++++++++++++')
    
    ln1 = len(Gm1_c)
    ln2 = len(Gm2_c)


    print('ln1', len(Gm1_c))
    print('ln2', len(Gm2_c))
    for x in Gm1_c:
        print(x)
    print('')
    for x in Gm2_c:
        print(x)
    print('')


    W_middle_int_1 = arithmetic_string() 
    W_middle_int_2 = arithmetic_string() 
    if theory == 'cc3':
        W_middle_int_3 =arithmetic_string()

    W_middle_int = arithmetic_string()

    integrate_list = []
    Gm2_mu_list_big = []


    zz = 0
    for j1 in range(0, ln1):
        for j2 in range(0, ln2):
            ln1s = len(Gm1_c[j1])
            ln2s = len(Gm2_c[j2])
            kk = ln1s * ln2s 
            print('gkgk', Gm1_middle[j1]['mbpt'], Gm2_middle[j2]['mbpt'])
            if (Gm1_middle[j1]['mbpt'] + Gm2_middle[j2]['mbpt']) == mbpt:
                k = 0
                for i1 in range(0, ln1s):
                    for i2 in range(0, ln2s):
                        k = k + 1
                        zz += 1
                        Gm1copy = deepcopy(Gm1_c[j1][i1])
                        Gm2copy = deepcopy(Gm2_c[j2][i2])
                        print('przed dis', Gm1copy)
                        print('przed dis', Gm2copy)
                                                
                        disambiguate(Gm1copy, Gm2copy)
                        print('po dis', Gm1copy)
                        print('po dis', Gm2copy)


                        r1 = Gm2copy
                        r0 = r1.fromleft(Gm1copy)
                        integrate_list.append(r0)

    print('tyle bedzie zadan', len(integrate_list))
    
    start = time.time()
    #
    #Tutaj nie ma całkowania w bazie biortonormalnej, bo z obu stron stoi wektor R  
    #
    z = Parallel(n_jobs=30,verbose=100)(delayed(integrate)(integrate_list[i]) for i in range(0, len(integrate_list)))

    rsimp_all = []
    for i in range(0, len(z)):
        rint = deepcopy(z[i])
        
        rint.exec_delta()
        rsimp = simplify(rint)
        for x in rsimp:
            rsimp_all.append(x)

    for x in rsimp_all:
        print('dupa1', x)
    print('')


    rsimp = arithmetic_string()
    for x in rsimp_all:
        rsimp.append(x)


    return rsimp


def integrate_Gm_one(Gm_c_big, oneel_exch_dict, oneel_delta_dict, Gm_middle, theory, mbpt):

    print('na wstepie gm')
    for x in Gm_middle:
        print(x)
    print('++++++++++++++++++++')
    
    ln1 = len(Gm_c_big)

    # print('ln1', len(Gm_c))
    print('do przecalkowania mam kazdy z kazdym')
    k = 1
    # print('las', Gm_c_big['iabj'][0][1], len(Gm_c_big['iabj'][0]))
    # print('las', Gm_c_big['iabj'][1][0], len(Gm_c_big['iabj'][1]))
    for key in Gm_c_big:
        Gm_c = Gm_c_big[key]
        print('przypadek', key, k, 'z 16')
        k = 0
        for x in Gm_c:
            k += 1
            k1 = 0
            for y in x[0]:
                k1 += 1
                print(k, k1, y)
            print('')
            k2 = 0
            for y in x[1]:
                k2 += 1
                print(k, k2, y)
            print('')
            print('-------------')
        print('')
    
    print('lla')
    W_middle_int = arithmetic_string() 
    Gm2_mu_list_big = []

    big_integrate_list = {}
    for key in Gm_c_big:
        Gm_c = Gm_c_big[key]
        integrate_list = []
        keylist = []
        for x in key:
            keylist.append(x)
        print('jestem tut',key,  len(Gm_c))
        for j1 in range(0, len(Gm_c)):
            for tt in Gm_c[j1]:
                if len(tt) > 0:
                    print('a', tt)
                else:
                    print('0')
            for i in range(0, len(Gm_c[j1][0])):
                print('iiii', i, Gm_c[j1][0][i], len(Gm_c[j1][1]))
                for j in range(0, len(Gm_c[j1][1])):
                    print( j1, i, j)
                    Gm1_temp = deepcopy(Gm_c[j1][0][i])
                    Gm2_temp = deepcopy(Gm_c[j1][1][j])
                    print('przedGm1_temp, Gm2_temp', Gm1_temp, Gm2_temp)
                    print(keylist)
                    disambiguate(Gm1_temp, Gm2_temp)
#                    disambiguate_two_list(Gm1_temp, Gm2_temp, keylist)
                    print('poGm1_temp, Gm2_temp', Gm1_temp, Gm2_temp)
                    r1 = Gm2_temp
                    r0 = r1.fromleft(Gm1_temp)
                    integrate_list.append(r0)
                    print('dodaje r0', key, r0)
            print('')
        big_integrate_list[key] = integrate_list

    print('')
    print('to jest lista do calkowania')
    print('')

    # print('keylist')

    # gr1 =  Gm_c_big['abcd'][0][0][0]
    # gr2 =  Gm_c_big['abcd'][0][1][0]
    # print('grgr')
    # print(gr1)
    # print(gr2)
    # e2 = ejbias
    # e1 = eaibjs
    # print('e1', e1)
    # print('e2', e2)
    # disambiguate_with_list(gr1, ['a', 'b', 'c', 'd'])
    # disambiguate_with_list(gr2, ['a', 'b', 'c', 'd'])
    # disambiguate_with_list(e1, ['a', 'b', 'c', 'd'])
    # disambiguate_with_list(e2, ['a', 'b', 'c', 'd'])
    # print('grgr')
    # print(gr1)
    # print(gr2)
    # print('e1', e1)
    # print('e2', e2)
    # disambiguate_two_list(e1, gr1, ['a', 'b', 'c', 'd'])
    # disambiguate_two_list(e2, gr2, ['a', 'b', 'c', 'd'])
    # list2 = ['a', 'b', 'c', 'd'] + e1.summation
    # print(list2)
    # disambiguate_two_list(gr1, gr2, list2)
    
    # print('grgr')
    # print(gr1, e1)
    # print(gr2, e2)

    # r1 = e1.fromleft(gr1)
    # r2 = gr2.fromleft(e2)
    # print('r1 - ab', r1)
    # print('r2 - cd', r2)
    # rint1  = r1.integrate(nodelta=True)
    # print('raw ab')
    # for x in rint1:
    #     print(x)
    # # print('exec')
    # # rint1.exec_delta()
    # # for x in rint1:
    # #     print(x)
    # # print('exec')
    # # rint1  =simplify(rint1)
    # # for x in rint1:
    # #     print(x)
    # print('-----------------')
    # rint2  = r2.integrate(nodelta=True)
    # print('raw cd')
    # for x in rint2:
    #     print(x)
    # # print('exec')
    # # rint2.exec_delta()
    # # for x in rint2:
    # #     print(x)
    # # print('exec')
    # # rint2  =simplify(rint2)
    # # for x in rint2:
    # #     print(x)
    # print('-----------------')
    # rint = arithmetic_string()
    # for x in rint1:
    #     for y in rint2:
    #         rint.append(y.fromleft(x))
    # print('poloczone')
    # for i in range(0, len(rint)):
    #     s1 = set(rint[i].summation)
    #     rint[i].summation = list(s1)
    #     rint[i].num_factor *= 0.5
    #     rint[i].remove_duplicate_summation
    #     print(rint[i])
    
    # rint.exec_delta()
    # print('exec')
    # for  x in rint:
    #     print(x)
    # rint  =simplify(rint)
    # print('rint')
    # for  x in rint:
    #     print(x)
#-------------------------------------------------------------------------

    # rint = arithmetic_string()
    # for x in rint1:
    #     for y in rint2:
    #         print('xyxy', x,y)
    #         disambiguate(x, y )
    #         rint.append(y.fromleft(x))
    # # rint = rint1 + rint2
    # print('oto jest rint')
    # for x in range(0, len(rint)):
    #     rint[x].exec_delta()
    #     rint[x].num_factor *= 0.5
    #     print(x, rint[x])
    # rsimp = simplify(rint)
    # print('rsimp')
    # for x in rsimp:
    #     print(x)
    # for x in range(0, len(rint2)):
    #     rint2[x].num_factor *= 0.5
    #     print(x)

    # print(gr1)
    # print(gr2)


#-----------------------------------------------------------------------------------
    # print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    # big_integrate_list['abcd'] = [big_integrate_list['abcd'][0]]
    # test0 = deepcopy(arithmetic_string(big_integrate_list['abcd'][0]))
    # print('test0', test0)
    # print('')
    # rint  = test0.integrate()
    # for x in rint:
    #     x.exec_delta()
    #     print(x)
    # sys.exit(0)
    
    for key in big_integrate_list:
        integrate_list = big_integrate_list[key]
        for y in range(0, len(integrate_list)):
            print(key, y, integrate_list[y])
        print('')
    print('to koniec tej listy')
    
    # for j1 in range(0, ln1):
    #     for j2 in range(0, ln2):
    #         ln1s = len(Gm1_c[j1])
    #         ln2s = len(Gm2_c[j2])
    #         kk = ln1s * ln2s 
    #         print('gkgk', Gm1_middle[j1]['mbpt'], Gm2_middle[j2]['mbpt'])
    #         if (Gm1_middle[j1]['mbpt'] + Gm2_middle[j2]['mbpt']) == mbpt:
    #             k = 0
    #             for i1 in range(0, ln1s):
    #                 for i2 in range(0, ln2s):
    #                     k = k + 1
    #                     zz += 1
    #                     Gm1copy = deepcopy(Gm1_c[j1][i1])
    #                     Gm2copy = deepcopy(Gm2_c[j2][i2])
    #                     print('przed dis', Gm1copy)
    #                     print('przed dis', Gm2copy)
                                                
    #                     disambiguate(Gm1copy, Gm2copy)
    #                     print('po dis', Gm1copy)
    #                     print('po dis', Gm2copy)


    #                     r1 = Gm2copy
    #                     r0 = r1.fromleft(Gm1copy)
    #                     integrate_list.append(r0)

    # print('tyle bedzie zadan', len(integrate_list))

    rsimp_big = {}
    for key  in big_integrate_list:
        print('CALKOWANIE DLA KEY', key)
        integrate_list = big_integrate_list[key]
        start = time.time()
        #
        #Tutaj nie ma całkowania w bazie biortonormalnej, bo z obu stron stoi wektor R  
        #

        for d in integrate_list:
            print(d)
        print('lne', len(integrate_list))
        z = Parallel(n_jobs=30,verbose=100)(delayed(integrate)(integrate_list[i]) for i in range(0, len(integrate_list)))

        print('bedzie', len(z), 'grup')

        rsimp_all = []
        for i in range(0, len(z)):
            rint = deepcopy(z[i])
            if len(rint) == 0:
                print('z tego calkowania wyszlo zero')
            else:
                print('')
                print('TO WYNIK CALKOWANIA')
                for k in z[i]:
                    print(k)
                print('')
                rint.exec_delta()
                rsimp = simplify(rint)
                for x in rsimp:
                    rsimp_all.append(x)

        # for x in rsimp_all:
        #     print('dupa1', x)
        # print('')


        rsimp = arithmetic_string()
        print('wynik całkowania nieuproszczonego dla ', key)
        for x in rsimp_all:
            print(x)
            rsimp.append(x)
        print('')
        # print(oneel_exch_dict)
        # print(oneel_delta_dict)

        rsimp = simplify(rsimp)
        print('wynik całkowania uproszczonego dla ', key)
        for x in rsimp:
            print(x)
        print('')

        disc = find_disconnected(rsimp)

        print('te wyrazy sa disconnected')
        for x in disc:
            print(x)

        print('')

        rs_exch = arithmetic_string()
        print('z macierzy exch wychodzi', key, oneel_exch_dict.keys())
        if key in oneel_exch_dict.keys():
            rsimp_oneel_exch = oneel_exch_dict[key]
            # print('oneel_exch_dict[key]', oneel_exch_dict[key])
            for x in rsimp_oneel_exch:
                print('oneeloneel', x)
                rsimp.append(x)
                rs_exch.append(x)
            print('')

        rss_exch = simplify(rs_exch)
        print('z macierzy exch simp wychodzi', key, oneel_exch_dict.keys())
        for x in rss_exch:
            print(x)
        print('')

        
        print('z macierzy delta wychodzi', key, oneel_delta_dict.keys())
        if key in oneel_delta_dict.keys():
            rsimp_oneel_delta = oneel_delta_dict[key]
            # print('oneel_exch_dict[key]', oneel_exch_dict[key])
            for x in rsimp_oneel_delta:
                print('oneeldel', x)
                rsimp.append(x)
            print('')

        rsimp_all = simplify(rsimp)
        # else:
        #     rsimp_all = rsimp

        print('a po uproszczeniu mamy', key)
        for x in rsimp_all:
            print(x)
        print('')
        print('-----------------------------------------------------------------------------------')

        disc = find_disconnected(rsimp_all)

        print('te wyrazy sa disconnected222')
        for x in disc:
            print(x)
        print('koniec wyrazow disconnected')

        print('')
        print('i to wpisuje do key', key)
        for x in rsimp_all:
            print(x)
        print('koniec wpisywania do key')
        rsimp_big[key] = rsimp_all

        print('lla')

    return rsimp_big

def find_disconnected(rsimp):

    disc = arithmetic_string()
    # for i in range(0, len(rsimp)):
    #     left_i = []
    #     for j in rsimp[i].coefficient_idx[0]:
    #         left_i.append(j)
    #     for j in rsimp[i].coefficient_idx[1]:
    #         left_i.append(j)
    #     right_i = []
    #     for j in rsimp[i].coefficient_idx[2]:
    #         right_i.append(j)
    #     for j in rsimp[i].coefficient_idx[3]:
    #         right_i.append(j)

    #     if len(list(set(left_i) & set(right_i))) == 0:
    #         disc.append(rsimp[i])

    return disc
        

def integrate_quadra_Wm(Wm1_c, Wm2_c, Wm3_c, Wm2_mu_list, Wm3_mu_list, theory):
    
    
    ln = len(Wm1_c) # len(Wm1_c)=len(Wm2_c)=len(Wm3_c)
    print(ln)


    print('l', len(Wm1_c))
    print('l2', len(Wm2_c))
    print('l3', len(Wm3_c))
    for x in Wm1_c:
        print(x)
    print('')
    for x in Wm2_c:
        print(x)
    print('')
    for x in Wm3_c:
        print(x)
    print('')


    wywalam = 0
    for i in range(0, ln):
        Wm1_c[i].transpose()
        Wm2_c[i].transpose()


    W_middle_int_1 = arithmetic_string() 
    W_middle_int_2 = arithmetic_string() 
    if theory == 'cc3':
        W_middle_int_3 =arithmetic_string()

    W_middle_int = arithmetic_string()

    integrate_list = []
    Wm2_mu_list_big = []
    Wm3_mu_list_big = []

    zz = 0
    for i in range(0, ln):
        ln1 = len(Wm1_c[i])
        ln2 = len(Wm2_c[i])
        ln3 = len(Wm3_c[i])
        kk = ln1 * ln2 * ln3
        
#        print('z wyrazu numer:', i, 'bedzie do policzenia', kk, 'calek', ln1, ln2, ln3)
        
        k = 0
        for i1 in range(0, ln1):
            for i2 in range(0, ln2):
                for i3 in range(0, ln3):
                    k = k + 1
                    zz += 1
                    Wm1copy = deepcopy(Wm1_c[i][i1])
                    Wm2copy = deepcopy(Wm2_c[i][i2])
                    Wm3copy = deepcopy(Wm3_c[i][i3])
                    # print('')
                    # print(zz, Wm1_c[i][i1], Wm2_c[i][i2], Wm3_c[i][i3])
                    disambiguate(Wm1copy, Wm2copy, Wm3copy)
                    # print(zz, 'po', Wm1copy, Wm2copy, Wm3copy)

                    r1 = Wm3copy
                    r0 = r1.fromleft(Wm1copy)
                    r  = r0.fromleft(Wm2copy)
                    Wm2_mu_list_big.append(Wm2_mu_list[i])
                    Wm3_mu_list_big.append(Wm3_mu_list[i])

                    integrate_list.append(r)
                    # print(zz, r)


    # print(len(integrate_list))
    # integrate_list = [integrate_list[147]]
#    sys.exit(0)


    # print('')
    # integrate_list = []
    # gg = flukt_potential
    # print(len(gg))
    # for t in range(0, len(gg)):
    #     disambiguate(gg[t], deepcopy(t2))
    #     ww = deepcopy(gg[t].fromright(t2))
    #     disambiguate(ww, deepcopy(t2))
    #     zz = deepcopy(ww.fromright(t2))
    #     disambiguate(zz, deepcopy(t2))
    #     uu = deepcopy(zz.fromleft(t2c))
    #     integrate_list.append(deepcopy(uu))
        
    # for t in range(0, len(integrate_list)):
    #     print(integrate_list[t])


#     integerate_list = []
#     ff = 1
#     for t in range(0, len(gg)):
#         disambiguate(t2, gg[t])
#         dd = gg[t].fromleft(t2)
#         disambiguate(t2, gg[t])                                                                                                                         
#         dd2 = dd.fromleft(t2) 
#         integrate_list.append(deepcopy(dd2))
#         print('dupa', ff, dd2)
#         ff += 1
#     print('teraz')
# #    sys.exit(0)


    # for i in range(0, len(integrate_list)):
    #     print(i)
    #     start = time.time()
    #     rint = integrate_list[i].integrate()
    #     # for x in rint:
    #     #     print('intint', x)
    #     end = time.time()
    #     print('time', start-end, k)

    # for x in rint:
    #     print(x)
    # rsimp = simplify(rint)
    # print('rsimp')
    # for x in rsimp:
    #     print(x)
    # # # print(Wm2_mu_list_big)
    # sys.exit(0)

    print('tyle bedzie zadan', len(integrate_list))
    
    start = time.time()
    #
    #Tutaj nie ma całkowania w bazie biortonormalnej, bo z obu stron stoi wektor R  
    #
    z = Parallel(n_jobs=30,verbose=100)(delayed(integrate)(integrate_list[i]) for i in range(0, len(integrate_list)))
    # for i in range(0, len(z)):
    #     print('w wyrazie', i)
    #     if len(z[i]) != 0:
    #         print(z[i])
    end = time.time()
    time_new = abs(end-start)
    print('czas caly', time_new)

    # for i in range(0, len(z)):
    #     rint = deepcopy(z[i])
    #     rint.exec_delta()
    #     rsimp = simplify(rint)
    #     print(rsimp)
    # sys.exit(0)
    # print('dupynik')
    # start = time.time()
    # for i in range(0, len(z)):
    #     rint = deepcopy(z[i])
    #     rint.exec_delta()
    #     rsimp = simplify(rint)
    #     for xx in rsimp:                                                                                                                               
    #         print(xx)   

    # print('koniec dupynika')
    # print(Wm2_mu_list_big)
#    sys.exit(0)

    for i in range(0, len(z)):
        rint = deepcopy(z[i])
        
        rint.exec_delta()
        rsimp = simplify(rint)

        if Wm2_mu_list_big[i] == 1:
            for x in rsimp:
                W_middle_int_1.append(deepcopy(x))
                print('1', x)
        elif Wm2_mu_list_big[i] == 2:
            for x in rsimp:
                W_middle_int_2.append(deepcopy(x))
                print('2', x)
        elif Wm2_mu_list_big[i] == 3:
            for x in rsimp:
                W_middle_int_3.append(deepcopy(x))
                print('3', x)


    end = time.time()
    time_new = abs(end-start)
    print('czas na simp 1', time_new)
    print('len wszystkich przed', len(W_middle_int_1), len(W_middle_int_2))

    for x in W_middle_int_1:
        print('dupa1', x)
    print('')
    for x in W_middle_int_2:
        print('dupa1', x)
    print('')

    for x in W_middle_int_3:
        print('dupa1', x)
    print('')


    if theory == 'ccsd':
        start = time.time()
        W_middle_int_1_simp = simplify_W_middle(W_middle_int_1)
        W_middle_int_2_simp = simplify_W_middle(W_middle_int_2)
        end = time.time()
        time_new = abs(end-start)
        print('czas simp 2', time_new)
        print('len wszystkich p0', len(W_middle_int_1_simp), len(W_middle_int_2_simp))
        return W_middle_int_1_simp, W_middle_int_2_simp

    elif theory == 'cc3':

        start = time.time()
        W_middle_int_1_simp = simplify_W_middle(W_middle_int_1)
        W_middle_int_2_simp = simplify_W_middle(W_middle_int_2)
        W_middle_int_3_simp = simplify_W_middle(W_middle_int_3)

        end = time.time()
        time_new = abs(end-start)
        print('czas simp 2', time_new)

        print('wmwmwm1')
        for x in W_middle_int_1_simp:
            print(x)
        print('wmwmwm2')
        for x in W_middle_int_2_simp:
            print(x)
        print('wmwmwm3')
        for x in W_middle_int_3_simp:
            print(x)


        return W_middle_int_1_simp, W_middle_int_2_simp, W_middle_int_3_simp

def integrate_quadra_Wm_triplet(Wm1_c, Wm2_c, Wm3_c, Wm2_mu_list, Wm3_mu_list, theory, multiplicity):
    
    ln = len(Wm1_c) # len(Wm1_c)=len(Wm2_c)=len(Wm3_c)
    #print(len(Wm1_c))
    #print(len(Wm2_c))
    #print(len(Wm3_c))

    for i in range(0, ln):
        Wm1_c[i].transpose()
        Wm2_c[i].transpose()

    W_middle_int_1 = arithmetic_string() 
    W_middle_int_2 = arithmetic_string() 
    W_middle_int_2_plus = arithmetic_string() 
    W_middle_int_2_minus = arithmetic_string() 
    if theory == 'cc3':
        W_middle_int_3 = arithmetic_string()

    # for i in range(0, ln):
    #     ln1 = len(Wm1_c[i])
    #     ln2 = len(Wm2_c[i])
    #     ln3 = len(Wm3_c[i])
    #     kk = ln1 * ln2 * ln3
    #     print('z wyrazu numer:', i, 'bedzie do policzenia', kk, 'calek', ln1, ln2, ln3)

    # sys.exit(0)

    # Prepare one list of all terms to integrate
    integrate_list = []
    Wm2_mu_list_big = []
    for i in range(0, ln):
        ln1 = len(Wm1_c[i])
        ln2 = len(Wm2_c[i])
        ln3 = len(Wm3_c[i])
        kk = ln1 * ln2 * ln3

        print('z wyrazu numer:', i, 'bedzie do policzenia', kk, 'calek')#, Wm1_c[i], Wm2_c[i], Wm3_c[i])                                                      
        k = 0
        time_old = 0
        for i1 in range(0, ln1):
            for i2 in range(0, ln2):
                for i3 in range(0, ln3):
                    k = k + 1

                    Wm1copy = deepcopy(Wm1_c[i][i1])
                    Wm2copy = deepcopy(Wm2_c[i][i2])
                    Wm3copy = deepcopy(Wm3_c[i][i3])
                    disambiguate(Wm1copy, Wm2copy, Wm3copy)

                    r1 = Wm3copy
                    r0 = r1.fromleft(Wm1copy)
                    r  = r0.fromleft(Wm2copy)
                    Wm2_mu_list_big.append(Wm2_mu_list[i])

                    integrate_list.append(r)

#     print('integrate_list')
#     start = time.time()
#     k = 0
#     for i in range(0, 10):#len(integrate_list)):#len(integrate_list)):
#         #print(i, integrate_list[i])
#         print(i)
#         rint = integrate_list[i].integrate()
#         for x in rint:
#             k+= 1
#             print(k, x)
# #   [integrate_list[i].integrate() for i in range(0, 10)]
#     end = time.time()
#     time_new = abs(end-start)
#     print('czas', time_new)
#     sys.exit(0)

    # print('integrate_list', len(integrate_list))
    # for i in range(0, len(integrate_list)):
    #     vidx = 0
    #     oidx = 0
    #     for j in integrate_list[i].operator_idx:
    #         if j[0] in virtual:
    #             vidx = vidx + 1
    #         elif j[0] in occupied:
    #             oidx = oidx + 1
    #         if j[1] in virtual:
    #             vidx = vidx + 1
    #         elif j[1] in occupied:
    #             oidx = oidx + 1
    #     print(i, integrate_list[i], len(integrate_list[i].operator_type), vidx, oidx)
    #     start = time.time()
    #     rint = integrate_list[i].integrate()
    #     end = time.time()
    #     print('time', start-end)

    # sys.exit(0)
        
    start = time.time()
#    Parallel(n_jobs=2)(delayed(integrate_list[i].integrate()) for i in range(0, 10))
 
    print('tyle bedzie zadan', len(integrate_list))
    z = Parallel(n_jobs=30,verbose=50)(delayed(integrate)(integrate_list[i]) for i in range(0, len(integrate_list)))    
    k = 0
    #for i in range(0, len(z)):
        #print('w wyrazie', i)
        #if len(z[i]) != 0:
         #   print(z[i])
        # for j in range(0, len(z[i])):
        #     k += 1
        #     print(k, i, j, z[i][j])
#    sys.exit(0)
    # for i in z:
    #     if (i != NoneType):
    #         print(i)
    # res, i = zip(*z)
    # for s in range(0, len(res)):
    #     print(i[s], res[s])
    # for i in range(0, 10):#len(integrate_list)):
    #     z = delayed(integrate)(integrate_list[i])
    #     print(z[1][0])
    #     print(i, integrate_list[i])
    #     rint = integrate_list[i].integrate()
    end = time.time()
    time_new = abs(end-start)
    print('czas', time_new)


    for i in range(0, len(z)):
        rint = deepcopy(z[i])

        rint.exec_delta()
        rsimp = simplify(rint)
        if Wm2_mu_list_big[i] == 1:
            for x in rsimp:
                W_middle_int_1.append(deepcopy(x))
        elif Wm2_mu_list_big[i] == '2p':
            for x in rsimp:
                W_middle_int_2_plus.append(deepcopy(x))
        elif Wm2_mu_list_big[i] == '2m':
            for x in rsimp:
                W_middle_int_2_minus.append(deepcopy(x))
        elif Wm2_mu_list_big[i] == 2:
            for x in rsimp:
                W_middle_int_2.append(deepcopy(x))
        elif Wm2_mu_list_big[i] == 3:
            for x in rsimp:
                W_middle_int_3.append(deepcopy(x))
        else:
            print(Wm2_mu_list_big[i])
            print('NIGDZIE')
            sys.exit(1)


    start = time.time()           
    W_middle_int_1_simp = simplify_W_middle(W_middle_int_1)
    W_middle_int_2_plus_simp = simplify_W_middle(W_middle_int_2_plus)
    W_middle_int_2_minus_simp = simplify_W_middle(W_middle_int_2_minus)
    W_middle_int_2_simp = simplify_W_middle(W_middle_int_2)
    end = time.time()
    time_new = abs(end-start)
    print('czas simp2', time_new)

    
    if theory == 'cc3':
        start = time.time()
        W_middle_int_3_simp = simplify_W_middle(W_middle_int_3)
        end = time.time()
        time_new = abs(end-start)
        print('czas simp2', time_new)


    if theory == 'ccsd':
        if multiplicity == 3:
            return  W_middle_int_1_simp, W_middle_int_2_plus_simp, W_middle_int_2_minus_simp
        elif multiplicity == 13:
            return  W_middle_int_1_simp, W_middle_int_2_simp
        elif multiplicity == 31:
            return  W_middle_int_1_simp, W_middle_int_2_plus_simp, W_middle_int_2_minus_simp
    elif theory == 'cc3':
        if multiplicity == 3:
            return  W_middle_int_1_simp, W_middle_int_2_plus_simp, W_middle_int_2_minus_simp, W_middle_int_3_simp
        elif multiplicity == 13:
            return  W_middle_int_1_simp, W_middle_int_2_simp, W_middle_int_3_simp
        elif multiplicity == 31:
            return  W_middle_int_1_simp, W_middle_int_2_plus_simp, W_middle_int_2_minus_simp, W_middle_int_3_simp

def integrate_quadra_Um(Um1_c, Um2_c, Um3_c, Um2_mu_list, Um3_mu_list):
    
    
    ln = len(Um1_c) # len(Wm1_c)=len(Wm2_c)=len(Wm3_c)

    for i in range(0, ln):
        Um3_c[i].transpose()
    

    U_middle_int_11 = arithmetic_string()
    U_middle_int_12 = arithmetic_string()
    U_middle_int_21 = arithmetic_string()
    U_middle_int_22 = arithmetic_string()

    U_middle_int = arithmetic_string()

    print(len(Um2_mu_list), len(Um3_mu_list))

    
    for i in range(0, ln):

        ln1 = len(Um1_c[i])
        ln2 = len(Um2_c[i])
        ln3 = len(Um3_c[i])
        kk = ln1 * ln2 * ln3
        print('z wyrazu numer:', i, 'bedzie do policzenia', kk, 'calek')
        k = 0
        for i1 in range(0, ln1):
            for i2 in range(0, ln2):
                for i3 in range(0, ln3):
                    k = k + 1
                    Um1copy = deepcopy(Um1_c[i][i1])
                    Um2copy = deepcopy(Um2_c[i][i2])
                    Um3copy = deepcopy(Um3_c[i][i3])
                    disambiguate(Um1copy, Um2copy, Um3copy)

                    
                    r1 = Um1copy
                    r0 = r1.fromleft(Um3copy)
                    r  = r0.fromright(Um2copy)
                    print('licze calke nr:', k, 'z', kk, '        ', len(r.operator_idx), '- operatorow')
                    print('----------------------------------------------------------------------------')
                    print(r)
                    print('----------------------------------------------------------------------------')

                    # rintr = r.integrate()
                    start = time.time()
                    rint = r.integrate()
                    end = time.time()
                    print('time elapsed=', end-start)

                    # # for x in rintr:
                    # #     print(x)
                    # # if (Wm2_mu_list[i] == 1):
                    # #     rint = rintr.fromleft(eomr1_amp)
                    # # elif (Wm2_mu_list[i] == 2):
                    # #     rint = rintr.fromleft(eomr2_amp)
                    # # elif (Wm2_mu_list[i] == 3):
                    # #     rint = rintr.fromleft(eomr3_amp)
                    # # for x in rint:
                    # #     print(x)

                    rint.exec_delta()
                    rsimp = simplify(rint)
                    print(len(rsimp), 'LEN(RSIMP)')

                    if Um2_mu_list[i] == 1 and Um3_mu_list[i] == 1:
                        print('')
                        print('DODAJE DO 11')
                        print('')
                        for x in rsimp:
                            U_middle_int_11.append(deepcopy(x))
                            print('11', x)
                    elif Um2_mu_list[i] == 1 and Um3_mu_list[i] == 2:
                        print('')
                        print('DODAJE DO 12')
                        print('')
                        for x in rsimp:
                            U_middle_int_12.append(deepcopy(x))
                            print('12', x)
                    elif Um2_mu_list[i] == 2 and Um3_mu_list[i] == 1:
                        print('')
                        print('DODAJE DO 21')
                        print('')
                        for x in rsimp:
                            U_middle_int_21.append(deepcopy(x))
                            print('21', x)
                    elif Um2_mu_list[i] == 2 and Um3_mu_list[i] == 2:
                        print('')
                        print('DODAJE DO 22')
                        print('')
                        for x in rsimp:
                            U_middle_int_22.append(deepcopy(x))
                            print('22', x)
                    else:
                        print('nie dodaje nigdzie')
                    print('')
                    for x in rsimp:
                        U_middle_int.append(x)
    

    U_middle_int_simp = simplify_W_middle(U_middle_int)

    U_middle_int_11_simp = simplify_W_middle(U_middle_int_11)
    U_middle_int_12_simp = simplify_W_middle(U_middle_int_12)
    U_middle_int_21_simp = simplify_W_middle(U_middle_int_21)
    U_middle_int_22_simp = simplify_W_middle(U_middle_int_22)

                        
    return U_middle_int_simp, U_middle_int_11_simp, U_middle_int_12_simp, U_middle_int_21_simp, U_middle_int_22_simp


def latex_W_middle_big(W_middle, Wm1_cc, Wm2_cc, Wm3_cc, maxmbpt, cumulative):       

    Wm1_c = deepcopy(Wm1_cc)
    Wm2_c = deepcopy(Wm2_cc)
    Wm3_c = deepcopy(Wm3_cc)
    

    ln = len(Wm1_c) # len(Wm1_c)=len(Wm2_c)=len(Wm3_c)


    for i in range(0, ln):
        
        Wm1_c[i].transpose()
        Wm2_c[i].transpose()

    n_int = []
    n_op = []

    for i in range(0, ln):
        ln1 = len(Wm1_c[i])
        ln2 = len(Wm2_c[i])
        ln3 = len(Wm3_c[i])
        kk = ln1 * ln2 * ln3

        # n_int.append(kk)
        # n_op.append(kk)

        n_int.append(1)
        n_op.append(1)

        
    latex_W_middle(W_middle, n_int, n_op, maxmbpt, cumulative)

def generate_descriptor(x):

    # [particle rank of S, particle rank of T]

    part_rank_t = 0
    part_rank_s = 0

    for i in range(0, len(x.coefficient)):
        if x.coefficient[i] == CC_AMPLITUDE:
            part_rank_t = len(x.coefficient_idx[i])//2
        if x.coefficient[i] == S_AMPLITUDE:
            part_rank_s = len(x.coefficient_idx[i])//2

    descriptor = [part_rank_s, part_rank_t]

    return descriptor

def generate_intermediate_fixed(interm_temp):

    idx_fx_sort = []
    for i in range(0, len(interm_temp.coefficient_idx)):
        for j in interm_temp.coefficient_idx[i]:
            if j not in interm_temp.summation:
                if j not in idx_fx_sort:
                    idx_fx_sort.append(j)

    return idx_fx_sort
    

def generate_intermediate_descriptor(coefa, coefb, interm_temp):
    
    int_desc = [[coefa, coefb]]
    for i in range(0, len(interm_temp.coefficient)):
        lst = []
        for j in interm_temp.coefficient_idx[i]:
            if j in interm_temp.summation:
                lst.append(1)
            else:
                lst.append(0)
        int_desc.append(lst)
    return int_desc

def simplify_intermediates(intermediates, interm_list, noninterm, name, k, multiplicity, mbpt):

    int_temp = []
    
    if intermediates != []:
        for x in intermediates:
            int_temp.append(x['interm'])

    # for x in interm_list:
    #     print(x['noninterm'])
    #     print(x['interm1'], x['interm2'], find_fx(x['interm1']), find_fx(x['interm2']))
    #     print('')

    # sys.exit(0)
    originals = []
    for i in range(0, len(interm_list)):        
        a = 0
        interm = interm_list[i]['interm1']
        if interm.num_factor != 0:
            a = 1
            if interm not in int_temp:
                if multiplicity == 1:
                    int_name = "{name}_interm_{k}_pt{mbpt}".format(name=name, k=k, mbpt=mbpt)
                elif multiplicity ==3:
                    int_name = "{name}_interm_{k}_triplet_pt{mbpt}".format(name=name, k=k, mbpt=mbpt)
                elif multiplicity == 13:
                    int_name = "{name}_interm_{k}_so_pt{mbpt}".format(name=name, k=k, mbpt=mbpt)
                elif multiplicity == 31:
                    int_name = "{name}_interm_{k}_so_left_pt{mbpt}".format(name=name, k=k, mbpt=mbpt)

                k += 1
                minidict = {}
                minidict['interm'] = interm
                minidict['int_name'] = int_name
                idx_fx = find_fx(interm)
                minidict['idx_fx'] = idx_fx
                intermediates.append(minidict)
                int_temp.append(interm)
                print(i, 'interm', interm)
                print(int_name, idx_fx)
            else:
                l = int_temp.index(interm)
                int_name = intermediates[l]['int_name']

            for j in range(0, len(interm_list[i]['noninterm'].coefficient)):
                if interm_list[i]['noninterm'].coefficient[j] == 'Q0' or \
                   interm_list[i]['noninterm'].coefficient[j] == 'Q1':
                    interm_list[i]['noninterm'].coefficient[j] = deepcopy(int_name)
                    # print('nadaje nazwe 1', deepcopy(int_name))
            print(interm_list[i]['noninterm'])

        interm = interm_list[i]['interm2']
        if interm.num_factor != 0:
            a = 1
            if interm not in int_temp:
                if multiplicity == 1:
                    int_name = "{name}_interm_{k}_pt{mbpt}".format(name=name, k=k, mbpt=mbpt)
                elif multiplicity == 3:
                    int_name = "{name}_interm_{k}_triplet_pt{mbpt}".format(name=name, k=k, mbpt=mbpt)
                elif multiplicity == 13:
                    int_name = "{name}_interm_{k}_so_pt{mbpt}".format(name=name, k=k, mbpt=mbpt)
                elif multiplicity == 31:
                    int_name = "{name}_interm_{k}_so_left_pt{mbpt}".format(name=name, k=k, mbpt=mbpt)
                k += 1
                minidict = {}
                minidict['interm'] = interm
                minidict['int_name'] = int_name
                idx_fx = find_fx(interm)
                minidict['idx_fx'] = idx_fx
                intermediates.append(minidict)
                int_temp.append(interm)
                print(i, 'interm', interm)
                print(int_name, idx_fx)
            else:
                l = int_temp.index(interm)
                int_name = intermediates[l]['int_name']

            for j in range(0, len(interm_list[i]['noninterm'].coefficient)):
                if interm_list[i]['noninterm'].coefficient[j] == 'Q2':
                    interm_list[i]['noninterm'].coefficient[j] = deepcopy(int_name)
            print(interm_list[i]['noninterm'])
                    # print('nadaje nazwe 2', deepcopy(int_name))

        interm = interm_list[i]['interm3']
        if interm.num_factor != 0:
            a = 1
            if interm not in int_temp:
                if multiplicity == 1:
                    int_name = "{name}_interm_{k}_pt{mbpt}".format(name=name, k=k, mbpt=mbpt)
                elif multiplicity ==3:
                    int_name = "{name}_interm_{k}_triplet_pt{mbpt}".format(name=name, k=k, mbpt=mbpt)
                elif multiplicity == 13:
                    int_name = "{name}_interm_{k}_so_pt{mbpt}".format(name=name, k=k, mbpt=mbpt)
                elif multiplicity == 31:
                    int_name = "{name}_interm_{k}_so_left_pt{mbpt}".format(name=name, k=k, mbpt=mbpt)

                k += 1
                minidict = {}
                minidict['interm'] = interm
                minidict['int_name'] = int_name
                idx_fx = find_fx(interm)
                minidict['idx_fx'] = idx_fx
                intermediates.append(minidict)
                int_temp.append(interm)
                print(i, 'interm', interm)
                print(int_name, idx_fx)
            else:
                l = int_temp.index(interm)
                int_name = intermediates[l]['int_name']

            for j in range(0, len(interm_list[i]['noninterm'].coefficient)):
                if interm_list[i]['noninterm'].coefficient[j] == 'Q0' or \
                   interm_list[i]['noninterm'].coefficient[j] == 'Q1':
                    interm_list[i]['noninterm'].coefficient[j] = deepcopy(int_name)
                    # print('nadaje nazwe 1', deepcopy(int_name))                                                                                        
            print(interm_list[i]['noninterm'])



        if a == 1:
            originals.append(interm_list[i]['noninterm'])

    for x in noninterm:
        originals.append(x)

    return intermediates, originals, k

                

def generate_best_intermediates(Wm_int, name, idx_start):

    Wm_intermediates_list = []
    
    
    print('intermediates  tuz przed podzialem na intermediates')
    for x in Wm_int:
        print(x)

    print('')

    intermediates = []
    nonintermediates = []
    print('generate_best_intermediates')

    sing = 0
    doub = 1
    trip = 2
    

    for i in range(0, len(Wm_int)):
        
        Wm_int[i].clear_fixed()

        minidict = {}
        print('szukam intermediate dla',i,  Wm_int[i])
        print('')

        if len(Wm_int[i].coefficient) < 4:
            # <4 szukam pojedynczego i potrojnego

            add_single, interm1, noninterm1 = compute_single_intermediate(deepcopy(Wm_int[i]))
            # interm2 = ugg()
            # interm2.num_factor = 0
            # add_triple, interm3, noninterm3 = compute_triple_intermediates(deepcopy(Wm_int[i]))
            
            # add_number = select_best_intermediates([[interm1, noninterm1], [interm3, noninterm3]], [add_single, add_triple])
            # if add_number == 1:
            #     #To jest po to zeby interm1 mial 0, interm2 mial 1 a interm3 mial 2. A tutaj interm3 ma 1.
            #     add_number == 2

            if add_single == True:
                add_number = 0
            else:
                add_number = -1
            
        elif len(Wm_int[i].coefficient) >= 4:
            # >=4 szukam pojedynczego, podwojnego i potrojnego i patrze ktory najlepszy
            add_single, interm1s, noninterms = compute_single_intermediate(deepcopy(Wm_int[i]))
            add_double, interm1d, interm2d, nonintermd = compute_double_intermediate(deepcopy(Wm_int[i]))
            add_triple, interm3t, nonintermt = compute_triple_intermediates(deepcopy(Wm_int[i]))


            # if add_double == True or add_triple == True:
            #     print('truuu')
            #     print(interm1s, noninterms)
            #     print(interm1d, interm2d, nonintermd)
            #     print(interm3t, nonintermt)
            #     add_number = select_best_intermediates([[interm1s, noninterms], [interm1d, interm2d, nonintermd], [interm3t, nonintermt]], \
            #                                                [add_single, add_double, add_triple])
            #     print(add_number)
            #     sys.exit(0)
            add_number = select_best_intermediates([[interm1s, noninterms], [interm1d, interm2d, nonintermd], [interm3t, nonintermt]], \
                                                       [add_single, add_double, add_triple])


        if add_number == sing:

            interm1s.standarize_interm()
            intermediates.append(interm1s)
            
            interm2 = ugg()
            interm2.num_factor = 0
            
            interm3 = ugg()
            interm3.num_factor = 0

            minidict['interm1'] = interm1s
            minidict['interm2'] = interm2
            minidict['interm3'] = interm3
            minidict['noninterm'] = noninterms
            minidict['original'] = Wm_int[i]

            # z = compute_cost_direct_ugg(noninterms)
            # if (z == [4, 4]):
            #     print( minidict['interm1'],  minidict['noninterm'],  minidict['original'])
            #     print(add_single, add_double, add_triple)
            #     sys.exit(0)
            Wm_intermediates_list.append(minidict)
        elif add_number == doub:

            interm1d.standarize_interm()
            intermediates.append(interm1d)
                
            interm2d.standarize_interm()
            intermediates.append(interm2d)

            interm3 = ugg()
            interm3.num_factor = 0

                
            minidict['interm1'] = interm1d
            minidict['interm2'] = interm2d
            minidict['interm3'] = interm3
            minidict['noninterm'] = nonintermd
            minidict['original'] = Wm_int[i]
            Wm_intermediates_list.append(minidict)
        elif add_number == trip:

            interm1 = ugg()
            interm1.num_factor = 0

            interm2 = ugg()
            interm2.num_factor = 0

            minidict['interm1'] = interm1
            minidict['interm2'] = interm2
            minidict['interm3'] = interm3t
            minidict['noninterm'] = deepcopy(nonintermt)
            minidict['original'] = Wm_int[i]
            Wm_intermediates_list.append(minidict)
        elif add_number == -1:

            nonintermediates.append(Wm_int[i])
        else:
            print('add_number error = ', add_number)
            sys.exit(0)

    print('wm_list')
    for x in Wm_intermediates_list:
        print(x)
    print('nonintermediates')
    for x in nonintermediates:
        print(x)

    return Wm_intermediates_list, nonintermediates

def find_mem_int(interm):

    fx = []
    for i in range(0, len(interm.coefficient_idx)):
        for j in interm.coefficient_idx[i]:
            if j not in interm.summation and j not in fx:
                fx.append(j)

    vt = 0
    oc = 0
    for j in fx:
        if j in virtual:
            vt += 1
        if j in occupied:
            oc += 1

    return vt, oc

def find_fx(interm):

    fx = []

    fx_v = []
    fx_o = []

    for i in range(0, len(interm.coefficient_idx)):
        for j in interm.coefficient_idx[i]:
            if j not in interm.summation and j not in fx:
                if j in virtual:
                    fx_v.append(j)
                elif j in occupied:
                    fx_o.append(j)

    fx = fx_v + fx_o

    return fx
  

def compute_triple_intermediates(Wm_int):

    best_interm3 = ugg()
    best_interm3.num_factor = 0
    best_noninterm3 = ugg()
    best_noninterm3.num_factor = 0
    vt_org, oc_org = comp_cost_sum(Wm_int)
    max_cost_old = [vt_org, oc_org]
    fx_cost_old = [0, 0]
    add_triple = False

    coef_idx = []
    for j in range(0, len(Wm_int.coefficient)):
        if Wm_int.coefficient[j] != OBSERVABLE_X and Wm_int.coefficient[j] != OBSERVABLE_X_ASYM:
            coef_idx.append(j)
        
    coef_triples = coef_combinations_triples(coef_idx)

    all_check = False
    for j in range(0, len(coef_triples)):
        c_idx_11 = Wm_int.coefficient_idx[coef_triples[j][0]]
        c_idx_12 = Wm_int.coefficient_idx[coef_triples[j][1]]
        c_idx_13 = Wm_int.coefficient_idx[coef_triples[j][2]]

        idx_sum, idx_fx = find_sum_fx([c_idx_11, c_idx_12, c_idx_13], coef_triples[j], Wm_int)
        interm3, noninterm3 = divide_triple_interm(Wm_int, coef_triples[j], idx_sum, idx_fx)


        print('interm3', interm3)
        print('nonint3', noninterm3)

        vt1, oc1 = comp_cost_sum_fx(interm3)
        vtn, ocn = comp_cost_sum(noninterm3)
        vt_org, oc_org = comp_cost_sum(Wm_int)

        check2 = check_fx(idx_fx)

        vt, ot = find_mem_int(interm3)
        if (vt + ot == 0):
            print('FIND MEM ZEROOOO')
            print(interm3, noninterm3, Wm_int)
#            sys.exit(0)

        if len(idx_sum[0])+ len(idx_sum[1])==0:
            check3 = False
            print('nie dodaje1')
        else:
            check3 = True            

        if check2 == True and check3 == True:
            check = True
        else:
            check = False
            print('nie dodaje2')

        if check == True:
            max_cost_new = find_expensive([interm3, noninterm3])
            max_mem_new = find_mem([interm3, noninterm3])
            dodaje = False
            if (max_cost_new[0] < max_cost_old[0]):
                dodaje = True
                max_mem_old = [max_mem_new[0], max_mem_new[1]]
                max_cost_old = max_cost_new
                best_interm3 = deepcopy(interm3)
                best_noninterm3 = deepcopy(noninterm3)
                all_check = True
            elif(max_cost_new[0] == max_cost_old[0]):
                if(max_cost_new[1] < max_cost_old[1]):
                    dodaje = True
                    max_mem_old = [max_mem_new[0], max_mem_new[1]]
                    max_cost_old = max_cost_new
                    best_interm3 = deepcopy(interm3)
                    best_noninterm3= deepcopy(noninterm3)
                    all_check = True
                else:
                    if (max_mem_new[0] < max_mem_old[0]):
                        dodaje = True
                        max_mem_old = [max_mem_new[0], max_mem_new[1]]
                        max_cost_old = max_cost_new
                        best_interm3 = deepcopy(interm3)
                        best_noninterm3= deepcopy(noninterm3)
                        all_check = True
                    elif (max_mem_new[0] == max_mem_old[0]):
                        if (max_mem_new[1] < max_mem_old[1]):
                            dodaje = True
                            max_mem_old = [max_mem_new[0], max_mem_new[1]]
                            max_cost_old = max_cost_new
                            best_interm3 = deepcopy(interm3)
                            best_noninterm3= deepcopy(noninterm3)
                            all_check = True
                        else:
                            dodaje = False
                    else:
                        dodaje = False
            else:
                dodaje = False

            if dodaje == True:
                #                    print('dodaje NOWY')
                best_interm3 = deepcopy(interm3)
                best_noninterm3 = deepcopy(noninterm3)
                add_triple = True
            else:
                print('nie dodaje ABCD')
        else:
            print('nie dodaje check-false')


        
    return(add_triple, best_interm3, best_noninterm3)
            

def select_best_intermediates(int_lists, add_list):

    max_cost_list = []
    max_mem_list = []
    cost_list = []
    A = 10000
    B = 1000
    C = 100
    D = 10
    max_cost = 10**8+1

    print('-------------------')
    print(add_list)
    for x in int_lists:
        for y in x:
            print(y)
        print('')
    print('-------------------')

    for i in range(0, len(int_lists)):
        
        x = int_lists[i]
        print('x',i, x)

        if add_list[i] == False:
            cost_list.append(max_cost)
        else:
            mc = find_expensive(x)
            mm = find_mem(x)
            # max_cost_list.append(mc)
            # max_mem.append(mm)
            cost_list.append(A*mc[0] + B*mc[1] + C*mm[0] + D*mm[1])

    print(cost_list)
    add_number = cost_list.index(min(cost_list))
    if add_number >= max_cost:
        add_number = -1

    if add_number == 0:
        print('tak, add zero')
        if cost_list[0] == cost_list[1]:
            print('tak rowne')
            add_number = 1

    print('COSTCOST', cost_list)

    return add_number
        
def compute_double_intermediate(Wm_int):

    best_interm1 = ugg()
    best_interm1.num_factor = 0
    best_interm2 = ugg()
    best_interm2.num_factor = 0
    best_noninterm = ugg()
    best_noninterm.num_factor = 0

    vt_org, oc_org = comp_cost_sum(Wm_int)
    max_cost_old = [vt_org, oc_org]
    fx_cost_old = [0, 0]
    add_double = False
    coef_idx = []
    for j in range(0, len(Wm_int.coefficient)):
        if Wm_int.coefficient[j] != OBSERVABLE_X and Wm_int.coefficient[j] != OBSERVABLE_X_ASYM:
            coef_idx.append(j)


    coef_pairs = coef_combinations(coef_idx)

    for j in range(0, len(coef_pairs)):
        c_idx_11 = Wm_int.coefficient_idx[coef_pairs[j][0][0]]
        c_idx_12 = Wm_int.coefficient_idx[coef_pairs[j][0][1]]
        c_idx_21 = Wm_int.coefficient_idx[coef_pairs[j][1][0]]
        c_idx_22 = Wm_int.coefficient_idx[coef_pairs[j][1][1]]


        idx_sum_1, idx_fx_1 = find_sum_fx([c_idx_11, c_idx_12], coef_pairs[j][0], Wm_int)
        idx_sum_2, idx_fx_2 = find_sum_fx([c_idx_21, c_idx_22], coef_pairs[j][1], Wm_int)

        interm1, interm2, noninterm = divide_double_interm(Wm_int, coef_pairs[j], idx_sum_1, idx_sum_2, idx_fx_1, idx_fx_2)
        vt1, oc1 = comp_cost_sum_fx(interm1)
        vt2, oc2 = comp_cost_sum_fx(interm2)
        vt3, oc3 = comp_cost_sum(noninterm)
        vt_org, oc_org = comp_cost_sum(Wm_int)

        print(interm1)
        print(interm2)
        print(noninterm)
        print('')

            
        check2a = check_fx(idx_fx_1)
        check2b = check_fx(idx_fx_2)

        if check2a == True and check2b == True:
            check2 = True
        else:
            check2 = False
            print('compute_double check2 False')

        if len(idx_sum_1[0])+ len(idx_sum_1[1])==0 or len(idx_sum_2[0]) + len(idx_sum_2[1]) ==0:
            check3 = False
            print('compute_double check3 False')
        else:
            check3 = True

        if check2 == True and check3 == True:
            check = True
        else:
            check = False
            print('compute_double check False')


            
        if check == True:
            max_cost_new = find_expensive([interm1, interm2, noninterm])
            max_mem_new = find_mem([interm1, interm2, noninterm])
            dodaje = False
            if (max_cost_new[0] < max_cost_old[0]):
                dodaje = True
                max_mem_old = [max_mem_new[0], max_mem_new[1]]
                max_cost_old = max_cost_new
            elif(max_cost_new[0] == max_cost_old[0]):
                if(max_cost_new[1] < max_cost_old[1]):
                    dodaje = True
                    max_mem_old = [max_mem_new[0], max_mem_new[1]]
                    max_cost_old = max_cost_new
                else:
                    if (max_mem_new[0] < max_mem_old[0]):
                        dodaje = True
                        max_mem_old = [max_mem_new[0], max_mem_new[1]]
                        max_cost_old = max_cost_new
                    elif (max_mem_new[0] == max_mem_old[0]):
                        if (max_mem_new[1] < max_mem_old[1]):
                            dodaje = True
                            max_mem_old = [max_mem_new[0], max_mem_new[1]]
                            max_cost_old = max_cost_new
                        else:
                            dodaje = False
                            print('compute_double dodaje1 False')
                    else:
                        dodaje = False
                        print('compute_double DODAJE2 False')
            else:
                dodaje = False
                print('compute_double dodaje3 False')

            if dodaje == True:
                print('dodaje NOWY')
                best_interm1 = deepcopy(interm1)
                best_interm2 = deepcopy(interm2)
                best_noninterm = deepcopy(noninterm)
                add_double = True
            else:
                print('nie dodaje1')
                print('compute_double nieee')
        else:
            print('nie dodaje check-false333')


    return add_double, best_interm1, best_interm2, best_noninterm

def find_mem(ug_list):
    v_list = []
    o_list = []

    for x in ug_list:
        v_listl = []
        o_listl = []
        for y in x.coefficient_idx:
            for z in y:
                if z not in x.summation:
                    if z in virtual:
                        v_listl.append(z)
                    elif z in occupied:
                        o_listl.append(z)
        v_list.append(len(v_listl))
        o_list.append(len(o_listl))

    idx_max = 0
    max_mem = [v_list[0], o_list[0]]
    for i in range(1, len(ug_list)):
        if v_list[i] > v_list[idx_max]:
            idx_max = i
            max_mem = [v_list[i], o_list[i]]
        elif v_list[i] == v_list[idx_max]:
            if o_list[i] >= o_list[idx_max]:
                idx_max = i
                max_mem = [v_list[i], o_list[i]]

    return max_mem



def find_expensive(ug_list):

    v_list = []
    o_list = []

    for x in ug_list:
        v_listl = []
        o_listl = []
        for y in x.coefficient_idx:
            for z in y:
                if z not in x.summation:
                    if z in virtual:
                        v_listl.append(z)
                    elif z in occupied:
                        o_listl.append(z)
        v_list.append(len(v_listl))
        o_list.append(len(o_listl))

    
    idx_max = 0
    max_cost = [v_list[0], o_list[0]]
    for i in range(1, len(ug_list)-1):
        if v_list[i] > v_list[idx_max]:
            idx_max = i
            max_cost = [v_list[i], o_list[i]]
        elif v_list[i] == v_list[idx_max]:
            if o_list[i] >= o_list[idx_max]:
                idx_max = i
                max_cost = [v_list[i], o_list[i]]
    return max_cost

def compute_single_intermediate(Wm_int):

    best_interm = ugg()
    best_interm.num_factor = 0
    best_noninterm = ugg()
    best_noninterm.num_factor = 0

    vt, oc = comp_cost_sum(Wm_int)

    coef_idx = []
    for j in range(0, len(Wm_int.coefficient)):
        if Wm_int.coefficient[j] != OBSERVABLE_X and Wm_int.coefficient[j] != OBSERVABLE_X_ASYM:
            coef_idx.append(j)


    coef_pairs = list(combinations(coef_idx, 2))

    vt_int_best = 100
    oc_int_best = 100
    
    vt_nint_best = 100
    oc_nint_best = 100
    
    vt_mem_best = 100
    oc_mem_best = 100
    
#    print('LENNNC COEF', len(coef_pairs))
    for j in range(0, len(coef_pairs)):
 #       print('jjjjjjjjjjjjjjj', j, len(coef_pairs), '-------------------------------------------')
  #      print('')
        c_idx_1 = Wm_int.coefficient_idx[coef_pairs[j][0]]
        c_idx_2 = Wm_int.coefficient_idx[coef_pairs[j][1]]
   #     print(c_idx_1, c_idx_2)

        idx_sum, idx_fx = find_sum_fx([c_idx_1, c_idx_2], coef_pairs[j], Wm_int)

        interm_new, noninterm_new = divide_interm(deepcopy(Wm_int), idx_sum, idx_fx, coef_pairs[j])
    #    print('interm, noninterm', interm_new, noninterm_new)
     #   print('')

        vt_int, oc_int = comp_cost_sum_fx(interm_new)
        vt_nint, oc_nint = comp_cost_sum(noninterm_new)
        vt_org, oc_org = comp_cost_sum(Wm_int)
        # print('cost interm', vt_int, oc_int, 'cost noninterm', vt_nint, oc_nint, 'cost original', vt_org, oc_org)
        # print('')
        # print('compare cost',vt_int, oc_int, vt_org, oc_org)
        # print('compare cost',vt_nint, oc_nint, vt_org, oc_org)
        # print('')
        check_cost_int = compare_cost(vt_int, oc_int, vt_org, oc_org)
        check_cost_nint = compare_cost(vt_nint, oc_nint, vt_org, oc_org)
        
        if check_cost_int == True and check_cost_nint == True:
            check_cost = True
        else:
            check_cost = False

        check_mem = check_fx(idx_fx)
        vt_mem = len(idx_fx[0])
        oc_mem = len(idx_fx[1])

      #  print('check best')
        check_best = compare_interm(vt_int, oc_int, vt_nint, oc_nint, vt_mem, oc_mem, \
                                        vt_int_best, oc_int_best, vt_nint_best, oc_nint_best, vt_mem_best, oc_mem_best)

        check_all = False

        if check_cost == True and check_mem == True and check_best == True:
            check_all = True

        if check_all == True:
            best_interm = deepcopy(interm_new)
            best_noninterm = deepcopy(noninterm_new)
            vt_int_best = deepcopy(vt_int)
            oc_int_best = deepcopy(oc_int)
            
            vt_nint_best = deepcopy(vt_nint)
            oc_nint_best = deepcopy(oc_nint)

            vt_mem_best = deepcopy(vt_mem)
            oc_mem_best = deepcopy(oc_mem)

            
            print('dodaje', '-----------------------------------------------------------------------')
        else:
            print('nie dodaje', '----------------------------------------------------------------------')

    if best_interm.num_factor == 0 and best_noninterm.num_factor == 0:
        add_single = False
    else:
        add_single = True

    return add_single, best_interm, best_noninterm

def compare_interm(vt_int, oc_int, vt_nint, oc_nint, vt_mem, oc_mem, \
                       vt_int_best, oc_int_best, vt_nint_best, oc_nint_best, vt_mem_best, oc_mem_best):

    # sprawdz ktory wyraz jest bardziej kosztowny w nowym intermediate
    # i w starym intermediate
#    print(vt_int, oc_int, vt_nint, oc_nint)
#    print(vt_int_best, oc_int_best, vt_nint_best, oc_nint_best)
    if vt_int > vt_nint:
#        print('a1')
        new = 'int'
    elif vt_int == vt_nint:
        if oc_int > oc_nint:
 #           print('a2')
            new = 'int'
        else:
  #          print('a3')
            new = 'nint'
    elif vt_int < vt_nint:
   #     print('a4')
        new = 'nint'
        
    if vt_int_best > vt_nint_best:
        old = 'int'
    elif vt_int_best == vt_nint_best:
        if oc_int_best > oc_nint_best:
            old = 'int'
        else:
            old = 'nint'
    elif vt_int_best < vt_nint_best:
        old = 'nint'

        
    if new == 'int' and old == 'int':
        new_list = [vt_int, oc_int, vt_nint, oc_nint, vt_mem, oc_mem]
        old_list = [vt_int_best, oc_int_best, vt_nint_best, oc_nint_best, vt_mem_best, oc_mem_best]
    if new == 'nint' and old == 'int':
        new_list = [vt_nint, oc_nint, vt_int, oc_int, vt_mem, oc_mem]
        old_list = [vt_int_best, oc_int_best, vt_nint_best, oc_nint_best, vt_mem_best, oc_mem_best]
    if new == 'int' and old == 'nint':
        new_list = [vt_int, oc_int, vt_nint, oc_nint, vt_mem, oc_mem]
        old_list = [vt_nint_best, oc_nint_best, vt_int_best, oc_int_best, vt_mem_best, oc_mem_best]
    if new == 'nint' and old == 'nint':
        new_list = [vt_nint, oc_nint, vt_int, oc_int, vt_mem, oc_mem]
        old_list = [vt_nint_best, oc_nint_best, vt_int_best, oc_int_best, vt_mem_best, oc_mem_best]
    check = compare_cost_list(new_list, old_list)

    return check

def compare_cost_list(new_list, old_list):

    new_highest_cost = new_list[0]+new_list[1]
    old_highest_cost = old_list[0]+old_list[1]

    if new_highest_cost == old_highest_cost:
        if new_list[0] < old_list[0]:
            return True
        elif new_list[0] > old_list[0]:
            return False
        elif new_list[0] == old_list[0]:
            if new_list[1] < old_list[1]:
                return True
            elif new_list[1] > old_list[1]:
                return False
            elif new_list[1] == new_list[1]:
                if len(new_list) == 2:
                    return False
                else:
                    return compare_cost_list(new_list[2:len(new_list)], old_list[2:len(old_list)])
    elif new_highest_cost > old_highest_cost:
        return False
    elif new_highest_cost < old_highest_cost:
        return True
    
def check_fx(idx_fx):
    
    check = True

    vt = len(idx_fx[0])
    oc = len(idx_fx[1])
#    print('vtttt, occcc fx', vt, oc)

    if vt < 2:
        if oc <=5:
            check = True
        else:
            False
    elif vt == 2:
        if oc <=3:
            check = True
        else:
            check = False
    elif vt == 3:
        if oc == 0:
            check == True
        else:
            check = False
    else:
        check = False
        

    # if vt + oc >= 5:
    #     print('czf1')
    #     check = False

    # if vt >= 3:
    #     print('czf2')
    #     check = False
        
    # if fxx == True:
    #     if vt == 2 and oc == 2:
    #         print('czf3')
    #         check = False
#    print('check', check)
    return check
    

def compare_cost(v1, o1, v_org, o_org):

    check = True

    # if v1 < v_org and o1 < o_org:
    #     print('true1')
    #     check = True

    # if v1 == v_org and o1 < o_org:
    #     print('true2')
    #     check = True

    # if v1 < v_org and o1 == o_org:
    #     print('true3')
    #     check = True
    
    if v1 > v_org:
#        print('false')
        check = False
    
    # if v1 < v_org and o1 > o_org:
    #     print('true4')
    #     check = True

    if v1 == 0 and o1 ==0:
#        print('false2')
        check = False

    return check

def divide_interm(Wm, idx_sum, idx_fx, coef_pairs):
    
    interm = ugg()
    noninterm = ugg()

    interm.summation = idx_sum[0] + idx_sum[1]
    noninterm.summation = []
    noninterm.num_factor = Wm.num_factor

    for i in Wm.summation:
        if i not in interm.summation:
            noninterm.summation.append(i)
    
    for j in range(0, len(Wm.coefficient)):
        if j in coef_pairs:
            interm.coefficient.append(Wm.coefficient[j])
            interm.coefficient_idx.append(Wm.coefficient_idx[j])
        else:
            noninterm.coefficient.append(Wm.coefficient[j])
            noninterm.coefficient_idx.append(Wm.coefficient_idx[j])

    noninterm.coefficient.append('Q0')
    noninterm.coefficient_idx.append(idx_fx[0]+idx_fx[1])

    return interm, noninterm

def comp_cost_sum_fx(Wm):
    
    idx = []
    
    for k in range(0, len(Wm.coefficient_idx)):
        for j in Wm.coefficient_idx[k]:
            if j not in idx:
                if j not in Wm.summation:
                    idx.append(j)
    # print('idx', idx)
    vt = 0
    oc = 0
    for k in idx:
        if k in virtual:
            vt += 1
        elif k in occupied:
            oc += 1
    return vt, oc
        
    

def comp_cost_sum(Wm):

    vt = 0
    oc = 0
    for k in Wm.summation:
        if k in virtual:
            vt += 1
        elif k in occupied:
            oc += 1
    return vt, oc


# def generate_intermediates_doubles(Wm_int, name, idx_start):

#     print('laciaton')
#     Wm_intermediates_list = []

#     k = 0
#     for i in range(0, len(Wm_int)):

#         Wm_int[i].clear_fixed()


#         Wm_intermediates_dict = {}

#         if len(Wm_int[i].coefficient) < 5:
#             # Znajdz pojedynczy intermediate
#             a = 2
#             Wm_intermediates_dict = {}#['original'] = Wm_int[i]
#             # Wm_intermediates_dict['noninterm'] = ''
#         elif len(Wm_int[i].coefficient) >=5:     

#             vt = 0
#             oc = 0
#             for ks in Wm_int[i].summation:
#                 if ks in virtual:
#                     vt += 1
#                 elif ks in occupied:
#                     oc += 1

#             # print(vt, ',', oc)
            
#             if vt >= 4 and oc >= 4:
#                 # Znajdz podwojny intermediate
#                 # Zapisz indeksy wszystkiego co nie jest x
#                 coef_idx = []
#                 for j in range(0, len(Wm_int[i].coefficient)):
#                     if Wm_int[i].coefficient[j] != 'x':
#                         coef_idx.append(j)
#                 coef_pairs = coef_combinations(coef_idx)


#                 ABCD_old = 100000

#                 print('----------------------------------------------')
#                 # print(Wm_int[i])
#                 s1 = """v^{vt}o^{oc} """.format(vt=vt, oc=oc)
#                 print('sss1', s1)
#                 for z in range(0, len(coef_pairs)):
#                     c_idx_11 = Wm_int[i].coefficient_idx[coef_pairs[z][0][0]]
#                     c_idx_12 = Wm_int[i].coefficient_idx[coef_pairs[z][0][1]]
#                     c_idx_21 = Wm_int[i].coefficient_idx[coef_pairs[z][1][0]]
#                     c_idx_22 = Wm_int[i].coefficient_idx[coef_pairs[z][1][1]]

#                     idx_sum_1, idx_fx_1 = find_sum_fx(c_idx_11, c_idx_12, coef_pairs[z][0], Wm_int[i])
#                     idx_sum_2, idx_fx_2 = find_sum_fx(c_idx_21, c_idx_22, coef_pairs[z][1], Wm_int[i])

#                     AAA = 100 - (len(idx_sum_1[0] + idx_sum_2[0]))* 2
#                     BBB = (100 - len(idx_sum_1[1] + idx_sum_2[1]))
#                     CCC = len(idx_fx_1[0] + idx_fx_2[0]) 
#                     DDD = len(idx_fx_1[1] + idx_fx_2[1]) * 0.6

#                     if i == 716:
#                         print('sniezner', z)
#                         print(idx_sum_1, idx_fx_1)
#                         print(idx_sum_2, idx_fx_2)

#                     ABCD_new = AAA + BBB + CCC + DDD

#                     comp_int = True


#                     if comp_int == True:
#                         if ABCD_new < ABCD_old:
#                             ABCD_old = deepcopy(ABCD_new)
#                             best_intermediates = coef_pairs[z]
#                             zz = deepcopy(z)
#                             sum_1 = deepcopy(idx_sum_1)
#                             sum_2 = deepcopy(idx_sum_2)
#                             fx_1 = deepcopy(idx_fx_1)
#                             fx_2 = deepcopy(idx_fx_2)
#                             compint = True


#                     c_11 = Wm_int[i].coefficient[coef_pairs[z][0][0]]
#                     c_12 = Wm_int[i].coefficient[coef_pairs[z][0][1]]
#                     c_21 = Wm_int[i].coefficient[coef_pairs[z][1][0]]
#                     c_22 = Wm_int[i].coefficient[coef_pairs[z][1][1]]

#                     # print('kandydat', z)
#                     s2 = """v^{vt1}o^{oc1} """.format(vt1 = len(idx_sum_1[0]), oc1 = len(idx_sum_1[1]))
#                     s3 = """v^{vt1}o^{oc1} """.format(vt1 = len(idx_sum_2[0]), oc1 = len(idx_sum_2[1]))
#                     vt2 = len(list(set(idx_fx_1[0] + idx_fx_2[0])))
#                     oc2 = len(list(set(idx_fx_1[1] + idx_fx_2[1])))
#                     s4 = """v^{vt2}o^{oc2} """.format(vt2 = vt2, oc2 = oc2)

#                     # print(c_idx_11, c_idx_12, c_idx_21, c_idx_22)
#                     # print(c_11, c_12, c_21, c_22)
#                     # print(s1, '--->', s2, s3, s4)
#                     # print(idx_sum_1, '|', idx_sum_2, '|',idx_fx_1, '|',idx_fx_2)
#                     # print(AAA, BBB, CCC, DDD, ABCD_new)
#                     # print('')
#                 interm1, interm2, noninterm = divide_double_interm(Wm_int[i], best_intermediates, sum_1, sum_2, fx_1, fx_2)
#                 # if i == 716:
#                 # if len(fx_1[0] + fx_1[1]) >4 or len(fx_2[0] + fx_2[1]) >4:
#                 #     print('laciaton')
#                 #     print(i, '111', interm1, interm2, noninterm)
#                 #     print(fx_1, fx_2)
#                 #     print(len(fx_1[0] + fx_1[1]), len(fx_2[0] + fx_2[1]))


#                 # print('wybranon kandydata', zz)
#                 # print(Wm_int[i])
#                 # print('', interm1, '', interm2, '', noninterm)
#                 interm11 = deepcopy(interm1)
#                 interm11.standarize()
#                 interm22 = deepcopy(interm2)
#                 interm22.standarize()
#                 # print('', interm11, '', interm22, '', noninterm)
#                 # print('')

#                 Wm_intermediates_dict['interm1'] = interm11
#                 Wm_intermediates_dict['interm2'] = interm22
#                 Wm_intermediates_dict['noninterm'] = noninterm
#                 Wm_intermediates_dict['original'] = Wm_int[i]
#             else:
#                 Wm_intermediates_dict = {}
#                 k += 1

#         Wm_intermediates_list.append(Wm_intermediates_dict)

#     print('itermediatesssss')
#     for i in range(0, len(Wm_intermediates_list)):        
#         if Wm_intermediates_list[i] != {}:
#             if Wm_intermediates_list[i]['interm1'] not in intermediates:
#                 intermediates.append(Wm_intermediates_list[i]['interm1'])

#                 v = 0
#                 o = 0
#                 fx = find_fx(Wm_intermediates_list[i]['interm1'])
#                 # if len(fx) == 4:
#                 #     print('1111', i, fx)
#                 for j in fx:
#                     if j in virtual:
#                         v += 1
#                     elif j in occupied:
#                         o += 1
#                 print(v, ',', o)


#             if Wm_intermediates_list[i]['interm2'] not in intermediates:
#                 intermediates.append(Wm_intermediates_list[i]['interm2'])
#                 v = 0
#                 o = 0
#                 fx = find_fx(Wm_intermediates_list[i]['interm2'])
#                 for j in fx:
#                     if j in virtual:
#                         v += 1
#                     elif j in occupied:
#                         o += 1
#                 print(v, ',', o)
#     print('ln', len(intermediates))

    
            
# def find_fx(Wm):

#     fx = []
#     for i in range(0, len(Wm.coefficient_idx)):
#         for l in Wm.coefficient_idx[i]:
#             if l not in Wm.summation and l not in fx:
#                 fx.append(l)
#     return fx


def divide_triple_interm(Wm_int, best_int, sum1, fx):

    interm = ugg()

    interm.summation = sum1[0] + sum1[1]

    for j in range(0, len(best_int)):
        for k in range(0, len(Wm_int.coefficient)):
            if k == best_int[j]:
                interm.coefficient.append(Wm_int.coefficient[k])
                interm.coefficient_idx.append(Wm_int.coefficient_idx[k])


    noninterm = ugg()
    noninterm.summation = list(set(fx[0] + fx[1]))
    noninterm.num_factor = Wm_int.num_factor

    for i in range(0, len(Wm_int.coefficient)):
        if i not in best_int:
            noninterm.coefficient.append(Wm_int.coefficient[i])
            noninterm.coefficient_idx.append(Wm_int.coefficient_idx[i])
            for j in Wm_int.coefficient_idx[i]:
                if j not in noninterm.summation:
                    noninterm.summation.append(j)

    noninterm.coefficient.append('Q1')
    noninterm.coefficient_idx.append(fx[0] + fx[1])


    return interm, noninterm

def divide_double_interm(Wm_int, best_int, sum_1, sum_2, fx_1, fx_2):

    
    interm1 = ugg()
    interm2 = ugg()
    
    interm1.summation = sum_1[0] + sum_1[1]
    interm2.summation = sum_2[0] + sum_2[1]
    
    for j in range(0, len(best_int[0])):
        for k in range(0, len(Wm_int.coefficient)):
            if k == best_int[0][j]:
                interm1.coefficient.append(Wm_int.coefficient[k])
                interm1.coefficient_idx.append(Wm_int.coefficient_idx[k])

    for j in range(0, len(best_int[1])):
        for k in range(0, len(Wm_int.coefficient)):
            if k == best_int[1][j]:
                interm2.coefficient.append(Wm_int.coefficient[k])
                interm2.coefficient_idx.append(Wm_int.coefficient_idx[k])

    noninterm = ugg()
    noninterm.summation = list(set(fx_1[0] + fx_1[1] + fx_2[0] + fx_2[1]))
    noninterm.num_factor = Wm_int.num_factor
    
    for i in range(0, len(Wm_int.coefficient)):
        if i not in best_int[0] and i not in best_int[1]:
            noninterm.coefficient.append(Wm_int.coefficient[i])
            noninterm.coefficient_idx.append(Wm_int.coefficient_idx[i])
            for j in Wm_int.coefficient_idx[i]:
                if j not in noninterm.summation:
                    noninterm.summation.append(j)

    noninterm.coefficient.append('Q1')
    noninterm.coefficient_idx.append(fx_1[0] + fx_1[1])

    noninterm.coefficient.append('Q2')
    noninterm.coefficient_idx.append(fx_2[0] + fx_2[1])
        
    return interm1, interm2, noninterm

def find_sum_fx(idx_list, coef_pair, Wm_int):

    sum_v_idx = []
    sum_o_idx = []
    fx_v_idx = []
    fx_o_idx = []

    # sprawdz czy indeks i wystepuje tylko wsrod 
    # tych indeksow z ktorych robimy intermediate
    all_other_idx = []
    for i in range(0, len(Wm_int.coefficient)):
        if len(coef_pair) == 2:            
            if i  != coef_pair[0] and i != coef_pair[1]:
                for j in Wm_int.coefficient_idx[i]:
                    all_other_idx.append(j)
        elif len(coef_pair) == 3:
             if i  != coef_pair[0] and i != coef_pair[1] and i != coef_pair[2]:
                for j in Wm_int.coefficient_idx[i]:
                        all_other_idx.append(j)


    
                
    idx = []

    for y in idx_list:
        for x in y:
            if x not in idx:
                idx.append(x)

    for i in idx:
        if i in all_other_idx:
            if i in virtual:
                fx_v_idx.append(i)
            elif i in occupied:
                fx_o_idx.append(i)
        else:
            if i in virtual:
                sum_v_idx.append(i)
            elif i in occupied:
                sum_o_idx.append(i)

    idx_sum = [sum_v_idx, sum_o_idx]
    idx_fx =  [fx_v_idx, fx_o_idx]

    return idx_sum, idx_fx
        
def coef_combinations(coef_idx):
    
    a = list(combinations(coef_idx, 2))
    
    coef_pairs = []

    for i in range(0, len(a)):
        cp = a[i]
        for j in range(i+1, len(a)):
            if a[j][0] not in a[i] and a[j][1] not in a[i]:
                cp2 = a[j]
                coef_pairs.append([cp, cp2])

    return coef_pairs

def coef_combinations_triples(coef_idx):

    a = list(combinations(coef_idx, 3))

    coef_triples = []

    for i in range(0, len(a)):
        coef_triples.append([a[i][0], a[i][1], a[i][2]])

    return coef_triples


def generate_intermediates(Wm_int, name, idx_start, offset):

    # indices_fixed_in_outer_intermediate  with r_{aibj}

    Wm_intermediates_list = []
    print('LEN WM_INT na poczatku')
    print(len(Wm_int))

    for i in range(0, len(Wm_int)):

        idx_sum_old = []
        idx_fx_old = ['a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a']

        # if 'b' in Wm_int[i].summation:
        #     idx_outer_fx = ['a', 'i']
        # else:
        #     idx_outer_fx = ['a', 'i', 'b', 'j']

        Wm_intermediates_dict = {}
        
        # idx_intermediate = set()
        print('----------------------', i, '---------------------------')
        print(Wm_int[i])
        for j in range(0, len(Wm_int[i].coefficient)-offset):
            for k in range(j+1, len(Wm_int[i].coefficient)-offset):
                if Wm_int[i].coefficient[j] !=  OBSERVABLE_X and Wm_int[i].coefficient[j] !=  OBSERVABLE_X_ASYM \
                        and Wm_int[i].coefficient[k] != OBSERVABLE_X and Wm_int[i].coefficient[k] != OBSERVABLE_X_ASYM:
                    idx_intermediate = Wm_int[i].coefficient_idx[j] + Wm_int[i].coefficient_idx[k]
                    idx_nonintermediate = set()
                    for l in range(0, len(Wm_int[i].coefficient)):
                        if l != j and l != k:
                            idx_nonintermediate = idx_nonintermediate |  set(Wm_int[i].coefficient_idx[l])

                    idx_fx  = []
                    idx_sum = []

                    for l in idx_intermediate:
                        if l in idx_nonintermediate:# or l in idx_outer_fx:
                            if l not in idx_fx:
                                idx_fx.append(l)
                        else:
                            if l not in idx_sum:
                                idx_sum.append(l)

                    #
                    # Choose the best intermediate (smallesd idx_fx, highest idx_sum)
                    #

                    choose_interm = False
                    print('wybieram pare', j, k)
                    if len(idx_fx) < len(idx_fx_old):
                        print('tutaj 1', idx_fx, len(idx_fx), idx_fx_old, len(idx_fx_old))
                        if len(idx_sum) > 0:
                            print('tutaj 2', idx_sum)

                            idx_sum_old = idx_sum
                            idx_fx_old = idx_fx
                            choose_interm = True

                    elif len(idx_fx) == len(idx_fx_old):
                        print('tutaj 3', idx_fx, len(idx_fx), idx_fx_old, len(idx_fx_old))

                        old_virt = 0
                        for idx_old in idx_sum_old:
                            if idx_old in virtual:
                                old_virt += 1

                        new_virt = 0
                        for idx_new in idx_sum:
                            if idx_new in virtual:
                                new_virt += 1

                        # if len(idx_sum) > len(idx_sum_old):
                        if new_virt > old_virt:
                            print('tutaj 4', idx_sum, len(idx_sum), idx_sum_old, len(idx_sum_old))
                            idx_sum_old = idx_sum
                            idx_fx_old = idx_fx
                            choose_interm = True
                        elif new_virt == old_virt:
                            if len(idx_sum) > len(idx_sum_old):
                                idx_sum_old = idx_sum
                                idx_fx_old = idx_fx
                                choose_interm = True

                    if choose_interm == True:
                        Wm_intermediates_dict['coef_idx'] = [j, k]
                        Wm_intermediates_dict['idx_sum'] = idx_sum
                        # Wm_intermediates_dict['idx_fx'] = idx_fx
                        Wm_intermediates_dict['original'] = Wm_int[i]

        if Wm_intermediates_dict != {}:
            
            # divide original ugg in nonintermediate (noninter) part and intermediate part (interm)
            
            # idxa and idxb are indices of coefficients in given ugg that form an intermediate
            idxa = Wm_intermediates_dict['coef_idx'][0]
            idxb = Wm_intermediates_dict['coef_idx'][1]

            
            interm = ugg()
            interm.summation = deepcopy(Wm_intermediates_dict['idx_sum'])
            interm.coefficient.append(deepcopy(Wm_int[i].coefficient[idxa]))
            interm.coefficient.append(deepcopy(Wm_int[i].coefficient[idxb]))
            interm.coefficient_idx.append(deepcopy(Wm_int[i].coefficient_idx[idxa]))
            interm.coefficient_idx.append(deepcopy(Wm_int[i].coefficient_idx[idxb]))

            interm.move_summation_to_the_left()
            int_desc = generate_intermediate_descriptor(idxa, idxb, interm)

            idx_fx_sort = generate_intermediate_fixed(interm)

            Wm_intermediates_dict['idx_fx'] = idx_fx_sort
            Wm_intermediates_dict['int_desc'] = int_desc
            Wm_intermediates_dict['interm'] = interm


            noninterm = ugg()
            for m in Wm_int[i].summation:
                if m not in Wm_intermediates_dict['idx_sum']:
                    noninterm.summation.append(m)
            for m in range(0, len(Wm_int[i].coefficient)):
                if m not in Wm_intermediates_dict['coef_idx']:
                    noninterm.coefficient.append(Wm_int[i].coefficient[m])
                    noninterm.coefficient_idx.append(Wm_int[i].coefficient_idx[m])
            noninterm.coefficient.append('Q')
            noninterm.coefficient_idx.append(Wm_intermediates_dict['idx_fx'])
            noninterm.num_factor = Wm_int[i].num_factor
            

            Wm_intermediates_dict['noninterm'] = noninterm
            
            # print(i) , 'original=', Wm_intermediates_dict['original'], 'interm=', Wm_intermediates_dict['interm'], 'noninterm=', Wm_intermediates_dict['noninterm'])
            print(Wm_intermediates_dict['original'], 'original')
            print(Wm_intermediates_dict['interm'], 'interm')
            print(Wm_intermediates_dict['noninterm'], 'noninterm')
            print('')
            
        Wm_intermediates_list.append(Wm_intermediates_dict)
    
    #
    # Some intermediates may be equivalent. In order to compute them only once,  
    # generate list only of different intermediates, and their unique descriptors.
    #
#    sys.exit(0)
    desc = []
    intermediates_dict = []
    k = idx_start
    for i in range(0, len(Wm_intermediates_list)):
        if Wm_intermediates_list[i] != {}:
            if Wm_intermediates_list[i]['int_desc'] not in desc:
                k += 1                
                interm_name = "{name}_interm_{k}".format(name=name, k=k)

                interm = Wm_intermediates_list[i]['interm']
                int_desc = Wm_intermediates_list[i]['int_desc']
                idx_fx = Wm_intermediates_list[i]['idx_fx']

                desc.append(int_desc)
                minidict = {}
                minidict['interm'] = interm
                minidict['int_desc'] = int_desc
                minidict['int_name'] = interm_name
                minidict['idx_fx'] = idx_fx
                intermediates_dict.append(minidict)

   
    #
    # Intermediates_dict contains only non-equivalent intermediates. 
    # Intermediates_dict=[{'interm', 'int_desc', 'int_name', 'idx_fx'}, ...]
    #

  

    #
    # First add to the Wm_intermediates all those arithemtic strings that
    # does not have intermediates, ergo their Wm_intermediates_list[i] == {}
    #


    sum_list = []
    Wm_intermediates = arithmetic_string()

    for j in range(0, len(Wm_intermediates_list)):
        if Wm_intermediates_list[j] == {}:
            Wm_intermediates.append(Wm_int[j])
            sum_list.append(Wm_int[j].summation)
        elif Wm_intermediates_list[j] != {}:
            for i in range(0, len(intermediates_dict)):
                if Wm_intermediates_list[j]['int_desc'] == intermediates_dict[i]['int_desc']:
                    for l in range(0, len(Wm_intermediates_list[j]['noninterm'].coefficient)):
                        if Wm_intermediates_list[j]['noninterm'].coefficient[l] == 'Q':
                            Wm_intermediates_list[j]['noninterm'].coefficient[l] = intermediates_dict[i]['int_name']
                            if(intermediates_dict[i]['int_name'] == 'wm_interm_13'):
                                print(Wm_intermediates_list[j]['original'], Wm_intermediates_list[j]['noninterm'], intermediates_dict[i]['int_name'])
                    Wm_intermediates.append(Wm_intermediates_list[j]['noninterm'])
                    sum_list.append(Wm_intermediates_list[j]['original'].summation)

    return intermediates_dict, Wm_intermediates, sum_list, k


def generate_intermediates_gaxi(gaxi_int, name, method, idxstart, block):

    gaxi_intermediates_list = []
    print('LEN WM_INT na poczatku')
    print(len(gaxi_int))

    if block == 'oo':
        idx_banned = ['i', 'j']
    elif block == 'ov' or block == 'vo':
        idx_banned = ['i', 'a']
    elif block == 'vv':
        idx_banned = ['a', 'b']


    for i in range(0, len(gaxi_int)):
        print(i, gaxi_int[i])
        idx_sum_old = []
        idx_fx_old = ['a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a']

        gaxi_intermediates_dict = {}
        if len(gaxi_int[i].coefficient) > 2:
            for j in range(0, len(gaxi_int[i].coefficient)-1):
                for k in range(j+1, len(gaxi_int[i].coefficient)):
                    idx_intermediate = gaxi_int[i].coefficient_idx[j] + gaxi_int[i].coefficient_idx[k]
                    idx_nonintermediate = set()
                    for l in range(0, len(gaxi_int[i].coefficient)):
                        if l != j and l != k:
                            idx_nonintermediate = idx_nonintermediate |  set(gaxi_int[i].coefficient_idx[l])

                    idx_fx  = []
                    idx_sum = []

                    for l in idx_intermediate:
                        if l in idx_nonintermediate:
                            if l not in idx_fx:
                                idx_fx.append(l)
                        else:
                            if l not in idx_sum:
                                if l not in general:
                                    if l not in idx_banned:
                                        idx_sum.append(l)

                    #
                    # Choose the best intermediate (smallesd idx_fx, highest idx_sum)
                    #

                    choose_interm = False
                    if len(idx_fx) < len(idx_fx_old):
                        if len(idx_sum) > 0:
                            idx_sum_old = idx_sum
                            idx_fx_old = idx_fx
                            choose_interm = True

                    elif len(idx_fx) == len(idx_fx_old):
                        if len(idx_sum) > len(idx_sum_old):
                            idx_sum_old = idx_sum
                            idx_fx_old = idx_fx
                            choose_interm = True

                    if choose_interm == True:
                        
                        sum_nonint = list(set(gaxi_int[i].summation) - set(idx_sum))
                        if block == 'oo':
                            oro = 2
                            orv = 0
                            nono = 2
                            nonv = 0
                            # into = 2
                            # intv = 0
                        elif block == 'ov' or block == 'vo':
                            oro = 1
                            orv = 1
                            nono = 1
                            nonv = 1
                            # into = 1
                            # intv = 1
                        elif block == 'vv':
                            oro = 0
                            orv = 2
                            nono = 0
                            nonv = 2
                            # into = 0
                            # intv = 2
                        into = 0
                        intv = 0
                        for q in gaxi_int[i].summation:
                            if q in virtual:
                                orv += 1
                            elif q in occupied:
                                oro += 1
                            
                        for q in sum_nonint:
                            if q in virtual:
                                nonv += 1
                            elif q in occupied:
                                nono += 1

                        for q in idx_fx:
                            if q in virtual:
                                intv += 1
                            elif q in occupied:
                                into += 1
                                
                        for q in idx_sum:
                            if q in virtual:
                                intv += 1
                            elif q in occupied:
                                into += 1
    
                        compare_org_nonint = [orv - nonv, oro - nono]
                        compare_org_int = [orv - intv, oro - into]
                        print(compare_org_nonint, compare_org_int, len(idx_fx))
                        nointerm = False
                        if compare_org_nonint[0] <= 0 and compare_org_nonint[1] <= 0:
                            gaxi_intermediates_dict = {}
                            nointerm = True
                        if compare_org_int[0] <= 0 and compare_org_int[1] <= 0:
                            gaxi_intermediates_dict = {}
                            nointerm = True
                        if nointerm == False:
                            gaxi_intermediates_dict['coef_idx'] = [j, k]
                            gaxi_intermediates_dict['idx_sum'] = idx_sum                            
                            gaxi_intermediates_dict['original'] = gaxi_int[i]


            if gaxi_intermediates_dict != {}:

                # divide original ugg in nonintermediate (noninter) part and intermediate part (interm)

                # idxa and idxb are indices of coefficients in given ugg that form an intermediate
                idxa = gaxi_intermediates_dict['coef_idx'][0]
                idxb = gaxi_intermediates_dict['coef_idx'][1]


                interm = ugg()
                interm.summation = deepcopy(gaxi_intermediates_dict['idx_sum'])
                interm.coefficient.append(deepcopy(gaxi_int[i].coefficient[idxa]))
                interm.coefficient.append(deepcopy(gaxi_int[i].coefficient[idxb]))
                interm.coefficient_idx.append(deepcopy(gaxi_int[i].coefficient_idx[idxa]))
                interm.coefficient_idx.append(deepcopy(gaxi_int[i].coefficient_idx[idxb]))

                interm.move_summation_to_the_left()
                int_desc = generate_intermediate_descriptor(idxa, idxb, interm)

                idx_fx_sort = generate_intermediate_fixed(interm)

                gaxi_intermediates_dict['idx_fx'] = idx_fx_sort
                gaxi_intermediates_dict['int_desc'] = int_desc
                gaxi_intermediates_dict['interm'] = interm

                noninterm = ugg()
                for m in gaxi_int[i].summation:
                    if m not in gaxi_intermediates_dict['idx_sum']:
                        noninterm.summation.append(m)
                for m in range(0, len(gaxi_int[i].coefficient)):
                    if m not in gaxi_intermediates_dict['coef_idx']:
                        noninterm.coefficient.append(gaxi_int[i].coefficient[m])
                        noninterm.coefficient_idx.append(gaxi_int[i].coefficient_idx[m])
                noninterm.coefficient.append('Q')
                noninterm.coefficient_idx.append(gaxi_intermediates_dict['idx_fx'])
                noninterm.num_factor = gaxi_int[i].num_factor


                gaxi_intermediates_dict['noninterm'] = noninterm

            gaxi_intermediates_list.append(gaxi_intermediates_dict)
        else:
            print(i, 'tu')
            gaxi_intermediates_list.append({})

    #
    # Some intermediates may be equivalent. In order to compute them only once,  
    # generate list only of different intermediates, and their unique descriptors.
    #

    desc = []
    intermediates_dict = []
    k = idxstart
    for i in range(0, len(gaxi_intermediates_list)):
        if gaxi_intermediates_list[i] != {}:
            if gaxi_intermediates_list[i]['int_desc'] not in desc:
                k += 1
                interm_name = "{name}_interm_{k}_{method}".format(name=name, k=k, method=method)

                interm = gaxi_intermediates_list[i]['interm']
                int_desc = gaxi_intermediates_list[i]['int_desc']
                idx_fx = gaxi_intermediates_list[i]['idx_fx']

                desc.append(int_desc)
                minidict = {}
                minidict['interm'] = interm
                minidict['int_desc'] = int_desc
                minidict['int_name'] = interm_name
                minidict['idx_fx'] = idx_fx
                minidict['original'] = gaxi_intermediates_list[i]['original']
                minidict['noninterm'] = gaxi_intermediates_list[i]['noninterm']
                minidict['idx_sum'] = gaxi_intermediates_list[i]['idx_sum']
                intermediates_dict.append(minidict)

   
    #
    # Intermediates_dict contains only non-equivalent intermediates. 
    # Intermediates_dict=[{'interm', 'int_desc', 'int_name', 'idx_fx'}, ...]
    #

    #
    # First add to the gaxi_intermediates all those arithemtic strings that
    # does not have intermediates, ergo their gaxi_intermediates_list[i] == {}
    #


    sum_list = []
    gaxi_intermediates = arithmetic_string()

    for j in range(0, len(gaxi_intermediates_list)):
        if gaxi_intermediates_list[j] == {}:
            gaxi_intermediates.append(gaxi_int[j])
            print('appenduje')
            sum_list.append(gaxi_int[j].summation)
        elif gaxi_intermediates_list[j] != {}:
            for i in range(0, len(intermediates_dict)):
                if gaxi_intermediates_list[j]['int_desc'] == intermediates_dict[i]['int_desc']:
                    for l in range(0, len(gaxi_intermediates_list[j]['noninterm'].coefficient)):
                        if gaxi_intermediates_list[j]['noninterm'].coefficient[l] == 'Q':
                            gaxi_intermediates_list[j]['noninterm'].coefficient[l] = intermediates_dict[i]['int_name']
                            if block == 'oo':
                                oro = 2
                                orv = 0
                                nono = 2
                                nonv = 0
                                into = 2
                                intv = 0
                            elif block == 'ov' or block == 'vo':
                                oro = 1
                                orv = 1
                                nono = 1
                                nonv = 1
                                into = 1
                                intv = 1
                            elif block == 'vv':
                                oro = 0
                                orv = 2
                                nono = 0
                                nonv = 2
                                into = 0
                                intv = 2
                                
                            for q in gaxi_intermediates_list[j]['original'].summation:
                                if q in virtual:
                                    orv += 1
                                elif q in occupied:
                                    oro += 1
                            
                            for q in gaxi_intermediates_list[j]['noninterm'].summation:
                                if q in virtual:
                                    nonv += 1
                                elif q in occupied:
                                    nono += 1

                            for q in gaxi_intermediates_list[j]['interm'].summation:
                                if q in virtual:
                                    intv += 1
                                elif q in occupied:
                                    into += 1

                            for q in gaxi_intermediates_list[j]['idx_fx']:
                                if q in virtual:
                                    intv += 1
                                elif q in occupied:
                                    into += 1
                    print(gaxi_intermediates_list[j]['original'])
                    print(gaxi_intermediates_list[j]['noninterm'])
                    print(gaxi_intermediates_list[j]['interm'])
                    print('')
                    gaxi_intermediates.append(gaxi_intermediates_list[j]['noninterm'])
                    sum_list.append(gaxi_intermediates_list[j]['original'].summation)


    return intermediates_dict, gaxi_intermediates, sum_list, k


def divide_to_r1r2_parts(Wm_with_intermediates, sum_list):

    
    Wm_with_intermediates_r1 = arithmetic_string()
    Wm_with_intermediates_r2 = arithmetic_string()

    print('')
    print('')
    print('PLUSZ')
    for i in range(0, len(Wm_with_intermediates)):
        print(Wm_with_intermediates[i])
        id = []
        if Wm_with_intermediates[i].summation == []:
            print('a00')
            for j in range(0, len(Wm_with_intermediates[i].coefficient)):
                for k in Wm_with_intermediates[i].coefficient_idx[j]:
                    if k not in id:
                        id.append(k)
            print('id', id)
            if 'b' not in id:
                Wm_with_intermediates_r1.append(Wm_with_intermediates[i])
            else:
                Wm_with_intermediates_r2.append(Wm_with_intermediates[i])
        elif 'b' in sum_list[i]:
            print('b00')
            Wm_with_intermediates_r1.append(Wm_with_intermediates[i])
        else:
            print('c00')
            Wm_with_intermediates_r2.append(Wm_with_intermediates[i])



    for i in range(0, len(Wm_with_intermediates_r1)):

        Wm_with_intermediates_r1[i].coefficient.append(EOM_CC_SINGLE_Rl)
        Wm_with_intermediates_r1[i].coefficient_idx.append(['a', 'i'])
        Wm_with_intermediates_r1[i].summation.append('a')
        Wm_with_intermediates_r1[i].summation.append('i')

    for i in range(0, len(Wm_with_intermediates_r2)):
        Wm_with_intermediates_r2[i].coefficient.append(EOM_CC_SINGLE_Rl)
        Wm_with_intermediates_r2[i].coefficient_idx.append(['a', 'i', 'b', 'j'])
        Wm_with_intermediates_r2[i].summation.append('a')
        Wm_with_intermediates_r2[i].summation.append('i')
        Wm_with_intermediates_r2[i].summation.append('b')
        Wm_with_intermediates_r2[i].summation.append('j')

    return Wm_with_intermediates_r1, Wm_with_intermediates_r2

def divide_to_r1r2_parts_2(Wm_with_intermediates_r1, Wm_with_intermediates_r2):

    for i in range(0, len(Wm_with_intermediates_r1)):
        Wm_with_intermediates_r1[i].coefficient.append(EOM_CC_SINGLE_Rl)
        Wm_with_intermediates_r1[i].coefficient_idx.append(['a', 'i'])
        Wm_with_intermediates_r1[i].summation.append('a')
        Wm_with_intermediates_r1[i].summation.append('i')

    for i in range(0, len(Wm_with_intermediates_r2)):
        Wm_with_intermediates_r2[i].coefficient.append(EOM_CC_SINGLE_Rl)
        Wm_with_intermediates_r2[i].coefficient_idx.append(['a', 'i', 'b', 'j'])
        Wm_with_intermediates_r2[i].summation.append('a')
        Wm_with_intermediates_r2[i].summation.append('i')
        Wm_with_intermediates_r2[i].summation.append('b')
        Wm_with_intermediates_r2[i].summation.append('j')

    return Wm_with_intermediates_r1, Wm_with_intermediates_r2

def append_r1r2_parts(Wm_with_intermediates_r1, Wm_with_intermediates_r2, Wm_with_intermediates_r3 = None):

    Wm_with_intermediates = arithmetic_string()
#    Teraz to juz nie dodaje r1 i r2 i r3 bo one sa liczone od poczatku
    for i in range(0, len(Wm_with_intermediates_r1)):
        Wm_with_intermediates.append(Wm_with_intermediates_r1[i])
    for i in range(0, len(Wm_with_intermediates_r2)):
        Wm_with_intermediates.append(Wm_with_intermediates_r2[i])
    if Wm_with_intermediates_r3:
        for i in range(0, len(Wm_with_intermediates_r3)):
            Wm_with_intermediates.append(Wm_with_intermediates_r3[i])


    # for i in range(0, len(Wm_with_intermediates_r1)):
    #     Wm_with_intermediates_r1[i].coefficient.append(EOM_CC_SINGLE_Rl)
    #     Wm_with_intermediates_r1[i].coefficient_idx.append(['a', 'i'])
    #     Wm_with_intermediates_r1[i].summation.append('a')
    #     Wm_with_intermediates_r1[i].summation.append('i')
    #     Wm_with_intermediates.append(Wm_with_intermediates_r1[i])


    # for i in range(0, len(Wm_with_intermediates_r2)):
    #     Wm_with_intermediates_r2[i].coefficient.append(EOM_CC_SINGLE_Rl)
    #     Wm_with_intermediates_r2[i].coefficient_idx.append(['a', 'i', 'b', 'j'])
    #     Wm_with_intermediates_r2[i].summation.append('a')
    #     Wm_with_intermediates_r2[i].summation.append('i')
    #     Wm_with_intermediates_r2[i].summation.append('b')
    #     Wm_with_intermediates_r2[i].summation.append('j')
    #     Wm_with_intermediates_r2[i].num_factor *= 0.5
    #     Wm_with_intermediates.append(Wm_with_intermediates_r2[i])

    # if Wm_with_intermediates_r3:
    #     for i in range(0, len(Wm_with_intermediates_r3)):
    #         Wm_with_intermediates_r3[i].coefficient.append(EOM_CC_SINGLE_Rl)
    #         Wm_with_intermediates_r3[i].coefficient_idx.append(['a', 'i', 'b', 'j', 'c', 'k'])
    #         Wm_with_intermediates_r3[i].summation.append('a')
    #         Wm_with_intermediates_r3[i].summation.append('i')
    #         Wm_with_intermediates_r3[i].summation.append('b')
    #         Wm_with_intermediates_r3[i].summation.append('j')
    #         Wm_with_intermediates_r3[i].summation.append('c')
    #         Wm_with_intermediates_r3[i].summation.append('k')
    #         Wm_with_intermediates_r3[i].num_factor *= 1./6.
    #         Wm_with_intermediates.append(Wm_with_intermediates_r3[i])


    # print('wmwmwm1')
    # for x in Wm_with_intermediates_r1:
    #     print(x)
    # print('wmwmwm2')
    # for x in Wm_with_intermediates_r2:
    #     print(x)
    # print('wmwmwm3')
    # for x in Wm_with_intermediates_r3:
    #     print(x)

    Wm_with_intermediates = simplify(Wm_with_intermediates)
    print('ilosc interm jak cal', len(Wm_with_intermediates))

    return Wm_with_intermediates

def append_r1r2_parts_triplet(Wm_with_intermediates_r1, Wm_with_intermediates_r2_plus, Wm_with_intermediates_r2_minus, \
                              Wm_with_intermediates_r3 = None):

    Wm_with_intermediates = arithmetic_string()

    print('wwwwww1')
    for i in range(0, len(Wm_with_intermediates_r1)):
        Wm_with_intermediates_r1[i].coefficient.append(EOM_CC_SINGLE_Rl)
        Wm_with_intermediates_r1[i].coefficient_idx.append(['a', 'i'])
        Wm_with_intermediates_r1[i].summation.append('a')
        Wm_with_intermediates_r1[i].summation.append('i')
        Wm_with_intermediates.append(Wm_with_intermediates_r1[i])
        print(Wm_with_intermediates_r1[i].summation)

    print('wwwwww plus')
    for i in range(0, len(Wm_with_intermediates_r2_plus)):
        Wm_with_intermediates_r2_plus[i].coefficient.append(EOM_CC_SINGLE_Rl_plus)
        Wm_with_intermediates_r2_plus[i].coefficient_idx.append(['a', 'i', 'b', 'j'])
        Wm_with_intermediates_r2_plus[i].summation.append('a')
        Wm_with_intermediates_r2_plus[i].summation.append('i')
        Wm_with_intermediates_r2_plus[i].summation.append('b')
        Wm_with_intermediates_r2_plus[i].summation.append('j')
        Wm_with_intermediates_r2_plus[i].num_factor *= 0.5
        Wm_with_intermediates.append(Wm_with_intermediates_r2_plus[i])
        print(Wm_with_intermediates_r2_plus[i].summation)

    print('wwwwww minus')
    for i in range(0, len(Wm_with_intermediates_r2_minus)):
        Wm_with_intermediates_r2_minus[i].coefficient.append(EOM_CC_SINGLE_Rl_minus)
        Wm_with_intermediates_r2_minus[i].coefficient_idx.append(['a', 'i', 'b', 'j'])
        Wm_with_intermediates_r2_minus[i].summation.append('a')
        Wm_with_intermediates_r2_minus[i].summation.append('i')
        Wm_with_intermediates_r2_minus[i].summation.append('b')
        Wm_with_intermediates_r2_minus[i].summation.append('j')
        Wm_with_intermediates.append(Wm_with_intermediates_r2_minus[i])
        print(Wm_with_intermediates_r2_minus[i].summation)

    if Wm_with_intermediates_r3:
        for i in range(0, len(Wm_with_intermediates_r3)):
            Wm_with_intermediates_r3[i].coefficient.append(EOM_CC_SINGLE_Rl)
            Wm_with_intermediates_r3[i].coefficient_idx.append(['a', 'i', 'b', 'j', 'c', 'k'])
            Wm_with_intermediates_r3[i].summation.append('a')
            Wm_with_intermediates_r3[i].summation.append('i')
            Wm_with_intermediates_r3[i].summation.append('b')
            Wm_with_intermediates_r3[i].summation.append('j')
            Wm_with_intermediates_r3[i].summation.append('c')
            Wm_with_intermediates_r3[i].summation.append('k')
            Wm_with_intermediates_r3[i].num_factor *= 1./2.
            Wm_with_intermediates.append(Wm_with_intermediates_r3[i])

    return Wm_with_intermediates
    
def select_cc3_terms_only(W_middle):

    W_middle_cc3 = []
    z = 0
    for x in W_middle:
        z += 1
        if 3 in x['Wm1_T_list'] or 3 in x['Wm1_S_list'] or 3 in x['Wm2_T_list'] \
           or 3 in x['Wm2_S_list'] or 3 in x['Wm3_S_list'] or x['Wm3_n'] == 3 or x['Wm2_n']==3:
            print('ZZZZZZZZZZZZZZZ' ,z)
            W_middle_cc3.append(x)
            
    return W_middle_cc3

def select_cc3_terms_overlap(Wm_overlap):

    Wm_overlap_cc3 = []

    for x in Wm_overlap:

        if 3 in x['Wm2_T_list'] or 3 in x['Wm2_S_list'] or 3 in x['Wm3_S_list']:
            Wm_overlap_cc3.append(x)
            
    return Wm_overlap_cc3

def dump_overlap(maxmbpt, theory, outo, multiplicity, cumulative):

    print('overlap')

    Wm_overlap = generate_commutators_overlap(maxmbpt, theory, cumulative)

    if theory == 'overlap_cc3':
        Wm_overlap = select_cc3_terms_overlap(Wm_overlap)

    print('latex')
    latex_Wm_overlap(Wm_overlap, cumulative, maxmbpt)


    if multiplicity == 1:
        Wm2_oc,  Wm3_oc, Wm2_omu_list, Wm3_omu_list = commutators_quadra_overlap(Wm_overlap)

    elif multiplicity == 3:

        Wm2_oc,  Wm3_oc, Wm2_omu_list, Wm3_omu_list = commutators_quadra_overlap_triplet(Wm_overlap)

    Wm_oint = integrate_quadra_Wm_overlap(Wm2_oc,  Wm3_oc, Wm2_omu_list, Wm3_omu_list)

    print('WYNIK CALKI')
    for x in Wm_oint:
        print(x)

    print('dumping Wm_overlap')
    pickle.dump(Wm_oint, outo)

    print('koniec')


def load_overlap(Wm_oint, theory, multiplicity, mbpt):


    for i in Wm_oint:
        print(i)


    Wm_oint = simplify(Wm_oint)
    print(len(Wm_oint))

    print('generate_best')
    Wm_oint_list, nonintermediates = generate_best_intermediates(Wm_oint, 'wmo', 0)


    # print('generate_best2')
    # Wm_oint_list_2, nonintermediates = generate_best_intermediates(nonintermediates, 'wmo', 0, False)

    # Wm_oint_list = []
    # for x in Wm_oint_list_1:
    #     Wm_oint_list.append(x)
    # for x in Wm_oint_list_2:
    #     Wm_oint_list.append(x)

    print('przed usunieciem drogich')
    for x in Wm_oint_list:
        print(x['interm1'], x['interm2'])

    Wm_oint_list_temp = deepcopy(Wm_oint_list)
    Wm_oint_list = []

    for x in Wm_oint_list_temp:
        vt1, oc1 = comp_cost_sum_fx(x['interm1'])
        vt2, oc2 = comp_cost_sum_fx(x['interm2'])
        vt3, oc3 = comp_cost_sum(x['noninterm'])
        print('adadada', vt1, oc1, ',', vt2, oc2, ',', vt3, oc3)
        if vt1 + oc1 < 7 and vt2 + oc2 < 7 and vt3 + oc3 < 7:
            Wm_oint_list.append(x)
            print('dodaje++')

    print('po usunieciu drogichch')
    for x in Wm_oint_list:
        print(x['interm1'], x['interm2'])


    intermediates, Wm_oint_intermediates, k = simplify_intermediates([], Wm_oint_list, nonintermediates, 'wmo', 0, multiplicity, mbpt)

    for x in Wm_oint_intermediates:
        if (len(x.summation) > 6):
            print('koszt', len(x.summation), x)

    if theory == 'overlap_ccsd':
        function_template_wm_intermediates(intermediates, 'ccsd', 'wmo', multiplicity, mbpt)
    elif theory == 'overlap_cc3':
        function_template_wm_intermediates(intermediates, 'cc3', 'wmo', multiplicity, mbpt)

    Wm_oint = arithmetic_string()
    for x in Wm_oint_intermediates:
        Wm_oint.append(x)

    function_template_wm_overlap(Wm_oint, theory, multiplicity, mbpt)


def execute_transition_exc(maxmbpt, maxcluster, theory, pick, multiplicity, cumulative):

    start11 = time.time()
    open_flags = ""
    if pick == "dump":
        open_flags = "wb"        
    elif pick == "load":
        open_flags = "rb"

    if theory == 'ccsd':
        
        s1 = './pickle/excexc/W_middle_ccsd_pt{pt}.pkl'.format(pt = maxmbpt)
        
        out = open(s1, open_flags)

        s2 = './pickle/excexc/Wm_int_1_ccsd_pt{pt}.pkl'.format(pt = maxmbpt)
        s3 = './pickle/excexc/Wm_int_2_pt{pt}.pkl'.format(pt = maxmbpt)
        s4 = './pickle/excexc/Wm_int_2_plus_pt{pt}.pkl'.format(pt = maxmbpt)
        s5 = './pickle/excexc/Wm_int_2_minus_pt{pt}.pkl'.format(pt = maxmbpt)

        out1 = open(s2, open_flags)
        out2 = open(s3, open_flags)
        out2_plus = open(s4, open_flags)
        out2_minus = open(s5, open_flags)
        
        s6 = './pickle/excexc/Wm_int_11_ccsd_pt{pt}.pkl'.format(pt = maxmbpt)
        s7 = './pickle/excexc/Wm_int_12_ccsd_pt{pt}.pkl'.format(pt = maxmbpt)
        s8 = './pickle/excexc/Wm_int_21_ccsd_pt{pt}.pkl'.format(pt = maxmbpt)
        s9 = './pickle/excexc/Wm_int_22_ccsd_pt{pt}.pkl'.format(pt = maxmbpt)

        s6a = './pickle/excexc/Wm_int_1_ccsd_so_pt{pt}.pkl'.format(pt = maxmbpt)
        s7a = './pickle/excexc/Wm_int_2_ccsd_so_pt{pt}.pkl'.format(pt = maxmbpt)

        out1so = open(s6a, open_flags)
        out2so = open(s7a, open_flags)


        out11 = open(s6, open_flags)
        out12 = open(s7, open_flags)
        out21 = open(s8, open_flags)
        out22 = open(s9, open_flags)

        s10 = './pickle/excexc/Um_int_11_ccsd_pt{pt}.pkl'.format(pt = maxmbpt)
        s11 = './pickle/excexc/Um_int_12_ccsd_pt{pt}.pkl'.format(pt = maxmbpt)
        s12 = './pickle/excexc/Um_int_21_ccsd_pt{pt}.pkl'.format(pt = maxmbpt)
        s13 = './pickle/excexc/Um_int_22_ccsd_pt{pt}.pkl'.format(pt = maxmbpt)

            
        outu11 = open(s10, open_flags)
        outu12 = open(s11, open_flags)
        outu21 = open(s12, open_flags)
        outu22 = open(s13, open_flags)
        
        s14 = './pickle/excexc/U_middle_ccsd_pt{pt}.pkl'.format(pt = maxmbpt)
        s15 = './pickle/excexc/Wl_ccsd_pt{pt}.pkl'.format(pt = maxmbpt)
        s16 = './pickle/excexc/Wr_ccsd_pt{pt}.pkl'.format(pt = maxmbpt)


        outu = open(s14, open_flags)
        outl = open(s15, open_flags)
        outr = open(s16, open_flags)

    elif theory == 'cc3':
        sss1 = """./pickle/excexc/Wm_int_1_cc3_pt{maxmbpt}.pkl""".format(maxmbpt = maxmbpt)
        sss2 = """./pickle/excexc/Wm_int_2_cc3_pt{maxmbpt}.pkl""".format(maxmbpt = maxmbpt)
        sss3 = """./pickle/excexc/Wm_int_3_cc3_pt{maxmbpt}.pkl""".format(maxmbpt = maxmbpt)
        cc3_out1 = open(sss1, open_flags)
        cc3_out2 = open(sss2, open_flags)
        cc3_out3 = open(sss3, open_flags)

        sss1so = """./pickle/excexc/Wm_int_1_cc3_so_pt{maxmbpt}.pkl""".format(maxmbpt = maxmbpt)
        sss2so = """./pickle/excexc/Wm_int_2_cc3_so_pt{maxmbpt}.pkl""".format(maxmbpt = maxmbpt)
        sss3so = """./pickle/excexc/Wm_int_3_cc3_so_pt{maxmbpt}.pkl""".format(maxmbpt = maxmbpt)
        cc3_out1so = open(sss1so, open_flags)
        cc3_out2so = open(sss2so, open_flags)
        cc3_out3so = open(sss3so, open_flags)

        # sss1 = """./pickle/excexc/Wm_int_1_cc3.pkl"""
        # sss2 = """./pickle/excexc/Wm_int_2_cc3.pkl"""
        # sss3 = """./pickle/excexc/Wm_int_3_cc3.pkl"""
        # cc3_out1 = open(sss1, open_flags)
        # cc3_out2 = open(sss2, open_flags)
        # cc3_out3 = open(sss3, open_flags)

        sss1trip = """./pickle/excexc/Wm_int_1_cc3_trip_pt{maxmbpt}.pkl""".format(maxmbpt = maxmbpt)
        sss2ptrip = """./pickle/excexc/Wm_int_2_cc3_plus_trip_pt{maxmbpt}.pkl""".format(maxmbpt = maxmbpt)
        sss2mtrip = """./pickle/excexc/Wm_int_2_cc3_minus_trip_pt{maxmbpt}.pkl""".format(maxmbpt = maxmbpt)
        sss3trip = """./pickle/excexc/Wm_int_3_cc3_trip_so_pt{maxmbpt}.pkl""".format(maxmbpt = maxmbpt)
        cc3_out1trip = open(sss1trip, open_flags)
        cc3_out2ptrip = open(sss2ptrip, open_flags)
        cc3_out2mtrip = open(sss2mtrip, open_flags)
        cc3_out3trip = open(sss3trip, open_flags)


        
        cc3_out1_trip = open('./pickle/excexc/Wm_int_1_cc3_trip.pkl', open_flags)
        cc3_out2_trip = open('./pickle/excexc/Wm_int_2_cc3_trip.pkl', open_flags)
        cc3_out2_plus_trip = open('./pickle/excexc/Wm_int_2_cc3_plus_trip.pkl', open_flags)
        cc3_out2_minus_trip = open('./pickle/excexc/Wm_int_2_cc3_minus_trip.pkl', open_flags)
        cc3_out3_trip = open('./pickle/excexc/Wm_int_3_cc3_trip.pkl', open_flags)

    elif theory == 'overlap_ccsd':
        if multiplicity == 1:
            outo = open('./pickle/excexc/Wm_int_ccsd_overlap.pkl', open_flags)
        elif multiplicity == 3:
            outo = open('./pickle/excexc/Wm_int_ccsd_overlap_triplet.pkl', open_flags)
    elif theory == 'overlap_cc3':
        if multiplicity == 1:
            outo = open('./pickle/excexc/Wm_int_cc3_overlap.pkl', open_flags)
        elif multiplicity == 3:
            outo = open('./pickle/excexc/Wm_int_cc3_overlap_triplet.pkl', open_flags)

    if pick == "dump":

        #-----------------------------------------------------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------------------------------------------------
        #
        # W_left, W_middle, W_right, U_middle are lists of little dictionaries. W_left and W_right contins lists of T operators
        # e.g. [T1,T1,T2...]. W_middle contains five lists Wm1_S_list, Wm1_T_list, Wm2_T_list, Wm2_S_list, Wm3_S_list.
        # The restrictions on the maximum excitation of T, s, mu, maximum number of nested commutators, and maximum number of MBPT
        # order are given by the user in the input.
        #
        # <e(-T)Ye(T)|\mu_m>    <e(-S)e(T*) X e(-T*) e(S) P(e(-S)e(T*) \mu_n e(-T*) e(S))| P(e(S*)\mu_l e(-S*))>    <mu_p|e(-T)Ye(T)>
        #-----Wl-----------#   #--------------------------------------Wm-----------------------------------#    #------Wr-------#
        #
        #                      Um = <  P(e(S*)\mu_l e(-S*)) | e(-S)e(T*) X e(-T*) e(S) P(e(-S)e(T*) \mu_n e(-T*) e(S))> 
        #                             --------Wm3---------# #--------Wm1-----------# #----------Wm2---------------#
        #
        #
        #-----------------------------------------------------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------------------------------------------------

        if theory == 'overlap_ccsd' or theory == 'overlap_cc3':
            dump_overlap(maxmbpt, theory, outo, multiplicity, cumulative)


        else:
            W_left, W_middle, W_right, U_middle = generate_commutators_for_trans_exc(maxmbpt, maxcluster, theory, cumulative)

            print('SAME KOMUTATORY')
            for x in W_middle:
                print('tra', x)


            ttt = []

            print('po gen')
            for i in range(0, len(W_middle)):
                print(i, W_middle[i])

            if theory == 'cc3':
                print('select cc3')
                W_middle = select_cc3_terms_only(W_middle)
                print('W-middle')
                for x in W_middle:
                    print(x)
            print(len(W_middle))
#            sys.exit(0)
            W_middle = cut_commutators_Wm(W_middle, maxmbpt)
            print('W-middle - po cut')
            print(len(W_middle))            
#            W_middle = [W_middle[0]]#!!!!!!!!!!!!!!!!
            k = 0
            for x in W_middle:
                print('kkk', k, x)
                k += 1
#            sys.exit(0)

            print(len(W_middle))


            if multiplicity == 1:
                W_middle_out, Wm1_c,  Wm2_c,  Wm3_c, Wm2_mu_list, Wm3_mu_list = commutators_quadra_Wm(W_middle)
            elif multiplicity == 3:
                W_middle_out, Wm1_c,  Wm2_c,  Wm3_c, Wm2_mu_list, Wm3_mu_list = commutators_quadra_Wm_triplet(deepcopy(W_middle))
            elif multiplicity == 13:
                print('jestem tu')
                W_middle_out, Wm1_c,  Wm2_c,  Wm3_c, Wm2_mu_list, Wm3_mu_list = commutators_quadra_Wm_singlet_triplet(deepcopy(W_middle))
                #Wm1_c,  Wm2_c,  Wm3_c, Wm2_mu_list, Wm3_mu_list = commutators_quadra_Wm_singlet_triplet_full(deepcopy(W_middle))
            elif multiplicity == 31:
                W_middle_out, Wm1_c,  Wm2_c,  Wm3_c, Wm2_mu_list, Wm3_mu_list = commutators_quadra_Wm_triplet_singlet(deepcopy(W_middle))


                
            print('W_middle')
            print(len(W_middle_out), len(Wm1_c), len(Wm2_c), len(Wm3_c), len(Wm2_mu_list))
            # Wm1cp  = deepcopy(Wm1_c)
            # Wm2cp = deepcopy(Wm2_c)
            # Wm3cp = deepcopy(Wm3_c)
            # Wm2listcp = deepcopy(Wm2_mu_list)
            # Wm3listcp = deepcopy(Wm3_mu_list)
            # nn = 1

            # Wm1_c = [Wm1cp[nn]]
            # Wm2_c = [Wm2cp[nn]]
            # Wm3_c = [Wm3cp[nn]]
            # Wm2_mu_list = [Wm2listcp[nn]]
            # Wm3_mu_list = [Wm3listcp[nn]]
            print('latex', theory, maxmbpt)
            print('przed calkowaniem')
            latex_W_middle_big(W_middle_out, Wm1_c,  Wm2_c,  Wm3_c, maxmbpt, cumulative)
            print('end')

            #pokeplusz
            
            
            print('Wm_integrate')
            if multiplicity == 1:

                if theory == 'ccsd':
                     Wm_int_1, Wm_int_2 = \
                        integrate_quadra_Wm(Wm1_c,  Wm2_c,  Wm3_c, Wm2_mu_list, Wm3_mu_list, theory)

                elif theory == 'cc3':
                    # print('llllllllll')
                    # print('')
                    # for x in range(0, len(Wm2_c)):
                    #     for y in Wm2_c[x]:
                    #         print('wm2', Wm2_mu_list[x], Wm3_mu_list[x], y)
                    # print('')
                    # for x in range(0, len(Wm2_c)):    

                    #     for y in Wm1_c[x]:
                    #         print('wm1', Wm2_mu_list[x], Wm3_mu_list[x], y)

                    # print('pikapika')

                    Wm_int_1, Wm_int_2, Wm_int_3 = \
                        integrate_quadra_Wm(Wm1_c,  Wm2_c,  Wm3_c, Wm2_mu_list, Wm3_mu_list, theory)

                    print('wmwmwm1')
                    for x in Wm_int_1:
                        print(x)
                    print('wmwmwm2')
                    for x in Wm_int_2:
                        print(x)
                    print('wmwmwm3')
                    for x in Wm_int_3:
                        print(x)

            elif multiplicity == 3:
                if theory == 'ccsd':
                    Wm_int_1, Wm_int_2_plus, Wm_int_2_minus = \
                       integrate_quadra_Wm_triplet(Wm1_c,  Wm2_c,  Wm3_c, Wm2_mu_list, \
                                                   Wm3_mu_list, theory, multiplicity=multiplicity)
                    print('dumping Wm, squirtle')

                    'W1'
                    for x in Wm_int_1:
                        print(x)
                        'W2plus'
                    for x in Wm_int_2_plus:
                        print(x)
                        'W2minus'
                    for x in Wm_int_2_minus:
                        print(x)
                    print('koniec')
                elif theory == 'cc3':
                    Wm_int_1, Wm_int_2_plus, Wm_int_2_minus, Wm_int_3 = \
                       integrate_quadra_Wm_triplet(Wm1_c,  Wm2_c,  Wm3_c, Wm2_mu_list, Wm3_mu_list, \
                                                   theory, multiplicity=multiplicity)
            elif multiplicity == 13:
                if theory == 'ccsd':
                    Wm_int_1, Wm_int_2= \
                       integrate_quadra_Wm_triplet(Wm1_c,  Wm2_c,  Wm3_c, Wm2_mu_list, Wm3_mu_list, \
                                                   theory, multiplicity=multiplicity)

                    print('to o co mi chodzi')
                    print('W1_int')
                    for x in Wm_int_1:
                        print(x)
                    print('')
                    print('W2_int')
                    for x in Wm_int_2:
                        print(x)
                    print('')


                elif theory == 'cc3':
                    Wm_int_1, Wm_int_2, Wm_int_3 = \
                       integrate_quadra_Wm_triplet(Wm1_c,  Wm2_c,  Wm3_c, Wm2_mu_list, \
                                                   Wm3_mu_list, theory, multiplicity=multiplicity)

            elif multiplicity == 31:
                if theory == 'ccsd':
                    Wm_int_1, Wm_int_2_plus, Wm_int_2_minus= \
                       integrate_quadra_Wm_triplet(Wm1_c,  Wm2_c,  Wm3_c, Wm2_mu_list, Wm3_mu_list, \
                                                   theory, multiplicity=multiplicity)


                elif theory == 'cc3 ':
                    Wm_int_1, Wm_int_2_plus, Wm_int_2_minus, Wm_int_3 = \
                       integrate_quadra_Wm_triplet(Wm1_c,  Wm2_c,  Wm3_c, Wm2_mu_list, \
                                                   Wm3_mu_list, theory, multiplicity=multiplicity)

            end11 = time.time()
            print('time elapsed NA CALE=',theory, end11-start11)


            if theory == 'ccsd':
            
                if multiplicity == 1:
                    
                    pickle.dump(Wm_int_1, out11)
                    pickle.dump(Wm_int_2, out22)
                    
                elif multiplicity == 3:
                    pickle.dump(Wm_int_1, out1)
                    pickle.dump(Wm_int_2_plus, out2_plus)
                    pickle.dump(Wm_int_2_minus, out2_minus)

                elif multiplicity == 13:
                    pickle.dump(Wm_int_1, out1so)
                    pickle.dump(Wm_int_2, out2so)

                elif multiplicity == 31:
                    pickle.dump(Wm_int_1, out1)
                    pickle.dump(Wm_int_2_plus, out2_plus)
                    pickle.dump(Wm_int_2_minus, out2_minus)


            elif theory == 'cc3':

                if multiplicity == 1:
                    pickle.dump(Wm_int_1, cc3_out1)
                    pickle.dump(Wm_int_2, cc3_out2)
                    pickle.dump(Wm_int_3, cc3_out3)
                    
                elif multiplicity == 3:
                    pickle.dump(Wm_int_1, cc3_out1trip)
                    pickle.dump(Wm_int_2_plus, cc3_out2ptrip)
                    pickle.dump(Wm_int_2_minus, cc3_out2mtrip)
                    pickle.dump(Wm_int_3, cc3_out3trip)

                elif multiplicity == 13:
                    pickle.dump(Wm_int_1, cc3_out1so)
                    pickle.dump(Wm_int_2, cc3_out2so)
                    pickle.dump(Wm_int_3, cc3_out3so)
                elif multiplicity == 31:
                    pickle.dump(Wm_int_1, cc3_out1_trip)
                    pickle.dump(Wm_int_2_plus, cc3_out2_plus_trip)
                    pickle.dump(Wm_int_2_minus, cc3_out2_minus_trip)
                    pickle.dump(Wm_int_3, cc3_out3_trip)

                        
                
    elif pick == 'load':

        if theory == 'overlap_ccsd' or theory == 'overlap_cc3':
            Wm_oint = pickle.load(outo)
            load_overlap(Wm_oint, theory, multiplicity, maxmbpt)

        else:
            if theory == 'ccsd':
                if multiplicity == 1:
                    Wm_int_1 = pickle.load(out11)
                    Wm_int_2 = pickle.load(out22)
                elif multiplicity == 3:
                    Wm_int_1 = pickle.load(out1)
                    Wm_int_2_plus = pickle.load(out2_plus)
                    Wm_int_2_minus = pickle.load(out2_minus)
                elif multiplicity == 13:
                    Wm_int_1 = pickle.load(out1so)
                    Wm_int_2 = pickle.load(out2so)
                elif multiplicity == 31:
                    Wm_int_1 = pickle.load(out1)
                    Wm_int_2_plus = pickle.load(out2_plus)
                    Wm_int_2_minus = pickle.load(out2_minus)

            elif theory == 'cc3':
                if multiplicity == 1:
                    Wm_int_1 = pickle.load(cc3_out1)
                    Wm_int_2 = pickle.load(cc3_out2)
                    Wm_int_3 = pickle.load(cc3_out3)

                elif multiplicity == 3:
                    Wm_int_1 = pickle.load(cc3_out1trip)
                    Wm_int_2_plus =pickle.load(cc3_out2ptrip)
                    Wm_int_2_minus = pickle.load(cc3_out2mtrip)
                    Wm_int_3 = pickle.load(cc3_out3trip)
                elif multiplicity == 13:
                    Wm_int_1 = pickle.load(cc3_out1so)
                    Wm_int_2 =pickle.load(cc3_out2so)
                    Wm_int_3 = pickle.load(cc3_out3so)
                elif multiplicity == 31:
                    Wm_int_1 = pickle.load(cc3_out1_trip)
                    Wm_int_2_plus = pickle.load(cc3_out2_plus_trip)
                    Wm_int_2_minus = pickle.load(cc3_out2_minus_trip)
                    Wm_int_3 = pickle.load(cc3_out3_trip)
                    

            if multiplicity == 1:
                if theory == 'ccsd':
                    Wm_int_1 = simplify(Wm_int_1)
                    Wm_int_2 = simplify(Wm_int_2)
                    Wm_with_intermediates = append_r1r2_parts(Wm_int_1, Wm_int_2)
                elif theory == 'cc3':

                    Wm_int_1 = simplify(Wm_int_1)
                    Wm_int_2 = simplify(Wm_int_2)
                    Wm_int_3 = simplify(Wm_int_3)

                    Wm_with_intermediates = append_r1r2_parts(Wm_int_1, Wm_int_2, Wm_int_3)

            elif multiplicity == 3:
                if theory == 'ccsd':
                    Wm_int_1 = simplify(Wm_int_1)
                    Wm_int_2_plus = simplify(Wm_int_2_plus)
                    Wm_int_2_minus = simplify(Wm_int_2_minus)
                    Wm_with_intermediates = append_r1r2_parts_triplet(Wm_int_1, Wm_int_2_plus, Wm_int_2_minus)
                elif theory == 'cc3':
                    Wm_int_1 = simplify(Wm_int_1)
                    Wm_int_2_plus = simplify(Wm_int_2_plus)
                    Wm_int_2_minus = simplify(Wm_int_2_minus)
                    Wm_int_3 = simplify(Wm_int_3)
                    Wm_with_intermediates = append_r1r2_parts_triplet(Wm_int_1, Wm_int_2_plus, Wm_int_2_minus, Wm_int_3)

            elif multiplicity == 13:

                if theory == 'ccsd':
                    Wm_int_1 = simplify(Wm_int_1)
                    Wm_int_2 = simplify(Wm_int_2)
                    Wm_with_intermediates = append_r1r2_parts(Wm_int_1, Wm_int_2)

                    print('to o co mi chodzi 2')
                    for x in Wm_with_intermediates:
                        print(x)
                    print('')
                    print('')


                elif theory == 'cc3':
                    Wm_int_1 = simplify(Wm_int_1)
                    Wm_int_2_plus = simplify(Wm_int_2)
                    Wm_int_3 = simplify(Wm_int_3)
                    Wm_with_intermediates = append_r1r2_parts(Wm_int_1, Wm_int_2, Wm_int_3)

            elif multiplicity == 31:

                if theory == 'ccsd':
                    Wm_int_1 = simplify(Wm_int_1)
                    Wm_int_2_plus = simplify(Wm_int_2_plus)
                    Wm_int_2_minus = simplify(Wm_int_2_minus)                    
                    Wm_with_intermediates = append_r1r2_parts_triplet(Wm_int_1, Wm_int_2_plus, Wm_int_2_minus)

                elif theory == 'cc3':
                    Wm_int_1 = simplify(Wm_int_1)
                    Wm_int_2_plus = simplify(Wm_int_2_plus)
                    Wm_int_2_minus = simplify(Wm_int_2_minus)
                    Wm_int_3 = simplify(Wm_int_3)
                    Wm_with_intermediates = append_r1r2_parts_triplet(Wm_int_1, Wm_int_2_plus, Wm_int_2_minus, Wm_int_3)


                    
            if (len(Wm_with_intermediates) == 0):
                print('ret')
                return

            Wm_with_intermediates.clear_fixed()

            print('Wm_with_intermediates', len(Wm_with_intermediates))
            kk = 1
            # for x in Wm_with_intermediates:
            #     print(kk, x, x.coefficient)
            #     kk +=1
            # print('')

            Wm_with_intermediates = simplify(Wm_with_intermediates)

            print('sniez', len(Wm_with_intermediates))
            for x in Wm_with_intermediates:
                print(x)


            # USUWAM WYRAZY disconnected X_ll i X_dd
            print('')
            print('USUWAM WYRAZY disconnected X_ll i X_dd ', len(Wm_with_intermediates))
            Wm_with_intermediates = remove_disconnected(Wm_with_intermediates)
            print('koniec', len(Wm_with_intermediates))
            print('')

            print('TERAZ FAKTORYZUJE')
            k = 1
            kkk = 0
            biglevel_super = []
            bigparameters_super = []
            basket_nointerm = []
            idx_interm = []
            idx_noninterm = []
            idxidx = -1
            for x in Wm_with_intermediates:

                biglevel_best, bigparameters_best, haveinterm = factorize(x)
                idxidx += 1
                if haveinterm == True:
                    idx_interm.append(idxidx)
                    biglevel_super.append(biglevel_best)
                    bigparameters_super.append(bigparameters_best)
                    print(k, x, biglevel_best)
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

            print('TERAZ CZYTAM PIRAMIDY')
            k = 0

            basket_super = []
            interm_dict = {}
            basket_outer_all = []
            basket_idx_fixed = []

            all_hash = {}
            xfx_dict = {}
            interm_fx = {}
            ni = 1
            lenall = 0


            for i in range(0, len(idx_interm)):
                               
                x = Wm_with_intermediates[idx_interm[i]]
                

                worm = []
                outer_worm = []
                idxfx_worm = []
                for minilevel in biglevel_super[i]:
                    ni = read_pyramid(minilevel, x, worm, outer_worm, idxfx_worm, interm_dict, all_hash, xfx_dict, interm_fx, ni)
                if worm == []:
                    print(i, idx_interm[i])
                    sys.exit(0)
                    print('')

                for hh in worm[-1]:
                    print(hh)
                lenall += len(worm[-1])

                basket_super.append(worm)
                basket_outer_all.append(outer_worm)
                print('idxfxworm', idxfx_worm)
                basket_idx_fixed.append(idxfx_worm)
                print(i, 'pirmidocz ', x)
                for gg in worm:
                    for ggg in gg:
                        print('worm', ggg)
                    print('')
                for gg in outer_worm:
                    print('outer_wrom', gg)
                print('')

            print('PIRAMIDY PRZECZYTANE')


            k = 1
            basket_interm = []
            basket_outer = []
            basket_xfx = []
            for x in range(0, len(basket_super)):
                basket_points = [] 
                k += 1
                for y in range(0, len(basket_super[x])):
                    points = 0
                    for z in range(0, len(basket_super[x][y])):
                        points += interm_dict[basket_super[x][y][z].binary_hash]
                        print('SPR POINT', basket_super[x][y][z], basket_idx_fixed[x][y])
                    basket_points.append(points)
                    
                    print('SPR POINT', points)
                print('')

                if basket_points != []:
                    
                    idx = basket_points.index(max(basket_points))
                    basket_interm.append(basket_super[x][idx])
                    basket_xfx.append(basket_idx_fixed[x][idx])
                    for jj in basket_super[x][idx]:
                        print('z jednego kosza dodaje', jj)
                    if len(basket_outer_all[x][idx].coefficient)==0:
                        print('la')
                        sys.exit(0)
                    basket_outer.append(basket_outer_all[x][idx])
                else:
                    print('Wm', Wm_with_intermediates[x])
                    basket_interm.append([])
                    basket_outer.append([])
                    print('nie ma')
                    sys.exit(0)


            n_max = 0
            for i in range(0, len(basket_interm)):
                # print('ten-len', len(basket_interm[i]))
                if (len(basket_interm[i])) > n_max:
                    n_max = deepcopy(len(basket_interm[i]))
                # for j in basket_interm[i]:
                #     print(j)
                # print('')
                    


            print('OSTATECZNE')
            list_of_int= []
            list_of_names = []
            list_of_int_mem = []
            list_of_int_disk = []
            list_of_int_mem_names = []
            list_of_int_disk_names = []
            for n in range(0, n_max):
                list_of_names.append([])
                list_of_int.append([])
                list_of_int_mem.append([])
                list_of_int_disk.append([])
                list_of_int_mem_names.append([])
                list_of_int_disk_names.append([])

            mem_dict = {}
            cost_dict = {}
            lenall2 = 0
            la = 0
            lb = 0
            print('pluszonek')
            for i in range(0, len(basket_interm)):
                whichlevel = 0
                for j in range(0, len(basket_interm[i])):
                    lenall2 += 1
                    print('teraz sprawdzam wyraz', basket_interm[i][j])
                    memidx, costidx = memcost(basket_interm[i][j])
                    memreal = (200**(memidx[0]) * 10**memidx[1])*8/(10**9)
                    print('COSTIDX', costidx)
                    if memreal > 10.0:
                        print(memidx[0], memidx[1])
                        print(basket_interm[i][j], memreal)
                        print('DUZE GOWNO')
                      #  sys.exit(0)
                    mem_dict[basket_interm[i][j].binary_hash] = memreal
                    cost_dict[basket_interm[i][j].binary_hash] = costidx
                    
                    # teraz sprawdzam ktorego rzedy jest dany intermediate
                    for coef in basket_interm[i][j].coefficient:
                        if 'interm' in coef:
                            whichlevel += 1
                            break
                    
                    if basket_interm[i][j] not in list_of_int[whichlevel]:
                        print(i, 'dodaje do', whichlevel, len(list_of_int[whichlevel]))
                        list_of_int[whichlevel].append(basket_interm[i][j])
                        list_of_names[whichlevel].append(all_hash[basket_interm[i][j].binary_hash])
                        print('dodaje', lenall2)
                        la += 1
                    else:
                        lb += 1

                        

            print(lenall2, la, lb)


            print('TERAZ STWORZE BACZE')
            for x in range(0, n_max):
                print('parampam', x, len(list_of_int[x]))
                k = 0
                for y in range(0, len(list_of_int[x])):
                    memidx, costidx = memcost(list_of_int[x][y])
#                    if sum(costidx) >= 7:
                    print(k, list_of_int[x][y], all_hash[list_of_int[x][y].binary_hash], costidx)
                    k+= 1

            print('')

            # extract intermediates of high memory
            # mem_disk= []
            # for x in range(0, n_max):
            #     mem_disk.append([])

            # for x in range(0, n_max):
            #     used_dict_temp = []

            #     print('len(list_of_int[n])', x, len(list_of_int[x]))
            #     mem_diskt = 0
            #     for l in range(0, len(list_of_int[x])):
            #         k = list_of_int[x][l]
            #         name_int = all_hash[k.binary_hash]
            #         mem_int = mem_dict[k.binary_hash]
            #         if k.binary_hash not in used_dict_temp:
            #             if (mem_int) > MEM_LITTLE_THRESH:
            #                 used_dict_temp.append(k.binary_hash)
            #                 list_of_int_disk[x].append(k)
            #                 list_of_int_disk_names[x].append(name_int)
            #                 mem_diskt += mem_int
            #     print('snol', len(list_of_int_disk[x]))
            #     mem_disk[x] = mem_diskt
            #     if mem_disk[x] < MEM_DISK_THRESH:
            #         allmem = True
            #     else:
            #         allmem = False

            # # jesli allmem == True to znaczy ze wszystkie duze zajmuja lacznie ponizej tresholdu i moga tez byc w pamiecie RAM
            # jesli allmem == False to musimy znalec wsystko co ich używa

            # allmem = False


            # large_interm_dict = []
            # if allmem == False:
            #     for x in range(0, n_max):
            #         # petla po wszystkich interm trzymanych na dysku
            #         l = -1
            #         for elem in list_of_int_disk[x]:
            #             names_to_add = []
            #             names_lv = []
            #             added = []
            #             added_lv = []
            #             minidict = {}
            #             l+= 1
            #             k = list_of_int_disk[x][l]
            #             big_name = all_hash[k.binary_hash]
            #             print('big_name', big_name, len(list_of_int_disk[x]))
            #             minidict['big_name'] = big_name
            #             minidict['big_level'] = x
            #             # petla po wszystkic pozostalych interm
                        
            #             for y in range(x+1, n_max):
            #                 ll = -1
            #                 for elem2 in list_of_int[y]:
            #                     ll += 1
            #                     kk = list_of_int[y][ll]
            #                     temp_int_name = all_hash[kk.binary_hash]
            #                     # sprawdzam czy interm z dysku jest uzywamy przez interm z tej petli
            #                     for coef in kk.coefficient:
            #                         # jesli jest uzywany to go dodaje                                    
            #                         if coef == big_name:
            #                             if temp_int_name not in names_to_add:
            #                                 nlv = check_level(temp_int_name,n_max, list_of_names)
            #                                 names_to_add.append(deepcopy(temp_int_name))
            #                                 names_lv.append(nlv)
            #                                 #list_of_int_disk[x].append(kk)
            #                             for coef2 in kk.coefficient:
            #                                 if 'interm' in coef2:
            #                                     if coef2 != big_name:
            #                                         if coef2 not in names_to_add:
            #                                             nlv = check_level(coef2, n_max, list_of_names)
            #                                             names_to_add.append(coef2)
            #                                             names_lv.append(nlv)
            #             minidict['added'] = added
            #             minidict['added_lv'] = added_lv
            #             minidict['names_to_add'] = names_to_add
            #             minidict['names_lv'] = names_lv
            #             large_interm_dict.append(minidict)
        

            # for x in range(0, len(large_interm_dict)):
            #     mini = large_interm_dict[x]
            #     print(mini['big_name'], mini['big_level'])
            #     print('addded')
            #     for y in range(0, len(mini['added'])):
            #         print(mini['added'][y], mini['added_lv'][y])
            #     print('names to add')
            #     for y in range(0, len(mini['names_to_add'])):
            #         print(mini['names_to_add'][y], mini['names_lv'][y])
            #     print('')

            # sys.exit(0)

            # separate low-memory 0 order intermediates
            
            all_interm_dict = {}
            mem_mem= []
            mem_disk = []
            for x in range(0, n_max):
                mem_mem.append([])
                mem_disk.append([])

            for x in range(0, n_max):            
                used_dict_temp = []

                print('len(list_of_int[n])', x, len(list_of_int[x]))
                mem_memt = 0
                mem_diskt = 0
                minidict = {}
                for l in range(0, len(list_of_int[x])):
                    minidict = {}
                    minidict['lv'] = x
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
                        if (mem_int) < MEM_LITTLE_THRESH:
                            minidict['mem'] = 'ram'
                            list_of_int_mem[x].append(k)                        
                            used_dict_temp.append(k.binary_hash)
                            list_of_int_mem_names[x].append(name_int)
                            mem_memt += mem_int
                        else:
                            minidict['mem'] = 'disk'
                            used_dict_temp.append(k.binary_hash)
                            list_of_int_disk[x].append(k)
                            list_of_int_disk_names[x].append(name_int)
                            mem_diskt += mem_int
                    all_interm_dict[name_int] = deepcopy(minidict)
                    all_interm_dict[name_int]['pt'] = maxmbpt
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

            print('dupa-all')
            sram = 1
            for x in all_interm_dict:
                if len(all_interm_dict[x]['xfx']) == 2:
                    sram += 1
                    if all_interm_dict[x]['xfx'][0] in virtual and all_interm_dict[x]['xfx'][1] in virtual:
                        print(x, all_interm_dict[x]['ugg'], all_interm_dict[x]['xfx'], all_interm_dict[x]['name'])


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


            for n in range(0, n_max):
                for y in big_batch_list[n]:
                    for x in y:
                        memidx, costidx = memcost(x)
                        print(n, x, 'DUPOST', costidx)
                    print('')


            print('lup basket outer')
            for x in basket_outer:
                print(x)

            print('')
            print('llub basetk no')
            for x in basket_nointerm:
                print(x)
                basket_outer.append(x)

            
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

            for n in range(0, n_max):
                print('nnn', n)
                unit_decl_temp = ""
                for l in range(0, len(list_of_int_disk[n])):
                    k = list_of_int_disk[n][l]
                    name_int = all_hash[k.binary_hash]
                    if (all_interm_dict[name_int]['mem'] == 'disk'):
                        unit_decl_temp += "integer :: u{name_int}_pt{pt}\n".format(name_int=name_int, pt=maxmbpt)
                unit_decl0.append(unit_decl_temp)
                print(unit_decl0)
                print( '')

            s_decl = ""
            for n in range(0, n_max):
                declname = "lev"+str(n)
                if theory == 'ccsd':
                    s_decl += """use decl_{declname}_interm_ccsd_pt{mbpt}           
                          """.format(declname=declname, mbpt=maxmbpt)
                elif theory == 'cc3':
                    s_decl += """use decl_{declname}_interm_cc3_pt{mbpt}                                                  
                          """.format(declname=declname, mbpt=maxmbpt)
                    

            dname = ''
            for n in range(0, n_max):
                new_decl = True
                declname = "lev"+str(n)
                for x in range(0, len(big_batch_list[n])):
                    if n == 0:
                        unit_decl = unit_decl0[x]
                    else:
                        unit_decl = ""
                    partname = "lev"+str(n)+"_batch"+str(x)
                    new_file = True
                    last_file = True
                    last_decl = False
                    if x ==  len(big_batch_list[n]) - 1:
                        last_decl = True
                    dname = function_template_batch_intermediates(big_batch_list[n][x], x, partname, declname, all_hash, all_interm_dict,\
                                                                      theory, maxmbpt, multiplicity, new_file, last_file, \
                                                                      new_decl, last_decl, s_decl, unit_decl)
                    new_decl = False


            for x in range(0, len(batch_outer)):
                partname = "batch"+str(x)
                new_file = True
                last_file = True
                print('sraka', s_decl)
                block_oo, block_ov, block_vo, block_vv = fill_dm_blocks_wm(batch_outer[x])
                function_template_batch(block_oo, block_ov, block_vo, block_vv, x, \
                                            partname, all_hash, theory, maxmbpt, multiplicity, new_file, last_file, dname, s_decl, all_interm_dict)


            print('KONCZE')
            sys.exit(0)

#             partname = "batch0"
#             new_file_list.append(0)
#             for x in range(0, len(batch_list)):
#                 if x in new_file_list:
#                     partname = "batch"+str(x)
#                     new_file = True
#                 else:
#                     new_file = False
#                 if (x+1) in new_file_list:
#                     last_file = True
#                 else:
#                     last_file = False
#                 print('przekazuje nazawe', x, partname)

#                 block_oo, block_ov, block_vo, block_vv = fill_dm_blocks_wm(batch_outer[x])
#                 dname = function_template_batch_intermediates(batch_list[x], x, partname, all_hash, interm_fx,\
#                                                                   theory, maxmbpt, multiplicity, new_file, last_file)
#                 function_template_batch(block_oo, block_ov, block_vo, block_vv, x, \
#                                             partname, all_hash, theory, maxmbpt, multiplicity, new_file, last_file, dname)
#             print('tralalal')
#             sys.exit(0)
            

#             block_oo, block_ov, block_vo, block_vv = fill_dm_blocks_wm(basket_outer)
            
#             print('')
#             print('block oo')
#             k = 1
#             for x in block_oo:
#                 print(k, x)
#                 k += 1
#             print('')

#             print('')
#             print('block ov')
#             k = 1
#             for x in block_ov:
#                 print(k, x)
#                 k += 1
#             print('')


#             print('')
#             print('block vo')
#             k = 1
#             for x in block_vo:
#                 print(k, x)
#                 k += 1
#             print('')


#             print('')
#             print('block vv')
#             k = 1
#             for x in block_vv:
#                 print(k, x)
#                 k += 1
#             print('')


#             sys.exit(0)


#             print('generate_best', len(Wm_with_intermediates))
#             Wm_intermediates_list, nonintermediates = generate_best_intermediates(Wm_with_intermediates, 'wm', 0)
            
#             print(len(Wm_intermediates_list))

#             k = 1
#             allbigmem = 0
#             allbigcost = 0
#             for x in Wm_intermediates_list:
#                 memidx, costidx = memcost(x['interm1'])
#                 print(k, x['interm1'], memidx, costidx)
#                 if memidx[0] >= 3:
#                     allbigmem += 1
#                     print('MEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEM')
#                 if sum(costidx) > 7:
#                     allbigcost += 1
#                 memidx, costidx = memcost(x['interm2'])
#                 print(k, x['interm2'], memidx, costidx)
#                 if memidx[0] >= 3:
#                     print('MEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEM')
#                     allbigmem += 1
#                 if sum(costidx) > 7:
#                     allbigcost += 1
#                 memidx, costidx = memcost(x['interm3'])
#                 print(k, x['interm3'],  memidx, costidx)
#                 if memidx[0] >= 3:
#                     print('MEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEM')
#                     allbigmem += 1
#                 if sum(costidx) > 7:
#                     allbigcost += 1

#                 k += 1
#             print('allbigmem', allbigmem)
#             print('allbigcost', allbigcost)





#             sys.exit(0)
            
#             # mem_of_all_intermediates(Wm_intermediates_list)

#             # print('generate_best2')
#             # Wm_intermediates_list_2, nonintermediates = generate_best_intermediates(nonintermediates, 'wm', 0, False)

#             print('nonintermediates po best2')
#             y = 0
#             for x in Wm_intermediates_list:
#                 #if x['interm3'].num_factor != 0:
#                 print(y, x['interm2'], x['interm1'], x['interm3'], x['noninterm'], x['original'])
#                 y+=1

#             print('nonintermediates')
#             for x in nonintermediates:
#                 print(x)

#             # Wm_intermediates_list = []
#             # for x in Wm_intermediates_list_1:
#             #     Wm_intermediates_list.append(x)
#             # for x in Wm_intermediates_list_2:
#             #     Wm_intermediates_list.append(x)

#             print('przed simp')
#             for x in Wm_intermediates_list:
#                 szczyna = compute_cost_direct_ugg(x['interm2'])
#                 szczyna1 = compute_cost_direct_ugg(x['interm1'])
#                 szczyna2 = compute_cost_direct_ugg(x['original'])
#                 szczyna3 = compute_cost_direct_ugg(x['noninterm'])
#                 szczyna4 = compute_cost_direct_ugg(x['interm3'])
#                 #if x['interm3'].num_factor !=0:
#                 print(szczyna2, szczyna, szczyna1, szczyna4, szczyna3, x['original'])
#                 #print(szczyna3, x['interm2'], x['interm1'], x['interm3'], x['noninterm'], x['original'])
                

#             print('nonintermediates3')
#             for x in nonintermediates:
#                 print(x)

#             compute_cost_with_intermediates(Wm_intermediates_list)	     

#             for x in Wm_intermediates_list:
#                 print(x['interm2'], x['interm1'], x['interm3'], x['noninterm'])
# #            sys.exit(0)
    

#             print('simplify')
#             intermediates, Wm_with_intermediates, k = \
#                         simplify_intermediates([], Wm_intermediates_list, nonintermediates, 'wm', 0, multiplicity, maxmbpt)

#             print('intermediates')
#             for x in Wm_with_intermediates:
#                 print(x)





# #            sys.exit(0)
#             # print('nonintermediates4')
#             # for x in nonintermediates:
#             #     print(x)

#                 # vegeta
#             print('vegeta')
#             print('Wm_with_intermediates')
#             for x in Wm_with_intermediates:
#                 print(x)
# #            sys.exit(0)
#             print('intermediates')
#             for x in intermediates:
#                 print(x['interm'])

#             mem_of_all_intermediates(intermediates)


#             block_oo, block_ov, block_vo, block_vv = fill_dm_blocks_wm(Wm_with_intermediates)
#             cost_list_inter = compute_cost_interm(intermediates)
#             memory_interm = []
#             for x in intermediates:
#                 vt, oc = find_mem_int(x['interm'])
#                 memory_interm.append([vt, oc])
#             cost_oo, cost_ov, cost_vo, cost_vv = compute_cost_ov(block_oo, block_ov, block_vo, block_vv)


#             print('intermediates')
#             print('')
#             print('\\begin{tabular}{lc}')
#             for x in range(0, len(intermediates)):
#                 print(x, '$',intermediates[x]['interm'],'$&', cost_list_inter[x], memory_interm[x],'\\\\')
#             print('\end{tabular}')
#             print('')
#             print('blok oo')
#             print('')
#             print('\\begin{tabular}{lc}')
#             for x in range(0, len(block_oo)):
#                 print(x, '$',block_oo[x],'$&', cost_oo[x], '\\\\')
#                 print('\end{tabular}')
#             print('')
#             print('blok ov')
#             print('')
#             print('\\begin{tabular}{lc}')
#             for x in range(0, len(block_ov)):
#                 print(x, '$',block_ov[x],'$&', cost_ov[x], '\\\\')
#             print('\end{tabular}')
#             print('')
#             print('blok vo')
#             print('')
#             print('\\begin{tabular}{lc}')
#             for x in range(0, len(block_vo)):
#                 print(x, '$',block_vo[x],'$&', cost_vo[x], '\\\\')
#             print('\end{tabular}')
#             print('')
#             print('blok vv')
#             print('')
#             print('\\begin{tabular}{lc}')
#             for x in range(0, len(block_vv)):
#                 print(x, '$',block_vv[x],'$&', cost_vv[x], '\\\\')
#             print('\end{tabular}')
#             print('')

#             if theory == 'ccsd':
#                 for x in intermediates:
#                     print(x['interm'])
                    
#                 function_template_wm_intermediates(intermediates, 'ccsd', 'wm', multiplicity, maxmbpt)

#                 function_template_wm2(block_oo, block_ov, block_vo, block_vv, 'ccsd', 'wm', multiplicity, maxmbpt)
#             elif theory == 'cc3':
#                 function_template_wm_intermediates(intermediates, 'cc3', 'wm', multiplicity, maxmbpt)
                
#                 function_template_wm2(block_oo, block_ov, block_vo, block_vv, 'cc3', 'wm', multiplicity, maxmbpt)

def extract_int_no(name):
    ns = '1234567890'
    no =""
    for s in name:
        if s in ns:
            no += s

    return no

def check_level(name, n_max, list_of_names):

    for n in range(0, n_max):
        for l in list_of_names[n]:
            if name == l:
                return n

    return -1

def add_terms(minibatch, used_dict, mem_dict, mem, n, n_max, name_int, list_of_int, all_hash, batch_all_names):

    for m in range(n+1, n_max):
        for ll in range(0, len(list_of_int[m])):
            kk = list_of_int[m][ll]
            for mm in range(0, len(kk.coefficient)):
                coef = kk.coefficient[mm]
#                print('coef', coef)
                if coef == name_int:
                    print('tak-coef', coef, name_int, kk)
                    if kk.binary_hash not in used_dict:
                        print('-~~~--')
                        print('dodaje do tego i owego', kk, name_int)
                        minibatch.append(kk)
                        used_dict.append(kk.binary_hash)
                        mem += mem_dict[kk.binary_hash]
                        if (n+1) < n_max:
                            name_int = all_hash[kk.binary_hash]
                            print('wrzucam z taka nazwa w miejscu a', name_int)
                            minibatch, used_dict, mem =  add_terms(minibatch, used_dict, mem_dict, mem, m, n_max, name_int, list_of_int, all_hash, batch_all_names)
                        for nn in range(0, len(kk.coefficient)):
                            coef2 = kk.coefficient[nn]
                            if mm != nn:
                                if 'interm' in coef2:
                                    if coef2 not in batch_all_names:
                                        print('musze jeszcze sprawdzic zatem', coef2)
                                        minibatch, used_dict, mem =  add_terms(minibatch, used_dict, mem_dict, mem, 0, n_max, coef2, list_of_int, all_hash, batch_all_names)
                                        print(kk)
                        print('i teraz jade dalej z', name_int)    
                print('i teraz jade dalej dwa', ll, name_int)
                

                                

    return minibatch, used_dict, mem


def remove_disconnected(Wm_with_intermediates):

    Wm_with_intermediatescp = deepcopy(Wm_with_intermediates)

    Wm_with_intermediates = []
    
    for i in range(0, len(Wm_with_intermediatescp)):
        for j in range(0, len(Wm_with_intermediatescp[i].coefficient)):
            if (Wm_with_intermediatescp[i].coefficient[j] == OBSERVABLE_X) or (Wm_with_intermediatescp[i].coefficient[j] == OBSERVABLE_X_ASYM):
                if Wm_with_intermediatescp[i].coefficient_idx[j][0] != Wm_with_intermediatescp[i].coefficient_idx[j][1]:
                    Wm_with_intermediates.append(Wm_with_intermediatescp[i])
                # else:
                #     print('usuwam disconnected', Wm_with_intermediatescp[i])
                    
    Wm_with_intermediatescp = deepcopy(Wm_with_intermediates)

    Wm_with_intermediates = []

    out = False
    for i in range(0, len(Wm_with_intermediatescp)):
        for j in range(0, len(Wm_with_intermediatescp[i].coefficient)):
            if out == True:
                break
            for k in range(0, len(Wm_with_intermediatescp[i].coefficient)):
                if j != k:
                    set1 = set(Wm_with_intermediatescp[i].coefficient_idx[j])
                    set2 = set(Wm_with_intermediatescp[i].coefficient_idx[k])
                    # print(Wm_with_intermediatescp[i])
                    # print(set1)
                    # print(set2)
                    if set1 == set2:
                        print('usuwam disconnected', Wm_with_intermediatescp[i])
                        out = True
                        break
        if out == False:
            # print('dodaje')                                                                                                                
            Wm_with_intermediates.append(Wm_with_intermediatescp[i]) 
        else:
            out = False

                    # else:
                    #     print('dodaje')
                    #     Wm_with_intermediates.append(Wm_with_intermediatescp[i])
                    #     out = True
                    #     break

    return Wm_with_intermediates


























































    
#--------------------------------------------------------------------------------------------------------------------------------------
#     intermediates_dict_1, Wm_with_intermediates_1, sum_list_1, k = generate_intermediates(Wm_with_intermediates_r1, 'wm', 0, 0)
#     intermediates_dict_2, Wm_with_intermediates_2, sum_list_2, k = generate_intermediates(Wm_with_intermediates_r2, 'wm', k, 0)

#     intermediates_dict_1b, Wm_with_intermediates_1, sum_list_1, k = generate_intermediates(Wm_with_intermediates_1, 'wm', k, 1)
#     intermediates_dict_2b, Wm_with_intermediates_2, sum_list_2, k = generate_intermediates(Wm_with_intermediates_2, 'wm', k, 1)

#     k = 0
#     l = 0
#     for x in intermediates_dict_1:
#         ll = 0
#         for y in x['idx_fx']:
#             if y in virtual:
#                 ll += 1
#         if ll == 2:
#             l += 1

#         if len(x['idx_fx']) >= 4:
#                print(x['int_name'], x['interm'], x['idx_fx'])
#                k += 1
#     for x in intermediates_dict_2:
#         ll = 0
#         for y in x['idx_fx']:
#             if y in virtual:
#                 ll += 1
#         if ll == 2:
#             l += 1
                
#         if len(x['idx_fx']) >= 4:
#                print(x['int_name'], x['interm'], x['idx_fx'])
#                k += 1

#     print('k', k, l)

#     k = 0
#     l = 0
#     for x in intermediates_dict_1b:
#         ll = 0
#         for y in x['idx_fx']:
#             if y in virtual:
#                 ll += 1
#         if ll == 2:
#             l += 1

#         if len(x['idx_fx']) >= 4:
#                print(x['int_name'], x['interm'], x['idx_fx'])
#                k += 1
#     for x in intermediates_dict_2b:
#         ll = 0
#         for y in x['idx_fx']:
#             if y in virtual:
#                 ll += 1
#         if ll == 2:
#             l += 1
                
#         if len(x['idx_fx']) >= 4:
#                print(x['int_name'], x['interm'], x['idx_fx'])
#                k += 1

#     print('k', k, l)


#     print('')
#     print('INTERMEDIATES')
#     print('')
#     cost_list = []
#     cost_old = 0
#     i = 0
#     for x in Wm_with_intermediates_1:
#         i += 1
#         vidx = 0
#         oidx = 0
#         for z in list(x.summation):
#             if z in virtual:
#                 vidx += 1
#             elif  z in occupied:
#                 oidx += 1
#         cost_list.append([vidx, oidx])
#         if vidx == 5:
#             print(i, [vidx, oidx], x)
#         cost_new = vidx * 20 + oidx * 0.1
#         if cost_new > cost_old:
#             cost_old = cost_new
#             cs = [vidx, oidx]
#             sniez = x
#             poln = i
#     print(cost_old, cs, sniez, poln)



#     print('')
#     cost_list = []
#     cost_old = 0
#     i = 0
#     for x in Wm_with_intermediates_2:
#         i += 1
#         vidx = 0
#         oidx = 0
#         for z in list(x.summation):
#             if z in virtual:
#                 vidx += 1
#             elif  z in occupied:
#                 oidx += 1
#         cost_list.append([vidx, oidx])
#         if vidx == 5:
#             print(i, [vidx, oidx], x)
#         cost_new = vidx * 20 + oidx * 0.1
#         if cost_new > cost_old:
#             cost_old = cost_new
#             cs = [vidx, oidx]
#             sniez = x
#             poln = i
#     print(cost_old, cs, sniez, poln)
#     sys.exit(0)


#     Wm_with_intermediates_r1 = arithmetic_string()
#     for x in Wm_with_intermediates_1:
#         Wm_with_intermediates_r1.append(x)

#     Wm_with_intermediates_r2 = arithmetic_string()
#     for x in Wm_with_intermediates_2:
#         Wm_with_intermediates_r2.append(x)

# #---------------------------------------------------------------------------------------------------------------------------

    
#     print('simpl1')
#     Wm_int_1 = simplify(Wm_int_11 + Wm_int_12)
#     print('simpl2')
#     Wm_int_2 = simplify(Wm_int_21 + Wm_int_22)

#     ln1 = len(Wm_int_1)
#     ln2 = len(Wm_int_2)

#     Wm_int_sum = Wm_int_1 + Wm_int_2

#     cost_list = []
#     cost_old = 0
#     i = 0
#     for x in Wm_int_sum:
#         i += 1
#         if 'b' in list(x.summation):
#             vidx = 1
#             oidx = 1
#         else:
#             vidx = 2
#             oidx = 2
#         for z in list(x.summation):
#             if z in virtual:
#                 vidx += 1
#             elif  z in occupied:
#                 oidx += 1
#         cost_list.append([vidx, oidx])
#         if vidx == 5:
#             print(i, [vidx, oidx], x)
#         cost_new = vidx * 20 + oidx * 0.1
#         if cost_new > cost_old:
#             cost_old = cost_new
#             cs = [vidx, oidx]
#             sniez = x
#             poln = i
#     print(cost_old, cs, sniez, poln)
#     #sys.exit(0)
   
#     intermediates_dict, Wm_with_intermediates, sum_list = generate_intermediates(Wm_int_sum, 'wm')

#     for x in intermediates_dict:
#         print(x['int_name'], x['interm'])

#     cost_list = []
#     cost_old = 0
#     i = 0
#     for x in Wm_with_intermediates:
#         i += 1
#         if 'b' in list(x.summation):
#             vidx = 1
#             oidx = 1
#         else:
#             vidx = 2
#             oidx = 2
#         for z in list(x.summation):
#             if z in virtual:
#                 vidx += 1
#             elif  z in occupied:
#                 oidx += 1
#         cost_list.append([vidx, oidx])
#         if vidx == 5:
#             print(i, [vidx, oidx], x)
#         cost_new = vidx * 20 + oidx * 0.1
#         if cost_new > cost_old:
#             cost_old = cost_new
#             cs = [vidx, oidx]
#             sniez = x
#     print(cost_old, cs, sniez)

# #    sys.exit(0)

#     print('Wm_int_1')
#     Wm_int_1 = arithmetic_string()
#     for x in range(0, ln1):
#         Wm_int_1.append(Wm_with_intermediates[x])
#         # vidx = 0
#         # oidx = 0
#         # z = list(Wm_with_intermediates[x].summation)
#         # for y in z:
#         #     if y in virtual:
#         #         vidx += 1
#         #     elif y in occupied:
#         #         oidx += 1

#         # vidx += 1
#         # oidx += 1
#         # s = """v^{vidx}o^{oidx}""".format(vidx=vidx, oidx=oidx)
#         # if vidx >= 4:
#         #     print(s, vidx+oidx, Wm_with_intermediates[x])

#     print('Wm_int_2')
#     Wm_int_2 = arithmetic_string()
#     for x in range(ln1, ln1+ln2):
#         Wm_int_2.append(Wm_with_intermediates[x])
#         # vidx = 0
#         # oidx = 0
#         # z = list(Wm_with_intermediates[x].summation)
#         # for y in z:
#         #     if y in virtual:
#         #         vidx += 1
#         #     elif y in occupied:
#         #         oidx += 1

#         # vidx += 1
#         # oidx += 1
#         # s = """v^{vidx}o^{oidx}""".format(vidx=vidx, oidx=oidx)
#         # if vidx >= 4:
#         #     print(s, vidx+oidx, Wm_with_intermediates[x])


# #    sys.exit(0)

#     # print('simplu1')
#     # Um_int_1 = simplify(Um_int_11 + Um_int_12)
#     # print('simplu2')
#     # Um_int_2 = simplify(Um_int_21 + Um_int_22)

#     # ln1u = len(Um_int_1)
#     # ln2u = len(Um_int_2)

#     # Um_int_sum = Um_int_1 + Um_int_2

#     # Um_intermediates_dict, Um_with_intermediates, Um_sum_list2 = generate_intermediates(Um_int_sum, 'um')

#     # Um_int_1 = arithmetic_string()
#     # for x in range(0, ln1u):
#     #     Um_int_1.append(Um_with_intermediates[x])

#     # Um_int_2 = arithmetic_string()
#     # for x in range(ln1u, ln1u+ln2u):
#     #     Um_int_2.append(Um_with_intermediates[x])

#     #
#     # Divide into r1 and r2 parts. R1 part is merged with R1 amplitudes and R2 part is merged with
#     # R2 amplitudes.

#     Wm_with_intermediates_r1, Wm_with_intermediates_r2 = divide_to_r1r2_parts_2(Wm_int_1, Wm_int_2)

#     for x in Wm_with_intermediates_r2:
#         print(x)
#     sys.exit(0)


    # Um_with_intermediates_r1, Um_with_intermediates_r2 = divide_to_r1r2_parts_2(Um_int_1, Um_int_2)


    # wl_block_oo_diag, wl_block_vv_diag, wl_block_oo, wl_block_ov, wl_block_vo, wl_block_vv = fill_dm_blocks_wm(Wl_int)
    # wr_block_oo_diag, wr_block_vv_diag, wr_block_oo, wr_block_ov, wr_block_vo, wr_block_vv = fill_dm_blocks_wm(Wr_int)



    # Um_r1_block_oo_diag, Um_r1_block_vv_diag, Um_r1_block_oo, \
    #     Um_r1_block_ov, Um_r1_block_vo, Um_r1_block_vv = fill_dm_blocks_wm(Um_with_intermediates_r1)
    # Um_r2_block_oo_diag, Um_r2_block_vv_diag, Um_r2_block_oo, \
    #     Um_r2_block_ov, Um_r2_block_vo, Um_r2_block_vv = fill_dm_blocks_wm(Um_with_intermediates_r2)


    # function_template_wl(wl_block_oo, wl_block_ov, wl_block_vo, wl_block_vv, 'ccsd')

    # r1_block_oo_diag, r1_block_vv_diag, r1_block_oo, r1_block_ov, r1_block_vo, r1_block_vv = fill_dm_blocks_wm(Wm_with_intermediates_r1)
    # r2_block_oo_diag, r2_block_vv_diag, r2_block_oo, r2_block_ov, r2_block_vo, r2_block_vv = fill_dm_blocks_wm(Wm_with_intermediates_r2)

    # function_template_wm_intermediates(intermediates_dict, 'ccsd', 'wm')

    # function_template_wm(r1_block_oo_diag, r1_block_oo, r1_block_ov, r1_block_vo, r1_block_vv, \
    #                          r2_block_oo_diag, r2_block_oo, r2_block_ov, r2_block_vo, r2_block_vv, 'ccsd', 'wm')


    # # function_template_wm_intermediates(Um_intermediates_dict_2, 'ccsd', 'um')

    # # function_template_wm(Um_r1_block_oo_diag, Um_r1_block_oo, Um_r1_block_ov, Um_r1_block_vo, Um_r1_block_vv, \
    # #                          Um_r2_block_oo_diag, Um_r2_block_oo, Um_r2_block_ov, Um_r2_block_vo, Um_r2_block_vv, 'ccsd', 'um')



    # sys.exit(0)

    
