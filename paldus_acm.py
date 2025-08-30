from params import *
import paldus_classes
from paldus_cas import cas
from paldus_cas import swap_count
from paldus_cas import cas_to_ugg
from paldus_cas import ugg_to_cas
from paldus_cas import add_spin_from_list
from paldus_cas import *
from paldus_basic import *
from paldus_classes import ugg
from copy import deepcopy
from paldus_classes import arithmetic_string
from paldus_classes import disambiguate
from collections import deque
from collections import Counter
from paldus_classes import pair_permutations
import math
from itertools import product
from itertools import permutations
from multiprocessing import Pool
from paldus_classes import integrate
from paldus_classes import virtual, occupied, general
from fortran_code import *
from templates import *
import sys
import time
import io
from eomccjac import jacobian_loop
from paldus_cas import are_spins_equal
import pickle
import time





def hh_rpa_koopmans():

    # for HH
    # r+ = AAr_cc
    # p = AAp_hh

    #KL = evaluate(AAr_cc, AAp_hh)

    # KL = evaluate(h1cas, AAp_hh, AAr_cc).scale(-1.0) +  evaluate(h2cas, AAp_hh, AAr_cc).scale(-1.0)
    
    # KL = evaluate(AAr_cc, h1cas, AAp_hh) +  evaluate(AAr_cc, h2cas, AAp_hh)
    

    KL = evaluate(AAr_cc, h1cas, AAp_hh).scale(0.5) +  evaluate(AAr_cc, h2cas, AAp_hh).scale(0.5) +\
        evaluate(h1cas, AAp_hh, AAr_cc).scale(-0.5) +  evaluate(h2cas, AAp_hh, AAr_cc).scale(-0.5)

    # kl = evaluate(AAp_hh, h1cas, AAr_cc).scale(0.5) +  evaluate(AAp_hh, h2cas, AAr_cc).scale(0.5) +\
    #     evaluate(h1cas, AAr_cc, AAp_hh).scale(-0.5) +  evaluate(h2cas, AAr_cc, AAp_hh).scale(-0.5)


    # KL = evaluate(AAp_hh, AAr_cc)

    print('wynik komutatorow KL')

    k = 0
    for x in KL:
        print(k, x)
        k += 1
    print('')

    print('Perform Wick')
    res = []
    for x in KL:
        result_list = x.wick_ca()
        res = res + result_list
    print('')

    
    print('Result po WICK:')    
    for x in res:
        x.exec_delta()
        if (x.num_factor) != 0.0:
            print(x)
        print(x)
    print('')

    print('Result po RENAME DENSITY:')    
    for x in res:
        x.rename_as_density()
        if (x.num_factor) != 0.0:
            print(x)
    print('')

    res2 = cas_to_ugg(res)
    print('po tranformacji')
    for x in res2:
        print(x)
    
    rsimp = simplify(res2, cas=True)
    k = 0
    print('po simp')
    for x in rsimp:
        k += 1
        # print("&", x, "\\\\")
        print(x)

    rsimpa = find_antysym_int(rsimp)
    k = 0
    print('')
    print('po simp antysym')
    print('')
    for x in rsimpa:
        k += 1
        # print("&", x, "\\\\")
        print(x)


def int_and_simp_cas(KL):

    print('wynik komutatorow KL')
    k = 0
    for x in KL:
        print(k, x)
        k += 1
    print('')

    print('Perform Wick')
    res = []
    for x in KL:
        result_list = x.wick_ca()
        res = res + result_list
    print('')

    
    print('Result po WICK:')    
    for x in res:
        x.exec_delta()
        if (x.num_factor) != 0.0:
            print('aa', x)

    print('')
#    sys.exit(0)
    # for x in res:
    #     x.exec_delta()

    print('Result po RENAME DENSITY:')    
    for x in res:
        x.rename_as_density()
        if (x.num_factor) != 0.0:
            print(x)
    print('')

    print('Result po RENAME DENSITY:')    
    for x in res:
        print(x)

    # tu bedzie przyblizenie Gamma_3 i Gamma4 przez 1-RDMy

    # dupa1 = res[0]
    # print(dupa1)
    # generate_1rdm_prod(dupa1, [1], [3])
    # sys.exit(0)


    
    res = approximate_rdm34_by_rdm1(res, True)
    print('teraz przyblizam rdm3 by rdm12')
    res = approximate_rdm3_by_rdm12(res)

    # res = approximate_rdm34_by_rdm12(res)

    print('po apporox')
    for x in res:
        print(x)
#    sys.exit(0)
    # for x in res:
    #     x.rename_as_density()

    print('cas')
    for x in res:
        print(x)
    print('ugg')
    res2 = cas_to_ugg(res)
    print('cas to ugg done')
    # for x in res2:
    #     print(x)
#    sys.exit(0)
    print('stand')
    for x in res2:
        x.standarize()
    print('')
    for x in res2:
        print(x)
#        print(x.standarize())
        

    #sys.exit(0)
    # res2kl1 = cas_to_ugg(reskl1)

    # print('po tranformacji')
    # for x in res2:
    #     print(x)

    print('przed simplify', len(res2))
    rsimp = simplify(res2, cas=True)
    print('po simplify', len(rsimp))
    for x in rsimp:
        print(x)
    print('')
    occup = ['i', 'j', 'k', 'l']
    for x in rsimp:
        nocc = 0
        isgam = False
        for i in range(0, len(x.coefficient)):
            coef = x.coefficient[i]
            if coef == DENS3 or coef == DENS4:
                isgam = True
                # print(x, 'ano')
                for j in x.coefficient_idx[i]:
                    # print('jjj', j, occup)
                    if j in occup:
                        nocc += 1
        if isgam == True:
            print(x, nocc)
        else:
            print(x)

    spin_dict1={}
    spin_dict1['i'] = '+'
    spin_dict1['j'] = '-'
    spin_dict1['r'] = '+'
    spin_dict1['s'] = '-'
    spin_dict1['q'] = '+'
    spin_dict1['p'] = '-'
    spin_dict1['k'] = '+'
    spin_dict1['l'] = '-'

    numf1 = 1.0

    # spin_dict2={}
    # spin_dict2['a'] = '+'
    # spin_dict2['b'] = '-'
    # spin_dict2['c'] = '-'
    # spin_dict2['d'] = '+'
    # numf2 = -1.0

    res3 = add_spin_driver(spin_dict1, numf1, rsimp)
    print('')
    print('i koniec')
    occup = ['i', 'j', 'k', 'l']
    for x in res3:
        nocc = 0
        isgam = False
        for i in range(0, len(x.coefficient)):
            coef = x.coefficient[i]
            if coef == DENS3 or coef == DENS4:
                isgam = True
                # print(x, 'ano')
                for j in x.coefficient_idx[i]:
                    # print('jjj', j, occup)
                    if j in occup:
                        nocc += 1
        if isgam == True:
            print(x, nocc)
        else:
            print(x)

        
    # rsimpkl1 = simplify(res2kl1, cas=True)

    return res3
#    return rsimp




def pphh_rpa_read():

    res = []

    f = open("wynik.dat", "r")
    lines = f.readlines()

    for line in lines:
#        print('llline', line)
        x = read_cas_from_str(line)
        res.append(x)
        old_stdout = sys.stdout
        new_stdout = io.StringIO()
        sys.stdout = new_stdout
        print(x)
        output = new_stdout.getvalue()
        sys.stdout = old_stdout
        if output != line:
            print('rozne')
            print('l', line)
            print('o', output)
            print('')
#        print('x', x)
        # print('')
    print('')
#    sys.exit(0)
    print('wynik')
    for x in res:
        print(x)
    print('')
    print('teraz inne')

    for i in range(0, len(res)):
        for j in range(0, len(res[i].coefficient)):
            if res[i].coefficient[j] == DENS2:
                w = are_spins_equal(res[i].coefficient_spin[j])
                if not w:
                    print('taak', res[i])
                    res[i].num_factor = 0.0
            if res[i].coefficient[j] == DENS1:
                if res[i].coefficient_spin[j][0] !=res[i].coefficient_spin[j][1]:
                    res[i].num_factor = 0.0
                if res[i].num_factor != 0.0:
                    res[i].new_delta(res[i].coefficient_idx[j][0], res[i].coefficient_idx[j][1], \
                                     res[i].coefficient_spin[j][0], res[i].coefficient_spin[j][1])
                    res[i].coefficient[j] = 'n'
                    res[i].coefficient_idx[j] = [res[i].coefficient_idx[j][0]]
        res[i].exec_delta()

    res = clean_cas(res)


    print('wynik')
    for x in res:
        print(x)

    print('only terms without Gamma2')
    res_1g, rest = res_1gam_only(res)

    res_pppp, res_pmpm, res_mpmp = res_divied_gammas(rest)

    print('lenlen', len(res_1g), len(res_pppp), len(res_pmpm))
    for x in res_1g:
        print(x)
    print('')
    print('only terms with Gamma_pppp')
    print('')
    for x in res_pppp:
        print(x)
    # print('only terms with Gamma_mmmm')
    # print('')
    # for x in res_mmmm:
    #     print(x)
    print('only terms with Gamma_pmpm')
    print('')
    for x in res_pmpm:
        print(x)
    print('only terms with Gamma_mpmp')
    print('')
    for x in res_mpmp:
        print(x)

    res_1g_ugg = cas_to_ugg(res_1g)
    res_pppp_ugg = cas_to_ugg(res_pppp)
#    res_mmmm_ugg = cas_to_ugg(res_mmmm)
    res_pmpm_ugg = cas_to_ugg(res_pmpm)


    res_1gs = simplify(res_1g_ugg, cas=True)
    res_pppps = simplify(res_pppp_ugg, cas=True)
#    res_mmmms = simplify(res_mmmm_ugg, cas=True)
    res_pmpms = simplify(res_pmpm_ugg, cas=True)


    print('')
    print('terms with DENS1', len(res_1gs))
    for x in res_1gs:
        print(x)
    print('')
    print('terms with DENS2_++++', len(res_pppps))
    for x in res_pppps:
        print(x)
    # print('')
    # print('terms with DENS2_----', len(res_mmmms))
    # for x in res_mmmms:
    #     print(x)
    print('')
    print('terms with DENS2_+-+-', len(res_pmpms))
    for x in res_pmpms:
        print(x)
    print('')


    sys.exit(0)
        
def pphh_rpa():

    KL = evaluate(AAtuvw, AAqprs)
    #KL = evaluate(AArs_aa, AAqp_pp)
    #KL = evaluate(AArp, AAqs)
        
    #KL = evaluate(AAtuvw, h1cas, AAqprs) +  evaluate(AAtuvw, h2cas, AAqprs)
    res = int_and_simp_cas(KL)
    

    for i in range(0, len(res)):
        for j in range(0, len(res[i].coefficient)):
            if res[i].coefficient[j] == DENS2:
                w = are_spins_equal(res[i].coefficient_spin[j])
                if not w:
                    print('taak', res[i])
                    res[i].num_factor = 0.0
            if res[i].coefficient[j] == DENS1:
                if res[i].coefficient_spin[j][0] !=res[i].coefficient_spin[j][1]:
                    res[i].num_factor = 0.0
                if res[i].num_factor != 0.0:
                    res[i].new_delta(res[i].coefficient_idx[j][0], res[i].coefficient_idx[j][1], \
                                     res[i].coefficient_spin[j][0], res[i].coefficient_spin[j][1])
                    res[i].coefficient[j] = 'n'
                    res[i].coefficient_idx[j] = [res[i].coefficient_idx[j][0]]
        res[i].exec_delta()

    res = clean_cas(res)
                                
    print('wynik')
    for x in res:
        print(x)


    
    sys.exit(0)
    k = 0
    print('po simp')
    for x in rsimp:
        k += 1
        # print("&", x, "\\\\")
        print(x)

def int_and_simp_cas_onedet(KL):

    
    print('wynik komutatorow KL')
    k = 0
    for x in KL:
        print(k, x)
        k += 1
    print('')

    print('Perform Wick')
    res = []
    for x in KL:
        result_list = x.wick_ca()
        res = res + result_list
    print('')

    
    print('Result po WICK:')
    k = 1
    for x in res:
        x.exec_delta()
        if (x.num_factor) != 0.0:
            print('aa', k,  x)
            k+=1

    print('')


    return res

def pphh_rpa_S_onedet():

    AA = cas()
    AA.operator_idx = ['c', 'd', 'k', 'l']
    AA.operator_type = ['+', '+', '0', '0']

    BB = cas()
    BB.operator_idx = ['i', 'j', 'a', 'b']
    BB.operator_type = ['+', '+', '0', '0']


    
    KL = evaluate(BB, AA)

    res = int_and_simp_cas_onedet(KL)
    
    res = zero_dens_matrix_onedet(res)
    res = clean_cas(res)

    print('Result po RENAME DENSITY:')    
    for x in res:
        x.rename_as_density()

    print('')


    print('')
    print('wynik')
    print('')
    for x in res:
        print(x)

    res2 = cas_to_ugg(res)    
    rsimp = simplify(res2, cas=True)

    print('')
    print('wynik', len(rsimp))
    print('')
    for x in rsimp:
        print(x)

    res = ugg_to_cas(rsimp)
    res = approximate_rdm3_by_rdm12(res)
    res = approximate_rdm2_by_rdm1_when_occ(res, True, True)

    # res = gamma_to_n(res)
    # res = clean_cas(res)

    # resugg = cas_to_ugg(res)
    # rsimp = simplify(resugg)
    # res = ugg_to_cas(rsimp)

    # print('')
    # print('hupla2', len(res))
    # print('')
    # for x in res:
    #     print(x)
    # print('')

    # sys.exit(0)


    #AA +-+-  +-+-
    #AB +-+-  +--+
    #AC +-+-  -++-
    #AD +-+-  -+-+

    #BA +--+  +-+-
    #BB +--+  +--+
    #BC +--+  -++-
    #BD +--+  -+-+

    spin_dict1={}
    spin_dict1['i'] = '+'
    spin_dict1['j'] = '-'
    spin_dict1['a'] = '-'
    spin_dict1['b'] = '+'
    spin_dict1['c'] = '-'
    spin_dict1['d'] = '+'
    spin_dict1['k'] = '-'
    spin_dict1['l'] = '+'

    numf1 = 1.0
    res3_abab = add_spin_driver(spin_dict1, numf1, res)
    # res4 = cas_to_ugg(res)
    # rsimp = simplify(res3_abab)

    # res2 = cas_to_ugg(res)    
    # rsimp = simplify(res2, cas=True)

    res = clean_cas(res3_abab)
    resugg = cas_to_ugg(res)
    rsimp = simplify(resugg)
    res = ugg_to_cas(rsimp)

    res = gamma_to_n(res, True)

    # for i in range(0, len(res)):
    #     for j in range(0, len(res[i].coefficient)):
    #         if res[i].coefficient[j] == DENS1:
    #             if res[i].num_factor != 0.0:
    #                 res[i].new_delta(res[i].coefficient_idx[j][0], res[i].coefficient_idx[j][1])
    #                 res[i].coefficient[j] = 'n'
    #                 res[i].coefficient_idx[j] = [res[i].coefficient_idx[j][0]]
    #                 res[i].num_factor *= 0.5
    #     res[i].exec_delta()

    # res = clean_cas(res)


    print('')
    print('wynik', len(res))
    print('')
    for x in res:
        print(x)


    sys.exit(0)


def test_rdm_approx():
    
    res = []

    f = open("rdm4.dat", "r")
    lines = f.readlines()

    for line in lines:
        x = read_cas_from_str(line, False)
        res.append(x)
        old_stdout = sys.stdout
        new_stdout = io.StringIO()
        sys.stdout = new_stdout
        print(x)
        output = new_stdout.getvalue()
        sys.stdout = old_stdout
        if output != line:
            print('rozne')
            print('l', line)
            print('o', output)
            print('')


    print('wynik')
    for x in res:
        print(x)
    print('')

    res_big = approximate_rdm4(res, False)
    print('wynik_big')
    for x in res_big:
        print(x)
    print('')
    nospin = True
    res3 = cum2_to_rdm2(res_big, nospin, False)
    print('wynik 3')
    for x in res3:
        print(x)
        # print(x.coefficient_idx)
        # print('')
    print('')

    res2 = cas_to_ugg(res3)
    #    sys.exit(0)
    print('ugg')
    for x in res2:
        print(x)
    print('')

    rsimp = simplify(res2, cas=True)

    print('')
    print('wyniker')
    print('')
    for x in rsimp:
        print(x)

    

    sys.exit(0)

    res6 = arithmetic_string()
    for x in res:
        res6 = generate_partition(x, 0, 4, 2, 2, False)
    print('')

    res2 = arithmetic_string()
    res2.append(res6[0])
    print('res2', res2)
    res3 = cum2_to_rdm2(res2, True, False)
    
    # res2 = approximate_rdm3_by_rdm12(res)
    print('wynik2res3')
    for x in res3:
        print(x)
    print('')

    sys.exit(0)

def TDA_pphh_ph_onedet(block, op, filename = '', read = False):

        
    if block == 'pp':
        AA = cas()
        AA.operator_idx = ['c', 'd']
        AA.operator_type = ['+', '+']
        
        BB = cas()
        BB.operator_idx = ['b', 'a']
        BB.operator_type = ['0', '0']
        

    if block == 11:
        # AA = cas()
        # AA.operator_idx = ['c', 'k']
        # AA.operator_type = ['+', '0']
        
        # BB = cas()
        # BB.operator_idx = ['i', 'a']
        # BB.operator_type = ['+', '0']

        AA = cas()
        AA.operator_idx = ['c', 'k']
        AA.operator_type = ['+', '0']
        
        BB = cas()
        BB.operator_idx = ['i', 'a']
        BB.operator_type = ['+', '0']

    elif block == 0:
        BB = cas()
        AA = cas()
    elif block == 2:
        BB = cas()

        AA = cas()
        AA.operator_idx = ['c', 'd', 'k', 'l']
        AA.operator_type = ['+', '+', '0', '0']
    elif block ==20:
        BB = cas()
        BB.operator_idx = ['l', 'k', 'd', 'c']
        BB.operator_type = ['+', '+', '0', '0']
        AA = cas()
    elif block ==12:
        BB = cas()
        BB.operator_idx = ['i', 'a']
        BB.operator_type = ['+', '0']
        
        AA = cas()
        AA.operator_idx = ['c', 'd', 'k', 'l']
        AA.operator_type = ['+', '+', '0', '0']
    elif block ==21:
        BB = cas()
        BB.operator_idx = ['l', 'k', 'd', 'c']
        BB.operator_type = ['+', '+', '0', '0']
        AA = cas()
        AA.operator_idx = ['a', 'i']
        AA.operator_type = ['+', '0']
    elif block ==22:

        AA = cas()
        AA.operator_idx = ['c', 'd', 'k', 'l']
        AA.operator_type = ['+', '+', '0', '0']

        BB = cas()
        BB.operator_idx = ['j', 'i', 'b', 'a']
        BB.operator_type = ['+', '+', '0', '0']


    if read:
        res = []

        f = open(filename, "r")
        lines = f.readlines()

        for line in lines:
            x = read_cas_from_str(line, False)
            res.append(x)
            old_stdout = sys.stdout
            new_stdout = io.StringIO()
            sys.stdout = new_stdout
            print(x)
            output = new_stdout.getvalue()
            sys.stdout = old_stdout
            if output != line:
                print('rozne')
                print('l', line)
                print('o', output)
                print('')
#            print('')

        print('wynik')
        for x in res:
            print(x)
        print('')

    else:
        
        # kod na S
        if op == 2:
            r1 = BB.fromright(AA)
            KL = arithmetic_string(r1)
        elif op == 1:
            # kod na H
            r1 = h1cas.fromleft(BB)
            r2 = h2cas.fromleft(BB)
            r1 = r1.fromright(AA)
            r2 = r2.fromright(AA)

            KL = arithmetic_string(r1)
            KL.append(r2)


        elif op ==3:
            KL = evaluate(BB, h1cas, AA).scale(0.5) +  evaluate(BB, h2cas, AA).scale(0.5) +\
                 evaluate(h1cas, AA, BB).scale(-0.5) +  evaluate(h2cas, AA, BB).scale(-0.5)
        elif op ==4:
            KL = evaluate(BB, AA)

        print('klkl')
        for x in KL:
            print(x)
        #    sys.exit(0)

        res = int_and_simp_cas_onedet(KL)
#        res = arithmetic_string(res[8])
        print('po int and simp')
        for x in res:
            print(x)
        print('')
        res = zero_dens_matrix_onedet(res)
        res = clean_cas(res)

        print('po zero dens onedet')
        for x in res:
            print(x)
        print('')


        print('Result po RENAME DENSITY:')    
        for x in res:
            x.rename_as_density()
            print(x)
        print('')

        res2 = cas_to_ugg(res)    
        rsimp = simplify(res2, cas=True)

        print('')
        print('wynik')
        print('')
        for x in rsimp:
            print(x)

        res = ugg_to_cas(rsimp)

        res = approximate_rdm3(res, True)
        res = approximate_rdm4(res, True)
        res = approximate_rdm2_by_rdm1_when_occ(res, True, True)


        
        # res = approximate_rdm34_by_rdm1(res, False)
        # print('hla')
        # for x in res:
        #     print(x)
        # res = approximate_rdm3_by_rdm12(res)
        # print( 'sra')
        # for x in res:
        #     print(x)

        # res = cum2_to_2rdm(res)
        # res = approximate_rdm2_by_rdm1_when_occ(res, True, True)

        print('hupla1')
        for x in res:
            print(x)

        print('')


    # res = gamma_to_n(res)
    
    # res = clean_cas(res)
    # res = approximate_rdm34_by_rdm1(res, True)

    print('wyn')
    for x in res:
        print(x)
    print('')
#    sys.exit(0)

    # print('approx')
    # i = 1
    # for x in res:
    #     print(i, x)
    #     i+=1
    # print('')





    resugg = cas_to_ugg(res)
    rsimp = simplify(resugg)
    res = ugg_to_cas(rsimp)
    for x in res:
        print('tra', x)
        x.exec_delta()
        print(x)



    print('')
    print('hupla2', len(res))
    print('')
    # res2 = arithmetic_string()
    # i = 0
    # for x in res:
    #     if (i==30):
    #         res2.append(x)
    #     if (i==39):
    #         res2.append(x)
    #     if (i==209):
    #         res2.append(x)
    #     if (i==215):
    #         res2.append(x)
                        
    #     i+=1
    #     print(x)
    # res = res2
    print( 'hupla3')
    for x in res:
        print(x)
    
    print('')

    # sys.exit(0)
    spin_dict_list=[]
    if block == 'pp':
        b = BB.operator_idx[0]
        a = BB.operator_idx[1]
        c = AA.operator_idx[0]
        d = AA.operator_idx[1]
        spin_dict_list.append({b: '-', a : '+', c:'+', d:'-'})

    if block == 2:

        p = AA.operator_idx[0]
        q = AA.operator_idx[1]
        k = AA.operator_idx[2]
        l = AA.operator_idx[3]

 #       spin_dict_list.append({p:'+', q:'+', k:'+', l:'+'})#0,AA
        spin_dict_list.append({p:'+', q:'-', k:'+', l:'-'})#0,CC
#        spin_dict_list.append({p:'-', q:'-', k:'-', l:'-'})#0,BB 
        spin_dict_list.append({p:'-', q:'+', k:'-', l:'+'})#0,DD



    if block == 20:
        i = BB.operator_idx[0]
        j = BB.operator_idx[1]
        r = BB.operator_idx[2]
        s = BB.operator_idx[3]

        #spin_dict_list.append({i: '+', j : '+', r:'+', s:'+'})#AA,A
        spin_dict_list.append({i: '+', j : '-', r:'+', s:'-'})#CC,A
        #spin_dict_list.append({i: '-', j : '-', r:'-', s:'-'})#BB,A
        spin_dict_list.append({i: '-', j : '+', r:'-', s:'+'})#DD,A

        
    if block == 11:
        i = BB.operator_idx[0]
        r = BB.operator_idx[1]
        p = AA.operator_idx[0]
        k = AA.operator_idx[1]

        spin_dict_list.append({i: '+', r : '+', p:'+', k:'+'})
        spin_dict_list.append({i: '+', r : '+', p:'-', k:'-'})

    elif block == 12:
        i = BB.operator_idx[0]
        r = BB.operator_idx[1]
        p = AA.operator_idx[0]
        q = AA.operator_idx[1]
        k = AA.operator_idx[2]
        l = AA.operator_idx[3]

        # spin_dict_list.append({i: '+', r : '+', p:'+', q:'+', k:'+', l:'+'})#A,AA
        spin_dict_list.append({i: '+', r : '+', p:'-', q:'-', k:'-', l:'-'})#A,BB 
        # spin_dict_list.append({i: '+', r : '+', p:'+', q:'-', k:'+', l:'-'})#A,CC 
        # spin_dict_list.append({i: '+', r : '+', p:'-', q:'+', k:'-', l:'+'})#A,DD


    elif block == 21:
        i = BB.operator_idx[0]
        j = BB.operator_idx[1]
        r = BB.operator_idx[2]
        s = BB.operator_idx[3]
        p = AA.operator_idx[0]
        k = AA.operator_idx[1]

        # spin_dict_list.append({i: '+', j : '+', r:'+', s:'+', p:'+', k:'+'})#AA,A
        spin_dict_list.append({i: '+', j : '-', r:'+', s:'-', p:'+', k:'+'})#CC,A
        spin_dict_list.append({i: '-', j : '+', r:'-', s:'+', p:'+', k:'+'})#DD,A
        #spin_dict_list.append({i: '-', j : '-', r:'-', s:'-', p:'+', k:'+'})#BB,A

    elif block == 22:
        
        j = BB.operator_idx[0]
        i = BB.operator_idx[1]
        b = BB.operator_idx[2]
        a = BB.operator_idx[3]
        
        c = AA.operator_idx[0]
        d = AA.operator_idx[1]
        k = AA.operator_idx[2]
        l = AA.operator_idx[3]
    
        # spin_dict_list.append({j: '+', i : '+', b:'+', a:'+', c:'+', d:'+', k:'+', l:'+'}) #AA_AA
        # spin_dict_list.append({j: '+', i : '+', b:'+', a:'+', c:'+', d:'-', k:'+', l:'-'}) #AA_CC
        # spin_dict_list.append({j: '+', i : '+', b:'+', a:'+', c:'-', d:'+', k:'-', l:'+'}) #AA_DD
        # spin_dict_list.append({j: '+', i : '+', b:'+', a:'+', c:'-', d:'-', k:'-', l:'-'}) #AA_BB

        # spin_dict_list.append({j: '+', i : '-', b:'+', a:'-', c:'+', d:'+', k:'+', l:'+'}) #CC_AA
        spin_dict_list.append({j: '+', i : '-', b:'+', a:'-', c:'+', d:'-', k:'+', l:'-'}) #CC_CC
        spin_dict_list.append({j: '+', i : '-', b:'+', a:'-', c:'-', d:'+', k:'-', l:'+'}) #CC_DD 
        # spin_dict_list.append({j: '+', i : '-', b:'+', a:'-', c:'-', d:'-', k:'-', l:'-'}) #CC_BB

        


        
        # spin_dict_list.append({j: '-', i : '+', b:'-', a:'+', c:'+', d:'+', k:'+', l:'+'}) #DD_AA
        # spin_dict_list.append({j: '-', i : '+', b:'-', a:'+', c:'+', d:'-', k:'+', l:'-'}) #DD_CC
        # spin_dict_list.append({j: '-', i : '+', b:'-', a:'+', c:'-', d:'+', k:'-', l:'+'}) #DD_DD
        # spin_dict_list.append({j: '-', i : '+', b:'-', a:'+', c:'-', d:'-', k:'-', l:'-'}) #DD_BB

        # spin_dict_list.append({j: '-', i : '-', b:'-', a:'-', c:'+', d:'+', k:'+', l:'+'}) #BB_AA
        # spin_dict_list.append({j: '-', i : '-', b:'-', a:'-', c:'+', d:'-', k:'+', l:'-'}) #BB_CC
        # spin_dict_list.append({j: '-', i : '-', b:'-', a:'-', c:'-', d:'+', k:'-', l:'+'}) #BB_DD
        # spin_dict_list.append({j: '-', i : '-', b:'-', a:'-', c:'-', d:'-', k:'-', l:'-'}) #BB_BB


        # spin_dict_list.append({j: '+', i : '+', b:'+', a:'+', c:'+', d:'+', k:'+', l:'+'})
        # spin_dict_list.append({j: '+', i : '+', b:'+', a:'+', c:'-', d:'-', k:'-', l:'-'})
        # spin_dict_list.append({j: '+', i : '+', b:'+', a:'+', c:'+', d:'+', k:'-', l:'-'})
        # spin_dict_list.append({j: '+', i : '+', b:'+', a:'+', c:'-', d:'-', k:'+', l:'+'})

    
    # spin_dict1['i'] = '+'
    # spin_dict1['r'] = '+'
    # spin_dict1['j'] = '+'
    # spin_dict1['s'] = '+'
    # spin_dict1['p'] = '-'
    # spin_dict1['k'] = '-'
    # spin_dict1['q'] = '-'    
    # spin_dict1['l'] = '-'
    start = time.time()
    reswielki = arithmetic_string()
    if len(spin_dict_list) > 0:
        for i in range(0, len(spin_dict_list)):
        
            numf1 = 1.0
            start = time.time()
            #        res3_abab = add_spin_driver(spin_dict1, numf1, res)
            res3_abab = add_spin_driver(spin_dict_list[i], numf1, deepcopy(res))
            for x in res3_abab:
                reswielki.append(x)
    else:
        reswielki = add_spin_driver({'pluszon':1}, 1.0, deepcopy(res))
        # reswielki = res
        print('jalala')
        for x in reswielki:
            print(x)
                
    end = time.time()
    print('spin time', end - start, 'lenn reswielki', len(reswielki))
    res3_abab = reswielki

    res = clean_cas(res3_abab)
    print('jalala2')
    for x in res:
        print(x)

    start = time.time()
    resugg = cas_to_ugg(res)    
    rsimp = simplify(resugg, cas=True)

    end = time.time()
    print('simpl2 time', end - start)
    res = ugg_to_cas(rsimp)
    print('huan1')
    for x in res:
        print(x)
    
    res = gamma_to_n(res, True)

    print('huan2')
    for x in res:
        print(x)

    
    res = clean_cas(res)
    print('huan3')
    for x in res:
        print(x)


    resugg = cas_to_ugg(res)
    rsimp = simplify(resugg, cas=True)
    res = ugg_to_cas(rsimp)
    print('huan1')
    for x in res:
        print(x)



    print('')
    print('wynik', len(res))
    print('')
    for x in res:
        print(x)

    res = dirac_to_coulomb(res)
    resugg = cas_to_ugg(res)
    rsimp = simplify(resugg, cas=True)
    res_coulomb = ugg_to_cas(rsimp)

    print('res_coulomb')
    for x in res_coulomb:
        print(x)
    print('')
    
    sys.exit(0)
#---------------------------------------------------------------------------------------
#    sys.exit(0)

    #AA +-+-  +-+-
    #AB +-+-  +--+
    #AC +-+-  -++-
    #AD +-+-  -+-+

    #BA +--+  +-+-
    #BB +--+  +--+
    #BC +--+  -++-
    #BD +--+  -+-+

    spin_dict1={}
    spin_dict1['i'] = '+'
    spin_dict1['j'] = '-'
    spin_dict1['a'] = '+'
    spin_dict1['b'] = '-'
    spin_dict1['c'] = '+'
    spin_dict1['d'] = '-'
    spin_dict1['k'] = '+'
    spin_dict1['l'] = '-'

    numf1 = 1.0
    res3_abab = add_spin_driver(spin_dict1, numf1, res)

    res = approximate_rdm2_by_rdm1_when_occ(res3_abab)
    res3_abab = deepcopy(res)
    # print('wynik')
    # for x in res3_abab:
    #     print(x)
    # sys.exit(0)

    # print('')
    # sys.exit(0)

    # res4 = cas_to_ugg(res)
    # rsimp = simplify(res3_abab)

    # res2 = cas_to_ugg(res)    
    # rsimp = simplify(res2, cas=True)

    res = clean_cas(res3_abab)
    resugg = cas_to_ugg(res)
    rsimp = simplify(resugg)
    res = ugg_to_cas(rsimp)
    
    res = gamma_to_n(res)
    
    res = clean_cas(res)

    resugg = cas_to_ugg(res)
    rsimp = simplify(resugg)
    res = ugg_to_cas(rsimp)



    print('')
    print('wynik', len(res))
    print('')
    for x in res:
        print(x)

        
    sys.exit(0)


def pphh_rpa_A(filename, op, read):

    AA = cas()
    AA.operator_idx = ['p', 'q', 'k', 'l']
    AA.operator_type = ['+', '+', '0', '0']

    BB = cas()
    BB.operator_idx = ['j', 'i', 's', 'r']
    BB.operator_type = ['+', '+', '0', '0']


    if read:
        
        f = open('.pickle_A.pkl', 'rb')
        res = pickle.load(f)
        f.close()

#         res = []

#         f = open(filename, "r")
#         lines = f.readlines()

#         for line in lines:
#             x = read_cas_from_str(line, False)
#             res.append(x)
#             old_stdout = sys.stdout
#             new_stdout = io.StringIO()
#             sys.stdout = new_stdout
#             print(x)
#             output = new_stdout.getvalue()
#             sys.stdout = old_stdout
#             if output != line:
#                 print('rozne')
#                 print('l', line)
#                 print('o', output)
#                 print('')
# #            print('')

#         print('wynik')
#         for x in res:
#             print(x)
#         print('')

    else:

        if op == 2:
            start = time.time()
            KL = evaluate(BB, AA)
            end = time.time()
            print('evaluate time', end - start)
        elif op == 1:
            start = time.time()
            KL = evaluate(BB, h1cas, AA).scale(0.5) +  evaluate(BB, h2cas, AA).scale(0.5) +\
                 evaluate(h1cas, AA, BB).scale(-0.5) +  evaluate(h2cas, AA, BB).scale(-0.5)
            end = time.time()
            print('evaluate time', end - start)

        # KL = evaluate(BB, h1cas, AA) +  evaluate(BB, h2cas, AA) 
        # KL = evaluate(h1cas, AA, BB).scale(-1.0) +  evaluate(h2cas, AA, BB).scale(-1.0)


        #   KL = KL[20:50]
        k = 0
        print('lenkl', len(KL))
        for x in KL:
            print(k, x)
            k+=1
        start = time.time()
        res = int_and_simp_cas_onedet(KL)
        end = time.time()
        print('int and simp', end - start)
        for x in res:
            print(x)
        start = time.time()
    
        #    res = zero_dens_matrix_onedet(res)
        end = time.time()
        print('zero_dens', end - start)

        res = clean_cas(res)

        start = time.time()
        print('Result po RENAME DENSITY:')    
        for x in res:
            x.rename_as_density()
        print('')
        end = time.time()
        print('rename density time', end - start)

        start = time.time()
        res2 = cas_to_ugg(res)
        end = time.time()
        print('cas to ug', end - start)

        start = time.time()
        rsimp = simplify(res2, cas=True)
        end = time.time()
        print('simplify', end - start)

        
        print('')
        print('wynik')
        print('')
        for x in rsimp:
            print(x)

        res = ugg_to_cas(rsimp)
        res = approximate_rdm3(res, True)
        res = approximate_rdm4(res, True)
        nospin = True
        # res = cum2_to_rdm2(res, nospin, False)

        # res = approximate_rdm2_by_rdm1_when_occ(res, True, True)
        print('hupla1', len(res))
        for x in res:
            print(x)

        print('')

        res2 = cas_to_ugg(res)
        rsimp = simplify(res2, cas=True)
        res = ugg_to_cas(rsimp)

        print('hupla2', len(res))
        for x in res:
            print(x)

        print('')

        if op == 1:
            f = open('.pickle_A.pkl', 'wb')
            pickle.dump(res, f)
            f.close()



    # delta = []
    # for x in range(0, len(res)):
    #     for i in range(0, len(res[x].delta)):
    #         if res[x].delta[i] not in delta:
    #                 delta.append(res[x].delta[i])

    # for x in delta:
    #     print(x)

    # sys.exit(0)

    spin_dict_list=[]

    j = BB.operator_idx[0]
    i = BB.operator_idx[1]
    b = BB.operator_idx[2]
    a = BB.operator_idx[3]
    
    c = AA.operator_idx[0]
    d = AA.operator_idx[1]
    k = AA.operator_idx[2]
    l = AA.operator_idx[3]
    
    spin_dict_list.append({j: '+', i : '+', b:'+', a:'+', c:'+', d:'+', k:'+', l:'+'}) #AA_AA
    # spin_dict_list.append({j: '+', i : '+', b:'+', a:'+', c:'+', d:'-', k:'+', l:'-'}) #AA_CC
    # spin_dict_list.append({j: '+', i : '+', b:'+', a:'+', c:'-', d:'+', k:'-', l:'+'}) #AA_DD
    
    # spin_dict_list.append({j: '+', i : '-', b:'+', a:'-', c:'+', d:'+', k:'+', l:'+'}) #CC_AA
    # spin_dict_list.append({j: '+', i : '-', b:'+', a:'-', c:'+', d:'-', k:'+', l:'-'}) #CC_CC
    # spin_dict_list.append({j: '+', i : '-', b:'+', a:'-', c:'-', d:'+', k:'-', l:'+'}) #CC_DD 
    # spin_dict_list.append({j: '+', i : '-', b:'+', a:'-', c:'-', d:'-', k:'-', l:'-'}) #CC_BB
    
    
    # spin_dict_list.append({j: '+', i : '+', b:'+', a:'+', c:'-', d:'-', k:'-', l:'-'}) #AA_BB
    
    # spin_dict_list.append({j: '-', i : '+', b:'-', a:'+', c:'+', d:'+', k:'+', l:'+'}) #DD_AA
    # spin_dict_list.append({j: '-', i : '+', b:'-', a:'+', c:'+', d:'-', k:'+', l:'-'}) #DD_CC
    # spin_dict_list.append({j: '-', i : '+', b:'-', a:'+', c:'-', d:'+', k:'-', l:'+'}) #DD_DD
    # spin_dict_list.append({j: '-', i : '+', b:'-', a:'+', c:'-', d:'-', k:'-', l:'-'}) #DD_BB
    
    # spin_dict_list.append({j: '-', i : '-', b:'-', a:'-', c:'+', d:'+', k:'+', l:'+'}) #BB_AA
    # spin_dict_list.append({j: '-', i : '-', b:'-', a:'-', c:'+', d:'-', k:'+', l:'-'}) #BB_CC
    # spin_dict_list.append({j: '-', i : '-', b:'-', a:'-', c:'-', d:'+', k:'-', l:'+'}) #BB_DD
    # spin_dict_list.append({j: '-', i : '-', b:'-', a:'-', c:'-', d:'-', k:'-', l:'-'}) #BB_BB

    reswielki = arithmetic_string()
    for i in range(0, len(spin_dict_list)):
        
        numf1 = 1.0
        start = time.time()
        #        res3_abab = add_spin_driver(spin_dict1, numf1, res)
        res3_abab = add_spin_driver(spin_dict_list[i], numf1, deepcopy(res))
        for x in res3_abab:
            reswielki.append(x)
                
    end = time.time()
    print('spin time', end - start, 'lenn reswielki', len(reswielki))
    res3_abab = reswielki

    res = clean_cas(res3_abab)
    start = time.time()
    resugg = cas_to_ugg(res)    
    rsimp = simplify(resugg, cas=True)

    end = time.time()
    print('simpl2 time', end - start)
    res = ugg_to_cas(rsimp)
    print('huan1')
    for x in res:
        print(x)
    
    res = gamma_to_n(res, True)

    print('huan2')
    for x in res:
        print(x)

    
    res = clean_cas(res)
    print('huan3')
    for x in res:
        print(x)


    resugg = cas_to_ugg(res)
    rsimp = simplify(resugg, cas=True)
    res = ugg_to_cas(rsimp)
    print('huan1')
    for x in res:
        print(x)



    print('')
    print('wynik', len(res))
    print('')
    for x in res:
        print(x)

    res = dirac_to_coulomb(res)
    resugg = cas_to_ugg(res)
    rsimp = simplify(resugg, cas=True)
    res_coulomb = ugg_to_cas(rsimp)

    print('res_coulomb', len(res_coulomb))
    for x in res_coulomb:
        print(x)
    print('')

    res = arithmetic_string()
    deltaps = ['p', 's']
    deltaqr = ['q', 'r']
    deltail = ['i', 'l']
    deltajk = ['j', 'k']

    for x in range(0, len(res_coulomb)):
        dl = res_coulomb[x].delta
        if ((deltaps in dl) and (deltaqr in dl)):
            res_coulomb[x].num_factor = 0.0

        if ((deltail in dl)) and ((deltajk in dl)):
            res_coulomb[x].num_factor = 0.0

    res = clean_cas(res_coulomb)
    print('wynik bez niedozwolonych delt', len(res))
    for x in res:
        print(x)

    
    sys.exit(0)

def test_idx_for_pphh():

    inlist = [1,2]
    aclist = [3,4,5]
    virtlist = [6,7]
    sel = [3,4]

    Ipp = []
    for p in range(1, 7):
        for q in range(1, p):
            if p>q:
                Ipp.append([p,q])

    Is = []
    for k in range(1, 7):
        for l in range(1, k):
            if (k in sel and l in sel):
                if k>l:
                    Is.append([k,l])

    Iaux = []
    for x in Ipp:
        p = x[0]
        q = x[1]
        for y in Is:
            k = y[0]
            l = y[1]

            if p!=k and q!=l:
                Iaux.append(x+y)
            if p==k and q==l:
                Iaux.append(x+y)


    for ii in range(0, len(Iaux)):
        x = Iaux[ii]
        p=x[0]
        q=x[1]
        k=x[2]
        l=x[3]
        for jj in range(0, len(Iaux)):
            y = Iaux[jj]
            r=y[0]
            s=y[1]
            i=y[2]
            j=y[3]
            
            if( (i==k) and (j==l) and (p==r) and (q==s)):
                print(r, s, i, j, p, q, k, l, ii, jj, -1)
            if( (i==k) and (j==l) and (p==r) and (q==s)):
                print(r, s, i, j, p, q, k, l, ii, jj, -1)
            if( (i==r) and (j==k) and (l==p) and (q==s)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==r) and (j==k) and (l==q) and (p==s)):
                print(r, s, i, j, p, q, k, l, ii, jj, -1)
            if( (i==s) and (j==k) and (l==p) and (q==r)):
                print(r, s, i, j, p, q, k, l, ii, jj, -1)
            if( (i==s) and (j==k) and (l==q) and (p==r)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==k) and (j==r) and (l==p) and (q==s)):
                print(r, s, i, j, p, q, k, l, ii, jj, -1)
            if( (i==k) and (j==r) and (l==q) and (p==s)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==k) and (j==s) and (l==p) and (q==r)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==k) and (j==s) and (l==q) and (p==r)):
                print(r, s, i, j, p, q, k, l, ii, jj, -1)
            if( (i==r) and (j==l) and (k==p) and (q==s)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==r) and (j==l) and (k==q) and (p==s)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==s) and (j==l) and (k==p) and (q==r)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==s) and (j==l) and (k==q) and (p==r)):
                print(r, s, i, j, p, q, k, l, ii, jj, -1)
            if( (i==l) and (j==r) and (k==p) and (q==s)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==l) and (j==r) and (k==q) and (p==s)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==l) and (j==s) and (k==p) and (q==r)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==l) and (j==s) and (k==q) and (p==r)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==k) and (j==r) and (l==q) and (p==s)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==l) and (j==r) and (k==q) and (p==s)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==r) and (j==k) and (l==q) and (p==s)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==r) and (j==l) and (k==q) and (p==s)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==k) and (j==r) and (l==p) and (q==s)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==l) and (j==r) and (k==p) and (q==s)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==r) and (j==k) and (l==p) and (q==s)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==r) and (j==l) and (k==p) and (q==s)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==k) and (j==s) and (l==q) and (p==r)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==l) and (j==s) and (k==q) and (p==r)):
                print(r, s, i, j, p, q, k, l, ii, jj, -1)
            if( (i==s) and (j==k) and (l==q) and (p==r)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==s) and (j==l) and (k==q) and (p==r)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==k) and (j==s) and (l==p) and (q==r)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==l) and (j==s) and (k==p) and (q==r)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==s) and (j==k) and (l==p) and (q==r)):
                print(r, s, i, j, p, q, k, l, ii, jj)
            if( (i==s) and (j==l) and (k==p) and (q==r)):
                print(r, s, i, j, p, q, k, l, ii, jj, -1)
            if( (p==r) and (q==s) and (i==k) and (j==l)):
                print(r, s, i, j, p, q, k, l, ii, jj, 1)                
            if( (p==r) and (q==s)):
                if (ii!=jj):
                    print(r, s, i, j, p, q, k, l, ii, jj, 'aa' )                    
            if( (i==k) and (j==l)):
                if (ii!=jj):
                    print(r, s, i, j, p, q, k, l, ii, jj, 'bb' )


    sys.exit(0)
              
                
        


    
    
def pphh_rpa_A_onedet():

    # AA = cas()
    # AA.operator_idx = ['c', 'k', 'd', 'l']
    # AA.operator_type = ['+', '0', '+', '0']

    # BB = cas()
    # BB.operator_idx = ['i', 'a', 'j', 'b']
    # BB.operator_type = ['+', '0', '+', '0']

    # BB = cas()
    # BB.operator_idx = ['k', 'l', 'c', 'd']
    # BB.operator_type = ['+', '+', '0', '0']

    # AA = cas()
    # AA.operator_idx = ['c', 'd']
    # AA.operator_type = ['+', '+']

    # AA = cas()
    # AA.operator_idx = ['a', 'i']
    # AA.operator_type = ['+', '0']

    # BB = cas()
    # BB.operator_idx = ['b', 'a']
    # BB.operator_type = ['0', '0']

    # BB = cas()
    # BB.operator_idx = ['i', 'a']
    # BB.operator_type = ['+', '0']

    # BB H AA
    # r1 = h1cas.fromleft(BB)
    # KL = [r1.fromright(AA)]
    # for x in KL:
    #     print(x)

    AA = cas()
    AA.operator_idx = ['c', 'd', 'k', 'l']
    AA.operator_type = ['+', '+', '0', '0']

    BB = cas()
    BB.operator_idx = ['j', 'i', 'b', 'a']
    BB.operator_type = ['+', '+', '0', '0']


    start = time.time()
    KL = evaluate(BB, h1cas, AA).scale(0.5) +  evaluate(BB, h2cas, AA).scale(0.5) +\
        evaluate(h1cas, AA, BB).scale(-0.5) +  evaluate(h2cas, AA, BB).scale(-0.5)
    end = time.time()
    print('evaluate time', end - start)

    # KL = evaluate(BB, h1cas, AA) +  evaluate(BB, h2cas, AA) 
    # KL = evaluate(h1cas, AA, BB).scale(-1.0) +  evaluate(h2cas, AA, BB).scale(-1.0)


#   KL = KL[20:50]
    k = 0
    print('lenkl', len(KL))
    for x in KL:
        print(k, x)
        k+=1
    start = time.time()
    res = int_and_simp_cas_onedet(KL)
    end = time.time()
    print('int and simp', end - start)
    for x in res:
        print(x)
    start = time.time()
    
    res = zero_dens_matrix_onedet(res)
    end = time.time()
    print('zero_dens', end - start)


    res = clean_cas(res)

    start = time.time()
    print('Result po RENAME DENSITY:')    
    for x in res:
        x.rename_as_density()
    print('')
    end = time.time()
    print('rename density time', end - start)

    start = time.time()
    res2 = cas_to_ugg(res)
    end = time.time()
    print('cas to ug', end - start)

    start = time.time()
    rsimp = simplify(res2, cas=True)
    end = time.time()
    print('simplify', end - start)


    # f = open('.pickle_A_onedet.pkl', 'wb')
    # pickle.dump(rsimp, f)
    # sys.exit(0)

    # f = open('.pickle_A_onedet.pkl', 'rb')
    # rsimp = pickle.load(f)
    # f.close()

    print('')
    print('wynik')
    print('')
    for x in rsimp:
        print(x)

    res = ugg_to_cas(rsimp)
    res = approximate_rdm3(res, True)
    res = approximate_rdm4(res, True)
    # res = approximate_rdm34_by_rdm1(res, False)

    res = approximate_rdm2_by_rdm1_when_occ(res, True, True)
    print('hupla1', len(res))
    for x in res:
        print(x)

    print('')

    res2 = cas_to_ugg(res)
    rsimp = simplify(res2, cas=True)
    res = ugg_to_cas(rsimp)

    print('hupla2', len(res))
    for x in res:
        print(x)

    print('')

    spin_dict_list=[]

    j = BB.operator_idx[0]
    i = BB.operator_idx[1]
    b = BB.operator_idx[2]
    a = BB.operator_idx[3]
    
    c = AA.operator_idx[0]
    d = AA.operator_idx[1]
    k = AA.operator_idx[2]
    l = AA.operator_idx[3]
    
    spin_dict_list.append({j: '+', i : '+', b:'+', a:'+', c:'+', d:'+', k:'+', l:'+'}) #AA_AA
    # spin_dict_list.append({j: '+', i : '+', b:'+', a:'+', c:'+', d:'-', k:'+', l:'-'}) #AA_CC
    # spin_dict_list.append({j: '+', i : '+', b:'+', a:'+', c:'-', d:'+', k:'-', l:'+'}) #AA_DD
    
    # spin_dict_list.append({j: '+', i : '-', b:'+', a:'-', c:'+', d:'+', k:'+', l:'+'}) #CC_AA
    # spin_dict_list.append({j: '+', i : '-', b:'+', a:'-', c:'+', d:'-', k:'+', l:'-'}) #CC_CC
    # spin_dict_list.append({j: '+', i : '-', b:'+', a:'-', c:'-', d:'+', k:'-', l:'+'}) #CC_DD 
    # spin_dict_list.append({j: '+', i : '-', b:'+', a:'-', c:'-', d:'-', k:'-', l:'-'}) #CC_BB
    
    
    # spin_dict_list.append({j: '+', i : '+', b:'+', a:'+', c:'-', d:'-', k:'-', l:'-'}) #AA_BB
    
    # spin_dict_list.append({j: '-', i : '+', b:'-', a:'+', c:'+', d:'+', k:'+', l:'+'}) #DD_AA
    # spin_dict_list.append({j: '-', i : '+', b:'-', a:'+', c:'+', d:'-', k:'+', l:'-'}) #DD_CC
    # spin_dict_list.append({j: '-', i : '+', b:'-', a:'+', c:'-', d:'+', k:'-', l:'+'}) #DD_DD
    # spin_dict_list.append({j: '-', i : '+', b:'-', a:'+', c:'-', d:'-', k:'-', l:'-'}) #DD_BB
    
    # spin_dict_list.append({j: '-', i : '-', b:'-', a:'-', c:'+', d:'+', k:'+', l:'+'}) #BB_AA
    # spin_dict_list.append({j: '-', i : '-', b:'-', a:'-', c:'+', d:'-', k:'+', l:'-'}) #BB_CC
    # spin_dict_list.append({j: '-', i : '-', b:'-', a:'-', c:'-', d:'+', k:'-', l:'+'}) #BB_DD
    # spin_dict_list.append({j: '-', i : '-', b:'-', a:'-', c:'-', d:'-', k:'-', l:'-'}) #BB_BB

    reswielki = arithmetic_string()
    for i in range(0, len(spin_dict_list)):
        
        numf1 = 1.0
        start = time.time()
        #        res3_abab = add_spin_driver(spin_dict1, numf1, res)
        res3_abab = add_spin_driver(spin_dict_list[i], numf1, deepcopy(res))
        for x in res3_abab:
            reswielki.append(x)
                
    end = time.time()
    print('spin time', end - start, 'lenn reswielki', len(reswielki))
    res3_abab = reswielki

    res = clean_cas(res3_abab)
    start = time.time()
    resugg = cas_to_ugg(res)    
    rsimp = simplify(resugg, cas=True)

    end = time.time()
    print('simpl2 time', end - start)
    res = ugg_to_cas(rsimp)
    print('huan1')
    for x in res:
        print(x)
    
    res = gamma_to_n(res, True)

    print('huan2')
    for x in res:
        print(x)

    
    res = clean_cas(res)
    print('huan3')
    for x in res:
        print(x)


    resugg = cas_to_ugg(res)
    rsimp = simplify(resugg, cas=True)
    res = ugg_to_cas(rsimp)
    print('huan1')
    for x in res:
        print(x)



    print('')
    print('wynik', len(res))
    print('')
    for x in res:
        print(x)

    res = dirac_to_coulomb(res)
    resugg = cas_to_ugg(res)
    rsimp = simplify(resugg, cas=True)
    res_coulomb = ugg_to_cas(rsimp)

    print('res_coulomb')
    for x in res_coulomb:
        print(x)
    print('')
    
    sys.exit(0)













#-----------------------------------------------------------------------------------------------------------------
    

#     res = approximate_rdm3_by_rdm12(res)
#     res = approximate_rdm34_by_rdm1(res, True)

# #    print('troc', len(res))    
#     res = approximate_rdm2_by_rdm1_when_occ(res, True, True)
# #    print('troc', len(res))
#     end = time.time()
#     print('apprix time', end - start)

    print('hupla1', len(res))
    k = 1
    for x in res:
        print(k, x)
        k+=1
    print('')

    # res = gamma_to_n(res)
    
    # res = clean_cas(res)
    # start = time.time()

    resugg = cas_to_ugg(res)
    rsimp = simplify(resugg, cas=True)
    end = time.time()
    print('simplifyprzed time', end - start)
    
    res = ugg_to_cas(rsimp)



    print('')
    print('hupla2', len(res))
    print('')
    for x in res:
        print(x)
    print('')


    # f2 = open('.pickle_A_onedet2.pkl', 'wb')
    # pickle.dump(res, f2)
    # f2.close()
    # sys.exit(0)

    # f2 = open('.pickle_A_onedet2.pkl', 'rb')
    # res = pickle.load(f2)
    # f2.close()


    #AA +-+-  +-+-
    #AB +-+-  +--+
    #AC +-+-  -++-
    #AD +-+-  -+-+

    #BA +--+  +-+-
    #BB +--+  +--+
    #BC +--+  -++-
    #BD +--+  -+-+

    # spin_dict1={}
    # spin_dict1['i'] = '+'
    # spin_dict1['j'] = '+'
    # spin_dict1['a'] = '-'
    # spin_dict1['b'] = '-'
    # spin_dict1['c'] = '+'
    # spin_dict1['d'] = '+'
    # spin_dict1['k'] = '-'
    # spin_dict1['l'] = '-'

    #AAAA
    spin_dict_list = []

    # AA.operator_idx = ['p', 'k', 'q', 'l']
    # BB.operator_idx = ['i', 'r', 'j', 's']
    i = BB.operator_idx[0]
    r = BB.operator_idx[1]
    j = BB.operator_idx[2]
    s = BB.operator_idx[3]

    p = AA.operator_idx[0]
    k = AA.operator_idx[1]
    q = AA.operator_idx[2]
    l = AA.operator_idx[3]

    
    spin_dict_list.append({i: '+', r : '+', j:'+', s:'+', p:'+', k:'+', q:'+', l:'+'})
    spin_dict_list.append({i: '+', r : '+', j:'+', s:'+', p:'+', k:'+', q:'-', l:'-'})
    spin_dict_list.append({i: '+', r : '+', j:'+', s:'+', p:'-', k:'-', q:'+', l:'+'})
    spin_dict_list.append({i: '+', r : '+', j:'+', s:'+', p:'-', k:'-', q:'-', l:'-'})
    spin_dict_list.append({i: '+', r : '+', j:'-', s:'-', p:'+', k:'+', q:'+', l:'+'})
    spin_dict_list.append({i: '+', r : '+', j:'-', s:'-', p:'+', k:'+', q:'-', l:'-'})
    spin_dict_list.append({i: '+', r : '+', j:'-', s:'-', p:'-', k:'-', q:'+', l:'+'})
    spin_dict_list.append({i: '+', r : '+', j:'-', s:'-', p:'-', k:'-', q:'-', l:'-'})
    
    # spin_dict1={}
    # spin_dict1['i'] = '+'
    # spin_dict1['r'] = '+'
    # spin_dict1['j'] = '+'
    # spin_dict1['s'] = '+'
    # spin_dict1['p'] = '-'
    # spin_dict1['k'] = '-'
    # spin_dict1['q'] = '-'
    # spin_dict1['l'] = '-'
    reswielki = arithmetic_string()
    for i in range(0, len(spin_dict_list)):
        
        numf1 = 1.0
        start = time.time()
        #        res3_abab = add_spin_driver(spin_dict1, numf1, res)
        res3_abab = add_spin_driver(spin_dict_list[i], numf1, deepcopy(res))
        for x in res3_abab:
            reswielki.append(x)
                
    end = time.time()
    print('spin time', end - start, 'lenn reswielki', len(reswielki))
    res3_abab = reswielki



    

    # print('przedclean', len(res3_abab))
    # for k in range(0, len(res3_abab)):
    #     print(k+1, res3_abab[k])
    # print('')


    # res = approximate_rdm2_by_rdm1_when_occ(res3_abab)
    #res3_abab = deepcopy(res)
    # print('przedclean', len(res3_abab))
    # for k in range(0, len(res3_abab)):
    #     print(k, res3_abab[k])
    # print('')

    # print('wynik')
    # for x in res3_abab:
    #     print(x)
    # sys.exit(0)

    # print('')
    # sys.exit(0)

    # res4 = cas_to_ugg(res)
    # rsimp = simplify(res3_abab)

    # res2 = cas_to_ugg(res)    
    # rsimp = simplify(res2, cas=True)

    res = clean_cas(res3_abab)

    # print('po res3_abab', len(res))
    # for k in range(0, len(res)):
    #     print(k+1, res[k])
    # print('')

    start = time.time()
    resugg = cas_to_ugg(res)    
    rsimp = simplify(resugg, cas=True)

    end = time.time()
    print('simpl2 time', end - start)


    res = ugg_to_cas(rsimp)
    # print('po rsimpx')
    # for k in range(0, len(res)):
    #     print(k+1, res[k])
    # print('')
    
    res = gamma_to_n(res, True)
    
    res = clean_cas(res)

    resugg = cas_to_ugg(res)
    rsimp = simplify(resugg, cas=True)
    res = ugg_to_cas(rsimp)



    print('')
    print('wynik', len(res))
    print('')
    for x in res:
        print(x)

        
    sys.exit(0)

def gamma_to_n(res, spin=False):

    for i in range(0, len(res)):
#        print(i, res[i])
        nnl = []
        for j in range(0, len(res[i].coefficient)):
#            print('teraz', j, res[i].coefficient[j])
            if res[i].coefficient[j] == DENS1:
                if res[i].num_factor != 0.0:
                    res[i].new_delta(res[i].coefficient_idx[j][0], res[i].coefficient_idx[j][1])
                    if res[i].coefficient_idx[j][0] in occupied and res[i].coefficient_idx[j][1] in general:
#                        print('zz1', res[i].coefficient[j], res[i].coefficient_idx[j])
                        nnl.append(j)
                        ind0 = deepcopy([res[i].coefficient_idx[j][0]])
                    elif res[i].coefficient_idx[j][1] in occupied and res[i].coefficient_idx[j][0] in general:
 #                       print('zz2', res[i].coefficient[j], res[i].coefficient_idx[j])
                        nnl.append(j)
                        ind0 = deepcopy([res[i].coefficient_idx[j][1]])
                    elif res[i].coefficient_idx[j][1] in occupied and res[i].coefficient_idx[j][0] in occupied:
  #                      print('zz3', res[i].coefficient[j], res[i].coefficient_idx[j])
                        nnl.append(j)                        
                        ind0 = deepcopy([res[i].coefficient_idx[j][0]])
                    elif res[i].coefficient_idx[j][1] in general and res[i].coefficient_idx[j][0] in general:
                        #                     print('zz3', res[i].coefficient[j], res[i].coefficient_idx[j])
                        #                        nnl.append(j)                        
                        ind0 = deepcopy([res[i].coefficient_idx[j][0],res[i].coefficient_idx[j][1]])
                    else:
                        ind0 = deepcopy([res[i].coefficient_idx[j][0]])
                        res[i].num_factor = 0.0
                                                
                                                
                    res[i].coefficient_idx[j] = ind0 
                    res[i].coefficient[j] = 'n'
                    # if not spin:
                    #     res[i].num_factor *= 2.0
        res[i].exec_delta()
        for k in range(0, len(nnl)):
            res[i].coefficient.pop(nnl[k])
            res[i].coefficient_idx.pop(nnl[k])
            for l in range(k+1, len(nnl)):
                nnl[l] -=1
        print('')
    print('ato', len(res))
    for x in res:
        print(x)

    return res


def twoel():

    cd = cas()
    cd.operator_idx = ['c', 'd']
    cd.operator_type = ['+', '+']

    ab = cas()
    ab.operator_idx = ['a', 'b']
    ab.operator_type = ['0', '0']

    print(cd, ab)
    KL1 = evaluate(h1cas, cd, ab).scale(-1.0)
    print('')
    print('wynik')
    for  x in KL1:
        print(x)

    reskl1 = []
    for x in KL1:
        result_list = x.wick_ca()
        reskl1 = reskl1 + result_list
    print('')



    for x in reskl1:
        x.exec_delta()

    print('Result po WICK:')    

    for x in reskl1:
        print(x)

    print('Result po RENAME DENSITY:')    
    
    for x in reskl1:
        x.rename_as_density()


    res2kl1 = cas_to_ugg(reskl1)

    rsimpkl1 = simplify(res2kl1, cas=True)
    k = 0
    print('po simp')
    for x in rsimpkl1:
        k += 1
        # print("&", x, "\\\\")
        print(x)

    sys.exit(0)

def pp_rpa():
    

    # -------------------------------------------------------------------------------------------------------
    # W scislym przypadku <0|[[dO, H],O+]|0> = <0|[dO,[H,O+]]|0>
    # <0|[dO, H, O+]|0>  = 1/2(<0|[[dO, H],O+]|0> + <0|[dO,[H,O+]]|0>)
    #
    # KL = <0|[[dO, H],O+]|0> = |jesli killer| =  <0|[dO, H],O+]0>
    # z tego <0|[[rs, H], q+p+]|0>  powinno wyjść to samo co z <0|[rs, H], q+p+|0>
    # KL = <0|[[rs, H], q+p+]|0>
    # KL_killer = <0|[rs, H], q+p+|0>  <--- tego nie bedziemy stosowac, bo ten powyzej jest leszpy - zwiazany
    #
    # KP = <0|[dO,[H,O+]]|0> = |jesli killer| = <0|dO,[H,O+]]0>
    # z tego  <0|[rs, [H, q+p+]|0>   powinno wyjść to samo co z   <0|rs, [H, p+q+]|0>  (w notacji komutatory na prawo)
    # z tego -<0|[[H, q+p+], rs|0>   powinno wyjść to samo co z   <0|rs, [H, q+p+]|0>  (w notacji komutatory na lewo)
    # KP = -<0|[[H, q+p+], rs|0>
    # KP_killer = <0|rs, [H, q+p+]|0>

    # wnioski
    # 1. Z komutatorow KL i KP powinno wyjsc to samo, jesli jest spelninony killer, czyli pq|0> = 0. Nigdzie pewnie tego nie zakladamy.
    # Zbadac czy takie wyrazy powstaja przed nazwaniem macierzy gestosci.
    

    # for PP
    # rs = AArs_aa
    # p+q+ = AAqp_pp

    # for HH
    # r+s+ = AArs_cc
    # pq = AAqp_hh

    # KL = <0|[[rs, H], q+p+]|0>

    #    KL = evaluate(AArs_aa, AAqp_pp)
#    KL = evaluate(AArs_aa, AApq_pp)
    #    KL = evaluate(AAsr_aa, AApq_pp)

    # KL = evaluate(h1cas, AAqp_pp, AArs_aa).scale(-1.0) +  evaluate(h2cas, AAqp_pp, AArs_aa).scale(-1.0)

    # KL = evaluate(AArs_aa, h1cas, AAqp_pp) +  evaluate(AArs_aa, h2cas, AAqp_pp)
    
    # KL = evaluate(AArs_aa, h1cas, AAqp_pp) +  evaluate(AArs_aa, h2cas, AAqp_pp) +\
    #     evaluate(h1cas, AAqp_pp, AArs_aa) +  evaluate(h2cas, AAqp_pp, AArs_aa)


    KL = evaluate(AArs_aa, h1cas, AAqp_pp).scale(0.5) +  evaluate(AArs_aa, h2cas, AAqp_pp).scale(0.5) +\
        evaluate(h1cas, AAqp_pp, AArs_aa).scale(-0.5) +  evaluate(h2cas, AAqp_pp, AArs_aa).scale(-0.5)

    KL1 = evaluate(AArs_aa, h1cas, AApq_pp).scale(0.5) +  evaluate(AArs_aa, h2cas, AApq_pp).scale(0.5) +\
        evaluate(h1cas, AApq_pp, AArs_aa).scale(-0.5) +  evaluate(h2cas, AApq_pp, AArs_aa).scale(-0.5)

    # KL = evaluate(AApq_pp, h1cas, AArs_aa).scale(0.5) +  evaluate(AApq_pp, h2cas, AArs_aa).scale(0.5) +\
    #     evaluate(h1cas, AArs_aa, AApq_pp).scale(-0.5) +  evaluate(h2cas, AArs_aa, AApq_pp).scale(-0.5)

#    KL1 = KL

    # KL = evaluate(AArs_aa, h1cas, AApq_pp)+  evaluate(AArs_aa, h2cas, AApq_pp)
    # KL = evaluate(h1cas, AApq_pp, AArs_aa) +  evaluate(h2cas, AApq_pp, AArs_aa).scale(-1.0)
    
    # KL = evaluate(AArs_aa, h1cas, AApq_pp).scale(0.5) +  evaluate(AArs_aa, h2cas, AApq_pp).scale(0.5) +\
    #     evaluate(h1cas, AApq_pp, AArs_aa).scale(-0.5) +  evaluate(h2cas, AApq_pp, AArs_aa).scale(-0.5)

    # KL = evaluate(AAsr_aa, h1cas, AAqp_pp).scale(0.5) +  evaluate(AAsr_aa, h2cas, AAqp_pp).scale(0.5) +\
    #     evaluate(h1cas, AAqp_pp, AAsr_aa).scale(-0.5) +  evaluate(h2cas, AAqp_pp, AAsr_aa).scale(-0.5)


    # KL = evaluate(AAsr_aa, h1cas, AApq_pp).scale(0.5) +  evaluate(AAsr_aa, h2cas, AApq_pp).scale(0.5) +\
    #     evaluate(h1cas, AApq_pp, AAsr_aa).scale(-0.5) +  evaluate(h2cas, AApq_pp, AAsr_aa).scale(-0.5)


    # KL = evaluate(AArs_cc, AAqp_hh)

    # KL = evaluate(h1cas, AAqp_hh, AArs_cc).scale(-1.0) +  evaluate(h2cas, AAqp_hh, AArs_cc).scale(-1.0)
    
    # KL = evaluate(AArs_cc, h1cas, AAqp_hh) +  evaluate(AArs_cc, h2cas, AAqp_hh)

    # KL = evaluate(AArs_cc, h1cas, AAqp_hh) +  evaluate(AArs_cc, h2cas, AAqp_hh) +\
    #     evaluate(h1cas, AAqp_hh, AArs_cc) +  evaluate(h2cas, AAqp_hh, AArs_cc)

    
    # KL = evaluate(AArs_cc, h1cas, AAqp_hh).scale(0.5) +  evaluate(AArs_cc, h2cas, AAqp_hh).scale(0.5) +\
    #     evaluate(h1cas, AAqp_hh, AArs_cc).scale(-0.5) +  evaluate(h2cas, AAqp_hh, AArs_cc).scale(-0.5)


    print('wynik komutatorow KL')

    k = 0
    for x in KL:
        print(k, x)
        k += 1
    print('')

    print('Perform Wick')
    res = []
    for x in KL:
        result_list = x.wick_ca()
        res = res + result_list
    print('')

    reskl1 = []
    for x in KL1:
        result_list = x.wick_ca()
        reskl1 = reskl1 + result_list
    print('')

    
    print('Result po WICK:')    
    for x in res:
        x.exec_delta()
        if (x.num_factor) != 0.0:
            print(x)
        print(x)
    print('')

    for x in reskl1:
        x.exec_delta()


    print('Result po RENAME DENSITY:')    
    for x in res:
        x.rename_as_density()
        if (x.num_factor) != 0.0:
            print(x)
    print('')
    
    for x in reskl1:
        x.rename_as_density()


    res2 = cas_to_ugg(res)
    res2kl1 = cas_to_ugg(reskl1)

    print('po tranformacji')
    for x in res2:
        print(x)
    
    rsimp = simplify(res2, cas=True)
    rsimpkl1 = simplify(res2kl1, cas=True)
    k = 0
    print('po simp')
    for x in rsimp:
        k += 1
        # print("&", x, "\\\\")
        print(x)

    res = cas_to_ugg(rsimp)
    reskl1 = cas_to_ugg(rsimpkl1)
    # print(
    # res = ari
    k = 0
    print('po tranformacji')
    for x in res:
        print(k, x)
        k += 1


    # sys.exit(0)

    spin_dict1={}
    spin_dict1['r'] = '+'
    spin_dict1['s'] = '-'
    spin_dict1['p'] = '+'
    spin_dict1['q'] = '-'
    numf1 = 1.0

    spin_dict2={}
    spin_dict2['r'] = '+'
    spin_dict2['s'] = '-'
    spin_dict2['p'] = '-'
    spin_dict2['q'] = '+'
    numf2 = -1.0


    spin_dict3={}
    spin_dict3['r'] = '+'
    spin_dict3['s'] = '+'
    spin_dict3['p'] = '+'
    spin_dict3['q'] = '+'

    spin_dict4={}
    spin_dict4['r'] = '-'
    spin_dict4['s'] = '-'
    spin_dict4['p'] = '-'
    spin_dict4['q'] = '-'


    res2 = deepcopy(res)
    res3_abab = add_spin_driver(spin_dict1, numf1, res)
    # rsimp = res3_abab
    res3_abba = add_spin_driver(spin_dict2, numf2, res2)
    rsimp = res3_abab + res3_abba
    # res3_abba = add_spin_driver(spin_dict2, numf1, res)
    # rsimp = res3_abba
    # res3_abab = res3_abba


    

    #-------------------------------------------------


    # res3_abab = add_spin_driver(spin_dict1, numf1, res)
    # res3_abba = add_spin_driver(spin_dict2, numf1, res)
    # rsimp = res3_abab + res3_abba
    
    
    # res3_aaaa = add_spin_driver(spin_dict3, numf1, res)
    # res3_bbbb = add_spin_driver(spin_dict4, numf1, res)
    # res3_abab = res3_aaaa + res3_bbbb

    # res3_aaaa = add_spin_driver(spin_dict3, numf1, res)
    # rsimp = res3_aaaa


    # res3_bbbb = add_spin_driver(spin_dict4, numf1, res)
    # rsimp = res3_bbbb

    print('')
    print('last result po spin abab')

  #  rsimp = res3_abab
        
    rsimp.cleanup()
    for x in rsimp:
        print(x)
#    sys.exit(0)
    # print('')
    # print('last result po spin abba')
    # res3_abba.cleanup()
    # for x in res3_abba:
    #     print(x)

    #    sys.exit(0)
    rgam_abab = simplify_for_mult_2(rsimp)
#    sys.exit(0)
    # rgam_abba = simplify_for_mult_2(res3_abba)

    # rgam = rgam_abab + rgam_abba
    # rsimp = simplify(rgam, cas=True)
    rsimp = rgam_abab
    rsimp = rename_gm1(rsimp)
    # rgam_abab = rename_gm1(rgam_abab)
    print('po simp and rename gamma')
    # for x in rgam_abab:
    #     print(x)
    for x in rsimp:
        print(x)
    sys.exit(0)


    res2 = cas_to_ugg(res)
    print('po tranformacji')
    for x in res2:
        print(x)
    
    rsimp = simplify(res2, cas=True)
    k = 0
    print('po simp')
    for x in rsimp:
        k += 1
        # print("&", x, "\\\\")
        print(x)

    rsimpa = find_antysym_int(rsimp)
    k = 0
    print('')
    print('po simp antysym')
    print('')
    for x in rsimpa:
        k += 1
        # print("&", x, "\\\\")
        print(x)

    sys.exit(0)



def process_single_block(args):

    key, Ablock_dict, rgam_original, myfixed = args
    
    Ablock_idx = key
    rgam = deepcopy(rgam_original)
    
    rgamSplit = arithmetic_string()
    for x in rgam:
        for fx in range(0, 4):
            x.orbital_type[myfixed[fx]] = Ablock_dict[Ablock_idx][fx]
            if len(x.summation) == 0:
                rgamSplit += arithmetic_string(x)
            else:
                new_elements = split_element(x)
                rgamSplit += new_elements
    
    remove_terms_with_virtual_in_any_DM(rgamSplit)
    remove_non_active_integrals_for_a0(rgamSplit)
    
    return key, rgamSplit

def process_blocks_parallel(Ablock_dict, rgam_original, myfixed=['p','q','r','s'], n_processes=None):


    args = [(key, Ablock_dict, rgam_original, myfixed) for key in Ablock_dict]
    
    with Pool(processes=n_processes) as pool:
        results = pool.map(process_single_block, args)
    
    A0_blocks = dict(results)
    
    return A0_blocks

def ph_rpa_plus_doubles_overlap():

    print(AArs)
    print(AAqp)

    #    KL = evaluate( AATqpqp, AArs)
    #KL = evaluate( AArsg, AApqpq)

    #M1
#    KL = evaluate( AArsrs, AApqpq)
    #M2
#    KL = evaluate( AArsrs, AAqpqp)
    #M3
#    KL = evaluate( AAsrss, AApqpq)
    #M4
 #   KL = evaluate( AAsrsr, AAqpqp)
        
 #KL = evaluate( AAsrsr, AAqp)
    #    KL = evaluate( AApq,  AArs)

    #KL = evaluate( AArs,  AApq)
    KL = evaluate( AArsrs,  AApq)

 #   KL = evaluate(AArsrs, AAqpqp)
   # KL = evaluate(AArsrs, AApq)

    print('kl')
    for x in KL:
        print(x)
    print()
    print('Perform Wick')
    res = []
    for x in KL:
        result_list = x.wick_ca()
        res = res + result_list
    print('Result po WICK:')    
    for x in res:
        x.exec_delta()
    for x in res:
        print(x)
    print()

    print()
    for x in res:
        print(x)

    for x in res:
        x.rename_as_density()
    print('Result po RENAME DENSITY:')
    for x in res:
        print(x)
    res2 = cas_to_ugg(res)
    print('typeres1', type(res2), type(res2[0]))
    for x in res2:
        print(x)
    rsimp = simplify(res2, cas=True)
    print('typeres2', type(rsimp), type(rsimp[0]))

    print('po simplify spinorb')
    for x in rsimp:
        print(x)
    print()


    spin_dict1={}
    spin_dict1['r'] = '+'
    spin_dict1['s'] = '+'
    spin_dict1['v'] = '-'
    spin_dict1['w'] = '-'

    spin_dict1['p'] = '+'
    spin_dict1['q'] = '+'
#    spin_dict1['t'] = '+'
#    spin_dict1['u'] = '+'

    # spin_dict1['r'] = '+'
    # spin_dict1['s'] = '+'

    # spin_dict1['p'] = '-'
    # spin_dict1['q'] = '-'


    # spin_dict1['r'] = '-'
    # spin_dict1['s'] = '-'

    # spin_dict1['p'] = '+'
    # spin_dict1['q'] = '+'

    # spin_dict1['r'] = '-'
    # spin_dict1['s'] = '-'

    # spin_dict1['p'] = '-'
    # spin_dict1['q'] = '-'

    numf1 = 1.0

        
    rsimp2 = deepcopy(rsimp)
    rgam = add_spin_driver(spin_dict1, numf1, rsimp)
    print('typeres3', type(rgam), type(rgam[0]))

        
    print('')
    print('last result po spin')
    for x in rgam:
        if x.num_factor != 0 and x.num_factor != 0.0:
            print(x)

#    rgams = simplify_for_mult_2(rgam)
    rgams = simplify_for_mult_234(rgam)
 
    

    print()
    print('po simp and rename gamma', len(rgams))
    print()
    print(type(rgams))
    for x in rgams:
        print(x)
    print()
     


def ph_rpa_plus_doubles(pick, cis=None):

    aaa = 1

    #if pick == "dump":
    if aaa == 1:

        # KL = evaluate(AArs, h1cas, AApq).scale(0.5) +  evaluate(AArs, h2cas, AApq).scale(0.5) +\
    #     evaluate(h1cas, AApq, AArs).scale(-0.5) +  evaluate(h2cas, AApq, AArs).scale(-0.5)

        # KL = evaluate( AArsrs, h1cas, AAqpqp).scale(0.5) + evaluate(h1cas, AAqpqp, AArsrs).scale(-0.5) 
        
        #E1
        # KL = evaluate( AArsrs, h1cas, AATqpqp).scale(0.5) +  evaluate( AArsrs, h2cas, AATqpqp).scale(0.5)+\
        #    evaluate(h1cas, AATqpqp, AArsrs).scale(-0.5) +  evaluate(h2cas, AATqpqp, AArsrs).scale(-0.5)
        #E2
        # KL = evaluate( AArsrs, h1cas, AApqpq).scale(0.5) +  evaluate( AArsrs, h2cas, AApqpq).scale(0.5)+\
        #    evaluate(h1cas, AApqpq, AArsrs).scale(-0.5) +  evaluate(h2cas, AApqpq, AArsrs).scale(-0.5)

        #P1
        # KL = evaluate(AArs, AATqpqp)
        #P2
        # KL = evaluate(AArs, AApqpq)
        #P3
 #       KL = evaluate(AArsg, AATqpqp)
        #P4
#        KL = evaluate(AArsg, AApqpq)

        #P5
        # KL = evaluate(AApqpq, AArsg)


#        KL = evaluate(AAsrsr, AApq)
 #       KL = evaluate(AAsrsr, AAqp)
  #      KL = evaluate(AArsrs, AApq)
#        KL = evaluate(AArsrs, AAqp)

        #M1
#        KL = evaluate(AAsrsr, AApqpq)
        #M2
        # KL = evaluate(AAsrsr, AATqpqp)
        #M3
        # KL = evaluate(AArsrs, AApqpq)
        #M4
        # KL = evaluate(AArsrs, AATqpqp)

#        KL = evaluate( AArs, AApqpq)
#        KL = evaluate( AArsrs, AApqpq)
        #KL = evaluate( AArsrs, AATqpqp)
#        KL = evaluate( AAsrsr, AApqpq)
#        KL = evaluate( AAsrsr, AATqpqp)

        #M1
        #KL = evaluate( AArsrs, AATqpqp)
        #M2
        #KL = evaluate( AArsrs, AApqpq)        
        #M3
        #KL = evaluate( AAsrsr, AATqpqp)
        #M4
        #KL = evaluate( AAsrsr, AApqpq)

        #KL = evaluate(AArs, AATqpqp)
#        KL = evaluate(AArs, AApqpq)
#        KL = evaluate(AArsg, AATqpqp)
#        KL = evaluate(AArsg, AApqpq)

        # basis
        # q*p   p*q   q*pq*p   p*qp*q
        #
        # r*s
        # s*r
        # r*sr*s
        # s*rs*r

        #C
        # KL = evaluate(h1cas, AATqpqp, AArs).scale(-1.0) +  evaluate(h2cas, AATqpqp, AArs).scale(-1.0)

        #C
#        KL = evaluate( AArs, h1cas, AATqpqp).scale(0.5) +  evaluate( AArs, h2cas, AATqpqp).scale(0.5)+\
 #          evaluate(h1cas, AATqpqp, AArs).scale(-0.5) +  evaluate(h2cas, AATqpqp, AArs).scale(-0.5)
        #C*
        # KL = evaluate( AArsg, h1cas, AApqpq).scale(0.5) +  evaluate( AArsg, h2cas, AApqpq).scale(0.5)+\
        #     evaluate(h1cas, AApqpq, AArsg).scale(-0.5) +  evaluate(h2cas, AApqpq, AArsg).scale(-0.5)

#        if cis == False:
#            KL = evaluate(h1cas, AAqp, AArs).scale(-1.0) +  evaluate(h2cas, AAqp, AArs).scale(-1.0)
        if cis == 'ph':
            #
            # this is for the CIS code NEVPTS *********************************
            #
#            ML = evaluate(h1cas, AApq, AAsr).scale(-1.0) +  evaluate(h2cas, AApq, AAsr).scale(-1.0)
#            ML = evaluate(h1cas, AApq, AArs).scale(-1.0) +  evaluate(h2cas, AApq, AArs).scale(-1.0)
#            ML =  evaluate(h2cas, AApq, AArs).scale(-1.0)
#            ML = evaluate(h1cas, AArs, AApq).scale(-1.0) +  evaluate(h2cas, AArs, AApq).scale(-1.0)
#            ML = evaluate(h1cas, AAsr, AAqp).scale(-1.0) +  evaluate(h2cas, AAsr, AAqp).scale(-1.0)
#            ML = evaluate(h1cas, AArs, AApq).scale(-1.0) +  evaluate(h2cas, AArs, AApq).scale(-1.0) B
            ML = evaluate(h1cas, AArs, AAqp).scale(-1.0) +  evaluate(h2cas, AArs, AAqp).scale(-1.0) #A1?
 #           ML1 = evaluate(h1cas, AAsr, AAqp).scale(-1.0) +  evaluate(h2cas, AAsr, AAqp).scale(-1.0)
#            ars H apq   [ars, [h, apq]]=-[[h, apq],ars]  | [[ars, h], apq]
            # ML = evaluate(h1cas, AAqp, AAsr).scale(-0.5) +evaluate(h2cas, AAqp, AAsr).scale(-0.5)   +\
            #     evaluate(AAsr, h1cas, AAqp).scale(0.5) +evaluate(AAsr, h2cas, AAqp).scale(0.5) 
            
#            ML = evaluate(h1cas, AArs, AApq).scale(-1.0) +  evaluate(h2cas, AArs, AApq).scale(-1.0)
#            ML = evaluate(h1cas, AAqp, AAsr).scale(-1.0) +  evaluate(h2cas, AAqp, AAsr).scale(-1.0)
            #ML = evaluate(h1cas, AAsr) +  evaluate(h2cas, AAsr)
#            ML = evaluate(h1cas, AAqp, AAsr).scale(-1.0) +  evaluate(h2cas, AAqp, AAsr).scale(-1.0)
#            KL = evaluate(h1cas, AApq) +  evaluate(h2cas, AApq)
            # print('kl')
            # for r in KL:
            #     print(r)
            # sys.exit(0)
#            KL = evaluate(h1cas, AAsr) +  evaluate(h2cas, AAsr)
#            KL = evaluate(h1cas, AArs) +  evaluate(h2cas, AArs)
        if cis == 'pp':
            KL = evaluate(h1cas, AAqp_pp) +  evaluate(h2cas, AAqp_pp)


        #C^T
        # KL = evaluate( AATqpqp, h1cas, AArs).scale(0.5) +  evaluate(AATqpqp, h2cas, AArs).scale(0.5)+\
        #     evaluate(h1cas, AArs, AATqpqp).scale(-0.5) +  evaluate(h2cas, AArs, AATqpqp).scale(-0.5)
        #C^T
        # KL = evaluate( AApqpq, h1cas, AArsg).scale(0.5) +  evaluate(AApqpq, h2cas, AArsg).scale(0.5)+\
        #     evaluate(h1cas, AArsg, AApqpq).scale(-0.5) +  evaluate(h2cas, AArsg, AApqpq).scale(-0.5)
        #D^T
        # KL = evaluate( AApqpq, h1cas, AArs).scale(0.5) +  evaluate(AApqpq, h2cas, AArs).scale(0.5)+\
        #     evaluate(h1cas, AArs, AApqpq).scale(-0.5) +  evaluate(h2cas, AArs, AApqpq).scale(-0.5)

        #D
        # KL = evaluate( AArs, h1cas, AApqpq).scale(0.5) +  evaluate( AArs, h2cas, AApqpq).scale(0.5)+\
        #     evaluate(h1cas, AApqpq, AArs).scale(-0.5) +  evaluate(h2cas, AApqpq, AArs).scale(-0.5)

        #D
        # KL = evaluate( AArs, h1cas, AApqpq).scale(0.5) +  evaluate( AArs, h2cas, AApqpq).scale(0.5)+\
        #     evaluate(h1cas, AApqpq, AArs).scale(-0.5) +  evaluate(h2cas, AApqpq, AArs).scale(-0.5)

        # KL = evaluate( AArsg, h1cas, AATqpqp).scale(0.5) +  evaluate( AArsg, h2cas, AATqpqp).scale(0.5)+\
        #     evaluate(h1cas, AATqpqp, AArsg).scale(-0.5) +  evaluate(h2cas, AATqpqp, AArsg).scale(-0.5)


        #A
        # KL = evaluate(AArs, h1cas, AAqp).scale(0.5) +  evaluate(AArs, h2cas, AAqp).scale(0.5) +\
        # evaluate(h1cas, AAqp, AArs).scale(-0.5) +  evaluate(h2cas, AAqp, AArs).scale(-0.5)
    
        
        # KL = evaluate(h1cas, AAqpqp, AArsrs).scale(-1.0) +  evaluate(h2cas, AAqpqp, AArsrs).scale(-1.0)

        print()
        print('------------')
        print()

        FLa = arithmetic_string(h1cas)
        FLa = FLa.fromleft(AAqp)
        FLa = FLa.fromright(AAsr)
        FL2a = arithmetic_string(h2cas)
        FL2a = FL2a.fromleft(AAqp)
        FL2a = FL2a.fromright(AAsr)


        FLb = arithmetic_string(h1cas)
        FLb = FLb.fromright(AAsr)
        FLb = FLb.fromright(AAqp)
        FL2b = arithmetic_string(h2cas)
        FL2b = FL2b.fromright(AAsr)
        FL2b = FL2b.fromright(AAqp)


        # FLc = arithmetic_string(h1cas)
        # FLc = FLc.fromleft(AAsr)
        # FLc = FLc.fromleft(AAqp)
        # FL2c = arithmetic_string(h2cas)
        # FL2c = FL2c.fromleft(AAsr)
        # FL2c = FL2c.fromleft(AAqp)


        FLc = arithmetic_string(h1cas)
        FLc = FLc.fromright(AArs)
        FLc = FLc.fromright(AApq)
        FL2c = arithmetic_string(h2cas)
        FL2c = FL2c.fromright(AArs)
        FL2c = FL2c.fromright(AApq)

        HLc = arithmetic_string(h1cas)
        HLc = HLc.fromright(AApq)
        HLc = HLc.fromright(AArs)
        HL2c = arithmetic_string(h2cas)
        HL2c = HL2c.fromright(AApq)
        HL2c = HL2c.fromright(AArs)

        


        FLd = arithmetic_string(h1cas)
        FLd = FLd.fromleft(AAsr)
        FLd = FLd.fromright(AAqp)
        FL2d = arithmetic_string(h2cas)
        FL2d = FL2d.fromleft(AAsr)
        FL2d = FL2d.fromright(AAqp)

        #        KL = FLa + FL2a #+ FLb + FL2b + FLc + FL2c + FLd + FL2d
#        KL =  FLc + FL2c
#        KL =  HLc + HL2c
#        KL = ML.fromleft(AAqp)
        # KL = arithmetic_string(h1cas)
        # KL = KL.fromright(AAps)
        # KL2 = arithmetic_string(h2cas)
        # KL2 = KL2.fromright(AAps)

        # KL = KL + KL2
        # for i in range(0, len(KL)):
        #     print(KL[i].delta)
        #     KL[i].delta.append(['q', 'r'])
        #     print(KL[i].delta)
        #     print(KL[i])
#        sys.exit(0)


#        Z1 = evaluate(h1cas, AApq) +  evaluate(h2cas, AApq)
#        Z1 = Z1.fromleft(AAsr)


        
        # testing commutators
        # KL = arithmetic_string(h1cas)
        # KL = KL.fromleft(AApq)
        # KL = KL.fromleft(AAsr)
        # KL2 = arithmetic_string(h2cas)
        # KL2 = KL2.fromleft(AApq)
        # KL2 = KL2.fromleft(AAsr)

        # Z2 = deepcopy(KL + KL2)

        # KL = arithmetic_string(h1cas)
        # KL = KL.fromright(AApq)
        # KL = KL.fromleft(AAsr)
        # KL2 = arithmetic_string(h2cas)
        # KL2 = KL2.fromright(AApq)
        # KL2 = KL2.fromleft(AAsr)


        # Z3 = KL.scale(-1.0) + KL2.scale(-1.0)
        
        # KL = Z1 + Z2 + Z3
 #       sys.exit(0)
        # KL = KL.fromleft(AAqp)
        # KL = KL.fromright(AArs)
#        KL2 = arithmetic_string(h2cas)
        #KL2 = KL2.fromleft(AAqp)
        #KL2 = KL2.fromright(AArs)
#        KL = KL + KL2
        #KL=KL2
#        print(KL)
        # sys.exit(0)
        # # N NEVPTS
        JL = arithmetic_string(h1cas)
        JL = JL.fromright(AApq)
        JL = JL.fromright(AArs)
        JL2 = arithmetic_string(h2cas)
        JL2 = JL2.fromright(AApq)
        JL2 = JL2.fromright(AArs)


        SL =  arithmetic_string(h1cas)
        SL = SL.fromleft(AArs)
        SL = SL.fromleft(AApq)
        SL2 = arithmetic_string(h2cas)
        SL2 = SL2.fromleft(AArs)
        SL2 = SL2.fromleft(AApq)

        
        SL =  arithmetic_string(h1cas)
        SL = SL.fromleft(AArs)
        SL = SL.fromright(AApq)
        SL2 = arithmetic_string(h2cas)
        SL2 = SL2.fromleft(AArs)
        SL2 = SL2.fromright(AApq)

#        KL =  JL + JL2
        KL = ML
#        KL = SL + SL2
        #ML#+ JL.scale(1.0) + JL2.scale(1.0)
        #        KL = KL.scale(-1.0) + KL2.scale(-1.0)
        #        KL = Z1 + KL
        for x in KL:
            print(x)
#        sys.exit(0)
#        if cis == 'ph':
#            KL = KL.fromleft(AAsr)    
#            KL = KL.fromleft(AAqp)    

            # print( 'wynik po fromleft')
            # k=0
            # for x in KL:
            #     print(k, x)
            #     k += 1
            # print('')

        if cis == 'pp':
            KL = KL.fromleft(AArs_aa)    

            print( 'wynik po fromleft')
            k=0
            for x in KL:
                print(k, x)
                k += 1
            print('')


        print('Perform Wick')
        res = []
        for x in KL:
            result_list = x.wick_ca()
            res = res + result_list
        print('Result po WICK:')    
        for x in res:
            x.exec_delta()
            print(x)


        print()
        print('Result po RENAME DENSITY:')
        print()

        for x in res:
            x.rename_as_density()
            print(x)

#        sys.exit(0)
        res2 = cas_to_ugg(res)
        print()
        print('afterc converting to cas')
        print()
        k = 0
        for x in res2:
            print(k, x)
            k+=1
        print()

        # print('simplified final')
        # rsimp = simplify(res2, cas=True)

        # print('')


        print('typeres1', type(res2), type(res2[0]))
#        sys.exit(0)
        # print()
        # rsimp_check = arithmetic_string()
        # for k in range(12,13):
        #     rsimp_check = rsimp_check + arithmetic_string(res2[k])
        #     print(res2[k], res2[k].coefficient)

        # res2 = rsimp_check
#        sys.exit(0)
#        res2 = dirac_to_coulomb(res2)
        print('res2')
        k = 0
        for x in res2:
            print(k, x)
            k+=1
#        sys.exit(0)
        rsimp = simplify(res2, cas=True)
        # print(len(rsimp))

        k = 0
        for x in rsimp:
            print(k, x)
            k+=1
#        sys.exit(0)


        # check
        # print()
        # rsimp_check = arithmetic_string()
        # for k in range(1,2):
        #     rsimp_check = rsimp_check + arithmetic_string(rsimp[k])
        #     print(rsimp[k])

        # rsimp = rsimp_check
        
#        sys.exit(0)
#        print('typeres2', type(rsimp), type(rsimp[0]))


        if cis  =='pp':
             spin_dict1={}
             spin_dict1['r'] = '+'
             spin_dict1['s'] = '-'
             spin_dict1['p'] = '+'
             spin_dict1['q'] = '-'
             numf1 = 1.0
        else:

            spin_dict1={}
            spin_dict1['p'] = '+'
            spin_dict1['q'] = '+'
            #        spin_dict1['t'] = '+'
            #        spin_dict1['u'] = '+'
            
            spin_dict1['r'] = '+'
            spin_dict1['s'] = '+'
            # spin_dict1['v'] = '+'
            # spin_dict1['w'] = '+'
            numf1 = 1.0

        

        rgam = add_spin_driver(spin_dict1, numf1, rsimp)

        rsimp2 = deepcopy(rsimp)

        print('rsimp2', type(rsimp2))
        k = 0
        for x in rsimp2:
            print(k, x)
            k+=1

        

        # check

            
        if cis == 'pp':
            spin_dict2={}
            spin_dict2['r'] = '+'
            spin_dict2['s'] = '-'
            spin_dict2['p'] = '-'
            spin_dict2['q'] = '+'
            numf2 = -1.0
        else:
            spin_dict1={}
            spin_dict1['p'] = '+'
            spin_dict1['q'] = '+'
            # spin_dict1['t'] = '+'
            # spin_dict1['u'] = '+'
            
            spin_dict1['r'] = '-'
            spin_dict1['s'] = '-'
            numf1 = 1.0
            
        rgam2 = add_spin_driver(spin_dict1, numf1, rsimp2)
        # print('typeres3', type(rgam), type(rgam[0]))
#        rgam = rgam + rgam2
        rgam = rgam
        print('')
        print('last result po spin')
        for x in rgam:
            if x.num_factor != 0:
                print(x)

#        rgam = simplify_for_mult_2(rgam)
        rgam = simplify_for_mult_234(rgam)


        print()
        print('po simp and reaname gamma', len(rgam))
        print()
        print(type(rgam))
        k = 0
        for x in rgam:
            print(k, x, x.coefficient_spin)
            k+=1
        print()


        rgam = simplify_final_touch(rgam)
        print('po simplify final touch', len(rgam))
        k = 0
        for x in rgam:
            print(k, x)
            k+=1
        print()


        

        #   rgam_original = deepcopy(rgam)
        #    rNotActive = add_spin_to_gamma(rNotActive)
        rgam_original = add_spin_to_gamma(rgam)
        print('po add spin to gamma')
        k = 0
        for x in rgam_original:
            print(k, x, x.coefficient_spin)
            k+=1
        print()

        A0_blocks = {}
        A1_blocks = {}
        start = time.time()
        for key in Ablock_dict:

            print('this_key', key)
            Ablock_idx = key        
            #            myfixed = ['r','s','p','q']

            myfixed = ['p','q', 'r', 's']
            rgam = deepcopy(rgam_original)
        
            rgamSplit = arithmetic_string()
            for x in rgam:
                for fx in range(0, 4):
                    x.orbital_type[myfixed[fx]] = Ablock_dict[Ablock_idx][fx]

                if len(x.summation) == 0:
                    rgamSplit += arithmetic_string(x)
                else:
                    new_elements = split_element(x)
                    rgamSplit += new_elements

                # for x in rgamSplit:
                #     print(x, x.orbital_type)

            print('len1', len(rgamSplit))
            remove_terms_with_virtual_in_any_DM(rgamSplit)
            print('len po virt', len(rgamSplit))
            remove_terms_with_delta_between_different_sets(rgamSplit)
            print('len po delta', len(rgamSplit))
            print()
            for x in rgamSplit:
                print(x)

            remove_terms_with_non_diag_h(rgamSplit)
            print('len po hdiag', len(rgamSplit))
            print()
            for x in rgamSplit:
                print(x)
            print()
#            remove_non_canonical_h(rgamSplit)
            print('popopopo')
            print()
            for x in rgamSplit:
                print(x)
            print()
            
            k = 0
            # for x in rgamSplit:
            #     print(k, x, x.orbital_type)
            #     k+=1
            # print()

            rgamSplit0, rgamSplit1 = divide_integrals_for_a0_a1(rgamSplit)
            print('len po spli0 split1', len(rgamSplit0), len(rgamSplit1))
            k = 0

            print('split0')
            for x in rgamSplit0:
                print(k, x, x.orbital_type)
                k+=1
            print()

            print('split1')
            for x in rgamSplit1:
                print(k, x, x.orbital_type)
                k+=1
            print()

            
            A0_blocks[key] = deepcopy(rgamSplit0)
            A1_blocks[key] = deepcopy(rgamSplit1)


            #        outA0 = open('./pickle/A0_blocks.pkl','wb')
            #        outA1 = open('./pickle/A1_blocks.pkl', 'wb')


            #        pickle.dump(A0_blocks, outA0)
            #        pickle.dump(A1_blocks, outA1)


            #    elif pick == 'dump2':
    if aaa == 1:
        #       outA0 = open('./pickle/A0_blocks.pkl','rb')
        #      outA1 = open('./pickle/A1_blocks.pkl', 'rb')

 
        #     A0_blocks = pickle.load(outA0)
        #    A1_blocks = pickle.load(outA1)
 

        print()
        print('A0 blocks')
        plusz = ['aaaa']
        A0_blocks_simp = {}
        #print(len(A0_blocks['aaaa']))
        #A0_blocks['aaaa'] = arithmetic_string(A0_blocks['aaaa'][101] )
        
        for key in A0_blocks:
            print()
            
            print('teraz robie klucz', key, len(A0_blocks[key]))
            # for x in A0_blocks[key]:
            #     print(x, x.orbital_type)
            # print()
            if len(A0_blocks[key]) > 0:
                print('ja')
                for x in A0_blocks[key]:
                    print(x)
                print('ja')
                res4 = decompose_gm4(A0_blocks[key])

                res3 = decompose_gm3(res4)

                
                res2 = decompose_gm2(res3, 0)
                # print('resssw2')
                # for x in res2:
                #     print(x)
                # sys.exit(0)

                
                ress = decompose_gm2(res2, 1)
                for x in ress:
                    print(x)
                print()
                ress = rename_gm1(ress)


                for x in ress:
                    #x.coefficient[:] = [coef if coef != DENS1 else DENSN for coef in x.coefficient]
                    x.exec_delta_cas_fixed()
                    print(x)
                rsimp = simplify(ress, cas=True)
                print('ostateczny wynik', len(ress), len(rsimp), key)
#                sys.exit(0)
                A0_blocks_simp[key] = deepcopy(rsimp)

                k = 0
                for xx in rsimp:
                    #print(k, xx)#, xx.orbital_type)
                    print(xx)
                    k+=1

                # if key == 'aoaa':
                #     sys.exit(0)
                
            print()
        end = time.time()

        print()
        print('A1 blocks')
        A1_blocks_simp = {}

        # print(len(A1_blocks['aoaa']))
        # A1_blocks['aoaa'] = arithmetic_string(A1_blocks['aoaa'][19] )
        # A1_blocks['aaao'] = arithmetic_string(A1_blocks['aaao'][19] )
        # # print(A1_blocks['aoaa'][0], A1_blocks['aoaa'][0].orbital_type)
        # plusz = ['aoaa', 'aaao']
        # for key in plusz:
        for key in A1_blocks:
            print()
            
            print('teraz robie klucz A1', key, len(A1_blocks[key]))

            if len(A1_blocks[key]) > 0:
                for gg in A1_blocks[key]:
                    print(len(A1_blocks[key]), 'elem', gg)
                #print('elem', A1_blocks[key])
                a1t = time.time()
                res4 = decompose_gm4(A1_blocks[key])
                a2t = time.time()
                # if (a2t-a1t) > 1.e-1:
                #     print(f"Czas res4: {a2t - a1t:.4f} sekund", key)

                print('res4 A1', 'len', len(res4), key)
                for xx in res4:
                    if len(xx.summation) == 0:
                        print(xx)
                print()

                
                res3 = decompose_gm3(res4)#A1_blocks[key])

                print('res3 A1', 'len', len(res3), key)
                for xx in res3:
                    print(xx)
#                    if len(xx.summation) == 0:
#                        print(xx)
                print()


                # print()
                # print('res3 A1', 'len', len(res3))
                # for xx in res3:
                #     print(xx)
                # print()

                res2 = decompose_gm2(res3, 0)
                print('res2 A1', 'len', len(res2), key)
                for xx in res2:
                    if len(xx.summation) == 0:
                        print(xx)
                print()

                # print()
                # print('res2 A1', 'len', len(res2))
                # for xx in res2:
                #     print(xx)
                # print()

                
                ress = decompose_gm2(res2, 1)
                print('res2 A1', 'len', len(ress), key)
                for xx in ress:
                    if len(xx.summation) == 0:
                        print(xx)
                print()

                # print('wynik przed simp A1', 'len', len(ress))
                # for xx in ress:
                #     print(xx)

                ress = rename_gm1(ress)
                for x in ress:
#                    x.coefficient[:] = [coef if coef != DENS1 else DENSN for coef in x.coefficient]                    
                    x.exec_delta_cas_fixed()

                print('pluszoprzed')
                for y in ress:
                    print(y)
                print()
                rsimp = simplify(ress, cas=True)
                print('pluszopo')
                for y in rsimp:
                    print(y)
                print()

                print('ostateczny wynik simp A1', len(ress), len(rsimp), key)
                A1_blocks_simp[key] = deepcopy(rsimp)

                k = 0
                for xx in rsimp:
                    #print(k, xx)#, xx.orbital_type)
                    print(xx)
                    k+=1

                
            print()
        end = time.time()

        # outA0 = open('./pickle/A0_blocks_decomposed2.pkl','wb')
        # outA1 = open('./pickle/A1_blocks_decomposed2.pkl', 'wb')


        # pickle.dump(A0_blocks_simp, outA0)
        # pickle.dump(A1_blocks_simp, outA1)


    #else:
    if aaa == 1:    

        print('pluszek')
        # outA0 = open('./pickle/A0_blocks_decomposed2.pkl','rb')
        # outA1 = open('./pickle/A1_blocks_decomposed2.pkl', 'rb')


        # A0_blocks = pickle.load(outA0)
        # A1_blocks = pickle.load(outA1)

        A0_blocks = A0_blocks_simp
        A1_blocks = A1_blocks_simp

        # for key in A0_blocks:
        #     if (len(A0_blocks[key]))> 0:
        #         for x in A0_blocks[key]:
        #             pr
        
       

        for key in A0_blocks:
            if (len(A0_blocks[key]))> 0:
                print(key, len(A0_blocks[key]))
                k = 0
                for x in A0_blocks[key]:
#                    print(k, x)
                    seen = set()
                    unique_delta = []
                    for item in x.delta:
                        item_str = str(item)
                        if item_str not in seen:
                            seen.add(item_str)
                            unique_delta.append(item)
                    x.delta = unique_delta
                    x.exec_delta_cas_fixed()

#                    print(k, x.delta)
#                    print(k, x)#, x.orbital_type)

 #               print()
                A0_blocks[key] = simplify(A0_blocks[key], cas=True)
#        print('lll', len(A0_blocks['aaaa']))
       
        for key in A1_blocks:
            if (len(A1_blocks[key]))> 0:
                print(key, len(A1_blocks[key]))
                k = 0
                for x in A1_blocks[key]:
                    seen = set()
                    unique_delta = []
                    for item in x.delta:
                        item_str = str(item)
                        if item_str not in seen:
                            seen.add(item_str)
                            unique_delta.append(item)
                    x.delta = unique_delta
                    x.exec_delta_cas_fixed()
                    
                    #print(k, x.delta)
                    k+=1
                print()
                A1_blocks[key] = simplify(A1_blocks[key], cas=True)

        end = time.time()

        print('final print')
        print()
        print('A0*******************************************************************************************************************************************')
        print('minipluszek')

        print('')
        print()
        
        for key in A0_blocks:
            if (len(A0_blocks[key]))> 0:
                print(f"+{'-'*30}+")
                print(f"| a0 block {key} {len(A0_blocks[key]):^10} |")
                print(f"+{'-'*30}+")
                all_g0 = []
                # for x in A0_blocks[key]:
                #     for j in range(0, len(x.coefficient)):
                #         if x.coefficient[j] ==   DENS3PM or x.coefficient[j] ==   DENS2P or x.coefficient[j] ==   DENS4PM or x.coefficient[j] ==   DENS4PPM:
                #             if len(set(x.coefficient_idx[j])) >1:                                                    
                #                 print(x, len(set(x.coefficient_idx[j])))
                #                 this = mark_repeated_indices(x.coefficient_idx[j])
                #                 if this not in all_g0:
                #                     all_g0.append(this)

                # print('required rmd')
                # for x in all_g0:
                #     print(x)
                #     print()

                print_struct(A0_blocks[key], with_ints = False)
                k = 0
                # for x in A0_blocks[key]:
                #     print(k, x)
                #     k+=1
                # print()
                print('================================================================================')

        print()
        print('A1*******************************************************************************************************************************************')
        print()
        
        for key in A1_blocks:
            if (len(A1_blocks[key]))> 0:
                print(f"+{'-'*30}+")
                print(f"| a1 block {key} {len(A1_blocks[key]):^10} |")
                print(f"+{'-'*30}+")
                # all_g1 = []
                # for x in A1_blocks[key]:
                #     for j in range(0, len(x.coefficient)):
                #         if x.coefficient[j] ==   DENS3PM or x.coefficient[j] ==   DENS2P or x.coefficient[j] ==   DENS4PM or x.coefficient[j] ==   DENS4PPM:
                #             if len(set(x.coefficient_idx[j])) >1:                                                    
                #                 print(x, len(set(x.coefficient_idx[j])))
                #                 this = mark_repeated_indices(x.coefficient_idx[j])
                #                 if this not in all_g1:
                #                     all_g1.append(this)

                # print('required rmd')
                # for x in all_g1:
                #     print(x)
                #     print()
                k = 0
                print_struct(A1_blocks[key], with_ints = False)
                # for x in A1_blocks[key]:
                #     print(k, x)
                #     k+=1
                print('================================================================================')

        
        #        print(f"Czas wykonania: {end - start:.2f} sekund")
        

def mark_repeated_indices(lst):
    counter = {}
    for item in lst:
        counter[item] = counter.get(item, 0) + 1
    
    result = []
    for item in lst:
        if counter[item] > 1:
            result.append(0)
        else:
            result.append(1)
    
    return result

    
def print_struct(r, with_ints = False):

    sumt = []
    sumtu = []
    sumtuv = []
    sumtuvw = []
    sumz = []
    counter_lstz = []
    counter_lstt = []
    counter_lsttu = []
    counter_lsttuv = []
    counter_lsttuvw = []
    printz = {}
    printt = {}
    printtu = {}
    printtuv = {}
    printtuvw = {}
    for x in r:
        counts = Counter(x.coefficient)
        counts_tuple = tuple(sorted(counts.items()))

        if len(x.summation) == 0:
            sumz.append(x)            
            if counts_tuple not in counter_lstz:
                counter_lstz.append(counts_tuple)
                printz[counts_tuple] = []
        elif len(x.summation) == 1:

            sumt.append(x)
            if counts_tuple not in counter_lstt:
                counter_lstt.append(counts_tuple)
                printt[counts_tuple] = []
        elif len(x.summation) == 2:
            sumtu.append(x)
            if counts_tuple not in counter_lsttu:
                counter_lsttu.append(counts_tuple)
                printtu[counts_tuple] = []
        elif len(x.summation) == 3:
            sumtuv.append(x)
            if counts_tuple not in counter_lsttuv:
                counter_lsttuv.append(counts_tuple)
                printtuv[counts_tuple] = []
        elif  len(x.summation) == 4:
            sumtuvw.append(x)
            if counts_tuple not in counter_lsttuvw:
                counter_lsttuvw.append(counts_tuple)
                printtuvw[counts_tuple] = []
        else:
            print('??', x.summation)
        
    
    for x in sumz:
        counts = Counter(x.coefficient)
        counts_tuple = tuple(sorted(counts.items()))
        printz[counts_tuple].append(x)

    for x in sumt:
        counts = Counter(x.coefficient)
        counts_tuple = tuple(sorted(counts.items()))
        printt[counts_tuple].append(x)
    for x in sumtu:
        counts = Counter(x.coefficient)
        counts_tuple = tuple(sorted(counts.items()))
        printtu[counts_tuple].append(x)

    for x in sumtuv:
        counts = Counter(x.coefficient)
        counts_tuple = tuple(sorted(counts.items()))
        printtuv[counts_tuple].append(x)

    for x in sumtuvw:
        counts = Counter(x.coefficient)
        counts_tuple = tuple(sorted(counts.items()))
        printtuvw[counts_tuple].append(x)

    print('{\\scriptsize')    
    print('\\begin{equation}')
    print('\\begin{split}\n&A_{rspq}^{(1)}=')
    ki = 0
#    print('zero summation idices')
#    print()
    for x  in printz.values():
        for elem in x:
            ki += 1
            print_mini(elem, ki, with_ints)

    #     print()
    # print()
    # print('one summation index')
    # print('---------------------------------------------------')
    # print()

    for x in printt.values():
        for t_val in ['a', 'o']:
            matching_elems = [elem for elem in x if elem.orbital_type['t'] == t_val]

            if matching_elems:  # Print only if there are elements in this category
                # if t_val=='o':
                #     print(f't = {t_val}')
                for elem in matching_elems:
                    if elem.num_factor.is_integer() == True:
                        if t_val=='o':
                            elem.summation = ["t'"]
                        ki +=1
                        print_mini(elem, ki, with_ints)

                            
                for elem in matching_elems:
                    if elem.num_factor.is_integer() == False:
                        #                        print(elem, int_type(elem))
                        if t_val=='o':
                            elem.summation = ["t'"]
                        ki +=1
                        print_mini(elem, ki, with_ints)


    # print('two summation idices')
    # print('---------------------------------------------------')
    # print()
    for x in printtu.values():
        for t_val, u_val in [('a', 'a'), ('a', 'o'), ('o', 'a'), ('o', 'o')]:
            matching_elems = [elem for elem in x if 
                              elem.orbital_type['t'] == t_val and 
                              elem.orbital_type['u'] == u_val]

            if matching_elems:  # Print only if there are elements in this category
                # if t_val=='o' or u_val == 'o':
                #     print(f'tu = {t_val}{u_val}')
                for elem in matching_elems:
                    if elem.num_factor.is_integer() == True:
                        if t_val=='o':
                            elem.summation[elem.summation.index("t")] = "t'"
                        if u_val == 'o':
                            elem.summation[elem.summation.index("u")] = "u'"
                        ki +=1
                        print_mini(elem, ki, with_ints)

                for elem in matching_elems:
                    if elem.num_factor.is_integer() == False:
                        if t_val=='o':
                            elem.summation[elem.summation.index("t")] = "t'"
                        if u_val == 'o':
                            elem.summation[elem.summation.index("u")] = "u'"

                        ki +=1
                        print_mini(elem, ki, with_ints)



#                        print(elem, int_type(elem))
    #                     print(elem)
    #             print()
    # print()

    # print('three summation idices', len(printtuv))
    # print('---------------------------------------------------')
    for x in printtuv.values():
        for t_val, u_val, v_val in [('a', 'a', 'a'), ('a', 'a', 'o'), ('a', 'o', 'a'), ('a', 'o', 'o'),
                                    ('o', 'a', 'a'), ('o', 'a', 'o'), ('o', 'o', 'a'), ('o', 'o', 'o')]:
            matching_elems = [elem for elem in x if 
                              elem.orbital_type['t'] == t_val and 
                              elem.orbital_type['u'] == u_val and 
                              elem.orbital_type['v'] == v_val]

            if matching_elems:  # Print only if there are elements in this category
                # if t_val=='o' or u_val == 'o' or v_val == 'o':                                    
                #     print(f'tuv = {t_val}{u_val}{v_val}')
                for elem in matching_elems:
                    if elem.num_factor.is_integer() == True:
#                        print(elem, int_type(elem))
                        if t_val=='o':
                            elem.summation[elem.summation.index("t")] = "t'"
                        if u_val == 'o':
                            elem.summation[elem.summation.index("u")] = "u'"
                        if v_val == 'o':
                            elem.summation[elem.summation.index("v")] = "v'"

                        ki +=1
                        print_mini(elem, ki, with_ints)

                        #print(elem)
                for elem in matching_elems:
                    if elem.num_factor.is_integer() == False:
                        if t_val=='o':
                            elem.summation[elem.summation.index("t")] = "t'"
                        if u_val == 'o':
                            elem.summation[elem.summation.index("u")] = "u'"
                        if v_val == 'o':
                            elem.summation[elem.summation.index("v")] = "v'"


                        ki +=1
                        print_mini(elem, ki, with_ints)

    for x in printtuvw.values():
        for t_val, u_val, v_val, w_val in [('a', 'a', 'a', 'a'), ('a', 'a', 'a', 'o'), ('a', 'a', 'o', 'a'), ('a', 'a', 'o', 'o'),
                                    ('a', 'o', 'a', 'a'), ('a', 'o', 'a', 'o'), ('a', 'o', 'o', 'a'), ('a', 'o', 'o', 'o'),
                                    ('o', 'a', 'a', 'a'), ('o', 'a', 'a', 'o'), ('o', 'a', 'o', 'a'), ('o', 'a', 'o', 'o'),
                                    ('o', 'o', 'a', 'a'), ('o', 'o', 'a', 'o'), ('o', 'o', 'o', 'a'), ('o', 'o', 'o', 'o')]:
            matching_elems = [elem for elem in x if 
                              elem.orbital_type['t'] == t_val and 
                              elem.orbital_type['u'] == u_val and 
                              elem.orbital_type['v'] == v_val and
                              elem.orbital_type['w'] == v_val]

            if matching_elems: 

                for elem in matching_elems:
                    if elem.num_factor.is_integer() == True:
                        if t_val=='o':
                            elem.summation[elem.summation.index("t")] = "t'"
                        if u_val == 'o':
                            elem.summation[elem.summation.index("u")] = "u'"
                        if v_val == 'o':
                            elem.summation[elem.summation.index("v")] = "v'"
                        if w_val == 'o':
                            elem.summation[elem.summation.index("w")] = "w'"

                        ki +=1
                        print_mini(elem, ki, with_ints)

                        #print(elem)
                for elem in matching_elems:
                    if elem.num_factor.is_integer() == False:
                        if t_val=='o':
                            elem.summation[elem.summation.index("t")] = "t'"
                        if u_val == 'o':
                            elem.summation[elem.summation.index("u")] = "u'"
                        if v_val == 'o':
                            elem.summation[elem.summation.index("v")] = "v'"
                        if w_val == 'o':
                            elem.summation[elem.summation.index("w")] = "w'"


                        ki +=1
                        print_mini(elem, ki, with_ints)
                        

                # if t_val=='o' or u_val == 'o' or v_val == 'o':
                #     print()

#                        print(elem, int_type(elem))
    #                     print(elem)
    #             print()
    # print()
    print('\\end{split}')
    print('\\end{equation}')
    print('}')

def print_mini(elem, ki, with_ints):

    if with_ints == True:
        if (ki%5==0):
            print(elem+"\\\&", int_type(elem))
        else:
            print(elem, int_type(elem))
    else:
        if (ki%5==0):
            print(elem+"\\\&")
        else:
            print(elem)

    
def int_type(elem):

    int_dict = {
    'IIII': ['oooo'],
    'IAII': ['oaoo', 'aooo', 'ooao', 'oooa'],
    'IAIA': ['oaoa', 'oaao', 'aooa', 'aoao'],
    'AAII': ['aaoo', 'ooaa'],
    'AAIA': ['aaoa', 'aaao', 'aoaa', 'oaaa'],
    'AAAA': ['aaaa'],

    'IVII': ['ovoo', 'vooo', 'ooov', 'oovo'],
    'IVIA': ['ovoa', 'ovao', 'vooa', 'voao', 'oaov', 'oavo', 'aoov', 'aovo'],
    'IVAA': ['ovaa', 'voaa', 'aaov', 'aavo'],
    'IVIV': ['ovov', 'ovvo', 'voov', 'vovo'],

    'VAII': ['vaoo', 'avoo', 'oova', 'ooav'],
    'VAIA': ['vaoa', 'vaao', 'avoa', 'avao', 'oava', 'oaav', 'aova', 'aoav'],
    'VAIV': ['vaov', 'vavo', 'avov', 'avvo', 'ovva', 'ovav', 'vova', 'voav'],
    'VAAA': ['vaaa', 'avaa', 'aava', 'aaav'],
    'VAVA': ['vava', 'vaav', 'avva', 'avav'],

    'VVII': ['vvoo', 'oovv'],
    'VVIA': ['vvoa', 'vvao', 'oavv', 'aovv'],
    'VVIV': ['vvov', 'vvvo', 'ovvv', 'vovv'],
    'VVAA': ['vvaa', 'aavv'],
    'VVVA': ['vvva', 'vvav', 'avvv', 'vavv'],
    'VVVV': ['vvvv'],
}

    
    if TWOEL_INT in elem.coefficient:
        for i, coef in enumerate(elem.coefficient):
            if coef == TWOEL_INT:
                intg = elem.coefficient_idx[i]

                ints = ""
                for j in intg:
                    ints += elem.orbital_type[j]
                for intb in int_dict:
                    if ints in int_dict[intb]:
                        return intb
                
        
    else:
        return ""

    return "bbb"
        
def decompose_gm2(r, tt):
#    print('BEFORE decomposing Gm2', len(r))

    # kk = 0
    # for x in r:
    #     print(kk, x)
    #     kk+=1
    # print()
    res = arithmetic_string()
    res2 = arithmetic_string()
    kk = -1
    for elem in r:
        kk+= 1

        #        has_gm2 = any('Gm2' in c for coefficient in elem.coefficient for c in elem.coefficient)

        gm2_idx = []
        for i, coef in enumerate(elem.coefficient):
            if 'Gm2' in coef:
                gm2_idx.append(i)

        
        
        if tt == 1:
            if len(gm2_idx)>1:
                has_gm2 = True
                gm2_idx = [gm2_idx[1]]
            else:
                has_gm2 = False
        elif tt ==0:
            if len(gm2_idx)==0:
                has_gm2 = False
            else:
                has_gm2 = True
                gm2_idx = [gm2_idx[0]]
        print(has_gm2, elem, elem.coefficient, elem.orbital_type)
        if has_gm2:
            # gm2_idx = []
            # for i, coef in enumerate(elem.coefficient):
            #     if 'Gm2' in coef:
            #         gm2_idx.append(i)

            # if tt == 0:

            # elif tt == 1:
            #     if len(gm2_idx) > 1:

            #     else
            # print('decomposing this', kk, elem, 'lenres', len(res))
            # print()
            for i in gm2_idx:
                #for i, coef in enumerate(elem.coefficient):
#                    print(kk, 'this_elem have G2 on place',i, elem, elem.coefficient_spin)
                    upper = elem.coefficient_idx[i][2:4]
                    lower = elem.coefficient_idx[i][0:2]

                    upper_spin = elem.coefficient_spin[i][2:4]
                    lower_spin = elem.coefficient_spin[i][0:2]
                    #print(upper, upper_spin)
                    for l in range(0, 2):
                        upper[l] = upper[l] + upper_spin[l]
                        lower[l] = lower[l] + lower_spin[l]

                    type_dict = {}
                    for l in range(0, 2):
                        key = upper[l]
                        type_dict[key] = elem.orbital_type[key[0]]
                        key = lower[l]
                        type_dict[key] = elem.orbital_type[key[0]]

                    are_active = check_active(upper, lower, type_dict)
                    if are_active:
                        # print('all indices active, leave G2')
                        # print(upper, lower, type_dict, elem)
                        
                        res += arithmetic_string(elem)
                    else:
                        print('upper', upper)
                        print('lower', lower)
                        print('type', type_dict)
                        res_this = decompose_g2_ul(upper, lower, type_dict)
                        
                        print('final g2', elem)
                        print('result of len', len(res_this))
                        for d in res_this:
                            print(d)
                        
                        final = convert_to_ugg(elem, res_this, type_dict, [i])

                        print('ugg result for G2 this', elem, type_dict)
                        print(len(final))
                        print('***************')
                        print('finall', type_dict)
                        for x in final:
                            print(x, x.coefficient, x.coefficient_idx, x.orbital_type)
                        
                        rsimp = simplify(final, cas=True)
                        print()

                        print('po simplify', len(rsimp))
                        for x in rsimp:
                            print(x.coefficient_idx)
                            print(x, x.orbital_type)
                        print()
                        print()
                        print('---------------------------------------------------')
                        
                        rsimp.exec_delta_cas()
                        rsimp.cleanup()
                        # print()
                        # print('po cleanup')
                        # print('original elem', elem)

                        # print(type_dict)
                        # print()
                        # for x in rsimp:
                        #     print(x)

                        # print()

                        # print('ugg result for G2 this', elem, type_dict)
                        for x in final:
                            res = res + arithmetic_string(x)
                            # print(x, x.orbital_type)

        else:
            # print('not decomposing this', kk, elem, 'lenres2', len(res2))

            res2 = res2 + arithmetic_string(elem)

    # print('decoplusz', len(res))
    # for x in res:
    #     print(x)
    # print()

    # print('deconieplusz', len(res2))
    # for x in res2:
    #     print(x)
    # print()

    res = res + res2

            
    rsimp = simplify(res, cas=True)
    # print('this is what we get g2', len(res), len(rsimp))
    # for x in rsimp:
    #     print(x)
    # print()
    # print()
    rsimp.exec_delta_cas()
    rsimp.cleanup()
    # print()
    # print('po cleanup')
    # print('this is what we get g2 po cleanup', len(res), len(rsimp))
    # for x in rsimp:
    #     print(x)

    return rsimp

def decompose_gm3(r):


    # print('BEFORE decomposing Gm3')

    # print()
    # k = 0
    # for x in r:
    #     print(k, x)
    #     k+=1
    # print()
    

    res = arithmetic_string()
    res2 = arithmetic_string()
    kk = -1
    for elem in r:
        kk+= 1
        #print(elem.coefficient)
        has_gm3 = any('Gm3' in c for coefficient in elem.coefficient for c in elem.coefficient)
        # print(has_gm3)
        if has_gm3:
            # print('decomposing this', kk, elem)
            for i, coef in enumerate(elem.coefficient):
                #if coef == DENS3PM or coef ==:
                if 'Gm3' in coef:
#                    print(kk, 'this_elem have G3', elem)
                    upper = elem.coefficient_idx[i][3:6]
                    lower = elem.coefficient_idx[i][0:3]

                    upper_spin = elem.coefficient_spin[i][3:6]
                    lower_spin = elem.coefficient_spin[i][0:3]

                    for l in range(0, 3):
                        upper[l] = upper[l] + upper_spin[l]
                        lower[l] = lower[l] + lower_spin[l]

                    type_dict = {}
                    for l in range(0, 3):

                        key = upper[l]
                        type_dict[key] = elem.orbital_type[key[0]]
                        key = lower[l]
                        type_dict[key] = elem.orbital_type[key[0]]
                        
                    #print(upper, lower, type_dict)


                    # spin_dict = {}
                    # for j, x in enumerate(elem.coefficient_idx[i]):
                    #     spin_dict[x] = elem.coefficient_spin[i][j]
                    # type_dict = elem.orbital_type

                    # print(spin_dict)
                    # print(type_dict)

                    # print()


                    # sprawdz czy wszystkie indeksy są active. jeśli tak to zostaw
                    are_active = check_active(upper, lower, type_dict)#elem.orbital_type)
                    if are_active:
#                        print('all indices active, leave G3')
                        
                        res += arithmetic_string(elem)
                    # jeśli nie, to za    
                    else:

                        res_this = decompose_g3_ul(upper, lower, type_dict)

                        # print('final g3', elem)
                        # print('result of len', len(res_this))
                        

                        final = convert_to_ugg(elem, res_this, type_dict, [i])

                        # print('ugg result for G3 this', elem, type_dict)
                        # print(len(final))
                        # print('***************')

                        
                        rsimp = simplify(final, cas=True)
                        # print()
                        # print('po simplify', len(rsimp))
                        # for x in rsimp:
                        #     print(x.coefficient_idx)
                        #     print(x, x.orbital_type)

                        
                        rsimp.exec_delta_cas()
                        rsimp.cleanup()
                        # print()
                        # print('po cleanup')
                        # print('original elem', elem)

                        # print(type_dict)
                        # print()
                        # for x in rsimp:
                        #     print(x)

                        # print()

                        # print('ugg result for G3 this', elem, type_dict)
                        for x in final:

                            res = res + arithmetic_string(x)
                            # print(x, x.orbital_type)

        else:
            # print('not decomposing this', kk, elem)
            res2 = res2 + arithmetic_string(elem)

    # print('decoplusz')
    # for x in res:
    #     print(x)
    # print()

    # print('deconieplusz')
    # for x in res2:
    #     print(x)
    # print()

    res = res + res2

            
    rsimp = simplify(res, cas=True)
    # print('this is what we get g3', len(res), len(rsimp))
    # for x in rsimp:
    #     print(x)
    # print()
    # print()
    rsimp.exec_delta_cas()
    rsimp.cleanup()
    # print()
    # print('po cleanup')
    # print('this is what we get g3 po cleanup', len(res), len(rsimp))
    # for x in rsimp:
    #     print(x)

    return rsimp


def decompose_g2_ul(upper, lower, type_dict):
    all_key_list = ['upper_g1', 'lower_g1', 'upper_g2', 'lower_g2', 'upper_g3', 'lower_g3', 'upper_c2', 'lower_c2', 'upper_c3', 'lower_c3', 'numf']
    
    base_dict = {key: [] for key in all_key_list}
    
    cases = [
        {
            'upper_g1': [[upper[0]], [upper[1]]],
            'lower_g1': [[lower[0]], [lower[1]]],
            'numf': 1.0
        },
        {
            'upper_g1': [[upper[0]], [upper[1]]],
            'lower_g1': [[lower[1]], [lower[0]]],
            'numf': -1.0
        }
    ]
    
    result = []
    for case in cases:
        print('case', case)
        dict_case = deepcopy(base_dict)
        dict_case.update(case)
        result.append(dict_case)
        print(dict_case)
        print()
    result = remove_zeros(result, type_dict)
        
    return result

def decompose_g3_ul(upper, lower, type_dict):
                    

    res_1111 = partition_indices_1111(upper, lower, 1.0)
    res_12 = partition_indices_12_13(upper, lower, 1.0)


    res_1111 = remove_zeros(res_1111, type_dict)
    res_12 = remove_zeros(res_12, type_dict)

    res_12 = cumulant2_to_gamma2(res_12, type_dict)
    res_12 = remove_zeros(res_12, type_dict)

    res_this = res_1111 + res_12
    res_this = cancel_out_dicts(res_this)
    return res_this


def decompose_gm4(r):

    res = arithmetic_string()
    kk = 0
    for elem in r:
        t1 = time.time()
        has_gm4 = any('Gm4' in c for coefficient in elem.coefficient for c in elem.coefficient)
        a1t_0 = time.time()
        if has_gm4:
            t3 = time.time()
            for i, coef in enumerate(elem.coefficient):
                if 'Gm4' in coef:
                    t9 = time.time()
                    t5 = time.time()

                    kk += 1
                    upper = elem.coefficient_idx[i][4:8]
                    lower = elem.coefficient_idx[i][0:4]

                    upper_spin = elem.coefficient_spin[i][4:8]
                    lower_spin = elem.coefficient_spin[i][0:4]
                    a2t = time.time()

                    for l in range(0, 4):
                        upper[l] = upper[l] + upper_spin[l]
                        lower[l] = lower[l] + lower_spin[l]

                    type_dict = {}
                    for l in range(0, 4):

                        key = upper[l]
                        type_dict[key] = elem.orbital_type[key[0]]
                        key = lower[l]
                        type_dict[key] = elem.orbital_type[key[0]]

                    a2t = time.time()

                    are_active = check_active(upper, lower, type_dict)
                    a2t = time.time()

                    t10 = time.time()

                    if are_active:
                 #       print('all indices active, leave G4')
                              
                        res += arithmetic_string(elem)
                    else:
                        t7 = time.time()
                        a1t = time.time()
                        res_this = decompose_g4_ul(upper, lower, type_dict)
                        a2t = time.time()

                        print(len(res_this))
                        for t in res_this:
                            if len(t['upper_g1']) == 4 and len(t['upper_g2']) == 0 and len(t['upper_c2']) == 0 and len(t['upper_g3']) == 0:
                                print(t['upper_g1'], t['numf'])
                                print(t['lower_g1'])

                        a1t = time.time()
                        final = convert_to_ugg(elem, res_this, type_dict, [i])
                        a2t = time.time()

                        #print('final', len(final))
                        
                        rsimp = simplify(final, cas=True)
                        
                        rsimp.exec_delta_cas()
                        rsimp.cleanup()
                        
                        for x in rsimp:
                            if len(x.summation) == 0:
                                if DENS2P not in x.coefficient and DENS2PM not in x.coefficient:
                                    print(x, elem)
                        print('---')
                            

                        
                        res = res + final

        else:
            res += arithmetic_string(elem)



    rsimp = simplify(res, cas=True)


    rsimp.exec_delta_cas()
    rsimp.cleanup()

    return rsimp


def decompose_g4_ul(upper, lower, type_dict):

    all_key_list = ['upper_g1', 'lower_g1', 'upper_g2', 'lower_g2', 'upper_g3', 'lower_g3', 'upper_c2', 'lower_c2', 'upper_c3', 'lower_c3']

    res_1111 = partition_indices_1111(upper, lower, 1.0)

    res_13 = partition_indices_12_13(upper, lower, 1.0)

    res_22 = partition_indices_2_2(upper, lower, 1.0)

    res_112 = partition_indices_112(upper, lower, 1.0)

    res_1111 = remove_zeros(res_1111, type_dict)

    for x in res_1111:
        if len(x['upper_g1'])==4:
            print('miniplusz 1111')
            print(x['upper_g1'])
            print(x['lower_g1'])
            print()
  #  print()
    res_13 = remove_zeros(res_13, type_dict)
    for x in res_13:
        if len(x['upper_g1'])==4:
            print('miniplusz 13')
            print(x['upper_g1'])
            print(x['lower_g1'])
            print()
  #  print()
    res_22 = remove_zeros(res_22,  type_dict)
    for x in res_22:
        if len(x['upper_g1'])==4:
            print('miniplusz 22')
            print(x['upper_g1'])
            print(x['lower_g1'])
            print()
  #  print()
    res_112 = remove_zeros(res_112, type_dict)
    for x in res_1111:
        if len(x['upper_g1'])==4:
            print('miniplusz 112')
            print(x['upper_g1'])
            print(x['lower_g1'])
            print()
  #  print()
    
    res_22 = cumulant2_to_gamma2(res_22, type_dict)
    for x in res_22:
        if len(x['upper_g1'])==4:
            print('miniplusz 22 cumul')
            print(x['upper_g1'])
            print(x['lower_g1'])
            print()
  #  print()
    res_112 = cumulant2_to_gamma2(res_112, type_dict)
    for x in res_112:
        if len(x['upper_g1'])==4:
            print('miniplusz 112 cumul')
            print(x['upper_g1'])
            print(x['lower_g1'])
            print()
  #  print()

    res_13 = cumulant3_to_gamma3(res_13, type_dict)
    for x in res_13:
        if len(x['upper_g1'])==4:
            print('miniplusz 13 cumul 33')
            print(x['upper_g1'])
            print(x['lower_g1'])
            print()
 #   print()


    res_13 = cumulant2_to_gamma2(res_13, type_dict)
    for x in res_13:
        if len(x['upper_g1'])==4:
            print('miniplusz 13 cumul 22')
            print(x['upper_g1'])
            print(x['lower_g1'])
            print()
#    print()

    res_13_decomp = []
    for g3elem in res_13:

        if len(g3elem['upper_g3']) == 1:

            upper_g3 = g3elem['upper_g3'][0]
            lower_g3 = g3elem['lower_g3'][0]

            res_this = decompose_g3_ul(upper_g3, lower_g3, type_dict)
            
                
            for x in res_this:
                new_elem = deepcopy(g3elem)
                new_elem['upper_g3'] = []
                new_elem['lower_g3'] = []
                for key in all_key_list:
                    new_elem[key].extend(x[key])
                #print()
                new_elem['numf'] *= x['numf']
                # if len(new_elem['upper_g1'])==4:
                #     print('miniplusz2 13')
                #     print(new_elem['upper_g1'])
                #     print(new_elem['lower_g1'])
                #     print()
                res_13_decomp.append(new_elem)
        else:
            res_13_decomp.append(g3elem)

    for x in res_13_decomp:
        if len(x['upper_g1'])==4:
            print('miniplusz 13 po decompose')
            print(x['upper_g1'])
            print(x['lower_g1'])
            print()
#    print()

    res_this = res_1111 + res_13_decomp + res_22 + res_112
    
#    print()
    for x in res_this:
        if x['upper_g3'] !=[]:
            sys.exit(0)
#    print()
    res_this = cancel_out_dicts(res_this)
    for x in res_this:
        if len(x['upper_g1'])==4:
            print('miniplusz res this fin')
            print(x['upper_g1'])
            print(x['lower_g1'])
            print()
    return res_this


def cancel_out_dicts(dict_list):

    new_dict_list = []
    for i, d in enumerate(dict_list):
        dict_without_numf = {k: v for k, v in d.items() if k != 'numf'}
        dict_key = str(dict_without_numf)

        mini_dict = {}
        mini_dict['hash'] = dict_key
        mini_dict['numf'] = d['numf']
        mini_dict['idx'] = i
        new_dict_list.append(mini_dict)

    for i in range(0, len(new_dict_list)):
        d = new_dict_list[i]
        if d['numf']!=0.0:

            for j in range(i+1, len(new_dict_list)):
                dd = new_dict_list[j]
                if dd['numf']!=0.0:
                    if d['hash'] == dd['hash']:
                        if d['numf'] * dd['numf'] < 0:
                            new_dict_list[i]['numf'] = 0.0
                            new_dict_list[j]['numf'] = 0.0
                            break

    result = []
    for i in range(0, len(new_dict_list)):
        if new_dict_list[i]['numf']!=0.0:
            result.append(dict_list[i])
            
    return result

def presort_gm2_gm3(res, type_dict):

    # name_g2 ++ +- -+ -- 
    # name_g3 +++ ++- +-+ -++ +-- -+- --+ ---

    for dct in res:
        dct['name_g2'] = []
        dct['name_g3'] = []
        if len(dct['upper_g2']) > 0:
          

            for i, g in enumerate(dct['upper_g2']):
                upper = dct['upper_g2'][i]
                lower = dct['lower_g2'][i]
                
                swap_upper = sort_spin_acorrding_to(upper, 'g2')
                swap_lower = sort_spin_acorrding_to(lower, 'g2')
                dct['name_g2'].append(swap_upper[0])
                
                dct['upper_g2'][i] = swap_upper[1]
                dct['lower_g2'][i] = swap_lower[1]
                
                dct['numf'] *= swap_upper[2] * swap_lower[2]

        if len(dct['upper_g3']) > 0:
            dct['name_g3'] = []

            for i, g in enumerate(dct['upper_g3']):
                upper = dct['upper_g3'][i]
                lower = dct['lower_g3'][i]
                
                swap_upper = sort_spin_acorrding_to(upper, 'g3')
                swap_lower = sort_spin_acorrding_to(lower,  'g3')
                dct['name_g3'].append(swap_upper[0])
                
                dct['upper_g3'][i] = swap_upper[1]
                dct['lower_g3'][i] = swap_lower[1]
                dct['numf'] *= swap_upper[2] * swap_lower[2]
                   
    return res
def sort_spin_acorrding_to(lst, g):

    pl = '+'
    mn = '-'
    sort_spin = True
    nosort_spin = False
    
    if g == 'g2':

        g2_dict = {'++': ['Gm2p', [lst[0], lst[1]], 1.0],
               '+-': ['Gm2pm', [lst[0], lst[1]], 1.0],
               '-+': ['Gm2pm', [lst[1], lst[0]], -1.0],
               '--': ['Gm2p', [lst[0], lst[1]], 1.0]
        }

        i0 = lst[0]
        i1 = lst[1]
        this = f"{i0[1]}{i1[1]}"
        return g2_dict[this]

    elif g == 'g3':
        g3_dict = {'+++': ['Gm3p',  [lst[0], lst[1], lst[2]],  1.0],
                   '++-': ['Gm3pm', [lst[0], lst[1], lst[2]],  1.0],
                   '+-+': ['Gm3pm', [lst[0], lst[2], lst[1]], -1.0],
                   '-++': ['Gm3pm', [lst[1], lst[2], lst[0]],  1.0],
                   '+--': ['Gm3pm', [lst[1], lst[2], lst[0]],  1.0],
                   '-+-': ['Gm3pm', [lst[0], lst[2], lst[1]], -1.0],
                   '--+': ['Gm3pm', [lst[0], lst[1], lst[2]],  1.0],
                   '---': ['Gm3pp', [lst[0], lst[1], lst[2]],  1.0]
                   }
    

        i0 = lst[0]
        i1 = lst[1]
        i1 = lst[2]
        this = f"{i0[1]}{i1[1]}{i2[1]}"
        return g3_dict[this]

    else:
        print('UnexpectedCaseError')
        sys.exit(1)
    


def convert_to_ugg(elem, res,  type_dict, del_list):

    res = presort_gm2_gm3(res, type_dict)
    

    result = arithmetic_string()

    new_elem_base = ugg()
    new_elem_base.num_factor = deepcopy(elem.num_factor)
    new_elem_base.orbital_type = deepcopy(elem.orbital_type)
    new_elem_base.summation = deepcopy(elem.summation)
    new_elem_base.delta = deepcopy(elem.delta)

    for i, coef in enumerate(elem.coefficient):
        if i not in del_list:
            new_elem_base.coefficient.append(elem.coefficient[i])
            new_elem_base.coefficient_idx.append(elem.coefficient_idx[i])
            #print(elem, elem.coefficient_spin, i, coef)
            new_elem_base.coefficient_spin.append(elem.coefficient_spin[i])
            

    for dct in res:
        #rint('this_dict', dct)
        new_elem = deepcopy(new_elem_base)
        
        ln_g1 = len(dct['upper_g1'])
        ln_g2 = len(dct['upper_g2'])
        ln_g3 = len(dct['upper_g3'])
        ln_c2 = len(dct['upper_c2'])
        ln_c3 = len(dct['upper_c3'])

        for up, low in zip(dct['upper_g1'], dct['lower_g1']):

            if type_dict[low[0]] == 'o':
                sortcoef = [low[0][0], up[0][0]]

                if (low[0][1] != up[0][1]):
                    print('rozne spiny w n?')
                    sys.exit(0)
                sortcoef.sort()
                if sortcoef not in new_elem.delta:
                    new_elem.delta.append(sortcoef)

            else:
                
                new_elem.coefficient.append(DENSN)
                sortcoef = [low[0][0], up[0][0]]
            
                if (low[0][1] != up[0][1]):
                    print('rozne spiny w n?')
                    sys.exit(0)
                sortcoef.sort()
                new_elem.coefficient_idx.append([sortcoef[0]])
                if sortcoef not in new_elem.delta:
                    new_elem.delta.append(sortcoef)
                new_elem.coefficient_spin.append([])


        for up, low, name in zip(dct['upper_g2'], dct['lower_g2'], dct['name_g2']):


            sl_list =[low[0][1], low[1][1]]
            su_list =[up[0][1], up[1][1]]

            
            new_elem.coefficient.append(name)
            new_elem.coefficient_idx.append([low[0][0], low[1][0], up[0][0], up[1][0]])
            new_elem.coefficient_spin.append([low[0][1], low[1][1], up[0][1], up[1][1]])

        for up, low, name in zip(dct['upper_g3'], dct['lower_g3'], dct['name_g3']):


            sl_list =[low[0][1], low[1][1], low[2][1]]
            su_list =[up[0][1], up[1][1], up[2][1]]

            
            new_elem.coefficient.append(name)
            new_elem.coefficient_idx.append([low[0][0], low[1][0], low[2][0], up[0][0], up[1][0], up[2][0]])
            new_elem.coefficient_spin.append([low[0][1], low[1][1], low[2][1], up[0][1], up[1][1], up[2][1]])


        for up, low in zip(dct['upper_c2'], dct['lower_c2']):

            
            sl_list =[low[0][1], low[1][1]]
            su_list =[up[0][1], up[1][1]]

            
            new_elem.coefficient.append(CL2)
            new_elem.coefficient_idx.append([low[0][0], low[1][0], up[0][0], up[1][0]])
            new_elem.coefficient_spin.append([low[0][1], low[1][1], up[0][1], up[1][1]])


        for up, low in zip(dct['upper_c3'], dct['lower_c3']):


            sl_list =[low[0][1], low[1][1], low[2][1]]
            su_list =[up[0][1], up[1][1], up[2][1]]

            
            new_elem.coefficient.append(CL3)
            new_elem.coefficient_idx.append([low[0][0], low[1][0], low[2][0], up[0][0], up[1][0], up[2][0]])
            new_elem.coefficient_spin.append([low[0][1], low[1][1], low[2][1], up[0][1], up[1][1], up[2][1]])



        new_elem.num_factor*= dct['numf']
        #print(new_elem)
        result = result + arithmetic_string(new_elem)

    if len(result) > 0:
        return result
    else:
        return arithmetic_string()

def remove_gx(elem, k):

    #    for i, coef in enumerate(elem.coefficient):
    #        if i == k:
    elem.coefficient.pop(k)
    elem.coefficient_idx.pop(k)
    elem.coefficient_spin.pop(k)
    return elem

def cumulant3_to_gamma3(res, type_dict):

    result = []
    all_key_list = ['upper_g1', 'lower_g1', 'upper_g2', 'lower_g2', 'upper_g3', 'lower_g3', 'upper_c2', 'lower_c2', 'upper_c3', 'lower_c3', 'numf']

    for elem in res:

        ln = len(elem['upper_c3'])

        if ln == 1:

            upper  = elem['upper_c3'][0]
            lower  = elem['lower_c3'][0]

            res_1111 = partition_indices_1111(upper, lower, -1.0)
            res_12 = partition_indices_12_13(upper, lower, -1.0)

            res_1111 = remove_zeros(res_1111, type_dict)
            res_12 = remove_zeros(res_12, type_dict)

            new_elem = deepcopy(elem)
            new_elem['upper_c3'] = []
            new_elem['lower_c3'] = []
            new_elem['upper_g3'] = [upper]
            new_elem['lower_g3'] = [lower]
            result.append(new_elem)

            for dct_1111 in res_1111:
                new_elem = deepcopy(elem)
                new_elem['upper_c3'] = []
                new_elem['lower_c3'] = []
                new_elem['upper_g1'].extend(dct_1111['upper_g1'])
                new_elem['lower_g1'].extend(dct_1111['lower_g1'])
                new_elem['numf'] *= dct_1111['numf']
                result.append(new_elem)

            for dct_12 in res_12:
            
                new_elem = deepcopy(elem)
                new_elem['upper_c3'] = []
                new_elem['lower_c3'] = []
                new_elem['upper_g1'].extend(dct_12['upper_g1'])
                new_elem['lower_g1'].extend(dct_12['lower_g1'])
                new_elem['upper_c2'].extend(dct_12['upper_c2'])
                new_elem['lower_c2'].extend(dct_12['lower_c2'])

                new_elem['numf'] *= dct_12['numf']
                result.append(new_elem)

    result = remove_zeros(result, type_dict)
    
    return result

def cumulant2_to_gamma2(res, type_dict):

    result = []
    all_key_list = ['upper_g1', 'lower_g1', 'upper_g2', 'lower_g2', 'upper_g3', 'lower_g3', 'upper_c2', 'lower_c2', 'upper_c3', 'lower_c3', 'numf']
    
    for elem in res:
        
        ln = len(elem['upper_c2'])

        if ln == 1: # mamy tylko jedna kumulante
            
            mappings = [
                ('upper_g2', 'lower_g2', 0),
                ('upper_g1', 'lower_g1', 1),
                ('upper_g1', 'lower_g1', 2)
            ]
            # print('jedna kumulanta')
            c2_dicts = c2_to_g2(elem['upper_c2'][0], elem['lower_c2'][0])

            for upper_key, lower_key, dict_idx in mappings:
                new_elem = deepcopy(elem)
                new_elem['upper_c2'] = []
                new_elem['lower_c2'] = []
                new_elem[upper_key].extend(c2_dicts[dict_idx][upper_key])
                new_elem[lower_key].extend(c2_dicts[dict_idx][lower_key])
                new_elem['numf'] *= c2_dicts[dict_idx]['numf']
                result.append(new_elem)
                # print('dodaje', new_elem)
               

        elif ln ==2:

            mappings = [
                ['upper_g2', 'lower_g2'],
                ['upper_g1', 'lower_g1'],
                ['upper_g1', 'lower_g1']
            ]        

            # print('dwie kumulanty', elem)
            c2_dicts1 = c2_to_g2(elem['upper_c2'][0], elem['lower_c2'][0])
            c2_dicts2 = c2_to_g2(elem['upper_c2'][1], elem['lower_c2'][1])
            

            for idx1, idx2 in product(range(3), range(3)):
                new_elem = deepcopy(elem)
                new_elem['upper_c2'] = []
                new_elem['lower_c2'] = []

                upper_key, lower_key = mappings[idx1]
                new_elem[upper_key].extend(c2_dicts1[idx1][upper_key])
                new_elem[lower_key].extend(c2_dicts1[idx1][lower_key])
                

                upper_key, lower_key = mappings[idx2]
                new_elem[upper_key].extend(c2_dicts2[idx2][upper_key])
                new_elem[lower_key].extend(c2_dicts2[idx2][lower_key])
                
                # Mnożymy numf z obu słowników
                new_elem['numf'] *= c2_dicts1[idx1]['numf'] * c2_dicts2[idx2]['numf']
                
                result.append(new_elem)
                
        else:
            result.append(elem)


    result = remove_zeros(result, type_dict)

    return result

            
def c2_to_g2(upper, lower):
    all_key_list = ['upper_g1', 'lower_g1', 'upper_g2', 'lower_g2', 'upper_g3', 'lower_g3', 'upper_c2', 'lower_c2', 'upper_c3', 'lower_c3', 'numf']
    
    base_dict = {key: [] for key in all_key_list}
    
    cases = [
        {
            'upper_g2': [upper],
            'lower_g2': [lower],
            'numf': 1.0
        },
        {
            'upper_g1': [[upper[0]], [upper[1]]],
            'lower_g1': [[lower[0]], [lower[1]]],
            'numf': -1.0
        },
        {
            'upper_g1': [[upper[0]], [upper[1]]],
            'lower_g1': [[lower[1]], [lower[0]]],
            'numf': 1.0
        }
    ]
    
    result = []
    for case in cases:
        dict_case = deepcopy(base_dict)
        dict_case.update(case)
        result.append(dict_case)
        
    return result



def remove_zeros_g1(res, type_dict):

    result = []

#    print('spin_dict', spin_dict)
#    print('type_dict', type_dict)
#    print()
    for elem in res:
     #   print()
      #  print(elem)
        add = True
        for i in range(0, len(elem['upper_g1'])):
            
            i1 = elem['upper_g1'][i][0]
            i2 = elem['lower_g1'][i][0]

            #if spin_dict[i1] != spin_dict[i2]:
            if i1[1] != i2[1]:   
#                print('nie dodaje rozne spiny', i1, i2, spin_dict[i1], spin_dict[i2])
                
                add = False
            if type_dict[i1] != type_dict[i2]:
 #               print('nie dodaje rozne typy', i1, i2)
                add = False
        if add:
  #          print('dodaje')
            result.append(elem)
   # print()
    #print('bylo', len(res), 'jest', len(result))
    return result


def remove_zeros(res, type_dict):

    res = remove_zeros_g1(res, type_dict)
    res = remove_zeros_c2c3g2g3(res, type_dict)
    res = remove_zeros_c2c3(res, type_dict)

    return res
    
def remove_zeros_c2c3g2g3(res, type_dict):

    result = []
#    print('spin_dict', spin_dict)
#    print('type_dict', type_dict)
    
    tensor_types = ['c2', 'g2', 'c3', 'g3']
    
    for elem in res:
        add = True
        
        for tensor_type in tensor_types:
            upper_key = f'upper_{tensor_type}'
            lower_key = f'lower_{tensor_type}'
            
            if upper_key in elem and lower_key in elem:
                for i in range(len(elem[upper_key])):
                    n_plus = 0
                    
                    for idx in elem[upper_key][i]:
                        if idx[1] == '+':
                            n_plus += 1
                            
                    for idx in elem[lower_key][i]:
                        if idx[1] == '+':
                            n_plus += 1
                            
                    if n_plus % 2 != 0:
#                        print(f'nie dodaje nplus%2 dla {tensor_type}')
                        add = False
                        break
                
                if not add:
                    break
        
        if add:
 #           print('dodaje')
            result.append(elem)
            
    return result


def remove_zeros_c2c3(res, type_dict):

    result = []
    
    tensor_types = ['c2', 'c3']
    
    for elem in res:
        add = True

        for tensor_type in tensor_types:
            upper_key = f'upper_{tensor_type}'
            lower_key = f'lower_{tensor_type}'
            
            for i in range(len(elem[upper_key])):
                n_plus = 0
                    
                for idx in elem[upper_key][i]:
                    if type_dict[idx] == 'o' or type_dict[idx] == 'v':
                        add = False
                        break
                for idx in elem[lower_key][i]:
                    if type_dict[idx] == 'o' or type_dict[idx] == 'v':
                        add = False
                        break
                
            if not add:
                break
        
        if add:
            result.append(elem)
            
    return result

    
def nswaps_lists(lst1, lst2):
    index_map = {v: i for i, v in enumerate(lst1)}
    
    lst2_indices = [index_map[v] for v in lst2]

    n_swaps = sum(1 for i in range(len(lst2_indices)) for j in range(i + 1, len(lst2_indices)) if lst2_indices[i] > lst2_indices[j])
    
    return n_swaps


def partition_indices_2_2(upper, lower, sign):
    results = []
    all_key_list = ['upper_g1', 'lower_g1', 'upper_g2', 'lower_g2', 'upper_g3', 'lower_g3', 'upper_c2', 'lower_c2', 'upper_c3', 'lower_c3', 'numf']
    
    n = len(upper)
    if n != 4:
        raise ValueError("Partition 2:2 only works for lists of length 4.")
    

    upper_combinations = [[0,1,2,3], [0,2,1,3], [0,3,1,2]]
    lower_combinations = [[0,1,2,3], [0,2,1,3], [0,3,1,2], [2,3, 0, 1], [1,3, 0,2], [1,2,0,3]]
    
    
    for upper_idx_pair, lower_idx_pair in product(upper_combinations,  lower_combinations):

        c2_1_upper = [upper[i] for i in upper_idx_pair[0:2]]
        c2_2_upper = [upper[i] for i in  upper_idx_pair[2:4]]
        
        c2_1_lower = [lower[i] for i in lower_idx_pair[0:2]]
        remaining_lower_idx = [i for i in range(n) if i not in lower_idx_pair]
        c2_2_lower = [lower[i] for i in lower_idx_pair[2:4]]
        
        n_swaps = (sum(upper_idx_pair) + sum(lower_idx_pair)) % 2
        numf = (-1)**n_swaps*sign
        
        mini_dict = {key: [] for key in all_key_list}
        mini_dict['upper_c2'].append(c2_1_upper)
        mini_dict['upper_c2'].append(c2_2_upper)
        mini_dict['lower_c2'].append(c2_1_lower)
        mini_dict['lower_c2'].append(c2_2_lower)
        mini_dict['numf'] = numf
        
        results.append(mini_dict)
            
    return results

def partition_indices_1111(upper, lower, sign):
    results = []
    all_key_list = ['upper_g1', 'lower_g1', 'upper_g2', 'lower_g2', 'upper_g3', 'lower_g3', 'upper_c2', 'lower_c2', 'upper_c3', 'lower_c3', 'numf']
    n = len(upper)
    
    if n < 3:
        raise ValueError("Partition 1:1:1:1 needs at least 3 elements.")
    
    upper_indices = list(range(n))
  #  print(upper_indices)
    for perm_upper_idx in permutations(upper_indices, n):
 #       print(perm_upper_idx)
        n_swaps = nswaps_lists(deepcopy(upper_indices), list(perm_upper_idx)) % 2
        numf = (-1)**n_swaps*sign

        mini_dict = {}
        for key in all_key_list:
            mini_dict[key] = []
        mini_dict['numf'] = numf

        for i in range(0, n):
        
            mini_dict['upper_g1'].append([upper[perm_upper_idx[i]]])
            mini_dict['lower_g1'].append([lower[i]])

        results.append(mini_dict)
#        print(mini_dict)

    # sys.exit(0)
#    print('la---')
 #   for x in results:
  #     print(x)
  #  print('la')
    return results
        
        

def partition_indices_12_13(upper, lower, sign):
    results = []
    all_key_list = ['upper_g1', 'lower_g1', 'upper_g2', 'lower_g2', 'upper_g3', 'lower_g3', 'upper_c2', 'lower_c2' , 'upper_c3', 'lower_c3', 'numf']
    n = len(upper)

    for i, j in product(range(n), range(n)):  
        g1 = [upper[i]]
        c2 = upper[:i] + upper[i+1:] 

        g1_lower = [lower[j]]
        c2_lower = lower[:j] + lower[j+1:]  

        nswap_upper = [upper[i]]
        nswap_upper.extend(c2)


        nswap_lower = [lower[j]]
        nswap_lower.extend(c2_lower)


        nu = nswaps_lists(deepcopy(upper), nswap_upper)
        nl = nswaps_lists(deepcopy(lower), nswap_lower)

        n_swaps = (nu+nl)%2
        

        numf = (-1)**n_swaps * sign

        mini_dict = {}
        for key in all_key_list:
            mini_dict[key] = []
        mini_dict['numf'] = numf

        mini_dict['upper_g1'].append(g1)
        mini_dict['lower_g1'].append( g1_lower )

        if n == 3:
            mini_dict['upper_c2'].append(c2)
            mini_dict['lower_c2'].append(c2_lower )
        elif n == 4:
            mini_dict['upper_c3'].append(c2)
            mini_dict['lower_c3'].append(c2_lower )

        results.append(mini_dict)
        
    return results


def partition_indices_112(upper, lower, sign):
    results = []
    all_key_list = ['upper_g1', 'lower_g1', 'upper_g2', 'lower_g2', 'upper_g3', 'lower_g3', 'upper_c2', 'lower_c2', 'upper_c3', 'lower_c3', 'numf']
    
    n = len(upper)
    if n != 4:
        raise ValueError("Partition 1:1:2 only works for lists of length 4.")
    
    for i1, i2 in combinations(range(n), 2):

        g1_1 = [upper[i1]]
        g1_2 = [upper[i2]]

        c2 = [x for i, x in enumerate(upper) if i not in [i1, i2]]
        

        for j1, j2 in permutations(range(n), 2):

            g1_1_lower = [lower[j1]]
            g1_2_lower = [lower[j2]]

            c2_lower = [x for i, x in enumerate(lower) if i not in [j1, j2]]

            nswap_upper = [upper[i1], upper[i2]]
            nswap_upper.extend(c2)


            nswap_lower = [lower[j1], lower[j2]]
            nswap_lower.extend(c2_lower)


            nu = nswaps_lists(deepcopy(upper), nswap_upper)
            nl = nswaps_lists(deepcopy(lower), nswap_lower)

            n_swaps = (nu+nl)%2
            
            numf = (-1)**n_swaps *sign

            mini_dict = {key: [] for key in all_key_list}
            mini_dict['upper_g1'].append(g1_1)
            mini_dict['upper_g1'].append(g1_2)
            mini_dict['lower_g1'].append(g1_1_lower)
            mini_dict['lower_g1'].append(g1_2_lower)
            mini_dict['upper_c2'].append(c2)
            mini_dict['lower_c2'].append(c2_lower)
            mini_dict['numf'] = numf
            
            results.append(mini_dict)
    
    return results
                
def check_active(lst1, lst2, otype):
    
    for x in lst1:
        if otype[x] != 'a':
            return False
    for x in lst2:
        if otype[x] != 'a':
            return False
        
    return True


def remove_terms_with_non_diag_h(r):

    for x in r:
        for i, coef in enumerate(x.coefficient):
            if coef == BARENUCL_HAM:
                i1 = x.coefficient_idx[i][0]
                i2 = x.coefficient_idx[i][1]
                if x.orbital_type[i1] != x.orbital_type[i2]:
                    x.num_factor = 0.0
    r.cleanup()

#    return(r)

def remove_non_canonical_h(r):

    for x in r:
        for i, coef in enumerate(x.coefficient):
            if coef == BARENUCL_HAM:
                i1 = x.coefficient_idx[i][0]
                i2 = x.coefficient_idx[i][1]
#                print('yuuuu', x, x.orbital_type[i1], x.orbital_type[i2])
                if x.orbital_type[i1] == x.orbital_type[i2]:
                    if  x.orbital_type[i1] =='o' or x.orbital_type[i1] =='v':
                        if i1 != i2:
                            xcoef = sorted(x.coefficient_idx[i])
                            if xcoef not in x.delta:
                                x.delta.append(xcoef)
                            x.coefficient_idx[i][0] = xcoef[0]
                            x.coefficient_idx[i][1] = xcoef[0]

def remove_terms_with_delta_between_different_sets(r):

    for x in r:
        for lst in x.delta:
            
            i1 = lst[0]
            i2 = lst[1]
            if x.orbital_type[i1] != x.orbital_type[i2]:
                x.num_factor = 0.0
    r.cleanup()
    k = 0
#    return(r)
    
def divide_integrals_for_a0_a1(r):

    r0 = arithmetic_string()
    r1 = arithmetic_string()
    for x in r:
        num_of_a = 0
        if TWOEL_INT in x.coefficient:
            for i, coef in enumerate(x.coefficient):
                if coef == TWOEL_INT:
                    for idx in x.coefficient_idx[i]:
                        if x.orbital_type[idx] == 'a':
                            num_of_a += 1
            if num_of_a == 4:
                r0 += arithmetic_string(x)
            else:
                r1 += arithmetic_string(x)
        else:
           r0 += arithmetic_string(x) 
           r1 += arithmetic_string(x)
    return(r0, r1)

def remove_terms_with_virtual_in_any_DM(r):


#    print('przed usunieciem virtualnych DM', len(r))
    for x in r:
        for i, coef in enumerate(x.coefficient):
            if 'Gm' in coef or 'gm' in coef:
                for idx in x.coefficient_idx[i]:
                
                    if x.orbital_type[idx] == 'v':
                        x.num_factor = 0.0
    r.cleanup()
    #print('po usunieciu wirtualnych_DM', len(r))
    #k = 0
    #for x in r:
    #    print(k, x, x.orbital_type)
    #    k+=1
#    return(r)
    
        

def split_element(elem1):
    """
    Split element based on summation indices into multiple elements with different orbital types.    
    Args:
        elem1: Element object with summation and orbital_type attributes    
    Returns:
        list: List of new elements with all possible combinations of orbital types
    """
    n_indices = len(elem1.summation)
    
    combinations = list(product(['a', 'o'], repeat=n_indices))


    new_elements = arithmetic_string()
    for comb in combinations:
        new_el = deepcopy(elem1)
        
        for idx, orb_type in zip(elem1.summation, comb):
            new_el.orbital_type[idx] = orb_type
            
        new_elements += arithmetic_string(new_el)

    return new_elements    
    
def ph_rpa(cis = False):
    
    # -------------------------------------------------------------------------------------------------------
    # W scislym przypadku <0|[[dO, H],O+]|0> = <0|[dO,[H,O+]]|0>
    # <0|[dO, H, O+]|0>  = 1/2(<0|[[dO, H],O+]|0> + <0|[dO,[H,O+]]|0>)
    #
    # KL = <0|[[dO, H],O+]|0> = |jesli killer| =  <0|[dO, H],O+]0>
    # z tego <0|[[rs, H], q+p+]|0>  powinno wyjść to samo co z <0|[rs, H], q+p+|0>
    # KL = <0|[[rs, H], q+p+]|0>
    # KL_killer = <0|[rs, H], q+p+|0>  <--- tego nie bedziemy stosowac, bo ten powyzej jest leszpy - zwiazany
    #
    # KP = <0|[dO,[H,O+]]|0> = |jesli killer| = <0|dO,[H,O+]]0>
    # z tego  <0|[rs, [H, q+p+]|0>   powinno wyjść to samo co z   <0|rs, [H, p+q+]|0>  (w notacji komutatory na prawo)
    # z tego -<0|[[H, q+p+], rs|0>   powinno wyjść to samo co z   <0|rs, [H, q+p+]|0>  (w notacji komutatory na lewo)
    # KP = -<0|[[H, q+p+], rs|0>
    # KP_killer = <0|rs, [H, q+p+]|0>

    # wnioski
    # 1. Z komutatorow KL i KP powinno wyjsc to samo, jesli jest spelninony killer, czyli pq|0> = 0. Nigdzie pewnie tego nie zakladamy.
    # Zbadac czy takie wyrazy powstaja przed nazwaniem macierzy gestosci.
    
    # for PH
    # r+s = AArs
    # p+q = AApq

    # KL = <0|[[r+s, H], p+q]|0>

    # KL = evaluate(AArs, AApq)

#    KL = evaluate(h1cas, AApq, AArs).scale(-1.0) +  evaluate(h2cas, AApq, AArs).scale(-1.0)

    # KL = evaluate(AArs, h1cas, AApq) +  evaluate(AArs, h2cas, AApq)
    #A

    #
    # this is main code *********************************
    #
    if cis == False:
#        KL = evaluate(h1cas, AAqp, AArs).scale(-1.0) +  evaluate(h2cas, AAqp, AArs).scale(-1.0)
        LK = evaluate(h1cas, AApq, AArs).scale(1.0) +  evaluate(h2cas, AApq, AArs).scale(1.0)
        
        KL = evaluate(AArs, h1cas, AApq).scale(0.5) +  evaluate(AArs, h2cas, AApq).scale(0.5) +\
            evaluate(h1cas, AApq, AArs).scale(-0.5) +  evaluate(h2cas, AApq, AArs).scale(-0.5)

        # LK = evaluate(AApq, h1cas, AArs).scale(-0.5) +  evaluate(AApq, h2cas, AArs).scale(-0.5) +\
        #     evaluate(h1cas, AArs, AApq).scale(0.5) +  evaluate(h2cas, AArs, AApq).scale(0.5)

        # KL = evaluate(AAsr, h1cas, AApq).scale(0.5) +  evaluate(AAsr, h2cas, AApq).scale(0.5) +\
        #     evaluate(h1cas, AApq, AAsr).scale(-0.5) +  evaluate(h2cas, AApq, AAsr).scale(-0.5)

        # LK = evaluate(AApq, h1cas, AAsr).scale(-0.5) +  evaluate(AApq, h2cas, AAsr).scale(-0.5) +\
        #     evaluate(h1cas, AAsr, AApq).scale(0.5) +  evaluate(h2cas, AAsr, AApq).scale(0.5)
#        LK  = KL
        # KL = evaluate(AArs, h1cas, AAqp).scale(0.5) +  evaluate(AArs, h2cas, AAqp).scale(0.5) +\
        #     evaluate(h1cas, AAqp, AArs).scale(-0.5) +  evaluate(h2cas, AAqp, AArs).scale(-0.5)

    elif cis ==True:
        #
        # this is for the CIS code NEVPTS *********************************
        #
        KL = evaluate(h1cas, AApq) +  evaluate(h2cas, AApq)
    
    # KL = evaluate(AArs, h1cas, AAqp).scale(0.5) +  evaluate(AArs, h2cas, AAqp).scale(0.5) +\
    #     evaluate(h1cas, AAqp, AArs).scale(-0.5) +  evaluate(h2cas, AAqp, AArs).scale(-0.5)

    #B
    # KL = evaluate(h1cas, AApq, AArs).scale(-1.0) +  evaluate(h2cas, AApq, AArs).scale(-1.0)

    # KL = evaluate(AArs, h1cas, AApq).scale(0.5) +  evaluate(AArs, h2cas, AApq).scale(0.5) +\
    #     evaluate(h1cas, AApq, AArs).scale(-0.5) +  evaluate(h2cas, AApq, AArs).scale(-0.5)

    # KL = evaluate(h1cas, AAqp, AArs).scale(-1.0) +  evaluate(h2cas, AAqp, AArs).scale(-1.0)

    # KL = evaluate(h1cas, AAqpqp, AArsrs).scale(-1.0) +  evaluate(h2cas, AAqpqp, AArsrs).scale(-1.0)

    # KL = evaluate(AArs, h1cas, AAqp) +  evaluate(AArs, h2cas, AAqp)



    # KL = evaluate(AArs, h1cas, AAqp) +  evaluate(AArs, h2cas, AAqp) +\
    #     evaluate(h1cas, AAqp, AArs) +  evaluate(h2cas, AAqp, AArs)
    
    # KL = evaluate(AArs,  h1cas, AAqp).scale(0.5) +  evaluate(AArs, h2cas, AAqp).scale(0.5) +\
    #     evaluate(h1cas, AAqp, AArs).scale(-0.5) +  evaluate(h2cas, AAqp, AArs).scale(-0.5)

    # KL = evaluate(h1cas, AArp) + evaluate(h2cas, AArp)

    # KL = evaluate(AArs, h1cas, AAqp).scale(0.5) +  evaluate(AArs, h2cas, AAqp).scale(0.5) +\
    #     evaluate(h1cas, AAqp, AArs).scale(-0.5) +  evaluate(h2cas, AAqp, AArs).scale(-0.5)


    #B_rspq*
    # KL = evaluate(AArsg, h1cas, AAqp).scale(0.5) +  evaluate(AArsg, h2cas, AAqp).scale(0.5) +\
    #     evaluate(h1cas, AAqp, AArsg).scale(-0.5) +  evaluate(h2cas, AAqp, AArsg).scale(-0.5)

    # A_rspa z pq a nie qp
    # KL = evaluate(AArs, h1cas, AApq).scale(0.5) +  evaluate(AArs, h2cas, AApq).scale(0.5) +\
    #     evaluate(h1cas, AApq, AArs).scale(-0.5) +  evaluate(h2cas, AApq, AArs).scale(-0.5)


    print('wynik komutatorow KL')

    k = 0
    for x in KL:
        print(k, x)
        k += 1
    print('')


    if cis ==True:
        KL = KL.fromleft(AAsr)    

        print( 'wynik po fromleft')
        for x in KL:
            print(k, x)
            k += 1
        print('')
        

    print('Perform Wick')
    res = []
    for x in KL:
        result_list = x.wick_ca()
        res = res + result_list

    print('Perform Wick')
    resLK = []
    for x in LK:
        result_list = x.wick_ca()
        resLK = resLK + result_list

    
    print('Result po WICK:')    
    for x in res:
        x.exec_delta()
        if (x.num_factor) != 0.0:
           print(x)

    print('Result po WICK:')    
    for x in resLK:
        x.exec_delta()
        if (x.num_factor) != 0.0:
           print(x)

    print()
    print('Result po RENAME DENSITY:')
    print()
    rres = []

    for x in res:
        x.rename_as_density()
        if (x.num_factor) != 0.0:

            print(x)
    for x in resLK:
        x.rename_as_density()
        if (x.num_factor) != 0.0:

            print(x)

    print('')

    #---------------------------
    res2 = cas_to_ugg(res)
    res2LK = cas_to_ugg(resLK)

    rsimp = simplify(res2, cas=True)
    rsimpLK = simplify(res2LK, cas=True)
    k = 0
    print('po simp')
    for x in rsimp:
       k += 1

       print(x)
    print()

    print('po simp LK')
    for x in rsimpLK:
       k += 1
       print(x)
    print()

    klkl =rsimp + rsimpLK
    rs = simplify(klkl, cas=True)

    print('po simp klkl')
    for x in rs:
       k += 1

       print(x)
    print()

    
    sys.exit(0)
#----------------------------ten kod znajdzie calki antysym---------------
    # rsimpa = find_antysym_int(rsimp)
    # k = 0
    # print()
    # print('po simp antysym')
    # print()
    # for x in rsimpa:
    #     k += 1
    #     # print("&", x, "\\\\")
    #     print(x)
    # sys.exit(0)



    # sys.exit(0)

    spin_dict1={}
    spin_dict1['r'] = '+'
    spin_dict1['s'] = '+'
    spin_dict1['p'] = '+'
    spin_dict1['q'] = '+'
    numf1 = 1.0

    spin_dict2={}
    spin_dict2['r'] = '+'
    spin_dict2['s'] = '+'
    spin_dict2['p'] = '-'
    spin_dict2['q'] = '-'
    numf2 = 1.0


    rsimp2 = deepcopy(rsimp)
    res3_aaaa = add_spin_driver(spin_dict1, numf1, rsimp)
    res3_aabb = add_spin_driver(spin_dict2, numf2, rsimp2)

    res3 = res3_aaaa + res3_aabb
    print('')
    print('last result po spin')
    k = 0
    res3.cleanup()
    for x in res3:
        print(k, x)
        k += 1

#    sys.exit(0)
    # print('')
    # print('last result po spin')
    # for x in res3_aabb:
    #     print(x)

    if cis == False:
        rgams = simplify_for_mult_2(res3)
    elif cis == True:
        rgams = simplify_for_mult_234(res3)

    

    print()
    print('po simp and rename gamma', len(rgams))
    print()
    print(type(rgams))
    for x in rgams:
        print(x)
    print()

    rgams = simplify_final_touch(rgams)
    print()
    print('po final touch', len(rgams))
    print()
    print(type(rgams))
    for x in rgams:
        print(x, x.coefficient)
    print()

    rgam_original = add_spin_to_gamma(rgams)
    
    print('len rgam_original', len(rgam_original))

    A0_blocks = {}
    A1_blocks = {}

    for key in Ablock_dict:

        print('this_key', key)
        Ablock_idx = key        
        #myfixed = ['p','q','r','s']
        myfixed = ['r','s', 'p', 'q']

        rgam = deepcopy(rgam_original)
        
        rgamSplit = arithmetic_string()
        for x in rgam:
            for fx in range(0, 4):
                x.orbital_type[myfixed[fx]] = Ablock_dict[Ablock_idx][fx]

            if len(x.summation) == 0:
                rgamSplit += arithmetic_string(x)
            else:
                new_elements = split_element(x)
                rgamSplit += new_elements

            
        print('len1', len(rgamSplit))
        remove_terms_with_virtual_in_any_DM(rgamSplit)
        

        print('len po virt', len(rgamSplit))
        remove_terms_with_delta_between_different_sets(rgamSplit)
        print('len po delta', len(rgamSplit))
        remove_terms_with_non_diag_h(rgamSplit)
        print('len po hdiag', len(rgamSplit))
        for x in rgamSplit:
            x.coefficient[:] = [coef if coef != DENS1 else DENSN for coef in x.coefficient]

        # k = 0
        # for x in rgamSplit:
        #     print(k, x, x.orbital_type)
        #     k+=1
        # print()

        rgamSplit0, rgamSplit1 = divide_integrals_for_a0_a1(rgamSplit)
        print('len po split0 split1', len(rgamSplit0), len(rgamSplit1))
        k = 0
        for x in rgamSplit1:
            print(k, x, x.orbital_type)
            k+=1
        print()


        A0_blocks[key] = deepcopy(rgamSplit0)
        A1_blocks[key] = deepcopy(rgamSplit1)



    A0_blocks_simp = {}
    print('A0 blocks')
    for key in A0_blocks:
        print()
        print('teraz robie klucz', key, len(A0_blocks[key]))
        if len(A0_blocks[key]) > 0:

            ress = decompose_gm2(A0_blocks[key], 0)


            for x in ress:
                x.exec_delta_cas_fixed()
            rsimp = simplify(ress, cas=True)
            print('ostateczny wynik', len(ress), len(rsimp), key)
            A0_blocks_simp[key] = deepcopy(rsimp)
            k = 0
            for xx in rsimp:
                print(k, xx, xx.orbital_type)
            k+=1

            print()

    print()
    print('A1 blocks')
    A1_blocks_simp = {}

    for key in A1_blocks:
        print()
        print('teraz robie klucz A1', key, len(A1_blocks[key]))

        if len(A1_blocks[key]) > 0:
            ress = decompose_gm2(A1_blocks[key], 0)

            for x in ress:
                x.exec_delta_cas_fixed()
            rsimp = simplify(ress, cas=True)
            print('ostateczny wynik simp A1', len(ress), len(rsimp), key)
            A1_blocks_simp[key] = deepcopy(rsimp)

            k = 0
            for xx in rsimp:
                print(k, xx, xx.orbital_type)
                k+=1
                
            print()


    A0_blocks = A0_blocks_simp
    A1_blocks = A1_blocks_simp

    for key in A0_blocks:
        if (len(A0_blocks[key]))> 0:
            print(key, len(A0_blocks[key]))
            k = 0
            for x in A0_blocks[key]:
                seen = set()
                unique_delta = []
                for item in x.delta:
                    item_str = str(item)
                    if item_str not in seen:
                        seen.add(item_str)
                        unique_delta.append(item)
                x.delta = unique_delta
                x.exec_delta_cas_fixed()

            A0_blocks[key] = simplify(A0_blocks[key], cas=True)
                
       
    for key in A1_blocks:
        if (len(A1_blocks[key]))> 0:
            print(key, len(A1_blocks[key]))
            k = 0
            for x in A1_blocks[key]:
                seen = set()
                unique_delta = []
                for item in x.delta:
                    item_str = str(item)
                    if item_str not in seen:
                        seen.add(item_str)
                        unique_delta.append(item)
                x.delta = unique_delta
                x.exec_delta_cas_fixed()
                    
            A1_blocks[key] = simplify(A1_blocks[key], cas=True)
                
    print('final print minipluszek')
    print()
    print('A0*******************************************************************************************************************************************')
    print('')
    print()
    for key in A0_blocks:
        if (len(A0_blocks[key]))> 0:
            print(f"+{'-'*30}+")
            print(f"| a0 block {key} {len(A0_blocks[key]):^10} |")
            print(f"+{'-'*30}+")
            print_struct(A0_blocks[key])
            print('================================================================================')

            # for x in A0_blocks[key]:
            #     print(x)
            print('================================================================================')
            print('================================================================================')
            print()

    print()
    print('A1*******************************************************************************************************************************************')
    print()
    for key in A1_blocks:
        if (len(A1_blocks[key]))> 0:
            print(f"+{'-'*30}+")
            print(f"| a1 block {key} {len(A1_blocks[key]):^10} |")
            print(f"+{'-'*30}+")

            print_struct(A1_blocks[key])
            print('================================================================================')
            # for x in A1_blocks[key]:
            #     print(x)
            print('================================================================================')
            print('================================================================================')
            print()
            
    sys.exit(0)

    

def add_spin_to_gamma(r):

    for x in range(0, len(r)):
        for i, y, in enumerate(r[x].coefficient):
            coef = r[x].coefficient[i]
            if 'Gm' in coef:
                print(coef)
                if coef == DENS2PM:
                    r[x].coefficient_spin[i] = ['+', '-', '+', '-']
                elif coef == DENS2P:
                    r[x].coefficient_spin[i] = ['+', '+', '+', '+']
                elif coef == DENS2M:
                    r[x].coefficient_spin[i] = ['-', '-', '-', '-']
                elif coef == DENS3PM:
                    r[x].coefficient_spin[i] = ['+', '+', '-', '+', '+', '-']
                elif coef == DENS3P:
                    r[x].coefficient_spin[i] = ['+', '+', '+', '+', '+', '+']
                elif coef == DENS3M:
                    r[x].coefficient_spin[i] = ['-', '-', '-', '-', '-', '-']
                elif coef == DENS4PM:
                    r[x].coefficient_spin[i] = ['+', '+', '-', '-', '+', '+', '-', '-']
                elif coef == DENS4P:
                    r[x].coefficient_spin[i] = ['+', '+', '+', '+', '+', '+', '+', '+']
                elif coef == DENS4M:
                    r[x].coefficient_spin[i] = ['-', '-', '-', '-', '-', '-', '-', '-']
                elif coef == DENS4PPM:
                    r[x].coefficient_spin[i] = ['+', '+', '+', '-', '+', '+', '+', '-']
                else:
                    print('other gm', r[x].coefficient[i])
                    sys.exit(0)
            else:
                r[x].coefficient_spin[i] = []

    return r



def decompose_rdm2(r):

    res_new = arithmetic_string()
    res_old = arithmetic_string()
    
    for x in range(0, len(r)):
        elem = r[x]
        print('rzzzz', elem.coefficient)
        idx = [i for i, x in enumerate(elem.coefficient) if x.startswith('Gm2')]
        print(idx)
        if len(idx)  > 0:
            res = decompose_rdm2_mini(elem, idx[0])
            
            res_new += res
        else:
            res_old += arithmetic_string(elem)

    return res_old, res_new

def decompose_rdm2_mini(elem, i):


    print('elem', elem, i)
    coef_idx = elem.coefficient_idx[i]
    coef_spin = elem.coefficient_spin[i]
    
    g1_idx1 = [coef_idx[0], coef_idx[2]]
    g1_idx2 = [coef_idx[1], coef_idx[3]]

    if coef_spin[0] == coef_spin[2] and coef_spin[1] == coef_spin[3]:
        g1_add = True
    else:
        g1_add = False

    g2_idx1 = [coef_idx[0], coef_idx[3]]
    g2_idx2 = [coef_idx[1], coef_idx[2]]

    if coef_spin[0] == coef_spin[3] and coef_spin[1] == coef_spin[2]:
        g2_add = True
    else:
        g2_add = False

    res = arithmetic_string(elem)
    if g1_add: 
        elem1 = deepcopy(elem)

        elem1.coefficient[i] = 'Gm1'
        elem1.coefficient_idx[i] = g1_idx2
        elem1.coefficient.insert(i, 'Gm1')
        elem1.coefficient_idx.insert(i, g1_idx1)
        elem1.delta.append(g1_idx2)
        elem1.delta.append(g1_idx1)
        res = arithmetic_string(elem1)
    if g2_add:
        elem2 = deepcopy(elem)
        elem2.coefficient[i] = 'Gm1'
        elem2.coefficient_idx[i] = g2_idx2
        elem2.coefficient.insert(i, 'Gm1')
        elem2.coefficient_idx.insert(i, g2_idx1)
        elem2.delta.append(g2_idx2)
        elem2.delta.append(g2_idx1)
        
        elem2.num_factor *= -1.0
        if g1_add:
            res += arithmetic_string(elem2)
        else:
            res = arithmetic_string(elem2)
        

    print()
    print(elem)
    print(res)

    return res

def zero_dens_matrix_onedet(res):

    #----------------------------------kod do zerowania macierzy gestosci dla jednowyznacznikowych-----------------

    resal = arithmetic_string()
    dd = 0
    for x in res:
        dd += 1
 #       print('sprawdzam', x, x.operator_idx)
        if len(x.operator_idx) > 0:

            if (len(list(set(x.operator_idx).intersection(virtual)))>0):
                x.num_factor = 0
            resal.append(x)
        else:
            resal.append(x)
            
            
            # opi0 = x.operator_idx[0]
            # opi1 = x.operator_idx[-1]
            # opt0 = x.operator_type[0]
            # opt1 = x.operator_type[-1]
#            print(opi0, opi1, opt0, opt1)
            
        #     if opi0 in virtual and opt0 == CRE:
        #         x.num_factor = 0
        #         continue
        #     if opi1 in virtual and opt1 == ANI:
        #         x.num_factor = 0
        #         continue
        #     if opi0 in occupied and opt0 == ANI:
        #         x.num_factor = 0
        #         continue
        #     if opi1 in occupied and opt1 == CRE:
        #         x.num_factor = 0
        #         continue
        #     resal.append(x)
        # else:
        #     resal.append(x)
                
            
            # if x.operator_idx[0] in virtual:
            #     print('virtual')
            #     if x.operator_type[0] == CRE:
            #         print('i kreator')
            #         x.num_factor = 0
            # if x.operator_idx[-1] in virtual:
            #     if x.operator_type[-1] == ANI:
            #         print('i anihil')
            #         x.num_factor = 0
            # if x.operator_idx[0] in occupied:
            #     print('occupaj')
            #     if x.operator_type[0] == ANI:
            #         print('i anihil2')                                    
            #         ## print('taaaaak')
            #         x.num_factor = 0
            # if x.operator_idx[-1] in occupied:
            #     if x.operator_type[-1] == CRE:
            #         print('i kreator2')
            #         x.num_factor = 0
            # if len(set(x.operator_idx).intersection(virtual)) != 0:
            #     x.num_factor = 0.0
            #     print('jest zero bo   ?')
            #     x.operator_idx = []

        # # print('witam1')

        # if x.num_factor != 0:
        # #     d = 1
        # #     # print('aa')
        # #     # print('jest zero')
        # # else:
        #     # print('nie ma')
        #     # if x.num_factor != -100:
        #     # print('dodaje', x)
        #     resal = resal + arithmetic_string(x)
            # else:
            #     print('tak dodaje te dwa', res1, res2, x)
            #     resal = resal + arithmetic_string(res1)
            #     resal = resal + arithmetic_string(res2)
        # if dd == 1:
        #     sys.exit(0)

#    res = deepcopy(resal)
    print('lenreszero', len(resal))

    #----------------------------------koniec kodu do zerowania macierzy gestosci dla jednowyznacznikowych-----------------
    return resal


    
def pp_rpa_hf_det():
    
    # -----sprawdzenie wersji jednowyznacznikowej
    # - utworz wszystkie mozliwe pary creator, anihilatr, virt, occ


    ovlist = [['i', 'j']]#, ['c', 'k'], ['k', 'c'], ['k','l']]
    print(ovlist)
    # ovlist = [['k', 'l']]
    crealist = [['0', '0']]
    Axx = cas()
    Axx.operator_idx = ['c', 'd']
    Axx.operator_type =['+', '+']
    # Axx.coefficient = ['n']
    # Axx.coefficient_idx.append(['a', 'b'])
    # Axx.coefficient_spin.append(['+', '+'])
    # print(Axx)
    # sys.exit(0)

    Ayy = cas()
    Ayy.operator_idx =   ['l', 'k']
    Ayy.operator_type =['+', '+']

    KL = arithmetic_string()
    for v1 in range(0, len(ovlist)):
        for v2 in range(0, len(crealist)):
            Ars = cas()
            Ars.operator_idx = ovlist[v1]
            Ars.operator_type = crealist[v2]

            # XX = evaluate(h1cas, Axx, Ars).scale(-1.0) +  evaluate(h2cas, Axx, Ars).scale(-1.0)
            XX = arithmetic_string()
            # YY = arithmetic_string()
            # XX = evaluate(Ars, h1cas, Axx).scale(0.5) +  evaluate(Ars, h2cas, Axx).scale(0.5) +\
            #     evaluate(h1cas, Axx, Ars).scale(0.5) +  evaluate(h2cas, Axx, Ars).scale(0.5)
            
            # XX = evaluate(Ars, h1cas, Axx) +  evaluate(Ars, h2cas, Axx)

            YY = evaluate(h1cas, Ayy, Ars ).scale(-1.0) +  evaluate(h2cas, Ayy, Ars).scale(-1.0)
            # YY = evaluate(Ars, h1cas, Ayy ) +  evaluate(Ars, h2cas, Ayy)

            # YY = evaluate(Ars, h1cas, Ayy ).scale(0.5) +  evaluate(Ars, h2cas, Ayy).scale(0.5) + \
            #     evaluate(h1cas, Ayy, Ars ).scale(-0.5) +  evaluate(h2cas, Ayy, Ars).scale(-0.5)

            
             # YY = evaluate(AArs, h1cas, Ayy) +  evaluate(Ars, h2cas, AAqp)

            # YY =  evaluate(Ars, h2cas, Ayy).scale(-1.0)
            KL = KL + XX
            KL = KL + YY
            
            print('')
            print('wynik komutatorow XX')
            k = 0
            for x in XX:
                print(k, x)
                k += 1
            print('')
            print('')
            print('wynik komutatorow YY')
            k = 0
            for x in YY:
                print(k, x)
                k += 1
            print('')


    print('wynik komutatorow KL')
    k = 0
    for x in KL:
        print(k, x)
        k += 1
    print('')

    k = 0
    for x in KL:
        print(k, x)
        k += 1
    print('')


    print('Perform Wick')
    res = []
    for x in KL:
        result_list = x.wick_ca()
        res = res + result_list
    print('')

    print('Result po WICK:')    
    for x in res:
        x.exec_delta()
        print(x)
    print('')


    #----------------------------------kod do zerowania macierzy gestosci dla jednowyznacznikowych-----------------


    resal = arithmetic_string()
    dd = 0
    for x in res:
        dd += 1
        print('sprawdzam', x, x.operator_idx)
        if len(x.operator_idx) > 0:
            if x.operator_idx[0] in virtual:
                if x.operator_type[0] == CRE:
                    x.num_factor = 0
            if x.operator_idx[-1] in virtual:
                if x.operator_type[-1] == ANI:
                    x.num_factor = 0
            if x.operator_idx[0] in occupied:
                if x.operator_type[0] == ANI:
                    #print('taaaaak')
                    x.num_factor = 0
            if x.operator_idx[-1] in occupied:
                if x.operator_type[-1] == CRE:
                    x.num_factor = 0
            if len(set(x.operator_idx).intersection(virtual)) != 0:
                x.num_factor = 0.0
                print('jest zero')

            # if len(x.operator_idx) == 2:
                
            #     x.coefficient.append('n')
            #     x.coefficient_idx.append([x.operator_idx[0], x.operator_idx[1]])
            #     #print(x)
            #     x.new_delta(x.operator_idx[0], x.operator_idx[1])
                #print(x)                            
                x.operator_idx = []
                #print(x)
            # if len(x.operator_idx) == 4:
            #     res1 = deepcopy(x)
            #     res2 = deepcopy(x)
            #     # res1.operator_idx = [x.operator_idx[0], x.operator_idx[2]]
            #     # res1.operator_idx.append([x.operator_idx[1], x.operator_idx[3]])
            #     res1.new_delta(x.operator_idx[0], x.operator_idx[2])
            #     res1.new_delta(x.operator_idx[1], x.operator_idx[3])
            #     res1.num_factor *= -1.0
            #     res2.new_delta(x.operator_idx[0], x.operator_idx[3])
            #     res2.new_delta(x.operator_idx[1], x.operator_idx[2])
            #     res1.operator_idx = []
            #     res2.operator_idx = []
            #     x.num_factor = -100

        print('wtiam1')

        if x.num_factor == 0:
            print('aa')
            print('jest zero')
        else:
            print('nie ma')
            # if x.num_factor != -100:
            print('dodaje', x)
            resal = resal + arithmetic_string(x)
            # else:
            #     print('tak dodaje te dwa', res1, res2, x)
            #     resal = resal + arithmetic_string(res1)
            #     resal = resal + arithmetic_string(res2)
        # if dd == 1:
        #     sys.exit(0)

    res = deepcopy(resal)
    #----------------------------------koniec kodu do zerowania macierzy gestosci dla jednowyznacznikowych-----------------

    # print(x)
    # print('')
    print('')

    print('Result po WICK:')    
    for x in res:
        x.exec_delta()
        if (x.num_factor) != 0.0:
            print(x)
        print(x)
    print('')


    print('Result po RENAME DENSITY:')    
    for x in res:
        x.rename_as_density()
        if (x.num_factor) != 0.0:
            print(x)
    print('')
    print('')
    print('')

    print('Result temp:')    
    for x in res:
        if (x.num_factor) != 0.0:
            print(x)
        print(x)
    print('')
    
    res2 = cas_to_ugg(res)
    print('po tranformacji')
    for x in res2:
        print(x)
    
    rsimp = simplify(res2, cas=True)
    k = 0
    print('po simp')
    for x in rsimp:
        k += 1
        x.exec_delta()
        # print("&", x, "\\\\")
        print(x)

        
    res_abab = deepcopy(res)
    res_abba = deepcopy(res)
        
    spin_dict1={}
    spin_dict1['i'] = '+'
    spin_dict1['j'] = '-'
    spin_dict1['k'] = '+'
    spin_dict1['l'] = '-'
    numf1 = 1.0

    spin_dict2={}
    spin_dict2['i'] = '+'
    spin_dict2['j'] = '-'
    spin_dict2['k'] = '-'
    spin_dict2['l'] = '+'
    numf2 = -1.0

    res3_abab = add_spin_driver(spin_dict1, numf1, res_abab)
    print('')
    print('res_abab-----dup')
    
    rgam_abab = simplify_for_hf(res3_abab)

    # rgam_abab = dirac_to_coulomb(rgam_abab)
    # print('kla')
    # for x in rgam_abab:
    #     print(x)
    # sys.exit(0)

    res3_abba = add_spin_driver(spin_dict2, numf2, res_abba)
    rgam_abba = simplify_for_hf(res3_abba)

    print('')
    print('res_abab')
    for x in rgam_abab:
        print(x)

    print('')
    print('res_abba')
    for x in rgam_abba:
        print(x)
    print('--------------')

    res_all = rgam_abab+rgam_abba
    rsimp = simplify(res_all)
    print('')
    print('resimp all', len(rgam_abab), len(rgam_abba), len(res_all))
    for x in rsimp:
        print(x)

    sys.exit(0)


def simplify_for_mult_2(res3):

    res5 = simplify_1rdm_mult(res3)
    print('')
    print('')
    print('po simp rdm1', len(res5))
    k = 0
    for x in res5:
        print(k, x)
        k += 1

    res5a = simplify_2rdm_mult_act(res5)
    print('')
    print('')
    print('po simp rdm2', len(res5a))
    for x in res5a:
        if x.num_factor !=0:
            print(x)
#    sys.exit(0)
    res6 = remove_spin(res5a)
    print('')
    print('')
    print('po remove spin', len(res6))
    for x in res6:
        print(x)
#    sys.exit(0)
    res2 = cas_to_ugg(res6)
    print('po tranformacji')
    k = 0
    for x in res2:
        print(k, x)
        k+=1
#    sys.exit(0)
    
    rsimp = simplify(res2, cas=True)
    rsimp = dirac_to_coulomb(rsimp)
    k = 0
    print('po simp')
    for x in rsimp:
        x.exec_delta(general_cond=True)
        k += 1
        # print("&", x, "\\\\")
        print(x)

#    rgam = rename_and_rescale_gamma(rsimp)
    # print('rename_and_rescale_gamma')
    # for x in rgam:
    #     print(x)
    # return rgam
    return rsimp

def simplify_for_mult_234(res3):

    res5 = simplify_1rdm_mult(res3)
    print('')
    print('')
    print('po simp rdm1', len(res5))
    k = 0
    for x in res5:
        print(k, x)
        k += 1

    res5a = simplify_2rdm_mult_act(res5)
    print('')
    print('')
    print('po simp rdm2', len(res5a))
    for x in res5a:
        if x.num_factor !=0:
            print(x)


    res5a = simplify_3rdm_mult_act(res5a)
    res5a = simplify_4rdm_mult_act(res5a)

    print()
    print('resresres')
    print()
    res6a = []
    for x in res5a:
        #print(type(x.num_factor), x)
        if x.num_factor != 0:
            res6a.append(x)
            print('le------', x)
        #else:
        #    print('else', x)
        #sys.exit(0)
#    sys.exit(0)
    res5a = res6a
    k = 0



    print('po tranformacji')
    for x in res5a:
        # print()
        print(k, x)
        #x.exec_delta_fixed(['p', 't'])
        #x.exec_delta_fixed(['q', 'u'])
        # x.exec_delta_fixed(['r', 'v'])
        # x.exec_delta_fixed(['s', 'w'])
        print(k, x)
        print()
        k+=1


    res6 = remove_spin(res5a)
    print('')
    print('')
    print('po remove spin', len(res6))
    for x in res6:
        print(x)

    print('typeres, typeres', type(res6), type(res6[0]))

    res2 = cas_to_ugg(res6)
    print('typeres, typeres', type(res2), type(res2[0]))
    print('po tranformacji')
    k = 0
    for x in res2:
        print(k, x)
        k+=1
#    sys.exit(0)
    
    rsimp = simplify(res2, cas=True)
    rgam = dirac_to_coulomb(rsimp)
    k = 0
    print('po simp')
    for x in rgam:
        x.exec_delta(general_cond=True)
        k += 1
        # print("&", x, "\\\\")
        print(k, x)

    #rgam = rename_and_rescale_gamma(rsimp)
    # print('rename_and_rescale_gamma')
    # for x in rgam:
    #     print(x)
    return rgam


def simplify_for_hf(res3):
    res5 = simplify_2rdm_occ(res3)
    print('')
    print('')
    print('po simp rdm', len(res5))
    for x in res5:
        print(x)

    res6 = remove_spin(res5)
    print('')
    print('')
    print('po simp rdm', len(res6))
    for x in res5:
        print(x)

    res2 = cas_to_ugg(res6)
    print('po tranformacji')
    for x in res2:
        print(x)
#    sys.exit(0)
    
    rsimp = simplify(res2, cas=True)
    k = 0
    print('po simp')
    for x in rsimp:
        x.exec_delta()
        k += 1
        # print("&", x, "\\\\")
        print(x)

    rgam = rename_and_rescale_gamma(rsimp)
    print('rename_and_rescale_gamma')
    for x in rgam:
        print(x)
    return rgam

def add_spin_driver(spin_dict, numf, rsimp, onedet = False):

    start = time.time()

    for x in rsimp:
        print(x)
    print('')
    res = ugg_to_cas(rsimp)

    print('cascas')
    for x in res:
        print(x)
    print('srala')

    print('przedwogule', len(res))
    for k in range(0, len(res)):
        print(k+1, res[k])
    print('')
    end = time.time()
    print('przedwogule time', end - start)
    
    start = time.time()
    res3 = add_spin_from_list(res, spin_dict, numf)
    print('gowno')
    end = time.time()
    print('add spin_l time', end - start)
    
    print('po add spin', len(res3))
    for x in res3:
        print(x)
    print('koniec')
    print()
    print('po add spinlist', len(res3))
    for k in range(0, len(res3)):
        print(k+1, res3[k])
    print('')

    for x in res3:
        if (x.num_factor) != 0.0:
            print(x)

    start = time.time()
    print('llk1', len(res3))

    res3b = simplify_spin_initially(res3)
    print('po simp initially')
    for x in res3b:
        print(x)
    print('koniec')
    print()

    print('llk2', len(res3b))
    end = time.time()
    print('simp spin_init time', end - start)

#    res3b = deepcopy(res3)
    print('')
    print('')
    # print('przed add sum spin', len(res3b))    
    # for x in res3b:
    #     if (x.num_factor) != 0.0:
    #         print(x)

    start = time.time()
    res3a = add_spin_summ(res3b)
    end = time.time()
    print('add spin sum time', end - start)
    print('po add spin')
    for x in res3a:
        print(x)
    print('koniec')
    print()




    # print('')
    # print('')
    # print('po add sum spin', len(res3a))    
    # for x in res3a:
    #     if (x.num_factor) != 0.0:
    #         print(x)

    print('przedsimpspin', len(res3a))
    for k in range(0, len(res3a)):
        print(k+1, res3a[k])
    print('')

    start = time.time()
    if onedet:
        res4 = simplify_spin_onedet(res3a)
    else:
        res4 = simplify_spin(res3a)
    end = time.time()
    print('simpl spin time', end - start)

    print('')
    print('')
    print('po simp spin-rra', len(res4))
    for x in res4:
        if (x.num_factor) != 0.0:
            print(x)

#    sys.exit(0)
    return res4

def test_rpa():

    # e1 = cas()
    # e1.operator_idx = ['r', 's', 'q', 'p', 't','u']
    # e1.operator_type = ['0', '0', '+','+', '+','0']

    # e1 = cas()
    # e1.operator_idx = ['t', 'u', 'q', 'p', 'r','s']
    # e1.operator_type = ['+', '0', '+','+', '0','0']

    # e1 = cas()
    # e1.operator_idx = ['q', 'p', 'r', 's', 't','u']
    # e1.operator_type = ['+', '+', '0','0', '+','0']

    e1 = cas()
    e1.operator_idx = ['q', 's', 't', 'u', 'v','w']
    e1.operator_type = ['+', '0', '+','+', '0','0']

    print(e1)
    res  = e1.wick_ca()

    print('wynik')
    for x in res:
        print(x)
    sys.exit(0)

        

def erpa_ph():

    # epq = cas()
    # epq.operator_idx = ['k', 'l']
    # epq.operator_type = ['0', '0']
    # epq1 = cas()
    # epq1.operator_idx = ['i', 'a', 'j', 'b', 'c','d','l','k']
    # epq1.operator_type = ['+', '0', '+','0', '+','+','0','0']
    # epq1.wick_ca()
    # ress = evaluate(epq, epq1)

    # res = []
    # for x in ress:
    #     ress = x.wick_ca()
    #     res = res + ress
    # print('')

    # print('otowick')
    # for x in res:
    #     print(x)
    # sys.exit(0)

    # epq = cas()
    # epq.operator_idx = ['s', 'u', 'q', 'p']
    # epq.operator_type = ['0', '0', '+', '+']
    # res = epq.wick_ca()
    # for x in res:
    #     print(x)
    # sys.exit(0)
    # print('')    
    # print('[', h2cas, ' , ', AApq, ']')
    # print(AArs)
    # print(h1cas)
    # print(AApq)
    # print('')

    # l1 = evaluate(AArs, h1cas, AApq)
    # l1 = evaluate(AArs, h2cas, AApq)
    # l1 = evaluate(AArs, h2cas, AApq)
    # l1 = evaluate(h2cas, AApq, AArs).scale(-1.0)
    # l1 = evaluate(AArs, h2cas, AApq).scale(0.5) + evaluate(h2cas, AApq, AArs).scale(-0.5)
    #l1 = evaluate(h2cas, AApq, AArs).scale(-0.5)
    # l1 = evaluate(h2cas, AApq, AArs)

    # l1 = evaluate(AArs, h1cas, AApq).scale(0.5) + evaluate(h1cas, AApq, AArs).scale(-0.5) \
    #     + evaluate(AArs, h2cas, AApq).scale(0.5) + evaluate(h2cas, AApq, AArs).scale(-0.5)

    # l1 = evaluate(h1cas, AApq, AArs).scale(-1.0) + evaluate(h2cas, AApq, AArs).scale(-1.0)
    
    # l1 = evaluate(AArs, h1cas, AApq).scale(0.5) + evaluate(AArs, h2cas, AApq).scale(0.5) \ 
    #     + evaluate(h1cas, AApq, AArs).scale(-0.5) + evaluate(h2cas, AApq, AArs).scale(-0.5)

    #-----------------------------------------------PP-------------------------------------------------
    # l1-pp
    # l1 = evaluate(AArs_aa, AAqp_pp)

    # l1 = evaluate(h1cas, AAqp_pp, AArs_aa).scale(-1.0) + evaluate(h2cas, AAqp_pp, AArs_aa).scale(-1.0)

    # l1 = evaluate(AArs_aa, h1cas, AAqp_pp).scale(0.5) + evaluate(AArs_aa, h2cas, AAqp_pp).scale(0.5) \
    #     + evaluate(h1cas, AAqp_pp, AArs_aa).scale(-0.5) + evaluate(h2cas, AAqp_pp, AArs_aa).scale(-0.5)

    #-----------------------------------------------HH-------------------------------------------------
    # l1-hh
    # l1 = evaluate(AArs_aa, AAqp_pp)

    # l1 = evaluate(h1cas, AAqp_hh, AArs_cc).scale(-1.0) + evaluate(h2cas, AAqp_hh, AArs_cc).scale(-1.0)

    # l1 = evaluate(AArs_cc, h1cas, AAqp_hh).scale(0.5) + evaluate(AArs_cc, h2cas, AAqp_hh).scale(0.5) \
    #         + evaluate(h1cas, AAqp_hh, AArs_cc).scale(-0.5) + evaluate(h2cas, AAqp_hh, AArs_cc).scale(-0.5)
    # l1 = evaluate(AArs_cc, AAqp_hh)


    # -------------------------------------------------------------------------------------------------------
    # W scislym przypadku <0|[[dO, H],O+]|0> = <0|[dO,[H,O+]]|0>
    # <0|[dO, H, O+]|0>  = 1/2(<0|[[dO, H],O+]|0> + <0|[dO,[H,O+]]|0>)
    #
    # KL = <0|[[dO, H],O+]|0> = |jesli killer| =  <0|[dO, H],O+]0>
    # z tego <0|[[rs, H], q+p+]|0>  powinno wyjść to samo co z <0|[rs, H], q+p+|0>
    # KL = <0|[[rs, H], q+p+]|0>
    # KL_killer = <0|[rs, H], q+p+|0>  <--- tego nie bedziemy stosowac, bo ten powyzej jest leszpy - zwiazany
    #
    # KP = <0|[dO,[H,O+]]|0> = |jesli killer| = <0|dO,[H,O+]]0>
    # z tego  <0|[rs, [H, q+p+]|0>   powinno wyjść to samo co z   <0|rs, [H, p+q+]|0>  (w notacji komutatory na prawo)
    # z tego -<0|[[H, q+p+], rs|0>   powinno wyjść to samo co z   <0|rs, [H, q+p+]|0>  (w notacji komutatory na lewo)
    # KP = -<0|[[H, q+p+], rs|0>
    # KP_killer = <0|rs, [H, q+p+]|0>

    # wnioski
    # 1. Z komutatorow KL i KP powinno wyjsc to samo, jesli jest spelninony killer, czyli pq|0> = 0. Nigdzie pewnie tego nie zakladamy.
    # Zbadac czy takie wyrazy powstaja przed nazwaniem macierzy gestosci.
    

    # for PP
    # rs = AArs_aa
    # p+q+ = AAqp_pp

    # for HH
    # r+s+ = AArs_cc
    # pq = AAqp_hh

    # KL = <0|[[rs, H], q+p+]|0>
    # KL = evaluate(AArs_aa, h2cas, AAqp_pp)

    # KL = evaluate(AArs_aa, AAqp_pp)

    # KL = evaluate(h1cas, AAqp_pp, AArs_aa).scale(-1.0) +  evaluate(h2cas, AAqp_pp, AArs_aa).scale(-1.0)

    # KL = evaluate(AArs_aa, h1cas, AAqp_pp) +  evaluate(AArs_aa, h2cas, AAqp_pp)
    
    # KL = evaluate(AArs_aa, h1cas, AAqp_pp).scale(0.5) +  evaluate(AArs_aa, h2cas, AAqp_pp).scale(0.5) +\
    #     evaluate(h1cas, AAqp_pp, AArs_aa).scale(-0.5) +  evaluate(h2cas, AAqp_pp, AArs_aa).scale(-0.5)


    # KL = evaluate(AArs_cc, AAqp_hh)

    # KL = evaluate(h1cas, AAqp_hh, AArs_cc).scale(-1.0) +  evaluate(h2cas, AAqp_hh, AArs_cc).scale(-1.0)

    # KL = evaluate(AArs_cc, h1cas, AAqp_hh) +  evaluate(AArs_cc, h2cas, AAqp_hh)

    # KL = evaluate(AArs_cc, h1cas, AAqp_hh).scale(0.5) +  evaluate(AArs_cc, h2cas, AAqp_hh).scale(0.5) +\
    #     evaluate(h1cas, AAqp_hh, AArs_cc).scale(-0.5) +  evaluate(h2cas, AAqp_hh, AArs_cc).scale(-0.5)


    # -----sprawdzenie wersji jednowyznacznikowej
    # - utworz wszystkie mozliwe pary creator, anihilatr, virt, occ

    # ovlist = [['c', 'd'], ['c', 'k'], ['k', 'c'], ['k','l']]

    # ovlist = [['k', 'l']]
    # crealist = [['0', '0']]
    # Axx = cas()
    # Axx.operator_idx = ['a', 'b']
    # Axx.operator_type =['+', '+']

    # Ayy = cas()
    # Ayy.operator_idx = ['j', 'i']
    # Ayy.operator_type =['+', '+']

    # KL = arithmetic_string()
    # for v1 in range(0, len(ovlist)):
    #     for v2 in range(0, len(crealist)):
    #         Ars = cas()
    #         Ars.operator_idx = ovlist[v1]
    #         Ars.operator_type = crealist[v2]

    #         # XX = evaluate(h1cas, Axx, Ars).scale(-1.0) +  evaluate(h2cas, Axx, Ars).scale(-1.0)
    #         XX = arithmetic_string()
    #         YY = arithmetic_string()
    #         # XX = evaluate(Ars, h1cas, Axx) +  evaluate(Ars, h2cas, Axx)

    #         YY = evaluate(Ars, h1cas, Ayy).scale(-1.0) +  evaluate(Ars, h2cas, Ayy).scale(-1.0)
    #         # YY =  evaluate(Ars, h2cas, Ayy).scale(-1.0)
    #         KL = KL + XX
    #         KL = KL + YY
            
    #         print('')
    #         print('wynik komutatorow XX')
    #         k = 0
    #         for x in XX:
    #             print(k, x)
    #             k += 1
    #         print('')
    #         print('')
    #         print('wynik komutatorow YY')
    #         k = 0
    #         for x in YY:
    #             print(k, x)
    #             k += 1
    #         print('')
#    sys.exit(0)
#-----------------------------------------------------pojedynczy 
    # print('wynik komutatorow KL')
    # k = 0
    # for x in KL:
    #     print(k, x)
    #     k += 1
    # print('')

    # k = 0
    # for x in KL:
    #     print(k, x)
    #     k += 1
    # print('')


    # print('Perform Wick')
    # res = []
    # for x in KL:
    #     result_list = x.wick_ca()
    #     res = res + result_list
    # print('')

    # print('Result po WICK:')    
    # for x in res:
    #     x.exec_delta()
    #     print(x)
    # print('')

#--------------------------------------------------------------------------

    epq1 = cas()
    epq1.operator_idx = ['i', 'a', 'j', 'b', 'c','d','l','k']
    epq1.operator_type = ['+', '0', '+','0', '+','+','0','0']
    res  = epq1.wick_ca()


    #----------------------------------kod do zerowania macierzy gestosci dla jednowyznacznikowych----------
    resal = arithmetic_string()
    for x in res:
        print('sprawdzam', x, x.operator_idx)
        if len(x.operator_idx) > 0:
            if x.operator_idx[0] in virtual:
                if x.operator_type[0] == CRE:
                    x.num_factor = 0
            if x.operator_idx[-1] in virtual:
                if x.operator_type[-1] == ANI:
                    x.num_factor = 0
            if x.operator_idx[0] in occupied:
                if x.operator_type[0] == ANI:
                    #print('taaaaak')
                    x.num_factor = 0
            if x.operator_idx[-1] in occupied:
                if x.operator_type[-1] == CRE:
                    x.num_factor = 0
            if len(set(x.operator_idx).intersection(virtual)) != 0:
                x.num_factor = 0.0
                print('jest zero')

            if len(x.operator_idx) == 2:
                x.new_delta(x.operator_idx[0], x.operator_idx[1])
                x.operator_idx = []
            if len(x.operator_idx) == 4:
                res1 = deepcopy(x)
                res2 = deepcopy(x)
                res1.new_delta(x.operator_idx[0], x.operator_idx[2])
                res1.new_delta(x.operator_idx[1], x.operator_idx[3])
                res1.num_factor *= -1.0
                res2.new_delta(x.operator_idx[0], x.operator_idx[3])
                res2.new_delta(x.operator_idx[1], x.operator_idx[2])
                res1.operator_idx = []
                res2.operator_idx = []
                x.num_factor = -100


        if x.num_factor == 0:
            print('jest zero')
        else:
            print('nie ma')
            if x.num_factor != -100:
                resal = resal + arithmetic_string(x)
            else:
                print('tak dodaje te dwa', res1, res2, x)
                resal = resal + arithmetic_string(res1)
                resal = resal + arithmetic_string(res2)

    res = deepcopy(resal)

        # print(x)
        # print('')
    print('')
    #----------------------------------koniec kodu do zerowania macierzy gestosci dla jednowyznacznikowych----------

    print('Result po WICK:')    
    for x in res:
        x.exec_delta()
        if (x.num_factor) != 0.0:
            print(x)
        print(x)
    print('')


    print('Result po RENAME DENSITY:')    
    for x in res:
        x.rename_as_density()
        if (x.num_factor) != 0.0:
            print(x)
    print('')

    res2 = cas_to_ugg(res)
    print('po tranformacji')
    for x in res2:
        print(x)
    
    rsimp = simplify(res2, cas=True)
    k = 0
    print('po simp')
    for x in rsimp:
        k += 1
        # print("&", x, "\\\\")
        print(x)
    sys.exit(0)

    

    print('')
    print('wynik komutatorow KL')

    k = 0
    for x in KL:
        print(k, x)
        k += 1
    print('')

    print('Perform Wick')
    res = []
    for x in KL:
        result_list = x.wick_ca()
        res = res + result_list
    print('')

    print('Result po WICK:')    
    for x in res:
        x.exec_delta()
        print(x)
    print('')

    print('Result po RENAME DENSITY:')    
    for x in res:
        x.rename_as_density()
        print(x)
    print('')

    res2 = cas_to_ugg(res)
    print('po tranformacji')
    for x in res2:
        print(x)
    
    rsimp = simplify(res2, cas=True)
    k = 0
    print('po simp')
    for x in rsimp:
        k += 1
        # print("&", x, "\\\\")
        print(x)
    sys.exit(0)
    #---------------------------------------------------------------

    epq = cas()
    epq.operator_idx = ['k', 'l', 't', 'u', 'v', 'w']
    epq.operator_type = ['0', '0', '+', '+', '0', '0']
    res = epq.wick_ca()
    for x in res:
        print(res)
    epq2 = cas()
    epq2.operator_idx = ['t', 'u', 'v', 's']
    epq2.operator_type = ['+', '+', '0', '0']


    # l1 = evaluate(epq, AApq)
    # l1 = evaluate(epq2, AApq)

    # epq = cas()
    # epq.operator_idx =  ['p', 'q', 'r', 's']
    # epq.operator_type = ['0', '0', '0', '+']

    # epq = cas()
    # epq.operator_idx =  ['r', 's', 'p', 'q']
    # epq.operator_type = ['0', '+', '0', '0']


    for x in l1:
        print(x)
    print('')
    print('Perform Wick')
    res = []
    for x in l1:
        result_list = x.wick_ca()
        res = res + result_list
    print('')

    
    # print('Perform Wick Theorem')
    # result_list = l1.wick_ca()
    # print('')
    # print('Result:')    
    # for x in result_list:
    #     print(x)
    # print('')
    sys.exit(0)
    
    epq = cas()
    epq.coefficient = ['h']
    epq.coefficient_idx = [['r', 's', 't', 'u']]
    epq.summation = ['r', 's', 't', 'u']
    epq.operator_idx = ['s', 't', 'r', 'u']
    epq.operator_type = ['+', '0', '+', '0']

    print(epq)
    temp1, e2_r = epq.right_split()
    print('temp', temp1)
    print('e2r', e2_r)
    sys.exit(0)
    print('Creator-Anihilator-string:', epq)
    print('')

    print('Perform Wick Theorem')
    result_list = epq.wick_ca()
    print('')
    print('Result:')    
    for x in result_list:
        x.rename_as_density()
        print(x)
    print('')

    # epq.normal_order()
    # epq.contraction([[1,2]])
    # contr = epq.all_contractions()
    # epq.remove_null_contractions(contr)
    # print(epq)
    # print(epq)


def tda_comp_idx():

    

    l1 = [1,2,3]#,4,5,6]#,3,4]
    vir = [3, 4, 5, 6]
    occ = [1,2]
    l2=[3]

    N = len(vir)
    Mp = len(occ)
    Pp = len(vir)

    A = int(N*(N+1)/2)
    B = int(Mp*(Mp+1)/2)
    C = int(Pp*(Pp+1)/2)


    for k in range(occ[0], occ[-1]+1):
        for l in range(k, occ[-1]+1):
            kk = k - occ[0] +1
            ll = l - occ[0] + 1
            
            kl = int((2*Mp-l+2)*(l-1)/2 +k -l +1)

            for qq in range(vir[0], vir[-1]+1):
                for pp in range(qq,vir[-1]+1):
                    p = pp - vir[0] +1
                    q = qq - vir[0] + 1
            
            
                    pq = int((2*N-q+2)*(q-1)/2 +p -q +1)
                    print(k, l, pp, qq, pq+ A*(kl-1))
            print('--------------------')
    sys.exit(0)

        
def comp_idx():

#    l1 = [1,2,3,4]
#    l1 = [1,2]
    
#    l2 = [1,2]

    l1 = [1,2,3]#,4,5,6]#,3,4]
    vir = [3, 4]
    occ = [1,2]
    l2=[3]

    N = len(vir)
    Mp = len(occ)
    Pp = len(vir)

    A = N*(N+1)/2
    B = Mp*(Mp+1)/2
    C = Pp*(Pp+1)/2


    dct = {}
    # for p in range(l1[0],l1[-1]+1):
    for pp in range(vir[0],vir[-1]+1):
        for qq in range(vir[0], pp+1):
            p = pp - vir[0] +1
            q = qq - vir[0] + 1
            
            pq = int((2*N-q+2)*(q-1)/2 +p -q +1)

            for k in range(occ[0],occ[-1]+1):
                for l in range(1, k+1):

                    klp = int((2*Mp-l+2)*(l-1)/2 +k -l +1)

                    pqkl = int(klp + (pq-1)*B)
                                        
                    dct[pqkl] = [pp, qq, k, l]
#    pqkl = 18
#    start = pqkl
    # print('pqkl', pqkl)
#     for p in range(l1[0],l1[-1]+1):
#         for q in range(1, p+1):
#             pq = int((2*N-q+2)*(q-1)/2 +p -q +1)
#             for k in range(vir[0],vir[-1]+1):
#                 for l in range(vir[0], k+1):

#                     kp = k - vir[0] +1
#                     lp = l - vir[0] + 1
# #                    print(kp, lp)
#                     klp = int((2*Pp-lp+2)*(lp-1)/2 +kp -lp +1)
#  #                   print(pq, klp)
#                     pqkl = start + int(klp + (pq-1)*C)
#                     dct[pqkl] = [p, q, k, l]
#                     # print(pqkl, dct[pqkl])

    dct2 = {}
    for i in range(0, len(dct)):
        for key in dct:
            if key == i:
                dct2[i] = dct[key]
                exit
            
    for key1 in dct2:
        print('key', key1, dct2[key1])
    print(len(dct2))
#    sys.exit(0)

    
#     dct = {}
#     for p in range(l1[0],l1[-1]+1):
#         for q in range(1, p+1):
#             pq = int((2*N-q+2)*(q-1)/2 +p -q +1)
#             #print(p, q, pq)
#             # print(l2[0], len(l2))
#             for k in range(l1[0],l1[-1]+1):
#                 for l in range(1, k+1):
# #                    klp = int((2*N-l+2)*(l-1)/2 +k -l +1)
#                     klp = int((2*N-l+2)*(l-1)/2 +k -l +1)

#             # for k in range(l2[0],l2[-1]+1):
#             #     for l in range(l2[0], k+1):
#             #         kp = k - l2[0] +1
#             #         lp = l - l2[0] + 1
#             #         kl = int((2*Mp-l+2)*(l-1)/2 +k -l +1)
                    
#             #         klp = int((2*Mp-lp+2)*(lp-1)/2 +kp -lp +1)
#                     #                    print(k, l, kp, lp, kl, klp)
                                        
# #                    pqkl = int(pq + (klp-1)*A)
#                     pqkl = int(klp + (pq-1)*A)
                                        
#                     dct[pqkl] = [p, q, k, l]

#                     # print(p, q, k, l, pqkl)

#     dct2 = {}
#     for i in range(0, len(dct)):
#         for key in dct:
#             if key == i:
#                 dct2[i] = dct[key]
#                 exit
            

    dct = deepcopy(dct2)
    tab = []
#    print(len(dct))
    for key1 in dct:
 #       print('key', key1)
        tabmini = []
        for key2 in dct:
            tabmini.append(0)
        tab.append(tabmini)
    # for key in dct:
    #     print(key)
    #     print(key, tab[key-1])
    
    # sys.exit(0)
    # print('')
    # print('start')

    for key1 in dct:
        for key2 in dct:
            w1 = dct[key1]
            w2 = dct[key2]

            #w1 = [p, q, k, l]
            #w2 = [i, j, r, s]
            # print(key1, dct[key1], key2, dct[key2])
            cc = w1[0]
            dd = w1[1]
            kk = w1[2]
            ll = w1[3]

            aa = w2[0]
            bb = w2[1]
            ii = w2[2]
            jj = w2[3]


            pp = w1[0]
            qq = w1[1]
            kk = w1[2]
            ll = w1[3]

            rr = w2[0]
            ss = w2[1]
            ii = w2[2]
            jj = w2[3]

            # print(w1, w2)
            # print('dup0', dct[key2-1], dct[key1-1], aa, bb, ii, jj, cc, dd, kk, ll)

            # if ii in occ and jj in occ and aa in vir and bb in vir:
            #     if ll == ii and aa ==cc:
            #         tab[key2-1][key1-1] = 1
            #         print('dd', dct[key2], dct[key1], key2, key1)
                    

            #onedet
            if ii in occ and jj in occ and ii == kk and jj == ll and aa in vir and bb in vir:
                if aa == cc and bb == dd:
                    tab[key2-1][key1-1] = 1
                    print('kk', dct[key2], dct[key1], key2, key1)
                elif aa == dd and bb == cc:
                    tab[key2-1][key1-1] = 1
                    print('kk', dct[key2], dct[key1], key2, key1)


            if ii in occ and jj in occ and ii == ll and jj == kk and aa in vir and bb in vir:
                if aa == cc and bb == dd:
                    tab[key2-1][key1-1] = 1
                    print('kk', dct[key2], dct[key1], key2, key1)
                elif aa == dd and bb == cc:
                    tab[key2-1][key1-1] = 1
                    print('kk', dct[key2], dct[key1], key2, key1)
            

            # if qq in vir and pp in vir and ss in vir and rr in vir:

            #     if qq == rr and kk == ii and ll == pp and ss == jj:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if qq == rr and kk == ii and ll == jj and ss == pp:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if pp == ss and kk == qq and ll == jj and rr == ii:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if pp == ss and kk == ii and ll == jj and rr == qq:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if ii == kk and ll == pp and rr == qq and ss == jj:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if ii == kk and ll == jj and rr == qq and ss == pp:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if jj == ll and kk == qq and rr == ii and ss == pp:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if jj == ll and kk == ii and rr == qq and ss == pp:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])

            #     if qq == rr and ll == jj:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if qq == rr and ss == pp:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if qq == rr and kk == ii:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if qq == rr and ss == jj:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if qq == rr and ll == pp:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if pp == ss and rr == ii:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if pp == ss and ll == jj:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if pp == ss and kk == qq:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if pp == ss and rr == qq:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if pp == ss and kk == ii:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if ii == kk and ss == jj:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if ii == kk and ll == pp:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if ii == kk and ll == jj:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if ii == kk and ss == pp:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if ii == kk and rr == qq:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if jj == ll and kk == qq:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if jj == ll and rr == ii:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if jj == ll and kk == ii:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if jj == ll and ss == pp:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if jj == ll and rr == qq:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if pp == ss and qq == rr:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])
            #     if ii == kk and jj == ll:
            #         tab[key2-1][key1-1] = 1
            #         print(key2, key1, dct[key2], dct[key1])


                
            # if w1[0]==w2[3] and w1[1]==w2[2]:
              #  print(key2, key1)

    print('stop')
                
        #    for key in dct:
    for i in range(0, len(dct)):
        print(tab[i])
                

            
    sys.exit(0)



# def indeksy_s():

#     # 0,1,2 - occupied
#     # 3,4,5 - virtual
#     # 0 - 5 - general


    

def test1():
    
    # AA = cas()
    # AA.operator_idx = ['c', 'k', 'd', 'l']
    # AA.operator_type = ['+', '0', '+', '0']

    AA = cas()
    AA.operator_idx = ['c', 'l']
    AA.operator_type = ['+','0']

    # BB = cas()
    # BB.operator_idx = ['i', 'b', 'j', 'a']
    # BB.operator_type = ['+', '0', '+', '0']

    BB = cas()
    BB.operator_idx = ['i', 'a']
    BB.operator_type = ['+', '0']

    
    KL = evaluate(BB, AA)

    res = int_and_simp_cas_onedet(KL)
    
    res = zero_dens_matrix_onedet(res)
    res = clean_cas(res)

    print('Result po RENAME DENSITY:')    
    for x in res:
        x.rename_as_density()
    print('')


    print('')
    print('wynik')
    print('')
    for x in res:
        print(x)

    res2 = cas_to_ugg(res)    
    rsimp = simplify(res2, cas=True)

    print('')
    print('wynik')
    print('')
    for x in rsimp:
        print(x)

    res = ugg_to_cas(rsimp)
    res = approximate_rdm3_by_rdm12(res)
    res = approximate_rdm2_by_rdm1_when_occ(res)
    print('wynik')
    print('')
    for x in res:
        print(x)
