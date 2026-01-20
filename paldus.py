
# ------------------------------------------------------
# Author: Aleksandra Tucholska, University of Warsaw
# ------------------------------------------------------
from params import *
import paldus_classes
from paldus_cc import *
from ccsd_f12 import *
from paldus_eom import *
from paldus_eom_mem import *
from paldus_eom_triplet import *
from paldus_basic import *
from paldus_density import *
from paldus_commutators import *
from paldus_omega import *
from eomccjac import jacobian_loop
from paldus_classes import ugg
from copy import deepcopy
from paldus_classes import arithmetic_string
from paldus_classes import disambiguate
from collections import deque
from paldus_classes import pair_permutations
from paldus_cas import *
from paldus_acm import erpa_ph
from paldus_acm import pp_rpa
from paldus_acm import pphh_rpa
from paldus_acm import test1
from paldus_acm import comp_idx
from paldus_acm import tda_comp_idx

from paldus_acm import pphh_rpa_S_onedet
from paldus_acm import pphh_rpa_A_onedet
from paldus_acm import pphh_rpa_A
from paldus_acm import TDA_pphh_ph_onedet
from paldus_acm import test_idx_for_pphh
from paldus_acm import test_rdm_approx
from paldus_acm import ph_rpa
from paldus_acm import hh_rpa_koopmans
from paldus_acm import test_rpa
from paldus_acm import pp_rpa_hf_det
from paldus_acm import *
import math
from itertools import product
from itertools import permutations
from paldus_classes import integrate
from paldus_basic import simplify
from paldus_classes import virtual, occupied, general
from slater_rules import *
from factor import *
from paldus_f12 import *
#from fortran_code import *
from templates import *
from templates_factor import *
#from paldus_diagrams import simple_diag
import sys
import time
import pickle
#from paldus_male_testy import *
from paldus_torun import Tlatex_read
from paldus_polaritonic import polarit_test


EPSILON = 10**(-14) 

##
# @mainpage PALDUS
#
# @section description_main Description
#The Paldus code is  designed to derive, simplify, and automatically implement expressions of the type
# \f[\braket{[V_1,\mu_n]_{k_1}|[V_2, V_3]_{k_2}|[V_4,\nu_m]_{k_3}}\f]
#where \f$k_1, k_2, k_3\f$ denote \f$k\f$-tuply nested commutators. The operators \f$V_1-V_4\f$ could be any excitation, de-excitation,
#or general operators that are represented by the products of the \f$E_{pq} \f$ and/or \f$T_{pq}\f$ operators.
#Each of the integrals is approximated within the requested level 
#of theory and integrated using the Wick's theorem.\cite{wick1950evaluation}
# 
#
# @section notes_main Notes
# - Add special project notes here that you want to communicate to the user.
#
# Copyright (c) 2021 Aleksandra Tucholska.  All rights reserved.
##
# """! lala \f$ f(x) = e^x \f$ witam"""


# r = evaluate(t2,t1c)
# for x in r:
#      print(x)
# sys.exit(0)

# zobs = ugg()
# zobs.summation = ['p','q']
# zobs.coefficient = [OBSERVABLE_X]
# zobs.coefficient_idx.append(['p','q'])
# zobs.operator_idx.append(['p','q'])
# zobs.operator_type.append("s")

# r = evaluate(obs, t1)
# for x in r:
#      print(x)


# d = ugg()
# d.summation = ['b', 'j']
# d.operator_idx.append(["j","b"])
# d.operator_type.append("s")

# print('asd')
# for x in range(0, len(r)):
#      r[x] = r[x].fromleft(d)
#      print(x)
# print('rrrrrrrrr')
# print(r[1])
# rint = r.integrate().scale(0.5)
# print('integrate')
# for y in rint:
#      print(y)
# sys.exit(0)

#----------------------------------------

# NEVPT

#nevpt2()

#sys.exit(0)

# CCSD                                                                                                                                                 
#----------------------------------------
# print('Execute CCSD')
# execute_ccsd()
# sys.exit(0)

#---------------------------------------
# -----Równania dla Michala
#execute_michal(3)
# execute_michal2()
# sys.exit(0)
#----------------------------------------


#---------------------------------------
# -----Macierz gęstości dla Polnona
#----------------------------------------


#execute_polnon_72_73()
#execute_ccsd()
# generate_t1_in_2_mbpt()
# generate_t1_in_3_mbpt()
# generate_t2_in_1_mbpt()
# generate_t2_in_2_mbpt()
# generate_t2_in_3_mbpt()
# generate_density_Wm_ground_driver('dump', 'ccd')
#sys.exit(0)
#-------------------------------------


#--------------------------------GENERATE DIAGRAMS---------------------------------------
# simple_diag()
# sys.exit(0)


#--------------------------------GENERATE S WITH F12---------------------------------------
# test_p_operator(4)
# sys.exit(0)

# test_bazy_biort()
#evaluate_s_f12_operator_separate_functions("load")
#sys.exit(0)
#------------------------------------s_OPERATORS----------------------------------------------------------

# test_p_operators()
# sys.exit(0)

# evaluate_s_operator_separate_functions()
# # generate_s_operators() 
# sys.exit(0)

#------------------------------------Density Matrix Ground F12----------------------------------------------------------                                                
# print('generate_density_matrix_ground_f12')
#generate_density_matrix_ground_f12_driver('dump', 'ccsd')
#generate_density_matrix_ground_f12('load', 'ccsd', 5)
#sys.exit(0)
#--------------------------------------Cumulant Ground CCD----------------------------------------------------
# generate_cumulant_ground_driver('dump', 'f12')
# sys.exit(0)
#------------------------------------Cumulant Ground F12----------------------------------------------------------                            
# generate_cumulant_ground_f12_driver('load', 'f12', only_quadratic=True)
# sys.exit(0)


#--------------------------------------------TOURUN-----------------------------------------------
# Tlatex_read()
# sys.exit(0)
#--------------------------------------------ERPA-----------------------------------------------
# twoel()
# sys.exit(0)
# hh_rpa_koopmans()
#pphh_rpa()
#comp_idx()
# test_macierzy_S_onedet()
#test1()
# test_komutatora()
#pphh_rpa_S_onedet()
#pphh_rpa_A_onedet()
#pphh_rpa_A("wynik-pphh-og.dat", 1, True)
#test_idx_for_pphh()
#test_nowy()
#tda_comp_idx()
#TDA_pphh_ph_onedet(22, 3, "wynik-h22-popr2.dat", False)
#test_rdm_approx()
#execute_cisd(2,2)
#execute_rowe(2,1, 1, 1)
#execute_cisd_rowe_triplet(2, 2, True, 1, -0.5)
# TDA_pphh_ph_onedet(22, 1, "rdm4.dat", True)

#pp_rpa(diagonal=True)
#hh_rpa_koopmans()

#pp_rpa_hf_det()
# r = evaluate(XXp, YYq)
# r = evaluate(XXp, AArs)
# r = evaluate(AArs, AApq)
# r = evaluate(AArs_aa, AApq_pp)
# for x in r:
#     print(x)
# sys.exit(0)     
#ph_rpa_plus_doubles('dump', cis='ph')
#ph_rpa_plus_doubles('dump2')
#ph_rpa_plus_doubles('load')
#ph_rpa_plus_doubles_overlap()
#ph_rpa(cis = False, spinres=True)
polarit_test()

# test_rpa()
sys.exit(0)








# r =nukl
# rint=r.integrate(bra=['i', 'j'], braspin=['s'])
# #rint = r.integrate()
# print(rint)
# sys.exit(0)


# aa = ugg()
# aa.summation = ["p1","a","r","a1"]
# aa.coefficient = [TWOEL_INT]
# aa.coefficient_idx.append(["p1","a","r","a1"])
# aa.delta.append(["p1","a1"])
# aa.num_factor = -1./2.

# print(aa)
# aa.exec_delta()
# print(aa)
# sys.exit(0)


# r = arithmetic_string(t2fa)
# rsimp = simplify(r)
# print(rsimp)
# sys.exit(0)

# r = nuaibjck
# print(r)
# rint = r.integrate(bra = ['d', 'l', 'e', 'm', 'f', 'n'], braspin = ['s', 's', 's'])
# rsimp = simplify(rint)
# for x in rsimp:
#      print(x)
# sys.exit(0)


# r = t1
# r.coefficient.append('t')
# r.coefficient_idx.append(['b', 'j'])
# r.operator_idx.append(["b", "j"])
# r.operator_type.append("s")
# r.summation = ["a","i", 'b', 'j']

# print(r)
# print('')
# rr = evaluate(hamiltoniant, r)

# for x in rr:
#      print(x)
# print('')
# rint = rr.integrate()
# rsimp = simplify(rint)
# print('wynik')
# for x in rsimp:
#      print(x)

# sys.exit()



#r1 = hamiltoniant + evaluate(hamiltoniant, t2)
#r1 = hamiltonian + evaluate(hamiltonian, t1) + evaluate(hamiltonian, t2) +evaluate(hamiltonian, t1, t1).scale(0.5) + evaluate(hamiltonian, t1, t2)
#r1 = evaluate(hamiltoniant, t2fa) + evaluate(hamiltonian, t2fb) 


#-------------------------------------------------------------GENERATE CCSD WITH F12---------------------------------
"""! main paldus file beginning"""
# sladkisz
#execute_ccsd_f12("ccsd", "dump", "energy", True)
# execute_ccsd_f12("ccsd", "load", "energy", True)
# sys.exit(0)
# execute_ccsd_f12("ccsd", "dump", "t1", True)
# execute_ccsd_f12("ccsd", "load", "t1", True)
# sys.exit(0)
# execute_ccsd_f12("ccsd", "dump", "t2", True)
# execute_ccsd_f12("ccsd", "load", "t2", True)
# sys.exit(0)
# execute_ccsd_f12("ccsd", "dump", "tf", True)
execute_ccsd_f12("ccsd", "load", "tf", True)
sys.exit(0)
# print(t2fa)

# print(hamiltonian_comp)
# print('')

#r1 = evaluate(hamiltonian_comp_fock, t2fa)
# r1 = hamiltonian_comp + evaluate(hamiltonian_comp, t1)+ evaluate(hamiltonian_comp, t2)+ evaluate(hamiltonian_comp, t2fa)\
#     + evaluate(hamiltonian_comp, t1, t2) + evaluate(hamiltonian_comp, t1, t2fa)+ evaluate(hamiltonian_comp, t1, t1).scale(0.5)+\
#     evaluate(hamiltonian_comp, t1, t1, t1).scale(0.16666666666666) \
#     + evaluate(hamiltonian_comp, t1, t1, t2).scale(0.5) + evaluate(hamiltonian_comp, t1, t1, t2fa).scale(0.5) + evaluate(hamiltonian_comp, t2, t2fa).scale(0.5)\
#     + evaluate(hamiltonian_comp, t2, t2).scale(0.5) + evaluate(hamiltonian_comp, t2fa, t2fa).scale(0.5) \
#     + evaluate(hamiltonian_comp, t1, t1, t1, t1).scale(0.0416666666)
# r1 = evaluate(hamiltonian_comp, t1, t1, t2fa).scale(0.5)
# r1 = evaluate(hamiltonian_comp, t2fa)
# r1 = evaluate(hamiltonian_comp, t2fa)
# r1 = arithmetic_string(hamiltonian_comp) + evaluate(hamiltonian_comp, t2)+ evaluate(hamiltonian_comp, t2fa)\
#     + evaluate(hamiltonian_comp, t2, t2fa).scale(0.5)\
#     + evaluate(hamiltonian_comp, t2, t2).scale(0.5) + evaluate(hamiltonian_comp, t2fa, t2fa).scale(0.5) 
# sys.exit(0)


# r2 = evaluate(hamiltoniant_comp, t2) + evaluate(hamiltoniant_comp, t2fa)\
#     + evaluate(hamiltoniant_comp, t2, t2fa)\
#     + evaluate(hamiltoniant_comp, t2, t2).scale(0.5) + evaluate(hamiltoniant_comp, t2fa, t2fa).scale(0.5) 
# r2 = evaluate(hamiltoniant_comp, t2fa)

# print('przed calkowaniem')
# for x in r1:
#      print(x)
# print('po calkowaniu')


# r1 = r1.fromleft(gemi)
#print('oto', gemi)

# k =1
# r3 = arithmetic_string()
# for x in r1:
#      print(k, x)
#      disambiguate(x, gemi)
#      y = x.fromleft(gemi)
#      r3 = r3 + arithmetic_string(y)
#      print(k,  y)
#      print('')
#      k += 1
# print(len(r3))

# print('oto gemi')
# for x in r3:
#      print(x)
# print('')

# sys.exit(0)

# #r1 = evaluate(hamiltonian_comp, t2fa) + evaluate(hamiltonian_comp, t1, t2fa)
# #r1 = arithmetic_string(hamiltonian_comp) + evaluate(hamiltonian_comp, t1)+ evaluate(hamiltonian_comp, t2)+ evaluate(hamiltonian_comp, t1, t1).scale(0.5) 
# #     evaluate(hamiltonian_comp, t1, t1, t1).scale(0.16666666666666)


# #r1 = evaluate(hamiltonian_comp, t2fa)
# # print('t2', t2)
# # r1 = evaluate(hamiltonian, t2)
# # print('przed calkowaniem')
# # r2 = arithmetic_string()
# # for t in range(0, 10):
# #      r2 += arithmetic_string(r1[t])
# r2 = r1
# #print(r1[0])
# for x in r2:
#      print(x)
# print('')

rint1 = r1.integrate()
#rint1 = r1.integrate(bra = ['a', 'i'], braspin = ['s']).scale(0.5)
# start1 = time.time()
#rint1 = r1.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s', 's']).scale(0.5)
# rint1 = r3.integrate().scale(0.5)
# end = time.time()
# #rint1 = r2.integrate()
# print('po calkowaniu', end - start1)
# print('')
for x in rint1:
     print(x)
print('')
rint1.exec_delta()
print('po exec delta')
for x in rint1:
     print(x)
print('')
print('ZACZYNAM RSIMP')

rsimp1 = simplify(rint1)
print('po rsimp1')
for x in rsimp1:
     print(x)
print('')
sys.exit(0)
# print('regroup')
# rgroup = find_L(rsimp1)
# print('po find L')
# for x in rgroup:
#      print(x)
# print('')


# rsimp2 = simplify(rgroup)
# print('po rsimp ost', len(rsimp2))
# for x in rsimp2:
#      print(x)
# print('')

res = find_permutations(rsimp1)
#res = rsimp1
print('po find permutations')
for p in res:
     print(p)
print('')
print(' a teraz identify interm')
res2 = identify_interm_V_f12(res)
print('po find V')
for p in res2:
     print(p)
print('')
res3 = identify_interm_X_f12(res2)
print('po find X')
for p in res3:
     print(p)
print('')
res4 = identify_interm_B_f12(res3)
print('po find B')
for p in res4:
     print(p)
print('')
res5 = identify_interm_P_f12(res4)
print('po find P')
for p in res5:
     print(p)
print('')
s = arithmetic_string()
for x in res5:
     s = s + arithmetic_string(x)

#sys.exit(0)
print('')

print('')
ss = s.cabstransform()

print('po cabs transform')
for x in ss:
     print(x)

res_ost = identify_interm_Ft_f12(ss)
print('i to jest ostateczny wyniczor')
for x in res_ost:
     print(x)
# res6 = cabssplit(res5)

sys.exit(0)


# for x in rsimp1:
#      x.optimize()
# rsimp1.cleanup()

# print('wynik')
# for x in rsimp1:
#      print(x)
# sys.exit(0)



# print(t2fa)
# print(t2fb)
# print('')
# #r = evaluate(hamiltoniant, t2fb)
# r = evaluate(hamiltonian, t1) + evaluate(hamiltonian, t2) + evaluate(hamiltonian, t1, t1).scale(0.5)
# #r = evaluate(hamiltoniant, t2) + evaluate(hamiltoniant, t2fa) + evaluate(hamiltoniant, t2fb)
# for x in r:
#      print(x)
# rint = r.integrate()
# rsimp = simplify(rint)
# print('wynik')
# for x in rsimp:
#      print(x)

# sys.exit(0)

# #r = ham +  evaluate(ham, t1) + evaluate(ham, t2) + evaluate(ham, t1, t1)
# r = hamiltonian +  evaluate(hamiltonian, t1) + evaluate(hamiltonian, t2) + evaluate(hamiltonian, t1, t1)

# rint = r.integrate()
# rsimp = simplify(rint)
# for x in rsimp:
#      print(x)

# sys.exit(0)
#---------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------FACTORIZATION--------------------

# lista1 = [[[1], [2, 3]], [[2], [1, 4]], [[3], [1, 4]], [[4], [2, 3]]]

lista2 = [[[1], [3, 7]], [[2], [3, 6, 7]], [[3], [1, 2, 4]], [[4], [3, 5, 7]], \
              [[5],[4, 6]], [[6], [2,5]], [[7], [1, 2,4]]]

#e = bb13
#e = e.fromright(bb12)
#e = e.fromright(bb11)
e = bb11
e = e.fromright(bb10)
e = e.fromright(bb9)
e = e.fromright(bb8)
e = e.fromright(bb7)
e = e.fromright(bb6)
e = e.fromright(bb5)
e = e.fromright(bb4)
e = e.fromright(bb3)
e = e.fromright(bb2)
e = e.fromright(bb1)
print(e)
#sys.exit(0)
# e = bb7
# e = e.fromright(bb5)
# e = e.fromright(bb8)
# e = e.fromright(bb9)
# e = e.fromright(bb10)

# print(e)

# e = b1
# e = e.fromright(b2)
# e = e.fromright(b3)
# e = e.fromright(b4)
# e = e.fromright(b5)
# print(e)
# # e = ww1
# # e = e.fromright(ww2)
# # e = e.fromright(ww3)
# # e = e.fromright(ww4)


# print(e)
# #sys.exit(0)

e_idx  = find_all_diff_idx(e)
e_cost = idx_to_cost(e_idx)
print(e_idx)

# factorize(e)
list1, list2, dupa  = factorize(e)
sys.exit(0)

# idx_matrix, idx_sum, len_idx, len_coef = construct_idx_matrix(e, e_idx)
# factorize(idx_matrix, idx_sum, len_idx, len_coef, e_cost)



# #pairs_list = find_pairs(idx_matrix, idx_sum, len_idx, len_coef)
# permut = find_permutations(lista2)
# print('')
# print('wynik')
# print('')
# k = 1
# for  x in permut:
#     print(k, x)
#     k+=1

# sys.exit(0)
# # permut = find_permutations(lista2)
# print('aaaa')
# k = 1
# for  x in permut:
#     print(k, x)
#     k+=1
#sys.exit(0)


#0------------------------------------SLATER CONDON RULES-------------------------------------------------------------------
# slater_rules_driver(4, 4)
# sys.exit(0)
#--------------------------------------------------------------

#r = ham

# rl = ugg()
# rl.summation = ['a', 'b', 'i', 'i']
# rl.coefficient = [CC_AMPLITUDE]
# rl.coeficient_idx = [['a','i','b', 'j']]

# rr = ugg()
# rr.summation = ['c', 'k', 'd', 'l']
# rr.coefficient = [CC_AMPLITUDE]
# r.coeficient_idx = [['c','k','d', 'l']]


# c2 = ugg()
# c2.summation = ["a","i","b","j"]
# c2.coefficient = [CC_AMPLITUDE]
# c2.coefficient_idx.append(["a","i", "b", "j"])
# c2.operator_idx.append(["a", "i"])
# c2.operator_idx.append(["b", "j"])
# c2.operator_type.append("s")
# c2.operator_type.append("s")
# c2.num_factor = 1./2.

# d2 = ugg()
# d2.summation = ["c","k","d","l"]
# d2.coefficient = [CC_AMPLITUDE]
# d2.coefficient_idx.append(["c","k", "d", "l"])
# d2.operator_idx.append(["k", "c"])
# d2.operator_idx.append(["l", "d"])
# d2.operator_type.append("s")
# d2.operator_type.append("s")
# d2.num_factor = 1./2.


# r1 = r.fromright(c2)
# r2 = r1.fromleft(d2)
# for x in r2:
#     print(x)
# #sys.exit(0)
# rint1 = r2.integrate()#bra = [], ket = [], braspin = [], ketspin = [])
# print('int')
# for x in rint1:
#     x.exec_delta()
# print('a')
# for x in rint1:
#     print(x)
# rsimp = simplify(rint1)
# print('rsim')
# for x in rsimp:
#     print(x)
# sys.exit(0)

#r =ugg()
#r.summation = ['a','b','c','d','i','j','k','l']
#r.coefficient = [EOM_CC_SINGLE_Rr, EOM_CC_SINGLE_Rr_plus, S_AMPLITUDE, CC_AMPLITUDE]
#r.coefficient_idx=[['a','p','b','i','c','j'],['b','k','c','l'],['d','k','q','l'],['a','i','d','j']]


# mu1 = ugg()
# mu1.operator_idx.append(["a", "i"])
# mu1.operator_type.append("s")

# mu2 = ugg()
# mu2.operator_idx.append(["a", "i"])
# mu2.operator_idx.append(["b", "j"])
# mu2.operator_type.append("s")
# mu2.operator_type.append("s")
# mu2.operator_type.append("s")

# mu3 = ugg()
# mu3.operator_idx.append(["a", "i"])
# mu3.operator_idx.append(["b", "j"])
# mu3.operator_idx.append(["c", "k"])
# mu3.operator_type.append("s")
# mu3.operator_type.append("s")
# mu3.operator_type.append("s")


# r = arithmetic_string(obst)

# r = ugg()
# r.operator_idx.append(['i', 'a'])
# r.operator_idx.append(['b', 'j'])
# #r.operator_idx.append(['m', 'n'])
# r.operator_type.append("t")
# r.operator_type.append("t")
# #r.operator_type.append("s")
# r.num_factor = 1.0

# r = arithmetic_string(r)
# rint = r.integrate()
# #rint = r.integrate(bra= ['i', 'j'], braspin = ['t'], ket = ['k', 'l', 'm', 'n'], ketspin = ['t', 's'])

# for x in rint:
#     print(x)
# sys.exit(0)

# r = evaluate(obsab, t1)

# rint = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s', 's'])
# for x in r:
#     print(x)
# sys.exit(0)


# r = arithmetic_string(t2c)
# disambiguate(t2, obsij, obsia, obsai, obsab)
# r2 = r.fromright(obsij) + r.fromright(obsia) + r.fromright(obsai) + r.fromright(obsab)
# for x in r2:
#     print(x)

# rint = r2.integrate(ket = ['a', 'i', 'b', 'j'], ketspin = ["s", "s"])#.scale(0.5)

# rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], braspin = ["s", "s"])#.scale(0.5)
# rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
# rsimp = simplify(rint)
# print('calka')
# for x in rsimp:
#     print(x)

# sys.exit(0)

# sys.exit(0)

# r1 = ugg()
# r1.summation = ["i", "a" , "b", "j"]
# r1.operator_idx.append(["i", "a"])
# r1.operator_idx.append(["j", "b"])
# r1.operator_type.append("s")
# r1.operator_type.append("s")
# r1.num_factor = 0.5

# r10 = ugg()
# r10.summation = ["i", "a" , "b", "j"]
# r10.operator_idx.append(["i", "a"])
# r10.operator_idx.append(["j", "b"])
# r10.operator_type.append("s")
# r10.operator_type.append("s")

# r2 = ugg()
# r2.summation = ["c", "k", "d", "l"]
# r2.operator_idx.append(["b", "j"])
# r2.operator_idx.append(["c", "k"])
# r2.operator_idx.append(["d", "l"])
# r2.operator_type.append("s")
# r2.operator_type.append("s")
# r2.operator_type.append("s")

# r = evaluate(r2, r1, r1, r1, r1, r1)
# print('komutator')
# for x in r:
#     print(x)
# sys.exit(0)

# rr1 = r2.fromleft(r1)
# rint = rr1.integrate()
# r = r2
# rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ["s", "s"]).scale(0.5)
# rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], braspin = ["s", "s"])#.scale(0.5)
# rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
# rsimp = simplify(rint1)
# print('calka')
# for x in rsimp:
#     print(x)

# sys.exit(0)
# x = ugg()
# x.summation = ['a', 'i', 'b', 'j']
# x.operator_idx.append(['j','b'])
# x.operator_idx.append(['i','a'])
# x.coefficient = [EOM_CC_AMPLITUDE_L]
# x.coefficient_idx.append(['a','i', 'b', 'j'])
# x.operator_type.append('s')
# x.operator_type.append('s')
# disambiguate(x, t2)
# print(t2)
# r = evaluate(obs, t2)
# r1 = r.fromleft(x)
# for a in r1:
#     print(a)
# rint = r1.integrate()
# for a in rint:
#     print(a)
# rsimp = simplify(rint)
# print('')
# for i in rsimp:
#     print(i)
# sys.exit(0)

# r = arithmetic_string(obs)
# # # r = evaluate(obs, t3)
# # # r = evaluate(fock_imp, t3)
# # # x = ugg()
# # # x.operator_idx.append(['c','k'])
# # # x.operator_type.append('s')
# # # basic_disambiguate(t2c, t2)
# # r1 = r.fromleft(t2c)
# r2 = r.fromright(t2)
# rint1 = r2.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ["s", "s"])
# rint2 = r2.integrate(bra = ['a', 'j', 'b', 'i'], braspin = ["s", "s"])
# rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
# # print(r2)
# # rint = r2.integrate(ket = ['a', 'i'], ketspin = ["s"])
# # # rint.exec_delta()
# rsimp = simplify(rint)
# rint = deepcopy(rsimp)

# for x in rint:
#     print(x)


# r2 = deepcopy(rint)
# for x in r2:
#     x.new_delta("a", "b")
#     x.new_delta("i", "j")
# rint += r2.scale(-1./2.)
# print('')
# for x in rint:
#     print(x)
# sys.exit(0)
# rsipm = simplify(rint)
# print('calka')
# for x in rsimp:
#     print(x)


# sys.exit(0)
# for i in r:
#     print(i)

# rint = r.integrate(bra = ['a', 'i'], braspin = ["s"]).scale(0.5)
# rint = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ["s", "s", "s"])
# rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], braspin = ["s", "s"])
# rint  = rint1.scale(1./3.) + rint2.scale(1./6.)

# rsimp = simplify(rint)
# print('wyn')
# print(fixed)
# for  i in rsimp:
#     print(i)
# sys.exit(0)

# r = evaluate(t2c,t1,t2)
# rint = r.integrate(bra = ['a','i'], braspin = ['s']).scale(0.5)

# rint.exec_delta()
# rsimp = simplify(rint)
# rsimp.cleanup()

# print('wynik')
# for x in rsimp:
#     print(x)
# sys.exit(0)

# r = ugg()
# r.operator_idx.append(["a", "i"])
# r.operator_idx.append(["a", "i"])
# r.operator_type.append("s")
# r.operator_type.append("s")

# rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ["s", "s"])
# rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], braspin = ["s", "s"])
# rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
    
# rsimp = simplify(rint)
# rint = deepcopy(rsimp)

# r2 = deepcopy(rint)
# for x in r2:
#     x.new_delta("a", "b")
#     x.new_delta("i", "j")
# rint += r2.scale(-1./2.)

# rsimp2 = simplify(rint)
# print('wynik')
# for x in rsimp2:
#     print(x)
# sys.exit(0)

# def bra_ort(r, a, i, b, j, c, k):

#     rint1 = triplet_bra_part(r, a, i, b, j, c, k)
#     rint2 = triplet_bra_part(r, a, k, b, j, c, i)
#     rint3 = triplet_bra_part(r, a, j, b, i, c, k)
    
#     rint4 = triplet_bra_part(r, b, i, a, j, c, k)
#     rint5 = triplet_bra_part(r, b, k, a, j, c, i)
#     rint6 = triplet_bra_part(r, b, j, a, i, c, k)
    
#     rint7 = triplet_bra_part(r, c, i, b, j, a, k)
#     rint8 = triplet_bra_part(r, c, k, b, j, a, i)
#     rint9 = triplet_bra_part(r, c, j, b, i, a, k)
    
#     rint0 = rint1.scale(9.0) + rint2 + rint3 +  \
#         rint4 + rint5.scale(-1.0) + rint6.scale(-1.0) + \
#         rint7 + rint8.scale(-1.0) + rint9.scale(-1.0) 
    
#     rint = rint0.scale(1.0/10.0)
#     print('')
#     print('wynik rint')
#     print('')
#     for x in rint:
#         print(x)
#     rsimp = simplify(rint)
#     return rsimp


# def triplet_bra_part(r, a, i, b, j, c, k):
    
#     rinta1 = r.integrate(bra = [a, i, b, j, c, k], braspin = ["s", "t", "s"])
#     rinta2 = r.integrate(bra = [a, i, c, k, b, j], braspin = ["s", "t", "s"])
#     rinta  = rinta1.scale(1.0/8.0) + rinta2.scale(1.0/8.0)
#     # print('AAAAAAAAAAa')
#     # print('wynik rint1')
#     # print('')
#     # for x in rint1:
#     #     x.exec_delta()
#     #     x.standarize_delta()
#     #     if x.num_factor != 0:
#     #         print(x)
#     # sys.exit(0)
#     # print('wynik rint2')
#     # print('')
#     # for x in rint1:
#     #     x.exec_delta()
#     #     x.standarize_delta()
#     #     if x.num_factor != 0:
#     #         print(x)

#     rsimpa = simplify(rinta)

#     return rsimpa

# def nu_triplet(a, i, b, j, c, k):

#     nu1 = operat3([[a, i], [b, j], [c, k]], ["s", "t", "s"])
#     nu2 = operat3([[a, i], [c, k], [b, j]], ["s", "t", "s"])

#     r = arithmetic_string(nu1) + arithmetic_string(nu2)
#     return r
#---------------------------------------------------------------------------------------------


# nu1 = operat2([["a", "i"], ["b", "j"]], ["t", "s"])
# nu2 = operat2([["b", "j"], ["a", "i"]], ["t", "s"])
# #nu2.scale(-1.0)
# r = arithmetic_string(nu1) #+ arithmetic_string(nu2)
# rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ["t", "s"]).scale(1.0/8.0)
# #rint2 = r.integrate(bra = ['b', 'j', 'a', 'i'], braspin = ["t", "s"]).scale(-1.0/8.0)
# rint  = rint1# + rint2
# rsimp = simplify(rint)
# print('wynik')
# for x in rsimp:
#     print(x)
# sys.exit(0)
    

# nu1 = operat3([["d", "l"], ["e", "m"], ["f", "n"]], ["s", "t", "s"])
# nu2 = operat3([["d", "l"], ["f", "n"], ["e", "m"]], ["s", "t", "s"])

# r = evaluate(fock_imp, nu1) + evaluate(fock_imp, nu2)

# rint1 = triplet_bra_part(r, 'a', 'i', 'b', 'j', 'c', 'k')
# rint2 = triplet_bra_part(r, 'a', 'k', 'b', 'j', 'c', 'i')
# rint3 = triplet_bra_part(r, 'a', 'j', 'b', 'i', 'c', 'k')

# rint4 = triplet_bra_part(r, 'b', 'i', 'a', 'j', 'c', 'k')
# rint5 = triplet_bra_part(r, 'b', 'k', 'a', 'j', 'c', 'i')
# rint6 = triplet_bra_part(r, 'b', 'j', 'a', 'i', 'c', 'k')

# rint7 = triplet_bra_part(r, 'c', 'i', 'b', 'j', 'a', 'k')
# rint8 = triplet_bra_part(r, 'c', 'k', 'b', 'j', 'a', 'i')
# rint9 = triplet_bra_part(r, 'c', 'j', 'b', 'i', 'a', 'k')


# rint0 = rint1.scale(9.0) + rint2 + rint3 +  \
#     rint4 + rint5.scale(-1.0) + rint6.scale(-1.0) + \
#     rint7 + rint8.scale(-1.0) + rint9.scale(-1.0) 

# rint = rint0.scale(1.0/10.0)

# rsimp = simplify(rint)

# print('')
# for x in rsimp:
#     print(x)

# sys.exit(0)

#r = evaluate(hamiltoniant, pluszR2p)

# nu1 = operat3([["d", "l"], ["e", "m"], ["f", "n"]], ["s", "t", "s"])
# nu2 = operat3([["d", "l"], ["f", "n"], ["e", "m"]], ["s", "t", "s"])

# r = arithmetic_string(nu1) + arithmetic_string(nu2)


# r = nu_triplet('a', 'i', 'a', 'j', 'c', 'k')
 
# #wynik = bra_ort9(r, 'e', 'n', 'f', 'l', 'g', 'm')
# wynik = bra_ort3_virt(r, 'e', 'n', 'f', 'l', 'g', 'm')
# print('wynik')
# for x in wynik:
#     print(x)
# sys.exit(0)

#r = nu_triplet('a', 'k', 'a', 'j', 'c', 'i')
#r = nu_triplet('a', 'j', 'a', 'i', 'c', 'k')


#r = nu_triplet('c', 'i', 'a', 'j', 'a', 'k')
#r = nu_triplet('c', 'k', 'a', 'j', 'a', 'i')
#r = nu_triplet('c', 'j', 'a', 'i', 'a', 'k')



# rint1 = triplet_bra_part(r, 'a', 'i', 'b', 'j', 'c', 'k')
# rint2 = triplet_bra_part(r, 'a', 'k', 'b', 'j', 'c', 'i')
# rint3 = triplet_bra_part(r, 'a', 'j', 'b', 'i', 'c', 'k')

# rint4 = triplet_bra_part(r, 'b', 'i', 'a', 'j', 'c', 'k')
# rint5 = triplet_bra_part(r, 'b', 'k', 'a', 'j', 'c', 'i')
# rint6 = triplet_bra_part(r, 'b', 'j', 'a', 'i', 'c', 'k')

# rint7 = triplet_bra_part(r, 'c', 'i', 'b', 'j', 'a', 'k')
# rint8 = triplet_bra_part(r, 'c', 'k', 'b', 'j', 'a', 'i')
# rint9 = triplet_bra_part(r, 'c', 'j', 'b', 'i', 'a', 'k')


#rint = triplet_bra_part(r, 'a', 'i', 'a', 'j', 'c', 'k')
#rint = triplet_bra_part(r, 'a', 'k', 'a', 'j', 'c', 'i')
#rint = triplet_bra_part(r, 'a', 'j', 'a', 'i', 'c', 'k')

#rint = triplet_bra_part(r, 'c', 'i', 'a', 'j', 'a', 'k')
#rint = triplet_bra_part(r, 'c', 'k', 'a', 'j', 'a', 'i')
#rint = triplet_bra_part(r, 'c', 'j', 'a', 'i', 'a', 'k')




# print('wynik')
# for x in rint:
#     print(x)
# sys.exit(0)

# rr = []
# rr.append

# rr.append(nu_triplet('a', 'i', 'b', 'j', 'c', 'k'))
# rr.append(nu_triplet('a', 'k', 'b', 'j', 'c', 'i'))
# rr.append(nu_triplet('a', 'j', 'b', 'i', 'c', 'k'))
# rr.append(nu_triplet('b', 'i', 'a', 'j', 'c', 'k'))
# rr.append(nu_triplet('b', 'k', 'a', 'j', 'c', 'i'))
# rr.append(nu_triplet('b', 'j', 'a', 'i', 'c', 'k'))
# rr.append(nu_triplet('c', 'i', 'b', 'j', 'a', 'k'))
# rr.append(nu_triplet('c', 'k', 'b', 'j', 'a', 'i'))
# rr.append(nu_triplet('c', 'j', 'b', 'i', 'a', 'k'))

# for i in range(0, 9):
    
#     r = rr[i]

#     rint1 = triplet_bra_part(r, 'a', 'i', 'b', 'j', 'c', 'k')

#     rsimp1 = simplify(rint1)
#     print('WYNIK', '1', i+1)
#     for x in rsimp1:
#         if len(x.delta) == 0:
#             print(x)
#     print('')

#     rint2 = triplet_bra_part(r, 'a', 'k', 'b', 'j', 'c', 'i')
#     rsimp2 = simplify(rint2)
#     print('WYNIK', '2', i+1)
#     for x in rsimp2:
#         if len(x.delta) == 0:
#             print(x)
#     print('')

#     rint3 = triplet_bra_part(r, 'a', 'j', 'b', 'i', 'c', 'k')
#     rsimp3 = simplify(rint3)
#     print('WYNIK', '3', i+1)
#     for x in rsimp3:
#         if len(x.delta) == 0:
#             print(x)
#     print('')

#     rint4 = triplet_bra_part(r, 'b', 'i', 'a', 'j', 'c', 'k')
#     rsimp4 = simplify(rint4)
#     print('WYNIK', '4', i+1)
#     for x in rsimp4:
#         if len(x.delta) == 0:
#             print(x)
#     print('')

#     rint5 = triplet_bra_part(r, 'b', 'k', 'a', 'j', 'c', 'i')
#     rsimp5 = simplify(rint5)
#     print('WYNIK', '5', i+1)
#     for x in rsimp5:
#         if len(x.delta) == 0:
#             print(x)
#     print('')

#     rint6 = triplet_bra_part(r, 'b', 'j', 'a', 'i', 'c', 'k')
#     rsimp6 = simplify(rint6)
#     print('WYNIK', '6', i+1)
#     for x in rsimp6:
#         if len(x.delta) == 0:
#             print(x)
#     print('')
    
#     rint7 = triplet_bra_part(r, 'c', 'i', 'b', 'j', 'a', 'k')
#     rsimp7 = simplify(rint7)
#     print('WYNIK', '7', i+1)
#     for x in rsimp7:
#         if len(x.delta) == 0:
#             print(x)
#     print('')
#     rint8 = triplet_bra_part(r, 'c', 'k', 'b', 'j', 'a', 'i')
#     rsimp8 = simplify(rint8)
#     print('WYNIK', '8', i+1)
#     for x in rsimp8:
#         if len(x.delta) == 0:
#             print(x)
#     print('')
#     rint9 = triplet_bra_part(r, 'c', 'j', 'b', 'i', 'a', 'k')
#     rsimp9 = simplify(rint9)
#     print('WYNIK', '9', i+1)
#     for x in rsimp9:
#         if len(x.delta) == 0:
#             print(x)
#     print('')

# sys.exit(0)
#---------------------------------------------------------------------------------------------

# r = nu_triplet('a', 'k', 'b', 'j', 'c', 'i')
# rint1 = triplet_bra_part(r, 'a', 'i', 'b', 'j', 'c', 'k')    
# rint2 = triplet_bra_part(r, 'a', 'k', 'b', 'j', 'c', 'i')    
# rint3 = triplet_bra_part(r, 'a', 'j', 'b', 'i', 'c', 'k')    

# rint4 = triplet_bra_part(r, 'b', 'i', 'a', 'j', 'c', 'k')    
# rint5 = triplet_bra_part(r, 'b', 'k', 'a', 'j', 'c', 'i')    
# rint6 = triplet_bra_part(r, 'b', 'j', 'a', 'i', 'c', 'k')    

# rint7 = triplet_bra_part(r, 'c', 'i', 'b', 'j', 'a', 'k')    
# rint8 = triplet_bra_part(r, 'c', 'k', 'b', 'j', 'a', 'i')    
# rint9 = triplet_bra_part(r, 'c', 'j', 'b', 'i', 'a', 'k')    

# rint = rint1.scale(9.0) + rint2 + rint3 + \
#     rint4 + rint5.scale(-1.0) + rint6.scale(-1.0)  + \
#     rint7 + rint8.scale(-1.0) + rint9.scale(-1.0)

# rint0 = rint.scale(1.0/10.0)

# rsimp = simplify(rint0)
# print('wynik')
# for x in rsimp:
#     if len(x.delta) == 0:
#         print(x)
# sys.exit(0)
#---------------------------------------------------------------------------------------------
# rr = []
# rr.append

# rr.append(nu_triplet('a', 'i', 'b', 'i', 'c', 'k'))
# rr.append(nu_triplet('b', 'i', 'a', 'i', 'c', 'k'))
# rr.append(nu_triplet('c', 'i', 'b', 'i', 'a', 'k'))
# # rr.append(nu_triplet('c', 'i', 'a', 'j', 'a', 'k'))
# # rr.append(nu_triplet('c', 'k', 'a', 'j', 'a', 'i'))
# # rr.append(nu_triplet('c', 'j', 'a', 'i', 'a', 'k'))

# for i in range(0, 3):
    
#     r = rr[i]

#     #rint1 = triplet_bra_part(r, 'c', 'k', 'a', 'j', 'a', 'i')    
#     rint1 = triplet_bra_part(r, 'a', 'i', 'b', 'i', 'c', 'k')
#     # rint1 = triplet_bra_part(r, 'a', 'k', 'a', 'j', 'c', 'i')
#     # rint1 = triplet_bra_part(r, 'a', 'i', 'a', 'j', 'c', 'k')
#     # rint1 = triplet_bra_part(r, 'a', 'j', 'a', 'k', 'c', 'i')

#     rsimp1 = simplify(rint1)
#     print('WYNIK', '1', i+1)
#     for x in rsimp1:
#         if len(x.delta) == 0:
#             print(x)
#     print('')

# #    rint2 = triplet_bra_part(r, 'a', 'k', 'a', 'j', 'c', 'i')
# #    rint2 = triplet_bra_part(r, 'c', 'i', 'a', 'j', 'a', 'k')
#     rint2 = triplet_bra_part(r, 'b', 'i', 'a', 'i', 'c', 'k')
#     rsimp2 = simplify(rint2)
#     print('WYNIK', '2', i+1)
#     for x in rsimp2:
#         if len(x.delta) == 0:
#             print(x)
#     print('')

# #    rint3 = triplet_bra_part(r, 'a', 'j', 'a', 'i', 'c', 'k')
# #    rint3 = triplet_bra_part(r, 'c', 'j', 'a', 'k', 'a', 'i')
#     rint3 = triplet_bra_part(r, 'c', 'i', 'b', 'i', 'a', 'k')

#     rsimp3 = simplify(rint3)
#     print('WYNIK', '3', i+1)
#     for x in rsimp3:
#         if len(x.delta) == 0:
#             print(x)
#     print('')
    
##     rint4 = triplet_bra_part(r, 'c', 'i', 'a', 'j', 'a', 'k')
#     rint4 = triplet_bra_part(r, 'a', 'k', 'c', 'j', 'a', 'i')

#     rsimp4 = simplify(rint4)
#     print('WYNIK', '4', i+1)
#     for x in rsimp4:
#         if len(x.delta) == 0:
#             print(x)
#     print('')
# #    rint5 = triplet_bra_part(r, 'c', 'k', 'a', 'j', 'a', 'i')
#     rint5 = triplet_bra_part(r, 'a', 'i', 'c', 'j', 'a', 'k')

#     rsimp5 = simplify(rint5)
#     print('WYNIK', '5', i+1)
#     for x in rsimp5:
#         if len(x.delta) == 0:
#             print(x)
#     print('')
# #    rint6 = triplet_bra_part(r, 'c', 'j', 'a', 'i', 'a', 'k')
#     rint6 = triplet_bra_part(r, 'a', 'j', 'c', 'k', 'a', 'i')
#     rsimp6 = simplify(rint6)
#     print('WYNIK', '6', i+1)
#     for x in rsimp6:
#         if len(x.delta) == 0:
#             print(x)
#     print('')

# sys.exit(0)
#------------------------------------------------------------------------------------------


# rint1 = triplet_bra_part(r, 'a', 'i', 'b', 'j', 'c', 'k')
# print('rint1 completed')
# rint2 = triplet_bra_part(r, 'a', 'k', 'b', 'j', 'c', 'i')
# print('rint2 completed')
# rint3 = triplet_bra_part(r, 'a', 'j', 'b', 'i', 'c', 'k')
# print('rint3 completed')
# rint4 = triplet_bra_part(r, 'b', 'i', 'a', 'j', 'c', 'k')
# print('rint4 completed')
# rint5 = triplet_bra_part(r, 'b', 'k', 'a', 'j', 'c', 'i')
# print('rint5 completed')
# rint6 = triplet_bra_part(r, 'b', 'j', 'a', 'i', 'c', 'k')
# print('rint6 completed')

# rint7 = triplet_bra_part(r, 'c', 'i', 'b', 'j', 'a', 'k')
# print('rint7 completed')
# rint8 = triplet_bra_part(r, 'c', 'k', 'b', 'j', 'a', 'i')
# print('rint8 completed')
# rint9 = triplet_bra_part(r, 'c', 'j', 'b', 'i', 'a', 'k')
# print('rint9 completed')


# print('')
# print('wynik rint1')
# print('')
# for x in rint1:
#     print(x)
# print('')
# print('wynik rint2')
# print('')
# for x in rint2:
#     print(x)
# print('')
# print('wynik rint3')
# print('')
# for x in rint3:
#     print(x)
# #sys.exit(0)

# rint0 = rint1#.scale(11.0) + rint2.scale(-1.0) + rint3.scale(-1.0)

# rint0 = rint1.scale(9.0) + rint2 + rint3 +  \
#     rint4 + rint5.scale(-1.0) + rint6.scale(-1.0) + \
#     rint7 + rint8.scale(-1.0) + rint9.scale(-1.0) 

# rint = rint0.scale(1.0/10.0)
# print('')
# print('wynik rint')
# print('')
# for x in rint:
#     print(x)
# rsimp = simplify(rint)

# print('')
# print('wynik')
# print('')
# for x in rsimp:
#     print(x)

# sys.exit(0)


#------------------------------------TEST 21 cc3----------------------------------------------------------

# r = ugg()
# #r.summation = ["b", "j"]
# #r.coefficient = [CC_AMPLITUDE]
# #r.coefficient_idx.append(["b", "j"])
# r.operator_idx.append(["b", "j"])
# r.operator_type.append("s")
# r.num_factor = 1.0

# rr = arithmetic_string(r)
# rint = rr.integrate(bra = ['a', 'i'], braspin = ["s"])
# rsimp = simplify(rint)
# print('')
# print('wynik')
# print('')
# for x in rsimp:
#     print(x)

# sys.exit(0)


# R2p1 = ugg()
# R2p1.summation = ["a", "i", "b", "j"]
# R2p1.coefficient = [EOM_TRIPLET_R2m]
# R2p1.coefficient_idx.append(["a", "i", "b", "j"])
# R2p1.operator_idx.append(["a", "i"])
# R2p1.operator_idx.append(["b", "j"])
# R2p1.operator_type.append("t")
# R2p1.operator_type.append("s")
# R2p1.num_factor = 1.0#/2.0


# R2p2 = ugg()
# R2p2.summation = ["a", "i", "b", "j"]
# R2p2.coefficient = [EOM_TRIPLET_R2m]
# R2p2.coefficient_idx.append(["a", "i", "b", "j"])
# R2p2.operator_idx.append(["b", "j"])
# R2p2.operator_idx.append(["a", "i"])
# R2p2.operator_type.append("t")
# R2p2.operator_type.append("s")
# R2p2.num_factor = -1.0#/2.0


# r = evaluate(hamiltoniant, R2p1) + evaluate(hamiltoniant, R2p2)
# r = evaluate(hamiltonian, plusz)

# print('')
# print('po evaluate')
# for x in r:
#     print(x)
# print('')


# op = ugg()
# op.operator_idx.append(['i', 'a'])
# op.operator_idx.append(['b', 'j'])
# op.operator_idx.append(['k', 'c'])
# op.operator_idx.append(['d', 'l'])
# op.operator_type.append('t')
# op.operator_type.append('t')
# op.operator_type.append('t')
# op.operator_type.append('t')

# r = arithmetic_string(op)
# rint = r.integrate()
# rsimp = simplify(rint)
# print('wynik')
# for x in rsimp:
#     print(x)
# sys.exit(0)

#1
# nu1 = operat3([["b", "i"], ["a", "j"], ["c", "k"]], ["s", "t", "s"])
# nu2 = operat3([["b", "i"], ["c", "k"], ["a", "j"]], ["s", "t", "s"])
# #2
# nu1 = operat3([["b", "i"], ["a", "k"], ["c", "j"]], ["s", "t", "s"])
# nu2 = operat3([["b", "i"], ["c", "j"], ["a", "k"]], ["s", "t", "s"])
# #3
# nu1 = operat3([["c", "k"], ["b", "j"], ["a", "i"]], ["s", "t", "s"])
# nu2 = operat3([["c", "k"], ["a", "i"], ["b", "j"]], ["s", "t", "s"])
# #4
# nu1 = operat3([["b", "j"], ["a", "k"], ["c", "i"]], ["s", "t", "s"])
# nu2 = operat3([["b", "j"], ["c", "i"], ["a", "k"]], ["s", "t", "s"])
# #5
# nu1 = operat3([["b", "k"], ["a", "i"], ["c", "j"]], ["s", "t", "s"])
# nu2 = operat3([["b", "k"], ["c", "j"], ["a", "i"]], ["s", "t", "s"])
# # #6
# nu1 = operat3([["b", "k"], ["a", "j"], ["c", "i"]], ["s", "t", "s"])
# nu2 = operat3([["b", "k"], ["c", "i"], ["a", "j"]], ["s", "t", "s"])


# r = arithmetic_string(nu1) + arithmetic_string(nu2) 




# nu1 = operat2([["b", "j"], ["c", "k"]], ["t", "s"])
# nu2 = operat2([["c", "k"], ["b", "j"]], ["t", "s"])
#nu2.scale(-1.0)
# r = evaluate(hamiltoniant, nu1) +  evaluate(hamiltoniant, nu2)     



# r = operat2([["c", "k"], ["d", "l"]], ["s", "t"])


# print('po komutatorach')
# print('')
# for x in r:
#     print(x)
# print('')

#rint = r.integrate(bra = ['a', 'i'], braspin = ["t"]).scale(1.0/2.0)
# rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ["t", "s"]).scale(1.0/8.0)
# rint2 = r.integrate(bra = ['b', 'j', 'a', 'i'], braspin = ["t", "s"]).scale(-1.0/8.0)
# rint  = rint1 + rint2


# rint1 = triplet_bra_part(r, 'e', 'n', 'f', 'm', 'g', 'l')
# rint2 = triplet_bra_part(r, 'e', 'l', 'f', 'm', 'g', 'n')
# rint3 = triplet_bra_part(r, 'e', 'm', 'f', 'n', 'g', 'l')

# rint4 = triplet_bra_part(r, 'f', 'n', 'e', 'm', 'g', 'l')
# rint5 = triplet_bra_part(r, 'f', 'l', 'e', 'm', 'g', 'n')
# rint6 = triplet_bra_part(r, 'f', 'm', 'e', 'n', 'g', 'l')

# rint7 = triplet_bra_part(r, 'g', 'n', 'f', 'm', 'e', 'l')
# rint8 = triplet_bra_part(r, 'g', 'l', 'f', 'm', 'e', 'n')
# rint9 = triplet_bra_part(r, 'g', 'm', 'f', 'n', 'e', 'l')
# print('po calkowaniu rint1')
# # print('')
# # s = 0.0
# for x in rint1:
# #     s += x.num_factor()
#     print(x)
# print('')
# # print(s)

# # rint2 = triplet_bra_part(r, 'a', 'k', 'b', 'j', 'c', 'i')
# print('po calkowaniu rint2')
# print('')
# for x in rint2:
#     print(x)
# print('')
# # rint3 = triplet_bra_part(r, 'a', 'j', 'b', 'i', 'c', 'k')
# print('po calkowaniu rint3')
# print('')
# for x in rint3:
#     print(x)
# print('')

# # rint4 = triplet_bra_part(r, 'b', 'i', 'a', 'j', 'c', 'k')
# print('po calkowaniu rint4')
# print('')
# for x in rint4:
#     print(x)
# # print('')
# # rint5 = triplet_bra_part(r, 'b', 'k', 'a', 'j', 'c', 'i')
# print('po calkowaniu rint5')
# print('')
# for x in rint5:
#     print(x)
# print('')
# # rint6 = triplet_bra_part(r, 'b', 'j', 'a', 'i', 'c', 'k')
# print('po calkowaniu rint6')
# print('')
# for x in rint6:
#     print(x)
# print('')

# # rint7 = triplet_bra_part(r, 'c', 'i', 'b', 'j', 'a', 'k')
# print('po calkowaniu rint7')
# print('')
# for x in rint7:
#     print(x)
# print('')
# # rint8 = triplet_bra_part(r, 'c', 'k', 'b', 'j', 'a', 'i')
# print('po calkowaniu rint8')
# print('')
# for x in rint8:
#     print(x)
# print('')
# # rint9 = triplet_bra_part(r, 'c', 'j', 'b', 'i', 'a', 'k')
# print('po calkowaniu rint9')
# print('')
# for x in rint9:
#     print(x)
# print('')


# rint0 = rint1.scale(9.0) + rint2 + rint3 +  \
#     rint4 + rint5.scale(-1.0) + rint6.scale(-1.0) + \
#     rint7 + rint8.scale(-1.0) + rint9.scale(-1.0) 

# print('rint0')
# print('')
# for x in rint0:
#     print(x)
# print('')

# rint = rint0.scale(1.0/80.0)
# print('rint')
# print('')
# for x in rint:
#     print(x)
# print('')

# rsimp = simplify(rint)

# print('')
# print('wynik')
# print('')
# for x in rsimp:
#     print(x)
# sys.exit(0)

# print(rsimp[0].delta[0][1])
# for x in range(0, len(rsimp)):
#     print(rsimp[x].delta[0][1]+rsimp[x].delta[1][1]+rsimp[x].delta[2][1], \
#               rsimp[x].delta[3][1]+rsimp[x].delta[4][1]+rsimp[x].delta[5][1])
#    # print(x)


# r = evaluate(mu1, t1c)

# for x in r:
#     print(x)

# rint = r.integrate(bra = ['a', 'i'], braspin = ['s', 's'])
# for x in rint:
#     print(x)
# sys.exit(0)
# r2 = evaluate(obs, s2c)

# print('r2')
# for x in r2:
#     print(x)

# rint = arithmetic_string()

# for i in range(0, len(r)):
#     for j in range(0, len(r2)):
#         disambiguate(r[i], r2[j], eomr2)
        
#         print('dismabiguate')
#         print(r[i])
#         print(r2[j])
#         print(eomr2)
#         print('')
        

#         r3 = r2[j].fromleft(r[i])
#         r4 = r3.fromright(eomr2)

#         rt = r4.integrate()
#         rt.exec_delta()
#         if len(rt) > 0 :
#             for x in rt:
#                 print(x)
#                 rint.append(x)

# # for x in rint:
# #     print(x)

# r = arithmetic_string(obs)
# r2 = r.fromright(eomr1)
# rint = r2.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s', 's'])
# for i in range(0, len(rint)):
#     rint[i].summation.append('a')
#     rint[i].summation.append('i')
#     rint[i].summation.append('b')
#     rint[i].summation.append('j')
# # rint.clear_fixed()
# rsimp = simplify(rint)
# for x in rsimp:
#     print(x)
# sys.exit(0)


# print(s1c)
# print(eomr2)
# r = evaluate(eomr2, s1c)
# print('r')
# for x in r:
#     print(x)
# rr = r.fromleft(obs)
# eoml1.transpose()
# for x in rr:
#     disambiguate(eoml1, x)
# rr = rr.fromleft(eoml1)
# print('rr')
# for x in rr:
#     print(x)
# print('')
# rint = rr.integrate()#bra  = ['a', 'j'], braspin = ['s', 's'])
# rsimp = simplify(rint)
# print('')
# for x in rsimp:
#     print(x)

# xxx = eij
# r = arithmetic_string(xxx)
# print(xxx)
# disambiguate(xxx, eoml1, eomr1)
# print(xxx, eoml1, eomr1)
# rr = r.fromright(eomr1)
# eoml1.transpose()
# rrr = rr.fromleft(eoml1)
# print(rrr)

# rint = rrr.integrate()

# #rint = r.integrate(bra  = ['a', 'i'], braspin = ['s', 's'], ket  = ['b', 'j'], ketspin = ['s', 's'])
# rsimp = simplify(rint)
# print('')
# for x in rsimp:
#     print(x)


# def gen_ug(a, b, c, i, j, k):

#     ug1 = ugg()
#     ug1.operator_idx.append([k, c])
#     ug1.operator_idx.append([j, b])
#     ug1.operator_idx.append([i, a])
#     ug1.operator_idx.append(['g', 'm'])
#     ug1.operator_idx.append(['f', 'l'])
#     ug1.operator_idx.append(['e', 'n'])
#     ug1.operator_type.append('t')
#     ug1.operator_type.append('s')
#     ug1.operator_type.append('s')
#     ug1.operator_type.append('s')
#     ug1.operator_type.append('t')
#     ug1.operator_type.append('s')


#     ug2 = ugg()
#     ug2.operator_idx.append([j, b])
#     ug2.operator_idx.append([k, c])
#     ug2.operator_idx.append([i, a])
#     ug2.operator_idx.append(['g', 'm'])
#     ug2.operator_idx.append(['f', 'l'])
#     ug2.operator_idx.append(['e', 'n'])
#     ug2.operator_type.append('t')
#     ug2.operator_type.append('s')
#     ug2.operator_type.append('s')
#     ug2.operator_type.append('s')
#     ug2.operator_type.append('t')
#     ug2.operator_type.append('s')

#     ug3 = ugg()
#     ug3.operator_idx.append([k, c])
#     ug3.operator_idx.append([j, b])
#     ug3.operator_idx.append([i, a])
#     ug3.operator_idx.append(['g', 'm'])
#     ug3.operator_idx.append(['e', 'n'])
#     ug3.operator_idx.append(['f', 'l'])
#     ug3.operator_type.append('t')
#     ug3.operator_type.append('s')
#     ug3.operator_type.append('s')
#     ug3.operator_type.append('s')
#     ug3.operator_type.append('t')
#     ug3.operator_type.append('s')


#     ug4 = ugg()
#     ug4.operator_idx.append([j, b])
#     ug4.operator_idx.append([k, c])
#     ug4.operator_idx.append([i, a])
#     ug4.operator_idx.append(['g', 'm'])
#     ug4.operator_idx.append(['e', 'n'])
#     ug4.operator_idx.append(['f', 'l'])
#     ug4.operator_type.append('t')
#     ug4.operator_type.append('s')
#     ug4.operator_type.append('s')
#     ug4.operator_type.append('s')
#     ug4.operator_type.append('t')
#     ug4.operator_type.append('s')

#     return ug1, ug2, ug3, ug4


# ug1, ug2, ug3, ug4 = gen_ug('a', 'a', 'c', 'i', 'j', 'k')
# ug5, ug6, ug7, ug8= gen_ug('a', 'a', 'c', 'k', 'j', 'i')
# ug9, ug10, ug11, ug12 = gen_ug('a', 'a', 'c', 'j', 'i', 'k')

# yg1, yg2, yg3, yg4 = gen_ug('a', 'a', 'c', 'i', 'j', 'k')
# yg5, yg6, yg7, yg8 = gen_ug('a', 'a', 'c', 'k', 'j', 'i')
# yg9, yg10, yg11, yg12 = gen_ug('a', 'a', 'c', 'j', 'i', 'k')

# ig1, ig2, ig3, ig4 = gen_ug('c', 'a', 'a', 'i', 'j', 'k')
# ig5, ig6, ig7, ig8 = gen_ug('c', 'a', 'a', 'k', 'j', 'i')
# ig9, ig10, ig11, ig12 = gen_ug('c', 'a', 'a', 'j', 'i', 'k')

# # ar = arithmetic_string(ug1, ug2, ug3, ug4).scale(9.0) + \
# #      arithmetic_string(ig1, ig2, ig3, ig4) + arithmetic_string(ig9, ig10, ig11, ig12).scale(-1.0) 

# # ar = arithmetic_string(ug1, ug2, ug3, ug4).scale(10.0) + \
# #      arithmetic_string(ig1, ig2, ig3, ig4) + arithmetic_string(ig1, ig2, ig3, ig4).scale(-1.0) \
# #      + arithmetic_string(ig9, ig10, ig11, ig12).scale(-1.0)


# ar = arithmetic_string(ug1, ug2, ug3, ug4) + arithmetic_string(ug5, ug6, ug7, ug8) +arithmetic_string(ug9, ug10, ug11, ug12) 
#      # arithmetic_string(ig1, ig2, ig3, ig4) + arithmetic_string(ig1, ig2, ig3, ig4).scale(-1.0) \
#      # + arithmetic_string(ig9, ig10, ig11, ig12).scale(-1.0)


# # ar = arithmetic_string(ug1, ug2, ug3, ug4).scale(9.0) + arithmetic_string(ug5, ug6, ug7, ug8, ug9, ug10, ug11, ug12) + \
# #      arithmetic_string(yg1, yg2, yg3, yg4) + arithmetic_string(yg5, yg6, yg7, yg8, yg9, yg10, yg11, yg12).scale(-1.0) + \
# #      arithmetic_string(ig1, ig2, ig3, ig4) + arithmetic_string(ig5, ig6, ig7, ig8, ig9, ig10, ig11, ig12).scale(-1.0)

# #ar = ar.scale(0.0125)
# # ug1 = ugg()
# # ug1.operator_idx.append(['k', 'c'])
# # ug1.operator_idx.append(['j', 'b'])
# # ug1.operator_idx.append(['i', 'a'])
# # ug1.operator_idx.append(['d', 'l'])
# # ug1.operator_idx.append(['e', 'm'])
# # ug1.operator_idx.append(['f', 'n'])
# # ug1.operator_type.append('t')
# # ug1.operator_type.append('s')
# # ug1.operator_type.append('s')
# # ug1.operator_type.append('s')
# # ug1.operator_type.append('t')
# # ug1.operator_type.append('s')

# # print(ug1)


# ug1 = ugg()
# ug1.operator_idx.append(['p', 'q'])
# ug1.operator_type.append('t0')

# ug2 = ugg()
# ug2.operator_idx.append(['a', 'i'])
# ug2.operator_type.append('s')


# epqrs1 = ugg()
# epqrs1.operator_idx.append(['p', 'q'])
# epqrs1.operator_idx.append(['r', 's'])
# epqrs1.operator_type.append('s')
# epqrs1.operator_type.append('s')

# epqrs2 = ugg()
# epqrs2.operator_idx.append(['p', 's'])
# epqrs2.operator_type.append('s')
# epqrs2.delta = [['q', 'r']]
# epqrs2.num_factor = -1.0

# tpqrs1t = ugg()
# tpqrs1t.operator_idx.append(['p', 'q'])
# tpqrs1t.operator_idx.append(['r', 's'])
# tpqrs1t.operator_type.append('t0')
# tpqrs1t.operator_type.append('s')

# tpqrs2t = ugg()
# tpqrs2t.operator_idx.append(['p', 's'])
# tpqrs2t.operator_type.append('t0')
# tpqrs2t.delta = [['q', 'r']]
# tpqrs2t.num_factor = -1.0

# tpqrs1p = ugg()
# tpqrs1p.operator_idx.append(['p', 'q'])
# tpqrs1p.operator_idx.append(['r', 's'])
# tpqrs1p.operator_type.append('s')
# tpqrs1p.operator_type.append('t0')

# tpqrs2p = ugg()
# tpqrs2p.operator_idx.append(['p', 's'])
# tpqrs2p.operator_type.append('t0')
# tpqrs2p.delta = [['q', 'r']]
# tpqrs2p.num_factor = -1.0

# tpqrs1 = ugg()
# tpqrs1.operator_idx.append(['p', 'q'])
# tpqrs1.operator_idx.append(['r', 's'])
# tpqrs1.operator_type.append('t0')
# tpqrs1.operator_type.append('t0')

# tpqrs2 = ugg()
# tpqrs2.operator_idx.append(['p', 's'])
# tpqrs2.operator_type.append('s')
# tpqrs2.delta = [['q', 'r']]
# tpqrs2.num_factor = -1.0

# ug3 = ugg()
# ug3.operator_idx.append(['a', 'i'])
# ug3.operator_type.append('s')



# r = evaluate(tpqrs1p, ug3) + evaluate(tpqrs2p, ug3)
# for x in r:
#     print(x)
# sys.exit(0)


# ug3 = ugg()
# ug3.operator_idx.append(['j', 'b'])
# ug3.operator_idx.append(['i', 'a'])
# ug3.operator_idx.append(['d', 'l'])
# ug3.operator_idx.append(['c', 'k'])
# ug3.operator_type.append('s')
# ug3.operator_type.append('t')
# ug3.operator_type.append('t')
# ug3.operator_type.append('s')
# ug3.num_factor = -1.0


# ug4 = ugg()
# ug4.operator_idx.append(['i', 'a'])
# ug4.operator_idx.append(['j', 'b'])
# ug4.operator_idx.append(['d', 'l'])
# ug4.operator_idx.append(['c', 'k'])
# ug4.operator_type.append('s')
# ug4.operator_type.append('t')
# ug4.operator_type.append('t')
# ug4.operator_type.append('s')
# ug4.num_factor = -1.0


# ar = arithmetic_string(ug1, ug2, ug3, ug4).scale(1.0/8.0)
# for x in ar:
#     print('ar', x)

# print('integrate')
# rint = ar.integrate()
# print('po integrate')
# for x in rint:
#     print(x)
# print('simplify')
# rsimp = simplify(rint)
# print('po simplify')
# for x in rsimp:
#     print('a', x)


# r = eomr3
# r = r.fromleft(obs)
# mu2.transpose()
# disambiguate(mu2, r)
# r = r.fromleft(mu2)
# print(r)
# rint = r.integrate()
# rs = simplify(rint)
# print('wynik')
# for x in rs:
#     if x.num_factor != 0:
#         print(x)


# r = evaluate(t2c, t2, t2)

# tplus = ugg()
# tplus.operator_idx.append(['j','a'])
# tplus.operator_type.append("t1")

# tminus = ugg()
# tminus.operator_idx.append(['i','b'])
# tminus.operator_type.append("tm1")

# print(tminus)
# for x in r:
#     disambiguate(tminus, tplus, x)
#     print(tminus, tplus, x)
# print('')

# r2 = r.fromleft(tminus)
# r3 = r2.fromleft(tplus)

# # for x in r3:
# #     print(x)
# # sys.exit(0)
# rint = r3.integrate()
# rint = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s','s'])
# rsimp = simplify(rint)
# print('wynik')
# for x in rsimp:
#     print(x)
# sys.exit(0)

# nu1 = operat2([["c", "k"], ["d", "l"]], ["s", "t0"])
# print(nu1)
# nu2 = operat2([["d", "l"], ["c", "k"]], ["t0", "s"])
# r = arithmetic_string(nu1, nu2)
# rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ["t0", "s"]).scale(1.0/8.0)
# rint2 = r.integrate(bra = ['b', 'j', 'a', 'i'], braspin = ["t0", "s"]).scale(1.0/8.0)
# rint  = rint1 + rint2
# print('pluszo')
# rsimp = simplify(rint)
# for x in rsimp:
#     print(x)
# sys.exit(0)


#nieb = ugg()
# nieb.summation = ['l','m','n','d','i','e','j','c','k']
# nieb.coefficient = [TWOEL_INT_TRANS]
# nieb.coefficient.append(EOM_TRIPLET_R3)
# nieb.coefficient_idx.append(['l','m','n','e'])
# nieb.coefficient_idx.append(['d', 'i', 'e', 'j', 'c', 'k' ])
#nieb.operator_idx.append(['a', 'i'])
# nieb.operator_idx.append(['i>', 'a'])
# nieb.operator_idx.append(['l', 'm'])
# nieb.operator_idx.append(['d', 'i'])
# nieb.operator_idx.append(['n', 'j'])
# nieb.operator_idx.append(['c', 'k'])
#nieb.operator_type.append('s')
# nieb.operator_type.append('t0')
# nieb.operator_type.append('s')
# nieb.operator_type.append('s')
# nieb.operator_type.append('t0')
# nieb.operator_type.append('s')
# nieb.num_factor = 0.25
# print('nieb', nieb)
# rint = nieb.integrate()
# print(rint)
# print('nieb', nieb)
# print('wynik')
# for x in rint:
#     print(x)

# sys.exit(0)

#nieb1 = ugg()
# nieb1.operator_idx.append(['j', 'b'])                                                                                                                                 
# nieb1.operator_idx.append(['i', 'a'])                                                                                                                                 
# nieb1.operator_idx.append(['c', 'k'])
#nieb1.operator_idx.append(['a', 'i'])
#nieb1.operator_type.append('s')
# nieb1.operator_type.append('t0')
# nieb1.operator_type.append('t0')    
# nieb1.operator_type.append('s')   

#nieb2 = ugg()
#nieb2.operator_idx.append(['b', 'j'])
#nieb2.operator_type.append('t0')


#r = evaluate(hamiltoniant, nieb1, nieb2)
#print('wynik')
#for x in r:
#    print(x)
#sys.exit(0)

# nieb2 = ugg()
# nieb2.operator_idx.append(['i', 'a'])
# nieb2.operator_idx.append(['j', 'b'])
# nieb2.operator_idx.append(['c', 'k'])
# nieb2.operator_idx.append(['d', 'l'])
# nieb2.operator_type.append('s')
# nieb2.operator_type.append('t0')
# nieb2.operator_type.append('t0')
# nieb2.operator_type.append('s')

#nieb2.summation = ['l','e', 'f']
# nieb2.coefficient = [TWOEL_INT_TRANS]
# nieb2.coefficient.append(EOM_TRIPLET_R3)
# nieb2.coefficient_idx.append(["l", "e", 'e', 'f'])
# nieb2.coefficient_idx.append(['f', 'i', 'b', 'l', 'a', 'j' ])
# nieb2.num_factor = 0.625

# r = arithmetic_string(nieb, nieb2)

# for x in r:
#     print(x)

# rsimp = simplify(r)
# print(r)
# a, b = nieb.swap(2)
# print('a', a)
# print('b', b)
#sys.exit(0)
#print('nieb', nieb)
#print('')
#r = nieb.integrate()
#print('po integrate')
#for x in r:
#    print(x)
# a, b = nieb.ovsplit()
# print(a, b)
# print('nieb', nieb)

# rint1 = nieb1.integrate().scale(1.0/8.0)#bra = ['a', 'i', 'b', 'j'], braspin = ["t0", "s"]).scale(1.0/8.0)
# rint2 = nieb2.integrate().scale(-1.0/8.0)#bra = ['b', 'j', 'a', 'i'], braspin = ["t0", "s"]).scale(-1.0/8.0)
# print('rint1')
# for x in rint1:
#      print(x)
# print('rint2')
# for x in rint2:
#      print(x)

# rint  = rint1 + rint2
# print('wynik')
# for x in rint:
#     print(x)
# rsimp = simplify(rint)
# print('ssimp')
# for x in rsimp:
#     print(x)
# sys.exit(0)


#------------------------------------CI SINGLET----------------------------------------------------------
#execute_cisd(2, 2)
# sys.exit(0)


#------------------------------------Density Matrix Ground----------------------------------------------------------
# print('EXECUTplusz density matrix ground')
# execute_density_matrix_ground_ground('ccsd', 0)
# sys.exit(0)


# ------------------------------------EOM TRIPLET----------------------------------------------------------
# print('eom_triplet pika')
# task_eom_triplet(3, 2,'cc3','trans', 'load',  lmfold = "", rmfold = "m")
# sys.exit(0)

# ------------------------------------ci TRIPLET----------------------------------------------------------

# r = arithmetic_string(hamiltonian)
# rint = r.integrate(bra = ['a', 'i'], braspin = ["t"], ket = ['b', 'j'], ketspin = ["t"]).scale(0.5)
# rsimp = simplify(rint)

# print('')
# for x in rsimp:
#     print(x)
# sys.exit(0)

# r = ugg()
# r.operator_idx.append(['m','n'])
# r.operator_type.append("s")
# r.operator_idx.append(['p','q'])
# r.operator_type.append("t0")
# print('', r)
# a, b = r.swap(1)
# print(a, b)
# sys.exit(0)
# print(fixed)
# r = evaluate(hamiltonian, nubjckdl)
# rint = r.integrate(bra = ['a', 'i'], braspin=['s']).scale(0.5)
# ars = preprep_for_fortran(rint)
# print('lal')
# for x in rint:
#     print(x)
# sys.exit(0)

# print(fixed)
# for x in range(0, len(ars)):
#     ars[x].coefficient.append(EOM_CC_SINGLE_Rl)
#     ars[x].coefficient_idx.append(['a', 'i'])
#     ars[x].summation.append('a')
#     ars[x].summation.append('i')
#     ars[x].exec_delta()
#     print(ars[x])
#     # rint[x].summation.append('a')
#     # rint[x].summation.apppend('i')
#     # print(rint[x])

# sys.exit(0)
#r = evaluate(hamiltoniant, t2, nudl)
# r = evaluate(hamiltoniant, nudlem)
# rint = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's'])
# rsimp = simplify(rint)
# delta_subst = deepcopy(delta_to_dict(['j', 'k']))
# rsimp.exec_delta_fixed(delta_subst)
# rsimp2 = simplify_fort(rsimp)
# for x in rsimp2:
#     print(x)


# r = evaluate(hamiltoniant, nuck)
# rint = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s', 's'])
# rsimp = simplify(rint)
# for x in rsimp:
#     print(x)
# sys.exit(0)

#temper()
#sys.exit(0)
#pluszaplusza
#------------------------------------EOM----------------------------------------------------------
# print('------------------------EOM-------------------------------------')
# task_eom(1, 3,'cc3', 'trans', 'dump')
# sys.exit(0)
#------------------------------------eom-MEM------------------------------------------------------                                                                     
# pika
# task_eom_mem(3, 2,'cc3', 'trans', 'load', noR = False, restricted = True)
# sys.exit(0)

#------------------------------------TRANSITION MOMENTS----------------------------------------------------------
# maxcluster = 10
# execute_transition_gamma_dm('ccsd', 'dump', mbpt = 3, maxcluster=maxcluster, multiplicity = 1, cumulative = False)
# #execute_transition_xi_dm('ccsd', 'load')
# sys.exit(0)
#------------------------------------TRANSITION MOMENTS EXCITED----------------------------------------------------------
# REMEMBER THAT MAXMBPT is the ONLY MBPT order included. If you want lower orders, you need to compute them separately.
# pikaczu
# execute_transition_exc(maxexc=2, maxmu=2, maxcom=3, maxmbpt = 3, theory='ccsd', pick='load', multiplicity = 1)
#execute_transition_exc(maxexc=3, maxmu=3, maxcom=3, maxmbpt = 4, theory='cc3', maxop = 6, pick='load', multiplicity = 1)
# pika

#ccsd triplet
#execute_transition_exc(maxexc=2, maxmu=2, maxcom=3, maxmbpt=3, theory='ccsd', maxop = 6, pick='dump', multiplicity = 13, cumulative=True)
#execute_transition_exc(maxexc=2, maxmu=2, maxcom=3, maxmbpt=3, theory='ccsd', maxop = 6, pick='dump', multiplicity = 3, cumulative=True)
# execute_transition_exc(maxexc=2, maxmu=2, maxcom=3, maxmbpt=3, theory='ccsd', maxop = 6, pick='load', multiplicity = 3, cumulative=True)
#sys.exit(0)
#poke

#maximum sum of excitations of all cluster operators
#------------------------------------------------------------------------------------------------------------------
maxcluster = 10
#pikapikapika
#for pt in range(2, 3):
     #   print('pt = ', pt)
     # print('cc3')
#     execute_transition_exc(maxmbpt=pt, maxcluster=maxcluster, theory='ccsd', pick='dump', multiplicity = 13, cumulative=False)
     # print('cc3')
#     execute_transition_exc(maxmbpt=pt, maxcluster=maxcluster, theory='ccsd', pick='load', multiplicity = 13, cumulative=False)

    # print('ccsd trip')
#     execute_transition_exc(maxmbpt=pt, maxcluster=maxcluster, theory='cc3', pick='dump', multiplicity = 3, cumulative=False)
   # print('load trip')
#     execute_transition_exc(maxmbpt=pt, maxcluster=maxcluster, theory='cc3', pick='load', multiplicity = 3, cumulative=False)


   #  print('ccsd triplet')
   #  execute_transition_exc(maxexc=2, maxmu=2, maxcom=maxcom, maxmbpt=pt, theory='ccsd', maxcluster = 6, pick='dump', multiplicity = 3, cumulative=False)
   #  print('load')
   #  execute_transition_exc(maxexc=2, maxmu=2, maxcom=maxcom, maxmbpt=pt, theory='ccsd', maxcluster = 6, pick='load', multiplicity = 3, cumulative=False)
   #  print('cc3')
   #  start = time.time()
   #  execute_transition_exc(maxmbpt=pt, maxcluster=maxcluster, theory='cc3', pick='dump', multiplicity = 13, cumulative=False)
   #  end = time.time()
   #  print('time elapsed NA CALE CC3=', end-start)
   #  print('load')
   #  execute_transition_exc(maxmbpt=pt, maxcluster=maxcluster, theory='cc3', pick='load', multiplicity = 13, cumulative=False)
   #  print('cc3 triplet')
   #  execute_transition_exc(maxexc=3, maxmu=3, maxcom=maxcom, maxmbpt=pt, theory='cc3', maxcluster = 6, pick='dump', multiplicity = 3, cumulative=False)
   #  print('load')
   #  execute_transition_exc(maxexc=3, maxmu=3, maxcom=maxcom, maxmbpt=pt, theory='cc3', maxcluster = 6, pick='load', multiplicity = 3, cumulative=False)
#---------------------------------------------------------------------------------------
#sys.exit(0)
# print('dump overlap_ccsd...')
# execute_transition_exc(maxexc=2, maxmu=2, maxcom=2, maxmbpt=40, theory='overlap_ccsd', maxop = 6, pick='dump', multiplicity = 1, cumulative=True)
# sys.exit(0)
# print('load overlap_ccsd...')
# execute_transition_exc(maxexc=2, maxmu=2, maxcom=2, maxmbpt=40, theory='overlap_ccsd', maxop = 6, pick='load', multiplicity = 1, cumulative=True)
# sys.exit(0)
# print('dump overlap_ccsd triplet...')
# execute_transition_exc(maxexc=2, maxmu=2, maxcom=2, maxmbpt=45, theory='overlap_ccsd', maxop = 6, pick='dump', multiplicity = 3, cumulative=True)

# print('load overlap_ccsd triplet...')
# execute_transition_exc(maxexc=2, maxmu=2, maxcom=2, maxmbpt=45, theory='overlap_ccsd', maxop = 6, pick='load', multiplicity = 3, cumulative=True)
# sys.exit(0)
# print('dump overlap_cc3...')
# execute_transition_exc(maxexc=3, maxmu=3, maxcom=2, maxmbpt=40, theory='overlap_cc3', maxop = 6, pick='dump', multiplicity = 1, cumulative = True)
# sys.exit(0)
# print('load overlap_cc3...')
# execute_transition_exc(maxexc=3, maxmu=3, maxcom=2, maxmbpt=40, theory='overlap_cc3', maxop = 6, pick='load', multiplicity = 1, cumulative = True)
# sys.exit(0)
# print('dump overlap_cc3 triplet...')
# execute_transition_exc(maxexc=3, maxmu=3, maxcom=2, maxmbpt=40, theory='overlap_cc3', maxop = 6, pick='dump', multiplicity = 3, cumulative = True)
# print('load overlap_cc3 triplet...')
# execute_transition_exc(maxexc=3, maxmu=3, maxcom=2, maxmbpt=40, theory='overlap_cc3', maxop = 6, pick='load', multiplicity = 3, cumulative = True)
# sys.exit(0)

        # execute_transition_exc(maxexc=2, maxmu=2, maxcom, maxmbpt, theory='ccsd', maxop = 6, pick='dump', multiplicity = 3)
        # execute_transition_exc(maxexc=2, maxmu=2, maxcom, maxmbpt, theory='ccsd', maxop = 6, pick='load', multiplicity = 3)
# execute_transition_exc(maxexc=3, maxmu=3, maxcom=3, maxmbpt = 4, theory='cc3', maxop = 6, pick='load', multiplicity = 3)
#sys.exit(0)
#------------------------------------OMEGA----------------------------------------------------------    
#generate_omega_operators('ccsd', 2)  
#sys.exit(0)  
#------------------------------------TESTY----------------------------------------------------------
# r = evaluate(obs, t2)
# rint = r.integrate(bra = ['a','i', 'b', 'j'])
# for x in rint:
#     print(x)
# sys.exit(0)
# tauaibjck = ugg()
# tauaibjck.operator_idx = [["a", "i"], ["b", "j"], ['c', 'k']]
# r2 = evaluate(s2c, tauaibjck)
# print('[S2c, tau3]')
# for x in r2:
#     print(x)
# print('')
# r = arithmetic_string(obs) * r2
# for x in r:
#     print(x)
# rint = r.integrate()
# rint.exec_delta()
# rsimp = simplify(rint)
# print('')
# print('po symplifikacji')
# for x in rsimp:
#     print(x)
# sys.exit(0)

# #r = evaluate(t1c, t1, t3)
# #r = evaluate(t1c, t2, t2)
# #r = evaluate(t1c, t3, t1)
# #r1 = evaluate(obs, t1, t2)
# r = evaluate(obs, s3c)
# r = arithmetic_string(r[0])
# print('po evaluate')
# for x in r:
#     print(x)
# print('')
# #sys.exit(0)

# rint = r.integrate(ket = ['a','i', 'b', 'j', 'c', 'k'])
# print('po integrate')
# for x in rint:
#     print(x)
# rsimp = simplify(rint)
# for x in rsimp:
#     print(x)
# rsimp.cleanup()
# for x in rsimp:
#     print(x)
# sys.exit(0)

# rint.exec_delta()
# rsimp = simplify(rint)
# rsimp.cleanup()
# rsimp_clear = rsimp.clear_deltas()

# for x in rsimp:
#     print(x)
# sys.exit(0)



#----------------------------------------                                                                                                               
# CCSD-F12                                                                                                                    
#----------------------------------------                                                                                          
# print('Execute CCSD-F12')
# execute_ccsd_f12()
# sys.exit(0)

#----------------------------------------                                                                                                               
# CCSD                                                                                                                                                 
#----------------------------------------
print('Execute CCSD')
execute_ccsd()
sys.exit(0)

#----------------------------------------
# CC3
#----------------------------------------

#r = evaluate(hamiltonian, t3)
#execute_cc3(r, 100)
#r = evaluate(hamiltonian, t3) + evaluate(hamiltonian, t1, t3)
#execute_cc3(r, 101)

#r = evaluate(flukt_potential, t2) + evaluate(flukt_potential, t1, t2) + evaluate(flukt_potential, t1, t1, t2).scale(0.5) \
#     + evaluate(flukt_potential, t1, t1, t1, t2).scale(1/6)

# r = evaluate(flukt_potentialt, t2)
# execute_cc3(r, 86)
# sys.exit(0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TESTY2~~~~~~~~~~~~~~~~~~~~~~~~

# Equations from "Unitary Group Approach to Spin-Adapted Open-Shell
# Coupled Cluster Theory" Jeziorski, Jankowski, Paldus
# Y. 1995, 56, p. 135 eq. (15-19)
# 1 - E(ijk, abc )
# 2 - E(jki, abc )
# 3 - E(kji, abc )
# 4 - E(jik, abc )
# 5 - E(ikj, abc )
# 6 - E(kij, abc )


def intg2(r, n):
    if n == 1:
        rint1 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k']).scale(1./4.)
        rint2 = r.integrate(bra = ['a', 'k', 'b', 'i', 'c', 'j']).scale(1./12.)
        rint3 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i']).scale(1./6.)
        rint4 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k']).scale(1./6)
        rint5 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i']).scale(1./12.)
        r1 = rint1 + rint2 + rint3 + rint4 + rint5
        return r1

def intg(r, n):
    
    if n == 1:
        rint1 = r.integrate(bra = ['a','i', 'b', 'k', 'c', 'j']).scale(1./4.)
        rint2 = r.integrate(bra = ['a','j', 'b', 'i', 'c', 'k']).scale(1./12.)
        rint3 = r.integrate(bra = ['a','j', 'b', 'k', 'c', 'i']).scale(1./6.)
        rint4 = r.integrate(bra = ['a','k', 'b', 'i', 'c', 'j']).scale(1./6.)
        rint5 = r.integrate(bra = ['a','k', 'b', 'j', 'c', 'i']).scale(1./12.)
        rint = rint1 + rint2 + rint3 + rint4 + rint5
    elif n == 2:
        rint1 = r.integrate(bra = ['a','i', 'b', 'k', 'c', 'j']).scale(1./12.)
        rint2 = r.integrate(bra = ['a','j', 'b', 'i', 'c', 'k']).scale(1./4.)
        rint3 = r.integrate(bra = ['a','j', 'b', 'k', 'c', 'i']).scale(1./6.)
        rint4 = r.integrate(bra = ['a','k', 'b', 'i', 'c', 'j']).scale(1./6.)
        rint5 = r.integrate(bra = ['a','k', 'b', 'j', 'c', 'i']).scale(1./12.)
        rint = rint1 + rint2 + rint3 + rint4 + rint5
    elif n == 3:
        rint1 = r.integrate(bra = ['a','i', 'b', 'k', 'c', 'j']).scale(1./6.)
        rint2 = r.integrate(bra = ['a','j', 'b', 'i', 'c', 'k']).scale(1./6.)
        rint3 = r.integrate(bra = ['a','j', 'b', 'k', 'c', 'i']).scale(1./3.)
        rint4 = r.integrate(bra = ['a','k', 'b', 'i', 'c', 'j']).scale(1./6.)
        rint5 = r.integrate(bra = ['a','k', 'b', 'j', 'c', 'i']).scale(1./6.)
        rint = rint1 + rint2 + rint3 + rint4 + rint5
    elif n == 4:
        rint1 = r.integrate(bra = ['a','i', 'b', 'k', 'c', 'j']).scale(1./6.)
        rint2 = r.integrate(bra = ['a','j', 'b', 'i', 'c', 'k']).scale(1./6.)
        rint3 = r.integrate(bra = ['a','j', 'b', 'k', 'c', 'i']).scale(1./6.)
        rint4 = r.integrate(bra = ['a','k', 'b', 'i', 'c', 'j']).scale(1./3.)
        rint5 = r.integrate(bra = ['a','k', 'b', 'j', 'c', 'i']).scale(1./6.)
        rint = rint1 + rint2 + rint3 + rint4 + rint5
    elif n == 5: 
        rint1 = r.integrate(bra = ['a','i', 'b', 'k', 'c', 'j']).scale(1./12.)
        rint2 = r.integrate(bra = ['a','j', 'b', 'i', 'c', 'k']).scale(1./12.)
        rint3 = r.integrate(bra = ['a','j', 'b', 'k', 'c', 'i']).scale(1./6.)
        rint4 = r.integrate(bra = ['a','k', 'b', 'i', 'c', 'j']).scale(1./6.)
        rint5 = r.integrate(bra = ['a','k', 'b', 'j', 'c', 'i']).scale(1./4.)
        rint = rint1 + rint2 + rint3 + rint4 + rint5
        
    return rint



# rsimp2 = deepcopy(rsimp)
# for x in range(0,len(rsimp)):
#     rsimp2[x].coefficient = []
#     rsimp2[x].coefficient_idx = []
#     for i in range (0, len(rsimp[x].coefficient)):
#         if rsimp[x].coefficient[i] == TEMP1:
#             a = rsimp[x].coefficient_idx[i][0]
#             c = rsimp[x].coefficient_idx[i][1]
#         elif rsimp[x].coefficient[i] == TEMP2:
#             b = rsimp[x].coefficient_idx[i][0]
#             d = rsimp[x].coefficient_idx[i][1]
#         else:
#             rsimp2[x].coefficient.append(rsimp[x].coefficient[i])
#             rsimp2[x].coefficient_idx.append(rsimp[x].coefficient_idx[i])
#     rsimp2[x].coefficient.append('w')
#     rsimp2[x].coefficient_idx.append([a, b, c, d])

print('rsimp2')
for x in rsimp2:
    print(x)
print('')
rsimp3 = simplify(rsimp2)
for x in rsimp3:
    print(x)
print('')
sys.exit(0)
rsimp4 = arithmetic_string()

for x in rsimp3:
    if x.summation == ['a','b','c','i','j','k']:
        rsimp4.append(x)

for x in rsimp4:
    print(x)
print('')

sys.exit(0)


#----------------------BLOKI <SHS>, <DHS>, <DHD>, <S[H, T2]S>-----------
#r = arithmetic_string(hamiltonian)
# r = evaluate(hamiltonian, t2)
# rint = r.integrate(bra = ['a', 'i'], ket=['b', 'j'], braspin = ["s"], ketspin = ["s"]).scale(0.5)
# rsimp = simplify(rint)
# print('wynik')
# for x in rsimp:
#     print(x)
# sys.exit(0)

# r = arithmetic_string(hamiltonian)
# rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], ket=['c', 'k'], braspin = ["s", "s"], ketspin = ["s"])
# rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], ket=['c', 'k'], braspin = ["s", "s"], ketspin = ["s"])
# rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
# rsimp = simplify(rint)
# print('wynik')
# for x in rsimp:
#     print(x)
# sys.exit(0)


# r = arithmetic_string(hamiltonian)
# rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], ket=['c', 'k', 'd', 'l'], braspin = ["s", "s"], ketspin = ["s", "s"])
# rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], ket=['c', 'k', 'd', 'l'], braspin = ["s", "s"], ketspin = ["s", "s"])
# rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
# rsimp = simplify(rint)
# print('wynik')
# for x in rsimp:
#     print(x)
# sys.exit(0)

# r = arithmetic_string(hamiltonian)
# rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], ket=['a', 'i', 'b', 'j'], braspin = ["s", "s"], ketspin = ["s", "s"])
# rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], ket=['a', 'j', 'b', 'i'], braspin = ["s", "s"], ketspin = ["s", "s"])
# rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
# rsimp = simplify(rint)
# print('wynik')
# for x in rsimp:
#     print(x)
# sys.exit(0)

#----------------------------------------------------------------------------------------------------------------------------------------
