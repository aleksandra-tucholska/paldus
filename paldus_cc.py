from params import *
import paldus_classes
from paldus_basic import *
from paldus_classes import ugg
from copy import deepcopy
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
from eomccjac import jacobian_loop


# def nevpt2():

#     # r = ham
#     # rint  = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s', 's'])

#     r = evaluate(ham, ebj, eai)

#     rint = integrate(r)

    
#     print('wynik2 z ')
#     for x in rint:
#         print(x)
#     print('')



def execute_michal2():

    obs1 = ugg()
    obs1.summation = ['p','q']
    obs1.coefficient = [OBSERVABLE_X]
    obs1.coefficient_idx.append(['p','q'])
    obs1.operator_idx.append(['p','q'])
    obs1.operator_type.append("s")

    obs2 = deepcopy(obs1)
    obs3 = deepcopy(obs1)
    obs4 = deepcopy(obs1)


    rint = integrate(obs1)
    print('wynik z', obs1)
    print(rint)
    print('')
    print(obs1)
    print(obs2)

    disambiguate(obs1, obs2)
    print(obs2)

    r = obs1
    r = r.fromright(obs2)
    r2 = deepcopy(r)
    print(r)
    rint = integrate(r)
    print('wynik2 z ', r)
    print(rint)
    print('')

    disambiguate(r, obs3)
    print(r)
    print(obs3)


    r2 = r2.fromright(obs3)
    r3 = deepcopy(r2)
    print(r)
    rint = integrate(r2)
    print('wynik3 z ', r2)
    print(rint)
    print('')

    disambiguate(r2, obs4)
    print(r2)
    print(obs4)

    r3 = r3.fromright(obs4)
    print(r)
    rint = integrate(r2)
    print('wynik4 z ',r3)
    print(rint)
    print('')

    sys.exit(0)

    

def execute_michal(n):

    # ket
    ket = t1
    for i in range(2, n+1):
        a = deepcopy(t1)
        disambiguate(a, ket)
        ket = ket.fromright(a)

    bra = t1c
    for i in range(2, n+1):
        b = deepcopy(t1c)
        disambiguate(b, bra)
        bra = bra.fromright(b)

    #    for x in bra:
    disambiguate(bra, ket)

    rr = arithmetic_string()
    r = hamiltonian.fromright(ket)
    for x in r:
        rr.append(bra.fromright(x))
    for x in rr:
        print(x)
    r = rr



    
    # a = deepcopy(t1c)
    # disambiguate(a, t1)
    # p = a.fromright(t1)

    # r = arithmetic_string(p)
    rint = r.integrate()
    rint.exec_delta()
    rsimp  = simplify(rint)
    rsimp2  = simplify(rsimp)
    print('')
    print('wynik <T1^n|T1^n>, n=', n)
    for x in rsimp2:
        print(x)
    print('')
    sys.exit(0)


def execute_polnon_72_73():

    # rownanie 72 a wrtykule jeziorski explicitly
    p1 = deepcopy(g_1).fromright(t2)
    p2 = deepcopy(g_2).fromright(t2)
    p3 = deepcopy(g_1).fromleft(t2c)
    p4 = deepcopy(g_2).fromleft(t2c)

    p6 = evaluate(fock_imp, t2)
    for x in p6:
        disambiguate(x, t2c)

    p7 = deepcopy(p6).fromleft(t2c)

    r = arithmetic_string(p1) + arithmetic_string(p2) + arithmetic_string(p3) + arithmetic_string(p4) + arithmetic_string(p7)
    rint = r.integrate()
    rint.exec_delta()
    rsimp  = simplify(rint)
    rsimp2  = simplify(rsimp)
    print('')
    print('wynik rownania 72')
    for x in rsimp2:
        print(x)
    print('')
    sys.exit(0)
    # rownanie 73 a wrtykule jeziorski explicitly

    w1 = evaluate(g_1, t2)
    w2 = evaluate(g_2, t2)
    for x in w1:
        disambiguate(x, t2c)
    for x in w2:
        disambiguate(x, t2c)

    ww1 = deepcopy(w1).fromleft(t2c)
    ww2 = deepcopy(w2).fromleft(t2c)
        
    r = arithmetic_string(ww1) + arithmetic_string(ww2)
    rint = r.integrate()
    rint.exec_delta()
    rsimp  = simplify(rint)
    rsimp2  = simplify(rsimp)
    rsimp3  = simplify(rsimp2)
    print('')
    print('wynik rownania 73', len(rsimp3))

    for x in rsimp3:
        print(x)
    print('')
    

    sys.exit(0)


def execute_ccsd():
    
    # part for Polnon
    # r = hamiltonian + evaluate(hamiltonian, t1) + evaluate(hamiltonian, t2) + evaluate(hamiltonian, t1, t1).scale(0.5)
    # rint = r.integrate()
    # rint.exec_delta()
    # rsimp  = simplify(rint)
    # print('wynik')
    # for x in rsimp:
    #     print(x)
    # sys.exit(0)

    r1 = hamiltoniant + evaluate(hamiltoniant, t2)
    r2 = hamiltoniant + evaluate(hamiltoniant, t2) +evaluate(hamiltoniant, t2, t2).scale(0.5)

    # part for Polnon
    # r2 = hamiltonian  + evaluate(hamiltonian, t2) +evaluate(hamiltonian, t2, t2).scale(0.5)
    
    #
    # 'Molecular electronic structure theory' vol 2,  Trygve Helgaker
    # p. 692 (170 evince) eq 13.7.58 and 13.7.59
    # Omega vector function in biorthogonal basis.
    #
    # -------------------------------------------------
    print('r1')
    rint1 = r1.integrate(bra = ['a', 'i'], braspin = ['s']).scale(0.5)

    # part for Polnon
    # print('r20')
    # rint2  = r2.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s', 's'])

    
    print('r21')
    rint21  = r2.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s', 's']).scale(1.0/3.0)
    print('r22')
    rint22 = r2.integrate(bra = ['a', 'j', 'b', 'i'], braspin = ['s', 's']).scale(1.0/6.0)
    rint2 = (rint21 + rint22)

    rint1.exec_delta()
    rint2.exec_delta()
    
    rsimp1 = simplify(rint1)
    rsimp2 = simplify(rint2)

    for x in rsimp1:
        x.optimize()
    for x in rsimp2:
        x.optimize()
    #
    # Fixme: Add tunable elimination of non-diagonal fock elements
    # #
    # # DELETE ALL TERMS WITH NONDIAGONAL FOCK MATRIX ELEMENTS
    # # 
    for x in rsimp:
        for y in range(0, len(x.coefficient)):
            if x.coefficient[y] == FOCK_MATRIX:
                    x.num_factor = 0.0

    rsimp1.cleanup()
    rsimp2.cleanup()

    print('wynik')
    for x in rsimp1:
        print(x)

    function_template_ccsd(rsimp1, 1)
    function_template_ccsd(rsimp2, 2)

def execute_cc3(r, neq):
# The CC3 model: An iterative coupled cluster approach including connected triples
# Koch, Henrik Christiansen, Ove Jorgensen, Poul Sanchez De Meras, Alfredo M.
# The Journal of Chemical Physics 1997 v.106 issue 5 p. 1808

    if neq == 100:
        rint = r.integrate(bra = ['a', 'i']).scale(0.5)
    elif neq == 101:
        rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'])
        rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'])
        rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
    elif neq == 86:
        rint1 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k']).scale(1.0/6)
        rint1 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i']).scale(-1.0/6)
        # rint4 = deepcopy(rint1).scale(-1.0)
        # rint4.multisubst(['i','j','k'],['j','k','i'])
        rint = rint1 + rint4
    
#    print('po RINT')
#    for x in rint:
#        print(x)

    rint.exec_delta()
    rsimp = simplify(rint)
    rsimp.cleanup()

    if neq == 86:
        rsimp_perm = delete_permutations(rsimp)
    rsimp = deepcopy(rsimp_perm)
    rsimp.cleanup()

    for x in rsimp:
        x.optimize()

    rsimp.cleanup()
        
    for x in range(0, len(rsimp_perm)):
        print(rsimp_perm[x])

    function_template_cc3(rsimp_perm, neq)


def execute_cisd(BRA, KET):
    
    r = hamiltonian

    if BRA == 0 and KET == 0:
        rsimp = block_hf_hf(r)
    elif BRA == 2 and KET ==0:
        rsimp = block_d_hf(r)
    elif BRA == 1 and KET ==1:
        rsimp = block_s_s(r)
    elif BRA == 2 and KET ==1:
        rsimp = block_d_s(r)
    elif BRA == 0 and KET ==2:
        rsimp = block_hf_d(r)
    elif BRA == 1 and KET ==2:
        rsimp = block_s_d(r)
    elif BRA == 2 and KET ==2:
        rsimp = block_d_d(r)

    print('rsimp')
    for x in rsimp:
        print(x)
    sys.exit(0)
    delta_list, name1, name_list, arg_list = eom_func_cisd(rsimp, BRA, KET)
    
    print('name')
    for x in name_list:
        print(x)
    print('')
    jacobian_loop('cisd', arg_list, 'cisd', (BRA,KET), delta_list, name1, name_list)

def block_hf_hf(r):

    #BLOCK HF-HF
    rint = r.integrate()
    rint.exec_delta()
    rsimp = simplify(rint)
    
    return rsimp

def block_d_hf(r):
    
    #BLOCK D-HF
    # rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s', 's'])
    # rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], braspin = ['s', 's'])
    # rint = rint1.scale(1./3.) + rint2.scale(1./6.)

    # r = eaibj
    # rint = r.integrate(bra = ['c', 'k', 'd', 'l'], braspin = ['s', 's'])
    
    rint1 = r.integrate(bra = ['c', 'k', 'd', 'l'], braspin = ['s', 's'])
    rint2 = r.integrate(bra = ['c', 'l', 'd', 'k'], braspin = ['s', 's'])
    rint = rint1.scale(1./3.) + rint2.scale(1./6.)

    rint.exec_delta()
    rsimp = simplify(rint)

    rint = deepcopy(rsimp)

    # r2 = deepcopy(rint)
    # for x in r2:
    #     x.new_delta("a", "b")
    #     x.new_delta("i", "j")
    # rint += r2.scale(-1./2.)

    return rint
    
def block_s_s(r):

    #BLOCK S-S
    rint = r.integrate(bra = ['a', 'i'], ket = ['c', 'k'], braspin = ['s'], ketspin = ['s']).scale(0.5)
    rint.exec_delta()
    rsimp = simplify(rint)

    return rsimp


def block_d_s(r):
    
    #BLOCK D-S
    # rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], ket = ['c', 'k'], braspin = ['s', 's'], ketspin = ['s'])
    # rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], ket = ['c', 'k'], braspin = ['s', 's'], ketspin = ['s'])

    # rint = r.integrate(bra = ['c', 'k', 'd', 'l'], ket = ['a', 'i'], braspin = ['s', 's'], ketspin = ['s']).scale(0.5)

    rint1 = r.integrate(bra = ['c', 'k', 'd', 'l'], ket = ['a', 'i'], braspin = ['s', 's'], ketspin = ['s'])
    rint2 = r.integrate(bra = ['c', 'l', 'd', 'k'], ket = ['a', 'i'], braspin = ['s', 's'], ketspin = ['s'])

    rint = rint1.scale(1./3.) + rint2.scale(1./6.)

    rint.exec_delta()
    rsimp = simplify(rint)
    rint = deepcopy(rsimp)

    # r2 = deepcopy(rint)
    # for x in r2:
    #     x.new_delta("a", "b")
    #     x.new_delta("i", "j")
    # rint += r2.scale(-1./2.)
    
    return rint
    

def block_hf_d(r):

    #BLOCK HF-D
    # rint = r.integrate(ket = ['a', 'i', 'b', 'j'], ketspin = ['s', 's'])
    rint = r.integrate(ket = ['c', 'k', 'd', 'l'], ketspin = ['s', 's'])
    rint.exec_delta()
    rsimp = simplify(rint)
    
    return rsimp

def block_s_d(r):
    
    #BLOCK S-D
#    r = arithmetic_string(h)
#    rint = r.integrate(bra = ['a', 'i'], ket = ['b', 'j', 'c', 'k'], braspin = ['s'], ketspin = ['s', 's']).scale(0.5)
    rint = r.integrate(bra = ['a', 'i'], ket = ['c', 'k', 'd', 'l'], braspin = ['s'], ketspin = ['s', 's']).scale(0.5)

    rint.exec_delta()
    rsimp = simplify(rint)
    
    return rsimp

def block_d_d(r):
    
    # BLOCK D-D

    # rint = r.integrate(bra = ['a', 'i', 'b', 'j'], ket = ['c', 'k', 'd', 'l'], braspin = ['s', 's'], ketspin = ['s', 's'])


    rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], ket = ['a', 'i', 'b', 'j'], braspin = ['s', 's'], ketspin = ['s', 's'])
    rint2 = r.integrate(bra = ['a', 'i', 'b', 'j'], ket = ['a', 'j', 'b', 'i'], braspin = ['s', 's'], ketspin = ['s', 's'])
    rint3 = r.integrate(bra = ['a', 'j', 'b', 'i'], ket = ['a', 'i', 'b', 'j'], braspin = ['s', 's'], ketspin = ['s', 's'])
    rint4 = r.integrate(bra = ['a', 'j', 'b', 'i'], ket = ['a', 'j', 'b', 'i'], braspin = ['s', 's'], ketspin = ['s', 's'])
    rint  = rint1.scale(1./18.) + rint2.scale(1./36.)+ rint3.scale(1./9.) + rint4.scale(1./18.)

    # rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], ket = ['c', 'k', 'd', 'l'], braspin = ['s', 's'], ketspin = ['s', 's'])
    # rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], ket = ['c', 'k', 'd', 'l'], braspin = ['s', 's'], ketspin = ['s', 's'])

    
    # rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
    rint.exec_delta()
    delta_subst = delta_to_dict([])
    #rint.exec_delta_fixed(delta_subst)
    rsimp = simplify(rint)

    rint = deepcopy(rsimp)

    # r2 = deepcopy(rint)
    # for x in r2:
    #     x.new_delta("a", "b")
    #     x.new_delta("i", "j")
    # rint += r2.scale(-1./2.)

    return rint


def execute_rowe(BRA, KET, v, op):

    if BRA == 2 and KET ==0:
        rsimp = rowe_block_d_hf(op, v)
    elif BRA == 1 and KET ==1:
        rsimp = rowe_block_s_s(op, v)
    elif BRA == 2 and KET ==1:
        rsimp = rowe_block_d_s(op, v)
    elif BRA == 0 and KET ==2:
        rsimp = rowe_block_hf_d(op, v)
    elif BRA == 1 and KET ==2:
        rsimp = rowe_block_s_d(op, v)
    elif BRA == 2 and KET ==2:
        rsimp = rowe_block_d_d(op, v)

    print('rsimp')
    for x in rsimp:
        print(x)
    sys.exit(0)


def rowe_block_d_hf(r, v):

    # AA = ugg()
    # AA.operator_idx.append(['i', 'j'])
    # AA.operator_idx.append(['m', 'n'])
    # AA.operator_type.append("s")
    # AA.operator_type.append("s")


    # [BB, H, AA]
    if v == 2:
        BB1 = ebjai
        BB2 = ebiaj
    elif v ==1:
        BB1 = ejbia
        BB2 = eibja
        
    #biorthonormal basis

    r = evaluate(BB1,hamiltonian).scale(1./3.) + evaluate(BB2, hamiltonian).scale(1./6.)

    # r1 = evaluate(BB1,hamiltonian, AA).scale(1./3.) + evaluate(BB2, hamiltonian, AA).scale(1./6.)
    # r2 = evaluate(hamiltonian, AA, BB1).scale(1./3.) + evaluate(hamiltonian, AA, BB2).scale(1./6.)
    # r = r1.scale(0.5) + r2.scale(-0.5)


    # # [BB, H, AA]
    # if v == 1:
    #     BB1 = ebjai
    #     BB2 = ebiaj
    # elif v ==2:
    #     BB1 = ejbia
    #     BB2 = eibja
        

    # #biorthonormal basis
    # r = evaluate(BB1, hamiltonian).scale(1./3.) + evaluate(BB2, hamiltonian).scale(1./6.)
    rint = r.integrate()
    rint.exec_delta()
    rsimp = simplify(rint)

    return rsimp

def rowe_block_hf_d(r, v):

#     AA = ugg()
#     AA.operator_idx.append(['i', 'j'])
#     AA.operator_idx.append(['k', 'l'])
#     AA.operator_type.append("s")
#     AA.operator_type.append("s")

#     BB = ugg()
#     BB.operator_idx.append(['i', 'j'])
# #    BB.operator_idx.append(['m', 'n'])
#     BB.operator_type.append("s")
# #    BB.operator_type.append("s")
#     AA = eckdl
# #    AA = ekcld
# #    AA = ekcld

    if v == 1:
        AA = eckdl
    elif v ==2:
        AA = ekcld

#    r = evaluate(BB, hamiltonian, AA).scale(0.5) + evaluate(hamiltonian, AA, BB).scale(-0.5)
    r = evaluate(hamiltonian, AA)
    print('ha')
    for x in r:
        print(x)
        
#    sys.exit(0)

    # # [BB, H, AA]
    # if v == 1:
    #     AA = eckdl
    # elif v ==2:
    #     AA = ekcld

    # r = evaluate(hamiltonian, AA)

    rint = r.integrate()
    rint.exec_delta()
    rsimp = simplify(rint)
    print('wyn')
    for x in rsimp:
        print(x)


    return rsimp

def rowe_block_s_s(op, v):

    # [BB, H, AA]

    if v == 3:
        AA = eckr
        BB = eai
    elif v ==4:
        AA = ekcr
        BB = eai
    elif v==1:
        AA = eckr
        BB = eia
    elif v==2:
        AA = ekcr
        BB = eia

    if op == 1:
        r = evaluate(BB, hamiltonian, AA).scale(0.5) + evaluate(hamiltonian, AA, BB).scale(-0.5)
    else :
        r = evaluate(BB, AA)

    rint = r.integrate().scale(0.5)
    rint.exec_delta()
    rsimp = simplify(rint)

    return rsimp

def rowe_block_s_d(op, v):


    # [BB, H, AA]
    if v == 3:
        BB = eai
        AA = eckdl
    elif v ==4:
        BB = eai
        AA = ekcld
    elif v ==1:
        BB = eia
        AA = eckdl
    elif v ==2:
        BB = eia
        AA = ekcld

    #biorthonormal basis
    if op == 1:
        r = evaluate(BB, hamiltonian, AA).scale(0.25) + evaluate(hamiltonian, AA, BB).scale(-0.25)
    else:
        r = evaluate(BB, AA)

    rint = r.integrate()
    rint.exec_delta()
    rsimp = simplify(rint)

    return rsimp


def rowe_block_d_s(op, v):


    # [BB, H, AA]
    if v == 3:
        AA = eckr
        BB1 = ebjai
        BB2 = ebiaj

    elif v ==4:
        AA = ekcr
        BB1 = ebjai
        BB2 = ebiaj
    elif v ==1:
        AA = eckr
        BB1 = ejbia
        BB2 = eibja

    elif v ==2:
        AA = ekcr
        BB1 = ejbia
        BB2 = eibja

        
    #biorthonormal basis

    if op == 1:
        r1 = evaluate(BB1,hamiltonian, AA).scale(1./3.) + evaluate(BB2, hamiltonian, AA).scale(1./6.)
        r2 = evaluate(hamiltonian, AA, BB1).scale(1./3.) + evaluate(hamiltonian, AA, BB2).scale(1./6.)
        r = r1.scale(0.5) + r2.scale(-0.5)
    else :
        r = evaluate(BB1, AA).scale(1./3.) + evaluate(BB2, AA).scale(1./6.)

    print('wynik')
    for x in r:
        print(x)

    rint = r.integrate()
    rint.exec_delta()
    rsimp = simplify(rint)

    return rsimp
    

def rowe_block_d_d(op, v):

    # BB1 = ebjai
    # BB2 = ebiaj

    # r = evaluate(hamiltonian, eai, ekcld)
    # for x in r:
    #     print(x)
    # sys.exit(0)
    # t1 = evaluate(ebjai, deepcopy(r)).scale(1./3.)
    # # t2 = evaluate(eiajb, deepcopy(r)).scale(1./3.)
    # # t3 = evaluate(ejaib, deepcopy(r)).scale(1./6.)
    # # t4 = evaluate(ejaib, deepcopy(r)).scale(1./6.)
    # t = t1#+t2+t3+t4
    # print('laaa wynik ttt')
    # for x in t:
    #     print(x)
    
    # rint = t.integrate()
    # rint.exec_delta()
    # rsimp = simplify(rint)

    # print('laaa wynik')
    # for x in rsimp:
    #     print(x)


#    sys.exit(0)

              
    # [BB, H, AA]
    if v == 3:
        AA = eckdl
        BB1 = ebjai
        BB2 = ebiaj
    elif v ==4:
        AA = ekcld
        BB1 = ebjai
        BB2 = ebiaj
    elif v ==1:
        AA = eckdl
        BB1 = ejbia
        BB2 = eibja
    elif v ==2:
        AA = ekcld
        BB1 = ejbia
        BB2 = eibja
        

    #biorthonormal basis

    if op == 1:
        r1 = evaluate(BB1, hamiltonian, AA).scale(1./3.) + evaluate(BB2, hamiltonian, AA).scale(1./6.)
        r2 = evaluate(hamiltonian, AA, BB1).scale(1./3.) + evaluate(hamiltonian, AA, BB2).scale(1./6.)
        r = r1.scale(0.5) + r2.scale(-0.5)
    else:
        r = evaluate(BB1, AA).scale(1./3.) + evaluate(BB2, AA).scale(1./6.)
        # r2 = evaluate(AA, BB1).scale(1./3.) + evaluate(AA, BB2).scale(1./6.)
        # r = r1.scale(0.5) + r2.scale(-0.5)
    
    rint = r.integrate()
    rint.exec_delta()
    rsimp = simplify(rint)

    # print('laaa wynik')
    # for x in rsimp:
    #     print(x)



    return rsimp


def execute_cisd_rowe_triplet(BRA, KET, rowe, v, pm):
    
    if BRA == 0 and KET == 0:
        rsimp = block_00_trip(rowe, v, pm)
    elif BRA == 2 and KET ==0:
        rsimp = block_20_trip(pm)
    elif BRA == 1 and KET ==1:
        rsimp = block_11_trip(rowe, v, pm)
    elif BRA == 2 and KET ==1:
        rsimp = block_21_trip(rowe, v, pm)
    elif BRA == 0 and KET ==2:
        rsimp = block_02_trip(pm)
    elif BRA == 1 and KET ==2:
        rsimp = block_12_trip(rowe, v, pm)
    elif BRA == 2 and KET ==2:
        rsimp = block_22_trip(rowe, v, pm)

    print('rsimp')
    for x in rsimp:
        print(x)
    sys.exit(0)


def block_20_trip(pm):
    
    LL1 = etjbia
    LL2 = etiajb
    LL1.num_factor = (1.0/8.0)
    LL2.num_factor = (1.0/8.0)

    if pm == -1:
        LL2.num_factor = LL2.num_factor * (-1.0)


    r1 = deepcopy(hamiltonian).fromleft(LL1)
    r2 = deepcopy(hamiltonian).fromleft(LL2)

    rint = r1.integrate() + r2.integrate()

    rint.exec_delta()
    rsimp = simplify(rint)
    rint = deepcopy(rsimp)

    return rint


def block_02_trip(pm):

    RR1 = etckdl
    RR2 = etdlck

    if pm == -1:
        RR2.num_factor = RR2.num_factor * (-1.0)

        
    r1 = deepcopy(hamiltonian).fromright(RR1)
    r2 = deepcopy(hamiltonian).fromright(RR2)
    rint = r1.integrate() + r2.integrate()

    rint.exec_delta()
    rsimp = simplify(rint)
    rint = deepcopy(rsimp)

    return rint


def block_12_trip(rowe, v, pm):

    if v == 1:
        LL = ttia
        RR1 = etckdl
        RR2 = etdlck
    elif v ==2:
        if rowe == False:
            print('inadequate version for cisd')
            sys.exit()
        LL = ttia
        RR1 = etkcld
        RR2 = etldkc
    elif v == 3:
        if rowe == False:
            print('inadequate version for cisd')
            sys.exit()
        LL = ttai
        RR1 = etckdl
        RR2 = etdlck
    elif v ==4:
        if rowe == False:
            print('inadequate version for cisd')
            sys.exit()
        LL = ttai
        RR1 = etkcld
        RR2 = etldkc

    if pm == -1:
        RR2.num_factor = RR2.num_factor * (-1.0)

    LL.num_factor *= 0.5
    
    if rowe == False:
        r1 = arithmetic_string()
        for x in hamiltonian:
            disambiguate(x, LL, RR1)
            rt = deepcopy(x).fromleft(LL)
            rt = rt.fromright(RR1)
            r1.append(rt)
        r2 = arithmetic_string()
        for x in hamiltonian:
            disambiguate(x, LL, RR2)
            rt = deepcopy(x).fromleft(LL)
            rt = rt.fromright(RR2)
            r2.append(rt)

        # r1 = deepcopy(hamiltonian).fromleft(LL)
        # r1 = r1.fromright(RR1)
        # r2 = deepcopy(hamiltonian).fromleft(LL)
        # r2 = r2.fromright(RR2)
    else:
        r1 = evaluate(LL, hamiltonian, RR1).scale(0.5) + evaluate(hamiltonian, RR1, LL).scale(-0.5)
        r2 = evaluate(LL, hamiltonian, RR2).scale(0.5) + evaluate(hamiltonian, RR2, LL).scale(-0.5)
    
    rint = r1.integrate() + r2.integrate()


    rint.exec_delta()
    rsimp = simplify(rint)
    rint = deepcopy(rsimp)

    return rint


def block_11_trip(rowe, v,  pm):

    if v == 1:
        LL = ttia
        RR = ttck
    elif v ==2:
        if rowe == False:
            print('inadequate version for cisd')
            sys.exit()
        LL = ttia
        RR = ttkc
    elif v == 3:
        if rowe == False:
            print('inadequate version for cisd')
            sys.exit()
        LL = ttai
        RR = ttck
    elif v ==4:
        if rowe == False:
            print('inadequate version for cisd')
            sys.exit()
        LL = ttai
        RR = ttkc

    RR.num_factor *= 0.5
        
    
    if rowe == False:
        r = arithmetic_string()
        for x in hamiltonian:
            disambiguate(x, LL, RR)
            r2 = deepcopy(x).fromleft(LL)
            r2 = r2.fromright(RR)
            r.append(r2)
    else:
        r = evaluate(LL, hamiltonian, RR).scale(0.5) + evaluate(hamiltonian, RR, LL).scale(-0.5)
    
    print('rr')
    for x in r:
        print(x)

    rint = r.integrate()

    print('rr')
    for x in rint:
        print(x)

    rint.exec_delta()
    rsimp = simplify(rint)
    rint = deepcopy(rsimp)

    return rint

def block_21_trip(rowe, v, pm):

    if v == 1:
        RR = ttai
        LL1 = letldkc
        LL2 = letkcld
        LL1.num_factor = (1.0/8.0)
        LL2.num_factor = (1.0/8.0)

    elif v == 2:
        if rowe == False:
            print('inadequate version for cisd')
            sys.exit()

        RR = ttia
        LL1 = letldkc
        LL2 = letkcld
        LL1.num_factor = (1.0/8.0)
        LL2.num_factor = (1.0/8.0)

    elif v ==3:
        if rowe == False:
            print('inadequate version for cisd')
            sys.exit()
        RR = ttai
        LL1 = letdlck
        LL2 = letckdl
        LL1.num_factor = (1.0/8.0)
        LL2.num_factor = (1.0/8.0)
    elif v ==4:
        if rowe == False:
            print('inadequate version for cisd')
            sys.exit()
        RR = ttia
        LL1 = letdlck
        LL2 = letckdl
        LL1.num_factor = (1.0/8.0)
        LL2.num_factor = (1.0/8.0)


    if pm == -1:
        print('SRA')
        print(LL2)
        LL2.num_factor = LL2.num_factor * (-1.0)
        print(LL2)
        print('')


    if rowe == False:
        r1 = arithmetic_string()
        for x in hamiltonian:
            print(LL1, x, RR)
            disambiguate(x, LL1, RR)
            print(LL1, x, RR)
            print('-----')
            rt = deepcopy(x).fromleft(LL1)
            rt = rt.fromright(RR)
            r1.append(rt)
        r2 = arithmetic_string()
        for x in hamiltonian:
            print(LL2, x, RR)
            disambiguate(x, LL2, RR)
            print(LL2, x, RR)
            print('-------')
            rt = deepcopy(x).fromleft(LL2)
            rt = rt.fromright(RR)
            r2.append(rt)

        # r1 = deepcopy(hamiltonian).fromleft(LL1)
        # r1 = r1.fromright(RR)
        # r2 = deepcopy(hamiltonian).fromleft(LL2)
        # r2 = r1.fromright(RR)
    else:
        r1 = evaluate(LL1, hamiltonian, RR).scale(0.5) + evaluate(hamiltonian, RR, LL1).scale(-0.5)
        r2 = evaluate(LL2, hamiltonian, RR).scale(0.5) + evaluate(hamiltonian, RR, LL2).scale(-0.5)



    rint = r1.integrate() + r2.integrate()

    # for xz in rint:
    #     print(xz)
    # sys.exit(0)

    rint.exec_delta()
    rsimp = simplify(rint)
    rint = deepcopy(rsimp)

    return rint


def block_22_trip(rowe, v, pm):
    
    if v == 1:
        RR1 = etckdl
        RR2 = etdlck
        LL1 = etjbia
        LL2 = etiajb
        LL1.num_factor = (1.0/8.0)
        LL2.num_factor = (1.0/8.0)

    elif v == 2:
        if rowe == False:
            print('inadequate version for cisd')
            sys.exit()
        RR1 = etkcld
        RR2 = etldkc
        LL1 = etjbia
        LL2 = etiajb
        LL1.num_factor = (1.0/8.0)
        LL2.num_factor = (1.0/8.0)

    elif v ==3:
        if rowe == False:
            print('inadequate version for cisd')
            sys.exit()
        RR1 = etckdl
        RR2 = etdlck
        LL1 = etbjai
        LL2 = etaibj
        LL1.num_factor = (1.0/8.0)
        LL2.num_factor = (1.0/8.0)
    elif v ==4:
        if rowe == False:
            print('inadequate version for cisd')
            sys.exit()
        RR1 = etkcld
        RR2 = etldkc
        LL1 = etbjai
        LL2 = etaibj
        LL1.num_factor = (1.0/8.0)
        LL2.num_factor = (1.0/8.0)


    if pm == -1.0:
        LL2.num_factor = LL2.num_factor * (-1.0)
        RR2.num_factor = RR2.num_factor * (-1.0)
    elif pm == -0.5:
            LL2.num_factor = LL2.num_factor * (-1.0)
    elif pm == 0.5:
            RR2.num_factor = RR2.num_factor * (-1.0)

    if rowe == False:
        r1 = arithmetic_string()
        for x in hamiltonian:
            disambiguate(x, LL1, RR1)
            rt = deepcopy(x).fromleft(LL1)
            rt = rt.fromright(RR1)
            r1.append(rt)
        r2 = arithmetic_string()
        for x in hamiltonian:
            disambiguate(x, LL2, RR1)
            rt = deepcopy(x).fromleft(LL2)
            rt = rt.fromright(RR1)
            r2.append(rt)
        r3 = arithmetic_string()
        for x in hamiltonian:
            disambiguate(x, LL1, RR2)
            rt = deepcopy(x).fromleft(LL1)
            rt = rt.fromright(RR2)
            r3.append(rt)
        r4 = arithmetic_string()
        for x in hamiltonian:
            disambiguate(x, LL2, RR2)
            rt = deepcopy(x).fromleft(LL2)
            rt = rt.fromright(RR2)
            r4.append(rt)

        # r1 = deepcopy(hamiltonian).fromleft(LL1)
        # r1 = r1.fromright(RR1)
        # r2 = deepcopy(hamiltonian).fromleft(LL2)
        # r2 = r1.fromright(RR2)
        # r3 = deepcopy(hamiltonian).fromleft(LL1)
        # r3 = r1.fromright(RR2)
        # r4 = deepcopy(hamiltonian).fromleft(LL2)
        # r4 = r1.fromright(RR1)



    else:
        r1 = evaluate(LL1, hamiltonian, RR1).scale(0.5) + evaluate(hamiltonian, RR1, LL1).scale(-0.5)
        r2 = evaluate(LL1, hamiltonian, RR2).scale(0.5) + evaluate(hamiltonian, RR2, LL1).scale(-0.5)
        r3 = evaluate(LL2, hamiltonian, RR1).scale(0.5) + evaluate(hamiltonian, RR1, LL2).scale(-0.5)
        r4 = evaluate(LL2, hamiltonian, RR2).scale(0.5) + evaluate(hamiltonian, RR2, LL2).scale(-0.5)

        # r1 = evaluate(LL1, RR1).scale(0.5) + evaluate( RR1, LL1).scale(-0.5)
        # r2 = evaluate(LL1, RR2).scale(0.5) + evaluate( RR2, LL1).scale(-0.5)
        # r3 = evaluate(LL2, RR1).scale(0.5) + evaluate( RR1, LL2).scale(-0.5)
        # r4 = evaluate(LL2, RR2).scale(0.5) + evaluate( RR2, LL2).scale(-0.5)

    rint = r1.integrate() + r2.integrate() + r3.integrate() + r4.integrate()

    rint.exec_delta()
    rsimp = simplify(rint)
    rint = deepcopy(rsimp)

    return rint


    
    
    
