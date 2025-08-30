from params import *
import paldus_classes
# from paldus_classes import ugg
# from paldus_classes import arithmetic_string
from copy import deepcopy
# from paldus_classes import disambiguate
# from paldus_classes import pair_permutations
import math
# from paldus_classes import integrate
# from paldus_classes import virtual, occupied, general
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


def test_nowy():

#     z = ugg()
#     z.operator_idx.append(["c", "k"])
# #    z.operator_idx.append(["d", "l"])
#     z.operator_type.append("s")
#  #   z.operator_type.append("s")

#     z2 = ugg()
#     z2.operator_idx.append(["i", "a"])
#   #  z2.operator_idx.append(["j", "b"])
#     z2.operator_type.append("s")
#     # z2.operator_type.append("s")


#     z = ugg()
#     z.operator_idx.append(["l", "k"])
# #    z.operator_idx.append(["d", "l"])
#     z.operator_type.append("s")
#  #   z.operator_type.append("s")

#     z2 = ugg()
#     z2.operator_idx.append(["i", "j"])
#   #  z2.operator_idx.append(["j", "b"])
#     z2.operator_type.append("s")
#    # z2.operator_type.append("s")


#     r = z2.fromright(z)
#     print(r)
#     rint = r.integrate()
#     print('')
#     for x in rint:
#         print(x)

#     sys.exit(0)

    z = ugg()
    z.summation = ["p","q", 'k', 'l']
    z.coefficient = ["z1"]
    z.coefficient.append("z2")
    z.coefficient.append("z3")
    z.coefficient.append("z4")
    z.coefficient_idx.append(["p"])
    z.coefficient_idx.append(["q"])
    z.coefficient_idx.append(["k"])
    z.coefficient_idx.append(["l"])
    z.operator_idx.append(["p", "k"])
    z.operator_idx.append(["q", "l"])
    z.operator_type.append("s")
    z.operator_type.append("s")


    z1 = ugg()
    z1.summation = ["r","s", 'j', 'i']
    z1.coefficient = ["x1"]
    z1.coefficient.append("x2")
    z1.coefficient.append("x3")
    z1.coefficient.append("x4")
    z1.coefficient_idx.append(["r"])
    z1.coefficient_idx.append(["s"])
    z1.coefficient_idx.append(["i"])
    z1.coefficient_idx.append(["j"])
    # z1.summation = ["r","s"]
    # z1.coefficient = ["z2"]
    # z1.coefficient_idx.append(["r", "s"])
    z1.operator_idx.append(["i", "r"])
    z1.operator_idx.append(["j", "s"])
    z1.operator_type.append("s")
    z1.operator_type.append("s")
    print(z1)
#    sys.exit(0)

    KL = evaluate(z1, hamiltonian, z).scale(0.5)+  evaluate(hamiltonian, z, z1).scale(-0.5)

    for k in KL:
        print(k)
    print('')
#    sys.exit(0)
    rint = KL.integrate()
    rsimp = simplify(rint)
    print('wyn')
    for x in rsimp:
        print(x)
    sys.exit(0)






def test_macierzy_S_onedet():

    eia = ugg()
    eia.summetion = ["i","a"]
    eia.coefficient = ["z"]
    eia.coefficient_idx.append(["i","a"])
    eia.operator_idx.append(["i","a"])
    eia.operator_type.append("s")


    ejb = ugg()
    ejb.summetion = ["j","b"]
    ejb.coefficient = ["z"]
    ejb.coefficient_idx.append(["j","b"])
    ejb.operator_idx.append(["j","b"])
    ejb.operator_type.append("s")

    eck = ugg()
    eck.summetion = ["c","k"]
    eck.coefficient = ["z"]
    eck.coefficient_idx.append(["c","k"])
    eck.operator_idx.append(["c","k"])
    eck.operator_type.append("s")

    edl = ugg()
    edl.summetion = ["d","l"]
    edl.coefficient = ["z"]
    edl.coefficient_idx.append(["d","l"])
    edl.operator_idx.append(["d","l"])
    edl.operator_type.append("s")


    bra = eia.fromright(ejb)
    ket = eck.fromright(edl)
    
    r = evaluate(bra, ket)
    rint = r.integrate()
    rsimp = simplify(rint)

    print('')
    print('wynik')
    for x in rsimp:
        print(x)

    sys.exit(0)
def test_komutatora():

    eia = ugg()
    eia.operator_idx.append(["i","a"])
    eia.operator_type.append("s")

    ejb = ugg()
    ejb.operator_idx.append(["j","b"])
    ejb.operator_type.append("s")

    eck = ugg()
    eck.operator_idx.append(["c","k"])
    eck.operator_type.append("s")

    ff = hamiltonian
    print('przed')
    for x in ff:
        print(x)

    # sys.exit(0)
    
    # w = ff.integrate()
    # rsimp = simplify(w)

    # print('wynik')
    # for x in rsimp:
    #     print(x)

    # sys.exit(0)
    
    print(ff)
    print(eia)
    print('')

    for x in ff:
        disambiguate(eia, ejb, ekc, x)
    print('')
    print('po')
    for x in ff:
        print(x)

    print(eia)

    e = eia.fromright(ejb)


    # w = evaluate(e, hamiltonian, eck)
    w = evaluate(e, eck)
    print('')
    print('sym')
    for x in w:
        print(x)

    
    rint = w.integrate()
    rsimp = simplify(rint)


    print('wynik')
    for x in rsimp:
        print(x)

    sys.exit(0)


def test_bazy_biort():

    #  BAZA DLA SINGLI
    # r = ugg()
    # r.operator_idx.append(["a","i"])
    # r.operator_type.append("s")


    # rint = r.integrate(bra = ['c', 'k'], braspin = ['s']).scale(0.5)
    # for x in rint:
    #     print(x)

    # print('wwynika')
    # for x in rint:
    #     print(x)
    # sys.exit(0)


    # BAZA DLA DUBLI

    # r = ugg()
    # r.operator_idx.append(["c","k"])
    # r.operator_idx.append(["d","l"])
    # r.operator_type.append("s")
    # r.operator_type.append("s")



    # rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s', 's'])
    # rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], braspin = ['s', 's'])
    # rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
    
    # rsimp = simplify(rint)
    # rint = deepcopy(rsimp)
    
    # r2 = deepcopy(rint)
    # for x in r2:
    #     x.new_delta("a", "b")
    #     x.new_delta("i", "j")
    # rint += r2.scale(-1./2.)


    # print('wwynika')
    # for x in rint:
    #     print(x)
    # sys.exit(0)


    r = ugg()
    r.operator_idx.append(["d","l"])
    r.operator_idx.append(["e","m"])
    r.operator_idx.append(["f","n"])
    r.operator_type.append("s")
    r.operator_type.append("s")
    r.operator_type.append("s")

    
    # print('r1 integration')
    # rint1 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./4.)
    # print('wwynika')
    # for x in rint1:
    #     print(x)
    # rsimp  = simplify(rint1)

    # print('wwynikss')
    # for x in rsimp:
    #     print(x)
    # sys.exit(0)

    # sys.exit(0)

    rint1 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./4.)
    rint2 = r.integrate(bra = ['a', 'k', 'b', 'i', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./12.)
    rint3 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./6.)
    rint4 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./6)
    rint5 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./12.)
    r1 = rint1 + rint2 + rint3 + rint4 + rint5
    # print('wwynik')
    # for x in r1:
    #     print(x)
    # sys.exit(0)
    
    rsimp  = simplify(r1)

    print('wwynik')
    for x in rsimp:
        print(x)
    sys.exit(0)





