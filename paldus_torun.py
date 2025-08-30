from params import *
import paldus_classes
from paldus_basic import *
from paldus_classes import ugg
from copy import deepcopy
from collections import deque
import math
from itertools import product
from itertools import permutations
from templates import *
import sys
from paldus_f12 import *
import time
import pickle
import io
import pickle
import re

coef_pool = set()


def Tlatex_read():

    res = []

#    f = open("dipole_moment-hlf.out", "r")
#    f = open("dipole_moment.out", "r")
#    f = open("dpp.out", "r")
    f = open("dipole_moment13.out", "r")
    #f = open("dmh.out", "r")
#    f = open("exch_disp20_s4.out", "r")
    lines = f.readlines()

    res = arithmetic_string()
    for line in lines:

        x = read_ugg_from_str(line)
        res.append(x)


    print('resres')
    k=0
    for x in res:
        print(k, x)
        k+=1
        print('')
    rsimp = simplify(res)
    for x in rsimp:
        print(x)
    print('dlugosc po simplify', len(res), len(rsimp))

    
    start_time = time.time()
    excl = []
    # cProfile.runctx('factorize_driver(rsimp, excl, COST=6, \
    # filel=[], s_ampli=False, disconnected = True)', {'rsimp':rsimp, 'excl':excl, 'COST':6, 'filel':[],\
    #                                                  's_ampli':False, 'disconnected':True, 'factorize_driver':factorize_driver}, {})
    basket_outer, basket_nointerm, list_of_int, n_max, all_hash, mem_dict, \
        interm_fx, xfx_dict, list_of_int_xfx  = factorize_driver(rsimp, [], COST=6, filel=[], s_ampli=False, disconnected = True)

    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"DISC Elapsed time for all: {elapsed_time} seconds")
    sys.exit(0)
    print('baskett', len(basket_outer))
    for x in basket_outer:
        print(x)

    print('interm length')
    k = 0
    for x in list_of_int:
        print('oto poziom', k, len(x))
        k+=1
        for y in x:
            print(y)
        print('')
        


def read_ugg_from_str(line):

    # read coefficients
    symbols = re.findall(r'\b(\w+)\b(?![^{]*})',line)
    coef = [s.replace('_', '') for s in symbols]
    new_lst = []
    for item in coef:
        new_lst.extend(item.split())
    coefs = [s for s in new_lst if s.isalpha()]

    # add new coefficients to coefficients_set
    coef_pool_mini = coef_pool| set(coefs)

    # read coef idx
    coef_idx = re.findall(r'{([^\d{}][\w]*)}', line)

    # read numerical factor
    z = re.split('\s', line)

    no = float(z[2])
    if z[1] == '-':
        no = -no

    for i in range(0, len(coefs)):
        coefs[i] = 'TT'+coefs[i]

    offset = 0
    sumset = set()
    u_list = []
    l_list = []
    for i in range(0, len(coefs)):
        k = i + offset

        # upper indices for this coeffs are in coef_idx[k]
        # lower indices for this coeffs are in coef_idx[k+1]
        u_list_mini = str_to_list(coef_idx[k])
        l_list_mini = str_to_list(coef_idx[k+1])
        sumset  = sumset | set(u_list_mini)
        sumset  = sumset | set(l_list_mini)

        u_list.append(u_list_mini)
        l_list.append(l_list_mini)
        offset += 1

    this = ugg()
    this.summation = list(sumset)
    # this is parameter for printing ugg for this type of code
    this.coefficient = coefs

    this.coefficient_idx = []

    # in this part I change the names from the original files to working indices
    # in paldus code. so a,b, c, i, k ... remain unchanged
    # but j1 -> ja, i1->ia, i2->ib etc
    for i in range(0, len(u_list)):
        this.coefficient_idx.append(u_list[i] + l_list[i])
    this.num_factor = no

    this.standarize()
#    print('this po', this)

    return this


def str_to_list(s):

    # this procedure takes the string of indices and turns it into a list
    # e.g. 'j_3abi_41' --> ['j3', 'a', 'b', 'i41']
    
    s = s.replace("_", "")

    alp_idx = []
    dig_idx = []

    dig_lst = re.findall(r'\d+', s)

    for  i in range(0, len(s)):
        if s[i].isdigit():
            dig_idx.append(i)
        elif s[i].isalpha():
            alp_idx.append(i)

    lst = []
    l = 0
    for k in alp_idx:
        mini_s = s[k]
        if k+1 in dig_idx:
            mini_s += dig_lst[l]
            l += 1
        lst.append(mini_s)

    return lst
        
