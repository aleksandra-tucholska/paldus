from params import *
import paldus_classes
from paldus_basic import *
from paldus_classes import ugg
from paldus_commutators import *
from copy import deepcopy
from copy import copy
from eomccjac import jacobian_loop
from paldus_classes import arithmetic_string
from paldus_classes import disambiguate
from collections import deque
from paldus_classes import pair_permutations
import math
from itertools import product
from itertools import permutations
from itertools import combinations
from paldus_classes import integrate
from paldus_classes import virtual, occupied, general
from fortran_code import *
from templates import  *
from joblib import Parallel, delayed
import sys
import time
import pickle
import numpy as np
from joblib import Parallel, delayed

sys.setrecursionlimit(900000)
COST_THRESH_DEFAULT = 6

MEM_THRESH = 3
ORIGINAL_COST = 0


dct_idx_cost = {'a': 30, 'b':30, 'c':30, 'd':30, 'e':30, 'f':30, 'g':30, 'i':2, 'j':2, 'k':2, 'l':2, 'm':2, 'n':2, 'i>':2, 'j>':2, 'k>':2, 'l>':2,\
                    'Q':100, "S":1000, "P":1000, "R":1000, 'X':10, 'Y':10, 'Z':10, 'V':10, 'W':10, 'U':10,'T':10, 'X>':10, 'Y':10, 'Z>':10, "α":100, "β":100}


mem_real_cost = {'a': 232, 'b':232, 'c':232, 'd':232, 'e':232, 'f':232, 'g':232, 'i':20, 'j':20, 'k':20, 'l':20, 'm':20, 'n':20, 'i>':20, 'j>':20, 'k>':20, 'l>':20,\
                    'Q':100, "S":1000, "P":1000, "R":1000, 'X':10, 'Y':10, 'Z':10, 'V':10, 'W':10, 'U':10, 'T':10, 'X>':10, 'Y':10, 'Z>':10, "α":100, "β":100, \
                     'A':888, 'B':888, 'C':888, 'D':888, 'E':888}



dct_idx_cost = {}
mem_real_cost = {}
for x in virtual:
    dct_idx_cost[x]=30
    mem_real_cost[x]=232
for x in occupied:
    dct_idx_cost[x] = 2
    mem_real_cost[x]=20
for x in CABS:
    dct_idx_cost[x] = 100
    mem_real_cost[x]=888
for x in completev:
    dct_idx_cost[x] = 1000
    mem_real_cost[x]=1000

for i in range(0, len(mona_occup)):
    dct_idx_cost[mona_occup[i]] = 2
    dct_idx_cost[mona_virt[i]] = 30

    dct_idx_cost[monb_occup[i]] = 2
    dct_idx_cost[monb_virt[i]] = 30

    mem_real_cost[mona_occup[i]] = 20
    mem_real_cost[mona_virt[i]] = 232

    mem_real_cost[monb_occup[i]] = 20
    mem_real_cost[monb_virt[i]] = 232




# mem_real_cost = {'a': 200, 'b':200, 'c':200, 'd':200, 'e':200, 'f':200, 'g':200, 'i':10, 'j':10, 'k':10, 'l':10, 'm':10, 'n':10, \
#                      'i>':2, 'j>':2, 'k>':2, 'l>':2,"α":100, "β":100}


def idx_to_cost(e_idx):

    e_cost = []
    for i in e_idx:
        e_cost.append(dct_idx_cost[i])
        
    return e_cost

def find_all_diff_idx(e):

    e_idx = []
#    print('e', e)
    for x in e.coefficient_idx:
        for y in x:
            if y not in e_idx:
                e_idx.append(y)

    e_idx.sort(key=str.casefold)

    return e_idx

def construct_excluded_dict(e, e_idx_in, excl = ''):


    
    for i in range(0, len(e.coefficient)):
        vt = 0
        oc = 0
        for x in e.coefficient_idx[i]:
            if x in excl:
                if x in virtual:
                    vt += 1
                elif x in occupied:
                    oc += 1
        excl_dict.append([vt, oc])

    return excl_dict


def construct_idx_matrix(e, e_idx_in, excl = []):

    # print('z construct_idx_matrix')
    # print('ee', e)
    # print('e_idx_in', e_idx_in)
    # print('excl', excl)
    
    e_idx = []
    if excl != []:
        for i in e_idx_in:
            if i not in excl:
                e_idx.append(i)

    else:
        e_idx = e_idx_in

    # print('excluded', e_idx_in, e_idx, excl)

    len_idx = len(e_idx)
    len_coef = len(e.coefficient)

    # print('len_idx', len_idx)
    # print('len_coef', len_coef)
        
    idx_matrix = np.zeros((len_coef, len_idx+len(excl)), dtype=int)

    for i in range(0, len_coef):
        for x in e.coefficient_idx[i]:
            # print(e, ' ', i, x)
            if x in e_idx:
                place_idx = e_idx.index(x)
                # print('miejsce', i, place_idx)
                idx_matrix[i][place_idx] = 1

    for i in range(0, len_coef):
        for x in e.coefficient_idx[i]:
#            print(e, ' ', i, x)
            if x in excl:
                place_idx = excl.index(x) + len(e_idx)
                # print('miejsce fixed', i, place_idx)
                idx_matrix[i][place_idx] = 1

    # for i in range(0, len_idx):
    #     k = 0
    #     for j in range(0, len_coef):
    #         if idx_matrix[j][i] == 1:
    #             k+=1

    # for x in idx_matrix:
    #     print(x)

    e_idx_out = e_idx + excl
    # print('e_idx_out', e_idx_out)

    return idx_matrix, len_idx, len_coef, e_idx_out

def factorize(nn, e, excl_idx = [], COST = 0, is_x_excluded = False, disconnected=False):


    
    e_idx_in  = find_all_diff_idx(e)
    e_cost = idx_to_cost(e_idx_in)

    print('')
    print('robie wyraz', e)
    print('')

    # print('e_idx na poczatku', e_idx_in)
    # print('excl na poczatku', excl_idx)
    # print('e_cost to:', e_cost)
    # print('')
    omem, ocost = memcost(e)
    ORIGINAL_COST = sum(ocost)
  #  print('omem, ocost', omem, ocost)
    # print('ORIGINAL_COST', ORIGINAL_COST)

    if excl_idx != []:
#        print('sa indeksy excluded, tworze macierz idx_matrix')
        idx_matrix, len_idx, len_coef, e_idx = construct_idx_matrix(e, e_idx_in, excl_idx)
        # excl_list = construct_excluded_dict(e, e_idx, excl)
    else:
 #       print('nie ma indeksow excluded, tworze macierz idx_matrix')
        idx_matrix, len_idx, len_coef, e_idx = construct_idx_matrix(e, e_idx_in)
        # excl_list = construct_excluded_dict(e, e_idx)
        
    k = 0

    haveinterm = True
    if disconnected:
        max_coef = 1
    else:
        max_coef = 2
    if len(e.coefficient) <= max_coef:
#        print('Tutaj jest za malo wyrazow na intermediate ',e)
        haveinterm = False

        biglevel_best = []
        bigparameters_best = []
        return biglevel_best, bigparameters_best, haveinterm

    # if len(e.summation) >= 8:
    #     print('dużo indeksow sumacyjnych, >=8')
    idx_dict = {}
    pair_restr = [[],[]]


    if is_x_excluded:
        excl_coef_level = []
        for i in range(0, len(e.coefficient)):
            if (e.coefficient[i] == OBSERVABLE_X) or (e.coefficient[i] == OBSERVABLE_X_ASYM):
                excl_coef_level.append(i+1)
    else:
        excl_coef_level = []    

    # print('')
    # print('na poczatku dodaje pair_restr - czyli wyrazy z ktorych moge robic pary')
    for x in range(0, len_coef):
        idx_dict[x+1] = x+1
        if excl_coef_level != []:
            if (x+1) not in excl_coef_level:
                pair_restr[0].append(x)
                pair_restr[1].append(x)
        else:
            pair_restr[0].append(x)
            pair_restr[1].append(x)

    # print('pair_restr na poczatku', pair_restr)
    # print('idx_dict', idx_dict)

    minilevel = []
#    excl_coef_lvel = []
    biglevel = []
    minicost = []
    miniparameters = []
    bigcost = []
    bigparameters = []
    pair_status = []
    big_idx_matrix = [idx_matrix]#[deepcopy(idx_matrix)]
    big_idx_dict = [idx_dict]#[deepcopy(idx_dict)]

    pair_status_idx = 0
    idx_matrix_old = idx_matrix#deepcopy(idx_matrix)
    idx_dict_old = idx_dict#deepcopy(idx_dict)
    nn[1] += 1


    # excluded_idx = []
    # for i in range(0, len(e.coefficient)):
    #     if (e.coefficient[i] == "Pp"):
    #         excluded_idx.append(i+1)
    #     if (e.coefficient[i] == "Pm"):
    #         excluded_idx.append(i+1)

    # print('')
    # print('Do excl_coef, zostaly dodane wspolczynniki', excl_coef_level)

            

    if COST != 0:
        COST_THRESH = COST
    else:
        COST_THRESH = COST_THRESH_DEFAULT


    # print('COST_THRESHHOLD IS', COST_THRESH)


    
    # print('')
    # print('zaczynam factorize')
    # print('macierze idx_matrix na samym poczatku')
    # for ii in range(0, len_coef):
    #     print(idx_matrix[ii][:])

    for r1 in range(0, len_idx):
        pp = 0
        for r2 in range(0, len_coef):
            pp+=idx_matrix[r2][r1]
        if pp>2:
            print('tak, ten indeks wystepuje w trzech wyrazach')
            sys.exit(0)
            
    t1 = time.time()
    factorize3(nn, k, pair_status, pair_status_idx, pair_restr, biglevel, minilevel, bigparameters, miniparameters, \
                   bigcost, minicost, idx_dict, idx_dict_old, idx_matrix, idx_matrix_old, len_idx, len_coef, e_idx,\
                   big_idx_matrix, big_idx_dict, excl_coef_level, excl_idx, ORIGINAL_COST, COST_THRESH, disconnected=disconnected)
    t2 = time.time()
    print('factofacto3', t2-t1)

    # for x in biglevel:
    #     print('biglevelpo', x)
    # print('')


    haveinterm = True
    if len(biglevel) == 0:
#        print('dla COST_THRESH = ', COST_THRESH, 'nie ma intermediate bo', ORIGINAL_COST)
        haveinterm = False

        biglevel_best = []
        bigparameters_best = []

    else:
        
        group_idx = []
        first_elem = biglevel[0][0]
        minil = [0]
        for i in range(1, len(biglevel)):
            # print(i, 'group', biglevel[i])
            # print(biglevel[i][0], first_elem)
            if biglevel[i][0]==first_elem:
                # print('tak jest w tej grupie', minil)
                minil.append(i)
            else:
                group_idx.append(deepcopy(minil))
                nn[2]+=1
                minil = []
                first_elem = biglevel[i][0]
                minil.append(i)
                # print('nie, jest w nowej', minil)
        group_idx.append(deepcopy(minil))
        nn[3]+=1
        # print('group_idx', group_idx)
        if len(biglevel) == 1:
            group_idx = [[0]]


        if  haveinterm == True:
            cost_list = []
            num = []

            lowest_cost = [0] * len(group_idx)
            for l in range(0, len(group_idx)):
                lst = group_idx[l]
                lowest_cost[l] = 100
                for i in lst:

            # for i in range(0, len(bigparameters)):
                    this_cost = 0
                    this_mem = 0   
                    for j in range(0, len(biglevel[i])-1):
                        for k in range(0, len(bigparameters[i][j][0])):
                            this_cost += bigparameters[i][j][0][k]
                    for k in range(0, len(bigparameters[i][-1][0])):
                        this_cost += bigparameters[i][-1][0][k]

                    this_cost += bigparameters[i][-1][4][0] # cost of last row
                    if this_cost < lowest_cost[l]:
                        lowest_cost[l] = deepcopy(this_cost)
                        nn[4]+=1
                    # cost_list.append(this_cost)
                    # num.append(i)
            # k = 1
            # cost_list = []
            # num = []
            # for i in range(0, len(bigparameters)):
            #     # print('big gib big', bigparameters[i])
            #     cost_list.append(bigparameters[i][-1][5][1])
            #     # print(bigparameters[i][-1][5])
            #     num.append(i)

            # cost_sort, num_sort = (list(x) for x in zip(*sorted(zip(cost_list, num))))
            print(lowest_cost)


            # lowest_cost = cost_sort[0]
            biglevel_best = []
            bigparameters_best = []
            # s = """N^{a}""".format(a = bigparameters[i][-1][5][0]) 




            for l in range(0, len(group_idx)):
                lst = group_idx[l]
                print('robie grupe', lst)
                for i in lst:



            # for i in range(0, len(biglevel)):

                    this_cost = 0
                    print('')
                    print('--------------------------')
                    for j in range(0, len(biglevel[i])-1):
                        print('max', lowest_cost)
                        print('jj', i, j, bigparameters[i])
                        print(bigparameters[i][j][0][0])
                        #                print(bigparameters[i][-1][0])
                        for k in range(0, len(bigparameters[i][j][0])):
                            this_cost += bigparameters[i][j][0][k]
                        print('this_cost', this_cost)
                    for k in range(0, len(bigparameters[i][-1][0])):
                        this_cost += bigparameters[i][-1][0][k]
                    # this_cost += bigparameters[i][-1][0][0] # cost of last row
                    # print('this_cost-last1', this_cost)
                    this_cost += bigparameters[i][-1][4][0] # cost of last row
                    print('this_cost-last2', this_cost)

                    # if bigparameters[i][-1][5][1]== lowest_cost:
                    print('czy dodaje?', l, this_cost, lowest_cost)
                    if this_cost == lowest_cost[l]:
                        print('dodaje', len(biglevel[i]), e)
                        print(biglevel[i])
                        print(bigparameters[i],this_cost, this_mem)
                        # print(bigparameters[i], len(bigparameters[i]))
                        # print('')
                        # for bigp in range(0, len(bigparameters[i])):
                        #     print('bigp', bigp, bigparameters[i][bigp])
                        biglevel_best.append(biglevel[i])
                        bigparameters_best.append(bigparameters[i])
                        break
                    elif bigparameters[i][-1][5][1] < lowest_cost[l]:
                        qq=0
                        # print('WTF? lowest cost')
                        # print(biglevel[i])
                        # print('')
                    else:
                        print('nie dodaje', e)
                        print(biglevel[i])
                        print(bigparameters[i],this_cost, this_mem)
                        # print(bigparameters[i])
                        # print('')

                        continue
        print(len(biglevel_best), '$$$$$$$', len(biglevel), len(biglevel_best))
        print('===========================================================================================')

    return biglevel_best, bigparameters_best, haveinterm

            

def factorize3(nn, p, pair_status, pair_status_idx, pair_restr, biglevel, minilevel, bigparameters, miniparameters, bigcost, minicost, \
                   idx_dict, idx_dict_old, idx_matrix, idx_matrix_old, len_idx, len_coef, e_idx, \
                   big_idx_matrix, big_idx_dict, excl_coef_level, excl_idx, ORIGINAL_COST, COST_THRESH, disconnected=False):

#    print('', nn, nn[9])
    # print('znow zaczynam factorize3')
    # print('')
    # print('macierze idx_matrix na wejsciu')
    # for ii in range(0, len_coef):
    #     print(idx_matrix[ii][:])
        

    
    # print('minilevel na wejsciu', minilevel)
    # print('')
    # print('p', p)
    # print('pair_status', pair_status)
    # print('pair_status_idx', pair_status_idx)
    # print('')
    #print('excl_coef_level na wejsciu', excl_coef_level)
    #print('')
    # print('zaczynam od znalezienia pair-restr')
    pair_restr = find_restr(minilevel, len_coef, excl_coef_level)
    # print('minilevel na wejsciu', minilevel)
    # print('PAIR RESTRRRR', pair_restr)
    # print('')
    e_cost = idx_to_cost(e_idx)
    #print('e_cost', e_cost)
    #print('e_idx', e_idx)
    #print('')
    # print('teraz szukam listy mozliwych par dla kolejnego lewelu')
    t1 = time.time()
    pair_list, pair_cost, pair_contr = find_pairs(e_idx, pair_restr, idx_matrix, len_idx, len_coef, e_cost, excl_coef_level, excl_idx, COST_THRESH)
    t2 = time.time()
    print('factopair', t2-t1)
    # print('pair_list', pair_list)


    # print('teraz szukam leveli')
    level_pairs = find_level_pairs(pair_list, len_coef)
    # print('oto wynik w postaci tabelki')
    # for x in level_pairs:
    #     print(x)
    #     print('')
    # print('koniec tabelki')

   #print('level pairs', level_pairs)

#    print('minilevel na wejsciu', minilevel)
    # print('pair_status', pair_status_idx, pair_status)

    # pair status mowi ktora piramidke z level pairs aktualnie sprawdzamy

    if (p == 0):
        pair_status.append([0, len(level_pairs)])
        big_idx_matrix.append(deepcopy(idx_matrix))
        big_idx_dict.append(deepcopy(idx_dict))
        nn[5]+=1
        # print('TO JEST PAIR STATUS', pair_status)

    # print('jestem na poziomie p=', p, len(big_idx_matrix))
    # for nn in range(0, len(big_idx_matrix)):
    #     print(big_idx_matrix[nn])
    #     print(big_idx_dict[nn])
    # vspace(1)


    if (len(level_pairs)) == 0:    #LABEL1, no ELSEIF
        # print('')
        # print('tak len jest zero')
        # print('')
        # print('pair status', pair_status)
        if pair_status[-1][0]==pair_status[-1][1]:
            # print('pair_status[-1][0]==pair_status[-1][1]',pair_status[-1][0], pair_status[-1][1])
            if len(pair_status) == 1:

                hhh = 1
                # print('hhhhhhaaaaaa')
            else:
                # print('elsse')
                big_idx_matrix = big_idx_matrix[:-1]
                
                big_idx_dict = big_idx_dict[:-1]
                d = last_row_cost(idx_matrix)
                # print(e_idx)
                # print(idx_matrix)
                lrp = last_row_parameters(idx_matrix, e_idx)
                miniparameters[-1].append(lrp)
                # print('lrp', lrp)
                # print('miniparam0', miniparameters)
                mp = max_parameters(miniparameters)
                # print('mmmppp', mp)
                miniparameters[-1].append(mp)

                minicost[-1].append(d)
                #print('o tutaj minicost', minicost)
                maxcost = max_cost(minicost)
                minicost[-1].append(maxcost)
                maxmem = mp[3]
                # print('maxmem', maxmem, MEM_THRESH)
                # print('larampim', maxcost, COST_THRESH, ORIGINAL_COST, minicost, minilevel)
                if disconnected:
                    condition = (maxcost <= ORIGINAL_COST)
                else:
                    condition = (maxcost < ORIGINAL_COST)                                                    

                if (maxcost <= COST_THRESH) and condition:
                    # print('tak mniejsze od cost')
                    if maxmem <= MEM_THRESH:
                        # print('tak mniejsze od mem')
                        #print(maxcost, COST_THRESH, maxcost<=COST_THRESH)
#                    if 1==1:
                        #print('append1', minilevel)
                        biglevel.append(deepcopy(minilevel))
                        bigcost.append(deepcopy(minicost))
                        nn[6]+=1
                        # print('miniparam1', miniparameters)
                        bigparameters.append(miniparameters)
                        # print('bigl', biglevel)
                        # print('bigp', bigparameters)

                        #print('macierze idx_matrix 1')
                        #for ii in range(0, len_coef):
                            #print(idx_matrix[ii][:])

                    else:
                        qq=0
#                        print('nie append1 a')
                else:
                    qq=0
 #                   print('nie append1 b', maxcost, COST_THRESH, ORIGINAL_COST)


                #print('minilevel przed', minilevel)
                minilevel = minilevel[:-1]
                excl_coef_level = excl_coef_level[:-1]
                #print('minilevel po', minilevel)
                minicost = minicost[:-1]
                miniparameters = miniparameters[:-1]
                #print('pair status przed', pair_status)
                pair_status = pair_status[:-1]
                #print('pair status po', pair_status)
                p = pair_status[pair_status_idx-1][0]
                new_len_coef = big_idx_matrix[-1].shape[0]
                #print('zaczynam factorize znow')
                #print('macierze idx_matrix 2')
#                for ii in range(0, len_coef):
                    #print(idx_matrix[ii][:])

                return factorize3(nn, p, pair_status, pair_status_idx-1,  pair_restr, biglevel, minilevel, bigparameters, miniparameters,\
                                      bigcost, minicost, \
                                      big_idx_dict[-1], idx_dict, big_idx_matrix[-1], \
                                      idx_matrix, len_idx, new_len_coef, e_idx, big_idx_matrix, big_idx_dict, excl_coef_level, excl_idx, ORIGINAL_COST, COST_THRESH, disconnected=disconnected)

        else:
            pair_status[-1][0] +=1
            new_len_coef = idx_matrix.shape[0]
            minilevel = minilevel[:-1]
            excl_coef_level = excl_coef_level[:-1]
            minicost = minicost[:-1]
            miniparameters = miniparameters[:-1]
            new_len_coef = big_idx_matrix[-1].shape[0]
            p = pair_status[-1][0]
            #print('macierze idx_matrix 3')
#            for ii in range(0, len_coef):
                #print(idx_matrix[ii][:])

            return factorize3(nn, p, pair_status, pair_status_idx-1,  pair_restr, biglevel, minilevel, bigparameters, miniparameters, \
                                  bigcost, minicost, \
                                  big_idx_dict[-1], idx_dict, big_idx_matrix[-1], \
                                  idx_matrix, len_idx, new_len_coef, e_idx, big_idx_matrix, big_idx_dict, excl_coef_level, excl_idx, ORIGINAL_COST, COST_THRESH, disconnected=disconnected)



    if p == len(level_pairs): #LABEL2
        
        if pair_status_idx == 0:
            hhh = 1
            # print('kuaaaa')
        else:
            # print('sraa')
            big_idx_matrix = big_idx_matrix[:-1]
            big_idx_dict = big_idx_dict[:-1]
            minilevel = minilevel[:-1]
            excl_coef_level = excl_coef_level[:-1]
            minicost = minicost[:-1]
            miniparameters = miniparameters[:-1]
            pair_status = pair_status[:-1]
            p = pair_status[pair_status_idx-1][0]
            new_len_coef = big_idx_matrix[-1].shape[0]
            #print('factordupa3')
            #print('macierze idx_matrix 4')
#            for ii in range(0, len_coef):
                #print(idx_matrix[ii][:])

            return factorize3(nn, p, pair_status, pair_status_idx-1,  pair_restr, biglevel, minilevel, bigparameters, miniparameters, \
                                  bigcost, minicost, \
                                  big_idx_dict[-1], idx_dict, big_idx_matrix[-1], \
                                  idx_matrix, len_idx, new_len_coef, e_idx, big_idx_matrix, big_idx_dict, excl_coef_level, excl_idx, ORIGINAL_COST, COST_THRESH, disconnected=disconnected)
            

    else:   #LABEL2
        # print('atu')

        for i in range(p, len(level_pairs)):
            # print('dla lewelu w zakresie', p, len(level_pairs))
            level_cost = []
            level_parameters = [[], [], [],[]]
            # print('minilevel', minilevel)
            # print('pair_status', i, pair_status)
            
            # print('level_pairs', level_pairs)
            # print('level_pairs[i], pair_status_idx', level_pairs[i], pair_status_idx)
            pair_status[pair_status_idx] = [i+1, len(level_pairs)]
            #print('pair status juz jest', pair_status)
            #print('minil przed append', minilevel)

            minilevel.append(level_pairs[i])
            # print('minil po append', minilevel)

            excl_coef_level = append_excl_coef_level(minilevel, excl_coef_level)
            #print('excl_coef_level po append', excl_coef_level)
            idx_matrix_temp = deepcopy(idx_matrix)
            nn[7]+=1

            del_list = []
            level_pairs_cp = deepcopy(level_pairs)
            nn[8]+=1
            # #print('len(level_pairs_cp)', len(level_pairs_cp))
            # #print('czy wchodze')
            #print('macierze idx_matrix 5')
#            for ii in range(0, len_coef):
                #print(idx_matrix[ii][:])

            

            idx_dict_work = deepcopy(idx_dict)
#            print('nnnn', nn)
            nn[9]+=1
            t2 = time.time()
            for enum in range(0, len(level_pairs_cp[i])):
                # print('to miejscw w ktorym chyba zmieniam macierz idx')
                # print('tak', enum)
                elem = level_pairs_cp[i][enum]
                #print('czesc level_pairs', elem)
                if len(elem) == 2:
                    #print('dlugosc elem =2')
                    ii = elem[0]-1
                    jj = elem[1]-1
                    #print('biore elementy ii oraz jj')

                    d = sum(np.bitwise_or(idx_matrix_temp[ii][:], idx_matrix_temp[jj][:]))
                    #print('biore wektory')
                    #print('minilevel', minilevel)
                    #print(idx_matrix_temp[ii][:])
                    #print(idx_matrix_temp[jj][:])
                    # #print('tego szukasz', e_idx, excl)
                    cost_n, cost_real, mem_n, mem_real = detailed_cost_mem(idx_matrix_temp[ii][:], idx_matrix_temp[jj][:], e_idx, excl_idx)
                    # print('detailed cost', cost_n, cost_real, mem_n, mem_real)
                    level_cost.append(d)
                    # #print('')
                    # #print('level_parameeters')
                    # #print('')
                    level_parameters[0].append(cost_n)
                    level_parameters[1].append(cost_real)
                    level_parameters[2].append(mem_n)
                    level_parameters[3].append(mem_real)
                    c = np.bitwise_xor(idx_matrix_temp[ii][:], idx_matrix_temp[jj][:])
                    #print('wektor c jest suma ich', c)
                    idx_matrix_temp[ii][:] = deepcopy(c)
                    nn[10]+=1
                    del_list.append(jj)
                    idx_matrix_temp = np.delete(idx_matrix_temp, jj, 0)

                    new_len_coef = idx_matrix_temp.shape[0]
                    #print('czy zmieniam slownik', idx_dict_work, 'ELEM', elem, ii, jj)
                    #print('jj', jj)
                    idx_dict_temp = {}
                    #print('macierze idx_matrix original 6')
#                    for ii in range(0, len_coef):
                        #print(idx_matrix[ii][:])
                    #print('')
                    #print('macierze idx_matrix with deleted row 6')
#                    for ii in range(0, new_len_coef):
                        #print(idx_matrix_temp[ii][:])

                    for enumrest in range(enum+1, len(level_pairs_cp[i])):
                        rest_elem = level_pairs_cp[i][enumrest]
                        #print('rest_elem', rest_elem)
                        if len(rest_elem) == 2:
                            if rest_elem[0]-1 >= jj:
                                #print('rest elem zero zmieniam')
                                level_pairs_cp[i][enumrest][0] -=1
                            if rest_elem[1]-1 >=jj:
                                #print('rest elem jeden zmieniam')
                                level_pairs_cp[i][enumrest][1] -=1
                    #print('zmieniam? oto temp', idx_dict_temp)
                    #print('oraz teraz bedzie nowy elem')
                    idx_dict = deepcopy(idx_dict_temp)
                    nn[11]+=1
            t3 = time.time()
#            print(t3-t2, 'time nn10')
            minicost.append(level_cost)
            # print('levvvel parameters', level_parameters)
            miniparameters.append(level_parameters)
            #print('minicost', minicost)
            # print('zminiparameters', miniparameters)
            # #print('tutaj-srrrrrrr')
            #print('macierze idx_matrix original 7')
 #           for ii in range(0, len_coef):
                #print(idx_matrix[ii][:])
            #print('')
            #print('macierze idx_matrix deleted row 7')
#            for ii in range(0, new_len_coef):
                #print(idx_matrix_temp[ii][:])


            if new_len_coef <= 2:
                # print('mniejsze od 2 i pair status', pair_status)
                if pair_status[-1][0] ==pair_status[-1][1]:
                    #print('i jeszcze tego', idx_matrix_temp)
                    lrp = last_row_parameters(idx_matrix_temp, e_idx)
                    miniparameters[-1].append(lrp)
                    #print('miniparam', miniparameters)
                    mp = max_parameters(miniparameters)
                    miniparameters[-1].append(mp)
                    #print('minil', minilevel)
                    #print('idxmt', idx_matrix_temp)
                    
                    d = last_row_cost(idx_matrix_temp)
                    minicost[-1].append(d)
                    maxcost = max_cost(minicost)
                    minicost[-1].append(maxcost)
                    maxmem = mp[3]
                    #print('macierze idx_matrix 8')
#                    for ii in range(0, len_coef):
                        #print(idx_matrix[ii][:])

                    #print('larampim2', maxcost, COST_THRESH, ORIGINAL_COST, minicost, minilevel)
                    #print('maxmem', maxmem)
                    #print('d', d)

                    
                    if disconnected:
                        condition = (maxcost <= ORIGINAL_COST)
                    else:
                        condition = (maxcost < ORIGINAL_COST)                            
                        
                    if (maxcost <= COST_THRESH)  and condition:
                        # print('tak cost')
                        if maxmem <= MEM_THRESH:
                            # print('tak mem')
                            #                        if 1==1:
                            #print('append2', minilevel)
                            biglevel.append(deepcopy(minilevel))
                            bigcost.append(deepcopy(minicost))
                            nn[11]+=1
                            # print('miniparam2', miniparameters)
                            bigparameters.append(miniparameters)
                            # print('bigl1', biglevel)
                            #print('macierze idx_matrix 9')
                            #for ii in range(0, len_coef):
                                #print(idx_matrix[ii][:])

                        else:
                            qq=0
#                            print('nie append2 a', maxmem, MEM_THRESH)
                    else:
                        qq=0
 #                       print('nie append2 b', maxcost, COST_THRESH, ORIGINAL_COST)   


                    big_idx_matrix = big_idx_matrix[:-1]
                    big_idx_dict = big_idx_dict[:-1]

                    minilevel = minilevel[:-2]
                    excl_coef_level = excl_coef_level[:-2]
                    minicost = minicost[:-2]
                    miniparameters = miniparameters[:-2]
                    pair_status = pair_status[:-1]

                    #print('macierze idx_matrix 10')
#                    for ii in range(0, len_coef):
                        #print(idx_matrix[ii][:])

                    if pair_status_idx == 0:
                        hhh = 1
                    else:
                        p = pair_status[pair_status_idx-1][0]
                        new_len_coef = big_idx_matrix[-1].shape[0]
                        #print('new_len_coef', new_len_coef)
                        #print('factordupa4')
                        #print('macierze idx_matrix 11')
 #                       for ii in range(0, len_coef):
                            #print(idx_matrix[ii][:])
                        # print('ku1')
                        return factorize3(nn, p, pair_status, pair_status_idx-1,  pair_restr, biglevel, minilevel, bigparameters, miniparameters,\
                                              bigcost, minicost, \
                                              big_idx_dict[-1], idx_dict, big_idx_matrix[-1], \
                                              idx_matrix, len_idx, new_len_coef, e_idx, big_idx_matrix, big_idx_dict, excl_coef_level, excl_idx, ORIGINAL_COST, COST_THRESH, disconnected=disconnected)

                else:
                    lrp = last_row_parameters(idx_matrix_temp, e_idx)
                    miniparameters[-1].append(lrp)
                    mp = max_parameters(miniparameters)

                    miniparameters[-1].append(mp)

                    d = last_row_cost(idx_matrix_temp)
                    minicost[-1].append(d)
                    maxcost = max_cost(minicost)
                    maxmem = mp[3]
                    minicost[-1].append(maxcost)
                    #print('larampim3', maxcost, COST_THRESH, ORIGINAL_COST, minicost, minilevel)
                    if disconnected:
                        condition = (maxcost <= ORIGINAL_COST)
                    else:
                        condition = (maxcost < ORIGINAL_COST)                                                    

                    if (maxcost <= COST_THRESH)  and condition:
                        if maxmem <= MEM_THRESH:
                            #print('append3', minilevel)
                            biglevel.append(deepcopy(minilevel))
                            bigcost.append(deepcopy(minicost))
                            nn[12]+=1
                            # print('miniparam3',miniparameters)
                            bigparameters.append(miniparameters)
                            #print('macierze idx_matrix 12')
  #                          for ii in range(0, len_coef):
                                #print(idx_matrix[ii][:])

                        else:
                            qq=0
                            # print('nie append3 a')
                            # print('maxmem, > MEM_THRESH', maxmem, '>', MEM_THRESH)
                    else:
                        qq=0
                        # print('nie append3 b')   
                            
                    minilevel = minilevel[:-1]
                    excl_coef_level = excl_coef_level[:-1]
                    minicost = minicost[:-1]
                    miniparameters = miniparameters[:-1]
            else:
                pair_status_idx += 1

                #print('macierze idx_matrix 13')
 #               for ii in range(0, len_coef):
                    #print(idx_matrix[ii][:])
                #print('')

#                for ii in range(0, new_len_coef):
                    #print(idx_matrix_temp[ii][:])

                # print('ku2')
                return factorize3(nn, 0, pair_status, pair_status_idx,  pair_restr, biglevel, minilevel, bigparameters, miniparameters, \
                                      bigcost, minicost, \
                                      idx_dict_temp, idx_dict, idx_matrix_temp, \
                               idx_matrix, len_idx, new_len_coef, e_idx, big_idx_matrix, big_idx_dict, excl_coef_level, excl_idx, ORIGINAL_COST, COST_THRESH, disconnected=disconnected)



def detailed_cost_mem(a, b, e_idx, excl_idx):

    e_cost = idx_to_cost(e_idx)
    e_cost2 = idx_to_cost(excl_idx)

    e_cost = e_cost+e_cost2
    
    cost_n = sum(np.bitwise_or(a, b))
    c = np.bitwise_or(a, b)

    cost_real = 1
    # print('witaram', c)
    # print(e_cost)
    for h in range(0, len(c)):
        if c[h] == 1:
            cost_real *= e_cost[h]                                       

    mem_n = sum(np.bitwise_xor(a,b))
    d = np.bitwise_xor(a,b)
    
    mem_real = 1
    for h in range(0, len(d)):
        if d[h] == 1:
            mem_real *= mem_real_cost[e_idx[h]]

    mem_real = mem_real * 8/(10**9)
#    if mem_real > 1.0:
        #print('mem-wieksze', mem_real)

    return cost_n, cost_real, mem_n, mem_real

def append_excl_coef_level(minilevel, excl_coef_level):

    #print('')
    #print('wewnatrz appedn')
    if len(excl_coef_level)>=1:
        #print('tak len coef >=1')
        for i in range(0, len(minilevel[-1])):
            #print(minilevel[-1][i], excl_coef_level[-1])
            if minilevel[-1][i][0]==excl_coef_level[-1]:
                excl_coef_level.append(i+1)

    return excl_coef_level

def find_restr(minilevel, len_coef, excl_coef_level):

    idx_dict = {}
    idx_dict_len = 0
    del_list = []
    if len(minilevel) > 0:
        minil = minilevel[-1]
        for j in range(0, len(minil)):
            idx_dict_len += len(minil[j])
            if len(minil[j])==2:
                del_list.append(minil[j][1])
                

    for i in range(1, idx_dict_len+1):
        idx_dict[i] = i
        ##print('testuje', i)
        if idx_dict[i] not in del_list:
            ##print('del list', del_list)
            for j in del_list:
                ##print('sprawdzam j', j)
                if i > j:
                    ##print('zmniejszm')
                    idx_dict[i] = idx_dict[i] -1
                    ##print(idx_dict[i])
                ##print('')

                    
    pair_restr = [[],[]]
    ##print('')
    ##print('minilevel z find_restr')
    ##print(minilevel)
    ##print('len_coef', len_coef)
    ##print(idx_dict)
    ##print('')
    if len(minilevel) > 0:
        minil = minilevel[-1]
        ##print('******dctdct******')
        ##print(len(excl_coef_level), excl_coef_level)
        ##print(minilevel[-1])
        ##print(idx_dict)
        ##print('minil lala', minil)
        for j in range(0, len(minil)):
            x = minil[j]
            if len(x) == 2:
                ##print('dodaje pierwszy z podwojnego', idx_dict[x[0]])
                # #print('dodaje', x, x[0], idx_dict[x[0]])
                # #print('nowe dodaje', j+1, x[0], idx_dict[x[0]])
                # pair_restr.append(x[0])
                # pair_restr[0].append(j+1)
                pair_restr[0].append(idx_dict[x[0]])
                #pair_restr[0].append(j+1)                                
            else:                
                # if (x[0] != excl_coef_level[-1]):
                if excl_coef_level != []:
                    if (idx_dict[x[0]] != excl_coef_level[-1]):
                        ##print('dodaje pojedynczy', idx_dict[x[0]])                        
                        pair_restr[1].append(idx_dict[x[0]])
                    # if ((j+1) != excl_coef_level[-1]):
                    #     pair_restr[1].append(j+1)#idx_dict[x[0]])
                else:
                    ##print('dodaje pojedynczy', idx_dict[x[0]])                        
                    pair_restr[1].append(idx_dict[x[0]])
                    # pair_restr[1].append(j+1)
                # pair_restr[1].append(j+1)

    else:
        ##print('maly1')
        for i in range(0, len_coef):
            ##print('maly2')
            if excl_coef_level != []:
                ##print('tak rozne od 0', i, i+1, excl_coef_level, excl_coef_level[0])
                if (i+1) != excl_coef_level[0]:
                    ##print('dodaje pojedynczy excl', i+1)                        
                    pair_restr[0].append(i+1)
                    pair_restr[1].append(i+1)
            else:
                ##print('dodaje pojedynczy nie excl', i+1)                        
                pair_restr[0].append(i+1)
                pair_restr[1].append(i+1)
        ##print('else', pair_restr)
    #print('ostatecznie', pair_restr)
    #print('')
    return pair_restr


def last_row_parameters(idx_matrix, e_idx):

    e_cost = idx_to_cost(e_idx)
    c = []
    d = []
    v1 = idx_matrix[00][:]
    v2 = idx_matrix[00][:]
    #print(idx_matrix)
    for tt in range(1, idx_matrix.shape[0]):
        c = np.bitwise_or(v1, idx_matrix[tt][:])
        d = np.bitwise_xor(v2, idx_matrix[tt][:])
        v1 = deepcopy(c)
        v2 = deepcopy(d)
        
        


    cost_n = sum(c)
    mem_n = sum(d)

    cost_real = 1
    for h in range(0, len(c)):
        if c[h] == 1:
            cost_real *= e_cost[h]


    mem_real = 1
    for h in range(0, len(d)):
        if d[h] == 1:
            mem_real *= mem_real_cost[e_idx[h]]#e_cost[h]                                                                                                

    mem_real = mem_real * 64/(10**9)
#    #print('memreallast', mem_real)


    return [cost_n, cost_real, mem_n, mem_real]

def max_parameters(miniparameters):
#    print('szukam max')
    max_cost_n = 0
    for x in miniparameters:
#        print('')
#        print('xx', x, x[0])
        for y in x[0]:
            if y > max_cost_n:
                max_cost_n = y

    max_cost_real = 0
    for x in miniparameters:
        for y in x[1]:
            if y > max_cost_real:
                max_cost_real = y

    max_mem_n = 0
    for x in miniparameters:
        for y in x[2]:
            if y > max_mem_n:
                max_mem_n = y

    max_mem_real = 0
    for x in miniparameters:
        for y in x[3]:
            if y > max_mem_real:
                max_mem_real = y


    # teraz jeszcze musze sprawdzic last row parameters ktore zawsze sa zapisane w
    # maxparameters
#    print('musze sprawdzic last row parameters', miniparameters[-1][-1][0])
    if miniparameters[-1][-1][0] > max_cost_n:
        max_cost_n = miniparameters[-1][-1][0]
    if miniparameters[-1][-1][1] > max_cost_real:
        max_cost_real = miniparameters[-1][-1][1]
    if miniparameters[-1][-1][2] > max_mem_n:
        max_mem_n = miniparameters[-1][-1][2]
    if miniparameters[-1][-1][3] > max_mem_real:
        max_mem_real = miniparameters[-1][-1][3]


    return [max_cost_n, max_cost_real, max_mem_n, max_mem_real]

    
def last_row_cost(idx_matrix):

    v1 = idx_matrix[00][:]
    # #print('LAST row cost')
    # #print(v1)
    for tt in range(1, idx_matrix.shape[0]):
        # #print(v1)
        # #print(idx_matrix[tt][:])
        v1 = np.bitwise_or(v1, idx_matrix[tt][:])
        # #print('v1po', v1)
    d = sum(v1)
    # #print('i taki wychodzi cost', d)

    return(d)

def max_cost(minicost):

    maxcost = 0
    for x in minicost:
        for y in x:
            if y > maxcost:
                maxcost = y

    return(maxcost)
                    


def find_level_pairs(pair_list, len_coef):


    level_pairs = []
    n = len_coef//2
    npair = len(pair_list)
    big_lp = []
    # #print('pair_list', pair_list)
    # #print('robie sublevel', len_coef)
    # #print('range to jest', 1, n+1)
    for i in range(1, n+1):
        # #print('')
        # #print('--------------', i, ' i teraz juz jest TO JE ONO', i)
        # #print('')
        lp = sublevel(pair_list, npair, i, len_coef)
        # #print('lp poooo', lp)
        for x in lp:
            big_lp.append(x)
    # #print('big_lp', big_lp)
    return big_lp
                        
def sublevel(pair_list, npair, n, len_coef):

    big_lp = []

    # #print('pair_list', pair_list, n)
    comb_set = list(combinations(pair_list, n))
    # #print('comb_set', comb_set)
    comb_list = []
    for elem in comb_set:
        comb_list.append(list(elem))
    # #print('kombinacja wszystkich', comb_list)
    
    comb_list_unique = []
    for elem in comb_list:
        add = True
        elemlist  = []
        for minielem in elem:
            elemlist = elemlist + minielem

        elemlist.sort()
        for e in elemlist:
            count_e = elemlist.count(e)
            if count_e > 1:
                add = False

        if add == True:
            comb_list_unique.append(elem)
    # #print('kombinacja unique', comb_list_unique)
    for lp in comb_list_unique:        
        elemset  = set()
        for minielem in lp:
            # #print('minielem', minielem, set(minielem))        
            elemset = elemset |  set(minielem)
        # #print('elemset', elemset)
        for sing_elem in range(1, len_coef+1):
            # #print('sprawdzam czy', sing_elem, 'jest w ', elemset)
            if sing_elem not in elemset:
                # #print('dodaje sing', sing_elem)
                lp.append([sing_elem])
        # #print('oto lp', lp)
        big_lp.append(lp)


    # #print('pair_list w sublevel', pair_list)
    # for k1 in range(0, npair):
    #     c = pair_list[k1][0]
    #     d = pair_list[k1][1]
    #     lp = [pair_list[k1]]
    #     #print('lp dla k1k1=', k1, lp)
    #     idx_set = {c, d}
    #     #print('idx_set', idx_set)
    #     #print('nnn', len(lp), n)
    #     pair_counter[0] += 1
    #     mm = 1
    #     if len(lp)<n:
    #         #print('if')
    #         for k2 in range(k1+1, npair):
    #             a = pair_list[k2][0]
    #             b = pair_list[k2][1]
    #             #print('ab', a, b)
    #             #print('sprawdzam z para', pair_list[k2])
    #             #print('idx_set, a, b', idx_set, a, b)
    #             if a not in idx_set and b not in idx_set:
    #                 #print('oba nie sa')
    #                 lp.append(pair_list[k2])
    #                 pair_counter[mm] = k2
    #                 mm += 1
    #                 #print('lp', lp)
    #                 idx_set.add(a)
    #                 idx_set.add(b)
    #                 #print('sprawdzam', len(lp), n)
    #                 if len(lp) == n:
    #                     #print('tak, juz max')
    #                     lpcp = deepcopy(lp)
    #                     for sing_elem in range(1, len_coef+1):
    #                         if sing_elem not in idx_set:
    #                             lp.append([sing_elem])
    #                     big_lp.append(lp)
    #                     lp = lpcp[:-1]
    #                     idx_set.remove(a)
    #                     idx_set.remove(b)
    #                     pair_counter[n] = 0
    #                     if pair_counter[n-1] == npair:
    #                         pair_counter[n-1] = 0
    #                     #print('lp tutaj po usunieciu ostatniego', lp)
    #             else:
    #                 #print('ktorys jednak jest')
    #     else:
    #         #print('else')
    #         for sing_elem in range(1, len_coef+1):
    #             if sing_elem not in idx_set:
    #                 lp.append([sing_elem])
    #                 #print('lp w elese', lp)
    #         big_lp.append(lp)


    return big_lp

def real_contractions(i, j, n, k, idx_matrix, c):

    add = True
    vec = []
    for l in range(0, n):
        vec.append(0)

    for l in range(0, k):
        if (l != i and l!= j):
            b = np.bitwise_and(c, idx_matrix[l][0:n])
            if sum(b) > 0:
                vec = deepcopy(np.bitwise_or(vec, b))
    if np.array_equal(vec, c):
        add = False

    return add


def find_pairs(e_idx, restr, idx_matrix, n, k, e_cost, excl_coef_level = [], excl_idx = [], COST_THRESH = 0):


    excl_idx_to_no = []
    for i in range(0, len(e_idx)):
        if e_idx[i] in excl_idx:
            excl_idx_to_no.append(i+1)

    
    #print('')
    #print('zaczynam szukanie par')
    perm_list = []
    perm_cost = []
    perm_contr = []
    #print('')
    #print('find pairs')
    #print('restr', restr)
    #print('excluded coefs', excl_coef_level)
    #print('excluded_idx', excl_idx)
    #print('excluded_idx_to_no', excl_idx_to_no)
    for i in range(0, k):
        for j in range(i+1, k):
            #print('----------sprawdzam parę------------------', i+1, j+1)
            if i != j:
                a = idx_matrix[i][0:n]
                b = idx_matrix[j][0:n]
                a_excl = idx_matrix[i][:]
                b_excl = idx_matrix[j][:]
                zz = np.bitwise_and(a,b)
#                vspace(0)
                #print('bitwise')
                #print(a, b)
                #print(zz)
 #               vspace(0)
                c = sum(np.bitwise_and(a,b))
                #print('paraa', a)
                #print('parab', b)
                if c > 0:
                    d = sum(np.bitwise_or(a, b))
                    dd = np.bitwise_or(a, b)
                    #print('bitor', dd)
                    d_excl = sum(np.bitwise_or(a_excl, b_excl))
                    if d_excl <= COST_THRESH: 
                        # sprawdz czy indeksy sumacyjne nie wystepuja w zadnych innych wyrazach
                        cc = np.bitwise_and(a, b)
                        add = real_contractions(i, j, n, k, idx_matrix, cc)
                        if add == True:
                            #print('wymogi formalne sa')
                            dd = np.bitwise_xor(a,b)
                            # #print('restrrrrr', restr, i+1, j+1)
                            #print('czy oba ind naleza do restr?', restr)
                            if ((i+1) in restr[0] and (j+1) in restr[1]) or ((i+1) in restr[1] and (j+1) in restr[0]):
                                #print('tak')
                                # #print('restr', restr, i+1, j+1, excl)
                                # #print(i+1, excl)
                                # #print(j+1, excl)
                                # #print('czy oba sa rozne od exlc_coef?')
                                # if (i+1) not in excl_coef_level and (j+1) not in excl_coef_level:
                                    # #print('tak')
                                    # if (i+1) not in excluded_idx and (j+1) not in excluded_idx:
                                    # if 1 == 1:
                                perm_list.append([i+1,j+1])
                                perm_cost.append([d_excl])
                                mem = 1
                                for h in range(0, len(dd)):
                                    if dd[h] == 1:
                                        # #print('e_cost[h]', e_cost[h])
                                        mem *= e_cost[h]
                                perm_contr.append([mem])

    #print('')
    return perm_list, perm_cost, perm_contr




        
#-----------------------------------------Part for translating ugg with intermediates-----------------------------------------------


def read_pyramid(minilevel, e, worm, outer_worm, idxfx_worm, interm_dict,
                 all_hash, xfx_dict, interm_fx, interm_hash, ni, disc_mini, \
                 excl = [], COST_THRESH = COST_THRESH_DEFAULT, disconnected = False):

    #disconnected out says whether there is disconnected term in this pyramid

#    print('to jest e', e)
    disconn_out = False
    # print('to jest minilevel', minilevel)
    # print('to jest worm')
    # for x in worm:
    #     print(x)
    # print('to jest outer_worm')
    # for x in outer_worm:
    #     print(x)
    # print('idxfx_worm', idxfx_worm)
    # print('interm_dict', interm_dict)
    # print('all hash', all_hash)
    # print('xfx_dict', xfx_dict)
    # print('interm_fx', interm_fx)
    # print('ni', ni)
    
    e_work = deepcopy(e)
    miniworm = []
    miniworm_xfx = []

    for x in range(0, len(e.coefficient)):
        if e.coefficient[x] == OBSERVABLE_X or e.coefficient[x] == OBSERVABLE_X_ASYM:
            xfx = e.coefficient_idx[x]
            break
        else:
            xfx = []

    used = []
    for i in range(0, len(minilevel)):
        s_idx = []

        minicp = deepcopy(minilevel[i])
#        print('minicp', minicp)
        for enum in range(0, len(minicp)):
            elem = minicp[enum]
            if len(elem) == 2:
 #               print('elem', elem)
                pair = [deepcopy(elem[0])-1, deepcopy(elem[1])-1]

    #            print('bede extractowac z ', e_work)
                e_work.clear_fixed()
                e_work.establish_fixed()
                interm1, new_e, cn, ni, interm_fixed_old, interm_fixed, xfx_out = \
                    extract_interm(pair, e_work, xfx, all_hash, xfx_dict, ni, s_idx, excl, disconnected)
#                print('i wynik to interm1', interm1)
                # print('new_e', new_e)
                # print('cn', cn)
                # print('ni', ni)
                # print('interm_fixed_old', interm_fixed_old)
                # print('interm_fixed', interm_fixed)
                # print('xfx_out', xfx_out)
                
                s_idx = s_idx + interm_fixed_old
                if interm_fixed != []:
                    new_e.coefficient[pair[0]] = cn
                    new_e.coefficient_idx[pair[0]] = deepcopy(interm_fixed_old)
                else:
                    if disconnected:
#                        print('ten jest disconnected')
                        disconn_out = True
                        disc_mini[i] += 1
                        new_e.coefficient[pair[0]] = cn
                        new_e.coefficient_idx[pair[0]] = []
                    else:
                        new_e = ugg()

                interm_fx[interm1.binary_hash] = deepcopy(interm_fixed)

                
                e_work = deepcopy(new_e)
                miniworm.append(interm1)
  #              print('appenduje do miniworma', interm1)

                miniworm_xfx.append(interm_fixed)#xfx_out)
   #             print('appenduje do miniworma_xfx', interm_fixed)
                
                for enumrest in range(enum+1, len(minicp)):
                    rest_elem = minicp[enumrest]
                    if len(rest_elem) == 2:
                        if rest_elem[0]-1 >= pair[1]:
                            minicp[enumrest][0] -=1
                        if rest_elem[1]-1 >=pair[1]:
                            minicp[enumrest][1] -=1

    # print('miniworm')
    # for x in miniworm:
    #     print(x)
    # print('')
    add = True
    for d in miniworm:
        # print('d in miniworm', d)
        temp = []
        cost2 = 0
        for dd in d.coefficient:
            if dd == OBSERVABLE_X or dd == OBSERVABLE_X_ASYM:
                sys.exit(0)
                add = False
        for dd in d.coefficient_idx:
            for ddd in dd:
                if ddd not in temp:
                    cost2 += 1
                    temp.append(ddd)
        if cost2 > COST_THRESH:
            print(cost2, COST_THRESH)
            add = False
    if add == True:
        for interm1 in miniworm:
            if interm1.binary_hash in interm_dict:
#                vspace(0)
                # print('wyrazz e ', e)
                # print('ma TEŻ ten interm', interm1, all_hash[interm1.binary_hash])
                interm_dict[interm1.binary_hash] += 1
                # print('i interm_dict[interm1.binary_hash]', interm_dict[interm1.binary_hash])
                # print('interm is here', all_hash[interm1.binary_hash],  interm_hash[interm1.binary_hash])
 #               vspace(0)
            else:
  #              vspace(0)
                # print('wyrazz e ', e)
                # print('ma ten interm', interm1, all_hash[interm1.binary_hash])
                interm_dict[interm1.binary_hash] = 1
                # print('i interm_dict[interm1.binary_hash]', interm_dict[interm1.binary_hash])
                interm_hash[interm1.binary_hash] = interm1
                #             vspace(0)

        
    if add == True:
#        print('appenduje do worma', miniworm)
        worm.append(miniworm)
#        print('appenduje do outer worma', e_work)
        outer_worm.append(deepcopy(e_work))
#        print('appenduje do idxfx', miniworm_xfx)
        idxfx_worm.append(miniworm_xfx)
        if len(e_work.coefficient) == 0:
            sys.exit(0)
    else:
        print('add bylo false')

            
    return ni, disconn_out


def extract_interm(pair, e, xfx, all_hash, xfx_dict, ni, s_idx = [], excl = [], disconnected = False):

    
    interm = ugg()    
    new_e = ugg()
    new_e.num_factor = e.num_factor

    for i in range(0, len(e.coefficient)):
        if i == pair[0]:
            interm.coefficient.append(e.coefficient[i])
            interm.coefficient_idx.append(e.coefficient_idx[i])
            new_e.coefficient.append('')
            new_e.coefficient_idx.append([])
        elif i == pair[1]:
            interm.coefficient.append(e.coefficient[i])
            interm.coefficient_idx.append(e.coefficient_idx[i])
        else:
            new_e.coefficient.append(e.coefficient[i])
            new_e.coefficient_idx.append(e.coefficient_idx[i])


    # #print('e', e, xfx)
    #print('interm-yyy', interm)

    interm_fixed = find_sum_interm_rest(interm, new_e, s_idx, excl)


    
    interm_fixed_old = deepcopy(interm_fixed)
#    print('interm przedppp', interm)

    xfx_interm = []
    for elem in interm.coefficient_idx:
#        print('elem', elem)
        for i in elem:
 #           print('i', i)
            if i in xfx:
                if i not in xfx_interm:
                    xfx_interm.append(i)
                
    # print('interm_fixed', interm_fixed)
    # print('interm', interm)
    # print('xfx_interm', xfx_interm)
    if disconnected:
        if len(interm_fixed) == 0:
            interm.standarize()
    else:
        interm.standarize_drag_list(interm_fixed, xfx_interm)
    
    # print('interm po standarize', interm)
    # print('interm_fixexd', interm_fixed)
    # print('xfx_interm', xfx_interm)
    # print('all_hash', all_hash)
    interm.binary_hash_gen()
    # print(interm.binary_hash)

    # all hash is a dictionary that has an interm hash
    # and corresponding name which is in cn

    # xfx_dict is a dictionary that has an interm hash
    # and corresponding fixed indices
    
    if interm.binary_hash not in all_hash:
        if disconnected:
            cn = 'TTinterm'+str(ni)
        else:
            cn = 'interm'+str(ni)
        all_hash[interm.binary_hash] =  cn
        xfx_dict[interm.binary_hash] = xfx_interm
        ni += 1
    else:
        cn = all_hash[interm.binary_hash]

        a1 = xfx_dict[interm.binary_hash]

    for idx_lst in interm.coefficient_idx:
        for idx in idx_lst:
            if idx in excl and idx not in interm_fixed_old:
                interm_fixed_old.append(idx)
            if idx in excl and idx not in interm_fixed:
                interm_fixed.append(idx)
            

    return interm, new_e, cn, ni, interm_fixed_old, interm_fixed, xfx_interm



def find_fixed_for_interm2(e):

    fx = []
    for i in range(0, len(e.coefficient_idx)):
        for k in e.coefficient_idx[i]:
            if k not in fx and k not in e.summation:
                fx.append(k)
    return fx


def find_sum_interm_rest(interm, rest, s_idx = [], excl = []):

#    #print('a to bedzie dopiero plusz')
#    print('interm', interm)
#    print('rest', rest)
#    print('excl', excl)
    interm_idx = []
    interm_fixed = []

    for x in interm.coefficient_idx:
        for y in x:
            if y not in interm_idx:
                if y not in excl:
                    interm_idx.append(y)
 #   print('if one', interm_idx)


    for x in rest.coefficient_idx:
        for y in x:
            if y not in rest.summation:
                if y not in excl:
                    rest.summation.append(y)

#    print('ii', interm_idx)
#    print('rs', rest.summation)
#    print('s_idx', s_idx)
    forb_idx = s_idx + rest.summation
#    print('forb idx', forb_idx)
    int_occ = []
    int_virt = []

    for x in interm_idx:
        if x not in forb_idx:
            interm.summation.append(x)
        else:
            if x in occupied or x in mona_occup  or x in monb_occup :
                int_occ.append(x)
            elif x in virtual or x in mona_virt or x in monb_virt:
                int_virt.append(x)
                
    interm_fixed = int_virt + int_occ


    # print('interm summ', interm.summation)
    #vspace(0)
    return interm_fixed


def find_fixed_for_interm(e):

    fx = []
    for i in range(0, len(e.coefficient_idx)):
        for k in e.coefficient_idx[i]:
            if k not in fx and k not in e.summation:
                fx.append(k)
    return fx


def latex_with_cost(f, n, rsimp_ost_in, name=False, all_hash = [], f12 = False):

    # grupuję, żeby były najpierw male sumowania, pozniej duze
    all_summation_list = []
    for x in rsimp_ost_in:
        x.summation.sort()
        if x.summation not in all_summation_list:
            all_summation_list.append(x.summation)
    all_summation_list = sorted(all_summation_list, key=lambda i: len(i))
    rsimp_ost_sort = []
    for s in all_summation_list:
        for x in rsimp_ost_in:
            if x.summation == s:
                rsimp_ost_sort.append(x)


    
    # tutaj zamieniam indeksy i> na i1 dla latexa
    rsimp_ost = deepcopy(rsimp_ost_in)    
    if f12 == True:
        for k in range(0, len(rsimp_ost)):

            subst_list_old = []
            subst_list_new = []
            for i in rsimp_ost[k].summation:
                if ">" in i:
                    if i in virtual:
                        ii = free_idx(latex_virtual, subst_list_new)
                        subst_list_old.append(deepcopy(i))
                        subst_list_new.append(deepcopy(ii))
                    if i in occupied:
                        ii = free_idx(latex_occupied, subst_list_new)
                        subst_list_old.append(deepcopy(i))
                        subst_list_new.append(deepcopy(ii))
                # if i in CABS:
                #     ii = free_idx(latex_CABS, subst_list_new)
                #     subst_list_old.append(deepcopy(i))
                #     subst_list_new.append(deepcopy(ii))

            #print("subst1", k, rsimp_ost[k])
            #print(subst_list_old)
            #print(subst_list_new)
            for j in range(0, len(subst_list_old)):
                rsimp_ost[k].substitute(subst_list_old[j], subst_list_new[j])
            #print("subst2", rsimp_ost[k])


              

    int_name = ""
    
    k = 1
    kk = 1
    #print("\equag{")

    lx = "\\begin{minipage}[t]{0.5\\textwidth}\n"
    lx += "\\begin{alignat*}{4}\n"
    # lx = "\equag{"
    for z in rsimp_ost:
        # lx += f'teraz idzie k i kk {k}, {kk}'
        fx = ""
        if name:
            namestr = all_hash[z.binary_hash]
            fx = find_fixed_for_interm(z)
            fxstr = ""
            for x in fx:
                fxstr += x
            lenl = 0
            for s in namestr:
                if s == "l":
                    lenl += 1
            no = namestr[lenl:len(namestr)]
            if lenl == 1:
                namestr = f"(\Theta_{{{no}}})" 
            elif lenl == 2:
                namestr = f"(\Xi_{{{no}}})" 
            elif lenl == 3:
                namestr = f"(\\vartheta_{{{no}}})"
            elif lenl == 4:
                namestr = f"(\\upsilon_{{{no}}})" 
            else:
                #print('implement name for higher level intermediates')
                sys.exit(0)

            int_name = f"{namestr}_{{{fxstr}}}\qquad "
        memidx, costidx = memcost(z)
        # lx+= f'w tym momencie jest {kk},{len(rsimp_ost)} {kk%33} '
        s = ""
        if (k-1) == 16:
            if (kk == len(rsimp_ost)):
                s = """"
                \\end{alignat*}
                \\end{minipage}
                \\newpage
                """
            else:
                if (kk-1)%32 == 0:#33:
                    s = """"
                    \\end{alignat*}
                    \\end{minipage}
                    \\newpage
                    \\begin{minipage}[t]{0.5\\textwidth}
                    \\begin{alignat*}{4}
                    """
                else:
                    s = """
                    \\end{alignat*}
                    \\end{minipage}
                    \\begin{minipage}[t]{0.5\\textwidth}
                    \\begin{alignat*}{4}
                    """
            k = 1
            # if kk%33 == 0:#33:
            # lx += f'SRAM1,  {len(rsimp_ost)}'
            # if (kk == len(rsimp_ost)):
            #     lx += 'SRAM2'
            #     s = """"
            #     \\end{alignat*}
            #     \\end{minipage}
            #     \\newpage
            #     """
            # else:
            #     s = """szczyna
            #     \\end{alignat*}
            #     \\end{minipage}
            #     \\newpage
            #     \\begin{minipage}[t]{0.5\\textwidth}
            #     \\begin{alignat*}{4}

            #     """
            # k = 1
        if k%n == 0:
            printpref = "& "
        else:
            printpref = ""

        printend = ""
        if k%16 !=0:
            if (kk != len(rsimp_ost)):                
                printend = "\\\\"
                
        #     else:
        #         #printend = ''
        # else:
        #     #printend = 'haaaa'
            
        if n > 1:
            cost = "{int_name}V^{vv}v^{v}o^{o} &\\quad N^{{{s}}}".format(int_name = int_name, vv=costidx[0], v=costidx[1], o=costidx[2], s = sum(costidx))
            over = "^{{{cost}}}".format(cost=cost)
            st = s+printpref+"\overbracket[0.05ex]{"+str(z)+"}"+over+printend
        else:
            inm = f"{int_name}"
            cs = "V^{vv}v^{v}o^{o} &\\quad N^{{{s}}}".format(vv=costidx[0], v=costidx[1], o=costidx[2], s = sum(costidx))
            if all_hash == []:
                st = s+ printpref+inm + str(z) + "&\\quad "+cs + printend
            else:
                st = s+printpref+ inm +"=" + str(z) + "&\\quad "+cs + printend 
        #print(st)
        lx += st
        # lx += f"{k}, {len(rsimp_ost)}"
        if (kk != len(rsimp_ost)):
            if (k%16 !=0):
                lx += "\n"
        k += 1
        kk += 1
    lx+= "\\end{alignat*}\n \\end{minipage}"
    lx+="\n"
    if f != []:
        f.write(lx)


def just_latex(f, n, rsimp_ost_in, name=False, all_hash = [], f12 = False, idx_str = '', interm_all_xfx =[]):

    print('JUST LATEX')
    print('innnn')
    for x in rsimp_ost_in :
        print(x)
    #print('')
    
    # grupuję, żeby były najpierw male sumowania, pozniej duze
    all_summation_list = []
    for x in rsimp_ost_in:
        x.summation.sort()
        if x.summation not in all_summation_list:
            all_summation_list.append(x.summation)
    all_summation_list = sorted(all_summation_list, key=lambda i: len(i))
    rsimp_ost_sort = []
    for s in all_summation_list:
        for x in rsimp_ost_in:
            if x.summation == s:
                rsimp_ost_sort.append(x)

    #print('soooor')
 #   for x in rsimp_ost_in :
        #print(x)
    #print('')
                
    
    # tutaj zamieniam indeksy i> na i1 dla latexa
    rsimp_ost = deepcopy(rsimp_ost_in)    
    if f12 == True:
        for k in range(0, len(rsimp_ost)):

            subst_list_old = []
            subst_list_new = []
            for i in rsimp_ost[k].summation:
                if ">" in i:
                    if i in virtual:
                        ii = free_idx(latex_virtual, subst_list_new)
                        subst_list_old.append(deepcopy(i))
                        subst_list_new.append(deepcopy(ii))
                    if i in occupied:
                        ii = free_idx(latex_occupied, subst_list_new)
                        subst_list_old.append(deepcopy(i))
                        subst_list_new.append(deepcopy(ii))
                # JESLI CHCESZ MIEC INDEKSY X_1, X2 w latexu to odkomentuj to
                # if i in CABS:
                #     ii = free_idx(latex_CABS, subst_list_new)
                #     subst_list_old.append(deepcopy(i))
                #     subst_list_new.append(deepcopy(ii))

            #print("subst1", k, rsimp_ost[k])
            #print(subst_list_old)
            #print(subst_list_new)
            for j in range(0, len(subst_list_old)):
                rsimp_ost[k].substitute(subst_list_old[j], subst_list_new[j])
            #print("subst2", rsimp_ost[k])


              

    int_name = ""
    
    k = 1
    kk = 1
    #    #print("\equag{")
    t = -1
    lx = "\equa{"
    for z in rsimp_ost:
        t += 1

        fx_ret = z.establish_fixed_and_return()
        # lx += f'teraz idzie k i kk {k}, {kk}'
        fx = ""
        if name:
            namestr = all_hash[z.binary_hash]
            #print('namestr', z, namestr)
            if len(interm_all_xfx) > 0:
                fx = interm_all_xfx[t]
                # fx = find_fixed_for_interm(z)
            
            fxstr = ""
            for x in fx:
                fxstr += x
            lenl = 0
            for s in namestr:
                if s == "l":
                    lenl += 1
            no = namestr[lenl:len(namestr)]
            if lenl == 1:
                namestr = f"(\Theta_{{{no}}})" 
            elif lenl == 2:
                namestr = f"(\Xi_{{{no}}})" 
            elif lenl == 3:
                namestr = f"(\\vartheta_{{{no}}})"
            elif lenl == 4:
                namestr = f"(\\upsilon_{{{no}}})" 
            else:
                #print('implement name for higher level intermediates')
                sys.exit(0)

            int_name = f"{namestr}_{{{fxstr}}}\qquad "
        memidx, costidx = memcost(z)
        # lx+= f'w tym momencie jest {kk},{len(rsimp_ost)} {kk%33} '
        s = ""
        if (k-1) == 16:
            if (kk == len(rsimp_ost)):
                s = """"
                }
                """
            else:
                if (kk-1)%32 == 0:#33:
                    s = """"
                    }
                    \\newpage
                    \equa{
                    """
                else:
                    s = """
                    }
                    \equa{
                    """
            k = 1

        if k%n == 0:
            printpref = "& "
        else:
            printpref = ""

        printend = ""
        if k%16 !=0:
            if (kk != len(rsimp_ost)):
                printend = "\\\\"
            
        if n > 1:
            cost = "{int_name}V^{vv}v^{v}o^{o} &\\quad N^{{{s}}}".format(int_name = int_name, vv=costidx[0], v=costidx[1], o=costidx[2], s = sum(costidx))
            over = "^{{{cost}}}".format(cost=cost)
            st = s+printpref+"\overbracket[0.05ex]{"+str(z)+"}"+over+printend
        else:
            inm = f"{int_name}"
            cs = "V^{vv}v^{v}o^{o} &\\quad N^{{{s}}}".format(vv=costidx[0], v=costidx[1], o=costidx[2], s = sum(costidx))
            if all_hash == []:
                st = s+ printpref+inm + str(z)+ printend
            else:
                st = s+printpref+ inm +"=" + str(z) + printend 
        #print(st)
        lx += st

        if (kk != len(rsimp_ost)):
            if (k%16 !=0):
                lx += "\n"
        k += 1
        kk += 1
    lx+= "}"
    lx+="\n"
    if f != []:
        f.write(lx)
        print('INTERM LATEX')
        print(lx)
        vspace(0)

        
def memcost(k):

#    #print('sprawdzam koszt kxk', k)
    memc = [0, 0, 0]
    cc = [0, 0, 0]

    memidx = []
    cidx = []

    for f in k.coefficient_idx:
        for ff in f:
            if ff not in memidx:
                if ff not in k.summation:
                    if ff in CABS or ff in fortran_CABS:
                        memc[0] += 1
                    elif ff in virtual or ff in fortran_virtual or ff in mona_virt or ff in monb_virt:
                        memc[1] += 1
                    elif ff in occupied or ff in fortran_occupied  or ff in mona_occup or ff in monb_occup:
                        memc[2] += 1
                    memidx.append(ff)
            if ff not in cidx:
#                #print('sprawdzam ff', ff)
                if ff in CABS or ff in fortran_CABS or ff in latex_CABS:
 #                   #print('cabs')
                    cc[0] += 1
                elif ff in virtual or ff in fortran_virtual or ff in latex_virtual or  ff in mona_virt or ff in monb_virt:
  #                  #print('virt')
                    cc[1] += 1
                elif ff in occupied or ff in fortran_occupied or ff in latex_occupied or ff in mona_occup or ff in monb_occup:
   #                 #print('occ')
                    cc[2] += 1
                cidx.append(ff)
   # #print('cc', cc)
    return memc, cc
