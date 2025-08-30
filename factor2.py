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
from templates import  *
from joblib import Parallel, delayed
import sys
import time
import pickle
import numpy as np
from joblib import Parallel, delayed


COST_THRESH = 7


dct_idx_cost = {'a': 10, 'b':10, 'c':10, 'd':10, 'e':10, 'f':10, 'g':10, 'i':1, 'j':1, 'k':1, 'l':1, 'm':1, 'n':1, \
                    'Q':100, 'X':10, 'Y':10, 'Z':10, 'Z>':10}


def idx_to_cost(e_idx):

    e_cost = []
    for i in e_idx:
        e_cost.append(dct_idx_cost[i])

    return e_cost

def find_all_diff_idx(e):

    e_idx = []

    for x in e.coefficient_idx:
        for y in x:
            if y not in e_idx:
                e_idx.append(y)

    e_idx.sort(key=str.casefold)

    return e_idx

def construct_idx_matrix(e, e_idx):
    
    len_idx = len(e_idx)
    len_coef = len(e.coefficient)


    idx_matrix = np.zeros((len_coef, len_idx), dtype=int)

#    print(len_idx, len_coef)
    for i in range(0, len_coef):
        for x in e.coefficient_idx[i]:
            place_idx = e_idx.index(x)
            idx_matrix[i][place_idx] = 1

#    print(idx_matrix)
    for i in range(0, len_idx):
        k = 0
        for j in range(0, len_coef):
            if idx_matrix[j][i] == 1:
                k+=1
            
    return idx_matrix, len_idx, len_coef

def pre_factorize2(e):
    e_idx  = find_all_diff_idx(e)
    e_cost = idx_to_cost(e_idx)
    idx_matrix, len_idx, len_coef = construct_idx_matrix(e, e_idx)
#    factorize2(idx_matrix, len_idx, len_coef, e_cost)
    k = 0

    idx_dict = {}
    for x in range(0, len_coef):
        idx_dict[x+1] = x+1

    minilevel = []
    biglevel = []
    pair_status = []
    big_idx_matrix = [deepcopy(idx_matrix)]
    big_idx_dict = [deepcopy(idx_dict)]

    pair_status_idx = 0
    idx_matrix_old = deepcopy(idx_matrix)
    idx_dict_old = deepcopy(idx_dict)
    factorize3(k, pair_status, pair_status_idx, biglevel, minilevel, idx_dict, idx_dict_old, idx_matrix, idx_matrix_old, len_idx, len_coef, e_cost, \
                   big_idx_matrix, big_idx_dict)

    # do ustalenia przed factorize 3

#     tm_list_big = []
#     tm_cost_big = []
#     tm_small = np.zeros((len_coef-1, len_coef-1), dtype=int)

#     factorize3(0, 0, tm_small, tm_list_big, tm_cost_big, pair_status_idx, pair_status,\
#                    idx_dict, idx_matrix, len_idx, new_len_coef, e_cost)

# def factorize3(ti, tj, tm_small, tm_list_big, tm_cost_big, pair_status_idx, pair_status, 
#                                idx_dict, idx_matrix, len_idx, new_len_coef, e_cost):


#     pair_list, pair_cost, pair_contr = find_pairs(idx_matrix, len_idx, len_coef, e_cost)

    
#     for y in range(0, len(pair_list)):
#         pair_status[pair_status_idx] = [y, len(pair_list)-1]
#         pair = pair_list[y] 
#         i = pair[0]-1
#         j = pair[1]-1

#         idx_matrix_temp = deepcopy(idx_matrix)

#         # Dodaj do siebie rzedy i oraz j zeby skontraktowac wyrazy                                                    
#         c = np.bitwise_xor(idx_matrix_temp[i][:], idx_matrix_temp[j][:])
#         # wpisz skontraktowane wyrazy na miejsce i                                                                                        
#         idx_matrix_temp[i][:] = deepcopy(c)
#         # nowy wymiar macierzy bedzie o 1 mniejszy                                                            
#         new_len_coef = len_coef - 1
#         print('na miejscu', i+1, 'jest para', i+1, j+1)
#         # usun wiersz j-ty                                                                                                                             
#         idx_matrix_temp = np.delete(idx_matrix_temp, j, 0)
#         # nadpisz slownik        

#         idx_dict_temp = {}
#         for m in range(0, j):
#             idx_dict_temp[m+1] = idx_dict[m+1]

#         for m in range(j, len_coef-1):
#             idx_dict_temp[m+1] = idx_dict[m+1] + 1


#         if i+1 in pair: # jesli w macierzy jest indeks ktory jest skontraktowany juz to nowy indeks bedzie dobudowany w prawo         

#             print('tak, i+1=', i+1, 'jest w perm', perm)
#             new_level_v+=1
#             # sprawdz na którym miejscu w perm jest indeks perm = [a, b]                                                                                 
#             if perm[0] == i+1:
#                 dict_matrix_temp[new_level_v][new_level_h] = new_perm[1]
#                 print('bede usuwac wyraz', perm[1])
#                 ii = perm[0]-1
#                 jj = perm[1]-1
#                 print('ktory jest tak naprawde', idx_dict[perm[1]], 'i zapisywac na miejsce', ii)
#                 print('')
#             else:
#                 dict_matrix_temp[new_level_v][new_level_h] = new_perm[0]
#                 print('bede usuwac wyraz', perm[0])
#                 ii = perm[1]-1
#                 jj = perm[0]-1
#                 print('ktory jest tak naprawde', idx_dict[perm[0]], 'i zapisywac na miejsce', ii)
#                         print('')
#         else:
#             # nie, to jest swieza para
#             tm_small[0][tj] = pair 


def factorize3(p, pair_status, pair_status_idx, biglevel, minilevel, idx_dict, idx_dict_old, idx_matrix, idx_matrix_old, len_idx, len_coef, e_cost, \
                   big_idx_matrix, big_idx_dict):
    
    print('na wejsciu do factorize', pair_status)
    print('len_coef', len_coef)
    print(idx_matrix)
    pair_list, pair_cost, pair_contr = find_pairs(idx_matrix, len_idx, len_coef, e_cost)
    print(pair_list)
    level_pairs = find_level_pairs(pair_list, len_coef)
    if (len(level_pairs)) == 0:
        print('len zero')
        sys.exit(0)

    print('level_pairs')
    kk = 0
    for x in level_pairs:
        print(kk, x)
        kk+=1
    print('')
    # pikapika
    #p = pair_status[pair_status_idx][0]

    print('zaczynam moja petle', p)
    print('pair status_idx', pair_status_idx)
    if (p == 0):
        pair_status.append([0, len(level_pairs)-1])
        big_idx_matrix.append(deepcopy(idx_matrix))
        big_idx_dict.append(deepcopy(idx_dict))
    print(pair_status)


    if p == len(level_pairs):
        if pair_stauts_idx == 0:
            print('koniec')
        else:
            print('cofnij sie o jeszcze jeden krok')


        for i in range(p, len(level_pairs)+1):
            print('i', i)
            print('jestem na zestawie,  i', i, 'z', len(level_pairs), level_pairs[i])
            pair_status[pair_status_idx] = [i+1, len(level_pairs)]
            print('pair popopopopopop', pair_status)
            print('idx_matrix')
            print(idx_matrix)
            print('idx matrx pair status')
            print(big_idx_matrix[pair_status_idx])
            minilevel.append(level_pairs[i])
            idx_matrix_temp = deepcopy(idx_matrix)

            print('i', i, level_pairs[i])
            del_list = []
            for enum in range(0, len(level_pairs[i])):
                elem = level_pairs[i][enum] 
                print('elem', elem)
                if len(elem) == 2:
                    print(elem)
                    ii = elem[0]-1
                    jj = elem[1]-1
                    c = np.bitwise_xor(idx_matrix_temp[ii][:], idx_matrix_temp[jj][:])
                    idx_matrix_temp[ii][:] = deepcopy(c)
                    del_list.append(jj)
                    idx_matrix_temp = np.delete(idx_matrix_temp, jj, 0)

                    new_len_coef = idx_matrix_temp.shape[0]

                    print('idx_dict', idx_dict)
                    print('lacze', elem)
                    idx_dict_temp = {}
                    for m in range(0, jj):
                        idx_dict_temp[m+1] = idx_dict[m+1]
     #                   idx_dict_temp[m+1] = big_idx_dict[pair_status_idx][m+1]
                    print(idx_dict_temp)
                    for m in range(jj, len_coef-1):
                        print('m+1', m+1, big_idx_dict[pair_status_idx][m+1] + 1)
                        idx_dict_temp[m+1] = idx_dict[m+2] 
    #                    idx_dict_temp[m+1] = big_idx_dict[pair_status_idx][m+2]
                    print('idx_dict_temp', idx_dict_temp)
                    print('idx_dict', idx_dict)
                    # inv_idx_dict_temp = {v: k for k, v in deepcopy(idx_dict_temp).items()}
                    # print('inv', inv_idx_dict_temp)
                    for enumrest in range(enum+1, len(level_pairs[i])):
                        rest_elem = level_pairs[i][enumrest]
                        if len(rest_elem) == 2:
                            print(rest_elem)
                            if rest_elem[0] >= jj:
                                level_pairs[i][enumrest][0] -=1
                            if rest_elem[1] >=jj:
                                level_pairs[i][enumrest][1] -=1


            if new_len_coef <= 2:
                print('pair', pair_status)
                print('LENNNN', new_len_coef)
                print(pair_status[-1])
                if pair_status[-1][0] ==pair_status[-1][1]:
                    biglevel.append(minilevel)
                    print('dodaje minilevel')
                    print('')
                    for ff in minilevel:
                        print(ff)
                    print('')
                    # koniec tej listy par, cofnij sie o jeden duzy
                    print(pair_status)
                    print('koniec tej listy par, cofnij sie o jeden duzy - wroc do poprzedniego factorize')

                    big_idx_matrix = big_idx_matrix[:-1]
                    big_idx_dict = big_idx_dict[:-1]
                    print(len(big_idx_matrix))
                    # for ff in big_idx_matrix:
                    #     print(ff)
                    #     print('')
                    print('big idx')
                    print(big_idx_matrix[-1])
                    print('oto co mam')
                    rr = 0
                    for r in biglevel:
                        print(rr, r)
                        rr +=1
                    minilevel = minilevel[:-2]
                    print('minilevel', minilevel)
                    pair_status = pair_status[:-1]
                    print('pair status', pair_status)
                    print('idx_dict', idx_dict)
                    print('idx_dict_temp', idx_dict_temp)

                    p = pair_status[pair_status_idx-1][0]
                    new_len_coef = big_idx_matrix[-1].shape[0]
    #                big_idx_matrix = big_idx_matrix[:-2]

                    return factorize3(p, pair_status, pair_status_idx-1,  biglevel, minilevel, big_idx_dict[-1], idx_dict, big_idx_matrix[-1], \
                                          idx_matrix, len_idx, new_len_coef, e_cost, big_idx_matrix, big_idx_dict)

                else:
                    # nowy element z tej listy par, cofnij sie o maly krok
                    biglevel.append(minilevel)
                    print('przed', minilevel)
                    minilevel = minilevel[:-1]
                    print('po   ', minilevel)

                    print('idx matrix temp')
                    print(idx_matrix_temp)
                    print('idx matrix')
                    print(idx_matrix)
                    print('nowy element z tej listy par, cofnij sie o maly krok')
                    print(idx_dict)
                    print(pair_status)
                    print(pair_status_idx)
                    print(i, len(level_pairs))

            else:
                print('new len coef')
                print('NIE COFAM SIE IDE DALEJ')
                print(pair_status_idx)
                print(' status', pair_status)
                print(idx_matrix_temp)
                print(idx_dict_temp)
                print('')
                print('obecne minilvel to', new_len_coef)
                print('minilevel')
                for g in minilevel:
                    print(g)
                pair_status_idx += 1
                print('wywoluje nowe factorize', pair_status)
    #            p = pair_status[pair_status_idx][0]
                return factorize3(0, pair_status, pair_status_idx,  biglevel, minilevel, idx_dict_temp, idx_dict, idx_matrix_temp, \
                               idx_matrix, len_idx, new_len_coef, e_cost, big_idx_matrix, big_idx_dict)

    


def find_level_pairs(pair_list, len_coef):


    level_pairs = []
    n = len_coef//2
    npair = len(pair_list)
    big_lp = []

    for i in range(1, n+1):
        lp = sublevel(pair_list, npair, i, len_coef)
        for x in lp:
            big_lp.append(x)
#        big_lp.append(lp)


    # r = 1
    # for x in big_lp:
    #     print(r, x)
    #     r+=1
        # k = 1
        # for y in x:
        #     print(r, k, y)
        #     k+=1
        # print('')
        # r+=1

    return big_lp
                        
def sublevel(pair_list, npair, n, len_coef):

    big_lp = []

    for k1 in range(0, npair):
        c = pair_list[k1][0]
        d = pair_list[k1][1]
        lp = [pair_list[k1]]
        idx_set = {c, d}
        if len(lp)<n:
            
            for k2 in range(k1+1, npair):
                a = pair_list[k2][0]
                b = pair_list[k2][1]
                if a not in idx_set and b not in idx_set:
                    lp.append(pair_list[k2])
                    idx_set.add(a)
                    idx_set.add(b)
                    if len(lp) == n:
                        lpcp = deepcopy(lp)
                        for sing_elem in range(1, len_coef+1):
                            if sing_elem not in idx_set:
                                lp.append([sing_elem])
                        big_lp.append(lp)
                        lp = lpcp[:-1]
                        idx_set.remove(a)
                        idx_set.remove(b)
        else:
            for sing_elem in range(1, len_coef+1):
                if sing_elem not in idx_set:
                    lp.append([sing_elem])
            big_lp.append(lp)


    return big_lp





def factorize2(idx_matrix, len_idx, len_coef, e_cost):

    perm_list, perm_cost, perm_contr = find_pairs(idx_matrix, len_idx, len_coef, e_cost)
    
    idx_dict = {}
    for x in range(0, len_coef):
        idx_dict[x+1] = x+1

    tm_list_big = []
    tm_cost_big = []

    print(perm_list)

    big_dict_matrix = []
    big_dict_idx = []

    perm_status = []
    perm_status.append([0, len(perm_list)-1])
    perm_status_idx = 0

    for y in range(0, len(perm_list)):
        tm_small = []
        perm_status[perm_status_idx] = [y, len(perm_list)-1]

        print('I KOLEJNA PERMUTACJA')
#        x = perm_list[y]
        x = [4,5]
        tm_small.append([x])
        # utworz macierz slownikow
        dict_matrix = np.zeros((2, len_coef-1), dtype=int)
        big_dict_matrix.append(dict_matrix)

        i = x[0]-1
        j = x[1]-1
        idx_matrix_temp = deepcopy(idx_matrix)

        # Dodaj do siebie rzedy i oraz j zeby skontraktowac wyrazy
        c = np.bitwise_xor(idx_matrix_temp[i][:], idx_matrix_temp[j][:])
        # wpisz skontraktowane wyrazy na miejsce i
        idx_matrix_temp[i][:] = deepcopy(c)
        # nowy wymiar macierzy bedzie o 1 mniejszy
        new_len_coef = len_coef - 1
        print('na miejscu', i+1, 'jest para', i+1, j+1)
        # usun wiersz j-ty
        idx_matrix_temp = np.delete(idx_matrix_temp, j, 0)
        # nadpisz slownik
        
        
        idx_dict_temp = {}
        for m in range(0, j):
            idx_dict_temp[m+1] = idx_dict[m+1]

        for m in range(j, len_coef-1):
            idx_dict_temp[m+1] = idx_dict[m+1] + 1

        big_dict_idx.append(idx_dict_temp)
        
        level_v = 0
        level_h = 0

        big_dict_matrix[y][level_v][level_h] = i+1
        big_dict_matrix[y][level_v+1][level_h] = j+1

        print('')
        print(dict_matrix)

        
        print('wywoluje update')
        print('')
        perm_status.append([0, 0])
        perm_status_idx += 1
        update_dict_matrix(tm_small, perm_status_idx, perm_status, i, level_v, level_h, \
                               big_dict_matrix[y], idx_dict_temp, idx_matrix_temp, len_idx, new_len_coef, e_cost)


def update_dict_matrix(perm_status_idx, perm_status, i, level_v, level_h, dict_matrix, idx_dict, idx_matrix, len_idx, len_coef, e_cost):
    
    # wyszukaj pary w nowym  ciągu
    perm_list, perm_cost, perm_contr = find_pairs(idx_matrix, len_idx, len_coef, e_cost)
    # przetlumacz na nowe
    level_h+=1
    print('perm list', perm_list)
    print('slownik wierszy', idx_dict)
    print('')
    print('ide po permutacjach')
    for y in range(0, len(perm_list)):
        perm_status[perm_status_idx] = [y, len(perm_list)-1]
        print('perm_status', perm_status)
        print('')
        print('tu bedzie tyle permutacji do sprawdzenia', len(perm_list))
        print('')
        perm = perm_list[y]
        print('stary lewel v jest rowny', level_v)
        print('stary lewel h jest rowny', level_h)
        print('biore permutacje', y, 'rowna', perm)
        new_level_v = deepcopy(level_v)
        new_level_h = deepcopy(level_h)
        dict_matrix_temp = deepcopy(dict_matrix)
        print('dict')
        print(dict_matrix)
        # musze przetlumaczyc nowe m na przesuniete indeksy                                                                                              
        new_perm = deepcopy(perm)
        new_perm[0] = idx_dict[perm[0]]
        new_perm[1] = idx_dict[perm[1]]
        print('jakie bylo i', i+1)

        if i+1 in perm: # jesli w macierzy jest indeks ktory jest skontraktowany juz to nowy indeks bedzie dobudowany w prawo
            print('tak, i+1=', i+1, 'jest w perm', perm)
            new_level_v+=1
            # sprawdz na którym miejscu w perm jest indeks perm = [a, b]
            if perm[0] == i+1:
                dict_matrix_temp[new_level_v][new_level_h] = new_perm[1]
                print('bede usuwac wyraz', perm[1])
                ii = perm[0]-1
                jj = perm[1]-1
                print('ktory jest tak naprawde', idx_dict[perm[1]], 'i zapisywac na miejsce', ii)
                print('')
            else:
                dict_matrix_temp[new_level_v][new_level_h] = new_perm[0]
                print('bede usuwac wyraz', perm[0])
                ii = perm[1]-1
                jj = perm[0]-1
                print('ktory jest tak naprawde', idx_dict[perm[0]], 'i zapisywac na miejsce', ii)
                print('')
        else:
            print('nie, i+1=', i+1, 'jnie jest w perm', perm)
            dict_matrix_temp[new_level_v][new_level_h] = new_perm[0]
            dict_matrix_temp[new_level_v+1][new_level_h] = new_perm[1]
            ii = new_perm[0]-1
            jj = new_perm[1]-1
            print('bede usuwac wyraz', perm[1])
            print('ktory jest tak naprawde', idx_dict[perm[1]], 'i zapisywac na miejsce', ii)
            print('')

        print('dcmatrix', perm)
        print(dict_matrix_temp)
        if len_coef >=3:
        
            print('stara macierz')
            print(idx_matrix)
            print('')

            # dodaj kolejne elelementy ciagu, zmniejsz macierz
            idx_matrix_temp = deepcopy(idx_matrix)
            c = np.bitwise_xor(idx_matrix_temp[ii][:], idx_matrix_temp[jj][:])
            idx_matrix_temp[ii][:] = deepcopy(c)
            new_len_coef = len_coef - 1
            print('na miejscu', ii+1, 'jest para', ii+1, jj+1)
            idx_matrix_temp = np.delete(idx_matrix_temp, jj, 0)

            print('nowa macierz', new_len_coef)
            print(idx_matrix_temp)
            print('')
            
            # zmien nazwy w slowniku
            idx_dict_temp = {}
            for m in range(0, jj):
                idx_dict_temp[m+1] = idx_dict[m+1]

            for m in range(jj, new_len_coef):
                idx_dict_temp[m+1] = idx_dict[m+1+1]
    
            print('idx_dict_temp', idx_dict_temp)


            if i == ii:
                new_level_v -=1

            
            # sprawdz czy teraz nie doszlismy do konca, jesli nie, to wywolaj update dict matrix
            if new_len_coef >= 3:
                print('jeszcze moge dokladac cegielki')
                perm_status.append([0,0])
                perm_status_idx += 1
                update_dict_matrix(perm_status_idx, perm_status, ii, new_level_v, \
                                       new_level_h, dict_matrix_temp, idx_dict_temp, idx_matrix_temp, len_idx, new_len_coef, e_cost)
            # doszlismy do konca tej galezi, zobacz czy sa jeszcze permutacje
            else:
                # zobacz czy mamy jeszcze permutacje z tej galezi
                print('koniec tej galezi, cofnij sie o krok', y, len(perm_list)-1)
                if y == len(perm_list)-1:
                    print('koniec listy permutacji')
                    print(new_level_v, new_level_h)
                    new_level_h -=1
                    big.append(dict_matrix_temp)
                    print('koniec listy permutacji')

                    update_dict_matrix(ii, new_level_v, new_level_h, dict_matrix_temp, idx_dict_temp, \
                                           idx_matrix_temp, len_idx, new_len_coef, e_cost)
                else:
                    
                    print('sa jeszcze permutacje, kontynuuj petle')
                
        else:
            if perm_status[perm_status_idx][0] == perm_status[perm_status_idx][1]:
                print('koniec tych permutacji')
            else:
                print('jeszcze do sprawdzenia inne permutacje')
            print('to koniec nie zostalo nic')
            sys.exit(0)


def factorize(e):
    e_idx  = find_all_diff_idx(e)
    e_cost = idx_to_cost(e_idx)
    idx_matrix, len_idx, len_coef = construct_idx_matrix(e, e_idx)
    

    perm_list, perm_cost, perm_contr = find_pairs(idx_matrix, len_idx, len_coef, e_cost)

    for x in perm_list:
        print(x)
#    sys.exit(0)

    for i in range(0, len_coef-2):
        print('i', i, 'rekur')
        perm_list, perm_cost, perm_contr = recur_idx_matrix(idx_matrix, perm_list, perm_cost, perm_contr, len_idx, len_coef, e_cost)
    sys.exit(0)

    big_list_all = []
    big_list_sum = []
    for x in range(0, len(perm_list)):
        big_list_all.append([perm_list[x], perm_cost[x], perm_contr[x]])
        cst = 0
        for y in perm_cost[x]:
            cst = cst + 28**y
        big_list_sum.append(cst)

    if len(perm_list) != 0:
        list1, list2 = (list(x) for x in zip(*sorted(zip(big_list_sum, big_list_all))))
        print('POSORTOWANE i tylko najtansze')
        v = np.array(list1)
        s = v[0]
        vv = v%s
        vvv = v/s
        ii = 0
        for i in range(0, len(v)):
            if (vv[i]==0 and vvv[i]==1):
                ii+= 1
        list1 = list1[0:ii]
        list2 = list2[0:ii]
        k = 1
        for x in range(0, len(list1)):
            print(k, list1[x], sum(list2[x][2]), list2[x])
            k+=1
        print(len(perm_list))
    else:
        print('NIE MA INTERMEDIATE!')
        list1 = []
        list2 = []
        print(e)
    

    return list1, list2


def rearrange_string(e, list2):

    e_out = deepcopy(e)
    e_out.coefficient = []
    e_out.coefficient_idx = []


    for j in reversed(list2):
        e_out.coefficient.append(e.coefficient[j-1])
        e_out.coefficient_idx.append(e.coefficient_idx[j-1])

    return e_out
        
    


def real_contractions(i, j, n, k, idx_matrix, c):

    add = True
    vec = []
    for l in range(0, n):
        vec.append(0)

    for l in range(0, k):
        if (l != i and l!= j):
            b = np.bitwise_and(c, idx_matrix[l][:])
            if sum(b) > 0:
                vec = deepcopy(np.bitwise_or(vec, b))
    if np.array_equal(vec, c):
        add = False

    return add

def find_pairs(idx_matrix, n, k, e_cost):
    
    perm_list = []
    perm_cost = []
    perm_contr = []
    for i in range(0, k):
        for j in range(i+1, k):
            if i != j:
                a = idx_matrix[i][:]
                b = idx_matrix[j][:]
                c = sum(np.bitwise_and(a,b))
                if c > 0:
                    d = sum(np.bitwise_or(a, b))
                    if d <= COST_THRESH:
                        # sprawdz czy indeksy sumacyjne nie wystepuja w zadnych innych wyrazach
                        cc = np.bitwise_and(a, b)
                        add = real_contractions(i, j, n, k, idx_matrix, cc)
                        if add == True:
                            dd = np.bitwise_xor(a,b)
                            perm_list.append([i+1,j+1])
                            perm_cost.append([d])
                            mem = 1
                            for h in range(0, len(dd)):
                                if dd[h] == 1:
                                    mem *= e_cost[h]
                            perm_contr.append([mem])

    return perm_list, perm_cost, perm_contr

def find_pairs_restricted(idx_matrix, i, n, k, pl, pcst, pcntr, perm_list_bigger, perm_cost_bigger, perm_contr_bigger, e_cost):

    perm_list = []
    perm_cost = []
    perm_contr = []
    
    a = idx_matrix[i][:]
    for j in range(0, k):
        pl_temp = deepcopy(pl)
        pcst_temp = deepcopy(pcst)
        pcntr_temp = deepcopy(pcntr)
        if j != i:            
            b = idx_matrix[j][:]
            c = sum(np.bitwise_and(a,b))
            if c > 0:
                d = sum(np.bitwise_or(a, b))
                if d <= COST_THRESH:
                    cc = np.bitwise_and(a, b)
                    add = real_contractions(i, j, n, k, idx_matrix, cc)
                    if add == True:
                        dd = np.bitwise_xor(a,b)
                        pl_temp.append(j+1)
                        pcst_temp.append(d)
                        mem = 1
                        for h in range(0, len(dd)):
                            if dd[h] == 1:
                                mem *= e_cost[h]

                        pcntr_temp.append(mem)
                        perm_list_bigger.append(pl_temp)
                        perm_cost_bigger.append(pcst_temp)
                        perm_contr_bigger.append(pcntr_temp)
                    
                    

def recur_idx_matrix(idx_matrix, perm_list, perm_cost, perm_contr, len_idx, len_coef, e_cost):
        
    perm_list_bigger = []
    perm_cost_bigger = []
    perm_contr_bigger = []
    # print('aaaaa')
    # print(perm_list)
    # print(perm_cost)

    for y in range(0, len(perm_list)):
        x = perm_list[y]
        print(x)
        i = x[0]-1

        idx_matrix_temp = deepcopy(idx_matrix)
        for k in range(1, len(x)):
            j = x[k]-1
            c = np.bitwise_xor(idx_matrix_temp[i][:], idx_matrix_temp[j][:])
            idx_matrix_temp[i][:] = deepcopy(c)
            idx_matrix_temp[j][:] = 0
    
        find_pairs_restricted(idx_matrix_temp, i, len_idx, len_coef, x, perm_cost[y], perm_contr[y], perm_list_bigger, \
                                  perm_cost_bigger, perm_contr_bigger, e_cost)

    return perm_list_bigger, perm_cost_bigger, perm_contr_bigger


############################################################


def find_permutations(pairs_list):#, idx_matrix, n, k):
    print(pairs_list, 'pairs_list')
    print('')
    permutation_list = []

    for i in range(0, len(pairs_list)):
        pairs_list[i].append(0)
        print(pairs_list[i])

    for i in range(0, len(pairs_list)):
        print('teraz sprawdzam gdy na 1 miejscu jest', pairs_list[i][0], )
        minip = []
        minip = deepcopy(pairs_list[i][0])
        print('minip - wchodze do basic permute', minip)
        basic_permute(i, 0, pairs_list, minip, permutation_list)
        for ll in range(0, len(pairs_list)):
            pairs_list[ll][2] = 0
        print(permutation_list, 'lissst')
        print('-----------------------------------')

    return permutation_list

def basic_permute(i, k, pairs_list, minip, permutation_list):

    print('zaczynam basic permute dla', pairs_list[i][0])
    print(pairs_list)
    lnp = len(pairs_list[i][1])
    print(pairs_list[i][2] , lnp)
    if (pairs_list[i][2] < lnp):        
        for j in range(pairs_list[i][2], lnp):
            jj = pairs_list[i][1][j] # pobieram element z listy sasiadujacych indeksow
            print(pairs_list[i][0],  'element', j, 'rowny', pairs_list[i][1][j], 'minip', minip)
            pairs_list[i][2] += 1
            if jj not in minip:
                print('dodaje')
                minip.append(jj)
                print('minip', minip, pairs_list[i])

                if len(minip) < len(pairs_list):
                    print('mniejsze wiec ide dalej')
                    return basic_permute(jj-1, 0, pairs_list, minip, permutation_list)
                elif len(minip) == len(pairs_list):
                    print('rowne')
                    permutation_list.append(minip)
                    print('wracam do poprzedniego kroku', j, lnp-1)
                    print('================================================')
                    print('')
                    print('================================================')
                    if (j==lnp-1):
                        minip_old = deepcopy(minip)
                        minip = minip_old[0:len(minip)-2]
                        minipdalsze = minip_old[len(minip_old)-2:len(minip_old)]
                        print('minipdalsze', minipdalsze)
                        for ll in minipdalsze:
                            pairs_list[ll-1][2] = 0
                        print('to bedzie nowe minip')
                        nowe_i = minip[len(minip)-2]
                        #print(pairs_list[nowe_i])
                        print(minip, minip[len(minip)-1]-1)
                        return basic_permute(minip[len(minip)-1]-1, 0, pairs_list, minip, permutation_list)
                    else:
                        minip_old = deepcopy(minip)
                        minip = minip_old[0:len(minip_old)-1]
                        minipdalsze = minip_old[len(minip_old)-1:len(minip_old)]
                        for ll in minipdalsze:
                            pairs_list[ll-1][2]= 0
            else:
                print('juz jest', minip)
                if (j==lnp-1):
                    print('laaaaaaaaaaaaa')
                    minip_old = deepcopy(minip)
                    minip = minip_old[0:len(minip_old)-1]
                    print('minip', minip)
                    minipdalsze = minip_old[len(minip_old)-1:len(minip_old)]
                    for ll in minipdalsze:
                        pairs_list[ll-1][2]= 0

                    print('to bedzieNOWE minip', len(minip)-1)
                    nowe_i = minip[len(minip)-1]
                    print('nowe_i', nowe_i)
                    #print(pairs_list[nowe_i])
                    print('minip')
                    print(minip, minip[len(minip)-1])
                    return basic_permute(minip[len(minip)-1]-1, 0, pairs_list, minip, permutation_list)

    else:
        print('dupa')
        minip_old = deepcopy(minip)
        minip = minip_old[0:len(minip_old)-1]
        if len(minip) > 0:
            minipdalsze = minip_old[len(minip_old)-1:len(minip_old)]
            for ll in minipdalsze:
                pairs_list[ll-1][2]= 0

            print('to bedzieNOWE minip', len(minip)-1)
            nowe_i = minip[len(minip)-1]
            print('nowe_i', nowe_i)
            print(pairs_list[nowe_i])
            print(minip, minip[len(minip)-1])
            return basic_permute(minip[len(minip)-1]-1, 0, pairs_list, minip, permutation_list)


        
