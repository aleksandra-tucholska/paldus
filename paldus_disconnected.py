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
import cProfile

import os, sys

def identify_disconnected(all_hash, basket_super, basket_outer_all, sorted_interms, interm_dict, interm_hash, \
                          disconnected_super, basket_original_only_with_interm):

    # basket_super_nodisc - basket super without disconnected terms
    # basket_outer_nodisc - basket outer without disconnected terms
    # basket_super_disconnected - basket super but only with disconnected terms
    # basket_outer_disconneted - basket outer but only with disconnected terms
    # disconnected_super_new - disconnected super but only for basket_super_disconnected
    # basket_original_disc - basket original but only for basket_super_disconnected

    basket_super_nodisc = []
    basket_outer_nodisc = []
    basket_super_disc = []
    basket_outer_disc = []
    disc_super_new = []
    basket_original_disc = []
    basket_original_nodisc = []

    basket_point_onlydisc = []
    disconnected_interms = {}
    nodisconnected_interms = {}

    n_of_disc = 0
    rest = 0
    print('huara')
    for x in range(0, len(basket_super)):
        print('basket nr', x, 'o dlugosci', len(basket_super[x]))
        print('disconnected dlan', disconnected_super[x])
        print('basket_original', basket_original_only_with_interm[x])
        disc_sum =0
        for z in disconnected_super[x]:
            disc_sum += sum(z)
        if disc_sum>0:
            print(x, 'do disc')
            n_of_disc += 1
            basket_super_disc.append(basket_super[x])
            basket_outer_disc.append(basket_outer_all[x])
            disc_super_new.append(disconnected_super[x])
            basket_original_disc.append(basket_original_only_with_interm[x])
        else:
            print(x, 'do nodisc')
            rest += 1
            basket_super_nodisc.append(basket_super[x])
            basket_outer_nodisc.append(basket_outer_all[x])
            basket_original_nodisc.append(basket_original_only_with_interm[x])

    print(len(basket_super_disc))
#    sys.exit(0)
    vspace(0)
    print('nadaje punkty disconnected')
    vspace(0)
    basket_points_super_disc = []
    basket_points_only_disc = []
    for x in range(0, len(basket_super_disc)):
        basket_points = []
        basket_points_disc = []

        for y in range(0, len(basket_super_disc[x])):
            points = 0

            for z in range(0, len(basket_super_disc[x][y])):
                
                if isdisc(interm_hash[basket_super_disc[x][y][z].binary_hash]):
                    basket_points_disc.append(interm_dict[basket_super_disc[x][y][z].binary_hash])
                    disconnected_interms[basket_super_disc[x][y][z].binary_hash] = interm_dict[basket_super_disc[x][y][z].binary_hash]
                points += interm_dict[basket_super_disc[x][y][z].binary_hash]
            basket_points.append(points)

        basket_points_super_disc.append(basket_points)

    print('len(basket_original_disc)', len(basket_original_disc))
    print('len(basket_original_nodisc)', len(basket_original_nodisc))

    sorted_disc_interms = sorted(disconnected_interms.items(), key=lambda x: x[1], reverse=True)

    print('sorteed disc')
    for key, value in sorted_disc_interms:
        print(value, all_hash[key], interm_hash[key], '                     ', isdisc(interm_hash[key]))
    print('')


    print('a to conncected')
    print(len(basket_super_nodisc))
    vspace(0)
    basket_points_super_nodisc = []
    basket_points_only_nodisc = []
    for x in range(0, len(basket_super_nodisc)):

        basket_points = []

        for y in range(0, len(basket_super_nodisc[x])):
            points = 0

            for z in range(0, len(basket_super_nodisc[x][y])):
                
                if isdisc(interm_hash[basket_super_nodisc[x][y][z].binary_hash]):
                    'WTF DISCONNECTED IN CONNECTED LIST?'
                else:
                    nodisconnected_interms[basket_super_nodisc[x][y][z].binary_hash] = interm_dict[basket_super_nodisc[x][y][z].binary_hash]
                points += interm_dict[basket_super_nodisc[x][y][z].binary_hash]
            basket_points.append(points)

        basket_points_super_nodisc.append(basket_points)

    print('len(basket_original_disc)', len(basket_original_disc))
    print('len(basket_original_nodisc)', len(basket_original_nodisc))

    sorted_nodisc_interms = sorted(nodisconnected_interms.items(), key=lambda x: x[1], reverse=True)

    print('sorteed nodisc')
    for key, value in sorted_nodisc_interms:
        print(value, all_hash[key], interm_hash[key], '                     ', isdisc(interm_hash[key]))
    print('')


    
    return sorted_disc_interms, basket_super_disc, basket_outer_disc, basket_original_disc, \
        sorted_nodisc_interms, basket_super_nodisc, basket_outer_nodisc, basket_original_nodisc, disc_super_new, basket_points_super_disc, \
        basket_points_super_nodisc
        

def pick_best_basket(sorted_interms, basket_super, basket_outer_all, basket_original, basket_points_super, interm_dict,  \
                     interm_hash, disconnected_super, interm_whichlevel_dict, all_hash, old_new_dict, n_interm, disconnected):


    print('poczatek pick_best_basekt')
    vspace(0)
    print()

    print('na wejsciu')
    k = 0
    for x in basket_original:
        print(k, x)
        k+= 1
    print('')

    outer_terms = arithmetic_string()

    # sorted interms
    # all intermediates sorted in descending orders according
    # to number of occuerence

    # all_hash_key
    # dictionary, takes binary_tag from ugg and gives interm name (e.g. inter23)
    #

    # basket_super
    # a list of baskets, one for each basket
    
    basket_super_used_idx = []

    # basket points zawiera punkty dla kazego zestawu
    # | 3 52   6   8   |, | 2 2  3 4|, |3 10  8 8 2  100|...
    #
    #
    start_time = time.time()
    print('teraz wybieranie intermsow')
    picked_interms = []
    picked_interms_dict = {}
    not_picked_interms = []
    numl = []



    picked_interms = set()
    picked_interms_dict = {}
    not_picked_interms = set()

    
    # Iterate through a sorted list of interms and check if the current most common interm key is present in the baskets of words
    for key, value in sorted_interms:        
        current_most_common_interm_key = all_hash[key]
        #        print(' current_most_common_interm_key', key,  current_most_common_interm_key, interm_hash[key])
        if (len(basket_super) - len(basket_super_used_idx)) ==0:
            break
        basket_super_not_used_idx = [x for x in range(0, len(basket_super)) if x not in basket_super_used_idx]

        # Iterate through baskets of words that contain the current interm
        minipicked = 0
        zostaje = 0
        for x in range(0, len(basket_super)):
            if x not in basket_super_used_idx:
                #                print('basket nr', x, 'o dlugosci', len(basket_super[x]))

                # Iterate through sets of words within the current basket
                is_present =  [0] * len(basket_super[x])
                for y in range(0, len(basket_super[x])):

                    # Iterate through intermediates in the current set of words
                    for z in range(0, len(basket_super[x][y])):
                        this_interm = basket_super[x][y][z]
                        this_interm_key =  basket_super[x][y][z].binary_hash
                        #
                        # check if this interm is equal to the current most common
                        #
                        if key == this_interm_key:
                            is_present[y] = 1
                            break

                # w is_present mam jedynki tam gdzie w ogole ten interm jest
                if sum(is_present) > 0:
                    if all_hash[key] in not_picked_interms:
                        picked_interms.add(all_hash[key])
                        picked_interms_dict[all_hash[key]] = deepcopy(interm_hash[key])
                        not_picked_interms.discard(all_hash[key])
                    basket_super_used_idx.append(x)
                    best_minibasket_idx = pick_best_minibasket(is_present, basket_points_super[x])
                    outer_terms.append(basket_outer_all[x][best_minibasket_idx])
                    numl.append(x)
                    # Iterate through the , intermediates in the selected minibasket and add them to the list of picked interms               
                    for z in range(0, len(basket_super[x][best_minibasket_idx])):
                        if all_hash[basket_super[x][best_minibasket_idx][z].binary_hash] not in picked_interms:
                            if all_hash[basket_super[x][best_minibasket_idx][z].binary_hash] in not_picked_interms:
                                not_picked_interms.discard(all_hash[basket_super[x][best_minibasket_idx][z].binary_hash])
                                picked_interms.add(all_hash[basket_super[x][best_minibasket_idx][z].binary_hash])
                                picked_interms_dict[all_hash[basket_super[x][best_minibasket_idx][z].binary_hash]] =\
                                    deepcopy(interm_hash[basket_super[x][best_minibasket_idx][z].binary_hash])
                            else:
                                picked_interms.add(all_hash[basket_super[x][best_minibasket_idx][z].binary_hash])
                                picked_interms_dict[all_hash[basket_super[x][best_minibasket_idx][z].binary_hash]] =\
                                    deepcopy(interm_hash[basket_super[x][best_minibasket_idx][z].binary_hash])
                                minipicked += 1

                else:
                    zostaje += 1

                    for y in range(0, len(basket_super[x])):
                        for z in range(0, len(basket_super[x][y])):
                            if all_hash[basket_super[x][y][z].binary_hash] not in not_picked_interms:
                                if all_hash[basket_super[x][y][z].binary_hash] in picked_interms:
                                    continue
                                else:
                                    not_picked_interms.add(all_hash[basket_super[x][y][z].binary_hash])


    k = 1

    
    print('finally', len(not_picked_interms), len(set(not_picked_interms)), len(picked_interms_dict))
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"DISC Elapsed time1: {elapsed_time} seconds")
#    sys.exit(0)
    start_time = time.time()



    for x in not_picked_interms:
        for y in picked_interms:
            if x == y:
                print(x, y)
                sys.exit(0)

    # interm_hash: klucz-binarny - intermediate
    # all_hash:    klucz-binarny - nazwa
    print('to sa nazwy intermediates', len(picked_interms_dict))
    for key in picked_interms_dict:
        print(key, picked_interms_dict[key])

    print('')
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"DISC Elapsed time2: {elapsed_time} seconds")
    start_time = time.time()


    print('na wyjsciu')
    k = 0
    for x in outer_terms:
        print(k, numl[k], x)
        k+=1 
    print('')


    print('na wyjsciu')
    k = 0
    for x in outer_terms:
        print(k, x)
        k+=1 
    print('')

    if disconnected:

        princp = deepcopy(picked_interms_dict)
        print('INTERMEDIATES I przed zmianami')
        k = 1
        for key in princp:
            num = int(''.join(filter(str.isdigit, key)))
            fx = find_fixed_for_interm2(princp[key])
            temp = ""
            for y in range(0, len(fx)):
                temp = temp + str(Tdict[fx[y]]) + str("")
            if len(fx)>0:
                temp = "_{"+temp+"}"
            princp[key].summation = []        
            if k%3 == 0:
                print('$', 'I^{{{}}}{temp}'.format(num, temp=temp),'$', '=', '$',princp[key],'$', '\\\\')
            else:
                print('$','I^{{{}}}{temp}'.format(num, temp=temp),'$', '=', '$',princp[key], '$', '&')
            k+=1
        print('')
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"DISC Elapsed time3: {elapsed_time} seconds")
        start_time = time.time()


        
        old_new_dict = {}
        outer_terms,  picked_interms_dict, old_new_dict, freq_dict = change_names(outer_terms, picked_interms_dict, old_new_dict, all_hash, interm_hash, True)
        for x in old_new_dict:
            print('stary', x, 'staje sie', old_new_dict[x])
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"DISC Elapsed time3: {elapsed_time} seconds")
        start_time = time.time()


        princp = deepcopy(picked_interms_dict)
        print('INTERMEDIATES I')
        k = 1
        for key in princp:
            num = int(''.join(filter(str.isdigit, key)))
            fx = find_fixed_for_interm2(princp[key])
            mem, cc = memcost(deepcopy(princp[key]))
            memstr = f"v^{{{mem[1]}}}o^{{{mem[2]}}}"
            ccstr = f"v^{{{cc[1]}}}o^{{{cc[2]}}}"
            temp = ""
            for y in range(0, len(fx)):
                temp = temp + str(Tdict[fx[y]]) + str("")
            if len(fx)>0:
                temp = "_{"+temp+"}"
            princp[key].summation = []        
            if k%3 == 0:
                print(key, '$', 'I^{{{}}}{temp}'.format(num, temp=temp),'$', '=', '$',princp[key],'$', 'mem=', memstr, 'cost=', ccstr,'\\\\')
            else:
                print(key, '$','I^{{{}}}{temp}'.format(num, temp=temp),'$', '=', '$',princp[key], '$', 'mem=', memstr, 'cost=', ccstr,'&')
            k+=1
        print('')
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"DISC Elapsed time3: {elapsed_time} seconds")
        start_time = time.time()

        no  =1
        disc_one_dict = {}
        for x in outer_terms:
            disc_one_dict, no = isdisc_one(x, disc_one_dict, no)

        print('na wyjsciu ost')
        k = 0
        for x in outer_terms:
            print(k, x)
            k+=1 
        print('')

    
        # rename all self disconnected so there are only unique ones
        # also rename the outer_terms
        outer_terms, disc_one_dict = identify_unique_self_disc(outer_terms, disc_one_dict)

        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"DISC Elapsed time4: {elapsed_time} seconds")
        start_time = time.time()


        # print('and original')
        # for x in range(0, len(outer_terms)):
        #     basket_original[numl[x]].summation = []
        #     print(basket_original[numl[x]])
        # print('')

        print('and last')
        for x in range(0, len(outer_terms)):
            print(outer_terms[x])
        print('')


        print('and last latex disconnected')
        for x in range(0, len(outer_terms)):
            orsum = deepcopy(len(basket_original[numl[x]].summation))
            outsum = deepcopy(len(outer_terms[x].summation))
            outer_terms[x].summation = []
            basket_original[numl[x]].summation = []
            print('$', basket_original[numl[x]],'$' , '&', '$', outer_terms[x],'$',  'cost przed:', orsum, 'cost po:', outsum,'\\\\')
        print('')


        n_interm = len(picked_interms_dict)

        
        print_interms(picked_interms_dict, outer_terms)


    else:
        print('conncected')
        for x in old_new_dict:
            print(x, old_new_dict[x])

        all_hash_swapped = {v: k for k, v in all_hash.items()}
            
        for x in old_new_dict:
            for key in picked_interms_dict:

                if all_hash_swapped[x] == picked_interms_dict[key].binary_hash:
                    print('old', x, old_new_dict[x])
                    print('picked', key, picked_interms_dict[key])
                    print('')
                    
        princp = deepcopy(picked_interms_dict)
        print('INTERMEDIATES przed zmianami nazw')
        k = 1
        for key in princp:
            num = int(''.join(filter(str.isdigit, key)))
            fx = find_fixed_for_interm2(princp[key])
            temp = ""
            for y in range(0, len(fx)):
                temp = temp + str(Tdict[fx[y]]) + str("")
            if len(fx)>0:
                temp = "_{"+temp+"}"
            princp[key].summation = []        
            if k%3 == 0:
                print('$', 'I^{{{}}}{temp}'.format(num, temp=temp),'$', '=', '$',princp[key],'$', '\\\\')
            else:
                print('$','I^{{{}}}{temp}'.format(num, temp=temp),'$', '=', '$',princp[key], '$', '&')
            k+=1
        print('')

        print('chchchange')
        outer_terms,  picked_interms_dict, old_new_dict,freq_dict = change_names(outer_terms, picked_interms_dict, old_new_dict, all_hash, interm_hash, False)
        
        princp = deepcopy(picked_interms_dict)
        print('INTERMEDIATES I')
            
        k = 1
        for key in princp:            
            num = int(''.join(filter(str.isdigit, key)))
            fx = find_fixed_for_interm2(princp[key])
            mem, cc = memcost(princp[key])
            memstr = f"v^{{{mem[1]}}}o^{{{mem[2]}}}"
            ccstr = f"v^{{{cc[1]}}}o^{{{cc[2]}}}"

            temp = ""
            for y in range(0, len(fx)):
                temp = temp + str(Tdict[fx[y]]) + str("")
            if len(fx)>0:
                temp = "_{"+temp+"}"
            princp[key].summation = []        
            if k%3 == 0:
                print('$', 'I^{{{}}}{temp}'.format(num, temp=temp),'$', '=', '$',princp[key],'$', 'mem=', memstr, 'cost=', ccstr, '\\\\')
            else:
                print('$','I^{{{}}}{temp}'.format(num, temp=temp),'$', '=', '$',princp[key], '$', 'mem=', memstr, 'cost=', ccstr,'&')
            k+=1
        print('')
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"DISC Elapsed time3: {elapsed_time} seconds")
        start_time = time.time()

        # print_interms(picked_interms_dict, outer_terms)

        
        print('and last latex connected')
        for x in range(0, len(outer_terms)):
            orsum = deepcopy(len(basket_original[numl[x]].summation))
            outsum = deepcopy(len(outer_terms[x].summation))
            outer_terms[x].summation = []
            basket_original[numl[x]].summation = []
            print('$', basket_original[numl[x]],'$' , '&', '$', outer_terms[x],'$',  'cost przed:', orsum, 'cost po:', outsum,'\\\\')
            for xx in outer_terms[x].coefficient:
                if 'interm' in xx:
                    print('xx',xx, picked_interms_dict[xx])
        print('')
        

    return old_new_dict, n_interm



def isdisc(res):

    for x in res.coefficient_idx:
        for y in x:
            if y not in res.summation:
                return False

    return True

def isdisc_one(res, disc_one_dict, no):

    
    for i in range(0, len(res.coefficient_idx)):
        
        disc_idx = res.coefficient_idx[i]
        all_other_idx = []
        for j in range(0, len(res.coefficient_idx)):
            if j!=i:
                all_other_idx += res.coefficient_idx[j]
        

        not_in_sum = set(disc_idx).difference(set(res.summation))
        if not_in_sum:
            print('are there some fixed indices?')
        else:
            not_in_others = set(disc_idx).difference(set(all_other_idx))
            if not_in_others:
                # print('this coef is self disconnected')
                # print(i, res, res.coefficient[i], res.coefficient_idx[i])
                no += 1
                cn = 'TTsinterm'
                selfdisc = ugg()
                selfdisc.summation = list(set(disc_idx))
                selfdisc.coefficient.append(deepcopy(res.coefficient[i]))
                selfdisc.coefficient_idx.append(deepcopy(res.coefficient_idx[i]))
                selfdisc.standarize()
                selfdisc.binary_hash_gen()
                disc_one_dict[deepcopy(cn+str(no))] = deepcopy(selfdisc)
                # print(selfdisc, cn+str(no))
                res.coefficient_idx[i] = []
                res.coefficient[i] = cn+str(no)
                res.summation = [elem for elem in res.summation if elem not in disc_idx]
                
                

    return disc_one_dict, no
    



def pick_best_minibasket(is_present, points):


    n = len(is_present)

    max_points = -1
    for i in range(0, n):
        if is_present[i] == 1:
            if points[i] > max_points:
                max_points = deepcopy(points[i])
                max_points_idx = i
                
    return max_points_idx


def identify_unique_self_disc(outer_terms, disc_one_dict):

    used_keys = []
    unique_dict = {}
    disc_one_dict_unique = {}

    for key in disc_one_dict:
        if key not in used_keys:
            used_keys.append(key)
            disc_one_dict_unique[key] = disc_one_dict[key]
            for key2 in disc_one_dict:
                if key2 != key:
                    if key2 not in used_keys:
                        if disc_one_dict[key2].binary_hash == disc_one_dict[key].binary_hash:
                            unique_dict[key2] = key
                            used_keys.append(key2)


    for i in range(0, len(outer_terms)):
        for k in range(0, len(outer_terms[i].coefficient)):
            if 'TTsinterm' in outer_terms[i].coefficient[k]:
                if outer_terms[i].coefficient[k] in list(unique_dict.keys()):
                    outer_terms[i].coefficient[k] = unique_dict[outer_terms[i].coefficient[k]]

    

    old_new_dict = {}    
    outer_terms, disc_one_dict_unique, old_new_dict,freq_dict = change_names(outer_terms, disc_one_dict_unique, old_new_dict, True)

    print('INTERMEDIATES SELF')
    k = 1
    for key in  disc_one_dict_unique:
        num = int(''.join(filter(str.isdigit, key)))
        disc_one_dict_unique[key].summation = []
        if k%3 == 0:
            print('$', 'K^{{{}}}'.format(num),'$', '=', '$', disc_one_dict_unique[key],'$', '\\\\')
        else:
            print('$','K^{{{}}}'.format(num),'$', '=', '$', disc_one_dict_unique[key], '$', '&')
        k+=1
    print('')

    for i in range(0, len(outer_terms)):
        print(outer_terms[i])

        
    return outer_terms, disc_one_dict_unique
    sys.exit(0)


def identify_square_terms(outer_terms):

    square_dict = {}
    no = 0
    for x in range(0, len(outer_terms)):
        resx = outer_terms[x]
        print('')
        print(resx)
        for y in range(0, len(resx.coefficient)):
            if 'interm' not in resx.coefficient[y]:
                miniugg1 = ugg()
                miniugg1.coefficient.append(resx.coefficient[y])
                miniugg1.coefficient_idx.append(resx.coefficient_idx[y])
                miniugg1.binary_hash_gen()

                for z in range(0, len(resx.coefficient)):
                    if z!= y:
                        if 'interm' not in resx.coefficient[z]:
                            miniugg2 = ugg()
                            miniugg2.coefficient.append(resx.coefficient[z])
                            miniugg2.coefficient_idx.append(resx.coefficient_idx[z])
                            miniugg2.binary_hash_gen()

                            if miniugg1.binary_hash == miniugg2.binary_hash:
                                miniugg3 = miniugg1.fromright(miniugg2)
                                miniugg3.summation = list(set(miniugg2.coefficient_idx[0]))
                                miniugg3.standarize()
                                miniugg3.binary_hash_gen()

                                break_out_flag = False
                                # check if miniugg3 is already in dict
                                for key in square_dict:
                                    if miniugg3.binary_hash ==square_dict[key].binary_hash:
                                        resx.summation = [elem for elem in resx.summation if elem not in miniugg3.summation]
                                        resx.coefficient[y]=key
                                        resx.coefficient_idx[y]=[]
                                        del resx.coefficient[z]
                                        del resx.coefficient_idx[z]
                                        break_out_flag = True
                                        break

                                if not break_out_flag:
                                    no += 1
                                    cn = 'TTqinterm'
                                    square_dict[deepcopy(cn+str(no))] = deepcopy(miniugg3)
                                    resx.summation = [elem for elem in resx.summation if elem not in miniugg3.summation]
                                    resx.coefficient[y]=cn+str(no)
                                    resx.coefficient_idx[y]=[]
                                    del resx.coefficient[z]
                                    del resx.coefficient_idx[z]
                                    break_out_flag2 = True
                                    break
                                else:
                                    break_out_flag2 = True
                                    break
                if break_out_flag2:
                    break

    for x in square_dict:
        print(x, square_dict[x])

    for x in outer_terms:
        x.standarize()
        print(x)

    return outer_terms



def change_names(outer_terms, picked_interms_dict, old_new_dict, all_hash={}, interm_hash={}, disconnected = False):

    no = 0
    cn = 'TTzinterm'

    pattern = r'[a-zA-Z]+'
    first_key = next(iter(picked_interms_dict))
    m = re.search(pattern, first_key)
    if m:
        int_name = m.group(0)

    freq_dict = {}
    if not disconnected:
        old_new_dict_swapped = {v: k for k, v in old_new_dict.items()}
        no = len(old_new_dict)
        
    picked_interms_dict_temp = {}
#    print(old_new_dict)
    for key in picked_interms_dict:

        no += 1
 #       print('keyy', key,no, picked_interms_dict[key])
        newname = cn + str(no)

        if disconnected:
            old_new_dict[key] = int_name+str(no)#newname
        no_minus = False
        for i in range(0, len(outer_terms)):
            for k in range(0, len(outer_terms[i].coefficient)):
                if outer_terms[i].coefficient[k] == key:
                    if not disconnected:
                        if outer_terms[i].coefficient[k] in old_new_dict.keys():
#                            print('tak ten', outer_terms[i].coefficient[k], 'ma sie nazywac', old_new_dict[key])
                            number = int(re.findall(r'\d+',  old_new_dict[key])[-1])
                            newname_from_disc = cn+str(number)
 #                           print('a chcialam go nazwiac', newname, cn+str(number))
                            outer_terms[i].coefficient[k] = newname_from_disc
                            
                            # value = picked_interms_dict_temp.pop(newname)
                            # print('value', value)
                            # picked_interms_dict_temp.update({newname_from_disc: value})
                            newname = newname_from_disc
                            no_minus = True

                            #print('ale teraz picked', newname, 'to', key)
                        else:
#                            print('a tego interm nie bylo wiec', newname)
                            outer_terms[i].coefficient[k] = newname
#                            old_new_dict[key] = int_name+str(no)#newname

                    else:
                        outer_terms[i].coefficient[k] = newname
                        old_new_dict[key] = int_name+str(no)#newname


                        
        picked_interms_dict_temp[newname] = picked_interms_dict[key]
        freq_dict[newname] = 1
        if no_minus:
            no -=1
   #     print('')
  #      print('jeszcze w slowniku', key, 'len', len(picked_interms_dict))
        for keyx, value in picked_interms_dict.items():
            for k in range(0, len(value.coefficient)):
                if value.coefficient[k] == key:
#                    print('zmieniam jeszcze w',  keyx, 'w', value, 'na', newname)
                    value.coefficient[k] = newname
                    freq_dict[newname] += 1
 #                   print('i to',  value)
                    #        print('')


    for i in range(0, len(outer_terms)):
        for k in range(0, len(outer_terms[i].coefficient)):
            if cn in outer_terms[i].coefficient[k]:
                outer_terms[i].coefficient[k] = outer_terms[i].coefficient[k].replace(cn, int_name)

    for value in picked_interms_dict_temp.values():
        for k in range(0, len(value.coefficient)):
            if cn in value.coefficient[k]:
                value.coefficient[k] = value.coefficient[k].replace(cn, int_name)

    picked_interms_dict = {}
    for old_key in picked_interms_dict_temp:
        value = picked_interms_dict_temp[old_key]
        new_key = old_key.replace(cn, int_name)
        picked_interms_dict[new_key] = value


    freq_dict_new = {}
    for old_key in freq_dict:
        value = freq_dict[old_key]
        new_key = old_key.replace(cn, int_name)
        freq_dict_new[new_key] = value

    freq_dict = dict(sorted(freq_dict_new.items(), key=lambda item: item[1], reverse=True))
    print('freq_dict')
    for x in freq_dict:
        print(x, freq_dict[x])


#    picked_interms_dict = picked_interms_dict_temp


    return outer_terms, picked_interms_dict, old_new_dict, freq_dict_new

# def free_no(old_new_dict, no):

#     numbers_list = [int(re.search(r'\d+', value).group()) for value in old_new_dict.values()]

#     print('number_list', number_list)


def change_names_drag_old_dict(outer_terms, picked_interms_dict, old_new_dict):

    no = 0
    cn = 'TTzinterm'

    pattern = r'[a-zA-Z]+'
    first_key = next(iter(picked_interms_dict))
    m = re.search(pattern, first_key)
    if m:
        int_name = m.group(0)


    picked_interms_dict_temp = {}
    old_new_dict = {}
    for key in picked_interms_dict:
        no += 1
        newname = cn + str(no)
        picked_interms_dict_temp[newname] = picked_interms_dict[key]
        old_new_dict[key] = newname
        for i in range(0, len(outer_terms)):
            for k in range(0, len(outer_terms[i].coefficient)):
                if outer_terms[i].coefficient[k] == key:
                    outer_terms[i].coefficient[k] = newname
                    old_new_dict[key] = newname

        for value in picked_interms_dict.values():
            for k in range(0, len(value.coefficient)):
                if value.coefficient[k] == key:
                    value.coefficient[k] = newname
#        print('')
        
    for i in range(0, len(outer_terms)):
        for k in range(0, len(outer_terms[i].coefficient)):
            if cn in outer_terms[i].coefficient[k]:
                outer_terms[i].coefficient[k] = outer_terms[i].coefficient[k].replace(cn, int_name)

    for value in picked_interms_dict_temp.values():
        for k in range(0, len(value.coefficient)):
            if cn in value.coefficient[k]:
                value.coefficient[k] = value.coefficient[k].replace(cn, int_name)



    picked_interms_dict = picked_interms_dict_temp


    return outer_terms, picked_interms_dict, old_new_dict



def print_interms(picked_interms_dict, outer_terms):


    interm_whichlevel_dict = {} # 0 - base interm, 1 - interm that has at least one interm0, 2 - interm that has at least one interm1...
    maxl = 10
    true_max = 0

    for key, value in picked_interms_dict.items():
        notzero = False
        for coef in value.coefficient:
            if 'interm' in coef:
                notzero = True
                break
        if not notzero:
            interm_whichlevel_dict[key] = 0

    for lev in range(1, maxl):
        for key, value in picked_interms_dict.items():
            addthis = False
            maxlev = 0
            for coef in value.coefficient:
                thislev = 0
                if 'interm' in coef:
                    if coef in interm_whichlevel_dict.keys():
                        addthis = True
                        thislev +=  interm_whichlevel_dict[coef] + 1
                        if thislev>maxlev:
                            maxlev = deepcopy(thislev)
                            if thislev>true_max:
                                true_max = deepcopy(thislev)

            if addthis:
                interm_whichlevel_dict[key] = maxlev


    # for key in interm_whichlevel_dict:
    #     print('interm', key, picked_interms_dict[key], 'is at level', interm_whichlevel_dict[key])


    lev2 = {}

    for key, value in interm_whichlevel_dict.items():
        if value not in lev2:
            lev2[value] = [key]
        else:
            lev2[value].append(key)

    # print('lev2', lev2)


    # for x in outer_terms:
    #     print(x)
    # zapiszmy wszystkie outer ktore maja interm 0 levelu
    # used_idx = []

    # for n in range(0, true_max):
    #     # jestem na levelu n
