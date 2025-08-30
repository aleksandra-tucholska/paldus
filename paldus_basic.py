from params import *
import paldus_classes
from paldus_classes import ugg
from paldus_cas import cas
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
from random import shuffle
import sys
import time
import pickle
from joblib import Parallel, delayed

# from paldus_density import cost_after_intermediates



def identify_equal_bin(cluster):

    for l in range(0, len(cluster)):
   #     print('cluster[l],', cluster[l])
        if abs(cluster[l].num_factor) > EPSILON :
            for k in range(l + 1, len(cluster)):
                if abs(cluster[k].num_factor) > EPSILON:
                    phase = cluster[l].binary_hash == cluster[k].binary_hash
  #                  print('cluster[k],', cluster[k], phase)
 #                   print()
                    if phase != 0:
#                        print('phase', phase, cluster[l])
                        cluster[l].num_factor = cluster[l].num_factor + phase * cluster[k].num_factor
                        cluster[k].num_factor = 0.0
                        # print(cluster[l])
                        # print('----------')

    return cluster


def identify_equal(cluster):
    """ In given cluster, identifies identical uqq elements.
    If identical elements are found, num_factor of uqq1 changes,
    and selector[place[uqq2]] is set to 0"""

    for l in range(0, len(cluster)):
        if abs(cluster[l].num_factor) > EPSILON :
            for k in range(l + 1, len(cluster)):
                if abs(cluster[k].num_factor) > EPSILON:
                    phase = cluster[l] == cluster[k]

                    if phase != 0:
                        cluster[l].num_factor = cluster[l].num_factor + phase * cluster[k].num_factor
                        cluster[k].num_factor = 0.0

    return cluster

def identify_equal_two_clusters(clist):
    
    cluster1 = clist[0]
    cluster2 = clist[1]
    
#    print('')
    for l in range(0, len(cluster1)):
        if abs(cluster1[l].num_factor) > EPSILON :
            for k in range(0, len(cluster2)):
 #               print('asdasd', cluster2[k])
                if abs(cluster2[k].num_factor) > EPSILON:
                    phase = cluster1[l] == cluster2[k]
                    if phase != 0:
                        cluster1[l].num_factor = cluster1[l].num_factor + phase * cluster2[k].num_factor
                        cluster2[k].num_factor = 0.0

  #  print('')
#    print(len(cluster1), len(cluster2))
    
    cluster1.cleanup()
    cluster2.cleanup()

#    print(len(cluster1), len(cluster2))
#    print('')

    clist = [cluster1, cluster2]

    return clist

def delete_permutations(cluster):

    len0 = len(cluster)
    len1 = 0
    for l in range(0, len(cluster)):
        if abs(cluster[l].num_factor) > EPSILON :

            perm_cluster = deepcopy(cluster[l])
            perm_cluster.multisubst(['i','j','k'],['j','k','i'])

            perm_cluster.standarize()

            
            for k in range(0, len(cluster)):
                if abs(cluster[k].num_factor) > EPSILON:
                    phase = perm_cluster == cluster[k]

                    if phase != 0:
                        cluster[k].num_factor = 0.0
                        len1 += 1

    if len0/2 != len1:
        print('NOT ALL PERMUTATIONS WERE DELETED')
        sys.exit(1)

    return cluster

def evaluate(*args0):
    """Evaluate nested commutator,
    evaluate(A, B, D, ...) <- [[...[A, B], C]...],
    where A, B, C, D are arithmetic_string instances
    or ugg instances.
    """
    args = []
    for x in args0:
        if isinstance(x, ugg):
            #print('oto x ugg ', x)
            args.append(arithmetic_string(x))
        elif isinstance(x, cas):
            args.append(arithmetic_string(x))
        else:
            #print('oto x arg ', x)
            args.append(x)

    b = arithmetic_string()

    a = product(*args)

    for x in a:
        # print('x', x)
        s = arithmetic_string(*x)
        #print('s', s)
        s.disambiguate()
        #print('po disambiguate w evaluate')
        #print(s)
        #sys.exit(0)
        b = b + nested_commutator(*s)
        # print('bbb', b)
        
    b.exec_delta()

    #print('----------------------------')
#    sys.exit(0)
    return b

def evaluate_gen(*args0):
    """Evaluate nested commutator,
    evaluate(A, B, D, ...) <- [[...[A, B], C]...],
    where A, B, C, D are arithmetic_string instances
    or ugg instances.
    """
    args = []
    for x in args0:
        if isinstance(x, ugg):
            args.append(arithmetic_string(x))
        elif isinstance(x, cas):
            args.append(arithmetic_string(x))
        else:
            args.append(x)

    b = arithmetic_string()

    a = product(*args)

    for x in a:
        s = arithmetic_string(*x)
        s.disambiguate()
        b = b + nested_commutator(*s)
    b.exec_delta(True)

    return b


def nested_commutator(*args):
    """
    nested_commutator(A, B, D, ...) <- [[...[A, B], C]...],
    where A, B, C, D are ugg instances.
    """
    nest_order = len(args) - 1
    left = args[0:len(args) - 1]
    right = args[len(args) - 1]


    if nest_order > 1:
        a = nested_commutator(*left)
        b = arithmetic_string()
        for x in a:
            if isinstance(x, ugg):
                #print('ugggg nest >1')
                b = b + commute(x, right)
            elif isinstance(x, cas):
                b = b + commute_cas(x, right)
        return b
    else:
        if isinstance(right, ugg):
            #print('uggg nest =1')
            z = commute(left[0], right)
        elif isinstance(right, cas):

            z = commute_cas(left[0], right)            
        return z


def commute(e1, e2):
    """
    Compute [E1, E2] commutator, where E1, E2 are
    ugg instances.
    """

    if len(e2.operator_idx) > 1:
        temp1, e2_r = e2.left_split()
        return commute(e1, temp1).fromright(e2_r) + commute(e1, e2_r).fromleft(temp1)

    if len(e1.operator_idx) > 1:
        temp1, e1_r = e1.left_split()
        return commute(e1_r, e2).fromleft(temp1) + commute(temp1, e2).fromright(e1_r)

    #print('op e0')
    #print(e1, e2)
    return basic_commute(e1, e2)


def basic_commute(e1, e2):

#    print('-------TUTAJ------')
#    print(e1, e2)
    res1 = ugg()
    res2 = ugg()

    res1.coefficient = e1.coefficient + e2.coefficient
    res2.coefficient = e1.coefficient + e2.coefficient

    res1.num_factor = e1.num_factor * e2.num_factor
    res2.num_factor = -e1.num_factor * e2.num_factor

    res1.summation = e1.summation + e2.summation
    res2.summation = e1.summation + e2.summation

    res1.operator_idx = [[e1.operator_idx[0][0] , e2.operator_idx[0][1]]]
    res2.operator_idx = [[e2.operator_idx[0][0], e1.operator_idx[0][1] ]]

    # if e1.operator_type[0] == "t" and e2.operator_type[0] == "t":
    #     res1.operator_type.append("s")
    #     res2.operator_type.append("s")
    # elif e1.operator_type[0] == "s" and e2.operator_type[0] == "t":
    #     res1.operator_type.append("t")
    #     res2.operator_type.append("t")
    # elif e1.operator_type[0] == "t" and e2.operator_type[0] == "s":
    #     res1.operator_type.append("t")
    #     res2.operator_type.append("t")
    # elif e1.operator_type[0] == "s" and e2.operator_type[0] == "s":
    #     res1.operator_type.append("s")
    #     res2.operator_type.append("s")


    if e1.operator_type[0] == "t0" and e2.operator_type[0] == "t0":
        res1.operator_type.append("s")
        res2.operator_type.append("s")
    elif e1.operator_type[0] == "s" and e2.operator_type[0] == "t0":
        res1.operator_type.append("t0")
        res2.operator_type.append("t0")
    elif e1.operator_type[0] == "t0" and e2.operator_type[0] == "s":
        res1.operator_type.append("t0")
        res2.operator_type.append("t0")
    elif e1.operator_type[0] == "s" and e2.operator_type[0] == "s":
        res1.operator_type.append("s")
        res2.operator_type.append("s")


    res1.coefficient_idx = e1.coefficient_idx + e2.coefficient_idx
    res2.coefficient_idx = e1.coefficient_idx + e2.coefficient_idx

    res1.new_delta(e2.operator_idx[0][0], e1.operator_idx[0][1])
    res2.new_delta(e2.operator_idx[0][1], e1.operator_idx[0][0])

    if e1.delta != []:
        for delta in e1.delta:
            res1.new_delta(delta[0], delta[1])
            res2.new_delta(delta[0], delta[1])

    if e2.delta != []:
        for delta in e2.delta:
            res1.new_delta(delta[0], delta[1])
            res2.new_delta(delta[0], delta[1])

 #   print('----I TUTAJ--------')
 #   print(res1, res2)
    result = arithmetic_string(res1, res2)

    return result


def res_compact(a):
    """
    Adds equal components, and adds permutation operator
    before permuted components. 
    """
    start = time.time()
    b = a.expand()
    b.cluster()
    result = arithmetic_string()
    end = time.time()
    start = time.time()
    k = 0
    for key in b.clusters:
        k += 1
        # print(k)
        cluster = b.clusters[key]
        if len(cluster) > 1:
            cluster = identify_equal(cluster)
            cluster.cleanup()

        for elem in cluster:
            result.append(elem)

    return result



def res_compact2(a):
    """                                                                                             
    Adds equal components, and adds permutation operator                                           
    before permuted components.                                
    """
#    t1 = time.time()

    b = a.expand()
#    t2 = time.time()
#    print('exp', t2-t1)
#    t1 = time.time()
    b.cluster()
#    t2 = time.time()
#    print('cluster', t2-t1)
#    t1 = time.time()

    result = arithmetic_string()
    
    
    k = 0
    for key in b.clusters:
        k += 1
        cluster = b.clusters[key]
 #       start1 = time.time()
        cluster = simplify_large_cluster2(cluster)

        for elem in cluster:
            result.append(elem)
        # end1 = time.time()
        # print('lenlen', len(cluster), end1-start1)
  #  t2 = time.time()
   # print('simplarge', t2-t1)
    #t1 = time.time()

    #end = time.time()
#    print('identify', end-start)

    return result

def simplify_large_cluster2(cluster):
    min_len = 20
    div_len = 200


    for x in cluster:
        x.binary_hash_gen()
    cluster = identify_equal_bin(cluster)

    return cluster


#     i = 0
#     if len(cluster) > min_len:
#         cluster_div = []
#         minicluster = arithmetic_string()
#         for x in cluster:
#             i += 1
#             minicluster.append(x)
#             if (i%div_len == 0):
#                 print(len(minicluster))
#                 cluster_div.append(minicluster)
#                 minicluster = arithmetic_string()

#         start = time.time()

#         z = Parallel(n_jobs=30,verbose=100)(delayed(identify_equal)(cluster_div[j]) for j in range(0, len(cluster_div)))

# #        end = time.time()
# #        print('czas na parallel POJEDYNCZEGO', end-start)

#         for j in range(0, len(z)):
#             z[j].cleanup()
            
#         ln = len(z)
#         if ln > 1:
#             if ln%2 != 0:
#                 for j in range(0, len(z[ln-1])):
#                     z[ln-2].append(z[ln-1][j])
#                 z = z[0:ln-1]

#             ln2 = len(z)

#             clust_even=[]
#             clust_odd=[]
#             cluster_list_big = []
#             for j in range(0, ln2):
#                 if j%2 ==0:
#                     clust_even = z[j]
#                 else:
#                     clust_odd = z[j]
#                     cluster_list_big.append([deepcopy(clust_even), deepcopy(clust_odd)])


#             print('PLUUUUUSZON', ln2)

#             for l in range(0, ln2//2):
#                 print('llllllllllllllllllllllllllll', l, 'z', ln2/2)
#                 start = time.time()
#                 if l != 0:
#                     clust_even = []
#                     clust_odd = []
#                     for ii in range(0, ln2//2):
#                         clust_even.append(deepcopy(z2[ii][0]))
#                         clust_odd.append(deepcopy(z2[ii][1]))

#                     first = deepcopy(clust_odd[0])
#                     clust_odd.remove(clust_odd[0])
#                     clust_odd.append(first)

#                     cluster_list_big = []
#                     for ii in range(0, ln2//2):
#                         cluster_list_big.append([clust_even[ii], clust_odd[ii]])
#                     stmal = time.time()
#                     z = Parallel(n_jobs=30,verbose=0)(delayed(identify_equal_two_clusters)(cluster_list_big[j]) for j in range(0, ln2//2))
#                     z2 = deepcopy(z)
#                     end = time.time()
#                     print('czas na pojedynczy chunk', end-start, end-stmal)
#                 else:
                    
#                     z = Parallel(n_jobs=30,verbose=0)(delayed(identify_equal_two_clusters)(cluster_list_big[j]) for j in range(0, ln2//2))
#                     z2 = deepcopy(z)

#     cluster_res = []
#     for i in z:
#         for j in i:
#             for k in j:
#                 cluster_res.append(k)

#     return cluster_res
        

def simplify(result, no_fx=None, last_idx='', fixed_fx = [], cas=False, drag_fixed = []):
    """ Simplifies result, i.e:
    1. Standarizes
    2. Compact: Identifies equal components and permuted components
    3. Findes the shortest version of result, by making all
    possible arrangemenst of E operators.
    """

    # t1 = time.time()
    if drag_fixed!=[]:
        result.clear_fixed()
    else:
        if not no_fx:
            result.clear_fixed()
            result.establish_fixed()
        
        if fixed_fx != []:
            result.clear_fixed()
            result.establish_fixed(fixed_fx)


    # t2 = time.time()
    # print()
    # print('ple1 time', t2-t1)
    # t1 = time.time()
    result.as_standarize(last_idx, fixed_fx, cas, drag_fixed)
    # t2 = time.time()
    # print()
    # print('standardize time', t2-t1)
    # print()
    # if (t2-t1)> 1:
    #     sys.exit(0)
    # print('mandalore2-wystandaryzowane')
    # k = 0
    # for x in result:
    #     print(k, x)
    #     k += 1
        
    # print('')
    #    sys.exit(0)
    # print('stand', result)
    # print('wynikissimo')
    # t1 = time.time()
    result.cleanup()
    # t2 = time.time()
    # print()
    # print('cleanup', t2-t1)
    # k  = 1
    # for x in result:
    #     print(k, x)
    #     k += 1
#    sys.exit(0)
#    result = res_compact(result)
    # t1 = time.time()

    result = res_compact2(result)

    # t2 = time.time()
    # print()
    # print('res_compact time', t2-t1)
    # t1 = time.time()
    # print('po res_compact2')
    result.cleanup()
    # k  = 1
    # for x in result:
    #     print(k, x)
    #     k += 1

    # print('tt', result)
    if not no_fx:
        free_fixed()
    # print('zz', result)


    n_g3 = 0
    n_g4 = 0
    
    # for x in result:
    #     for i, y in enumerate(x.coefficient):
    #         if y == DENS3:
    #             n_g3+=1
    #         elif y == DENS4:
    #             n_g4+=1

    # print('number of gamm3 after standarize', n_g3)
    # print('number of gamm4 after standarize', n_g4)
    # t2 = time.time()
    # print()
    # print('restt', t2-t1)

    return result

def simple_swap_count(lista):

    liczba_swapow = 0
    for i in range(len(lista)):
        for j in range(i + 1, len(lista)):
            if lista[i] > lista[j]:
                liczba_swapow += 1
    return liczba_swapow


def simplify_final_touch(res):

    for i, elem in enumerate(res):
        s = elem.summation
        p = list(permutations(s))

        new_elem = deepcopy(elem)
        print(elem)

        original_count = 0
        for j, coef_idx in enumerate(elem.coefficient_idx):
            original_count += simple_swap_count(coef_idx)
        print('original_count', original_count)
        if len(p)>1:
            original = list(p[0])

        for k in range(1, len(p)):
            print(k, p[k])
            swap = list(p[k])
            maps = {org: swp for org, swp in zip(original, swap)}

            for j, l in enumerate(new_elem.coefficient_idx):
               new_elem.coefficient_idx[j] = [maps.get(element, element) for element in l]
            new_elem.orbital_type = {maps.get(k, k): maps.get(v, v) for k, v in new_elem.orbital_type.items()}
#            new_elem.orbital_type = [maps.get(element, element) for element in l]
            print(new_elem.coefficient)
            new_elem.standarize_g()
            new_elem.standarize_v()
            new_elem.standarize_gs()
            new_elem.standarize_gm2()
            new_elem.standarize_gm34()
            new_elem.standarize_gm234spin_new()
            print(new_elem)
            this_elem_count = 0
            for j, coef_idx in enumerate(new_elem.coefficient_idx):
                this_elem_count += simple_swap_count(coef_idx)
            print('this_count', this_elem_count)
            if (this_elem_count < original_count):
                original_count = deepcopy(this_elem_count)
                res[i] = deepcopy(new_elem)

        print()

        # for r in p:
        #     print(list[r])
            #print(i, s, elem)
        print()
    print('finally')
    for x in res:
        print(x)
    print()
    res = res_compact2(res)
    res.cleanup()
    print('finally2')
    for x in res:
        print(x)
    return res


def find_antysym_int(res0):

    banned = []
    
    for i in range(0, len(res0)):
        # print('')
        # print('------------------------------------------------------------------------------------------')
        # print('')
        # print('teraz bedzie', res0[i])
        # remove the twoel int
        if i not in banned:
            l1 = -1
            if TWOEL_INT_DIRAC in res0[i].coefficient:
                term1 = ugg()
                term1.num_factor = res0[i].num_factor
                term1.delta = res0[i].delta
                term1.summation = res0[i].summation
                for k in range(0, len(res0[i].coefficient)):
                    if res0[i].coefficient[k] != TWOEL_INT_DIRAC:
                        term1.coefficient.append(res0[i].coefficient[k])
                        term1.coefficient_idx.append(res0[i].coefficient_idx[k])
                    else:
                        l1 = deepcopy(k)
                # print('taktuuuu', term1)

                for j in range(0, len(res0)):
                    # print('teraz bedzie res2', res0[j])
                    # remove the twoel int
                    if j not in banned and j != i:
                        l2 = -1
                        if TWOEL_INT_DIRAC in res0[j].coefficient:
                            term2 = ugg()
                            term2.num_factor = res0[j].num_factor
                            term2.delta = res0[j].delta
                            term2.summation = res0[j].summation
                            for k in range(0, len(res0[j].coefficient)):
                                if res0[j].coefficient[k] != TWOEL_INT_DIRAC:
                                    term2.coefficient.append(res0[j].coefficient[k])
                                    term2.coefficient_idx.append(res0[j].coefficient_idx[k])
                                else:
                                    l2 = deepcopy(k)

                            # print('termterm')
                            # print(term1)
                            # print(term2)
                            if term1 == term2:
                                if res0[i].num_factor > 0 and res0[j].num_factor < 0:
                                    # print('w tym pierwszym')
                                    x =find_twoel_dirac_asym(res0[i].coefficient_idx[l1], res0[j].coefficient_idx[l2])
                                    term1.coefficient.append(TWOEL_INT_DIRAC_A)
                                    term1.coefficient_idx.append(x)
                                    res0[j].num_factor = 0
                                    res0[i] = deepcopy(term1)
                                    banned.append(i)
                                    banned.append(j)
                                    # print(res0[i], # print(res0[j]))
                                    break
                                elif res0[i].num_factor < 0 and res0[j].num_factor > 0:
                                    # print('w tym drugim')
                                    x = find_twoel_dirac_asym(res0[j].coefficient_idx[l2], res0[i].coefficient_idx[l1])
                                    term2.coefficient.append(TWOEL_INT_DIRAC_A)
                                    term2.coefficient_idx.append(x)
                                    res0[i].num_factor = 0
                                    res0[j] = deepcopy(term2)                                                                    
                                    banned.append(i)
                                    banned.append(j)
                                    # print(# print(res0[i]), res0[j])
                                    break
                                else:
                                    # print('COS JEST NAROBIONE')
                                    sys.exit(0)

    res = arithmetic_string()
    for x in res0:
        if x.num_factor != 0:
            res.append(x)

    # print('wynik z  calkami asym')
    # for x in res:
    #     print(x)

    return res

                                

def find_twoel_dirac_asym(list1, list2):

    a1 = [list1[0], list1[3], list1[2], list1[1]]
    a2 = [list1[1], list1[0], list1[3], list1[2]]
    a3 = [list1[1], list1[2], list1[3], list1[0]]
    a4 = [list1[2], list1[3], list1[0], list1[1]]
    a5 = [list1[2], list1[1], list1[0], list1[3]]
    a6 = [list1[3], list1[2], list1[1], list1[0]]
    a7 = [list1[3], list1[0], list1[1], list1[2]]

    list1_big = [list1, a1, a2, a3, a4, a5, a6, a7]


    b1 = [list2[0], list2[3], list2[2], list2[1]]
    b2 = [list2[1], list2[0], list2[3], list2[2]]
    b3 = [list2[1], list2[2], list2[3], list2[0]]
    b4 = [list2[2], list2[3], list2[0], list2[1]]
    b5 = [list2[2], list2[1], list2[0], list2[3]]
    b6 = [list2[3], list2[2], list2[1], list2[0]]
    b7 = [list2[3], list2[0], list2[1], list2[2]]

    list2_big = [list2, b1, b2, b3, b4, b5, b6, b7]

    for x in list1_big:
        for y in list2_big:
            z = deepcopy(y)
            temp = deepcopy(z[3])
            z[3] = deepcopy(z[2])
            z[2] = deepcopy(temp)
            if x == z:
                return(x)


def find_L(result, typ):

    # print('przed find L')

    if typ == 0:
        coef = 2.0
    else:
        coef = 1.0

    result_group = []
    rg = arithmetic_string()
    # print(len(result))
    k = 1
    # print(len(result), 'przed as standarize w LLLL')
    result.as_standarize(TWOEL_INT_TRANS)
    # print(len(result), 'po as standarize w LLLL')
    rsimp = simplify(result, None, TWOEL_INT_TRANS)
    # print(len(rsimp), 'po simplify w LLLL')
    result = deepcopy(rsimp)
    k = 1
    for x in result:
        print(k, x)
        k += 1
    

    for x in result:
        newres = deepcopy(x)


        newres.coefficient = []
        newres.coefficient_idx = []

        for s in range(0, len(x.coefficient)):  
            if x.coefficient[s] != TWOEL_INT and x.coefficient[s] != TWOEL_INT_TRANS:           
                newres.coefficient.append(x.coefficient[s])
                newres.coefficient_idx.append(x.coefficient_idx[s])
        result_group.append(newres)
        rg += arithmetic_string(newres)
        # print(k, 'nwr', x, newres)
        k += 1

    # print(len(result_group), 'lenlnelnelnelne')
    for x in range(0, len(result)):
        print(result[x], result_group[x])

    print('wynik PRZED', len(result))
    k = 0
    for x in result:
        print(k, x)
        k += 1
    print('')

    print('wyniki')
    k = 1
    for x in range(0, len(result_group)):
        print('teraz sprawdzam ten', k, result_group[x])
        l = k + 1
        k += 1
        for y in range(x+1, len(result_group)):
            if len(result_group[x].summation)==len(result_group[y].summation):
                print('z tym sprawdzam ten', l, result_group[y])

            l += 1
            hash1 = result_group[x].hash
            hash2 = result_group[y].hash
            if hash1 == hash2:
                print('hashe sa rowne')
                print(result[x])
                print(result[y])
                listx = [152]
                listy = [21]
                for z in range(0, len(result[x].coefficient)):
                    if result[x].coefficient[z] == TWOEL_INT or result[x].coefficient[z] == TWOEL_INT_TRANS:
                        listx = deepcopy(result[x].coefficient_idx[z])
                for z in range(0, len(result[y].coefficient)):
                    if result[y].coefficient[z] == TWOEL_INT or result[y].coefficient[z] == TWOEL_INT_TRANS:
                        listy = deepcopy(result[y].coefficient_idx[z])

                listx.sort()
                listy.sort()
                print('listx', listx)
                print('listy', listy)
                if listx == listy:
                    print('listy sa rowne')
                    numx = result_group[x].num_factor
                    numy = result_group[y].num_factor
                    print(numx, numy)
                    if TWOEL_INT in result[x].coefficient or TWOEL_INT_TRANS in result[x].coefficient and \
                            TWOEL_INT in result[y].coefficient or TWOEL_INT_TRANS in result[y].coefficient:
                        if abs((coef * numx + numy)) < 0.0001:
                            print( 'i raz')
                            result[x].num_factor = 0.0
                            result[y].num_factor = result[y].num_factor /2.0
                            for z in range(0, len(result[y].coefficient)):
                                if result[y].coefficient[z] == TWOEL_INT or result[y].coefficient[z] == TWOEL_INT_TRANS:
                                    result[y].coefficient[z] = TWOEL_INT_COMBO
                                    print(result[y])
                        elif abs((coef * numy + numx))< 0.0001:
                            print('mamy 2 w x i -1 w y')
                            result[y].num_factor = 0.0
                            result[x].num_factor = result[x].num_factor /2.0
                            for z in range(0, len(result[x].coefficient)):
                                if result[x].coefficient[z] == TWOEL_INT or result[x].coefficient[z] == TWOEL_INT_TRANS:
                                    result[x].coefficient[z] = TWOEL_INT_COMBO
                                    print(result[x])
                        print('tak')
                        print(result[x])
                        print(result[y])
                        print('')
                break
        
    print(len(result))
    print('wynik przed cleanup', len(result))
    k = 1
    for x in result:
        print(k, x)
        k+= 1
    print('')

    result.cleanup()
    result.as_standarize(TWOEL_INT_TRANS)
    print('wynik ost', len(result))
    for x in result:
        print(x)
    return result
    # b = rg.expand()
    # b.cluster()
    # print('clusters')
    # for key in b.clusters:
    #     for f in b.clusters[key]:
    #         print(f)
    #     print('')


def find_F(result, typ):

    print('przed find F')

    if typ == 0:
        coef = 2.0
    else:
        coef =1.0


    result_group = []
    rg = arithmetic_string()
    print(len(result))
    k = 1
    result.as_standarize(F12_TWOEL)
    print(len(result), 'po as standarize w FFFF')
    for x in result:
        print(x)


    # rsimp = simplify(result, F12_TWOEL)
    # result = deepcopy(rsimp)

    for x in result:
        newres = deepcopy(x)
        newres.coefficient = []
        newres.coefficient_idx = []

        for s in range(0, len(x.coefficient)):  
            if x.coefficient[s] != F12_TWOEL: 
                newres.coefficient.append(x.coefficient[s])
                newres.coefficient_idx.append(x.coefficient_idx[s])
        result_group.append(newres)
        rg += arithmetic_string(newres)
        print(k, 'nwr', x, newres)
        k += 1

    print(len(result_group), 'lenlnelnelnelne')
    for x in range(0, len(result)):
        print(result[x], result_group[x])

    print('wynik PRZED', len(result))
    k = 0
    for x in result:
        print(k, x)
        k += 1
    print('')

    print('wyniki')
    k = 1
    for x in range(0, len(result_group)):
        print('teraz sprawdzam ten F', k, result_group[x])
        l = k + 1
        k += 1
        for y in range(x+1, len(result_group)):
            if len(result_group[x].summation)==len(result_group[y].summation):
                print('z tym sprawdzam ten F', l, result_group[y])

            l += 1
            hash1 = result_group[x].hash
            hash2 = result_group[y].hash
            if hash1 == hash2:
                print('hashe sa rowne')
                listx = [152]
                listy = [21]
                for z in range(0, len(result[x].coefficient)):
                    if result[x].coefficient[z] == F12_TWOEL:
                        listx = deepcopy(result[x].coefficient_idx[z])
                for z in range(0, len(result[y].coefficient)):
                    if result[y].coefficient[z] == F12_TWOEL:
                        listy = deepcopy(result[y].coefficient_idx[z])

                listx.sort()
                listy.sort()
                print('listx', listx)
                print('listy', listy)
                if listx == listy:
                    print('listy sa rowne')
                    numx = result_group[x].num_factor
                    numy = result_group[y].num_factor
                    print(numx, numy, abs((coef * numx + numy)), abs((coef * numy + numx)))
                    if abs((coef * numx + numy)) < 0.0001:
                        print( 'i raz')
                        result[x].num_factor = 0.0
                        result[y].num_factor = result[y].num_factor /2.0
                        for z in range(0, len(result[y].coefficient)):
                            if result[y].coefficient[z] == F12_TWOEL:
                                result[y].coefficient[z] = F12_TWOEL_COMBO
                    elif abs((coef * numy + numx))< 0.0001:
                        print('mamy 2 w x i -1 w y')
                        result[y].num_factor = 0.0
                        result[x].num_factor = result[x].num_factor /2.0
                        for z in range(0, len(result[x].coefficient)):
                            if result[x].coefficient[z] == F12_TWOEL:
                                result[x].coefficient[z] = F12_TWOEL_COMBO

                    print('tak')
                    print(result[x])
                    print(result[y])
                    print('')
                break
        
    print(len(result))
    print('wynik przed cleanup', len(result))
    k = 1
    for x in result:
        print(k, x)
        k+= 1
    print('')

    result.cleanup()
    result.as_standarize(F12_TWOEL)
    print('wynik ost', len(result))
    for x in result:
        print(x)
    return result
    # b = rg.expand()
    # b.cluster()
    # print('clusters')
    # for key in b.clusters:
    #     for f in b.clusters[key]:
    #         print(f)
    #     print('')


def identify_interm_V_f12(res):
    
    result = []
    
    
    for x in res:
        append = False
        excluded1 = []
        sx = []
        for i in range(0, len(x.coefficient)):
            y = x.coefficient[i]
            if y == F12_TWOEL or y == F12_TWOEL_COMBO:
                excluded1.append(i)
                n_of_completev = 0
                lst_of_completev = []
                rest_f_idx = []
                for idx in x.coefficient_idx[i]:
                    if idx in completev:
                        n_of_completev += 1
                        lst_of_completev.append(idx)
                    else:
                        rest_f_idx.append(idx)
                lst_of_completev.sort()
                if n_of_completev == 2:
                    for i2 in range(0, len(x.coefficient)):
#                        print(x.coefficient[i2])
                        if i2 not in excluded1:
                            y2 = x.coefficient[i2]
                            if y2 == TWOEL_INT or y2 == TWOEL_INT_TRANS or y2==TWOEL_INT_COMBO:
#                                print('mam y2', x, y2)
#                                print('a moze tu?')
                                lst_of_completev2 = []
                                n_of_completev2 = 0
                                rest_v_idx = []
                                for idx in x.coefficient_idx[i2]:
                                    if idx in completev:
                                        n_of_completev2 += 1
                                        lst_of_completev2.append(idx)
                                    else:
                                        rest_v_idx.append(idx)
                                lst_of_completev2.sort()
                                if n_of_completev2 ==2:
                                    if lst_of_completev == lst_of_completev2:
                                        sx.append(summed_idx([x.coefficient_idx[i], x.coefficient_idx[i2]], x.summation))
#                                        print('TU BEDZIE INTERMEDIATE V', x)
                                        resmini = x.copy_without([i, i2])
                                        if y == F12_TWOEL  and y2 in TWOEL:
                                            resmini.coefficient.append(INTERM_V_F12)
                                        elif y == F12_TWOEL  and y2 == TWOEL_INT_COMBO:
                                            resmini.coefficient.append(INTERM_V_F12+"_2")
                                        elif y == F12_TWOEL_COMBO  and y2 in TWOEL:
                                            resmini.coefficient.append(INTERM_V_F12+"_3")
                                        elif y == F12_TWOEL_COMBO  and y2 == TWOEL_INT_COMBO:
                                            resmini.coefficient.append(INTERM_V_F12+"_4")
                                        resmini.coefficient_idx.append(rest_f_idx+rest_v_idx)
                                        resmini.num_factor *= (2.0)
 #                                       print('resmini', resmini)
                                        append = True
        if append == True:
            sx_all = []
            for l in sx:
                for ll in l:
                    sx_all.append(ll)

            for s in x.summation:
                if s not in sx_all:
                    # print('dodaje do sum resmini', s, sx_all)
                    resmini.summation.append(s)

            # print('intermedaiate V = v* ff')
            # print('')
            # print(x)
            # print(resmini)
            # print('')
            result.append(resmini)
        else:
            result.append(x)
        k = 0

    # print('result z intermediate')
    # for x in result:
    #     print(x)

    return result


def identify_interm_X_f12(res):

    result = []
    for x in res:
        append = False
        excluded1 = []
        sx = []
        for i in range(0, len(x.coefficient)):
            y = x.coefficient[i]
            if y == F12_TWOEL or y == F12_TWOEL_COMBO:
                excluded1.append(i)
                n_of_completev = 0
                lst_of_completev = []
                rest_f_idx = []
                for idx in x.coefficient_idx[i]:
                    if idx in completev:
                        n_of_completev += 1
                        lst_of_completev.append(idx)
                    else:
                        rest_f_idx.append(idx)
                lst_of_completev_sort = deepcopy(lst_of_completev)
                lst_of_completev_sort.sort()
                if n_of_completev == 2:
                   for i2 in range(0, len(x.coefficient)):
                       if i2 not in excluded1:
                           y2 = x.coefficient[i2]
                           if y2 == F12_TWOEL or y2 == F12_TWOEL_COMBO:
                               lst_of_completev2 = []
                               n_of_completev2 = 0
                               rest_v_idx = []
                               for idx in x.coefficient_idx[i2]:
                                   if idx in completev:
                                       n_of_completev2 += 1
                                       lst_of_completev2.append(idx)
                                   else:
                                       rest_v_idx.append(idx)
                                       
                               lst_of_completev2_sort = deepcopy(lst_of_completev2)
                               lst_of_completev2_sort.sort()
                            
                               if n_of_completev2 ==2:
                                   if lst_of_completev == lst_of_completev2:
                                       sx.append(summed_idx([x.coefficient_idx[i], x.coefficient_idx[i2]], x.summation))

                                       resmini = x.copy_without([i, i2])
                                       if y == F12_TWOEL  and y2 == F12_TWOEL:
                                           resmini.coefficient.append(INTERM_X_F12)
                                       elif y == F12_TWOEL  and y2 == F12_TWOEL_COMBO:
                                           resmini.coefficient.append(INTERM_X_F12+"_2")
                                       elif y == F12_TWOEL_COMBO  and y2 == F12_TWOEL:
                                           resmini.coefficient.append(INTERM_X_F12+"_3")
                                       elif y == F12_TWOEL_COMBO  and y2 == F12_TWOEL_COMBO:
                                            resmini.coefficient.append(INTERM_X_F12+"_4")

                                       # resmini.coefficient.append(INTERM_X_F12)
                                       resmini.coefficient_idx.append(rest_f_idx+rest_v_idx)
                                       resmini.num_factor *= (2.0)
                                       append = True
                                   elif lst_of_completev_sort == lst_of_completev2_sort:
                                       sx.append(summed_idx([x.coefficient_idx[i], x.coefficient_idx[i2]], x.summation))

                                       resmini = x.copy_without([i, i2])
                                       if y == F12_TWOEL  and y2 == F12_TWOEL:
                                           resmini.coefficient.append(INTERM_X_F12)
                                       elif y == F12_TWOEL  and y2 == F12_TWOEL_COMBO:
                                           resmini.coefficient.append(INTERM_X_F12+"_2")
                                       elif y == F12_TWOEL_COMBO  and y2 == F12_TWOEL:
                                           resmini.coefficient.append(INTERM_X_F12+"_3")
                                       elif y == F12_TWOEL_COMBO  and y2 == F12_TWOEL_COMBO:
                                            resmini.coefficient.append(INTERM_X_F12+"_4")

                                       # resmini.coefficient.append(INTERM_X_F12)
                                       # tutaj sprawdzam czy posortowana jest pierwsza czy druga lista, zeby wiedziec
                                       # jak dodac wspolczynniki
                                       if lst_of_completev == lst_of_completev2_sort:
                                           # posortowana jest pierwsza wiec zmieniam indeksy dolne drugiego
                                           rest_v_idx = [rest_v_idx[1], rest_v_idx[0]]
                                       elif lst_of_completev_sort == lst_of_completev2:
                                           # posortowana jest druga wiec zmieniam indeksy dolne pierwszego
                                           rest_f_idx = [rest_f_idx[1], rest_f_idx[0]]
                                       resmini.coefficient_idx.append(rest_f_idx+rest_v_idx)
                                       resmini.num_factor *= (2.0)
                                       append = True

        if append == True:
            sx_all = []
            for l in sx:
                for ll in l:
                    sx_all.append(ll)
            for s in x.summation:
                if s not in sx_all:
                    resmini.summation.append(s)

            # print('intermedaiate X = ff* ff')
            # print('')
            # print(x)
            # print(resmini)
            # print('')

            result.append(resmini)
        else:
            result.append(x)


    # print('result z intermediate')
    # for x in result:
    #     print(x)

    return result


def identify_interm_Ft_f12(res):

#    print('ident interm Ft')
    # for x in res:
    #     print('innnt', x)

    result = []
    k = 0
    for x in res:
        #printt('xxxplusz', k, x )
        k += 1
        append = False
        append2 = False
        excluded1 = []
        resmini = deepcopy(x)
        pair = []
        pair_idx = []
        sign = []
        sx = []
        for i in range(0, len(x.coefficient)):
            y = x.coefficient[i]
            if y == F12_TWOEL or y == F12_TWOEL_COMBO:
                #printt(' i oto y', i, x, y)
                excluded1.append(i)
                n_of_occupied = 0
                lst_of_occupied = []
                rest_f_idx = []
                for idx in x.coefficient_idx[i]:
                    if idx in occupied:
                        n_of_occupied += 1
                        lst_of_occupied.append(idx)
                    else:
                        rest_f_idx.append(idx)
                lst_of_occupied_sort = deepcopy(lst_of_occupied)
                lst_of_occupied_sort.sort()
                #printt(lst_of_occupied, rest_f_idx)
                if n_of_occupied == 2:
                    excluded2 = []
                    for i2 in range(0, len(x.coefficient)):
                        if i2 not in excluded1:
                           y2 = x.coefficient[i2]
                           if y2 == F12_AMPLITUDE:
                               #printt(' i oto y2', i2, x, y2)
                               excluded2.append(i2)
                               lst_of_occupied2 = []
                               n_of_occupied2 = 0
                               rest_tf_idx = []
                               for idx in x.coefficient_idx[i2]:
                                   if idx in lst_of_occupied:
                                       lst_of_occupied2.append(idx)
                                       n_of_occupied2 += 1
                                   else:
                                       rest_tf_idx.append(idx)
                               #printt(lst_of_occupied2, rest_tf_idx)
                               lst_of_occupied2_sort = deepcopy(lst_of_occupied2)
                               lst_of_occupied2_sort.sort()
                               #printt('lsls', lst_of_occupied, lst_of_occupied2_sort)

                                #                              lst_of_occupied2.sort()
                               if n_of_occupied2 == 2:
                                   if lst_of_occupied == lst_of_occupied2:
                                       #printt('rowne listy')
                                       pair.append(deepcopy(i))
                                       #printt('pair1', pair)
                                       pair.append(deepcopy(i2))
                                       #printt('pair2', pair)
                                       # if (rest_f_idx[0] in virtual and rest_f_idx[1] in CABS):
                                       #     z = deepcopy(rest_f_idx[0])
                                       #     rest_f_idx[0] = deepcopy(rest_f_idx[1])
                                       #     rest_f_idx[1] = deepcopy(z)
                                       pair_idx.append(deepcopy(rest_f_idx)+deepcopy(rest_tf_idx))
                                       #printt('pair3tu', pair_idx)
                                       
                                       sx.append(summed_idx([x.coefficient_idx[i], x.coefficient_idx[i2]], x.summation))
                                       #printt('sx1', sx)
                                       sign.append(0)
                                       append = True
                                   elif lst_of_occupied_sort == lst_of_occupied2_sort :
                                       #printt('rowne listy sort')
                                       pair.append(deepcopy(i))
                                       #printt('pair1', pair)
                                       pair.append(deepcopy(i2))
                                       #printt('pair2', pair)
                                       # if (rest_f_idx[0] in virtual and rest_f_idx[1] in CABS):
                                       #     z = deepcopy(rest_f_idx[0])
                                       #     rest_f_idx[0] = deepcopy(rest_f_idx[1])
                                       #     rest_f_idx[1] = deepcopy(z)
                                       pair_idx.append(deepcopy(rest_f_idx) + deepcopy([rest_tf_idx[1], rest_tf_idx[0]]))
                                       #printt('pair3nietu', pair_idx)
                                       sx.append(summed_idx([x.coefficient_idx[i], x.coefficient_idx[i2]], x.summation))
                                       #printt('sx2', sx)
                                       # sign.append(1)

                                       append2 = True
                                   else:
                                       continue

        if append == True or append2 == True:
            sx_all = []
            for l in sx:
                for ll in l:
                    sx_all.append(ll)
                
            resmini = x.copy_without(pair)
            for z in range(0, len(pair_idx)):
                if append == True:
                    resmini.coefficient.append(INTERM_Ft_F12)
                elif append2 == True:
                    resmini.coefficient.append(INTERM_Ftt_F12)
                resmini.coefficient_idx.append(pair_idx[z])
#                resmini.num_factor *= (2.0)

            for s in x.summation:
                if s not in sx_all:
                    resmini.summation.append(s)

            # print('')
            # print('nowy x      ', x)
            # print('nowy resmini', resmini)
            # print('')
            result.append(resmini)
        else:
            result.append(x)


    return result

def identify_interm_Z_f12(res):

    result = []
    
    for x in res:
        # z += 1
        append = False
        excluded = []
        sx = []

        f12_indices = []
        obsx_indices = []
        for i, item in enumerate(x.coefficient):
            if item ==  F12_TWOEL or item == F12_TWOEL_COMBO :
                listf = x.coefficient_idx[i]
                if not any(elem in completev for elem in listf):
                    break
                else:
                    f12_indices.append(i)
            elif item == OBSERVABLE_X or item ==OBSERVABLE_AX:
                listx = x.coefficient_idx[i]
                if not all(elem in completev for elem in listx):
                    break
                else:
                    obsx_indices.append(i)

        if len(f12_indices)==2 and len(obsx_indices)==1:

            
            xidx = obsx_indices[0]
            listx = x.coefficient_idx[xidx]
            summin = set(listx)

            fidx1 = f12_indices[0]
            fidx2 = f12_indices[1]

            indicesf1 = [i for i, elem in enumerate(x.coefficient_idx[fidx1]) if elem in completev]
            if indicesf1[0] % 2 == 0:
                f1_even = True
                listf1_up = [x.coefficient_idx[fidx1][0], x.coefficient_idx[fidx1][2]]
                listf1_d = [x.coefficient_idx[fidx1][1], x.coefficient_idx[fidx1][3]]
            else:
                f1_even = False
                listf1_d = [x.coefficient_idx[fidx1][0], x.coefficient_idx[fidx1][2]]
                listf1_up = [x.coefficient_idx[fidx1][1], x.coefficient_idx[fidx1][3]]
            indicesf2 = [i for i, elem in enumerate(x.coefficient_idx[fidx2]) if elem in completev]
            if indicesf2[0] % 2 == 0:
                f2_even = True
                listf2_up = [x.coefficient_idx[fidx2][0], x.coefficient_idx[fidx2][2]]
                listf2_d = [x.coefficient_idx[fidx2][1], x.coefficient_idx[fidx2][3]]
            else:
                f2_even = False
                listf2_d = [x.coefficient_idx[fidx2][0], x.coefficient_idx[fidx2][2]]
                listf2_up = [x.coefficient_idx[fidx2][1], x.coefficient_idx[fidx2][3]]
            
            summin = list(set(listx) | set(listf1_up) |set(listf2_up))
            print(summin, x.summation)
            resmini = x.copy_without([xidx, fidx1, fidx2])                        
            resmini.coefficient.append(INTERM_Z_F12)
            z_idx = []

            if listx[0] in listf1_up and listx[1] in listf2_up:
                idxa = listf1_up.index(listx[0])
                idxb = listf2_up.index(listx[1])
                ff1 = deepcopy(listf1_d)
                ff2 = deepcopy(listf2_d)
            elif listx[0] in listf2_up and listx[1] in listf1_up:
                idxa = listf2_up.index(listx[0])
                idxb = listf1_up.index(listx[1])
                ff1 = deepcopy(listf2_d)
                ff2 = deepcopy(listf1_d)

            else:
                print('find Z something weird')

            if idxa == 0:
                z_idx.append(ff1[1])
                z_idx.append(ff1[0])
            elif idxa == 1:
                z_idx.append(ff1[0])
                z_idx.append(ff1[1])
            if idxb == 0:

                z_idx.append(ff2[0])
                z_idx.append(ff2[1])
            elif idxb == 1:
                z_idx.append(ff2[1])
                z_idx.append(ff2[0])

            resmini.coefficient_idx.append(z_idx)
            resmini.summation = list(set(x.summation) - set(summin))
            append = True


        if append:
            result.append(resmini)
        else:
            result.append(x)
            
    return result

def identify_interm_B_or_Z_f12(res, operat):

    if operat == OBSERVABLE_X:
        name = INTERM_Z_F12
    else:
        name = INTERM_B_F12

    result = []
    # z = 0
    for x in res:
        # z += 1
        append = False
        excluded = []
        sx = []
        
        for i in range(0, len(x.coefficient)):
            if (i not in excluded):
                y = x.coefficient[i]
                if y == F12_TWOEL or y == F12_TWOEL_COMBO :
                    # print('tak', y)
                    excluded.append(i)
                    n_of_completev = 0
                    lst_of_completev = []
                    rest_f_idx = []
                    for idx in x.coefficient_idx[i]:
                        if idx in completev:
                            n_of_completev += 1
                            lst_of_completev.append(idx)
                        else:
                            rest_f_idx.append(idx)
                    lst_of_completev.sort()
                    if n_of_completev == 2:
                        # print('n of complete', n_of_completev)
                        for i2 in range(0, len(x.coefficient)):
                            if (i2 not in excluded):
                                y2 = x.coefficient[i2]
                                if operat == OBSERVABLE_X:
                                    condition =  y2 == OBSERVABLE_X
                                else:
                                    condition =  y2 == FOCK_MATRIX or y2 == FOCK_MATRIX_TRANS
                                if condition:
                                    # if y2 == FOCK_MATRIX or y2 == FOCK_MATRIX_TRANS:
                                    # print('tak fock')
                                    excluded.append(i2)
                                    lst_of_completev2 = []
                                    n_of_completev2 = 0
                                    rest_2_idx = []
                                    for idx in x.coefficient_idx[i2]:
                                        if idx in completev:
                                            n_of_completev2 += 1
                                            lst_of_completev2.append(idx)
                                        else:
                                            rest_2_idx.append(idx)
                                    lst_of_completev2.sort()
                                    if n_of_completev2 ==2:
                                        # print('n of complete', n_of_completev2)                                                                
                                        if ((lst_of_completev2[0] in lst_of_completev) or (lst_of_completev2[1] in lst_of_completev)):
                                            for i3 in range(0, len(x.coefficient)):
                                                if (i3 not in excluded):
                                                    y3 = x.coefficient[i3]
                                                    if y3 == F12_TWOEL or y3 == F12_TWOEL_COMBO:
                                                        excluded.append(i3)
                                                        lst_of_completev3 = []
                                                        n_of_completev3 = 0
                                                        rest_3_idx = []
                                                        for idx in x.coefficient_idx[i3]:
                                                            if idx in completev:
                                                                n_of_completev3 += 1
                                                                lst_of_completev3.append(idx)
                                                            else:
                                                                rest_3_idx.append(idx)
                                                        lst_of_completev3.sort()
                                                        if n_of_completev3 == 2:
                                                            if (((lst_of_completev3[0] in lst_of_completev2)\
                                                                 and (lst_of_completev3[1] in lst_of_completev)) or \
                                                                    ((lst_of_completev3[1] in lst_of_completev2)\
                                                                     and (lst_of_completev3[0] in lst_of_completev))):
                                                                sx.append(summed_idx([x.coefficient_idx[i], \
                                                                                      x.coefficient_idx[i2], \
                                                                                      x.coefficient_idx[i3]], x.summation))
                                                                resmini = x.copy_without([i, i2, i3])
                                                                # print('yy', y, y2)
                                                                if y == F12_TWOEL  and y3 == F12_TWOEL:
                                                                    resmini.coefficient.append(name)
                                                                elif y == F12_TWOEL  and y3 == F12_TWOEL_COMBO:
                                                                    resmini.coefficient.append(name+"_2")
                                                                elif y == F12_TWOEL_COMBO  and y3 == F12_TWOEL:
                                                                    resmini.coefficient.append(name+"_3")
                                                                elif y == F12_TWOEL_COMBO  and y3 == F12_TWOEL_COMBO:
                                                                    resmini.coefficient.append(name+"_4")
                                                                    

                                                                # resmini.coefficient.append(INTERM_B_F12)
                                                                resmini.coefficient_idx.append(rest_f_idx+rest_3_idx)
                                                                append = True
        if append == True:
            sx_all = []
            for l in sx:
                for ll in l:
                    sx_all.append(ll)
            for s in x.summation:
                if s not in sx_all:
                    resmini.summation.append(s)

            # print('intermedaiate B = ff*v* ff')
            # print('')
            # print(x)
            # print(resmini)
            # print('')
            result.append(resmini)
        else:
            result.append(x)


    # print('result z intermediate')
    # for x in result:
    #     print(x)

    return result


def identify_interm_P_f12(res):

    result = []
    for x in res:
        append = False
        excluded = []
        sx = []
        for i in range(0, len(x.coefficient)):
            if (i not in excluded):
                y = x.coefficient[i]
                if y == F12_TWOEL or y == F12_TWOEL_COMBO:
                    excluded.append(i)
                    n_of_completev = 0
                    lst_of_completev = []
                    rest_f_idx = []
                    for idx in x.coefficient_idx[i]:
                        if idx in completev:
                            n_of_completev += 1
                            lst_of_completev.append(idx)
                        else:
                            rest_f_idx.append(idx)
                    lst_of_completev.sort()
                    if n_of_completev == 2:
                        for i2 in range(0, len(x.coefficient)):
                            if (i2 not in excluded):
                                y2 = x.coefficient[i2]
                                if y2 == TWOEL_INT or y2 == TWOEL_INT_TRANS or y2 == TWOEL_INT_COMBO:
                                    # print('mam v', y2, x.coefficient_idx[i2])
                                    excluded.append(i2)
                                    lst_of_completev2 = []
                                    n_of_completev2 = 0
                                    rest_2_idx = []
                                    for idx in x.coefficient_idx[i2]:
                                        if idx in completev:
                                            n_of_completev2 += 1
                                            lst_of_completev2.append(idx)
                                        else:
                                            rest_2_idx.append(idx)
                                    lst_of_completev2.sort()
                                    if n_of_completev2 == 4:
                                        if ((lst_of_completev[0] in lst_of_completev2) and (lst_of_completev[1] in lst_of_completev2)):
                                            for i3 in range(0, len(x.coefficient)):
                                                if (i3 not in excluded):
                                                    y3 = x.coefficient[i3]
                                                    if y3 == F12_TWOEL or y3 == F12_TWOEL_COMBO:
                                                        # print('i mam kolejne F', y3, x.coefficient_idx[i3])
                                                        excluded.append(i3)
                                                        lst_of_completev3 = []
                                                        n_of_completev3 = 0
                                                        rest_3_idx = []
                                                        for idx in x.coefficient_idx[i3]:
                                                            if idx in completev:
                                                                n_of_completev3 += 1
                                                                lst_of_completev3.append(idx)
                                                            else:
                                                                rest_3_idx.append(idx)
                                                        lst_of_completev3.sort()
                                                        # print('n_of_competev', n_of_completev3)
                                                        if n_of_completev3 == 2:
                                                            if (((lst_of_completev3[0] in lst_of_completev2)and (lst_of_completev3[1] in lst_of_completev2)) and \
                                                                    ((lst_of_completev3[0] not in lst_of_completev ) \
                                                                         and (lst_of_completev3[1] not in lst_of_completev))):
                                                                sx.append(summed_idx([x.coefficient_idx[i], x.coefficient_idx[i2], x.coefficient_idx[i3]], x.summation))

                                                                resmini = x.copy_without([i, i2, i3])
                                                                if y == F12_TWOEL  and y2 in TWOEL_INT and y3 ==F12_TWOEL:
                                                                    resmini.coefficient.append(INTERM_P_F12)
                                                                elif y == F12_TWOEL  and y2 in TWOEL_INT and y3 ==F12_TWOEL_COMBO:
                                                                    resmini.coefficient.append(INTERM_P_F12+"_2")
                                                                elif y == F12_TWOEL  and y2 == TWOEL_INT_COMBO and y3 ==F12_TWOEL:
                                                                    resmini.coefficient.append(INTERM_P_F12+"_3")
                                                                elif y == F12_TWOEL_COMOBO  and y2 in TWOEL_INT and y3 ==F12_TWOEL:
                                                                    resmini.coefficient.append(INTERM_P_F12+"_4")
                                                                elif y == F12_TWOEL  and y2 == TWOEL_INT_COMBO and y3 ==F12_TWOEL_COMBO:
                                                                    resmini.coefficient.append(INTERM_P_F12+"_5")
                                                                elif y == F12_TWOEL_COMBO  and y2 in TWOEL_INT and y3 ==F12_TWOEL_COMBO:
                                                                    resmini.coefficient.append(INTERM_P_F12+"_6")
                                                                elif y == F12_TWOEL_COMBO  and y2 == TWOEL_INT_COMBO and y3 ==F12_TWOEL:
                                                                    resmini.coefficient.append(INTERM_P_F12+"_7")
                                                                elif y == F12_TWOEL_COMBO  and y2 == TWOEL_INT_COMBO and y3 ==F12_TWOEL_COMBO:
                                                                    resmini.coefficient.append(INTERM_P_F12+"_8")

                                                                resmini.coefficient_idx.append(rest_f_idx+rest_3_idx)
                                                                resmini.num_factor *= 4.0
                                                                append = True
        if append == True:
            sx_all = []
            for l in sx:
                for ll in l:
                    sx_all.append(ll)
            for s in x.summation:
                if s not in sx_all:
                    resmini.summation.append(s)

            # print('intermedaiate P = ff v* ff')
            # print('')
            # print(x)
            # print(resmini)
            # print('')

            result.append(resmini)
        else:
            result.append(x)

    # print('result z intermediate')
    # for x in result:
    #     print(x)

    return result


def summed_idx(idx_lst_lst, sum_idx):

    # print(idx_lst_lst)
    # print(sum_idx)
    # print('hueee')
    summed_idx = []
    idx_lst = []
    for x in idx_lst_lst:
        for y in x:
            idx_lst.append(y)
    # print('idx_lst', idx_lst)

        
    idx_dict = {}
    for x in idx_lst:
        if x in idx_dict.keys():
            idx_dict[x]+=1
        else:
            idx_dict[x] = 1
    # print('idx_dict', idx_dict)
    # print('fixed', fixed)
    for x in idx_dict.keys():
        if idx_dict[x] == 2:
            #if x not in fixed:
            if x in sum_idx:
                summed_idx.append(x)
                    
        
    # for x in idx1:
    #     if x in idx2:
    #         if x not in fixed:
    #             if x in sum_idx:
                    # summed_idx.append(x)
                    
    return summed_idx



# def cabssplit(res):
#     result = arithmetic_string()

#     print('llllaaaa')
#     for x in res:
#         print('xxx', x)
#         result = result + cabssplit_ugg(x)

#     return result
        


# def cabssplit_ugg(res):
    
#     result = []
#     print('ressssss', res)
#     complidx = False
#     for p in res.summation:
#         if p in completev:
#             complidx = True
#             break
#     print('complidx', complidx)
#     if complidx:
#         s = res.cabssplit()
#         return cabssplit_ugg(s)

    



def find_permutations(result):
    # print('find permutations')
    # determine fixed
    fx_virt = []
    fx_occ = []
    for x in result:
        for y in x.coefficient_idx:
            for z in y:
                if z not in x.summation:
                    if z in virtual:
                        if z not in fx_virt:
                            fx_virt.append(z)
                    elif z in occupied:
                        if z not in fx_occ:
                            fx_occ.append(z)


    if len(fx_virt) == 0 and len(fx_occ) == 4:
        fx_a = ['i', 'j']
        fx_b = ['k', 'l']
    elif len(fx_virt) == 2 and len(fx_occ) ==2:
        fx_a = fx_virt
        fx_b = fx_occ
    elif len(fx_virt) == 1 and len(fx_occ) ==1:
        fx_a = fx_virt
        fx_b = fx_occ
    # print(fx_a)
    # print(fx_b)

    fx_virt.sort()
    fx_occ.sort()
    # print('zaczynam sprawdzanie')
    for i in range(0, len(result)):
        # print
        #tempi1 = Pab
        #tempi2 = PabPij
        #tempi3 = Pij
        tempi1 = deepcopy(result[i])
        tempi3 = deepcopy(result[i])
        # print(result[i])
        # print(fx_a)
        # print('tempi1')
        tempi1.multisubst([fx_a[0]], ['z'])
        tempi1.multisubst([fx_a[1]], [fx_a[0]])
        tempi1.multisubst(['z'], [fx_a[1]])
        tempi1.num_factor == 1.0

        # print('tempi2')
        tempi2 = deepcopy(tempi1)
        tempi2.multisubst([fx_b[0]], ['z'])
        tempi2.multisubst([fx_b[1]], [fx_b[0]])
        tempi2.multisubst(['z'], [fx_b[1]])
        tempi2.num_factor == 1.0

        # print('tempi3')
        tempi3.multisubst([fx_b[0]], ['z'])
        tempi3.multisubst([fx_b[1]], [fx_b[0]])
        tempi3.multisubst(['z'], [fx_b[1]])
        tempi3.num_factor == 1.0

        # print('sprawdzam', result[i])
        # print(tempi1)
        # print(tempi2)
        # print(tempi3)

        tempi1.standarize()
        tempi2.standarize()
        tempi3.standarize()

        # print('po standarize')
        # print(tempi1)
        # print(tempi2)
        # print(tempi3)
        # print('')
        for j in range(i+1, len(result)):
            tempj1 = deepcopy(result[j])
            tempj1.num_factor = 1.0

            # print('z', result[j])
            # print('tempj1', tempj1)
            # print('tempi1', tempi1)
            # print('tempi2', tempi2)
            # print('tempi3', tempi3)
            found = False
            if tempi1 == tempj1:
                # print('tak1')
                found = True
                if (result[i].num_factor < 0.0 and result[j].num_factor < 0.0) \
                        or (result[j].num_factor > 0.0 and result[j].num_factor > 0.0):
                    # print('a')
                    result[i].num_factor = 0.0
                    result[j].coefficient.append('Pp')
                    result[j].coefficient_idx.append(fx_a)
                elif result[i].num_factor < 0.0 and result[j].num_factor > 0.0:
                    # print('b')
                    result[i].num_factor = 0.0
                    result[j].coefficient.append('Pm')
                    result[j].coefficient_idx.append(fx_a)
                elif result[i].num_factor > 0.0 and result[j].num_factor < 0.0:
                    # print('c')
                    result[j].num_factor = 0.0
                    result[i].coefficient.append('Pm')
                    result[i].coefficient_idx.append(fx_a)
            if found:
                continue
            if tempi3 == tempj1:
                # print('tak4')
                found = True
                if (result[i].num_factor < 0.0 and result[j].num_factor < 0.0) \
                        or (result[j].num_factor > 0.0 and result[j].num_factor > 0.0):
                    # print('a')
                    result[i].num_factor = 0.0
                    result[j].coefficient.append('Pp')
                    result[j].coefficient_idx.append(fx_b)
                elif result[i].num_factor < 0.0 and result[j].num_factor > 0.0:
                    # print('b')
                    result[i].num_factor = 0.0
                    result[j].coefficient.append('Pm')
                    result[j].coefficient_idx.append(fx_b)
                elif result[i].num_factor > 0.0 and result[j].num_factor < 0.0:
                    # print('c')
                    result[j].num_factor = 0.0
                    result[i].coefficient.append('Pm')
                    result[i].coefficient_idx.append(fx_b)

            if found:
                continue
            if tempi2 == tempj1:
                # print('tak3')                
                found = True
                if (result[i].num_factor < 0.0 and result[j].num_factor < 0.0) \
                        or (result[j].num_factor > 0.0 and result[j].num_factor > 0.0):
                    # print('a')
                    result[i].num_factor = 0.0
                    result[j].coefficient.append('Pp')
                    result[j].coefficient_idx.append(fx_a+fx_b)
                    # result[j].coefficient.append('Pp')
                    # result[j].coefficient_idx.append(fx_b)
                elif result[i].num_factor < 0.0 and result[j].num_factor > 0.0:
                    # print('b')
                    result[i].num_factor = 0.0
                    result[j].coefficient.append('Pm')
                    result[j].coefficient_idx.append(fx_a+fx_b)
                    # result[j].coefficient.append('Pm')
                    # result[j].coefficient_idx.append(fx_b)
                elif result[i].num_factor > 0.0 and result[j].num_factor < 0.0:
                    # print('c')
                    result[j].num_factor = 0.0
                    result[i].coefficient.append('Pm')
                    result[i].coefficient_idx.append(fx_a+fx_b)
                    # result[i].coefficient.append('Pm')
                    # result[i].coefficient_idx.append(fx_b)


    result.cleanup()
    # print('to jest ostateczny wynik')

    # for x in result:
    #     print(x)
    return result
            


def preprep_for_fortran(rint, no_fx=None):

    #
    # Optimization is in fortran_code.py
    #
    rsimp = simplify(rint, no_fx)
    rint = deepcopy(rsimp)
    rint.exec_delta()
    rsimp = simplify(rint, no_fx)
#    print('x', no_fx)
#    for x in rsimp:
#        print(x)

    for x in rsimp:
        for y in range(0, len(x.coefficient)):
            if x.coefficient[y] == FOCK_MATRIX:
                if x.coefficient_idx[y][0] in virtual and x.coefficient_idx[y][1] in occupied:
                    x.num_factor = 0.0
                elif x.coefficient_idx[y][0] in occupied and x.coefficient_idx[y][1] in virtual:
                    x.num_factor = 0.0

    rsimp.cleanup()

    return rsimp

def evaluate_s_operator():
    """ All equations are in a/doc/cc/S_operators

    S_1 = T_1 + S_1[2] + S_1[3] = 
    T_1 +  P_1([T1c, T2]) + P_1(0.5[[T1c, T1], T1] + 0.5[[T2c, T2], T1]
    + 0.5[[T2c, T1], T2])

    S_2 = T_2 + S_2[3 ] = T_2 + P_2(0.5[[T1c, T2], T1]
    + [T1c, T3] + 0.5[[T1c, T1], T2])

    S_3 = T_3 + S_3[3] = T_3 + P_3(0.5[[T1c, T2], T2]
    + 0.5[[T1c, T1], T3] + [[T1c, T3], T1])
    """
    s11 = arithmetic_string(t1)
    
    s12 = evaluate(t1c, t2)
    
    s13 = evaluate(t1c, t1, t1).scale(0.5) + evaluate(t2c, t2, t1).scale(0.5) \
        + evaluate(t2c, t1, t2).scale(0.5)
    
    ss1 = s11 + s12 + s13
    sint1  = ss1.integrate(bra = ['a','i']).scale(0.5)
    sint1.exec_delta()
    ssimp1 = simplify(sint1)
    ssimp1.cleanup()
    for x in ssimp1:
        x.optimize()
        print(x)

    function_t_to_s(ssimp1, 1)

    s21 = arithmetic_string(t2)
    
    s22 = evaluate(t1c, t3)

    s23 = evaluate(t1c, t2, t1) \
        + evaluate(t1c, t1, t2)

    ss20 = s21 + s23
    ss2 = simplify(ss20)

    rint1 = ss2.integrate(bra = ['a', 'i', 'b', 'j'])
    rint2 = ss2.integrate(bra = ['a', 'j', 'b', 'i'])
    sint2  = rint1.scale(1./3.) + rint2.scale(1./6.)
    sint2.exec_delta()
    
    ssimp2 = simplify(sint2)
    ssimp2.cleanup()
    for x in ssimp2:
        print(x)
    for x in ssimp2:
        x.optimize()

    function_t_to_s(ssimp2, 2)

    s31 = arithmetic_string(t3)
    
    s33 = evaluate(t1c, t1, t3).scale(0.5) + evaluate(t1c, t2, t2).scale(0.5)\
        + evaluate(t1c, t3, t1)

    ss30 = s31 + s33

    for x in ss30:
        print(x)
    print('-----------')
    ss3 = simplify(ss30)
    ss3.cleanup()

    print('teraz')
    ssimp3 = nazwa_funkcji_ktory_wybiera_kawalek_wzoru(ss3) 

    for x in ssimp3:
        print(x)
    for x in ssimp3:
        x.optimize()

    function_t_to_s(ssimp3, 3)

def generate_s_f12():
    """All equations are in a/doc/cc/S_operators
    
    CCSD
    MBPT order
    S1 = T1                 # 2                                                                                                     
         + P1(T1c,T2)       # 3                                                                                       
         + P1(T2c,T1,T2)    # 4                                                                                               
         + 0.5P1(T1c,T1,T1) # 6                                                                       
    S2 = T2                 # 1                                                                                                 
         + 0.5P2(T2c,T2,T2) # 3     
         + P2(T1c,T1,T2)    # 5                                    
         """
    s1_ccsd = []
    s1_ccsd.append(arithmetic_string(t1))             #2 
    s1_ccsd.append(evaluate(t1c,t2))                  #3         
    s1_ccsd.append(evaluate(t1c,t2fa))                  #3
    #    s1_ccsd.append(evaluate(t2c,t1,t2))               #4     
    
    for x in s1_ccsd:
        for y in x:
            print(y)

    s2_ccsd = []
    s2_ccsd.append(arithmetic_string(t2))             #1                                                                                                             
    s2_ccsd.append(evaluate(t2c,t2,t2).scale(0.5))    #3   
    s2_ccsd.append(evaluate(t2c,t2,t2fa).scale(0.5))    #3
    s2_ccsd.append(evaluate(t2c,t2fa,t2).scale(0.5))    #3
    t2fac = deepcopy(t2fa)
    t2fac.transpose()
    print('a')
    print(t2fac)

    s2_ccsd.append(evaluate(t2fac,t2,t2).scale(0.5))    #3
    s2_ccsd.append(evaluate(t2fac,t2fa,t2).scale(0.5))    #3
    s2_ccsd.append(evaluate(t2fac,t2,t2fa).scale(0.5))    #3
    s2_ccsd.append(evaluate(t2fac,t2fa,t2fa).scale(0.5))    #3
    s2_ccsd.append(evaluate(t2c,t2fa,t2fa).scale(0.5))    #3
    
    for x in s2_ccsd:
        for y in x:
            y.excitation_rank()
#            print(y, y.excitation_rank())

    return s2_ccsd


def generate_s_operators():
    """ All equations are in a/doc/cc/S_operators

    CCSD
                            MBPT order
    S1 = T1                 # 2
         + P1(T1c,T2)       # 3
         + P1(T2c,T1,T2)    # 4
         + 0.5P1(T1c,T1,T1) # 6

    S2 = T2                 # 1
         + 0.5P2(T2c,T2,T2) # 3
         + P2(T1c,T1,T2)    # 5   

    S3 = 0.5P3(T1,T2,T2)    # 4

    CC3
                            MBPT order
    S1 = T1                 # 2
         + P1(T1c,T2)       # 3
         + P1(T2c,T3)
         + P1(T2c,T1,T2)    # 4
         + 0.5P1(T3c,T2,T2)
         + 0.5P1(T3,T2c,T2)
         + 0.5P1(T1c,T1,T1) # 6
         + 0.5P1(T3c,T1,T3) 
         + 0.5P1(T3,T1c,T3) 
         + 0.5P1(T3c,T3,T1) 

    S2 = T2                 # 1
         + 0.5P2(T2c,T2,T2) # 3
         + P2(T1c,T3)       # 4
         + P2(T1c,T1,T2)    # 5   
         + 0.5P2(T2c,T1,T3)   
         + 0.5P2(T2c,T3,T1)   
         + 0.5P2(T3c,T2,T3)   
         + 0.5P2(T3c,T3,T2)   

    S3 = T3                 # 2  
         + 0.5P3(T1c,T2,T2) # 4
         + P3(T2c,T2,T3) 
         + 0.5P3(T3c,T3,T3) 
         + P3(T1c,T1,T3)    # 6

    """

    f = open('s_functions.f90', 'w')

    s1_ccsd = []
    s1_ccsd.append(arithmetic_string(t1))             #2
    s1_ccsd.append(evaluate(t1c,t2))                  #3
    s1_ccsd.append(evaluate(t2c,t1,t2))               #4
    # s1_ccsd.append(evaluate(t1c,t1,t1).scale(0.5))


    s2_ccsd = []
    s2_ccsd.append(arithmetic_string(t2))             #1
    s2_ccsd.append(evaluate(t2c,t2,t2).scale(0.5))    #3
    # s2_ccsd.append(evaluate(t1c,t1,t2))
    
    # s3_ccsd = []                                      #4
    # s3_ccsd.append(evaluate(t1c,t2,t2).scale(0.5))

    s1_cc3 = []
    s1_cc3.append(arithmetic_string(t1))               #2
    s1_cc3.append(evaluate(t1c,t2) + evaluate(t2c,t3)) #3
    s1_cc3.append(evaluate(t2c,t1,t2) + evaluate(t3c,t2,t2).scale(0.5) \
                      + evaluate(t3,t2c,t2).scale(0.5))      #4
    # s1_cc3.append(evaluate(t1c,t1,t1).scale(0.5)     #6
    # s1_cc3.append(evaluate(t3c,t1,t3).scale(0.5))   
    # s1_cc3.append(evaluate(t3,t1c,t3).scale(0.5))   
    # s1_cc3.append(evaluate(t3c,t3,t1).scale(0.5))   

    s2_cc3 = []
    s2_cc3.append(arithmetic_string(t2))          #1
    s2_cc3.append(evaluate(t2c,t2,t2).scale(0.5)) #3
    s2_cc3.append(evaluate(t1c,t3))               #4
    # s2_cc3.append(evaluate(t1c,t1,t2))            #5  
    # s2_cc3.append(evaluate(t2c,t1,t3).scale(0.5)) #5  
    # s2_cc3.append(evaluate(t2c,t3,t1).scale(0.5)) #5  
    # s2_cc3.append(evaluate(t3c,t2,t3).scale(0.5)) #5  
    # s2_cc3.append(evaluate(t3c,t3,t2).scale(0.5)) #5  
    
    # s3_cc3 = []
    # s3_cc3.append(arithmetic_string(t3))          #2
    # s3_cc3.append(evaluate(t1c,t2,t2).scale(0.5)) #4
    # s3_cc3.append(evaluate(t2c,t2,t3))            #4
    # s3_cc3.append(evaluate(t3c,t3,t3).scale(0.5)) #4
    # s3_cc3.append(evaluate(t1c,t1,t3))            #6


    dc_s1_ccsd = deepcopy(s1_ccsd)
    for i in range(0, len(s1_ccsd)):
        j = i + 2
        s1_ccsd[i] = int_and_simp(dc_s1_ccsd[i], 1)
        # print('')
        # for x in s1_ccsd[i]:
        #     print(x)
        function_t_to_s(s1_ccsd[i], 1, 's1_ccsd_{j}'.format(j=j), f)
#    sys.exit(0)
    dc_s2_ccsd = deepcopy(s2_ccsd)
    for i in range(0, len(s2_ccsd)):
        if(i==0):
            j = i + 1
        else:
            j = i + 2        
        print('j', j)
        for k in s2_ccsd[i]:
            print(k)
#        s2_ccsd[i] = int_and_simp(dc_s2_ccsd[i], 2)

#        function_t_to_s(s2_ccsd[i], 2, 's2_ccsd_{j}'.format(j=j), f)

    # dc_s1_cc3 = deepcopy(s1_cc3)
    # for i in range(0, len(s1_cc3)):
    #     j = i + 2
    #     s1_cc3[i] = int_and_simp(dc_s1_cc3[i], 1)
    #     for x in s1_cc3[i]:
    #         s = []
    #         for y in x.summation:
    #             s.append(y)
    #         for y in x.coefficient_idx:
    #             for z in y:
    #                 if z not in s:
    #                     s.append(z)
    #         sv = []
    #         so = []
    #         for y in s:
    #             if y in occupied:
    #                 so.append(y)
    #             elif y in virtual:
    #                 sv.append(y)
    #         print(x)
    #         if len(x.coefficient) > 2:
    #             svi, svo = cost_after_intermediates(deepcopy(x))
    #             print(x, '         ', sv, so, len(sv), len(so),svi, svo, len(svi), len(svo))
    #         else:
    #             print(x, '         ', sv, so, len(sv), len(so))

    #     function_t_to_s(s1_cc3[i], 1, 's1_cc3_{j}'.format(j=j), f)

    # dc_s2_cc3 = deepcopy(s2_cc3)
    # for i in range(0, len(s2_cc3)):
    #     if(i==0):
    #         j = i + 1
    #     else:
    #         j = i + 2
    #     s2_cc3[i] = int_and_simp(dc_s2_cc3[i], 2)
    #     for x in s2_cc3[i]:
    #         s = []
    #         for y in x.summation:
    #             s.append(y)
    #         for y in x.coefficient_idx:
    #             for z in y:
    #                 if z not in s:
    #                     s.append(z)
    #         sv = []
    #         so = []
    #         for y in s:
    #             if y in occupied:
    #                 so.append(y)
    #             elif y in virtual:
    #                 sv.append(y)
    #         print(x)
    #         if len(x.coefficient) > 2:
    #             svi, svo = cost_after_intermediates(deepcopy(x))
    #             print(x, '         ', sv, so, len(sv), len(so),svi, svo, len(svi), len(svo))
    #         else:
    #             print(x, '         ', sv, so, len(sv), len(so))

    #     function_t_to_s(s2_cc3[i], 2, 's2_cc3_{j}'.format(j=j), f)
    sys.exit(0)
    f.close()

# def evaluate_s_f12_operator_separate_functions():

#     ev_s21 = arithmetic_string(t2) + arithmetic_string(t2fa)

#     t2fac = deepcopy(t2fa)
#     t2fac.transpose()
#     # print(t2fa)
#     # print(t2fac)
#     # print(t2)
#     # print(t2c)
#     # sys.exit(0)

#     ev_s22 = evaluate(t2c, t2, t2).scale(0.5) + evaluate(t2c, t2, t2fa).scale(0.5) + evaluate(t2c, t2fa, t2).scale(0.5) \
#         + evaluate(t2c, t2fa, t2fa).scale(0.5) + evaluate(t2fac, t2fa, t2).scale(0.5) + evaluate(t2fac, t2, t2fa).scale(0.5)\
#         + evaluate(t2fac, t2fa, t2fa).scale(0.5) + evaluate(t2fac, t2, t2).scale(0.5)

#     print('ev_s22 przed calkowaniem')
#     k = 1
#     for x in ev_s22:
#         print(k, x)
#         k += 1
#     print('')

#     s21 = int_and_simp(ev_s21, 2)

#     s22 = int_and_simp(ev_s22, 2)


#     for x in range(0, len(s21)):
#         for y in s21[x].coefficient:
#             at_least_one = True
#             if y == F12_TWOEL:
#                 at_least_one = False
#                 for i in s21[x].coefficient_idx:
#                     if i in completev:
#                         at_least_one = True
#                         break
#             if at_least_one == False:
#                 s21[x].num_factor == 0

#     s21.cleanup()

#     print('no to zaczynam s22')
#     print('len', len(s22))
#     for x in range(0, len(s22)):
#         for j in range(0, len(s22[x].coefficient)):
#             y = s22[x].coefficient[j]
#             at_least_one = True
#             if y == F12_TWOEL:
#                 print('tak jest f12_twoel')
#                 print(s22[x])
#                 at_least_one = False
#                 print('at_least_one', at_least_one)
#                 for i in s22[x].coefficient_idx[j]:
#                     print(i)
#                     if i in completev:
#                         at_least_one = True
#                         print('at_least_one', at_least_one)
#                         break
#             if at_least_one == False:
#                 print('tak, falsz')
#                 s22[x].num_factor = 0

#     s22.cleanup()
#     print('len2', len(s22))

                

#     print('s22 wynik')
#     for x in s22:
#         print(x)


#     print('TERAZ FAKTORYZUJE')
#     k = 1
#     kkk = 0
#     biglevel_super = []
#     bigparameters_super = []
#     basket_nointerm = []
#     idx_interm = []
#     idx_noninterm = []
#     idxidx = -1
#     for x in s22:
#         biglevel_best, bigparameters_best, haveinterm = factorize(x)
#         idxidx += 1
#         if haveinterm == True:
#             idx_interm.append(idxidx)
#             biglevel_super.append(biglevel_best)
#             bigparameters_super.append(bigparameters_best)
#             print(k, x, biglevel_best)
#             k += 1
#             s = 1
#             if len(biglevel_best) == 0:
#                 kkk += 1
#             for f in range(0, len(bigparameters_best)):
#                 z = bigparameters_best[f]
#                 zz = biglevel_best[f]
#             s+= 1
#         else:
#             basket_nointerm.append(x)
#             idx_noninterm.append(idxidx)


#     sys.exit(0)


def evaluate_s_operator_separate_functions():
    """ All equations are in a/doc/cc/S_operators

    S_1 = T_1 + S_1[2](3) = 
    T_1 +  P_1([T1c, T2]) +  P_1([T2c, T3])

    S_2 = T_2 + S_2[3](3) = T_2 + 0.5 P_2([[T2c, T2], T2]

    S_3 = T_3

    """

#-----------------------------TO BYL KOD DLA POLNONA DO SPRAWDZENIA S1
    # ev_s2 = evaluate(s1, t2c, t2)
    # s2 = int_and_simp(ev_s2, 1, 1)
    # print('to wynik tego kom')
    # for x in s2:
    #     print(x)
    # sys.exit(0)
    

    
    # ev_s2 = evaluate(t2c, t2, t2).scale(0.5)
    # s2 = int_and_simp(ev_s2, 2, 1)
    # print('to jest s2 pierwsza czesc')
    # for x in range(0, len(s2)):
    #     s2[x] = s2[x].fromright(eaibj)
    #     print(s2[x])
    # ars = arithmetic_string()
    # for x in s2:
    #     print('x', x, t2c, t2c)
    #     rr = evaluate(x, t2c)
    #     for x in rr:
    #         print('dfdf', x)
    #     print('rra', len(rr))
    #     for y in rr:
    #         ars = ars + deepcopy(arithmetic_string(y))

    # for x in ars:
    #       print(x)
    # s11 = int_and_simp(ars, 1)
    # print('wynik-ostateczny')
    # for x in s11:
    #     print(x)
            
    # sys.exit(0)


    # ev = evaluate(t2c, t2, t2, t2, t2)
    # for x in ev:
    #     print(x)
    # s11 = int_and_simp(ev, 1)
    # print('posimp')
    # for x in s11:
    #     print(x)

    # sys.exit(0)
#------------------------------------------------------    

    f = open('s_functions.f90', 'w')

    ev_s11 = arithmetic_string(t1)
    
    ev_s12a = evaluate(t1c, t2)
    ev_s12b = evaluate(t2c, t3)

    s11 = int_and_simp(ev_s11, 1)

    s12a = int_and_simp(ev_s12a, 1)
    s12b = int_and_simp(ev_s12b, 1)

    function_t_to_s(s11, 1, 's11', f)

    function_t_to_s(s12a, 1, 's12a', f)
    function_t_to_s(s12b, 1, 's12b', f)


    ev_s21 = arithmetic_string(t2)
    
    ev_s23 = evaluate(t2c, t2, t2).scale(0.5)

    #--------------------------------------dla POLNONA
    # ev_s23p = deepcopy(ev_s23)
    # s23p = int_and_simp(ev_s23p, 2)
    # print('')
    # for x in s23p:
    #     print(x)
    # sys.exit()


    # dla Polnona S_2^(5) = 1/4 P_2 [[[[T_2,T_2*],T_2*],T_2],T_2]
    
    # ev_s25 = evaluate(s2,t2c,t2).scale(-0.5)
    # s25 = int_and_simp(ev_s25, 2)
    # for x in ev_s25:
    #     print(x)

    # print('')
    # print('wynik')
    # for x in s25:
    #     print(x)

    # sys.exit(0)
    #--------------------------------------KONIEC dla POLNONA

 
    s21 = int_and_simp(ev_s21, 2)

    s23 = int_and_simp(ev_s23, 2)


    function_t_to_s(s21, 2, 's21', f)

    function_t_to_s(s23, 2, 's23', f)


    ev_s31 = arithmetic_string(t3)
    print('ev_s31')
    for x in ev_s31:
        print(x)
    print('')
    s31 = nazwa_funkcji_ktory_wybiera_kawalek_wzoru(ev_s31) 
    print('s31')
    for x in s31:
        print(x)
    print('')

    s31.optimize()

    function_t_to_s(s31, 3, 's31', f)

    f.close()


def int_and_simp(ars0, rank, complete=False):
    """ Integrates and simplifies and optimize
    arithmetic string 
    """
    # ars = simplify(ars0)
    # print('po simplify w int and simp')

    print('przed calkowaniem')
    
    ars = ars0    
    k = 1
    for x in ars:
        print(k, x)
        k += 1
    print('')


    if rank == 1:
        if complete:
            rint = ars.integrate(bra = ['α','i'], braspin = ['s']).scale(0.5)
        else:
            rint = ars.integrate(bra = ['a','i'], braspin = ['s']).scale(0.5)
    elif rank ==2:
        if complete:
            rint1 = ars.integrate(bra = ['α', 'i', 'β', 'j'], braspin =['s', 's'])
            rint2 = ars.integrate(bra = ['α', 'j', 'β', 'i'], braspin =['s', 's'])

            # rint1 = ars.integrate(bra = ['a', 'i', 'b', 'j'], braspin =['s', 's'])
            # rint2 = ars.integrate(bra = ['a', 'j', 'b', 'i'], braspin =['s', 's'])

            # rint3 = ars.integrate(bra = ['A', 'i', 'a', 'j'], braspin =['s', 's'])
            # rint4 = ars.integrate(bra = ['A', 'j', 'a', 'i'], braspin =['s', 's'])

            # rint5 = ars.integrate(bra = ['a', 'i', 'B', 'j'], braspin =['s', 's'])
            # rint6 = ars.integrate(bra = ['a', 'j', 'B', 'i'], braspin =['s', 's'])
            
            # rint7 = ars.integrate(bra = ['A', 'i', 'B', 'j'], braspin =['s', 's'])
            # rint8 = ars.integrate(bra = ['A', 'j', 'B', 'i'], braspin =['s', 's'])

            rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
            # rint  = rint1.scale(1./3.) + rint3.scale(1./3.) + rint5.scale(1./3.)+rint7.scale(1./3.) \
            #     + rint2.scale(1./6.) + rint4.scale(1./6.) + rint6.scale(1./6.) + rint8.scale(1./6.)
        else:
            rint1 = ars.integrate(bra = ['a', 'i', 'b', 'j'], braspin =['s', 's'])
            rint2 = ars.integrate(bra = ['a', 'j', 'b', 'i'], braspin =['s', 's'])
            rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
    elif rank ==3:
        rint = ars.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's'])

        
    print('to jest wynik tuz po calkowaniu')
    k =1
    for x in rint:
        print(k, x)
        k += 1
    print('')
    rint.exec_delta()
    print('to jest wynik tuz exec delta')
    k =1
    for x in rint:
        print(k, x)
        k += 1
    print('')

    if (len(rint) == 0):
        print('WYSZLO ZERO')
        return rint
    
    rsimp = simplify(rint)
    print('to jest wynik tuz po symplifikacji')
    k =1

    for x in rsimp:
        print(k, x)
        k += 1
    print('')

    rsimp.cleanup()
    print('to jest wynik tuz po cleanup')
    k =1
    for x in rsimp:
        print(k, x)
        k += 1
    print('')

#    rsimp.optimize()

    return rsimp


def int_and_simp_primitive(ars0, rank):
    """ Integrates and simplifies and optimize
    arithmetic string in primitive basis
    """
    # ars = simplify(ars0)
    # print('po simplify w int and simp')

    ars = ars0
    k = 1
    for x in ars:
        print(k, x)
        k += 1
    print('')


    if rank == 1:
        rint = ars.integrate(bra = ['a','i'], braspin = ['s'])
    elif rank ==2:
        rint = ars.integrate(bra = ['a', 'i', 'b', 'j'], braspin =['s', 's']).scale(0.5) #<----- to jest 1/2 z oepratora P
                    
    print('to jest wynik tuz po calkowaniu')
    k =1
    for x in rint:
        print(k, x)
        k += 1
    print('')
    rint.exec_delta()
    print('to jest wynik tuz exec delta')
    k =1
    for x in rint:
        print(k, x)
        k += 1
    print('')

    if (len(rint) == 0):
        print('WYSZLO ZERO')
        return rint
    
    rsimp = simplify(rint)
    print('to jest wynik tuz po symplifikacji')
    k =1

    for x in rsimp:
        print(k, x)
        k += 1
    print('')

    rsimp.cleanup()
    print('to jest wynik tuz po cleanup')
    k =1
    for x in rsimp:
        print(k, x)
        k += 1
    print('')

#    rsimp.optimize()

    return rsimp

def nazwa_funkcji_ktory_wybiera_kawalek_wzoru(ss3):

    listv = ['x','y','z']
    listo = ['u','v','w']
    ss3.multisubst(['a','b','c','i','j','k'],['x','y','z','u', 'v', 'w'])

    for x in ss3:
        idx = x.operator_idx
        x.multisubst([idx[0][0]],['a'])
        x.multisubst([idx[0][1]],['i'])
        x.multisubst([idx[1][0]],['b'])
        x.multisubst([idx[1][1]],['j'])
        x.multisubst([idx[2][0]],['c'])
        x.multisubst([idx[2][1]],['k'])
        for i in x.summation:
            if i in listv:
                x.multisubst([i],['d'])
            if i in listo:
                x.multisubst([i],['l'])
        
    listfx = ['a','b','c','i','j','k']

    for i in range(0, len(ss3)):
        temp = ugg()
        temp.summation = []
        for j in ss3[i].summation:
            if j not in listfx:
                temp.summation.append(j)

        temp.coefficient = ss3[i].coefficient
        temp.coefficient_idx = ss3[i].coefficient_idx
        temp.num_factor = ss3[i].num_factor
        ss3[i] = temp

    return ss3


def compare_num(num_a, num_b, sign):

    if sign == 'plus': 
        diff = abs(num_a-num_b)
        if diff < 1.e-05:
            return True
        else:
            return False

    if sign == 'minus':
        diff = abs(num_a+num_b)
        if diff < 1.e-05:
            return True
        else:
            return False


    
#--------------------------------------------------------FOR CAS CLASS-------------------------

def commute_cas(e1, e2):

    print('')
    print('commute_cas', e1, e2)
    print('')
    print(e1.operator_idx, e2.operator_idx)
    """
    Compute [A1, A2] commutator, where A1, A2 are
    cas instances.
    """
    if len(e2.operator_idx) >= 4:
        print('')
        print('split1')
        temp1, e2_r = e2.left_split()
        print('bede robic dwa commute', e1, 'oraz', temp1, 'i z prawej dodam', e2_r)
        print('bede robic dwa commute', e1, 'oraz', e2_r,' i z lewej dodam', temp1)        
        return commute_cas(e1, temp1).fromright(e2_r) + commute_cas(e1, e2_r).fromleft(temp1)

    if len(e1.operator_idx) >= 4:
        print('')
        print('split2')
        temp1, e1_r = e1.left_split()

        print('bede robic dwa commute', e1_r, 'oraz', e2, 'i z lewej dodam', temp1)
        print('bede robic dwa commute', temp1, 'oraz', e2,' i z prawej dodam', e1_r)        

        return commute_cas(e1_r, e2).fromleft(temp1) + commute_cas(temp1, e2).fromright(e1_r)


    return basic_commute_cas2(e1, e2)


def basic_commute_cas2(e1, e2):



    print('e1', e1)
    print('e2', e2)

        
    l1 = cas()
    l2 = cas()
    l1.operator_idx = deepcopy(e1.operator_idx)
    for x in e2.operator_idx:
        l1.operator_idx.append(x)
    l2.operator_idx = deepcopy(e2.operator_idx)
    for x in e1.operator_idx:
        l2.operator_idx.append(x)

    l1.operator_type = deepcopy(e1.operator_type)
    for x in e2.operator_type:
        l1.operator_type.append(x)
    l2.operator_type = deepcopy(e2.operator_type)
    for x in e1.operator_type:
        l2.operator_type.append(x)

    print('la', l1.operator_idx)
    print('la',l1.operator_type)
    print('la',l2.operator_idx)
    print('la',l2.operator_type)

    l1 = [l1]
    l2 = [l2]
    # l1.operator_type = deepcopy(e1.operator_type) +deepcopy(e2.operator_type)
    # l1.operator_idx = deepcopy(e2.operator_idx) + deepcopy(e1.operator_idx)
    # l1.operator_type = deepcopy(e2.operator_type) +deepcopy(e1.operator_type)
    print('l1 i l2', l1)
    print('l2', l2)
    print('')

    resmin1 = []
    for x in l1:
        result_list = x.wick_ca()
        resmin1 = resmin1 + result_list
    print('')
    resmin2 = []
    for x in l2:
        result_list = x.wick_ca()
        resmin2 = resmin2 + result_list
    print('')

    print('Result l1:')    
    for x in resmin1:
        print(x)
    print('')
    print('Result l2:')    
    for x in resmin2:
        x.num_factor *= -1.0
        print(x)
    print('')

    print('---------------------')


    for x in resmin1:
        for y in resmin2:
            res_temp = cas()
            if x.operator_idx == y.operator_idx and x.operator_type == y.operator_type :
                d1 = deepcopy(x.delta)
                d2 = deepcopy(y.delta)
                d1.sort()
                d2.sort()
                if d1== d2:
                    print('lllla', x, y)
                    x.num_factor = 0
                    y.num_factor = 0

    print('i finalnie mamy')
    for x in resmin1:
        if x.num_factor != 0:
            print(x)
    print('----------------')
    for x in resmin2:
        if x.num_factor != 0:
            print(x)

    print('')
    res_all = arithmetic_string()

    for x in resmin1:
        if x.num_factor != 0:

            e11 = deepcopy(e1)
            e22 = deepcopy(e2)
            res = cas()
            res.coefficient = e11.coefficient + e22.coefficient
            res.coefficient_idx = e11.coefficient_idx + e22.coefficient_idx
            print('lacze te     ', e11, '|', e22, '| w-->', x)
            res.num_factor = e11.num_factor * e22.num_factor * x.num_factor                
            res.summation = e11.summation + e22.summation
            res.operator_idx = x.operator_idx
            res.operator_type = x.operator_type
            res.delta = e11.delta + e22.delta
            if len(x.delta) >0:
                for i in range(0, len(x.delta)):
                    res.new_delta(x.delta[i][0], x.delta[i][1])
#            res.delta = x.delta
            print('res1', res)
            res.exec_delta()
            print('res1', res)
            res_all.append(deepcopy(res))
    print('')

    print('---------------------------------------------------------------------------')
    
    for x in resmin2:
        if x.num_factor != 0:
            e11 = deepcopy(e1)
            e22 = deepcopy(e2)            
            res = cas()
            res.coefficient = e11.coefficient + e22.coefficient
            res.coefficient_idx = e11.coefficient_idx + e22.coefficient_idx
            print('lacze te     ',e11, '|', e22, '|', x)
            # res.coefficient_idx = []
            # res.coefficient_idx.append(e11.coefficient_idx)
            # res.coefficient_idx.append(e22.coefficient_idx)
            res.num_factor = e11.num_factor * e22.num_factor * x.num_factor                
            res.summation = e11.summation + e22.summation
            res.operator_idx = x.operator_idx
            res.operator_type = x.operator_type
            res.delta = e11.delta + e22.delta
            if len(x.delta) >0:
                for i in range(0, len(x.delta)):
                    res.new_delta(x.delta[i][0], x.delta[i][1])
            # res.delta = x.delta
            res.exec_delta()
            print()
            print('res2', res)
            res_all.append(deepcopy(res))


    # for y in res_all:
            
    #     print('yyyyyyyy',y)


    return res_all
            
def basic_commute_cas(e1, e2):
    
    
    res1 = cas()
    res2 = cas()

    res1.coefficient = e1.coefficient + e2.coefficient
    res2.coefficient = e1.coefficient + e2.coefficient

    res1.num_factor = e1.num_factor * e2.num_factor
    res2.num_factor = e1.num_factor * e2.num_factor

    res1.summation = e1.summation + e2.summation
    res2.summation = e1.summation + e2.summation
    print('eeee', e1,'      ', e2)
    if e1.operator_type[0] == CRE and e1.operator_type[1] == ANI and \
       e2.operator_type[0] == CRE and e2.operator_type[1] == ANI:
        type = 1 
    elif e1.operator_type[0] == CRE and e1.operator_type[1] == ANI and \
       e2.operator_type[0] == ANI and e2.operator_type[1] == CRE:
        type = 2
    elif e1.operator_type[0] == ANI and e1.operator_type[1] == CRE and \
       e2.operator_type[0] == CRE and e2.operator_type[1] == ANI:
        type = 3
    elif e1.operator_type[0] == ANI and e1.operator_type[1] == CRE and \
       e2.operator_type[0] == ANI and e2.operator_type[1] == CRE:
        type = 4

    if type == 1:
        res1.operator_idx = [e1.operator_idx[0], e2.operator_idx[1]]
        res1.new_delta(e2.operator_idx[0], e1.operator_idx[1])
        res1.operator_type = ['+', '0']

        res2.operator_idx = [e2.operator_idx[0], e1.operator_idx[1]]
        res2.new_delta(e1.operator_idx[0], e2.operator_idx[1])
        res2.operator_type = ['+', '0']
        res2.num_factor *= -1.0
    if type == 2:
        res1.operator_idx = [e2.operator_idx[1], e1.operator_idx[1]]
        res1.new_delta(e2.operator_idx[0], e1.operator_idx[0])
        res1.operator_type = ['+', '0']
        res1.num_factor *= -1.0

        res2.operator_idx = [e1.operator_idx[0], e2.operator_idx[0]]
        res2.new_delta(e1.operator_idx[1], e2.operator_idx[1])
        res2.operator_type = ['+', '0']

    if type == 3:
        res1.operator_idx = [e1.operator_idx[1], e2.operator_idx[1]]
        res1.new_delta(e2.operator_idx[0], e1.operator_idx[0])
        res1.operator_type = ['+', '0']

        res2.operator_idx = [e2.operator_idx[0], e1.operator_idx[0]]
        res2.new_delta(e1.operator_idx[1], e2.operator_idx[1])
        res2.operator_type = ['+', '0']
        res2.num_factor *= -1.0
    if type == 4:
        res1.operator_idx = [e2.operator_idx[1], e1.operator_idx[0]]
        res1.new_delta(e2.operator_idx[0], e1.operator_idx[1])
        res1.operator_type = ['+', '0']
        res1.num_factor *= -1.0

        res2.operator_idx = [e1.operator_idx[1], e2.operator_idx[0]]
        res2.new_delta(e2.operator_idx[1], e1.operator_idx[0])
        res2.operator_type = ['+', '0']


    res1.coefficient_idx = e1.coefficient_idx + e2.coefficient_idx
    res2.coefficient_idx = e1.coefficient_idx + e2.coefficient_idx

    if e1.delta != []:
        for delta in e1.delta:
            res1.new_delta(delta[0], delta[1])
            res2.new_delta(delta[0], delta[1])

    if e2.delta != []:
        for delta in e2.delta:
            res1.new_delta(delta[0], delta[1])
            res2.new_delta(delta[0], delta[1])

    result = arithmetic_string(res1, res2)

    return result

def append_E_test_p(rint, d, m, e, n, f, l):
    
    for i in range(0, len(rint)):
        rint[i].operator_idx.append([d, m])
        rint[i].operator_idx.append([e, n])               
        rint[i].operator_idx.append([f, l])               
        rint[i].operator_type.append("s")
        rint[i].operator_type.append("s")
        rint[i].operator_type.append("s")
        
    return rint

def test_p_operator(u):

    if u==2:

        # Sprawdzenie bazy biortonormalnej dla S2
        # Tutaj sprawdzamy rownanie 14
        # <EiaEjb P_2(E_{em}E_{fn})>
        # <EiaEjb EckEdl>  (1/3<E_{ld}E_{kc} E_{em}E_{fn})> + 1/6 <E_{kd}E_{lc} E_{em}E_{fn})>)
        #       lewa                  prawa1                               prawa2
        
        # na poczatku jeszcze sprawdzamy rownanie 13

        qq = ugg()
        qq.summation = ['a', 'i', 'b', 'j']
        qq.coefficient.append('Y')
        qq.coefficient_idx.append(['a', 'b', 'i', 'j'])
        qq.operator_idx.append(["i","a"])
        qq.operator_idx.append(["j","b"])
        qq.operator_type.append("s")
        qq.operator_type.append("s")

        disambiguate(qq,t2)

        cala = qq.fromright(t2)
        print('sprawdzam rownanie 13', cala)
        rcala = cala.integrate()
        rscala = simplify(rcala)
        print('oto wynik')
        for x in rscala:
            print(x)
        print('a teraz sprawdzam rownanie 14')
            
        r0 = arithmetic_string(t2ef)
        print('r0', r0)
        prawa1 = r0.fromleft(sw1)
        prawa2 = r0.fromleft(sw2)
        
        print('')
        print('prawa1')
        print(prawa1)
        print('prawa2')
        print(prawa2)
        print('')
        rintpr1 = prawa1.integrate(nodelta=True)
        rintpr2 = prawa2.integrate(nodelta=True)
        
        print('')
        print('prawa1 strona po int')
        for x in rintpr1:
            print(x)
        print('')
        print('prawa2 strona po int')
        for x in rintpr2:
            print(x)
        print('')
        rintpr = arithmetic_string()
        
        for x in rintpr1:
            rintpr.append(x)
        for x in rintpr2:
            rintpr.append(x)

    
        r2 = arithmetic_string(ekclds)
        r =  r2.fromleft(eiajbs)
        print('')
        print('lewa strona')
        print(r)
        print('')
        rintlw = r.integrate(nodelta=True)
        print('lewa po int')
        
        print('')
        for x in rintlw:
            print(x)

        print('')
        print('teraz dodam sumowanie i wykonam delty')
        print('')
        r_all = arithmetic_string()
        k = 1
        for x in rintlw:
            for y in rintpr:
                temp = y.fromleft(x)
                temp.summation = temp.summation + ['a', 'b', 'c', 'd', 'i', 'j', 'k', 'l', 'e', 'f', 'm', 'n']
                print(k, temp)
                temp.exec_delta()
                r_all.append(temp)
                
                k += 1
        print('')
        rsimp = simplify(r_all)
        print('la')
        for  x in rsimp:
            print(x)
        sys.exit(0)
    elif u ==3:
        # Sprawdzenie P2(T2)
        #--------------------------------------------------------------------------------------------
        # P_2 (T2)
        t2ns = ugg()
        t2ns.coefficient = ["t"]
        t2ns.coefficient_idx.append(["c","k", "d", "l"])
        t2ns.operator_idx.append(["c", "k"])
        t2ns.operator_idx.append(["d", "l"])
        t2ns.operator_type.append("s")
        t2ns.operator_type.append("s")
        t2ns.num_factor = 0.5

        # r = arithmetic_string(t2ns)
        # print('r0', r)

        # rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s', 's'])
        # rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], braspin = ['s', 's'])
        # rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
        # rint = rint.scale(0.5)# to 0.5 pochodzi z rzutownika
        # # for i in range(0, len(rint)):
        # #     rint[i].operator_idx.append(["c","k"])
        # #     rint[i].operator_idx.append(["d","l"])               
        # #     rint[i].operator_type.append("s")
        # #     rint[i].operator_type.append("s")

        # for i in range(0, len(rint)):
        #     rint[i].summation = ['a', 'k', 'b', 'j', 'c', 'i','d', 'l']
            
        # rsimp  = simplify(rint)
        
        # print('rrr')
        # for x in rsimp:
        #     print(x)
        # sys.exit(0)

    elif u == 4:

        # P_3(T3)
        t3ns = ugg()
        t3ns.coefficient = ["t"]
        t3ns.coefficient_idx.append(["d","l", "e", "m", "f", "n"])
        t3ns.operator_idx.append(["d", "l"])
        t3ns.operator_idx.append(["e", "m"])
        t3ns.operator_idx.append(["f", "n"])
        t3ns.operator_type.append("s")
        t3ns.operator_type.append("s")
        t3ns.operator_type.append("s")
        t3ns.num_factor = 1.0/6.

        r = arithmetic_string(t3ns)
        print('r0', r)

        rint1 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./4.)
        rint2 = r.integrate(bra = ['a', 'k', 'b', 'i', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./12.)
        rint3 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./6.)
        rint4 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./6)
        rint5 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./12.)

        r1 = rint1 + rint2 + rint3 + rint4 + rint5
        r1 = r1.scale(1./3.)

        for i in range(0, len(r1)):
            r1[i].summation = ['a', 'k', 'b', 'j', 'c', 'i','d', 'l', 'e', 'f', 'm', 'n']
            
        rsimp  = simplify(r1)
        
        print('rrr')
        for x in rsimp:
            print(x)
        sys.exit(0)
        #--------------------------------------------------------------------------------------------

    elif u == 7:
        # sprawdzenie wektora wlasnego

        # P_3(T3)
        t3ns = ugg()
        # t3ns.coefficient = ["t"]
        # t3ns.coefficient_idx.append(["d","l", "e", "m", "f", "n"])
        t3ns.operator_idx.append(["d", "l"])
        t3ns.operator_idx.append(["e", "m"])
        t3ns.operator_idx.append(["f", "n"])

        # t3ns.operator_idx.append(["a", "i"])
        # t3ns.operator_idx.append(["b", "j"])
        # t3ns.operator_idx.append(["c", "k"])

        t3ns.operator_type.append("s")
        t3ns.operator_type.append("s")
        t3ns.operator_type.append("s")
        # t3ns.num_factor = 1.0/6.

        r = arithmetic_string(t3ns)
        print('r0', r)

        rint1 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's'])
        for x in rint1:
            print(x)
        sys.exit(0)

        rint1 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's']).scale(17.0)#.scale(5.0)
        rint2 = r.integrate(bra = ['a', 'i', 'b', 'k', 'c', 'j'], braspin = ['s', 's', 's']).scale(-1.0)#.scale(-1.0)
        rint3 = r.integrate(bra = ['a', 'k', 'b', 'i', 'c', 'j'], braspin = ['s', 's', 's']).scale(-1.0)#.scale(-1.0)
        rint4 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's']).scale(-7.0)#.scale(-1.0)
        rint5 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k'], braspin = ['s', 's', 's']).scale(-7.0)#.scale(-1.0)
        rint6 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i'], braspin = ['s', 's', 's']).scale(-1.0)#.scale(-1.0)


        r1 = rint1 + rint2 + rint3 + rint4 + rint5 + rint6



        # r1 = r1.scale(1./12.)
        # rsimp = simplify(r1)
        # print('wynik calki nakrywania')
        # for x in rsimp:
        #     print(x)
        # sys.exit(0)
        
        
        r1 = append_E_test_p(r1, 'a', 'i', 'b', 'j', 'c', 'k')
        print('fffff')
        for x in r1:
            print(x)

        r1 = r1.scale(1./12.)

        for i in range(0, len(r1)):
            r1[i].summation = ['a', 'k', 'b', 'j', 'c', 'i','d', 'l', 'e', 'f', 'm', 'n']
#            r1[i].summation = ['d', 'l', 'e', 'f', 'm', 'n']
            
        rsimp  = simplify(r1)
        
        print('P_3(T_3)=')
        for x in rsimp:
            print(x)
        sys.exit(0)
        r = rsimp
        for i in range(0, len(r)):
            r[i].summation = []

        rint1 = r.integrate(bra = ['d', 'l', 'e', 'm', 'f', 'n'], braspin = ['s', 's', 's']).scale(5.0)
        rint2 = r.integrate(bra = ['d', 'l', 'e', 'n', 'f', 'm'], braspin = ['s', 's', 's']).scale(-1.0)
        rint3 = r.integrate(bra = ['d', 'n', 'e', 'l', 'f', 'm'], braspin = ['s', 's', 's']).scale(-1.0)
        rint4 = r.integrate(bra = ['d', 'n', 'e', 'm', 'f', 'l'], braspin = ['s', 's', 's']).scale(-1.0)
        rint5 = r.integrate(bra = ['d', 'm', 'e', 'l', 'f', 'n'], braspin = ['s', 's', 's']).scale(-1.0)
        rint6 = r.integrate(bra = ['d', 'm', 'e', 'n', 'f', 'l'], braspin = ['s', 's', 's']).scale(-1.0)


        r1 = rint1 + rint2 + rint3 + rint4 + rint5 + rint6
        r1 = append_E_test_p(r1, 'd', 'l', 'e', 'm', 'f', 'n')
        print('fffff')
        for x in r1:
            print(x)

        r1 = r1.scale(1./12.)

        for i in range(0, len(r1)):
            r1[i].summation = ['a', 'k', 'b', 'j', 'c', 'i','d', 'l', 'e', 'f', 'm', 'n']
            
        rsimp  = simplify(r1)
        
        print('P_3(T_3)=')
        for x in rsimp:
            print(x)



            
        sys.exit(0)
        #--------------------------------------------------------------------------------------------

        
    elif u == 5:
        
        #--------------------------------------------------------------------------------------------

        # Sprawdzenie bazy biortonormalnej dla S3
        # Tutaj sprawdzamy rownanie 14
        # <EiaEjbEkc P_3(E_{dl}E_{em}E_{fn})>
        # <EiaEjbEkc EdlEemEfn>  ---???(1/3<E_{ld}E_{kc} E_{em}E_{fn})> + 1/6 <E_{kd}E_{lc} E_{em}E_{fn})>)
        #       lewa                  prawa1                               prawa2
        
        


        proj = ugg()
        proj.operator_idx.append(["l","d"])
        proj.operator_idx.append(["m","e"])
        proj.operator_idx.append(["n","f"])               
        proj.operator_type.append("s")
        proj.operator_type.append("s")
        proj.operator_type.append("s")
        #        proj.num_factor = 1./3.

        # sw2 = ugg()
        # sw2.operator_idx.append(["k","d"])
        # sw2.operator_idx.append(["l","c"])
        # sw2.operator_type.append("s")
        # sw2.operator_type.append("s")
        # sw2.num_factor = 1./6.

        # t3ns = ugg()
        # t3ns.coefficient = ["t"]
        # t3ns.coefficient_idx.append(["a","i", "b", "j", "c", "k"])
        # t3ns.operator_idx.append(["a", "i"])
        # t3ns.operator_idx.append(["b", "j"])
        # t3ns.operator_idx.append(["c", "k"])
        # t3ns.operator_type.append("s")
        # t3ns.operator_type.append("s")
        # t3ns.operator_type.append("s")
        # t3ns.num_factor = 1.0/6.

        t3ns = ugg()
        t3ns.coefficient = ["t"]
        t3ns.coefficient_idx.append(["d","l", "e", "m", "f", "n"])
        t3ns.operator_idx.append(["d", "l"])
        t3ns.operator_idx.append(["e", "m"])
        t3ns.operator_idx.append(["f", "n"])
        t3ns.operator_type.append("s")
        t3ns.operator_type.append("s")
        t3ns.operator_type.append("s")
        # t3ns.num_factor = 1.0/6.


        r = arithmetic_string(t3ns)
        print('r0', r)
        # rint1 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's'])
        # rsimp = simplify(rint1)
        # print('wynik')
        # for x in rsimp:
        #     print(x)
        # sys.exit(0)

        # rint0 = r.integrate(bra = ['d', 'l', 'e', 'm', 'f', 'n'], braspin = ['s', 's', 's'])
        # for i in range(0, len(rint0)):
        #     rint0[i].operator_idx.append(["d","l"])
        #     rint0[i].operator_idx.append(["e","m"])               
        #     rint0[i].operator_idx.append(["f","n"])               
        #     rint0[i].operator_type.append("s")
        #     rint0[i].operator_type.append("s")
        #     rint0[i].operator_type.append("s")
        
        # rint1 = r.integrate(bra = ['d', 'l', 'e', 'n', 'f', 'm'], braspin = ['s', 's', 's'])
        # for i in range(0, len(rint1)):
        #     rint1[i].operator_idx.append(["d","l"])
        #     rint1[i].operator_idx.append(["e","n"])               
        #     rint1[i].operator_idx.append(["f","m"])               
        #     rint1[i].operator_type.append("s")
        #     rint1[i].operator_type.append("s")
        #     rint1[i].operator_type.append("s")

        # rint2 = r.integrate(bra = ['d', 'n', 'e', 'l', 'f', 'm'], braspin = ['s', 's', 's'])
        # for i in range(0, len(rint2)):
        #     rint2[i].operator_idx.append(["d","n"])
        #     rint2[i].operator_idx.append(["e","l"])               
        #     rint2[i].operator_idx.append(["f","m"])               
        #     rint2[i].operator_type.append("s")
        #     rint2[i].operator_type.append("s")
        #     rint2[i].operator_type.append("s")

        # rint3 = r.integrate(bra = ['d', 'n', 'e', 'm', 'f', 'l'], braspin = ['s', 's', 's'])
        # for i in range(0, len(rint3)):
        #     rint3[i].operator_idx.append(["d","n"])
        #     rint3[i].operator_idx.append(["e","m"])               
        #     rint3[i].operator_idx.append(["f","l"])               
        #     rint3[i].operator_type.append("s")
        #     rint3[i].operator_type.append("s")
        #     rint3[i].operator_type.append("s")

        # rint4 = r.integrate(bra = ['d', 'm', 'e', 'l', 'f', 'n'], braspin = ['s', 's', 's'])
        # for i in range(0, len(rint4)):
        #     rint4[i].operator_idx.append(["d","m"])
        #     rint4[i].operator_idx.append(["e","l"])               
        #     rint4[i].operator_idx.append(["f","n"])               
        #     rint4[i].operator_type.append("s")
        #     rint4[i].operator_type.append("s")
        #     rint4[i].operator_type.append("s")

        # rint5 = r.integrate(bra = ['d', 'm', 'e', 'n', 'f', 'l'], braspin = ['s', 's', 's'])
        # for i in range(0, len(rint5)):
        #     rint5[i].operator_idx.append(["d","m"])
        #     rint5[i].operator_idx.append(["e","n"])               
        #     rint5[i].operator_idx.append(["f","l"])               
        #     rint5[i].operator_type.append("s")
        #     rint5[i].operator_type.append("s")
        #     rint5[i].operator_type.append("s")

        # r1 = rint0 + rint1 + rint2 + rint3 + rint4 + rint5
        # for i in range(0, len(r1)):
        #     r1[i].summation = ['a', 'k', 'b', 'j', 'c', 'i','d', 'l', 'e', 'm', 'f', 'n']
        # rsimp1 = simplify(r1)


        rint1 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./4.)
        rint2 = r.integrate(bra = ['a', 'k', 'b', 'i', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./12.)
        rint3 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./6.)
        rint4 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./6)
        rint5 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./12.)

        rint1 = append_E_test_p(rint1, 'a', 'i', 'b', 'j', 'c', 'k')
        rint2 = append_E_test_p(rint2,'a', 'k', 'b', 'i', 'c', 'j')
        rint3 = append_E_test_p(rint3,'a', 'k', 'b', 'j', 'c', 'i')
        rint4 = append_E_test_p(rint4,'a', 'j', 'b', 'i', 'c', 'k')
        rint5 = append_E_test_p(rint5,'a', 'j', 'b', 'k', 'c', 'i')

        r1 = rint1 + rint2 + rint3 + rint4 + rint5
        rsimp1 = simplify(r1)
        
        rint1 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./12.)
        rint2 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./4.)
        rint3 = r.integrate(bra = ['a', 'i', 'b', 'k', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./6.)
        rint4 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./6.)
        rint5 = r.integrate(bra = ['a', 'k', 'b', 'i', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./12.)

        rint1 = append_E_test_p(rint1, 'a', 'j', 'b', 'k', 'c', 'i')
        rint2 = append_E_test_p(rint2, 'a', 'i', 'b', 'j', 'c', 'k')
        rint3 = append_E_test_p(rint3, 'a', 'i', 'b', 'k', 'c', 'j')
        rint4 = append_E_test_p(rint4, 'a', 'k', 'b', 'j', 'c', 'i')
        rint5 = append_E_test_p(rint5, 'a', 'k', 'b', 'i', 'c', 'j')

        
        r2 = rint1 + rint2 + rint3 + rint4 + rint5
   #     rsimp2 = simplify(r2)

        rint1 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./6.)
        rint2 = r.integrate(bra = ['a', 'i', 'b', 'k', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./6.)
        rint3 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./3.)
        rint4 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./6.)
        rint5 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./6.)

        rint1 = append_E_test_p(rint1, 'a', 'k', 'b', 'j', 'c', 'i')
        rint2 = append_E_test_p(rint2, 'a', 'i', 'b', 'k', 'c', 'j')
        rint3 = append_E_test_p(rint3, 'a', 'i', 'b', 'j', 'c', 'k')
        rint4 = append_E_test_p(rint4, 'a', 'j', 'b', 'k', 'c', 'i')
        rint5 = append_E_test_p(rint5, 'a', 'j', 'b', 'i', 'c', 'k')
        
        r3 = rint1 + rint2 + rint3 + rint4 + rint5
   #     rsimp3 = simplify(r3)
        rint1 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./6.)
        rint2 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./6.)
        rint3 = r.integrate(bra = ['a', 'k', 'b', 'i', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./6.)
        rint4 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./3.)
        rint5 = r.integrate(bra = ['a', 'i', 'b', 'k', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./6.)

        rint1 = append_E_test_p(rint1, 'a', 'j', 'b', 'i', 'c', 'k')
        rint2 = append_E_test_p(rint2, 'a', 'k', 'b', 'j', 'c', 'i')
        rint3 = append_E_test_p(rint3, 'a', 'k', 'b', 'i', 'c', 'j')
        rint4 = append_E_test_p(rint4, 'a', 'i', 'b', 'j', 'c', 'k')
        rint5 = append_E_test_p(rint5, 'a', 'i', 'b', 'k', 'c', 'j')

        r4 = rint1 + rint2 + rint3 + rint4 + rint5
 #       rsimp4 = simplify(r4)

        rint1 = r.integrate(bra = ['a', 'k', 'b', 'i', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./12.)
        rint2 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i'], braspin = ['s', 's', 's']).scale(1./12.)
        rint3 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./6.)
        rint4 = r.integrate(bra = ['a', 'i', 'b', 'k', 'c', 'j'], braspin = ['s', 's', 's']).scale(1./6.)
        rint5 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's']).scale(1./4.)

        rint1 = append_E_test_p(rint1, 'a', 'k', 'b', 'i', 'c', 'j')
        rint2 = append_E_test_p(rint2, 'a', 'j', 'b', 'k', 'c', 'i')
        rint3 = append_E_test_p(rint3, 'a', 'j', 'b', 'i', 'c', 'k')
        rint4 = append_E_test_p(rint4, 'a', 'i', 'b', 'k', 'c', 'j')
        rint5 = append_E_test_p(rint5,'a', 'i', 'b', 'j', 'c', 'k')

        r5 = rint1 + rint2 + rint3 + rint4 + rint5
        #        rsimp5 = simplify(r5)
        
        rsimp = r1 + r2+r3+r4+r5
        #       rsimp = rsimp1+rsimp2+rsimp3+rsimp4+rsimp5
        for i in range(0, len(rsimp1)):
            rsimp1[i].summation = ['a', 'k', 'b', 'j', 'c', 'i','d', 'l', 'e', 'm', 'f', 'n']
            #            rsimp1[i].exec_delta()
            #           print(rsimp1[i])
            
            #        rintpr1  = simplify(rsimp)
        rsimp1  = simplify(rsimp1)
        
        print('rrr')
        for x in rsimp1:
            print(x)
            
    elif u == 6:

        # sprawdzenie rownania 13 i 14 dla P3
        
        qq = ugg()
        qq.summation = ['a', 'i', 'b', 'j', 'c', 'k']
        qq.coefficient.append('Y')
        qq.coefficient_idx.append(['a', 'b', 'c', 'i', 'j', 'k'])
        qq.operator_idx.append(["i","a"])
        qq.operator_idx.append(["j","b"])
        qq.operator_idx.append(["k","c"])
        qq.operator_type.append("s")
        qq.operator_type.append("s")
        qq.operator_type.append("s")
        
        disambiguate(qq,t3)
        
        cala = qq.fromright(t3)
        print('sprawdzam rownanie 13', cala)
        rcala = cala.integrate()
        rscala = simplify(rcala)
        print('oto wynik')
        for x in rscala:
            print(x)
        print('a teraz sprawdzam rownanie 14')

        
        t3ns = ugg()
        # t3ns.coefficient = ["t"]
        # t3ns.coefficient_idx.append(["d","l", "e", "m", "f", "n"])
        t3ns.operator_idx.append(["d", "l"])
        t3ns.operator_idx.append(["e", "m"])
        t3ns.operator_idx.append(["f", "n"])
        t3ns.operator_type.append("s")
        t3ns.operator_type.append("s")
        t3ns.operator_type.append("s")
        # t3ns.num_factor = 1.0/6.

        r = arithmetic_string(t3ns)
        print('r0', r)


        rint1 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's'], nodelta=True).scale(5.0)
        rint2 = r.integrate(bra = ['a', 'i', 'b', 'k', 'c', 'j'], braspin = ['s', 's', 's'], nodelta=True).scale(-1.0)
        rint3 = r.integrate(bra = ['a', 'k', 'b', 'i', 'c', 'j'], braspin = ['s', 's', 's'], nodelta=True).scale(-1.0)
        rint4 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's'], nodelta=True).scale(-1.0)
        rint5 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k'], braspin = ['s', 's', 's'], nodelta=True).scale(-1.0)
        rint6 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i'], braspin = ['s', 's', 's'], nodelta=True).scale(-1.0)

        # rint1 = r.integrate(bra = ['a', 'i', 'b', 'j', 'c', 'k'], braspin = ['s', 's', 's'], nodelta=True).scale(1./4.)
        # rint2 = r.integrate(bra = ['a', 'k', 'b', 'i', 'c', 'j'], braspin = ['s', 's', 's'], nodelta=True).scale(1./12.)
        # rint3 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's'], nodelta=True).scale(1./6.)
        # rint4 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k'], braspin = ['s', 's', 's'], nodelta=True).scale(1./6)
        # rint5 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i'], braspin = ['s', 's', 's'], nodelta=True).scale(1./12.)
#        rint6 = r.integrate(bra = ['a', 'i', 'b', 'k', 'c', 'j'], braspin = ['s', 's', 's'], nodelta=True)#.scale(1./4.)
        rintpr = rint1 + rint2 + rint3 + rint4 + rint5 + rint6
        # rsimp = simplify(rintpr)
        # print('wynik')
        # for x in rsimp:
        #     print(x)
        # sys.exit(0)
        
        r2 = arithmetic_string(e3s)

        r =  r2.fromleft(eiajbkcs)
        print('')
        print('lewa strona')
        print(r)
        print('')
        rintlw = r.integrate(nodelta=True)
        print('lewa po int')

        print('')
        for x in rintlw:
            print(x)
            
        print('')
        print('teraz dodam sumowanie i wykonam delty')
        print('')
        r_all = arithmetic_string()
        k = 1
        for x in rintlw:
            for y in rintpr:
                temp = y.fromleft(x)
                temp.summation = temp.summation + ['a', 'b', 'c', 'd', 'i', 'j', 'k', 'l', 'e', 'f', 'm', 'n', 'a>', 'b>', 'c>', 'i>', 'j>', 'k>']
                print(k, temp)
                temp.exec_delta()
                r_all.append(temp)
                
                k += 1
        print('')
        rsimp = simplify(r_all)
        print('la')
        for  x in rsimp:
            print(x)
        sys.exit(0)
    #--------------------------------------------------------------------------------------------

    # combo = ugg()
    # combo.summation = ['a', 'b', 'i', 'j', "e","m","f","n"]
    # combo.coefficient.append('Y')
    # combo.coefficient.append(CC_AMPLITUDE)
    # combo.coefficient_idx.append(["a",'i', 'b', 'j'])    
    # combo.coefficient_idx.append(["e","m", "f", "n"])
    # combo.operator_idx.append(["i","a"])
    # combo.operator_idx.append(["j","b"])
    # combo.operator_idx.append(["e", "m"])
    # combo.operator_idx.append(["f", "n"])
    # combo.operator_type.append("s")
    # combo.operator_type.append("s")
    # combo.operator_type.append("s")
    # combo.operator_type.append("s")
    # combo.num_factor = 1./2.

    # rcombo = arithmetic_string(combo)
    # rint = rcombo.integrate()
    
    # print(combo)
    # print('')
    # for x in rint:
    #     x.exec_delta()
    #     print(x)
                
    #     rsimp = simplify(rint)
    # print('la')
    # for  x in rsimp:
    #     print(x)

    # sys.exit(0)

    
    # r = evaluate(obs, t2)

    # print('disamb')
    # r2 = arithmetic_string()
    # for x in r:
    #     rr = deepcopy(eckdls)
    #     disambiguate(x, rr)
    #     print(x)
    #     print('disdis', x, rr)
    #     r2.append(x.fromleft(rr))
    #     print(x)
    # print('')
    # # r = r.fromleft(eckdls)
    # print('')
    # print('wynik fromleft')
    # k = 0
    # for x in r2:
    #     print(k, x)
    #     k+= 1

    
    # r = evaluate(t2c, t2, t2)
    # print('')
    # print('wynik')
    # for x in r:
    #     print(x)
    # sys.exit(0)
    
    # ekl = ugg()
    # ekl.summation = ['c', 'i', 'b', 'j']
    # ekl.operator_idx.append(["c","i"])
    # ekl.operator_idx.append(["b","j"])
    # ekl.operator_type.append("s")
    # ekl.operator_type.append("s")

    # r = arithmetic_string(ekl)
    # print(r)
    # rint = r.integrate()
    # for x in rint:
    #     print(x)
    # sys.exit(0)
    
    # r = evaluate(ekl, t2)#arithmetic_string(t2)
    # print('wynik')
    # for x in r:
    #     print(x)
    # sys.exit(0)

    # r = arithmetic_string(t2)

#     r = evaluate(obsab, t2)

#     # r = arithmetic_string(ekl)
    
#     print('wynik evaluat')
#     for x in r:
#         print(x)
#     print('')
# #    sys.exit(0)
#     # sys.exit(0)
#     # rint1 = r.integrate(bra = ['c', 'k', 'd', 'l'], braspin =['s', 's'])
#     # rint2 = r.integrate(bra = ['c', 'l', 'd', 'k'], braspin =['s', 's'])
#     # rint  = rint1.scale(1./3.) + rint2.scale(1./6.)

#     # r2 = deepcopy(r)
#     print('disamb')
#     r2 = arithmetic_string()
#     for x in r:
#         rr = deepcopy(eckdls)
#         disambiguate(x, rr)
#         print(x)
#         print('disdis', x, rr)
#         r2.append(x.fromleft(rr))
#         print(x)
#     print('')
#     # r = r.fromleft(eckdls)
#     print('')
#     print('wynik fromleft')
#     k = 0
#     for x in r2:
#         print(k, x)
#         k+= 1
# #    sys.exit(0)
#     # r2 = r.fromleft(eckdls2)
#     print(r2)
#     # rint = r.integrate(bra = ['c', 'k', 'd', 'l'], braspin =['s', 's'])
#     rint = r2.integrate()
#     print('integrate')
#     for x in rint:
#         print(x)
#     # rint2 = r2.integrate()
#     # rint  = rint1.scale(1./3.) + rint2.scale(1./6.)
#     rsimp = simplify(rint)
#     print('wynik')
#     for x in rsimp:
#         print(x)
#     sys.exit(0)


