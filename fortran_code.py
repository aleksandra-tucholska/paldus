# ------------------------------------------------------
# Author: Aleksandra Tucholska, University of Warsaw
# ------------------------------------------------------
from params import *
from copy import deepcopy
import paldus_classes
from paldus_classes import *
from paldus_classes import ugg
from paldus_classes import arithmetic_string
import math
from params import virtual, occupied, general, EPSILON
from itertools import combinations
import pickle
import datetime
from eomccjac import whichvk
from eomccjac import whichvk_triplet
from eomccjac import whichvk_triplet2
from paldus_classes import free_idx


#
# Treshold for davidson parallelization
#
MEM_THRESH_FOR_DISK_INTERM = 8
def is_redundant_triplet(BRA, KET, nm, lmfold, rmfold):

    a = nm[0]
    i = nm[1]
    b = nm[2]
    j = nm[3]
    
    if(BRA + KET) > 2:
        c = nm[4]
        k = nm[5]

    if (BRA + KET) > 3:
        d = nm[6]
        l = nm[7]
    
    if (BRA + KET) > 4:
        e = nm[8]
        m = nm[9]


    cond_virt = []
    cond_occ = []

    if BRA == 1:
        if KET == 2:
            if rmfold == 'p':
                if b == c or j == k:
                    return True
            elif rmfold == 'm':
                if b == c and j == k:
                    return True
        if KET == 3:
            if c == d or k == l:
                return True
    if BRA == 2:
        if lmfold == 'p':
            if a == b or i == j:
                return True
            cond_add([a, b], cond_virt)
            cond_add([i, j], cond_occ)
        elif lmfold == 'm':
            if a == b and i == j:
                return True
            elif a != b and i != j:
                cond_add([a, b], cond_virt)
            elif a != b and i == j:
                cond_add([a, b], cond_virt)
            elif a == b and i != j:
                cond_add([i, j], cond_occ)
        if KET == 2:
            if rmfold == 'p':
                if c == d or k == l:
                    return True
                cond_add([c, d], cond_virt)
                cond_add([k, l], cond_occ)
            elif rmfold == 'm':
                if c == d and k == l:
                    return True
                elif c != d and k != l:
                    cond_add([c, d], cond_virt)
                elif c != d and k == l:
                    cond_add([c, d], cond_virt)
                elif c == d and k != l:
                    cond_add([k, l], cond_occ)
        elif KET == 3:
            if d == e or l == m:
                return True
            cond_add([d, e], cond_virt)
            cond_add([l, m], cond_occ)

    elif BRA == 3:
        if b == c or j == k:
            print('b==c, j==k')
            return True
        cond_add([b, c], cond_virt)
        cond_add([j, k], cond_occ)
        
        if KET == 2:
            if rmfold == 'p':
                if d == e or l == m:
                    return True
                cond_add([d, e], cond_virt)
                cond_add([l, m], cond_occ)
            elif rmfold == 'm':
                if d == e and l == m:
                    return True
                elif d != e and l != m:
                    cond_add([d, e], cond_virt)
                elif d != e and l == m:
                    cond_add([d, e], cond_virt)
                elif d == e and l != m:
                    cond_add([l, m], cond_occ)

    is_contra_virt = is_contradictory(cond_virt)
    is_contra_occ = is_contradictory(cond_occ)

    if is_contra_virt or is_contra_occ:
#        print('odrzucam  nowe', nm)
        return True
    else:
 #       print('zostawiam nowe', nm)
        return False
    
                
def is_contradictory(cond_list):
    
    for i in range(0, len(cond_list)):
        p = cond_list[i][0]
        q = cond_list[i][1]
        cond_rev = [q, p]
        for j in range(i+1, len(cond_list)):
            if cond_rev == cond_list[j]:
                return True

    return False

def mem2(x):

    memc = [0, 0, 0]

    memidx = []

    for i in range(0, len(x.coefficient)):
        if 'interm' in x.coefficient[i]:
            for idx in x.coefficient_idx[i]:
                if idx in CABS:
                    memc[0] += 1
                elif idx in virtualall:
                    memc[1] += 1
                elif idx in occupied:
                    memc[2] += 1
                memidx.append(idx)

    return memc


def cond_add(cond, cond_list):
    """
    Add new cond operations.
    """
    if cond not in cond_list:
        cond_list.append(cond)


def is_redundant(BRA, KET, nm):
    """ For <BRA|X|KET> check if
    creation of function of name 'nm'
    is redundant
    """
    
    #
    # Create empty dictionaries
    #
    dict_gt = {}
    dict_lt = {}
    for i in range(0, len(nm)):
        dict_gt[nm[i]] = []
        dict_lt[nm[i]] = []
        
    #
    # Scan bra to find the condintions, for breaking loop
    # Break loop if conditions are broken
    #
    for i in range(0, 2*BRA-2, 2):
        if is_contradiction(i, nm, dict_gt, dict_lt):
            return True

    #
    # Scan ket to find the condintions, for breaking loop
    # Break loop if conditions are broken
    #
    for i in range(2*BRA, 2*BRA + 2*KET-2, 2):
        if is_contradiction(i, nm, dict_gt, dict_lt):
            return True
        
    return False
    
                
def is_contradiction(i, nm, dict_gt, dict_lt) :

    if nm[i] != nm[i+2]:
        if nm[i] not in dict_gt[nm[i+2]]:
            if nm[i+2] not in dict_lt[nm[i]]:
                dict_gt[nm[i]].append(nm[i+2])
                dict_lt[nm[i+2]].append(nm[i])
                dict_update(dict_gt)
                dict_update(dict_lt)
            else:
                return True
        else:
            return True
    else:
        if nm[i+1] != nm[i+3]:
            if nm[i+1] not in dict_gt[nm[i+3]]:
                if nm[i+3] not in dict_lt[nm[i+1]]:
                    dict_gt[nm[i+1]].append(nm[i+3])
                    dict_lt[nm[i+3]].append(nm[i+1])
                    dict_update(dict_gt)
                    dict_update(dict_lt)
                else:
                    return True
            else:
                return True
        
    return False
    
def dict_update(dct):
    
    for i in dct.keys():
        for j in dct[i]:
            if j in dct.keys():
                dct[i] += dct[j]

def delete_permutations(cluster):

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
    return cluster

def simplify_fort(ars):
    """ simplifies result
    """
    ars.as_standarize()

    ars_temp = deepcopy(ars.expand())
    ars_temp.cluster()
    result = arithmetic_string()

    for key in ars_temp.clusters:
        cluster = ars_temp.clusters[key]
        if len(cluster) > 1:
            cluster = identify_equal(cluster)
            cluster.cleanup()

        for elem in cluster:
            result.append(elem)

    return result


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

def create_omega(name):

    """ name = aibj | ck
"""
    s = ''

    if name[0] == name[2] and name[1] == name[3]: #aiaick
        if name[4] == name[0]:
            s += omega1_for_21(name[0], name[1], name[3], name[5])
        if name[5] == name[1]:
            s += omega2_for_21(name[0], name[2], name[4], name[3])
    else:
        if name[4] == name[0]:  #aibjak   c = a
            s += omega1_for_21(name[2], name[3], name[1], name[5])
        if name[4] == name[2]: #aibjbk  c = b
            s += omega1_for_21(name[0], name[1], name[3], name[5])
        if name[5] == name[1]: #aibjci  k = i
            s += omega2_for_21(name[0], name[2], name[4], name[3])
        if name[5] == name[3]: #aibjcj   k = j
            s += omega2_for_21(name[2], name[0], name[4], name[1])

    return s

def omega1_for_21(a, i, j, k):
    
    s = """
    triple_w1 = triple_w1 + 2.d+0 * w1({a}, {i}, {j}, {k})
    """.format(a = a, i = i, j = j, k = k)

    return s

def omega2_for_21(a, b, c, i):
    
    s = """
    triple_w2 = triple_w2 + 2.d+0 * w2({a}, {b}, {c}, {i})
    """.format(a = a, b = b, c = c, i = i)

    return s

def arstofort(ars, nterms=0, all_interm_dict = {}, f12=False, int_app = []):
    """ Translates arithmetic string 
    to fortran code.  Each ugg in arithmetic
    string is translated to set of loops. There is no
    declaration of variables. Returns string 's',
    containg fortran code.
    """
    print('ZACZYNAM')
    nocc = False
    nactive = False

    if f12:
        for x in ars:
            x.optimize_f12()
    print('popoopt')
    s = ''
    if f12:
        ars.cluster_for_fortran_f12(all_interm_dict)
    else:
        ars.cluster_for_fortran()

    if (nterms != 0):
        y = nterms
        y2 = y
    else:
        y = -1
        y2 = y
    print(len(ars.clusters))
    z = 0
    print(ars.clusters.keys())

    # for key in ars.clusters:
    #     cluster = deepcopy(ars.clusters[key])
    #     print(len(cluster))
    #     for y in cluster:
    #         print(y)
    #     print('')
    # sys.exit(0)
    
    if f12 == True:
        new_key = max(ars.clusters.keys())+1
        print('to jest lista kluczy', ars.clusters.keys())
        print('to jest  new key na samym poczatku', new_key)
        key_to_del = []
        new_keys_list = []
        new_clusters = {}
        qq = 0 
        for key in ars.clusters:
            cluster = deepcopy(ars.clusters[key])
            print('taki jest atkualny cluster dla ', qq, key, 'dla dlugosci', len(cluster))
            z = 1
            for i in range(0, len(cluster)):
                print(z, cluster[i])
                z += 1
                
            print('----------------------------------------------------')
            print('')

            mem_batch = 0
            key_to_del.append(key)
            new_clusters[new_key] = arithmetic_string()
            new_keys_list.append(new_key)

            for i in range(0, len(cluster)):
                x = cluster[i]
                print('rozwazam x i dodaje do new_clusters', i, len(cluster), x)
                new_clusters[new_key].append(deepcopy(x))

                for coef in x.coefficient:
                    if 'interm' in coef:
                        coef_name = extract_coef(coef)
                        if all_interm_dict[coef_name]['mem'] == 'disk':
                            memidx = mem2(x)
                            memreal = (888**(memidx[0]) * 232**(memidx[1]) * 20**memidx[2])*8/(10**9)
                            mem_batch += memreal
                            print('byl na liscie interm disk i licze pamiec', memreal, mem_batch)
                if mem_batch >= MEM_THRESH_FOR_DISK_INTERM:
                    print('mem jest wieksze i dziele po ', i, 'wyrazach',  mem_batch)
                    print('taki powstal cluster')
                    ww = 0
                    for kk in new_clusters[new_key]:
                        print(ww, kk)
                        ww += 1
                    print('')
                    print('zeruje mem, zwiekszam new_key', new_key)
                    new_keys_list.append(new_key)
                    new_key += 1
                    mem_batch = 0
                    new_clusters[new_key] = arithmetic_string()
                else:
                    print('nie dziele i dalej mam taki cluster do ktorego dalej bede dodawac')
                    xx = 0
                    for t in new_clusters[new_key]:
                        print(xx, new_key, t)
                        xx += 1
            new_key += 1

        rr = 0
        print('new_keys_list', new_keys_list)
        for ktd in key_to_del:
            ars.clusters.pop(ktd)
        for ktd in new_keys_list:
            ars.clusters[ktd] = arithmetic_string()
            for i in new_clusters[ktd]:
                ars.clusters[ktd].append(i)
                rr += 1

            
        print('to jest lista kluczy po dodaniu nowych', ars.clusters.keys(), rr)


        
        # ars.cluster_for_fortran_f12(all_interm_dict, MEM_THRESH_FOR_DISK_INTERM)            
    all_int_type = []
    for key in ars.clusters:
        print('jestem na', z, 'key', key)

        
        cluster = ars.clusters[key]

        print('len(cluster)', len(cluster))

        if f12 == True:
            gg = add_read_disk_f12(cluster, all_interm_dict, int_app)
            s += gg
            
        
        # for x in cluster:
        #     print(' i jak?')
        #     if all_interm_dict != {}:
        #         print('ano')
        #         if f12 == True:
                    
        #             gg = add_read_disk_f12(x, all_interm_dict)
        #             print('ggrr', gg, x)
        #             s += gg
        #     else:
        #         print('nie')    


        # tutaj zamieniam indeksy i> na i1 dla fortrana
        if f12 == True:
            print('jestem1')
            subst_list_old = []
            subst_list_new = []
            for i in cluster[0].summation:
                if ">" in i:
                    if i in virtual:
                        ii = free_idx(fortran_virtual, subst_list_new)
                        subst_list_old.append(deepcopy(i))
                        subst_list_new.append(deepcopy(ii))
                    if i in occupied:
                        ii = free_idx(fortran_occupied, subst_list_new)
                        subst_list_old.append(deepcopy(i))
                        subst_list_new.append(deepcopy(ii))
                if i in CABS:
                    ii = free_idx(fortran_CABS, subst_list_new)
                    subst_list_old.append(deepcopy(i))
                    subst_list_new.append(deepcopy(ii))
            for i in range(0, len(cluster)):
                print("subst1", i, cluster[i])
                print(subst_list_old)
                print(subst_list_new)
                for j in range(0, len(subst_list_old)):
                    cluster[i].substitute(subst_list_old[j], subst_list_new[j])
                print("subst2", cluster[i])
                                
                        
        for x in cluster:
            print('szukam B2_', x)

                        
        for i in cluster[0].summation:
            print('jestem2')
                        
            if i in virt_all:
                s += "do {i} = nocc + 1, nactive \n".format(i = i)
                nocc = True
                nactive = True
            if i in occ_all:
                s += "do {i} = 1, nocc \n".format(i = i)
                nocc = True

            if i in CABS_all:
                s += "do {i} = ncabs0, ncabs1 \n".format(i = i)
                ncabs= True
                
        y2 = y
        print('jestem3')
        print('key', key)
        for x in cluster:
            print('xrxrxrxrxr', x)
            for i in range(0, len(x.coefficient)):
                print('jestem420', x, x.coefficient[i])
            if all_interm_dict != {}:
                if f12 == False:
                    s += add_read_disk(x, all_interm_dict)
            print('xxxx', y+1, x)
            y += 1
            if(len(x.coefficient) == 0):
                s += "term({y}) = one".format(y = y)
            else:
                s += "term({y}) = term({y}) + ".format(y = y)
            for i in range(0, len(x.coefficient)):
                print('jestem4')

                if x.coefficient[i] == TWOEL_INT:
                    g_type = x.int_type("", i)
                    if(g_type not in all_int_type):
                        all_int_type.append(g_type)
                    print('all_int_type', all_int_type)

                    idx = x.coefficient_idx[i]
    
                    if (i + 1) < len(x.coefficient):
                        s += ("{g_type}({i0}, {i1}, {i2}, {i3}) * ".format(g_type = g_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                    else :
                        s += ("{g_type}({i0}, {i1}, {i2}, {i3})".format(g_type = g_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                    print('jestem5')

                elif x.coefficient[i] == TWOEL_INT_TRANS:
                    v_type = x.int_type("t", i)
                    if(v_type not in all_int_type):
                        all_int_type.append(v_type)
                    print('all_int_type', all_int_type)

                    idx = x.coefficient_idx[i]
    
                    if (i + 1) < len(x.coefficient):
                        s += ("{v_type}({i0}, {i1}, {i2}, {i3}) * ".format(v_type = v_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                    else :
                        s += ("{v_type}({i0}, {i1}, {i2}, {i3})".format(v_type = v_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                    print('jestem6')
                                        

                elif x.coefficient[i] == F12_TWOEL:
                    print('takjestem2o')
                    g_type = x.int_type("ff", i)
                    print(g_type)
                    if(g_type not in all_int_type):
                        all_int_type.append(g_type)
                    print('all_int_type', all_int_type)
                    idx = x.coefficient_idx[i]
    
                    if (i + 1) < len(x.coefficient):
                        s += ("{g_type}({i0}, {i1}, {i2}, {i3}) * ".format(g_type = g_type, 
                                                                             i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                    else :
                        s += ("{g_type}({i0}, {i1}, {i2}, {i3})".format(g_type = g_type, 
                                                                          i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))


                elif x.coefficient[i] == BARENUCL_HAM_TRANS:
                    idx = x.coefficient_idx[i]
                    k_type = x.int1_type_trans(i)
                    

                    if (i + 1) < len(x.coefficient):
                        s += "{k_type}({i0}, {i1}) * ".format(k_type = k_type, i0 = idx[0], i1 = idx[1])
                    else :
                        s += "{k_type}({i0}, {i1})".format(k_type = k_type, i0 = idx[0], i1 = idx[1])


                elif x.coefficient[i] == FOCK_MATRIX or x.coefficient[i] == FOCK_MATRIX_TRANS:
                    idx = x.coefficient_idx[i]
                    if idx[0] == idx[1]:
                        if (i + 1) < len(x.coefficient):
                            s += "eorb({i0}) * ".format(i0 = idx[0])
                        else :
                            s += "eorb({i0})".format(i0 = idx[0])
                            
                elif x.coefficient[i] in OBSERVABLE_MATRIX:
                    idx = x.coefficient_idx[i]
                    if (i + 1) < len(x.coefficient):
                        s += "Obs({i0}, {i1}) * ".format(i0 = idx[0], i1 = idx[1])
                    else :
                        s += "Obs({i0}, {i1})".format(i0 = idx[0], i1 = idx[1])

                elif x.coefficient[i] == CC_AMPLITUDE:
                    #print(x)

                    ln, at = x.amplitude_type(i)
                    if ln == '3' :
                        temp = "(nocc, nactive, "
                    else:
                        temp = "("
                    for j in x.coefficient_idx[i]:
                        temp = temp + str(j) + ","
                        #print(temp)

                    temp = temp[0:len(temp)-1]
                    if (i + 1) < len(x.coefficient):
                        s += "{at}{temp}) * ".format(at = at, temp = temp)
                    else:
                        s += "{at}{temp})".format(at = at, temp = temp)

                elif x.coefficient[i] == EOM_CC_AMPLITUDE_R:
                    ln, at = x.amplitude_type(i)
                    # if ln == '3' :
                    #     temp = "(w, vrdav, t2, t1, nocc, nactive, eorb, "
                    # else:
                    temp = "(vrdav(:, k1), "
                    for j in x.coefficient_idx[i]:
                        temp = temp + str(j) + ","
                    temp = temp[0:len(temp)-1]
                    if (i + 1) < len(x.coefficient):
                        s += "{at}{temp}) * ".format(at = at, temp = temp)
                    else:
                        s += "{at}{temp})".format(at = at, temp = temp)

                elif x.coefficient[i] == EOM_CC_SINGLE_Rl:
                    ln, at = x.amplitude_type(i)
                    # if ln == '3' :
                    #     temp = "(w, vrdav, t2, t1, nocc, nactive, eorb, "
                    # else:
                    if ln == '3':
                        temp = "(r_idxl, vrdav_Rl, "
                    else:
                        temp = "(vrdav_Rl, "
                    for j in x.coefficient_idx[i]:
                        temp = temp + str(j) + ","
                    temp = temp[0:len(temp)-1]
                    if (i + 1) < len(x.coefficient):
                        s += "{at}{temp}) * ".format(at = at, temp = temp)
                    else:
                        s += "{at}{temp})".format(at = at, temp = temp)

                elif x.coefficient[i] == EOM_CC_SINGLE_Rl_plus:
                    ln, at = x.amplitude_type(i)
                    at += 'p'
                    if ln == '3':
                        temp = "(r_idxl, vrdav_Rl, "
                    else:
                        temp = "(vrdav_Rl, "
                    for j in x.coefficient_idx[i]:
                        temp = temp + str(j) + ","
                    temp = temp[0:len(temp)-1]
                    if (i + 1) < len(x.coefficient):
                        s += "{at}{temp}) * ".format(at = at, temp = temp)
                    else:
                        s += "{at}{temp})".format(at = at, temp = temp)

                elif x.coefficient[i] == EOM_CC_SINGLE_Rl_minus:
                    ln, at = x.amplitude_type(i)
                    at += 'm'
                    if ln == '3':
                        temp = "(r_idxl, vrdav_Rl, "
                    else:
                        temp = "(vrdav_Rl, "
                    for j in x.coefficient_idx[i]:
                        temp = temp + str(j) + ","
                    temp = temp[0:len(temp)-1]
                    if (i + 1) < len(x.coefficient):
                        s += "{at}{temp}) * ".format(at = at, temp = temp)
                    else:
                        s += "{at}{temp})".format(at = at, temp = temp)

                elif x.coefficient[i] == EOM_CC_SINGLE_Rr_plus:
                    ln, at = x.amplitude_type(i)
                    at += 'p'
                    if ln == '3':
                        temp = "(r_idxr, vrdav_Rr, "
                    else:
                        temp = "(vrdav_Rr, "
                    for j in x.coefficient_idx[i]:
                        temp = temp + str(j) + ","
                    temp = temp[0:len(temp)-1]
                    if (i + 1) < len(x.coefficient):
                        s += "{at}{temp}) * ".format(at = at, temp = temp)
                    else:
                        s += "{at}{temp})".format(at = at, temp = temp)

                elif x.coefficient[i] == EOM_CC_SINGLE_Rr_minus:
                    ln, at = x.amplitude_type(i)
                    at += 'm'
                    if ln == '3':
                        temp = "(r_idxr, vrdav_Rr, "
                    else:
                        temp = "(vrdav_Rr, "
                    for j in x.coefficient_idx[i]:
                        temp = temp + str(j) + ","
                    temp = temp[0:len(temp)-1]
                    if (i + 1) < len(x.coefficient):
                        s += "{at}{temp}) * ".format(at = at, temp = temp)
                    else:
                        s += "{at}{temp})".format(at = at, temp = temp)


                elif x.coefficient[i] == EOM_CC_SINGLE_Rr:
                    ln, at = x.amplitude_type(i)
                    # if ln == '3' :
                    #     temp = "(w, vrdav, t2, t1, nocc, nactive, eorb, "
                    # else:
                    if ln == '3' :
                        temp = "(r_idxr, vrdav_Rr, "
                    else:
                        temp = "(vrdav_Rr, "
                    for j in x.coefficient_idx[i]:
                        temp = temp + str(j) + ","
                    temp = temp[0:len(temp)-1]
                    if (i + 1) < len(x.coefficient):
                        s += "{at}{temp}) * ".format(at = at, temp = temp)
                    else:
                        s += "{at}{temp})".format(at = at, temp = temp)

                elif x.coefficient[i] == EOM_CC_SINGLE_Ll:
                    ln, at = x.amplitude_type(i)
                    # if ln == '3' :
                    #     temp = "(w, vrdav, t2, t1, nocc, nactive, eorb, "
                    # else:
                    temp = "(vldav_Ll, "
                    for j in x.coefficient_idx[i]:
                        temp = temp + str(j) + ","
                    temp = temp[0:len(temp)-1]
                    if (i + 1) < len(x.coefficient):
                        s += "{at}{temp}) * ".format(at = at, temp = temp)
                    else:
                        s += "{at}{temp})".format(at = at, temp = temp)

                elif x.coefficient[i] == EOM_CC_SINGLE_Lr:
                    ln, at = x.amplitude_type(i)
                    # if ln == '3' :
                    #     temp = "(w, vrdav, t2, t1, nocc, nactive, eorb, "
                    # else:
                    temp = "(vldav_Lr, "
                    for j in x.coefficient_idx[i]:
                        temp = temp + str(j) + ","
                    temp = temp[0:len(temp)-1]
                    if (i + 1) < len(x.coefficient):
                        s += "{at}{temp}) * ".format(at = at, temp = temp)
                    else:
                        s += "{at}{temp})".format(at = at, temp = temp)

                elif x.coefficient[i] == EOM_CC_AMPLITUDE_L:
                    ln, at = x.amplitude_type(i)
                    # if ln == '3' :
                    #     temp = "(k1, vldav, t2, t1, nocc, nactive, eorb, "
                    # else:
                    temp = "(vldav(:, k1), "
                    for j in x.coefficient_idx[i]:
                        temp = temp + str(j) + ","
                    temp = temp[0:len(temp)-1]
                    if (i + 1) < len(x.coefficient):
                        s += "{at}{temp}) * ".format(at = at, temp = temp)
                    else:
                        s += "{at}{temp})".format(at = at, temp = temp)

                elif x.coefficient[i] == S_AMPLITUDE:
                    ln, at = x.amplitude_type(i)
                    if ln == '3' :
                        temp = "(nocc, nactive, "
                    else:
                        temp = "("
                    for j in x.coefficient_idx[i]:
                        temp = temp + str(j) + ","
                    temp = temp[0:len(temp)-1]
                    if (i + 1) < len(x.coefficient):
                        s += "{at}{temp}) * ".format(at = at, temp = temp)
                    else:
                        s += "{at}{temp})".format(at = at, temp = temp)

                elif x.coefficient[i] == INTERM_Ft_F12:
                    at = "Ft_interm"
                    temp = "("
                    for k in x.coefficient_idx[i]:
                        temp = temp + str(k) + ","
                    temp = temp[0:len(temp)-1]
                    if (i + 1) < len(x.coefficient):
                        s += "{at}{temp}) * ".format(at = at, temp = temp)
                    else:
                        s += "{at}{temp})".format(at = at, temp = temp)

                elif x.coefficient[i] == INTERM_V_F12:
                    v_type = x.int_type("vvf", i)
                    if(v_type not in all_int_type):
                        all_int_type.append(v_type)
                    print('all_int_type', all_int_type)


                    idx = x.coefficient_idx[i]
                    if (i + 1) < len(x.coefficient):
                        s += ("{v_type}({i0}, {i1}, {i2}, {i3}) * ".format(v_type = v_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                    else :
                        s += ("{v_type}({i0}, {i1}, {i2}, {i3})".format(v_type = v_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                    # at = "V_interm"                    
                    # temp = "("
                    # for k in x.coefficient_idx[i]:
                    #     temp = temp + str(k) + ","
                    # temp = temp[0:len(temp)-1]
                    # if (i + 1) < len(x.coefficient):
                    #     s += "{at}{temp}) * ".format(at = at, temp = temp)
                    # else:
                    #     s += "{at}{temp})".format(at = at, temp = temp)

                elif x.coefficient[i] == INTERM_X_F12:
                    v_type = x.int_type("xff", i)
                    if(v_type not in all_int_type):
                        all_int_type.append(v_type)

                    idx = x.coefficient_idx[i]
                    if (i + 1) < len(x.coefficient):
                        s += ("{v_type}({i0}, {i1}, {i2}, {i3}) * ".format(v_type = v_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                    else :
                        s += ("{v_type}({i0}, {i1}, {i2}, {i3})".format(v_type = v_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))

                    # at = "X_interm"
                    # temp = "("
                    # for k in x.coefficient_idx[i]:
                    #     temp = temp + str(k) + ","
                    # temp = temp[0:len(temp)-1]
                    # if (i + 1) < len(x.coefficient):
                    #     s += "{at}{temp}) * ".format(at = at, temp = temp)
                    # else:
                    #     s += "{at}{temp})".format(at = at, temp = temp)
                    # print('jestem7')


                elif x.coefficient[i] == INTERM_B_F12:
                    v_type = x.int_type("bff", i)
                    if(v_type not in all_int_type):
                        all_int_type.append(v_type)

                    idx = x.coefficient_idx[i]
                    if (i + 1) < len(x.coefficient):
                        s += ("{v_type}({i0}, {i1}, {i2}, {i3}) * ".format(v_type = v_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                    else :
                        s += ("{v_type}({i0}, {i1}, {i2}, {i3})".format(v_type = v_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))

                    # at = "B_interm"
                    # temp = "("
                    # for k in i.coefficient_idx[i]:
                    #     temp = temp + str(k) + ","
                    # temp = temp[0:len(temp)-1]
                    # if (i + 1) < len(x.coefficient):
                    #     s += "{at}{temp}) * ".format(at = at, temp = temp)
                    # else:
                    #     s += "{at}{temp})".format(at = at, temp = temp)
                    # print('jestem8')


                elif x.coefficient[i] == INTERM_P_F12:
                    v_type = x.int_type("pfvf", i)
                    if(v_type not in all_int_type):
                        all_int_type.append(v_type)

                    idx = x.coefficient_idx[i]
                    if (i + 1) < len(x.coefficient):
                        s += ("{v_type}({i0}, {i1}, {i2}, {i3}) * ".format(v_type = v_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                    else :
                        s += ("{v_type}({i0}, {i1}, {i2}, {i3})".format(v_type = v_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))

                    # at = "P_interm"
                    # temp = "("
                    # for k in x.coefficient_idx[i]:
                    #     temp = temp + str(k) + ","
                    # temp = temp[0:len(temp)-1]
                    # if (i + 1) < len(interm[i].coefficient):
                    #     s += "{at}{temp}) * ".format(at = at, temp = temp)
                    # else:
                    #     s += "{at}{temp})".format(at = at, temp = temp)
                    # print('jestem9')


                elif x.coefficient[i] == F12_AMPLITUDE:
                    at = "t2f"
                    temp = "("
                    for k in x.coefficient_idx[i]:
                        temp = temp + str(k) + ","
                    temp = temp[0:len(temp)-1]
                    if (i + 1) < len(x.coefficient):
                        s += "{at}{temp}) * ".format(at = at, temp = temp)
                    else:
                        s += "{at}{temp})".format(at = at, temp = temp)
                    print('jestem10')
                                        

                else:
                    if len(x.coefficient_idx[i]) != 0:
                        if f12 == True:
                            temp = "{name}_f12(".format(name=x.coefficient[i])
                        else:
                            temp = "{name}(".format(name=x.coefficient[i])
                        print(x.coefficient_idx[i])
                        for j in x.coefficient_idx[i]:
                            temp += "{j},".format(j=j)
                        temp = temp[0:len(temp)-1]
                        print(temp)
                        #temp += "({a} - nvirt0) * nocc + ({k} - nocc0) + 1".format(a = x.coefficient_idx[i][0], k = x.coefficient_idx[i][1])
                        # temp += "({a} - nvirt0) * nocc + ({i} - nocc0) + 1, ({b} - nvirt0) * nocc + ({j} - nocc0) + 1 , & \n \
                        # ({c} - nvirt0) * nocc + ({k} - nocc0) + 1 " \
                        #     .format(a = x.coefficient_idx[i][0], i = x.coefficient_idx[i][1], \
                        #             b = x.coefficient_idx[i][2], j = x.coefficient_idx[i][3], \
                        #             c = x.coefficient_idx[i][4], k = x.coefficient_idx[i][5])
                        #print(temp)
                        if (i + 1) < len(x.coefficient):
                            s += "{temp}) * ".format(temp=temp)
                        else:
                            s += "{temp})".format(temp=temp)
                    else:
                        temp = "{name}".format(name=x.coefficient[i])
                        if (i + 1) < len(x.coefficient):
                            s += "{temp} * ".format(temp=temp)
                        else:
                            s += "{temp}".format(temp=temp)


            s += "\n"
        print('jestem11')

        for i in x.summation:
            s += "end do \n"
        s += "\n"

        # for x in cluster:
        #     if all_interm_dict != {}:
        if f12 == True:
            s += add_close_disk_f12(cluster, all_interm_dict, int_app)

        if all_interm_dict != {}:
            if f12 == False:
                s += add_close_disk(x, all_interm_dict)
        print('jestem12')
        zz = 0
        print(len(cluster))
        for x in cluster:
            print(z, x)
            z += 1
            y2 += 1

            #
            # Alternatywny kod do przetestowania
            #
            #print('PRZETETSTUJ ALTERNATYWNY KOD ARSTOFORT')
            #sys.exit(0)
            if x.num_factor < 0:
                sss = "term({y2}) = term({y2}) * ({nf}d+0) \n".format(y2 = y2, nf = x.num_factor)
            elif x.num_factor > 0:
                if abs((abs(x.num_factor) - 1.0)) < 1.e-10:
                    sss = ""
                else:
                    sss = "term({y2}) = term({y2}) * ({nf}d+0) \n".format(y2 = y2, nf = x.num_factor)


            print('jestem13')
                
            # if x.num_factor < 0:
            #     print('mniejsze')
            #     if abs((abs(x.num_factor) - 1.0)) < 1.e-10:
            #         print('a')
            #         sss = "term({y2}) = -term({y2}) \n".format(y2 = y2)
            #     else:
            #         print('b')
            #         sss = "term({y2}) = term({y2}) * ({nf}d+0) \n".format(y2 = y2, nf = x.num_factor)
            # else:
            #     print('wieksze', x.num_factor - 1.0)
            #     if (abs(x.num_factor - 1.0)) > 1.e-10:
            #         print('a')
            #         sss = "term({y2}) = term({y2}) * {nf}d+0 \n".format(y2 = y2, nf = x.num_factor)
            #     elif (x.num_factor - 1.0) < 0.0:
            #         print('b')
            #         if(abs(x.num_factor - 1.0)) > 1.e-10:
            #             print('c')
           
            s += sss
        print('jestem14')

        s += "\n"
        print('jestem15')

        print('y2', y2)
        print('laaaaaaaaaaaaaaaaaaa')
        # print(s)

        print('KONCZE')
    print('all int types')
    for x in all_int_type:
        print(x)
    print('')

    print('koncze na dobre')
    return s, y2, nocc, nactive


def add_close_disk(x, all_interm_dict):
    s = ""

    if all_interm_dict != {}:
        for coef in x.coefficient:
            if 'interm' in coef:
                coef_name = extract_coef(coef)
                if all_interm_dict[coef_name]['mem'] == 'disk':
                    s += """close(u{coef_name}_pt{pt}) \n""".format(coef_name=coef_name, pt = all_interm_dict[coef_name]['pt'])
                    s += "deallocate ({coef_name}) \n".format(coef_name=coef_name)
    return s


def add_close_disk_f12(cluster, all_interm_dict, int_app):
    s = ""

    coef_name_list = []
    for x in cluster:
        for coef in x.coefficient:
            if 'l' in coef:
                coef_name = extract_coef(coef)
                if all_interm_dict[coef_name]['mem'] == 'disk':
                    if coef_name not in coef_name_list:
                        coef_name_list.append(coef_name)

    for coef_name in coef_name_list:
        s += """close(u{coef_name}_{int_app}) \n""".format(coef_name=coef_name, int_app=int_app)
        s += "deallocate ({coef_name}_{int_app}) \n".format(coef_name=coef_name, int_app=int_app)
    return s


# def add_close_disk_f12(x, all_interm_dict):
#     s = ""

#     if all_interm_dict != {}:
#         for coef in x.coefficient:
#             if 'l' in coef:
#                 coef_name = extract_coef(coef)
#                 if all_interm_dict[coef_name]['mem'] == 'disk':
#                     s += """close(u{coef_name}) \n""".format(coef_name=coef_name)
#                     s += "deallocate ({coef_name}) \n".format(coef_name=coef_name)
#     return s

def add_read_disk(x, all_interm_dict):


    s = ""
    
    if all_interm_dict != {}:
        nn = 0
        print('xx', x)
        for coef in x.coefficient:
            if 'interm' in coef:
                coef_name = extract_coef(coef)
                if all_interm_dict[coef_name]['mem'] == 'disk':
                    s += """open(newunit=u{coef_name}_pt{pt}, file='{filename}_pt{pt}', access='stream',&                                        
 form='unformatted', status='unknown') \n""".format(filename = all_interm_dict[coef_name]['filename'], coef_name=coef_name, pt = all_interm_dict[coef_name]['pt'])
                    s += "allocate ({coef_name}(".format(coef_name=coef_name)
                    for j in all_interm_dict[coef_name]['if']:
                        if j in virt_all:
                            s += "nocc+1: nactive, "
                        elif j in occ_all:
                            s += "1: nocc, "
                    s = s[0:len(s)-2]
                    s += "))\n"

                    hh = 0
                    gx =""
                    for idx in x.coefficient_idx[nn]:
                        # if all_interm_dict[coef_name]['if'][hh] in virtual:
                        #     s += "do {idx} = nocc + 1, nactive \n".format(idx = idx)
                        # elif all_interm_dict[coef_name]['if'][hh] in occupied:
                        #     s += "do {idx} = 1, nocc \n".format(idx = idx)
                        hh += 1
                        gx += idx + ','
                    gx = gx[0:len(gx)-2]
                    # s += 'read (u{coef_name}_pt{pt}) {coef_name}({gx})\n'.format(coef_name=coef_name, pt = all_interm_dict[coef_name][\
                    #         'pt'], gx=gx)
                    s += 'read (u{coef_name}_pt{pt}) {coef_name}\n'.format(coef_name=coef_name, pt = all_interm_dict[coef_name][\
                            'pt'], gx=gx)

                    # for idx in x.coefficient_idx[nn]:
                    #     s += "end do \n"
                    nn += 1

    return s


def add_read_disk_f12(cluster, all_interm_dict, int_app):

    
    s = ""

    coef_name_list = []
    coef_name_list2 = []
    for x in cluster:
        print('jadajadajad', x)
        j = 0
        for coef in x.coefficient:
            j += 1
            if 'l' in coef:
                coef_name = extract_coef(coef)

                if all_interm_dict[coef_name]['mem'] == 'disk':
                    print('name-disk', all_interm_dict[coef_name], x, coef_name, coef, j)
                    if coef_name not in coef_name_list:
                        coef_name_list.append(coef_name)
                        coef_name_list2.append(coef)

    for i in range(0, len(coef_name_list)):
        coef_name = coef_name_list[i]
        coef_name2 = coef_name_list2[i]

        s += """open(newunit=u{coef_name}_{int_app}, file='{filename}_f12', access='stream',&                                        
        form='unformatted', status='unknown') \n""".format(filename = all_interm_dict[coef_name]['filename'], coef_name=coef_name, int_app=int_app)
        s += "allocate ({coef_name}_{int_app}(".format(coef_name=coef_name, int_app=int_app)
        for j in all_interm_dict[coef_name]['if']:
            if j in virt_all:
                s += "nocc+1: nactive, "
            elif j in occ_all:
                s += "1: nocc, "
            elif j in CABS_all:
                s += "ncabs0:ncabs1, "
        s = s[0:len(s)-2]
        s += "))\n"
        s += 'read (u{coef_name}_{int_app}) {coef_name}_{int_app}\n'.format(coef_name=coef_name, int_app=int_app)
    return s

# def add_read_disk_f12(x, all_interm_dict):
#     print('all_interm_dict234')
#     print(all_interm_dict)

#     s = ""
    
#     if all_interm_dict != {}:
#         nn = 0
#         for coef in x.coefficient:
#             print(x.coefficient, coef)
#             if 'l' in coef:
#                 coef_name = extract_coef(coef)
#                 print('coeffff name', coef_name)
#                 if all_interm_dict[coef_name]['mem'] == 'disk':
#                     s += """open(newunit=u{coef_name}, file='{filename}', access='stream',&                                        
#  form='unformatted', status='unknown') \n""".format(filename = all_interm_dict[coef_name]['filename'], coef_name=coef_name)
#                     s += "allocate {coef_name}(".format(coef_name=coef_name)
#                     for j in all_interm_dict[coef_name]['if']:
#                         if j in virt_all:
#                             s += "nocc+1: nactive, "
#                         elif j in occ_all:
#                             s += "1: nocc, "
#                         elif j in CABS_all:
#                             s += "ncabs0:ncabs1, "
#                     s = s[0:len(s)-2]
#                     s += "))\n"

#                     hh = 0
#                     gx =""
#                     for idx in x.coefficient_idx[nn]:
#                         hh += 1
#                         gx += idx + ','
#                     gx = gx[0:len(gx)-2]
#                     s += 'read (u{coef_name}) {coef_name}\n'.format(coef_name=coef_name, gx=gx)

#                     nn += 1

#     return s

def extract_coef(coef):

    coef_name =""
    no = "0123456789"
    for s in coef:
        # if s in no:
        #     coef_name += s
        if s == "_":
            break
        else:
            coef_name += s

 

    # coef_name ="interm"
    # no = "0123456789"
    # for s in coef:
    #     if s in no:
    #         coef_name += s
    #     if s == "_":
    #         break

    # return coef
    return coef_name

def cluster_for_fortran_mem2(ars):

    print('CLUSTERY')
    banned = []
    clusters_mem = []
    cluster_name = []

    for x in range(0, len(ars)):
        mini_cluster = arithmetic_string()
        if x not in banned:
            mini_cluster.append(ars[x])
            cluster_name.append(ars[x].name)
            #cluster_name.append(name_list[x])
            #print('cluster name', name_list[x])
#            print('cluster name', ars[x].name)
            for y in range(x +1, len(ars)):
                if y not in banned:
                    if ars[y].summation == ars[x].summation:
                        #if name_list[y] == name_list[x]:
                        if ars[y].name == ars[x].name:
                            if ars[y].constraints  == ars[x].constraints:
                                mini_cluster.append(ars[y])
                                print('dodaje', ars[y])
                                banned.append(y)
            clusters_mem.append(mini_cluster)
    print('clustery')
    for x in clusters_mem:
        for g in x:
            print(g)
        print('zz')
            
    print('')
    print('asdfa')
    return clusters_mem, cluster_name


def cluster_for_fortran_mem3(ars):

    print('CLUSTERY')
    banned = []
    clusters_mem = []
    cluster_name = []

    for x in range(0, len(ars)):
#        print(x)
        mini_cluster = arithmetic_string()
        if x not in banned:
            mini_cluster.append(ars[x])
            cluster_name.append(ars[x].name)
            #cluster_name.append(name_list[x])
            #print('cluster name', name_list[x])
#            print('cluster name', ars[x].name)
            for y in range(x +1, len(ars)):
                if y not in banned:
                    if ars[y].delta == ars[x].delta:
                        #if name_list[y] == name_list[x]:
                        if ars[y].name == ars[x].name:
                            mini_cluster.append(ars[y])
                            print('dodaje', ars[y])
                            banned.append(y)
            clusters_mem.append(mini_cluster)
    print('clustery')
    for x in clusters_mem:
        for g in x:
            print(g)
        print('zz')
            
    print('')
    print('asdfa')
    return clusters_mem, cluster_name


def arstofort_mem(ars, dct, name_lst, BRA, KET):
    """ Translates arithmetic string 
    to fortran code.  Each ugg in arithmetic
    string is translated to set of loops. There is no
    declaration of variables. Returns string 's',
    containg fortran code.
    """
    nocc = False
    nactive = False
    print('arstofort_mem-PLUSZOPLUSZO')

    for x in ars:
        x.optimize()

    s = ''
    cluster_mem, cluster_name = cluster_for_fortran_mem2(ars)

    pluszyna = 1
    print('poczatekplusz')
    for x in range(0, len(cluster_mem)):
        cluster = cluster_mem[x]
        if BRA == 2 or KET == 2:
            for y in cluster:
                y.ibra = ['a', 'i', 'b', 'j']
                w = []
                for f in y.constraints.items():
                    w.append([f[0], list(f[1])])
                if BRA == 2 and KET == 3:
                    ibra0 = ['a', 'i', 'b', 'j']
                elif BRA == 3 and KET ==2:
                    ibra0 = ['d', 'l', 'e', 'm']
                for f in w:
                    for g in f[1]:
                        if g == ibra0[0]:
                            y.ibra[0] = f[0]
                        if g == ibra0[1]:
                            y.ibra[1] = f[0]
                        if g == ibra0[2]:
                            y.ibra[2] = f[0]
                        if g == ibra0[3]:
                            y.ibra[3] = f[0]
                print('2', y, y.loops, y.constraints, y.ibra)
        elif BRA == 1 or KET == 1:
            for y in cluster:
                y.ibra = ['a', 'i']
                w = []
                for f in y.constraints.items():
                    w.append([f[0], list(f[1])])
                if BRA == 1 and KET == 3:
                    ibra0 = ['a', 'i']
                elif BRA == 3 and KET == 1:
                    ibra0 = ['d', 'l']
                for f in w:
                    for g in f[1]:
                        if g == ibra0[0]:
                            y.ibra[0] = f[0]
                        if g == ibra0[1]:
                            y.ibra[1] = f[0]
                print(pluszyna,  y, y.loops, y.constraints, y.ibra)
                pluszyna += 1
        print('')

    print('koniecplusz')
    y2 = 1
    y = -1
    kk = -1
    kkk = 1
    for x in range(0, len(cluster_mem)):
        # cluster = ars.clusters[key]
        cluster = cluster_mem[x]
        n_if = 0
        kk += 1

        for key in cluster[0].constraints:
            if key not in cluster[0].summation:
                limit_lst = list(cluster[0].constraints[key])
                limit_lst.sort()
                limit = ''
                for l in limit_lst:
                    limit += l
                s += "if ({a} >= {a0} .and. {a} <= {a1}) then \n".format(a = key, a0 = 'n0' + limit, a1 = 'n1'+ limit)
                n_if += 1

        n = 1
        s1 = ""
        for j in cluster[0].summation:
            s1 += "{j}, ".format(j=j)
            n += 1
        s1 += "p, idx"
        s += "!{nm} \n".format(nm=cluster_name[kk])
#        s += "ai = ({n} - nvirt0) * nocc + ({m} - nocc0) + 1 \n".format(n=cluster_name[kk][6], m = cluster_name[kk][7])
#        s += "print*, ai \n"
#        s += "!$omp parallel private({s1})& \n".format(s1=s1)
#        s += "!$omp default(shared)\n"
#        if len(cluster[0].summation) == 0:
#            s += "!$omp do &\n".format(n=n)
#        else:
#            s += "!$omp do collapse({n})&\n".format(n=n)
           #s += "{z} \n".format(z=cluster[0].summation)
#        s += "!$omp reduction (+:{gname}_{sname})\n".format(gname = dct['gname'], sname = dct['sname'])
        nm = []

        # for l in range (dct['n0'], dct['n1']):
        #     if cluster_name[kk][l] in cluster[0].summation:
        #         s += " do {a} = {a0}, {a1} \n".format(a = cluster_name[kk][l], a0 = dct['bound'][l-offset][0], a1 = dct['bound'][l-offset][1])
        #         nm.append(cluster_name[kk][l])
        #         n_do += 1

#        s += "ai = ({n} - nvirt0) * nocc + ({m} - nocc0) + 1 \n".format(n=cluster_name[kk][6], m = cluster_name[kk][7])
#        s += "print*, ai \n"

        n_do = 0
        for i in cluster[0].summation:
            if i in cluster[0].constraints:
                limit_lst = list(cluster[0].constraints[i])
                limit_lst.sort()
                limit = ''
                for l in limit_lst:
                    limit += l
                s += " do {a} = {a0}, {a1} \n".format(a = i, a0 = 'n0' + limit, a1 = 'n1'+ limit)
                n_do += 1
        endif_for_doubles = False
        if BRA == 1 or KET == 1:
            for r in range(0, len(cluster)):
                print(cluster[r].ibra, cluster[0].ibra)
                if cluster[r].ibra != cluster[0].ibra:
                    print('ROZNE')
                    print(cluster[0], cluster[0].constraints, cluster[0].ibra)
                    print(cluster[r], cluster[r].constraints, cluster[r].ibra)
                    sys.exit(0)
            s += "idx = ai_mem({a}, {i})\n".format(a=cluster[0].ibra[0], i=cluster[0].ibra[1])
        elif BRA == 2:
            endif_for_doubles = True
            s+=""" ai = ({a} - nvirt0) * nocc + ({i} - nocc0) + 1 
            bj = ({b} - nvirt0) * nocc + ({j} - nocc0) + 1
            if (ai .ge. bj) then
            idx = npair + ((2 * npair - bj + 2) * (bj - 1)) / 2 + ai - bj + 1 \n""" \
                .format(a=cluster[0].ibra[0], 
                        i=cluster[0].ibra[1], b=cluster[0].ibra[2], j=cluster[0].ibra[3])
        elif KET == 2:
            endif_for_doubles =True
            s+=""" dl = ({d} - nvirt0) * nocc + ({l} - nocc0) + 1
            em = ({e} - nvirt0) * nocc + ({m} - nocc0) + 1
            if (dl .ge. em) then
            idx = npair + ((2 * npair - em + 2) * (em - 1)) / 2 + dl - em + 1 \n""" \
                .format(d = cluster[0].ibra[0], \
                        l=cluster[0].ibra[1], e=cluster[0].ibra[2], m=cluster[0].ibra[3])
        n_do2 = 0
        for i in cluster[0].summation:
            if i not in cluster[0].constraints:
                print('d')
                if i in virt_all and i not in nm:
                    print('d1')
                    s += "do {i} = nocc + 1, nactive \n".format(i = i)
                    nocc = True
                    nactive = True
                    n_do2 += 1
                elif i in occ_all and i not in nm:
                    print('d2')
                    s += "do {i} = 1, nocc \n".format(i = i)
                    nocc = True
                    n_do2 += 1

        s += "do p = 1, ntrial \n"
        y2 = y

        ai_pairs = []

        for x in cluster:
            print('totototot', x)
            s += "\n"
            y += 1
            ai = ""
            # for i in range(0, len(x.coefficient)):
            #     if x.coefficient[i] == EOM_CC_SINGLE_Rl or x.coefficient[i] == EOM_CC_SINGLE_Rl_plus \
            #        or x.coefficient[i] == EOM_CC_SINGLE_Rl_minus or x.coefficient[i] == EOM_CC_SINGLE_Rr_plus \
            #        or x.coefficient[i] == EOM_CC_SINGLE_Rr_minus or x.coefficient[i] == EOM_CC_SINGLE_Rr:
            #         ai = ""
            #         for j in x.coefficient_idx[i]:
            #             ai = ai + str(j)
            
            # if BRA == 1 or KET == 1:
            #     s += "idx = ai_mem({a}, {i})\n".format(a=ai[0], i=ai[1])
                        
            # if len(ai) ==0:
            #     s += ""
            # elif len(ai) == 2:
            #     s += "idx = ai_mem({a}, {i})\n".format(a=ai[0], i=ai[1])
            # elif len(ai) == 4:
            #     s+=" ai = ({a} - nvirt0) * nocc + ({i} - nocc0) + 1 \n \
            #     bj = ({b} - nvirt0) * nocc + ({j} - nocc0) + 1 \n \
            #     if (ai .ge. bj) then \n \
            #     ibra = npair + ((2 * npair - bj + 2) * (bj - 1)) / 2 + ai - bj + 1\n".format(a=ars[x].ibra[0], \
            #                                                                                  i=ars[x].ibra[1], b=ars[x].ibra[2], j=ars[x].ibra[3])

            #     s += "idx = aibj_mem({a}, {i}, {b}, {j})\n".format(a=ai[0], i=ai[1], b=ai[2], j=ai[3])
            #     #s+= "if(idx==26) print*, 'idx-32', {kkk}, idx \n".format(kkk=kkk)
            #     kkk += 1
            # else:
            #     print('FORTRAN_CODE.PY: NOT IMPLEMENTED NUMBER OF AI_PAIRS')

            nfa = abs(x.num_factor)-1.0
            st = "" 
            if(len(x.coefficient) == 0):
                s += "{gname}_{sname}(p) = {nf}_F64".format(gname = dct['gname'], sname = dct['sname'], nf = x.num_factor)
            else:
                if (abs(nfa)< 10.0**(-8)):
                    if x.num_factor > 0.0:
                        s += "{gname}_{sname}(p)= {gname}_{sname}(p) + ".format(gname = dct['gname'], sname = dct['sname'], nf = x.num_factor)
                        #st += "\n print*, "
                    elif x.num_factor < 0.0 :
                          s += "{gname}_{sname}(p)= {gname}_{sname}(p) - ".format(gname = dct['gname'], sname = dct['sname'], nf = x.num_factor)
                          #st += "\n print*,-"
                else:
                    if x.num_factor > 0.0:
                        s += "{gname}_{sname}(p)= {gname}_{sname}(p) + {nf}_F64 * ".format(gname = dct['gname'], sname = dct['sname'], nf = x.num_factor)
                        #st += "\n print*, {nf}_F64 * ".format(nf = x.num_factor)
                    elif x.num_factor < 0.0 :
                          s += "{gname}_{sname}(p)= {gname}_{sname}(p) {nf}_F64 * ".format(gname = dct['gname'], sname = dct['sname'], nf = x.num_factor)
                         # st += "\n print*, {nf}_F64 * ".format(nf = x.num_factor)
            s_in = ""
            for i in range(0, len(x.coefficient)):
                
                if x.coefficient[i] == TWOEL_INT:
                    g_type = x.int_type("", i)
                    if(g_type not in all_int_type):
                        all_int_type.append(g_type)

                    idx = x.coefficient_idx[i]
    
                    if (i + 1) < len(x.coefficient):
                        s_in += ("{g_type}({i0}, {i1}, {i2}, {i3}) * ".format(g_type = g_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                    else :
                        s_in += ("{g_type}({i0}, {i1}, {i2}, {i3})".format(g_type = g_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))

                elif x.coefficient[i] == TWOEL_INT_TRANS:
                    v_type = x.int_type("t", i)
                    if(v_type not in all_int_type):
                        all_int_type.append(v_type)

                    idx = x.coefficient_idx[i]
    
                    if (i + 1) < len(x.coefficient):
                        s_in += ("{v_type}({i0}, {i1}, {i2}, {i3}) * ".format(v_type = v_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                    else :
                        s_in += ("{v_type}({i0}, {i1}, {i2}, {i3})".format(v_type = v_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))

                elif x.coefficient[i] == BARENUCL_HAM_TRANS:
                    idx = x.coefficient_idx[i]
                    k_type = x.int1_type_trans(i)

                    if (i + 1) < len(x.coefficient):
                        s_in += "{k_type}({i0}, {i1}) * ".format(k_type = k_type, i0 = idx[0], i1 = idx[1])
                    else :
                        s_in += "{k_type}({i0}, {i1})".format(k_type = k_type, i0 = idx[0], i1 = idx[1])


                elif x.coefficient[i] == FOCK_MATRIX:
                    idx = x.coefficient_idx[i]
                    if idx[0] == idx[1]:
                        if (i + 1) < len(x.coefficient):
                            s_in += "eorb({i0}) * ".format(i0 = idx[0])
                        else :
                            s_in += "eorb({i0})".format(i0 = idx[0])
                            
                elif x.coefficient[i] in OBSERVABLE_MATRIX:
                    idx = x.coefficient_idx[i]
                    if (i + 1) < len(x.coefficient):
                        s_in += "Obs({i0}, {i1}) * ".format(i0 = idx[0], i1 = idx[1])
                    else :
                        s_in += "Obs({i0}, {i1})".format(i0 = idx[0], i1 = idx[1])

                elif x.coefficient[i] == CC_AMPLITUDE:
                    ln, at = x.amplitude_type(i)
                    if ln == '3' :
                        temp = "(nocc, nactive, "
                    else:
                        temp = "("
                    for j in x.coefficient_idx[i]:
                        temp = temp + str(j) + ","
                    temp = temp[0:len(temp)-1]
                    if (i + 1) < len(x.coefficient):
                        s_in += "{at}{temp}) * ".format(at = at, temp = temp)
                    else:
                        s_in += "{at}{temp})".format(at = at, temp = temp)

                elif x.coefficient[i] == EOM_CC_SINGLE_Rl or x.coefficient[i] == EOM_CC_SINGLE_Rl_plus \
                     or x.coefficient[i] == EOM_CC_SINGLE_Rl_minus or x.coefficient[i] == EOM_CC_SINGLE_Rr_plus \
                     or x.coefficient[i] == EOM_CC_SINGLE_Rr_minus or x.coefficient[i] == EOM_CC_SINGLE_Rr:
                    if (i + 1) < len(x.coefficient):
                        s_in += "vdav(p, idx) * "
                    else:
                        s_in += "vdav(p, idx)"

                else:
                    if len(x.coefficient_idx[i]) != 0:
                        temp = "{name}(".format(name=x.coefficient[i])
                        for j in x.coefficient_idx[i]:
                            temp += "{j},".format(j=j)
                        temp = temp[0:len(temp)-2]
                        if (i + 1) < len(x.coefficient):
                            s_in += "{temp}) * ".format(temp=temp)
                        else:
                            s_in += "{temp})".format(temp=temp)
                    else:
                        temp = "{name}".format(name=x.coefficient[i])
                        if (i + 1) < len(x.coefficient):
                            s_in += "{temp} * ".format(temp=temp)
                        else:
                            s_in += "{temp}".format(temp=temp)
            s += s_in
            #st += s_in 
            #st += "\n"
            #s += st
            s += "\n"
        s += "end do \n"
        for i in range(0, n_do2):
                s += "end do \n"
        if endif_for_doubles == True:
            s += "end if \n"
        for i in range(0, n_do):
            s += "end do \n"
        s += "\n"
#        s += "!$omp end do\n"
#        s += "!$omp end parallel \n"
        #        print('n_if', n_if)
        for i in range(0, n_if):
            s += "end if \n"
        s += "\n"

        print(s)


    return s, y2, nocc, nactive



def arstofort_mem_noR(ars, dct, BRA, KET):
    """ Translates arithmetic string 
    to fortran code.  Each ugg in arithmetic
    string is translated to set of loops. There is no
    declaration of variables. Returns string 's',
    containg fortran code.
    """
    nocc = False
    nactive = False

    for x in ars:
        x.optimize()

    s = ''
    # cluster_mem, cluster_name = cluster_for_fortran_mem2(ars)

    # for x in range(0, len(cluster_mem)):
    #     cluster = cluster_mem[x]
    #     for y in cluster:
    #         print(y)
    #     print('')

    y2 = 1
    y = -1
    kk = -1
    print('lenlen')
    print(len(ars))
    for x in range(0, len(ars)):
        print('totototototo', ars[x], ars[x].constraints, ars[x].loops, ars[x].ibra)
        n_if = 0
        kk += 1
        for key in ars[x].constraints:
            if key  in list(dct['zname']):
                limit_lst = list(ars[x].constraints[key])
                limit_lst.sort()
                limit = ''
                for l in limit_lst:
                    limit += l
                s += "if ({a} >= {a0} .and. {a} <= {a1}) then \n".format(a = key, a0 = 'n0' + limit, a1 = 'n1'+ limit)
                n_if += 1

        n = 1
        s1 = ""
        for j in ars[x].loops:
            s1 += "{j}, ".format(j=j)
            n += 1

        s1 += "tmp = zero \n"
        s1 += " do p1 = 1, ntrial"
#        s1 += "kk, ibra, idx"
#        s += "!{nm} \n".format(nm=cluster_name[kk])
#        s += "!$omp parallel private({s1})& \n".format(s1=s1)
#        s += "!$omp default(shared)\n"
#        if len(ars[x].loops) == 0:
#            s += "!$omp do &\n".format(n=n)
#        else:
#            s += "!$omp do collapse({n})&\n".format(n=n)
#        s += "!$omp reduction (+:SIGROWS_CC3)\n"
        nm = []
        n_do = 0
        print(s, ars[x].loops, ars[x].constraints)


        for i in ars[x].loops:
            # if i not in ars[x].constraints:
            #     print('gowno1', i)
            #     if i in virtual and i not in nm:
            #         s += "do {i} = nocc + 1, nactive \n".format(i = i)
            #         nocc = True
            #         nactive = True
            #         n_do += 1
            #     elif i in occupied and i not in nm:
            #         s += "do {i} = 1, nocc \n".format(i = i)
            #         nocc = True
            #         n_do += 1
            # else:
            if i in ars[x].constraints:
                limit_lst = list(ars[x].constraints[i])
                limit_lst.sort()
                limit = ''
                for l in limit_lst:
                    limit += l
                s += " do {a} = {a0}, {a1} \n".format(a = i, a0 = 'n0' + limit, a1 = 'n1'+ limit)
                n_do += 1
        for i in ars[x].summation:
            if i not in ars[x].loops:
                if i not in ars[x].constraints:
                    print(i, 'fff1')
                    if i in virt_all and i not in nm:
                        s += "do {i} = nocc + 1, nactive \n".format(i = i)
                        nocc = True
                        nactive = True
                        n_do += 1
                    elif i in occ_all and i not in nm:
                        s += "do {i} = 1, nocc \n".format(i = i)
                        nocc = True
                        n_do += 1
                else:
                    print(i, 'fff2')
                    limit_lst = list(ars[x].constraints[i])
                    limit_lst.sort()
                    limit = ''
                    for l in limit_lst:
                        limit += l
                    s += " do {a} = {a0}, {a1} \n".format(a = i, a0 = 'n0' + limit, a1 = 'n1'+ limit)
                    n_do += 1

        if BRA == 1:
            s+="ibra = ({a} - nvirt0) * nocc + ({i} - nocc0) + 1 \n".format(a=ars[x].ibra[0], i=ars[x].ibra[1])
        elif BRA == 2:
            s+=" ai = ({a} - nvirt0) * nocc + ({i} - nocc0) + 1 \n \
            bj = ({b} - nvirt0) * nocc + ({j} - nocc0) + 1 \n \
            if (ai .ge. bj) then \n \
            ibra = npair + ((2 * npair - bj + 2) * (bj - 1)) / 2 + ai - bj + 1\n".format(a=ars[x].ibra[0], \
                                                                                         i=ars[x].ibra[1], b=ars[x].ibra[2], j=ars[x].ibra[3])
        elif BRA == 3:
            if KET == 1:
                s += "iket = ({a} - nvirt0) * nocc + ({i} - nocc0) + 1 \n".format(a=ars[x].iket[0], i=ars[x].iket[1])
            elif KET == 2:
                print(ars[x])
                s += """ dl = ({a} - nvirt0) * nocc + ({i} - nocc0) + 1
                em = ({b} - nvirt0) * nocc + ({j} - nocc0) + 1
                if (dl .ge. em) then
                iket = npair + ((2 * npair - em + 2) * (em -1)) / 2 + dl - em +1
                """.format(a = ars[x].iket[0], i=ars[x].iket[1], \
                           b = ars[x].iket[2], j = ars[x].iket[3])

        s += "do p1 = 1, ntrial \n"
        s += ""
        y2 = y

        ai_pairs = []

        #        for x in ars:
        s += "\n"
        y += 1
        ai = ""
        for i in range(0, len(ars[x].coefficient)):
            if ars[x].coefficient[i] == EOM_CC_SINGLE_Rl or ars[x].coefficient[i] == EOM_CC_SINGLE_Rl_plus \
               or ars[x].coefficient[i] == EOM_CC_SINGLE_Rl_minus or ars[x].coefficient[i] == EOM_CC_SINGLE_Rr_plus \
               or ars[x].coefficient[i] == EOM_CC_SINGLE_Rr_minus or ars[x].coefficient[i] == EOM_CC_SINGLE_Rr:
                ai = ""
                for j in ars[x].coefficient_idx[i]:
                    ai = ai + str(j)
                        
        nfa = abs(ars[x].num_factor)-1.0
        st = "" 
        s += "! {x}\n ".format(x=x)
        if(len(ars[x].coefficient) == 0):
            s += "tmp(p1) = (one/(eorb-wr)) * {nf}_F64".format(gname = dct['gname'], sname = dct['sname'], nf = ars[x].num_factor)
        else:
            if (abs(nfa)< 10.0**(-8)):
                if ars[x].num_factor > 0.0:
                    s += "tmp(p1) = (one/(eorb-wr)) * sigma(p1) * ".format(gname = dct['gname'], sname = dct['sname'], nf = ars[x].num_factor)
                    #st += "\n print*, "
                elif ars[x].num_factor < 0.0 :
                    s += "tmp(p1) =  (one/(eorb-wr)) * (- sigma(p1)) * ".format(gname = dct['gname'], sname = dct['sname'], nf = ars[x].num_factor)
                    #st += "\n print*,-"
            else:
                if ars[x].num_factor > 0.0:
                    s += "tmp(p1) = (one/(eorb-wr)) * sigma(p1) * {nf}_F64 * ".format(gname = dct['gname'], sname = dct['sname'], nf = ars[x].num_factor)
                    #st += "\n print*, {nf}_F64 * ".format(nf = ars[x].num_factor)
                elif ars[x].num_factor < 0.0 :
                    s += "tmp(p1) = (one/(eorb-wr)) * sigma(p1) * ({nf}_F64) * ".format(gname = dct['gname'], sname = dct['sname'], nf = ars[x].num_factor)
                    #st += "\n print*, {nf}_F64 * ".format(nf = ars[x].num_factor)
        s_in = ""
        for i in range(0, len(ars[x].coefficient)):
                
            if ars[x].coefficient[i] == TWOEL_INT:
                g_type = ars[x].int_type("", i)
                
                idx = ars[x].coefficient_idx[i]
    
                if (i + 1) < len(ars[x].coefficient):
                    s_in += ("{g_type}({i0}, {i1}, {i2}, {i3}) * ".format(g_type = g_type, 
                              i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                else :
                    s_in += ("{g_type}({i0}, {i1}, {i2}, {i3})".format(g_type = g_type, 
                            i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))

            elif ars[x].coefficient[i] == TWOEL_INT_TRANS:
                v_type = ars[x].int_type("t", i)
                idx = ars[x].coefficient_idx[i]
    
                if (i + 1) < len(ars[x].coefficient):
                    s_in += ("{v_type}({i0}, {i1}, {i2}, {i3}) * ".format(v_type = v_type, 
                        i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))
                else :
                    s_in += ("{v_type}({i0}, {i1}, {i2}, {i3})".format(v_type = v_type, 
                        i0 = idx[0], i1 = idx[1], i2 = idx[2], i3 = idx[3]))

            elif ars[x].coefficient[i] == BARENUCL_HAM_TRANS:
                idx = ars[x].coefficient_idx[i]
                k_type = ars[x].int1_type_trans(i)

                if (i + 1) < len(ars[x].coefficient):
                    s_in += "{k_type}({i0}, {i1}) * ".format(k_type = k_type, i0 = idx[0], i1 = idx[1])
                else :
                    s_in += "{k_type}({i0}, {i1})".format(k_type = k_type, i0 = idx[0], i1 = idx[1])


            elif ars[x].coefficient[i] == FOCK_MATRIX:
                idx = ars[x].coefficient_idx[i]
                if idx[0] == idx[1]:
                    if (i + 1) < len(ars[x].coefficient):
                        s_in += "eorb({i0}) * ".format(i0 = idx[0])
                    else :
                        s_in += "eorb({i0})".format(i0 = idx[0])
                            
            elif ars[x].coefficient[i] in OBSERVABLE_MATRIX:
                idx = ars[x].coefficient_idx[i]
                if (i + 1) < len(ars[x].coefficient):
                    s_in += "Obs({i0}, {i1}) * ".format(i0 = idx[0], i1 = idx[1])
                else :
                    s_in += "Obs({i0}, {i1})".format(i0 = idx[0], i1 = idx[1])

            elif ars[x].coefficient[i] == CC_AMPLITUDE:
                ln, at = ars[x].amplitude_type(i)
                if ln == '3' :
                    temp = "(nocc, nactive, "
                else:
                    temp = "("
                for j in ars[x].coefficient_idx[i]:
                    temp = temp + str(j) + ","
                temp = temp[0:len(temp)-1]
                if (i + 1) < len(ars[x].coefficient):
                    s_in += "{at}{temp}) * ".format(at = at, temp = temp)
                else:
                    s_in += "{at}{temp})".format(at = at, temp = temp)

            elif ars[x].coefficient[i] == EOM_CC_SINGLE_Rl or ars[x].coefficient[i] == EOM_CC_SINGLE_Rl_plus \
                 or ars[x].coefficient[i] == EOM_CC_SINGLE_Rl_minus or ars[x].coefficient[i] == EOM_CC_SINGLE_Rr_plus \
                 or ars[x].coefficient[i] == EOM_CC_SINGLE_Rr_minus or ars[x].coefficient[i] == EOM_CC_SINGLE_Rr:
                if (i + 1) < len(ars[x].coefficient):
                    s_in += "vdav(p, idx) * "
                else:
                    s_in += "vdav(p, idx)"

            else:
                if len(ars[x].coefficient_idx[i]) != 0:
                    temp = "{name}(".format(name=ars[x].coefficient[i])
                    for j in ars[x].coefficient_idx[i]:
                        temp += "{j},".format(j=j)
                    temp = temp[0:len(temp)-2]
                    if (i + 1) < len(ars[x].coefficient):
                        s_in += "{temp}) * ".format(temp=temp)
                    else:
                        s_in += "{temp})".format(temp=temp)
                else:
                    temp = "{name}".format(name=ars[x].coefficient[i])
                    if (i + 1) < len(ars[x].coefficient):
                        s_in += "{temp} * ".format(temp=temp)
                    else:
                        s_in += "{temp}".format(temp=temp)
        s += s_in
        s += "\n"
        if KET == 3:
            s += """do p2 = 1, ntrial
            d_small(p2, p1) = d_small(p2, p1) + vdav(p2, ibra) * tmp(p1)
            end do
            SIGROWS_CC3(p1, ibra) = SIGROWS_CC3(p1, ibra) + tmp(p1)        
            """ 
        if BRA == 3:
            s += """do p2 = 1, ntrial
            d_small(p2, p1) = d_small(p2, p1) + vdav(p2, iket) * tmp(p1)
            end do
            SIGROWS_CC3(p1, iket) = SIGROWS_CC3(p1, iket) + tmp(p1)
            """

        s += "end do \n"
        if BRA == 2 or KET == 2:
            s += "end if \n"
        for i in range(0, n_do):
            s += "end do \n"
        s += "\n"
#        s += "!$omp end do\n"
#        s += "!$omp end parallel \n"
        #print('n_if', n_if)
        for i in range(0, n_if):
            s += "end if \n"
        s += "\n"

    print(s)

    return s, y2, nocc, nactive



def arstofort_matrix_wm(intermediates_dict, all_interm_dict, multiplicity, mbpt, int_app):
    """ Translates arithmetic string 
    to fortran code.  Each ugg in arithmetic
    string is translated to set of loops. There is no
    declaration of variables. Returns string 's',
    containg fortran code.
    """

    nocc = False
    nactive = False

    
    interm = arithmetic_string()
    for i in range(0, len(intermediates_dict)):
        interm.append(intermediates_dict[i]['interm'])

    for i in range(0, len(interm)):
        interm[i].optimize()

    s = ""
    s_init = ""
    s_free = ""

    for i in range(0, len(intermediates_dict)):
        # if multiplicity == 1:
        #     int_name = intermediates_dict[i]['int_name'].format(mbpt=mbpt)
        # elif multiplicity == 3:
        #     int_name = intermediates_dict[i]['int_name'].format(mbpt=mbpt)
        # elif multiplicity == 13:
        #     int_name = intermediates_dict[i]['int_name'].format(mbpt=mbpt)
        # elif multiplicity == 31:
        #     int_name = intermediates_dict[i]['int_name']+'_so_left'+'_pt{mbpt}'.format(mbpt=mbpt)
        int_name = intermediates_dict[i]['int_name']
        if len(intermediates_dict[i]['if']) != 0:
            if intermediates_dict[i]['mem'] == 'ram':
                s_init += "allocate({int_name}(".format(int_name=int_name)
                for j in intermediates_dict[i]['if']:
                    if j in virt_all:
                        s_init += "nocc+1: nactive, "
                    elif j in occ_all:
                        s_init += "1: nocc, "
                s_init = s_init[0:len(s_init)-2]
                s_init += "))\n"

    for i in range(0, len(intermediates_dict)):
        # if multiplicity == 1:
        #     int_name = intermediates_dict[i]['int_name']+'_pt{mbpt}'.format(mbpt=mbpt)
        # elif multiplicity == 3:
        #     int_name = intermediates_dict[i]['int_name']+'_triplet'+'_pt{mbpt}'.format(mbpt=mbpt)
        # elif multiplicity == 13:
        #     int_name = intermediates_dict[i]['int_name']+'_so'+'_pt{mbpt}'.format(mbpt=mbpt)
        # elif multiplicity == 31:
        #     int_name = intermediates_dict[i]['int_name']+'_so_left'+'_pt{mbpt}'.format(mbpt=mbpt)
        if intermediates_dict[i]['mem'] == 'ram':
            int_name = intermediates_dict[i]['int_name']
            s_init += "{int_name} = zero \n".format(int_name=int_name)

    for i in range(0, len(intermediates_dict)):
        # if multiplicity == 1:
        #     int_name = intermediates_dict[i]['int_name']+'_pt{mbpt}'.format(mbpt=mbpt)
        # elif multiplicity == 3:
        #     int_name = intermediates_dict[i]['int_name']+'_triplet'+'_pt{mbpt}'.format(mbpt=mbpt)
        # elif multiplicity == 13:
        #     int_name = intermediates_dict[i]['int_name']+'_so'+'_pt{mbpt}'.format(mbpt=mbpt)
        # elif multiplicity == 31:
        #     int_name = intermediates_dict[i]['int_name']+'_so_left'+'_pt{mbpt}'.format(mbpt=mbpt)
        int_name = intermediates_dict[i]['int_name']
        if intermediates_dict[i]['mem'] == 'ram':
            if len(intermediates_dict[i]['if'])!= 0:
                s_free += "deallocate({int_name})\n".format(int_name=int_name)

    for i in range(0, len(intermediates_dict)):

        int_name = intermediates_dict[i]['int_name'] 
        s1 = ""

        # current intermediate is large, open file to store it on disk for later
        if intermediates_dict[i]['mem'] == 'disk':
            coef_name = intermediates_dict[i]['coef_name']
            s += """open(newunit=u{coef_name}_pt{pt}, file='{filename}_pt{pt}', access='stream',&
 form='unformatted', status='unknown') \n""".format(i=i, \
                                                    no = all_interm_dict[coef_name]['no'], \
                                                    filename = all_interm_dict[coef_name]['filename'],\
                                                    coef_name=coef_name, pt = intermediates_dict[i]['pt'])
            s += "allocate({coef_name}_{int_app}(".format(coef_name=coef_name, int_app=int_app)
            for j in intermediates_dict[i]['if']:
                if j in virt_all:
                    s += "nocc+1: nactive, "
                elif j in occ_all:
                    s += "1: nocc, "
            s = s[0:len(s)-2]
            s += "))\n"

        # current intermediate uses large intermediates that need to read from disk Open file and read it to temp array
        iii = 0
        for gg in range(0, len(intermediates_dict[i]['interm'].coefficient)):
            coef = intermediates_dict[i]['interm'].coefficient[gg]
            if coef in all_interm_dict:
                if all_interm_dict[coef]['mem'] == 'disk':
                    iii += 1
                    s += """open(newunit=u{coef}_pt{pt}, file='{filename}_pt{pt}', access='stream',&                              
 form='unformatted', status='unknown') \n""".format(iii=iii, no = all_interm_dict[coef]['no'], \
                                                    filename = all_interm_dict[coef]['filename'], pt = intermediates_dict[i]['pt'], coef=coef)
                    s += "allocate({int_name}_{int_app}(".format(int_name=coef, int_app=int_app)
                    for j in all_interm_dict[coef]['if']:
                        if j in virt_all:
                            s += "nocc+1: nactive, "
                        elif j in occ_all:
                            s += "1: nocc, "
                    s = s[0:len(s)-2]
                    s += "))\n"
                    
                    gx = ""
                    for gi in intermediates_dict[i]['interm'].coefficient_idx[gg]:
                        # if gi in virtual:
                        #     s += "do {gi} = nocc + 1, nactive \n".format(gi = gi)
                        #     nocc = True
                        #     nactive = True
                        # if gi in occ_all:
                        #     s += "do {gi} = 1, nocc \n".format(gi = gi)
                        #     nocc = True
                        gx += "{gi}, ".format(gi=gi)
                    gx = gx[0:len(gx)-2]
                    
                    s += 'read (u{coef}_pt{pt}) {coef}_{int_app}\n'.format(coef = coef, gx=gx, pt=mbpt, int_app=int_app)
                    # for gi in range(0, len(intermediates_dict[i]['interm'].coefficient_idx[gg])):
                    #     s += "end do \n"

                    
#            s += "read file\n"

#        if intermediates_dict[i]['symmetry'] == False:
        for j in intermediates_dict[i]['interm'].summation:
            s1 += "{j}, ".format(j=j)
        for j in intermediates_dict[i]['idx_fx']:
            s1 += "{j}, ".format(j=j)
        s1 += "sum"
        if (len(intermediates_dict[i]['idx_fx'])!= 0):
            s += "!$omp parallel private({s1})& \n".format(s1=s1)
            s += "!$omp default(shared)\n"
            s += "!$omp do collapse({n})\n".format(n=len(intermediates_dict[i]['idx_fx']))
        for j in intermediates_dict[i]['idx_fx']:
            if j in virt_all:
                s += "do {j} = nocc + 1, nactive \n".format(j = j)
                nocc = True
                nactive = True
            if j in occ_all:
                s += "do {j} = 1, nocc \n".format(j = j)
                nocc = True


        s += "sum = zero \n"
        for j in interm[i].summation:
            if j in virt_all:
                s += "do {j} = nocc + 1, nactive \n".format(j = j)
                nocc = True
                nactive = True
            if j in occ_all:
                s += "do {j} = 1, nocc \n".format(j = j)
                nocc = True

        s += "sum = sum + "  
        for j in range(0, len(interm[i].coefficient)):

            if interm[i].coefficient[j] == CC_AMPLITUDE:
                ln, at = interm[i].amplitude_type(j)
                if ln == '3' :
                    temp = "(nocc, nactive, "
                else:
                    temp = "("
                for k in interm[i].coefficient_idx[j]:
                    temp = temp + str(k) + ","
                temp = temp[0:len(temp)-1]
                if (j + 1) < len(interm[i].coefficient):
                    s += "{at}{temp}) * ".format(at = at, temp = temp)
                else:
                    s += "{at}{temp})".format(at = at, temp = temp)

            elif interm[i].coefficient[j] == EOM_CC_SINGLE_Rl or\
                 interm[i].coefficient[j] == EOM_CC_SINGLE_Rr or\
                 interm[i].coefficient[j] == EOM_CC_AMPLITUDE_R or\
                 interm[i].coefficient[j] == EOM_CC_SINGLE_Ll or \
                 interm[i].coefficient[j] == EOM_CC_SINGLE_Lr or\
                 interm[i].coefficient[j] == EOM_CC_AMPLITUDE_L :

                ln, at = interm[i].amplitude_type(j)

                if interm[i].coefficient[j] == EOM_CC_SINGLE_Rl:
                    if ln == '3':
                        ss = "r_idxl, vrdav_Rl"
                    else:
                        ss = "vrdav_Rl"
                elif interm[i].coefficient[j] == EOM_CC_SINGLE_Rr:
                    if ln == '3':
                        ss = "r_idxr, vrdav_Rr"
                    else:
                        ss = "vrdav_Rr"

                elif interm[i].coefficient[j] == EOM_CC_SINGLE_Ll:
                    ss = "vldav_Ll"
                elif interm[i].coefficient[j] == EOM_CC_SINGLE_Lr:
                    ss = "vldav_Lr"
                elif interm[i].coefficient[j] == EOM_CC_AMPLITUDE_R:
                    ss = "r_idxr, vrdav"
                elif interm[i].coefficient[j] == EOM_CC_AMPLITUDE_L:
                    ss = "vldav"


                temp = """({ss}, """.format(ss = ss)
                for k in interm[i].coefficient_idx[j]:
                    temp = temp + str(k) + ","
                temp = temp[0:len(temp)-1]
                if (j + 1) < len(interm[i].coefficient):
                    s += "{at}{temp}) * ".format(at = at, temp = temp)
                else:
                    s += "{at}{temp})".format(at = at, temp = temp)

            elif interm[i].coefficient[j] == EOM_CC_SINGLE_Rl_plus:
                ln, at = interm[i].amplitude_type(j)
                at += 'p'
                temp = "(vrdav_Rl, "
                for k in interm[i].coefficient_idx[j]:
                    temp = temp + str(k) + ","
                temp = temp[0:len(temp)-1]
                if (j + 1) < len(interm[i].coefficient):
                    s += "{at}{temp}) * ".format(at = at, temp = temp)
                else:
                    s += "{at}{temp})".format(at = at, temp = temp)

            elif interm[i].coefficient[j] == EOM_CC_SINGLE_Rl_minus:
                ln, at = interm[i].amplitude_type(j)
                at += 'm'
                temp = "(vrdav_Rl, "
                for k in interm[i].coefficient_idx[j]:
                    temp = temp + str(k) + ","
                temp = temp[0:len(temp)-1]
                if (j + 1) < len(interm[i].coefficient):
                    s += "{at}{temp}) * ".format(at = at, temp = temp)
                else:
                    s += "{at}{temp})".format(at = at, temp = temp)

            elif interm[i].coefficient[j] == EOM_CC_SINGLE_Rr_plus:
                ln, at = interm[i].amplitude_type(j)
                at += 'p'
                temp = "(vrdav_Rr, "
                for k in interm[i].coefficient_idx[j]:
                    temp = temp + str(k) + ","
                temp = temp[0:len(temp)-1]
                if (j + 1) < len(interm[i].coefficient):
                    s += "{at}{temp}) * ".format(at = at, temp = temp)
                else:
                    s += "{at}{temp})".format(at = at, temp = temp)

            elif interm[i].coefficient[j] == EOM_CC_SINGLE_Rr_minus:
                ln, at = interm[i].amplitude_type(j)
                at += 'm'
                if ln == '3':
                    temp = "(r_idxr, vrdav_Rr, "
                else:
                    temp = "(vrdav_Rr, "
                for k in interm[i].coefficient_idx[j]:
                    temp = temp + str(k) + ","
                temp = temp[0:len(temp)-1]
                if (j + 1) < len(interm[i].coefficient):
                    s += "{at}{temp}) * ".format(at = at, temp = temp)
                else:
                    s += "{at}{temp})".format(at = at, temp = temp)

            elif interm[i].coefficient[j] == S_AMPLITUDE:
                ln, at = interm[i].amplitude_type(j)
                if ln == '3' :
                    temp = "(nocc, nactive, "
                else:
                    temp = "("
                for k in interm[i].coefficient_idx[j]:
                    temp = temp + str(k) + ","
                temp = temp[0:len(temp)-1]
                if (j + 1) < len(interm[i].coefficient):
                    s += "{at}{temp}) * ".format(at = at, temp = temp)
                else:
                    s += "{at}{temp})".format(at = at, temp = temp)
            else:
                temp = "("
                for k in interm[i].coefficient_idx[j]:
                    temp = temp + str(k) + ","
                temp = temp[0:len(temp)-1]
                if (j + 1) < len(interm[i].coefficient):
                    s += "{coef}_{int_app}{temp}) * ".format(coef = interm[i].coefficient[j], int_app = int_app, temp = temp)
                else:
                    s += "{coef}_{int_app}{temp})".format(coef = interm[i].coefficient[j], int_app = int_app, temp = temp)
                
        s += "\n"

        for j in range(0, len(interm[i].summation)):
            s += "end do \n"

        if len(intermediates_dict[i]['if']) != 0:
            print(1, int_name)
            m = "{int_name}(".format(int_name=int_name)
            for k in intermediates_dict[i]['if']:
                m += "{k}, ".format(k=k)
            m = m[0:len(m)-2]
            m += ")"
        else:
            print(2, int_name)
            m = "{int_name}".format(int_name=int_name)


        m += " = + sum \n"
        s += m
       

        for j in range(0, len(intermediates_dict[i]['if'])):
            s += "end do \n"
        if (len(intermediates_dict[i]['if'])!= 0):
            s += "!$omp end do nowait \n"
            s += "!$omp end parallel \n"
            s += "\n"

        if intermediates_dict[i]['mem'] == 'disk':

            # for j in intermediates_dict[i]['if']:
            #     if j in virtual:
            #         s += "do {j} = nocc + 1, nactive \n".format(j = j)
            #     if j in occ_all:
            #         s += "do {j} = 1, nocc \n".format(j = j)

            coef_name = intermediates_dict[i]['coef_name']
            gx = ""
            for rr in intermediates_dict[i]['if']:
                gx += rr + ", "
            gx = gx[0:len(gx)-2]
#            s+= "write(u{coef_name}_pt{pt}) {coef_name}_{int_app}({gx}) \n".format(coef_name=coef_name, gx=gx, pt=mbpt, int_app = int_app)
            s+= "write(u{coef_name}_pt{pt}) {coef_name}_{int_app} \n".format(coef_name=coef_name, gx=gx, pt=mbpt, int_app = int_app)
            # for j in range(0, len(intermediates_dict[i]['if'])):
            #     s += "end do \n"
            s += "deallocate({coef_name}_{int_app}) \n".format(coef_name=coef_name, int_app=int_app)
            s += 'close(u{coef_name}_pt{pt}) \n'.format(coef_name=coef_name, pt=mbpt)

        iii = 0
        for gg in range(0, len(intermediates_dict[i]['interm'].coefficient)):
            coef = intermediates_dict[i]['interm'].coefficient[gg]
            if coef in all_interm_dict:
                if all_interm_dict[coef]['mem'] == 'disk':
                    iii += 1
                    s += """close(u{coef}_pt{pt}) \n""".format(coef=coef, pt=mbpt)
                    s += "deallocate({coef}_{int_app}) \n".format(coef=coef, pt=mbpt, int_app=int_app)



    return s_init, s_free, s

def delta_to_dict(delta_list):
    """ Creates dct of delta 
    substitutions.
    """
    delta_subst = {}
    used = []
    for i in delta_list:
        if i[0] not in used:
            delta_subst[i[0]] = i[0]
            used.append(i[0])
        if i[1] not in used:
            delta_subst[i[1]] = i[0]
            used.append(i[1])
    return delta_subst
    
def delta_combinations(fixed):
    """ From arithmetic string, takes all delta's
    from all ugg's. From list of all delta's 
    e.g. [\delta_{ab}, delta{ij}, \delta_{cd}]
    creates all possible combinations of all possible 
    lenghts.
    """
    delta_list = []

    virt = []
    occ = []
    print('pp')
    for i in fixed:
        if i in virt_all:
            virt.append(i)
        elif i in occ_all:
            occ.append(i)
    print('pp')
    cvirt = list(combinations(virt,2))
    cocc  = list(combinations(occ,2))
    delta_list = cvirt + cocc
    for i in range(0, len(delta_list)):
        delta_list[i] = list(delta_list[i])

    delta_combination_list = []

    print('ppz')
    print(delta_list)

    for i in range(1, len(delta_list)+1):
        temp = combinations(delta_list, i)
        for j in temp:
            delta_combination_list.append(j)
    print('srsr')
    # for x in delta_combination_list:
    #     print(x)
    new_delta_combination_list = disambiguate_delta_combination_list(delta_combination_list)
    # print('ol')
    new_delta_combination_list.append([])
        
    return new_delta_combination_list


def disambiguate_delta_combination_list(delta_combination_list):
    """ Disambiguates delta_combination_list
    such that when there is element \delta_{ab}\delta{ad}\delta{bd}
    and also element \delta{ab}\delta{ad}, the list is
    shortened because these elements are the same.
    """
    new_delta_combination_list = []
    name_list_temp = []
    

    for i in range(0, len(delta_combination_list)):
        temp = deepcopy(delta_combination_list[i])
        new_delta_temp = deepcopy(execute_delta_in_delta(temp))

        if new_delta_temp not in new_delta_combination_list:

            new_delta_combination_list.append(deepcopy(new_delta_temp))

    return new_delta_combination_list

def execute_delta_in_delta(delta_old):

    """ Executes the delta function in given 
    delta string.
    """

    delta = []
    for i in range(0, len(delta_old)):
        delta_old[i].sort()
        delta.append(delta_old[i])
    delta.sort()

    new_delta = []
    for i in range(0, len(delta)):
        for j in range(0, len(delta)):
            if j != i:
                for k in range(0, 2):
                    if delta[j][k] == delta[i][1]:
                        delta[j][k] = delta[i][0]
                delta[j].sort()


    for i in range(0, len(delta)):
        if delta[i][0] != delta[i][1]:
            new_delta.append(delta[i])
    
        
    new_delta.sort()

    return new_delta
    
def rename_name(name_org, delta_list):
    """executes delta in name of
    the function i.e. name = 'aibjck'
    delta_list = ['b', 'a'],
    return aiajck
    """
    name = deepcopy(name_org)
    dl = deepcopy(delta_list)

    for k in range(0, len(dl)):
        dl[k].sort()
        for j in range(0, len(name)):
            if name[j] == dl[k][1]:
                name = name[0:j] + dl[k][0] + name[j + 1:]

    return name


def summation_indices(ars):
    """ Returns list of all summation indices
    in whole arithmetic string.
    """
    
    s_ind_lst = []
    
    for x in ars:
        for j in x.summation:
            if j not in s_ind_lst:
                s_ind_lst.append(j)

    if len(s_ind_lst) > 0:
        s_ind = ","
    else:
        s_ind = ""
    for x in range(0, len(s_ind_lst)):
        if (x + 1) < len(s_ind_lst):
            s_ind = s_ind + str(s_ind_lst[x]) + ','
        else:
            s_ind = s_ind + str(s_ind_lst[x])

    return s_ind

def loops_indices(ars):
    """ Returns list of all summation indices                                                                                                                          
    in whole arithmetic string.                                                                                                                                        
    """

    s_ind_lst = []

    for x in ars:
        for j in x.loops:
                if j not in s_ind_lst:
                    s_ind_lst.append(j)

    if len(s_ind_lst) > 0:
        s_ind = ","
    else:
        s_ind = ""
    for x in range(0, len(s_ind_lst)):
        if (x + 1) < len(s_ind_lst):
            s_ind = s_ind + str(s_ind_lst[x]) + ','
        else:
                s_ind = s_ind + str(s_ind_lst[x])

    return s_ind


def variables(name, ars, nocc, nactive):
    """Returns list of all variables - arguments
    in whole arithmetic string.
    """

    #
    # Find all fixed indices used in current function
    #
    f_ind_lst0 = []

    for x in ars:
        for y in range(0, len(x.coefficient_idx)):
            for z in x.coefficient_idx[y]:
                if z not in x.summation and z not in f_ind_lst0:
                    f_ind_lst0.append(z)
    
    #
    # Arrange used fixed indices in the same order
    # that is in the function name
    #
    f_ind_lst = []
    for x in name:
        if x in f_ind_lst0:
            if x not in f_ind_lst:
                f_ind_lst.append(x)

    
    t1t2_lst = []
    for x in ars:
        for i in range(0, len(x.coefficient)):
            if x.coefficient[i] == CC_AMPLITUDE:
                ln, at = x.amplitude_type(i)            
                if ln != '3':
                    if at not in t1t2_lst:
                        t1t2_lst.append(at)
    nn_lst = []
    if nocc == True:
        nn_lst.append('nocc')
    elif len(t1t2_lst) > 0:
        nn_lst.append('nocc')
    if nactive == True:
        nn_lst.append('nactive')
    elif len(t1t2_lst) > 0:
        nn_lst.append('nactive')

    arguments = ""
    for x in range(0, len(t1t2_lst)):
        if (x + 1) <= len(t1t2_lst):
            arguments = arguments + str(t1t2_lst[x]) + ', '

    for x in range(0, len(nn_lst)):
        if (x + 1) <= len(nn_lst):
            arguments = arguments + str(nn_lst[x]) + ', '
            
    f_ind = ""
    for x in range(0, len(f_ind_lst)):
        if (x + 1) < len(f_ind_lst):
            arguments = arguments + str(f_ind_lst[x]) + ', '
            f_ind = f_ind + str(f_ind_lst[x]) + ', '
        else:
            arguments = arguments + str(f_ind_lst[x])
            f_ind = f_ind + str(f_ind_lst[x])


    return arguments, f_ind, t1t2_lst, nn_lst


def fixed_indices(name, ars):
    """Returns list of all fixed indices
    in whole arithmetic string.
    """
    
    f_ind_lst0 = []

    for x in ars:
        for y in range(0, len(x.coefficient_idx)):
            for z in x.coefficient_idx[y]:
                if z not in x.summation and z not in f_ind_lst0:
                    f_ind_lst0.append(z)
    
    f_ind_lst = []
    for x in name:
        if x in f_ind_lst0:
            if x not in f_ind_lst:
                f_ind_lst.append(x)

    # for x in name:
    #     if x not in f_ind_lst:
    #         f_ind_lst.append(x)
            
    f_ind = ""
    for x in range(0, len(f_ind_lst)):
        if (x + 1) < len(f_ind_lst):
            f_ind = f_ind + str(f_ind_lst[x]) + ','
        else:
            f_ind = f_ind + str(f_ind_lst[x])
            
    return f_ind

def single_function_eom_cc3(dct, t1t2_lst, nn_lst):

    s0 = """
    function {gname}_{sname}({variables}) 
    real(F64) :: {gname}_{sname}  """.format(**dct)
    if len(nn_lst) == 2:
        s1 = """
    integer, intent(in) :: nocc, nactive"""
    elif len(nn_lst) == 1:
        if nn_lst[0] == 'nocc':
            s1 = """
    integer, intent(in) :: nocc"""     
        elif nn_lst['0'] == 'nactive':
            s1 = """
    integer, intent(in) :: nactive"""     
    else:
        s1 = """ """
    if len(t1t2_lst) == 2:
        s2 = """
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 """
    elif len(t1t2_lst) == 1:
        if t1t2_lst[0] == 't2':
            s2 = """
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 """
        elif t1t2_lst[0] == 't1':
            s2 = """
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 """
    else:
        s2 = """ """
    s3 = """
    integer, intent(in) ::  {f_ind} 
    integer :: s {s_ind} 
    real(F64) :: triple_w1
    real(F64) :: triple_w2
    real(F64), dimension(0:{D}) :: term 
    term = 0.d+0 
    triple_w1 = 0.d+0
    triple_w2 = 0.d+0
    {s1_inner}
    {s_inner}
    {gname}_{sname} = triple_w1 + triple_w2
    do s = 0, {D}
    {gname}_{sname} = {gname}_{sname} + term(s)
    end do

    end function {gname}_{sname}
    """.format(**dct)
    s = s0 + s1 + s2 + s3
    return s

def single_function_eom(dct, t1t2_lst, nn_lst):

    s0 = """
    function {gname}_{sname}({variables}) 
    real(F64) :: {gname}_{sname} """.format(**dct)
    if len(nn_lst) == 2:
        s1 = """
    integer, intent(in) :: nocc, nactive"""
    elif len(nn_lst) == 1:
        if nn_lst[0] == 'nocc':
            s1 = """
    integer, intent(in) :: nocc"""     
        elif nn_lst['0'] == 'nactive':
            s1 = """
    integer, intent(in) :: nactive"""     
    else:
        s1 = """ """
    if len(t1t2_lst) == 2:
        s2 = """
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 """
    elif len(t1t2_lst) == 1:
        if t1t2_lst[0] == 't2':
            s2 = """
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 """
        elif t1t2_lst[0] == 't1':
            s2 = """
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 """
    else:
        s2 = """ """
    s3 = """
    integer, intent(in) :: {f_ind} 
    integer :: s {s_ind} 
    real(F64), dimension(0:{D}) :: term 
    term = 0.d+0 
    {s_inner}
    {gname}_{sname} = 0.d+0
    do s = 0, {D}
    {gname}_{sname} = {gname}_{sname} + term(s)
    end do

    end function {gname}_{sname}""".format(**dct)

    s = s0 + s1 + s2 + s3
    return s


def single_function_eom_mem(dct, t1t2_lst, nn_lst, BRA, KET, noR = False):

    if noR == False:
        s0 = """        
        subroutine sub_{gname}_{sname}({gname}_{sname}, {variables}) 
        real(F64), dimension(:), intent(out) :: {gname}_{sname}
        real(F64), dimension(:, :), intent(in) :: vdav
        integer, intent(in) :: nocc0, nvirt0""".format(**dct)
        print('klklkl', dct['variables'])

    else:
        s0 = """                                                                                                     
        subroutine sub_{gname}_{sname}({variables})
        real(F64), dimension(:, :), intent(inout) :: SIGROWS_CC3
        real(F64), dimension(:,:), intent(inout) :: d_small
        real(F64), dimension(:), intent(in) :: sigma
        real(F64), intent(in) :: eorb, wr
        real(F64), dimension(:, :), intent(in) :: vdav
        integer, intent(in) :: nocc0, nvirt0
        real(F64), dimension(:), allocatable :: tmp""".format(**dct)
    s0+="""
        integer, intent(in) :: {bounds_str}
    """.format(**dct)

    if len(nn_lst) == 2:
        s1 = """
    integer, intent(in) :: nocc, nactive"""
    elif len(nn_lst) == 1:
        if nn_lst[0] == 'nocc':
            s1 = """
    integer, intent(in) :: nocc"""     
        elif nn_lst['0'] == 'nactive':
            s1 = """
    integer, intent(in) :: nactive"""     
    else:
        s1 = """ """
    if len(t1t2_lst) == 2:
        s2 = """
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 """
    elif len(t1t2_lst) == 1:
        if t1t2_lst[0] == 't2':
            s2 = """
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 """
        elif t1t2_lst[0] == 't1':
            s2 = """
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 """
    else:
        s2 = """ """
    if noR == False:
        if KET == 1:
            s3 = """                                                                                                                                              
            integer, intent(in) :: ntrial, {f_ind}
            integer :: p, idx {s_ind}
            {gname}_{sname} = ZERO""".format(**dct)
        elif KET == 2:
            s3 = """
            integer, intent(in) :: ntrial, {f_ind}
            integer :: p, idx {s_ind}
            integer :: n0de, n1de, n0lm, n1lm, d, l, e, m
            integer :: dl, em
            n0de = max(n0d, n0e)
            n1de = min(n1d, n1e)
            n0lm = max(n0l, n0m)
            n1lm = min(n1l, n1m)
            {gname}_{sname} = ZERO
                """.format(**dct)

        if BRA == 1:
            s3="""
            integer, intent(in) :: ntrial, {f_ind}
            integer :: p, idx {s_ind}
            {gname}_{sname} = ZERO""".format(**dct)
        elif BRA == 2:
            s3 = """
            integer, intent(in) :: ntrial, {f_ind}
            integer :: p, idx {s_ind}
            integer :: n0ab, n1ab, n0ij, n1ij, a, b, i, j
            integer :: ai, bj
            n0ab = max(n0a, n0b)
            n1ab = min(n1a, n1b)
            n0ij = max(n0i, n0j)
            n1ij = min(n1i, n1j)
            {gname}_{sname} = ZERO
                """.format(**dct)

    else:
        if BRA == 2:
            s3 = """
            integer, intent(in) :: ntrial, npair, {f_ind}
            integer :: p1, p2, ai, bj, ibra {s_ind}
            integer :: n0ab, n1ab, n0ij, n1ij, a, b, i, j
            allocate(tmp(ntrial))
            n0ab = max(n0a, n0b)
            n1ab = min(n1a, n1b)
            n0ij = max(n0i, n0j)
            n1ij = min(n1i, n1j)
            """.format(**dct)
        elif BRA == 1:
            s3 = """                                            
            integer, intent(in) :: ntrial, {f_ind}
            integer :: p1, p2, ai, bj, ibra {s_ind}
            allocate(tmp(ntrial))
            """.format(**dct)
        elif BRA == 3:
            if KET == 1: 
                s3 = """                                                                                                                   
                integer, intent(in) :: ntrial, {f_ind}                
                integer :: p1, p2, dl, em, iket  {s_ind}
                allocate(tmp(ntrial))
                """.format(**dct)
                
            elif KET == 2:
                s3 = """
                integer, intent(in) :: ntrial, npair, nocc, {f_ind}
                integer :: p1, p2, dl, em, iket {s_ind}
                integer :: n0de, n1de, n0lm, n1lm
                allocate(tmp(p1))
                n0de = max(n0d, n0e)
                n1de = min(n1d, n1e)
                n0lm = max(n0l, n0m)
                n1lm = min(n1l, n1m)
                """.format(**dct)

    s3 +="""
    {s_inner}
    deallocate(tmp)
    end subroutine sub_{gname}_{sname}""".format(**dct)

    s = s0 + s1 + s2 + s3
    return s


def single_function_ccsd(dct):

    s = """
    function {gname}(t2, t1, nocc, eorb, nactive, {variables}) 
    real(F64) :: {gname}
    integer, intent(in) :: nocc, nactive
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 
    real(F64), dimension(:), intent(in) :: eorb
    integer, intent(in) :: {variables} 
    integer :: s {s_ind} 
    real(F64), dimension(0:{D}) :: term 
    term = 0.d+0
    {s_inner}
    {gname} = 0.d+0 
    do s = 0, {D}
    {gname} = {gname} + term(s)
    end do

    end function {gname}
    """.format(**dct)
    return s


def single_function_ccsd_f12(dct):

    s = """
    function {name}(t2f, t2, nocc, eorb, nactive, {variables}) 
    {decl_name}
    real(F64) :: {name}
    integer, intent(in) :: nocc, nactive
    real(F64), dimension(nocc, nocc, nocc, nocc), intent(in) :: t2f
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 
    real(F64), dimension(:), intent(in) :: eorb
    integer, intent(in) :: {variables} 
    integer :: s {s_ind}



    real(F64), dimension(0:{D}) :: term 
    term = 0.d+0
    {s_inner}
    {name} = 0.d+0 
    do s = 0, {D}
    {name} = {name} + term(s)
    end do

    end function {name}
    """.format(**dct)
    return s

def single_function_cc3(dct):

    s = """
    function {gname}(t2, t1, eorb, nocc, nactive, {variables}) 
    real(F64) :: {gname}
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 
    real(F64), dimension(:), intent(in) :: eorb
    integer, intent(in) :: nocc, nactive
    integer, intent(in) :: {variables} 
    integer :: s {s_ind} 
    real(F64), dimension(0:{D}) :: term 
    term = 0.d+0 
    {s_inner}
    {gname} = 0.d+0 
    do s = 0, {D}
    {gname} = {gname} + term(s)
    end do

    end function {gname}
    """.format(**dct)
    return s

def single_function_cc3_86(dct):

    s = """
    function {gname}(t2, t1, eorb, nocc, nactive, {variables}) 
    real(F64) :: {gname}
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 
    real(F64), dimension(:), intent(in) :: eorb
    integer, intent(in) :: nocc, nactive
    integer, intent(in) :: {variables} 
    integer :: s {s_ind} 
    real(F64), dimension(0:{D}) :: term 
    term = 0.d+0 
    {s_inner}
    {gname} = 0.d+0 
    do s = 0, {D}
    {gname} = {gname} + term(s)
    end do
    t3 = -t3 / (eorb(a) + eorb(b) + eorb(c) - eorb(i) - eorb(j) - eorb(k))
    end function {gname}
    """.format(**dct)
    return s



def function_t_to_s(ars, rank, name, f):

    ars.establish_fixed()
    gname = '{name}'.format(name = name)
    dictionary = {}
    dictionary['gname'] = gname
    if rank == 1:
        dictionary['variables'] = 'a, i'
    elif rank == 2:
        dictionary['variables'] = 'a, i, b, j'
    elif rank == 3:
        dictionary['variables'] = 'a, i, b, j, c, k'
    dictionary['s_ind'] = summation_indices(ars)
    s, D, nocc, nactive = arstofort(ars)
    dictionary['D'] = D
    dictionary['s_inner'] = s

#    name_to_open =  './'+gname+'.f90'
    
#    f = open(name_to_open, 'w')
    string = single_function_ccsd(dictionary)
    f.write(string)

def single_function_transition_xi(dct):

    s = """
    function {name}(Obs, t2, t1,nocc, nactive, {variables}) 
    real(F64) :: {name}
    real(F64), dimension(:,:) :: Obs
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 
    integer, intent(in) :: nocc, nactive
    integer, intent(in) :: {variables} 
    integer :: s {s_ind} 
    real(F64), dimension(0:{D}) :: term 
    term = 0.d+0 
    {s_inner}
    {name} = 0.d+0 
    do s = 0, {D}
    {name} = {name} + term(s)
    end do

    end function {name}
    """.format(**dct)
    return s

def single_function_transition_xi_triples(dct):

    s = """
    function {name}(Obs, t2, t1,nocc, nactive, eorb, {variables}) 
    real(F64) :: {name}
    real(F64), dimension(:,:) :: Obs
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 
    integer, intent(in) :: nocc, nactive
    integer, intent(in) :: {variables} 
    real(F64), dimension(:), intent(in) :: eorb
    integer :: s {s_ind} 
    real(F64), dimension(0:{D}) :: term 
    term = 0.d+0 
    {s_inner}
    {name} = 0.d+0 
    do s = 0, {D}
    {name} = {name} + term(s)
    end do

    end function {name}
    """.format(**dct)
    return s

def single_function_transition_gamma(dct):

    s = """
    function {name}(Obs, t2, t1, s2, s1, nocc, nactive, {variables}) 
    real(F64) :: {name}
    real(F64), dimension(:,:) :: Obs
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: s2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: s1 
    integer, intent(in) :: nocc, nactive
    integer, intent(in) :: {variables} 
    integer :: s {s_ind} 
    real(F64), dimension(0:{D}) :: term 
    term = 0.d+0 
    {s_inner}
    {name} = 0.d+0 
    do s = 0, {D}
    {name} = {name} + term(s)
    end do

    end function {name}
    """.format(**dct)
    return s

def single_function_transition_overlap(dct):

    s = """
    function {name}(t2, t1, s2, s1, nocc, nactive, vrdav_Rl, vrdav_Rr)
    real(F64) :: {name}
    integer, intent(in) :: nocc, nactive
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: s2
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: s1
    real(F64), dimension(:), intent(in) :: vrdav_Rl                                                                       
    real(F64), dimension(:), intent(in) :: vrdav_Rr 
    integer :: s {s_ind}
    real(F64), dimension(0:{D}) :: term
    term = 0.d+0
    {s_inner}
    {name} = 0.d+0
    do s = 0, {D}
    {name} = {name} + term(s)
    end do                                                                                                  
    end function {name}
    """.format(**dct)
    return s


def single_function_transition_dm(dct):

    s = """
    function {name}(t2, t1, s2, s1, nocc, nactive, {variables}) 
    real(F64) :: {name}
    integer, intent(in) :: nocc, nactive
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: s2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: s1 
    integer, intent(in) :: {variables} 
    integer :: s {s_ind} 
    real(F64), dimension(0:{D}) :: term 
    term = 0.d+0 
    {s_inner}
    {name} = 0.d+0 
    do s = 0, {D}
    {name} = {name} + term(s)
    end do

    end function {name}
    """.format(**dct)
    return s

def single_function_transition_wl(dct):

    s = """
    function {name}(t2, t1, s2, s1, nocc, nactive, vldav, {variables}) 
    real(F64) :: {name}
    integer, intent(in) :: nocc, nactive
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: s2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: s1 
    real(F64), dimension(:), intent(in) :: vldav
    integer, intent(in) :: {variables} 
    integer :: s {s_ind} 
    real(F64), dimension(0:{D}) :: term 
    term = 0.d+0 
    {s_inner}
    {name} = 0.d+0 
    do s = 0, {D}
    {name} = {name} + term(s)
    end do

    end function {name}
    """.format(**dct)
    return s

def single_function_init(dct):
#subroutine wm_{mult}intermediates_{method}_init_pt{mbpt}(nocc, nactive)
    s = """
    subroutine {name}_init(nocc, nactive)
    {uses}
    integer, intent(in) :: nocc
    integer, intent(in) :: nactive
    {s_init}
    end subroutine {name}_init
    """.format(**dct)
    return s

def single_function_f12_init(dct):
#subroutine wm_{mult}intermediates_{method}_init_pt{mbpt}(nocc, nactive)
    s = """
    subroutine {name}_init(nocc, nvirt0, nvirt1, ncabs0, ncabs1)
    {uses}
    integer, intent(in) :: nocc, nvirt0, nvirt1, ncabs0, ncabs1
    {s_init}
    end subroutine {name}_init
    """.format(**dct)
    return s

def single_function_free(dct):

    s = """
    subroutine {name}_free
    {uses}
    {s_free}
    end subroutine {name}_free
    """.format(**dct)
    return s

def single_function_intermediates_wm(dct):

    s = """
    subroutine {name}(t2, t1, s2, s1, nocc, nactive, vrdav_Rl, vrdav_Rr, r_idxl, r_idxr)
!    use {name_decl}
    {uses}
    integer, intent(in) :: nocc, nactive
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: s2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: s1 
    real(F64), dimension(:), intent(in) :: vrdav_Rl
    real(F64), dimension(:), intent(in) :: vrdav_Rr
    integer, dimension(:), intent(in) :: r_idxl, r_idxr
    real(F64) :: sum
    integer :: {s_ind} 

    {s_inner}

    end subroutine {name}
    """.format(**dct)
    return s


def single_function_intermediates_f12_SP(dct):

    s = """
    subroutine {name}(t2f, t2, eorb, nocc, nvirt0, nvirt1, ncabs0, ncabs1)
    {uses}
    integer, intent(in) :: nocc, nvirt0, nvirt1, ncabs0, ncabs1
    real(F64), dimension(nocc, nocc, nocc,nocc), intent(in) :: t2f
    real(F64), dimension(nvirt0:nvirt1,nvirt0:nvirt1,nocc,nocc), intent(in) :: t2 
    real(F64), dimension(1:nvirt1), intent(in) :: eorb
    real(F64) :: sum
    integer :: {s_ind} 

    {s_inner}

    end subroutine {name}
    """.format(**dct)
    return s


def single_function_wm(dct):

    s = """
    function {name}(t2, t1, s2, s1, nocc, nactive, vrdav_Rl, vrdav_Rr, k1, k2, p, q, r_idxl, r_idxr) 
    real(F64) :: {name}
    integer, intent(in) :: nocc, nactive
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: s2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: s1 
    real(F64), dimension(:), intent(in) :: vrdav_Rl
    real(F64), dimension(:), intent(in) :: vrdav_Rr
    integer, dimension(:), intent(in) :: r_idxl, r_idxr
    integer, intent(in) :: k1, k2
    integer, intent(in) :: p, q
    integer :: s {s_ind} 
    real(F64), dimension(0:{D}) :: term 
    term = 0.d+0 

    {s_inner1}
    {name} = 0.d+0 
    do s = 0, {D1}
    {name} = {name} + term(s)
    end do

    term = 0.d+0 
    {s_inner2}
    do s = 0, {D2}
    {name} = {name} + term(s)
    end do

    end function {name}
    """.format(**dct)
    return s


def single_function_wm2(dct):

    s = """
    function {name}(t2, t1, s2, s1, nocc, nactive, vrdav_Rl, vrdav_Rr, k1, k2, p, q, r_idxl, r_idxr) 
    {decl_name}
    real(F64) :: {name}
    integer, intent(in) :: nocc, nactive
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: t2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: t1 
    real(F64), dimension(nocc+1:nactive,nocc+1:nactive,nocc,nocc), intent(in) :: s2 
    real(F64), dimension(nocc+1:nactive,nocc), intent(in)                  :: s1 
    real(F64), dimension(:), intent(in) :: vrdav_Rl
    real(F64), dimension(:), intent(in) :: vrdav_Rr
    integer, dimension(:), intent(in) :: r_idxl, r_idxr
    integer, intent(in) :: k1, k2
    integer, intent(in) :: p, q
    integer :: s {s_ind} 
    real(F64), dimension(0:{D}) :: term 
    term = 0.d+0 

    term = 0.d+0 
    {s_inner}
    {name} = zero
    do s = 0, {D}
    {name} = {name} + term(s)
    end do

    end function {name}
    """.format(**dct)
    return s


def function_template_cisd(ars):
    
    gname = './fortran_codes/cisd/'


def function_template_ccsd(ars, rank):

    gname = './fortran_codes/cc/automatic_t'+str(rank)+'_h'
    dct = {}
    dct['gname'] = gname
    
    if rank == 1:
        dct['variables'] = 'a, i'
    elif rank == 2:
        dct['variables'] = 'a, i, b, j'
    

    dct['s_ind'] = summation_indices(ars)
    s, D, nocc, nactive = arstofort(ars)
    dct['D'] = D
    dct['s_inner'] = s

    name_to_open =  './'+gname+'.f90'

    f = open(name_to_open, 'w')
    string = single_function_ccsd(dct)
    f.write(string)

def function_template_cc3(ars, neq):

    gname = 'equation_'+str(neq)
    dct = {}
    dct['gname'] = gname
    
    if neq == 100:
        dct['variables'] = 'a, i'
    elif neq == 101:
        dct['variables'] = 'a, i, b, j'
    elif neq == 86:
        dct['variables'] = 'a, i, b, j, c, k'
    

    dct['s_ind'] = summation_indices(ars)
    s, D, nocc, nactive = arstofort(ars)
    dct['D'] = D
    dct['s_inner'] = s

    name_to_open =  './'+gname+'.f90'
    
    f = open(name_to_open, 'w')
    if neq == 86:
        string = single_function_cc3_86(dct)
    else:
        string = single_function_cc3(dct)
    f.write(string)


def function_template_density_matrix_ground(mbpt, block_oo, block_ov, block_vo, block_vv, \
                                                block_oo_diag, block_vv_diag, method):

    dct = {}
    if method == 'ccsd':
        name = 'density_matrix_ground_ccsd{mbpt}'.format(mbpt=mbpt)
        dct['name'] = name
    elif method == 'cc3':
        name = 'density_matrix_ground_cc3{mbpt}'.format(mbpt=mbpt)
        dct['name'] = name

    name_to_open = './fortran_codes/transition/{name}.f90'.format(name=name)
    print(name_to_open)
    f = open(name_to_open, 'w')

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")  
    dct['date'] = datestr

    s_beginning = """module {name}

    use ccsd_transformed_integrals
    use cc3_intermediates
    use s_gen
    use basis
    
    implicit none
    !
    ! File generated automatically on {date}
    !
    contains
    """.format(**dct)
    f.write(s_beginning)

    generate_calc_D(block_oo_diag, 'oo_diag', f, 'ground_mbpt{mbpt}'.format(mbpt = mbpt), method)
    generate_calc_D(block_vv_diag, 'vv_diag', f, 'ground_mbpt{mbpt}'.format(mbpt = mbpt), method)

    generate_calc_D(block_oo, 'oo', f, 'ground_mbpt{mbpt}'.format(mbpt = mbpt), method)
    generate_calc_D(block_ov, 'ov', f, 'ground_mbpt{mbpt}'.format(mbpt = mbpt), method)
    generate_calc_D(block_vo, 'vo', f, 'ground_mbpt{mbpt}'.format(mbpt = mbpt), method)
    generate_calc_D(block_vv, 'vv', f, 'ground_mbpt{mbpt}'.format(mbpt = mbpt), method)
 
    s_end = """
    end module {name}
    """.format(**dct)
    f.write(s_end)
    f.close()


def function_template_transition_gamma_dm(block_oo, block_ov, block_vo, block_vv, method):

    dct = {}
    if method == 'ccsd':
        name = 'gamma_density_f_ccsd'
        dct['name'] = name
    elif method == 'cc3':
        name = 'gamma_density_f_cc3'
        dct['name'] = name
        
    name_to_open = './fortran_codes/transition/'+name+'.f90'

    f = open(name_to_open, 'w')

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")  
    dct['date'] = datestr

    s_beginning = """module {name}

    use ccsd_transformed_integrals
    use cc3_intermediates
    use s_gen
    use basis
    
    implicit none
    !
    ! File generated automatically on {date}
    !
    contains
    """.format(**dct)
    f.write(s_beginning)

    generate_calc_D(block_oo, 'oo', f, 'gamma', method)
    generate_calc_D(block_ov, 'ov', f, 'gamma', method)
    generate_calc_D(block_vo, 'vo', f, 'gamma', method)
    generate_calc_D(block_vv, 'vv', f, 'gamma', method)
 
    s_end = """
    end module {name}
    """.format(**dct)
    f.write(s_end)
    f.close()

def function_template_transition_gamma_dm_with_intermediates(block_oo, block_ov, block_vo, block_vv, method):

    dct = {}
    if method == 'ccsd':
        name = 'gamma_density_f_ccsd_with_intermediates'
        dct['name'] = name
    elif method == 'cc3':
        name = 'gamma_density_f_cc3_with_intermediates'
        dct['name'] = name
        
    name_to_open = './fortran_codes/transition/'+name+'.f90'

    f = open(name_to_open, 'w')

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")  
    dct['date'] = datestr

    s_beginning = """module {name}

    use ccsd_transformed_integrals
    use cc3_intermediates
    use s_gen
    use basis
    
    implicit none
    !
    ! File generated automatically on {date}
    !
    contains
    """.format(**dct)
    f.write(s_beginning)

    generate_calc_D(block_oo, 'oo', f, 'gamma_with_int', method)
    generate_calc_D(block_ov, 'ov', f, 'gamma_with_int', method)
    generate_calc_D(block_vo, 'vo', f, 'gamma_with_int', method)
    generate_calc_D(block_vv, 'vv', f, 'gamma_with_int', method)
 
    s_end = """
    end module {name}
    """.format(**dct)
    f.write(s_end)
    f.close()

def function_template_transition_xi_dm(block_oo, block_ov, block_vo, block_vv, method):

    print('lnnn', len(block_ov))
    sys.exit(0)

    dct = {}
    if method == 'ccsd':
        name = 'xi_density_f_ccsd'
        dct['name'] = name
    elif method == 'cc3':
        name = 'xi_density_f_cc3'
        dct['name'] = name
        
    name_to_open = './fortran_codes/transition/'+name+'.f90'
    f = open(name_to_open, 'w')

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")  
    dct['date'] = datestr

    s_beginning = """module {name}

    use ccsd_transformed_integrals
    use cc3_intermediates
    use s_gen
    use basis
    
    implicit none
    !
    ! File generated automatically on {date}
    !
    contains
    """.format(**dct)
    f.write(s_beginning)

    generate_calc_D(block_oo, 'oo', f, 'xi', method)
    generate_calc_D(block_ov, 'ov', f, 'xi', method)
    generate_calc_D(block_vo, 'vo', f, 'xi', method)
    generate_calc_D(block_vv, 'vv', f, 'xi', method)
 
    s_end = """
    end module {name}
    """.format(**dct)
    f.write(s_end)
    f.close()

def function_template_transition_xi_dm_with_intermediates(block_oo, block_ov, block_vo, block_vv, method):

    dct = {}
    if method == 'ccsd':
        name = 'xi_density_f_ccsd_with_intermediates'
        dct['name'] = name
    elif method == 'cc3':
        name = 'xi_density_f_cc3_with_intermediates'
        dct['name'] = name
        
    name_to_open = './fortran_codes/transition/'+name+'.f90'
    f = open(name_to_open, 'w')

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")  
    dct['date'] = datestr

    s_beginning = """module {name}

    use ccsd_transformed_integrals
    use cc3_intermediates
    use s_gen
    use basis
    
    implicit none
    !
    ! File generated automatically on {date}
    !
    contains
    """.format(**dct)
    f.write(s_beginning)

    generate_calc_D(block_oo, 'oo', f, 'xi_with_int', method)
    generate_calc_D(block_ov, 'ov', f, 'xi_with_int', method)
    generate_calc_D(block_vo, 'vo', f, 'xi_with_int', method)
    generate_calc_D(block_vv, 'vv', f, 'xi_with_int', method)
 
    s_end = """
    end module {name}
    """.format(**dct)
    f.write(s_end)
    f.close()

def generate_calc_D(block, typ, f, name, method):
    
    if method == 'ccsd':
        name = 'calc_D_{typ}_{name}'.format(typ = typ, name=name)
    elif method == 'cc3':
        name = 'calc_D_{typ}_{name}_cc3'.format(typ = typ, name=name)
    dct = {}
    dct['name'] = name

    for x in range(0, len(block)):
        block[x].optimize()
    block.cleanup()
    if typ == 'oo':
        dct['variables'] = fixed_indices(['i', 'j'], block)
    elif typ == 'oo_diag':
        dct['variables'] = fixed_indices(['i'], block)
    elif typ == 'vv_diag':
        dct['variables'] = fixed_indices(['a'], block)
    elif typ == 'ov':
        dct['variables'] = fixed_indices(['i', 'a'], block)
    elif typ == 'vo':
        dct['variables'] = fixed_indices(['a', 'i'], block)
    elif typ == 'vv':
        dct['variables'] = fixed_indices(['a', 'b'], block)
    else:
        print(typ)
        sys.exit(0)
        
#    dct['variables'] = 'p, q'
    dct['s_ind'] = summation_indices(block)
    s, D, nocc, nactive = arstofort(block)
    dct['D'] = D
    dct['s_inner'] = s
    string = single_function_transition_dm(dct)
    if len(block) != 0:
        f.write(string)

def generate_calc_D_wl(block, typ, f, name, method):
    
    if method == 'ccsd':
        name = 'calc_D_{typ}_{name}'.format(typ = typ, name=name)
    elif method == 'cc3':
        name = 'calc_D_{typ}_{name}_cc3'.format(typ = typ, name=name)
    dct = {}
    dct['name'] = name

    for x in range(0, len(block)):
        block[x].optimize()
    block.cleanup()
    dct['variables'] = fixed_indices(['p', 'q'], block)
    dct['s_ind'] = summation_indices(block)
    s, D, nocc, nactive = arstofort(block)
    dct['D'] = D
    dct['s_inner'] = s
    string = single_function_transition_wl(dct)
    if len(block) != 0:
        f.write(string)

def generate_calc_D_wm(r1_block, r2_block, typ, f, name, method):
    
    if method == 'ccsd':
        name = 'calc_D_{typ}_{name}'.format(typ = typ, name=name)
    elif method == 'cc3':
        name = 'calc_D_{typ}_{name}_cc3'.format(typ = typ, name=name)
    dct = {}
    dct['name'] = name

    for x in range(0, len(r1_block)):
        r1_block[x].optimize()
    r1_block.cleanup()

    for x in range(0, len(r2_block)):
        r2_block[x].optimize()
    r2_block.cleanup()


    sum_idx = []
    for i in range(0, len(r1_block)):
        for j in r1_block[i].summation:
            if j not in sum_idx:
                sum_idx.append(j)
    for i in range(0, len(r2_block)):
        for j in r2_block[i].summation:
            if j not in sum_idx:
                sum_idx.append(j)
    s = ['a', 'i', 'b', 'j']
    for i in s:
        if i not in sum_idx:
            sum_idx.append(i)
    print(sum_idx)

    dct['s_ind'] = ""
    for i in sum_idx:
        dct['s_ind'] += ", {i}".format(i=i)

    print('ars')

    s1, D1, nocc, nactive = arstofort(r1_block)
    s2, D2, nocc, nactive = arstofort(r2_block)
    print('s1')
    print(s1)
    print('s2')
    print(s2)

    dct['D1'] = D1
    dct['D2'] = D2
    dct['D'] = max(D1, D2)
    dct['s_inner1'] = s1
    dct['s_inner2'] = s2
    string = single_function_wm(dct)
    f.write(string)


def generate_calc_D_wm2(block, typ, f, fname, method, mbpt, dname, all_interm_dict):
    
    if method == 'ccsd':
        name = 'calc_D_{typ}_{fname}'.format(typ = typ, fname=fname, mbpt=mbpt)
    elif method == 'cc3':
        name = 'calc_D_{typ}_{fname}'.format(typ = typ, fname=fname, mbpt=mbpt)
    dct = {}
    dct['name'] = name
    dct['decl_name'] = dname
    for x in range(0, len(block)):
        block[x].optimize()
    block.cleanup()


    sum_idx = []
    for i in range(0, len(block)):
        print('bbblock', block[i])
        for j in block[i].summation:
            if j not in sum_idx:
                sum_idx.append(j)

    s = ['a', 'i', 'b', 'j']
    for i in s:
        if i not in sum_idx:
            sum_idx.append(i)
    print(sum_idx)

    dct['s_ind'] = ""
    for i in sum_idx:
        dct['s_ind'] += ", {i}".format(i=i)

    s, D, nocc, nactive = arstofort(block, 0, all_interm_dict)


    dct['D'] = D
    dct['s_inner'] = s
    string = single_function_wm2(dct)
    f.write(string)

def function_template_wm_intermediates(intermediates_dict, method, partname, multiplicity, mbpt):

    dct = {}
    if method == 'ccsd':
        name = """{partname}_intermediates_ccsd_pt{mbpt}""".format(partname=partname, mbpt=mbpt)
        if multiplicity == 3:
            name = """{partname}_triplet_intermediates_ccsd_pt{mbpt}""".format(partname=partname, mbpt=mbpt)
        if multiplicity == 13:
            name = """{partname}_so_intermediates_ccsd_pt{mbpt}""".format(partname=partname, mbpt=mbpt)
        if multiplicity == 31:
            name = """{partname}_so_left_intermediates_ccsd_pt{mbpt}""".format(partname=partname, mbpt=mbpt)
        dct['name'] = name
        print(name)
    elif method == 'cc3':
        name = """{partname}_intermediates_cc3_pt{mbpt}""".format(partname=partname, mbpt=mbpt)
        if multiplicity == 3:
            name = """{partname}_triplet_intermediates_cc3_pt{mbpt}""".format(partname=partname, mbpt=mbpt)
        if multiplicity == 13:
            name = """{partname}_so_intermediates_cc3_pt{mbpt}""".format(partname=partname, mbpt=mbpt)
        if multiplicity == 31:
            name = """{partname}_so_left_intermediates_cc3_pt{mbpt}""".format(partname=partname, mbpt=mbpt)
        dct['name'] = name
        
    name_to_open = './fortran_codes/intermediates/'+name+'.f90'
    f = open(name_to_open, 'w')

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")  
    dct['date'] = datestr


    s_decl = ""
    for i in range(0, len(intermediates_dict)):
        int_name = intermediates_dict[i]['int_name'].format(mbpt=mbpt)
        if multiplicity == 3:
            int_name = intermediates_dict[i]['int_name']
        if multiplicity == 13:
            int_name = intermediates_dict[i]['int_name']
        if multiplicity == 31:
            int_name = intermediates_dict[i]['int_name']+'_so_left'
        s_dim = ""
        for j in range(0, len(intermediates_dict[i]['idx_fx'])):
            s_dim += ":, "
        s_dim = s_dim[0:len(s_dim) - 2]
        if len(intermediates_dict[i]['idx_fx'])!= 0:
            s_decl += "real(F64), dimension({s_dim}), allocatable :: {int_name} \n".format(s_dim=s_dim, int_name=int_name)
        else:
            s_decl += "real(F64) :: {int_name} \n".format(s_dim=s_dim, int_name=int_name)
    
    dct['s_decl'] = s_decl

    s_beginning = """module {name}

    use ccsd_transformed_integrals
    use cc3_intermediates
    use arithmetic
    use s_gen
    use basis
    use eom_vectors
    
    implicit none
    !
    ! File generated automatically on {date}
    !
    {s_decl}
    contains
    """.format(**dct)
    f.write(s_beginning)
    
    s_init, s_free, s = arstofort_matrix_wm(intermediates_dict, multiplicity, mbpt)
    print(s)

    dct['s_inner'] = s
    dct['s_init'] = s_init
    dct['s_free'] = s_free
    dct['mbpt'] = mbpt
    dct['mult'] = ""
    if multiplicity == 3:
        dct['mult'] = "triplet_"
    if multiplicity == 13:
        dct['mult'] = "so_"
    if multiplicity == 31:
        dct['mult'] = "so_left_"


    dct['multiplicity'] = multiplicity
    dct['method'] = method

    s_ind_lst = []
    s_ind = ""
    for i in intermediates_dict:
        for j in i['interm'].summation:
            if j not in s_ind:
                s_ind_lst.append(j)
                s_ind += "{j}, ".format(j=j)
        for j in i['idx_fx']:
            if j not in s_ind:
                s_ind_lst.append(j)
                s_ind += "{j}, ".format(j=j)

    dct['s_ind'] = s_ind[0:len(s_ind) - 2]
    string = single_function_init(dct)
    string += single_function_free(dct)
    string += single_function_intermediates_wm(dct)

    print('str')
    print(string)

    f.write(string)
    
    s_end = """
    end module {name}
    """.format(**dct)
    f.write(s_end)
    f.close()

def function_template_wl(wl_block_oo, wl_block_ov, wl_block_vo, wl_block_vv, method):

    dct = {}
    if method == 'ccsd':
        name = 'wl_density_f_ccsd'
        dct['name'] = name
    elif method == 'cc3':
        name = 'wl_density_f_cc3'
        dct['name'] = name
        
    name_to_open = './fortran_codes/transition_excexc/'+name+'.f90'
    f = open(name_to_open, 'w')

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")  
    dct['date'] = datestr

    s_beginning = """module {name}

    use ccsd_transformed_integrals
    use cc3_intermediates
    use arithmetic
    use s_gen
    use basis
    use eom_vectors
    
    implicit none
    !
    ! File generated automatically on {date}
    !
    contains
    """.format(**dct)
    f.write(s_beginning)

    generate_calc_D_wl(wl_block_oo, 'oo', f, 'wl', method)
    generate_calc_D_wl(wl_block_ov, 'ov', f, 'wl', method)
    generate_calc_D_wl(wl_block_vo, 'vo', f, 'wl', method)
    generate_calc_D_wl(wl_block_vv, 'vv', f, 'wl', method)
 
    s_end = """
    end module {name}
    """.format(**dct)
    f.write(s_end)
    f.close()
    
def function_template_wm_overlap(Wm_oint, theory, multiplicity, mbpt):

    if multiplicity == 1:
        name_to_open = './fortran_codes/transition_excexc/{theory}_pt{mbpt}.f90'.format(theory = theory, mbpt=mbpt)
    elif multiplicity == 3:
        name_to_open = './fortran_codes/transition_excexc/{theory}_triplet_pt{mbpt}.f90'.format(theory=theory, mbpt=mbpt)
    f = open(name_to_open, 'w')

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")
    date = datestr

    s_beginning = """module {theory}_pt{mbpt}

    use ccsd_transformed_integrals
    use cc3_intermediates
    use arithmetic
    use s_gen
    use basis
    
    implicit none
    !                                                                                                                                                                                  
    ! File generated automatically on {date}                                                                                                                                          
    !  
    contains
    """.format(theory=theory, mbpt=mbpt, date=date)
    f.write(s_beginning)

   
    dct = {}
    if multiplicity == 1:
        dct['name'] = '{theory}_pt{mbpt}'.format(theory=theory, mbpt=mbpt)
    elif multiplicity == 3:
        dct['name'] = '{theory}_triplet_pt{mbpt}'.format(theory=theory, mbpt=mbpt)
    elif multiplicity == 13:
        dct['name'] = '{theory}_so_pt{mbpt}'.format(theory=theory, mbpt=mbpt)
    elif multiplicity == 31:
        dct['name'] = '{theory}_so_left_pt{mbpt}'.format(theory=theory, mbpt=mbpt)
    dct['variables'] = ''
    for x in range(0, len(Wm_oint)):
        Wm_oint[x].optimize()
    Wm_oint.cleanup()
    dct['s_ind'] = summation_indices(Wm_oint)
    s, D, nocc, nactive = arstofort(Wm_oint)
    dct['D'] = D
    dct['s_inner'] = s
    string = single_function_transition_overlap(dct)
    if len(Wm_oint) != 0:
        f.write(string)

    s_end = """
    end module {name}
    """.format(**dct)
    f.write(s_end)
    f.close()

# def function_template_wm(r1_block_oo_diag, r1_block_oo, r1_block_ov, r1_block_vo, r1_block_vv, \
#                              r2_block_oo_diag, r2_block_oo, r2_block_ov, r2_block_vo, r2_block_vv, method, wmum):
def function_template_wm(r1_block_oo, r1_block_ov, r1_block_vo, r1_block_vv, \
                         r2_block_oo, r2_block_ov, r2_block_vo, r2_block_vv, method, wmum):

    dct = {}
    if method == 'ccsd':
        if wmum == 'wm':
            name = 'wm_density_f_ccsd'
        elif wmum == 'um':
            name = 'um_density_f_ccsd'
        dct['name'] = name
    elif method == 'cc3':
        if wmum == 'wm':
            name = 'wm_density_f_cc3'
        elif wmum == 'um':
            name = 'um_density_f_cc3'
        dct['name'] = name
        
    name_to_open = './fortran_codes/transition_excexc/'+name+'.f90'
    f = open(name_to_open, 'w')

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")  
    dct['date'] = datestr

    s_beginning = """module {name}

    use ccsd_transformed_integrals
    use cc3_intermediates
    use arithmetic
    use s_gen
    use basis
    use eom_vectors
    
    implicit none
    !
    ! File generated automatically on {date}
    !
    contains
    """.format(**dct)
    f.write(s_beginning)

    # generate_calc_D_wm(r1_block_oo_diag, r2_block_oo_diag, 'oo_diag', f, 'wm', method)
    generate_calc_D_wm(r1_block_oo, r2_block_oo, 'oo', f, 'wm', method)
    generate_calc_D_wm(r1_block_ov, r2_block_ov, 'ov', f, 'wm', method)
    generate_calc_D_wm(r1_block_vo, r2_block_vo, 'vo', f, 'wm', method)
    generate_calc_D_wm(r1_block_vv, r2_block_vv, 'vv', f, 'wm', method)
 
    s_end = """
    end module {name}
    """.format(**dct)
    f.write(s_end)
    f.close()


def function_template_wm2(block_oo, block_ov, block_vo, block_vv, method, wmum, multiplicity, mbpt):

    dct = {}
    if method == 'ccsd':
        if wmum == 'wm':
            name = 'wm_density_f_ccsd_pt{mbpt}'.format(mbpt=mbpt)
            if multiplicity == 3:
                name = 'wm_density_f_ccsd_triplet_pt{mbpt}'.format(mbpt=mbpt)
            if multiplicity == 13:
                name = 'wm_density_f_ccsd_so_pt{mbpt}'.format(mbpt=mbpt)
            if multiplicity == 31:
                name = 'wm_density_f_ccsd_so_left_pt{mbpt}'.format(mbpt=mbpt)
        elif wmum == 'um':
            name = 'um_density_f_ccsd'
        dct['name'] = name
    elif method == 'cc3':
        if wmum == 'wm':
            name = 'wm_density_f_cc3_pt{mbpt}'.format(mbpt=mbpt)
            if multiplicity == 3:
                name = 'wm_density_f_cc3_triplet_pt{mbpt}'.format(mbpt=mbpt)
            if multiplicity == 13:
                name = 'wm_density_f_cc3_so_pt{mbpt}'.format(mbpt=mbpt)
            if multiplicity == 31:
                name = 'wm_density_f_cc3_so_left_pt{mbpt}'.format(mbpt=mbpt)
        elif wmum == 'um':
            name = 'um_density_f_cc3'
        dct['name'] = name
        
    name_to_open = './fortran_codes/transition_excexc/'+name+'.f90'

    f = open(name_to_open, 'w')

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")  
    dct['date'] = datestr

    s_beginning = """module {name}

    use ccsd_transformed_integrals
    use cc3_intermediates
    use arithmetic
    use s_gen
    use basis
    use eom_vectors
    
    implicit none
    !
    ! File generated automatically on {date}
    !
    contains
    """.format(**dct)
    f.write(s_beginning)

    if multiplicity == 1:
        generate_calc_D_wm2(block_oo, 'oo', f, 'wm', method, mbpt)
        generate_calc_D_wm2(block_ov, 'ov', f, 'wm', method, mbpt)
        generate_calc_D_wm2(block_vo, 'vo', f, 'wm', method, mbpt)
        generate_calc_D_wm2(block_vv, 'vv', f, 'wm', method, mbpt)
    elif multiplicity ==3:
        generate_calc_D_wm2(block_oo, 'oo', f, 'wm_triplet', method, mbpt)
        generate_calc_D_wm2(block_ov, 'ov', f, 'wm_triplet', method, mbpt)
        generate_calc_D_wm2(block_vo, 'vo', f, 'wm_triplet', method, mbpt)
        generate_calc_D_wm2(block_vv, 'vv', f, 'wm_triplet', method, mbpt)
    elif multiplicity == 13:
        print('lalala')
        generate_calc_D_wm2(block_oo, 'oo', f, 'wm_so', method, mbpt)
        generate_calc_D_wm2(block_ov, 'ov', f, 'wm_so', method, mbpt)
        generate_calc_D_wm2(block_vo, 'vo', f, 'wm_so', method, mbpt)
        generate_calc_D_wm2(block_vv, 'vv', f, 'wm_so', method, mbpt)
    elif multiplicity == 31:
        generate_calc_D_wm2(block_oo, 'oo', f, 'wm_so_left', method, mbpt)
        generate_calc_D_wm2(block_ov, 'ov', f, 'wm_so_left', method, mbpt)
        generate_calc_D_wm2(block_vo, 'vo', f, 'wm_so_left', method, mbpt)
        generate_calc_D_wm2(block_vv, 'vv', f, 'wm_so_left', method, mbpt)

 
    s_end = """
    end module {name}
    """.format(**dct)
    f.write(s_end)
    f.close()

def generate_ganu_triple(ars, f):

    ars.establish_fixed()

    gname = "ganu3"
    dct = {}
    dct['gname'] = gname

    base_sname = ''
    for i in range(0, 3):
        base_sname += virtual[i]
        base_sname += occupied[i]
    
    delta_list = deepcopy(delta_combinations(fixed))
    name_list = []

    for i in range(0, len(delta_list)):

        if delta_list[i] == []:
            sname = base_sname
            name_list.append(sname)
        else:
            sname = rename_name(base_sname, delta_list[i])
            name_list.append(sname)

    for i in range(0, len(delta_list)):
        """ Creating of dct for 
        this iteration.
        """
        delta_subst = deepcopy(delta_to_dict(delta_list[i]))
        ars_temp = deepcopy(ars)
        ars_temp.exec_delta_fixed(delta_subst)

        ars_temp_simp = simplify_fort(ars_temp)
        for x in range(0, len(ars_temp_simp)):
            ars_temp_simp[x].optimize()
        ars_temp_simp.cleanup()
        if len(ars_temp_simp) == 0:
            delta_list[i] = -1
        else:
            dct['name'] = "ganu3_{sname}".format(sname = name_list[i])
            dct['variables'] = fixed_indices(name_list[i])
            dct['s_ind'] = summation_indices(ars_temp_simp)
            s, D, nocc, nactive = arstofort(ars_temp_simp)
            s1 = ''
            dct['D'] = D
            dct['s1_inner'] = s1
            dct['s_inner'] = s
            
            string = single_function_transition_gamma(dct)
            f.write(string)

    delta_list_not_null = []    
    name_list_not_null = []
    for i in range(0, len(delta_list)):

        if delta_list[i] != -1:
            delta_list_not_null.append(delta_list[i])
            name_list_not_null.append(name_list[i])

            
    return delta_list_not_null, gname, name_list_not_null


def generate_gaxinu(ars, rank, nm, f, name0):
    
    sm_dict = {}
    sm_dict['name'] = "{name0}{rank}_{nm}".format\
        (name0 = name0, rank = rank, nm = nm)
    ars_temp = deepcopy(ars)
    if nm == 'ai':
        sm_dict['variables'] = fixed_indices('ai')
    elif nm == 'aibj':
        delta_subst = deepcopy(delta_to_dict([]))
        sm_dict['variables'] = fixed_indices('aibj')
    elif nm == 'aibi':
        delta_subst = deepcopy(delta_to_dict([['i','j']]))
        sm_dict['variables'] = fixed_indices('aib')
    elif nm == 'aiaj':
        delta_subst = deepcopy(delta_to_dict([['a', 'b']]))
        sm_dict['variables'] = fixed_indices('aij')
    elif nm == 'aiai':
        delta_subst = deepcopy(delta_to_dict([['a','b'],['i','j']]))
        sm_dict['variables'] = fixed_indices('ai')
    
    if rank == 1:
        ars_temp_simp = ars_temp
    else:
        ars_temp.exec_delta_fixed(delta_subst)
        ars_temp_simp = simplify_fort(ars_temp)
        
    ars_temp_simp.optimize()
    ars_temp_simp.cleanup()

    sm_dict['s_ind'] = summation_indices(ars_temp_simp)

    s, D, nocc, nactive = arstofort(ars_temp_simp)
    sm_dict['D'] = D
    sm_dict['s_inner'] = s

    if name0 == 'ganu':
        string = single_function_transition_gamma(sm_dict)
    elif name0 == 'ganu_ccsd':
        string = single_function_transition_gamma(sm_dict)
    elif name0 == 'xinu':
        string = single_function_transition_xi(sm_dict)
    elif name0 == 'xinu_ccsd':
        string = single_function_transition_xi(sm_dict)
    
    f.write(string)

def eom_func_triple_mem(ars_list_big, name_list_big, z_name_big, BRA, KET,theory, trans, noR):

    print(ars_list_big[2][2][10], ars_list_big[2][2][10].name)


    for z in range(0, len(ars_list_big)):
        print('start bloku', z)
        for x in range(0, len(ars_list_big[z])):
            print('polnon', x)
            print(len(ars_list_big[z][x]))
            for y in range(0, len(ars_list_big[z][x])):
                print(z, x, y, ars_list_big[z][x][y])
            print('')
    print('')

    # In all ars in ars_list fixed indices are 'a,b,c,i,j,k'.
    # In ars_list[0] fixed indices are some five of the
    # list above, so to avoid confusion, all 'abcijk', are
    # fixed.
    #
    dct = {}
    i = 0

    ars_list_big[3][0].establish_fixed()
    print(fixed)

    if noR == True:
        gname = "eom_{theory}_{dim1}{dim2}_mem_noR".format(theory = theory, dim1 = BRA, dim2 = KET, trans = trans)
    elif noR == False:
        gname = "eom_{theory}_{dim1}{dim2}_mem".format(theory = theory, dim1 = BRA, dim2 = KET, trans = trans)
    dct['gname'] = gname
    
    base_sname = ''
    for i in range(0, BRA + KET):
        base_sname += virtual[i]
        base_sname += occupied[i]


    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")

    if noR == True:
        n_to_open = "./fortran_codes/eom_functions/eom_cc3_{BRA}{KET}_mem_noR.f90".format(BRA=BRA, KET=KET)
    elif noR == False:
         n_to_open = "./fortran_codes/eom_functions/eom_cc3_{BRA}{KET}_mem.f90".format(BRA=BRA, KET=KET)
    f = open(n_to_open, 'w')
    s = "eom_cc3_{BRA}{KET}_mem".format(BRA=BRA, KET=KET)
    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")
    dct['date'] = datestr

    s_beginning = """module {gname}
    use ccsd_transformed_integrals
    use t1_transformed_int
    use basis
    use arithmetic
    use cc3_intermediates
    use cc_gparams                                                                                                                                                    
    implicit none

    contains
    """.format(**dct)
    f.write(s_beginning)

    if KET == 1:
        bound = [['n0d', 'n1d'], ['n0l', 'n1l']]
        dct['n0'] = 6
        dct['n1'] = 8
    elif KET == 2:
        bound = [['n0d', 'n1d'], ['n0l', 'n1l'], ['n0e', 'n1e'], ['n0m', 'n1m']]
        dct['n0'] = 6
        dct['n1'] = 10
    dct['bound'] = bound

    
    for y in range(0, len(ars_list_big[0])):

        dct['sname'] = 'v0_' + z_name_big[0][y][0]

        dct['zname'] = z_name_big[0][y][1]
        if noR == False:
            s, D, nocc, nactive = arstofort_mem(ars_list_big[0][y], dct, name_list_big[0][y], BRA, KET)
            dct['variables'], dct['f_ind'], t1t2_lst, nn_lst = variables(z_name_big[0][y][1], ars_list_big[0][y], nocc, nactive)
            if KET == 1:
                dct['variables'] = dct['variables'] + ', vdav, ntrial, nocc0, nvirt0, '
            elif KET == 2:
                 dct['variables'] = dct['variables'] + ',vdav, ntrial, npair, nocc0, nvirt0, '
        elif noR == True:
            s, D, nocc, nactive = arstofort_mem_noR(ars_list_big[0][y], dct, BRA, KET)
            dct['variables'], dct['f_ind'], t1t2_lst, nn_lst = variables(z_name_big[0][y][1], ars_list_big[0][y], nocc, nactive)
            if KET == 1:
                dct['variables'] = dct['variables'] + ', SIGROWS_CC3, d_small, sigma, vdav, ntrial, nocc0, nvirt0, '
            elif KET == 2:
                 dct['variables'] = dct['variables'] + ', SIGROWS_CC3, d_small, sigma, vdav, ntrial, npair, nocc0, nvirt0, nocc, '


        dct['bounds_str'] = ''
        for bd in range(0, len(dct['bound'])):
            print(dct['bound'][bd])
            for i in range(0, 2):
                print(i)
                if bd == len(dct['bound'])-1 and i == 1:
                    dct['bounds_str'] = dct['bounds_str'] + str(dct['bound'][bd][i])
                else:
                    dct['bounds_str'] = dct['bounds_str'] + str(dct['bound'][bd][i]) + ','
        dct['variables'] += dct['bounds_str']
        dct['variables'] += ', eorb, wr'
        dct['s_ind'] = summation_indices(ars_list_big[0][y])
        dct['s_ind'] += loops_indices(ars_list_big[0][y])

        s1 = ''

        dct['D'] = D
        dct['s1_inner'] = s1
        dct['s_inner'] = s
        print('zzzzzzzzzzzzzzzzzzzzz6')
        string = single_function_eom_mem(dct, t1t2_lst, nn_lst, BRA, KET, noR)
        f.write(string)


    for y in range(0, len(ars_list_big[1])):
        dct['sname'] = 'v6_' + z_name_big[1][y][0]
        dct['zname'] = z_name_big[1][y][1]
        print(dct['zname'])
        if noR == False:
            s, D, nocc, nactive = arstofort_mem(ars_list_big[1][y], dct, name_list_big[1][y], BRA, KET)
            dct['variables'], dct['f_ind'], t1t2_lst, nn_lst = variables(z_name_big[1][y][1], ars_list_big[1][y], nocc, nactive)
            if KET == 1:
                dct['variables'] = dct['variables'] + ', vdav, ntrial, nocc0, nvirt0, '
            elif KET == 2:
                 dct['variables'] = dct['variables'] + ',vdav, ntrial, npair, nocc0, nvirt0, '
        elif noR == True:
            print('pluszon', y)
            s, D, nocc, nactive = arstofort_mem_noR(ars_list_big[1][y], dct, BRA, KET)
            dct['variables'], dct['f_ind'], t1t2_lst, nn_lst = variables(z_name_big[1][y][1], ars_list_big[1][y], nocc, nactive)
            if KET == 1:
                dct['variables'] = dct['variables'] + ', SIGROWS_CC3, d_small, sigma, vdav, ntrial, nocc0, nvirt0, '
            elif KET == 2:
                dct['variables'] = dct['variables'] + ', SIGROWS_CC3, d_small, sigma, vdav, ntrial, npair, nocc0, nvirt0, nocc, '


        dct['bounds_str'] = ''
        for bd in range(0, len(dct['bound'])):
            print(dct['bound'][bd])
            for i in range(0, 2):
                print(i)
                if bd == len(dct['bound'])-1 and i == 1:
                    dct['bounds_str'] = dct['bounds_str'] + str(dct['bound'][bd][i])
                else:
                    dct['bounds_str'] = dct['bounds_str'] + str(dct['bound'][bd][i]) + ','
        dct['variables'] += dct['bounds_str']
        dct['variables'] += ', eorb, wr'
        dct['s_ind'] = summation_indices(ars_list_big[1][y])
        dct['s_ind'] += loops_indices(ars_list_big[1][y])

        s1 = ''

        dct['D'] = D
        dct['s1_inner'] = s1
        dct['s_inner'] = s
        print('zzzzzzzzzzzzzzzzzzzzz5')
        string = single_function_eom_mem(dct, t1t2_lst, nn_lst, BRA, KET, noR)
        f.write(string)


    for y in range(0, len(ars_list_big[2])):
        print('teraz bedzie obliczany blok', 2)
        dct['sname'] = 'v06_' + z_name_big[2][y][0]
        dct['zname'] = z_name_big[2][y][1]
        if noR == False:
            s, D, nocc, nactive = arstofort_mem(ars_list_big[2][y], dct, name_list_big[2][y], BRA, KET)
            dct['variables'], dct['f_ind'], t1t2_lst, nn_lst = variables(z_name_big[2][y][1], ars_list_big[2][y], nocc, nactive)            
            if KET == 1:
                dct['variables'] = dct['variables'] + ', vdav, ntrial, nocc0, nvirt0, '
            elif KET == 2:
                 dct['variables'] = dct['variables'] + ',vdav, ntrial, npair, nocc0, nvirt0, '

        elif noR == True:
            s, D, nocc, nactive = arstofort_mem_noR(ars_list_big[2][y], dct, BRA, KET)
            dct['variables'], dct['f_ind'], t1t2_lst, nn_lst = variables(z_name_big[2][y][1], ars_list_big[2][y], nocc, nactive)
            if KET == 1:
                dct['variables'] = dct['variables'] + ', SIGROWS_CC3, d_small, sigma, vdav, ntrial, nocc0, nvirt0, '
            elif KET == 2:
                dct['variables'] = dct['variables'] + ', SIGROWS_CC3, d_small, sigma, vdav, ntrial, npair, nocc0, nvirt0, nocc, '


        dct['bounds_str'] = ''
        for bd in range(0, len(dct['bound'])):
            print(dct['bound'][bd])
            for i in range(0, 2):
                print(i)
                if bd == len(dct['bound'])-1 and i == 1:
                    dct['bounds_str'] = dct['bounds_str'] + str(dct['bound'][bd][i])
                else:
                    dct['bounds_str'] = dct['bounds_str'] + str(dct['bound'][bd][i]) + ','
        dct['variables'] += dct['bounds_str']
        dct['variables'] += ', eorb, wr'
        dct['s_ind'] = summation_indices(ars_list_big[2][y])
        dct['s_ind'] += loops_indices(ars_list_big[2][y])

        s1 = ''

        dct['D'] = D
        dct['s1_inner'] = s1
        dct['s_inner'] = s
        print('zzzzzzzzzzzzzzzzzzzzz2')
        string = single_function_eom_mem(dct, t1t2_lst, nn_lst, BRA, KET, noR)
        f.write(string)


    for y in range(0, len(ars_list_big[3])):
        dct['sname'] = 'vp_' + z_name_big[3][y][0]
        dct['zname'] = z_name_big[3][y][1]
        if noR == False:
            s, D, nocc, nactive = arstofort_mem(ars_list_big[3][y], dct, name_list_big[3][y], BRA, KET)
            dct['variables'], dct['f_ind'], t1t2_lst, nn_lst = variables(z_name_big[3][y][1], ars_list_big[3][y], nocc, nactive)
            if KET == 1:
                dct['variables'] = dct['variables'] + ', vdav, ntrial, nocc0, nvirt0, '
            elif KET == 2:
                 dct['variables'] = dct['variables'] + ',vdav, ntrial, npair, nocc0, nvirt0, '
        elif noR == True:
            s, D, nocc, nactive = arstofort_mem_noR(ars_list_big[3][y], dct, BRA, KET)
            dct['variables'], dct['f_ind'], t1t2_lst, nn_lst = variables(z_name_big[3][y][1], ars_list_big[3][y], nocc, nactive)
            if KET == 1:
                dct['variables'] = dct['variables'] + ', SIGROWS_CC3, d_small, sigma, vdav, ntrial, nocc0, nvirt0, '
            elif KET == 2:
                 dct['variables'] = dct['variables'] + ', SIGROWS_CC3, d_small, sigma, vdav, ntrial, npair, nocc0, nvirt0, nocc, '

        print(z_name_big[3][y][1])

        dct['bounds_str'] = ''
        print('')
        print('')
        print(dct['bound'])
        for bd in range(0, len(dct['bound'])):
            print('b', bd, dct['bound'])
            for i in range(0, 2):
                if bd == len(dct['bound'])-1 and i == 1:
                    dct['bounds_str'] = dct['bounds_str'] + str(dct['bound'][bd][i])
                else:
                    print('else')
                    dct['bounds_str'] = dct['bounds_str'] + str(dct['bound'][bd][i]) + ','
        dct['variables'] += dct['bounds_str']
        dct['variables'] += ', eorb, wr'
        dct['s_ind'] = summation_indices(ars_list_big[3][y])
        dct['s_ind'] += loops_indices(ars_list_big[3][y])

        s1 = ''

        dct['D'] = D
        dct['s1_inner'] = s1
        dct['s_inner'] = s
        print('zzzzzzzzzzzzzzzzzzzzz1')
        string = single_function_eom_mem(dct, t1t2_lst, nn_lst, BRA, KET, noR)
        f.write(string)

    if noR == False:
        s = "eom_cc3_{BRA}{KET}_mem".format(BRA=BRA, KET=KET)
    else:
        s = "eom_cc3_{BRA}{KET}_mem_noR".format(BRA=BRA, KET=KET)
    f_end(s, f)


def eom_func_triple(ars_list, BRA, KET,theory, trans):
    
    #
    # In all ars in ars_list fixed indices are 'a,b,c,i,j,k'.
    # In ars_list[0] fixed indices are some five of the
    # list above, so to avoid confusion, all 'abcijk', are
    # fixed.
    #
    dct = {}

    ars_list[1].establish_fixed()

    gname = "eom_{theory}_{dim1}{dim2}_{trans}".format(theory = theory, dim1 = BRA, dim2 = KET, trans = trans)

    dct['gname'] = gname
    
    base_sname = ''
    for i in range(0, BRA + KET):
        base_sname += virtual[i]
        base_sname += occupied[i]

    delta_list_with_zero = deepcopy(delta_combinations(fixed))

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")

    f_list = []
    
    for i in range(0, 7):
        n_to_open = "./fortran_codes/eom_functions/v{i}_{gname}.f90".format(gname = gname, i = i)
        f_list.append(open(n_to_open, 'w'))
        s = "v{i}_{gname}".format(gname = gname, i  = i)
        f_beginning(s, f_list[i], datestr)

    
    name_list_with_zero = []

    for i in range(0, len(delta_list_with_zero)):
        """ Creating of dct for 
        this iteration.
        """
        if delta_list_with_zero[i] == []:
            sname = base_sname
            name_list_with_zero.append(sname)
        else:
            sname = rename_name(base_sname, delta_list_with_zero[i])
            name_list_with_zero.append(sname)
        if is_redundant(BRA, KET, name_list_with_zero[i]):
            delta_list_with_zero[i] = -1

    name_list = []
    delta_list = []

    for i in range(0, len(delta_list_with_zero)):
        if delta_list_with_zero[i] != -1:
            delta_list.append(delta_list_with_zero[i])
            name_list.append(name_list_with_zero[i])
        

    vanish_case1 = []
    arg_list = []
    print(len(delta_list))
    for i in range(0, len(delta_list)):
        """ Creating of dct for 
        this iteration.
        """
        delta_subst = deepcopy(delta_to_dict(delta_list[i]))

        vk = whichvk((BRA,KET), delta_list[i])

        if vk == -1:
            delta_list[i] = -1
            arg_list.append(-1)
        elif vk == 0:
            dl = create_module(0, ars_list[0], dct, delta_subst, name_list[i], delta_list[i], gname, f_list[0], 'excitation', \
                                   arg_list)
            if dl == -1:
                delta_list[i] = -1
                arg_list.append(-1)
        if vk == 6:
            dl = create_module(6, ars_list[6], dct, delta_subst, name_list[i], delta_list[i], gname, f_list[6], 'excitation', \
                                   arg_list)
            if dl == -1:
                delta_list[i] = -1
                arg_list.append(-1)
        elif vk in range(1,6):
            counter = 0
            arg_list.append([])
            ln = len(arg_list)

            for j in range(1, 6):
                dl = create_module(j, ars_list[j], dct, delta_subst, name_list[i], delta_list[i], gname, f_list[j], 'excitation', \
                                       arg_list[ln-1])

                print('dldldlldlldldlllllllll', dl)

                if dl == -1:
                    counter += 1
                    s = "v{j}_{gname}_{sname}".format(j = j, gname = gname, sname = name_list[i])
                    vanish_case1.append(s)
                    arg_list[ln-1].append(-1)
            if counter == 5:
                delta_list[i] = -1
                arg_list[ln-1] = -1
            print(arg_list[ln-1])


    print('lalalala')

    for i in range(0, 7):
        s = "v{i}_{gname}".format(gname = gname, i = i)
        f_end(s, f_list[i])

    delta_list_not_null = []    
    name_list_not_null  = []

    arg_list_not_null = []
    print(len(delta_list))
    print(len(arg_list))

    for i in range(0, len(delta_list)):

        if delta_list[i] != -1:
            delta_list_not_null.append(delta_list[i])
            name_list_not_null.append(name_list[i])
            arg_list_not_null.append(arg_list[i])

#    print(len(arg_list), len(delta_list_not_null))
#    for i in range (0, len(delta_list_not_null)):
#        print(name_list_not_null[i], arg_list_not_null[i])
            

    return delta_list_not_null, gname, name_list_not_null, vanish_case1, arg_list_not_null




def f_beginning(name, f, datestr):

    s_beginning = """module {gname}

    use ccsd_transformed_integrals                                                                                                                                   
    use t1_transformed_int                                                                                                                   
    use eom_vectors                                                                                                                               
    use basis                                                                                                                           
    use arithmetic                                                                                                                           
    use s_gen                                                                    
    use cc3_intermediates 
    use ccsd_transformed_integrals
    use t1_transformed_int
    use basis
    
    implicit none
    !
    ! File generated automatically on {datestr} UTC.
    !
    contains
    """.format(gname = name, datestr = datestr)
    f.write(s_beginning)


def f_beginning_trans(name, f, datestr):

    s_beginning = """module {name}

    use ccsd_transformed_integrals
    use cc3_intermediates
    use s_gen
    use basis
    
    implicit none
    !
    ! File generated automatcally on {datestr}
    !
    contains
    """.format(name = name, datestr = datestr)
    f.write(s_beginning)

def f_end(name, f):

    s_end = """
    end module {name}
    """.format(name = name)
    f.write(s_end)
    f.close()    

def create_module(j, ars, dct, delta_subst, name, delta_list, gname, f, propertie, arg_list):

    print('module', j, 'fx', fixed, name)
    dl = 1

    s = "v{j}_{gname}".format(j = j, gname = gname)
    print(s, 'przed')
    for x in ars:
        print(x)

    dct['gname'] = s
    ars_temp = deepcopy(ars)
    
    print('a')
    ars_temp.exec_delta_fixed(delta_subst)
    print('b', 'po')
    for x in ars_temp:
        print(x)
    ars_temp_simp = simplify_fort(ars_temp)
    print('mod po simp', name)
    for x in ars_temp_simp:
        print(x)
    print('')
    print('po optimize')
    for x in range(0, len(ars_temp_simp)):
        ars_temp_simp[x].optimize()
        print(ars_temp_simp[x])
    ars_temp_simp.cleanup()
    if len(ars_temp_simp) == 0:
        dl = -1
    else:
        dct['sname'] = name
        dct['s_ind'] = summation_indices(ars_temp_simp)
        s, D, nocc, nactive = arstofort(ars_temp_simp)
        print('adfgadgag')
        dct['variables'], dct['f_ind'], t1t2_lst, nn_lst = variables(name, ars_temp_simp, nocc, nactive)
        arg_list.append(dct['variables'])
        dct['D'] = D
        dct['s_inner'] = s
        
        if propertie == 'excitation':
            string = single_function_eom(dct, t1t2_lst, nn_lst)
        elif propertie == 'transition':
            dct['name'] = '{gname}_{sname}'.format(gname=dct['gname'],sname=dct['sname'])
            string = single_function_transition_xi_triples(dct)
        f.write(string)

    return dl

def eom_func_cisd(ars, BRA, KET, triplet = ''):
    
    ars.establish_fixed()

    gname = "cisd{dim1}{dim2}".format(dim1 = BRA, dim2 = KET)
    dct = {}
    dct['gname'] = gname

    #
    # Find base name eg.
    # BRA = 1, KET = 1, base_sname = aibj
    # BRA = 1, KET = 2, base_sname = aibjck
    #
    base_sname = ''
    for i in range(0, BRA + KET):
        base_sname += virtual[i]
        base_sname += occupied[i]
        
    #
    # Delta_list_with_zero is the list of all posible deltas from
    # fixed indices, including deltas that gives zero when applied
    # to arithmetic string
    #

    delta_list_with_zero = deepcopy(delta_combinations(fixed))

    n_to_open = "./fortran_codes/cisd/{gname}.f90".format(gname = gname)
    f = open(n_to_open, 'w')
    s_beginning = """module {gname}

    use ccsd_transformed_integrals
    use basis
    
    implicit none

    contains
    """.format(**dct)
    f.write(s_beginning)

    #
    # Create name list, with specified name - sname
    # changed according to current delta combination.
    # 
    # Check if for the current delta combinations
    # ars = 0, if positive, change delta_list to -1
    #
    name_list_with_zero = []
    k = 0
    l = 0
    for i in range(0, len(delta_list_with_zero)):
        """ Creation of dct for 
        this iteration.
        """
        print(delta_list_with_zero[i])

        if delta_list_with_zero[i] == []:
            sname = base_sname
            name_list_with_zero.append(sname)
        else:
            sname = rename_name(base_sname, delta_list_with_zero[i])
            name_list_with_zero.append(sname)

        if triplet != "":
           if is_redundant_triplet(BRA, KET, name_list_with_zero[i], lmfold, rmfold):
               delta_list_with_zero[i] = -1
        else:
            if is_redundant(BRA, KET, name_list_with_zero[i]):
                delta_list_with_zero[i] = -1
    
        if delta_list_with_zero[i] != -1:
            print('zostawiam',  i, name_list_with_zero[i])
            k += 1
        else:
            print('odrzucam',  i, name_list_with_zero[i])
            l += 1
    print('k', k, l, len(delta_list_with_zero))
    
    
    

    #
    # name_list and delta_list without "zero" elements
    #
    name_list  = []
    delta_list = []
    for i in range(0, len(delta_list_with_zero)):
        print(i)
        if delta_list_with_zero[i] != -1:
            delta_list.append(delta_list_with_zero[i])
            name_list.append(name_list_with_zero[i])
            print(name_list_with_zero[i])

        
    arg_list = []

    print('')
    for x in ars:
        print(x)
    print('')

    for i in range(0, len(delta_list)):
        """ Creating of dct for 
        this iteration.
        """
        print(name_list[i])
        delta_subst = deepcopy(delta_to_dict(delta_list[i]))

        ars_temp = deepcopy(ars)
        
        ars_temp.exec_delta_fixed(delta_subst)

        ars_simp = simplify_fort(ars_temp)


        for x in range(0, len(ars_simp)):
            ars_simp[x].optimize()
        ars_simp.cleanup()


        if len(ars_simp) == 0:
            delta_list[i] = -1
            arg_list.append(0)
        else:
            dct['sname'] = name_list[i]            
            print('iiiiiii', len(ars_simp))

            s, D, nocc, nactive = arstofort(ars_simp)
            print(s)

            dct['variables'], dct['f_ind'], t1t2_lst, nn_lst = variables(name_list[i], ars_simp, nocc, nactive)
            arg_list.append(dct['variables'])
            dct['s_ind'] = summation_indices(ars_simp)

            s1 = ''

            dct['D'] = D
            dct['s1_inner'] = s1
            dct['s_inner'] = s

            string = single_function_eom(dct, t1t2_lst, nn_lst)
            f.write(string)

    s_end = """
    end module {gname}
    """.format(**dct)

    f.write(s_end)
    f.close()

    delta_list_nn = []    
    name_list_nn = []
    arg_list_nn = []
    print('')
    for i in range(0, len(delta_list)):
        
        if delta_list[i] != -1:
            delta_list_nn.append(delta_list[i])
            name_list_nn.append(name_list[i])
            arg_list_nn.append(arg_list[i])

    return delta_list_nn, gname, name_list_nn, arg_list_nn

def eom_func_mem(ars_list, name_list, z_name, BRA, KET, theory, trans, triplet = "", lmfold = "", rmfold = "", noR = False):
    #
    # establish general name of the function
    #
    if noR == False:
        print('atttttu')
        gname = "eom_cc3_{dim1}{dim2}_mem{triplet}{lmfold}{rmfold}"\
                .format(theory = theory, dim1 = BRA, dim2 = KET, triplet=triplet, lmfold = lmfold, rmfold = rmfold)
        print(gname)
    else:
        gname = "eom_cc3_{dim1}{dim2}_mem_noR{triplet}{lmfold}{rmfold}"\
                .format(theory = theory, dim1 = BRA, dim2 = KET, triplet=triplet, lmfold = lmfold, rmfold = rmfold)

    dct = {}
    dct['gname'] = gname
            
    if triplet == "":
        print('dfasdfasdf', noR)
        if noR == False:
            n_to_open = "./fortran_codes/eom_functions/eom_cc3_{BRA}{KET}_mem.f90".format(BRA=BRA, KET=KET)
        else:
            n_to_open = "./fortran_codes/eom_functions/eom_cc3_{BRA}{KET}_mem_noR.f90".format(BRA=BRA, KET=KET)
    else:
        if noR == False:
            n_to_open = "./fortran_codes/eom_functions/triplet/eom_cc3_{BRA}{KET}_mem.f90".format(BRA=BRA, KET=KET)
        else:
            n_to_open = "./fortran_codes/eom_functions/triplet/eom_cc3_{BRA}{KET}_mem_noR.f90".format(BRA=BRA, KET=KET)
    f = open(n_to_open, 'w')
    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")
    dct['date'] = datestr
    s_beginning = """module {gname}

    use ccsd_transformed_integrals
    use t1_transformed_int
    use basis
    use arithmetic
    use cc3_intermediates
    use cc_gparams
    
    implicit none
    !               
    ! File generated automatically on {date}  
    !  
    contains
    """.format(**dct)
    f.write(s_beginning)

    if BRA == 1:
        bound = [['n0a', 'n1a'], ['n0i', 'n1i']]
        dct['n0'] = 0
        dct['n1'] = 2
    elif BRA == 2:
        bound = [['n0a', 'n1a'], ['n0i', 'n1i'], ['n0b', 'n1b'], ['n0j', 'n1j']]
        dct['n0'] = 0
        dct['n1'] = 4
    dct['bound'] = bound

    print('')
    arg_list = []

    for y in range(0, len(ars_list)):
        for x in range(0, len(ars_list[y])):
            ars_list[y][x].optimize()

    for y in range(0, len(ars_list)):
        print('arslistARSLIST', y, z_name[y][0])
        #print(ars_list[y], z_name[y][0])
        dct['sname'] = z_name[y][0]
        dct['zname'] = z_name[y][1]
        if noR == False:
            s, D, nocc, nactive = arstofort_mem(ars_list[y], dct, name_list[y], BRA, KET)
        else:
            s, D, nocc, nactive = arstofort_mem_noR(ars_list[y], dct, BRA, KET)

        # dct['s_ind'] = summation_indices(ars_list[y])
        # for x in dct['s_ind']:
            
        dct['variables'], dct['f_ind'], t1t2_lst, nn_lst = variables(z_name[y][1], ars_list[y], nocc, nactive)

        if noR == False:
            if KET == 1 or BRA == 1:
                dct['variables'] = dct['variables'] + ', vdav, ntrial, nocc0, nvirt0, '
            elif KET == 2 or BRA == 2:
                 dct['variables'] = dct['variables'] + ',vdav, ntrial, npair, nocc0, nvirt0, '
        else:
            if BRA == 2:
                dct['variables'] = dct['variables'] + ', SIGROWS_CC3, d_small, sigma, vdav, ntrial, npair, nocc0, nvirt0, nocc, '
            elif BRA == 1:
                dct['variables'] = dct['variables'] + ', SIGROWS_CC3, d_small, sigma, vdav, ntrial, nocc0, nvirt0, nocc, '
        dct['bounds_str'] = ''
        for bd in range(0, len(dct['bound'])):
            print(dct['bound'][bd]) 
            for i in range(0, 2):
                print(i)
                if bd == len(dct['bound'])-1 and i == 1:
                    dct['bounds_str'] = dct['bounds_str'] + str(dct['bound'][bd][i]) 
                else:
                    dct['bounds_str'] = dct['bounds_str'] + str(dct['bound'][bd][i]) + ','
        dct['variables'] += dct['bounds_str']
        if noR :
            dct['variables'] += ', eorb, wr'
        print(dct['variables'])

        arg_list.append(dct['variables'])
        dct['s_ind'] = summation_indices(ars_list[y])

        s1 = ''

        dct['D'] = D
        dct['s1_inner'] = s1
        dct['s_inner'] = s
        if noR == False:
            string = single_function_eom_mem(dct, t1t2_lst, nn_lst, BRA, KET)
        else:
            string = single_function_eom_mem(dct, t1t2_lst, nn_lst, BRA, KET, noR)
        f.write(string)

    s_end = """
    end module {gname}
    """.format(**dct)

    f.write(s_end)
    f.close()

def eom_func_easy(ars, BRA, KET, theory, trans, triplet = "", lmfold = "", rmfold = ""):
    
    #
    # In all ars in ars_list fixed indices are 'a,b,c,i,j,k'.
    #
    dct = {}

    ars.establish_fixed()
    print('ars')
    for x in range(0, len(ars)):
        print(ars[x])
        ars[x].optimize()
    print('')


    gname = "{triplet}_{dim1}{dim2}{lmfold}{rmfold}_{trans}"\
        .format(dim1 = BRA, dim2 = KET, trans = trans, \
                    triplet=triplet, lmfold = lmfold, rmfold = rmfold)

    dct['gname'] = gname
    
    base_sname = ''
    for i in range(0, BRA + KET):
        base_sname += virtual[i]
        base_sname += occupied[i]


    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")

    
    n_to_open = "./fortran_codes/eom_functions/triplet/{gname}.f90".format(gname = gname)

    f = open(n_to_open, 'w')

    delta_list = []
    for x in ars:
        if x.delta not in delta_list:
            delta_list.append(x.delta)

    ars0_list = []
    for x in delta_list:
        ars_mini = arithmetic_string()
        for y in ars:
            if y.delta == x:
                y.delta = []
                if len(ars_mini) == 0:
                    ars_mini = arithmetic_string(y)
                else:
                    ars_mini += arithmetic_string(y)
        print('dla delt', x)
        for z in ars_mini:
            print(z)
        print('')
        ars0_list.append(ars_mini)

    s = "term = 0.d+0\n"
    D = 0
    print('ppp')
    for x in range(0, len(delta_list)):
        print('x', x)
        sb, se = ifdelta(delta_list[x])
        s += sb
        print('D', D)
        ss, D, nocc, nactive = arstofort(ars0_list[x], D)
        print(ss)
        if (D==0):
            D = 1
        s+= ss
        s+= se


    f.write(s)



def eom_func(ars, BRA, KET, theory, trans, triplet = "", lmfold = "", rmfold = ""):

    for x in ars:
        print(x)
    print('')

    ars.establish_fixed()
    #
    # Establish general name of the function
    #
    gname = "eom_{theory}_{dim1}{dim2}{triplet}{lmfold}{rmfold}_{trans}"\
        .format(theory = theory, dim1 = BRA, dim2 = KET, trans = trans, triplet=triplet, lmfold = lmfold, rmfold = rmfold)
    dct = {}
    dct['gname'] = gname

    #
    # Find base name eg.
    # BRA = 1, KET = 1, base_sname = aibj
    # BRA = 1, KET = 2, base_sname = aibjck
    #
    base_sname = ''
    for i in range(0, BRA + KET):
        base_sname += virtual[i]
        base_sname += occupied[i]
    
        
    #
    # Delta_list_with_zero is the list of all posible deltas from
    # fixed indices, including deltas that gives zero when applied
    # to arithmetic string
    #
    delta_list_with_zero = deepcopy(delta_combinations(fixed))

    if triplet == "":
        n_to_open = "./fortran_codes/eom_functions/{gname}.f90".format(gname = gname)
    else:
        n_to_open = "./fortran_codes/eom_functions/triplet/{gname}.f90".format(gname = gname)
    f = open(n_to_open, 'w')

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")
    dct['date'] = datestr
    s_beginning = """module {gname}

    use ccsd_transformed_integrals
    use t1_transformed_int
    use cc3_intermediates_for_21
    use basis
    
    implicit none
    !               
    ! File generated automatically on {date}  
    !  
    contains
    """.format(**dct)
    f.write(s_beginning)

    #
    # Create name list, with specified name - sname
    # changed according to current delta combination.
    # 
    # Check if for the current delta combinations
    # ars = 0, if positive, change delta_list to -1
    #
    name_list_with_zero = []
    k = 0
    l = 0
    for i in range(0, len(delta_list_with_zero)):
        """ Creation of dct for 
        this iteration.
        """
#        print(delta_list_with_zero[i])

        if delta_list_with_zero[i] == []:
            sname = base_sname
            name_list_with_zero.append(sname)
        else:
            sname = rename_name(base_sname, delta_list_with_zero[i])
            name_list_with_zero.append(sname)

        if triplet != "":
           if is_redundant_triplet(BRA, KET, name_list_with_zero[i], lmfold, rmfold):
               delta_list_with_zero[i] = -1
        else:
            if is_redundant(BRA, KET, name_list_with_zero[i]):
                delta_list_with_zero[i] = -1
    
        if delta_list_with_zero[i] != -1:
            print('zostawiam',  i, name_list_with_zero[i])
            k += 1
        else:
            print('odrzucam',  i, name_list_with_zero[i])
            l += 1

    print('k', k, l, len(delta_list_with_zero))

    #
    # name_list and delta_list without "zero" elements
    #
    name_list  = []
    delta_list = []
    for i in range(0, len(delta_list_with_zero)):
        print(i)
        if delta_list_with_zero[i] != -1:
            delta_list.append(delta_list_with_zero[i])
            name_list.append(name_list_with_zero[i])
            print(name_list_with_zero[i])

    print('')
    arg_list = []

    print('po utworzeniu delta_list_without_zero')
    kk = 1
    for x in ars:
        print(kk, x)
        kk += 1
    print('')

    
    for i in range(0, len(delta_list)):
        """ Creating dct for 
        this iteration.
        """

        print('i', i)

        delta_subst = deepcopy(delta_to_dict(delta_list[i]))

        ars_temp = deepcopy(ars)
        ars_temp.exec_delta_fixed(delta_subst)

        ars_simp = simplify_fort(ars_temp)


        for x in range(0, len(ars_simp)):
            ars_simp[x].optimize()
        ars_simp.cleanup()
        
        # kk = 1
        # for x in ars_simp:
        #     pr = False
        #     for y in name_list[i]:
        #         idx_sm = []
        #         for x1 in x.coefficient_idx:
        #             idx_sm  = idx_sm + x1
        #         if y not in idx_sm:
        #             pr = True
        #     # if (len(x.summation) >= 3):
        #     #     vvv = 0
        #     #     for y in x.summation:
        #     #         if y in virtual:
        #     #             vvv += 1
        #     #     if (vvv >= 2):
        #     if (pr == True):
        #         print(i, kk, x, name_list[i])
        #     kk += 1
        # print('')


        if len(ars_simp) == 0:
            print('len=0')
            delta_list[i] = -1
            arg_list.append(0)
        else:
            print('name', name_list[i])
            if name_list[i] =='aibjcicj':
                print('pluszlallal')

            dct['sname'] = name_list[i]

            print('przed arstofort')
            for x in ars_simp:
                print('z', x)
                
                
            s, D, nocc, nactive = arstofort(ars_simp)
            print(s)

            if name_list[i] =='aibjcicj':
                sys.exit(0)
            dct['variables'], dct['f_ind'], t1t2_lst, nn_lst = variables(name_list[i], ars_simp, nocc, nactive)
            arg_list.append(dct['variables'])
            dct['s_ind'] = summation_indices(ars_simp)

            if triplet == "":
                if theory == 'cc3' and BRA == 2 and KET == 1:
                    s1 = create_omega(name_list[i])
                else:
                    s1 = ''
            else:
                s1 = ''

            dct['D'] = D
            dct['s1_inner'] = s1
            dct['s_inner'] = s

            if theory == 'cc3' and BRA == 2 and KET == 1:
                string = single_function_eom_cc3(dct, t1t2_lst, nn_lst)
            else:
                string = single_function_eom(dct, t1t2_lst, nn_lst)
            f.write(string)



    s_end = """
    end module {gname}
    """.format(**dct)

    f.write(s_end)
    f.close()

    delta_list_nn = []    
    name_list_nn = []
    arg_list_nn = []
    print('')
    for i in range(0, len(delta_list)):
        
        if delta_list[i] != -1:
            delta_list_nn.append(delta_list[i])
            name_list_nn.append(name_list[i])
            arg_list_nn.append(arg_list[i])


    return delta_list_nn, gname, name_list_nn, arg_list_nn
    
def eom_func_triple_triplet(ars_list, BRA, KET, theory, trans, triplet = "", lmfold = "", rmfold = ""):
    
    #
    # In all ars in ars_list fixed indices are 'a,b,c,i,j,k'.
    #
    dct = {}

    for x in range(0, len(ars_list)):
        for y in range(0, len(ars_list[x])):
            ars_list[x][y].constraints = {}

    ars_list[0].establish_fixed()
    # print('start fx')
    # for x in ars_list:
    #     print('fixed', fixed)
    #     x.establish_fixed()
    #     print('fixed', fixed)
    #     x.clear_fixed()
    #     print('')
    # sys.exit(0)


    gname = "eom_{theory}_{dim1}{dim2}{triplet}{lmfold}{rmfold}_{trans}"\
        .format(theory = theory, dim1 = BRA, dim2 = KET, trans = trans, \
                    triplet=triplet, lmfold = lmfold, rmfold = rmfold)

    dct['gname'] = gname
    
    base_sname = ''
    for i in range(0, BRA + KET):
        base_sname += virtual[i]
        base_sname += occupied[i]

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")
    f_list = []
    n_to_open = "./fortran_codes/eom_functions/triplet/v9_{gname}.f90".format(gname = gname)
    f_list.append(open(n_to_open, 'w'))
    s = "v9_{gname}".format(gname = gname)
    f_beginning(s, f_list[0], datestr)

    for i in range(1, 9):
        n_to_open = "./fortran_codes/eom_functions/triplet/v{i}_{gname}.f90".format(gname = gname, i = i)
        f_list.append(open(n_to_open, 'w'))
        s = "v{i}_{gname}".format(gname = gname, i  = i)
        f_beginning(s, f_list[i], datestr)



    # delta_list_with_zero_list = []
    # name_list_with_zero_list = []
    # for i in range(0, len(ars_list)):
    #     ars_list[i].establish_fixed()
    #     delta_list_with_zero_list.append(deepcopy(delta_combinations(fixed)))
    #     ars_list[i].clear_fixed()
    #     name_list_with_zero_list.append([])

    # for k in range(0, len(delta_list_with_zero_list)):
    #     for i in range(0, len(delta_list_with_zero_list[k])):
    #         if delta_list_with_zero_list[k][i] == []:
    #             sname = base_sname
    #             name_list_with_zero_list[k].append(sname)
    #         else:
    #             sname = rename_name(base_sname, delta_list_with_zero_list[k][i])
    #             name_list_with_zero_list[k].append(sname)
    #         if triplet != "":
    #             if is_reundant_triplet(BRA, KET, name_list_with_zero[k][i], lmfold, rmfold):
    #                 delta_list_with_zero_list[k][i] = -1
    #         else:
    #             if is_redundant(BRA, KET, name_list_with_zero_list[k][i]):
    #                 delta_list_with_zero_list[k][i] = -1

    # name_big_list = []
    # delta_big_list = []
    
    # for k in range(0, len(delta_list_with_zero_list)):
    #     mini_delta = []
    #     mini_name = []
    #     for i in range(0, len(delta_list_with_zero[k])):
    #         if delta_list_with_zero[k][i] != -1:
    #             mini_delta.append(delta_list_with_zero[k][i])
    #             mini_name.append(name_list_with_zero[k][i])
    #     name_big_list.apend(mini_name)
    #     delta_big_list.append(mini_delta)



    print('fx', fixed)
    delta_list_with_zero = deepcopy(delta_combinations(fixed))
    name_list_with_zero = []

    for i in range(0, len(delta_list_with_zero)):
        """ Creating of dct for 
        this iteration.
        """
        if delta_list_with_zero[i] == []:
            sname = base_sname
            name_list_with_zero.append(sname)
        else:
            sname = rename_name(base_sname, delta_list_with_zero[i])
            name_list_with_zero.append(sname)
        if triplet != "":
            if is_redundant_triplet(BRA, KET, name_list_with_zero[i], lmfold, rmfold):
#                print('sral', name_list_with_zero[i])
                delta_list_with_zero[i] = -1
                
        else:
            if is_redundant(BRA, KET, name_list_with_zero[i]):
                delta_list_with_zero[i] = -1

    name_list = []
    delta_list = []

    for i in range(0, len(delta_list_with_zero)):
        if delta_list_with_zero[i] != -1:
            delta_list.append(delta_list_with_zero[i])
            name_list.append(name_list_with_zero[i])
            
    # big_arg_list = []

    # for k in range(0, len(ars_list[0])):
    #     arg_list = []
    #     for i in range(0, len(delta_big_list[k])):
    #         delta_subst = deepcopy(delta_to_dict(delta_big_list[k][i]))
    #         vk = whichvk_triplet2((BRA, KET), delta_big_list[k][i])

            
    #         if vk == -1:
    #             print('vk = -1')
    #             sys.exit(1)
    #         elif vk == 0:
    #             dl = create_module(9, ars_list[0], dct, delta_subst, name_list[i], delta_list[i], gname, f_list[0], 'excitation', \
    #                                arg_list)
    #                               arg_list)
    #                 if dl == -1:
    #                     delta_big_list[k][i] = -1
    #                     arg_list.append(-1)

    #         elif vk == k and vk != 0:
    #             dl = create_module(k, ars_list[k], dct, delta_subst, name_list[i], delta_list[i], gname, f_list[0], 'excitation', \
    #                                arg_list)
    #             if dl == -1:
    #                 delta_big_list[k][i] = -1
    #                 arg_list.append(-1)
    #     big_arg_list.append(arg_list)



    vanish_case1 = []
    arg_list = []
    for i in range(0, len(delta_list)):
        """ Creating dct for 
        this iteration.
        """
        delta_subst = deepcopy(delta_to_dict(delta_list[i]))
        
        print('delta_subst przed whichvk', name_list[i])
        
        vk = whichvk_triplet2((BRA, KET), delta_list[i])
        print('vkkkkkkk', vk)


        if vk == -1:
            delta_list[i] = -1
            arg_list.append(-1)
        elif vk == 9:
            dl = create_module(9, ars_list[0], dct, delta_subst, name_list[i], delta_list[i], gname, f_list[0], 'excitation', \
                                   arg_list)
            if dl == -1:
                delta_list[i] = -1
                arg_list.append(-1)

        elif vk == 1:
            dl = create_module(1, ars_list[1], dct, delta_subst, name_list[i], delta_list[i], gname, f_list[1], 'excitation', \
                                   arg_list)
            if dl == -1:
                delta_list[i] = -1
                arg_list.append(-1)
        elif vk == 2:
            dl = create_module(2, ars_list[2], dct, delta_subst, name_list[i], delta_list[i], gname, f_list[2], 'excitation', \
                                   arg_list)
            if dl == -1:
                delta_list[i] = -1
                arg_list.append(-1)
        elif vk == 3:
            dl = create_module(3, ars_list[3], dct, delta_subst, name_list[i], delta_list[i], gname, f_list[3], 'excitation', \
                                   arg_list)
            if dl == -1:
                delta_list[i] = -1
                arg_list.append(-1)
        elif vk == 4:
            dl = create_module(4, ars_list[4], dct, delta_subst, name_list[i], delta_list[i], gname, f_list[4], 'excitation', \
                                   arg_list)
            if dl == -1:
                delta_list[i] = -1
                arg_list.append(-1)
        elif vk == 5:
            dl = create_module(5, ars_list[5], dct, delta_subst, name_list[i], delta_list[i], gname, f_list[5], 'excitation', \
                                   arg_list)
            if dl == -1:
                delta_list[i] = -1
                arg_list.append(-1)
        elif vk == 6:
            dl = create_module(6, ars_list[6], dct, delta_subst, name_list[i], delta_list[i], gname, f_list[6], 'excitation', \
                                   arg_list)
            if dl == -1:
                delta_list[i] = -1
                arg_list.append(-1)
        elif vk == 7:
            dl = create_module(7, ars_list[7], dct, delta_subst, name_list[i], delta_list[i], gname, f_list[7], 'excitation', \
                                   arg_list)
            if dl == -1:
                delta_list[i] = -1
                arg_list.append(-1)

        elif vk == 8:
            print('vk8', name_list[i])
            dl = create_module(8, ars_list[8], dct, delta_subst, name_list[i], delta_list[i], gname, f_list[8], 'excitation', \
                                   arg_list)
            if dl == -1:
                delta_list[i] = -1
                arg_list.append(-1)
                

    s = "v9_{gname}".format(gname = gname)
    f_end(s, f_list[0])
    for i in range(1, 8):#4):
        s = "v{i}_{gname}".format(gname = gname, i = i)
        f_end(s, f_list[i])

    delta_list_not_null = []    
    name_list_not_null  = []

    arg_list_not_null = []

    print('delt')
    print(delta_list)

    for i in range(0, len(delta_list)):

        if delta_list[i] != -1:
            delta_list_not_null.append(delta_list[i])
            name_list_not_null.append(name_list[i])
            print(name_list[i])

            arg_list_not_null.append(arg_list[i])
            print('a', arg_list[i])
            
    

    return delta_list_not_null, gname, name_list_not_null, vanish_case1, arg_list_not_null


def eom_func_triple_triplet_easy(ars_list, BRA, KET, theory, trans, triplet = "", lmfold = "", rmfold = ""):
    
    #
    # In all ars in ars_list fixed indices are 'a,b,c,i,j,k'.
    #
    dct = {}

    k = 0
    for x in range(0, len(ars_list)):
        print(' wyrazlalalall', k, fixed)
        k += 1
        for y in range(0, len(ars_list[x])):
 #           print(ars_list[x][y])
            ars_list[x][y].optimize()
#            print(ars_list[x][y])
        print('')
    sys.exit(0)


    ars_list[0].establish_fixed()
    k = 0
    for x in ars_list:
        print(k)
        k += 1
        for y in x:
            print(y)



    gname = "{triplet}_{dim1}{dim2}{lmfold}{rmfold}_{trans}"\
        .format(dim1 = BRA, dim2 = KET, trans = trans, \
                    triplet=triplet, lmfold = lmfold, rmfold = rmfold)

    dct['gname'] = gname
    
    base_sname = ''
    for i in range(0, BRA + KET):
        base_sname += virtual[i]
        base_sname += occupied[i]


    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")

    
    n_to_open = "./fortran_codes/eom_functions/triplet/{gname}.f90".format(gname = gname)

    f = open(n_to_open, 'w')

    bigs = ""
    for z in range(0, len(ars_list)):
        print('zzzzzzzzzzzz', z)
        s = 'ars_list{z}\n'.format(z=z) 
        delta_list = []
        for x in ars_list[z]:
            if x.delta not in delta_list:
                delta_list.append(x.delta)

        print('niebieski')
        ars0_list = []
        for x in delta_list:
            ars_mini = arithmetic_string()
            for y in ars_list[z]:
                if y.delta == x:
                    y.delta = []
                    if len(ars_mini) == 0:
                        ars_mini = arithmetic_string(y)
                    else:
                        ars_mini += arithmetic_string(y)
            for y in ars_mini:
                print(y)
            print('')
            ars0_list.append(ars_mini)
        print('')

    
        s += "term = 0.d+0\n"
        D = 0
        print('ppp')
        for x in range(0, len(delta_list)):
            
            sb, se = ifdelta(delta_list[x])
            s += sb
            print('D', D)
            ss, D, nocc, nactive = arstofort(ars0_list[x], D)
            print('ss', ss)
            if (D==0):
                D = 1
            s+= ss
            s+= se


        bigs += s

    f.write(bigs)


def ifdelta(delta_list):

    if len(delta_list) == 0:
        sb =""
        se = ""
        return sb, se

    else:
        se = "end if \n"
        ln = len(delta_list)
        sb = "if ("
        for x in range(0, ln):
            sb += "{a} == {b}".format(a=delta_list[x][0], b=delta_list[x][1])
            if (x+1==ln):
                sb+= ") then \n"
            else:
                sb += ".and."

    return sb, se
