from params import *
import paldus_classes
from paldus_basic import *
from paldus_classes import ugg
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
from templates import *
import sys
import time
import pickle


def check_name(name_org, nr, delta_list, BRA):
    """ for nam <aibj|ckdkcl>
    n0 and n1 indicate the boundaries of the
    right fixed side. The right part have indices
    intentionaly repeated. There is 10 different
    cases for triple indices: all different, two virt
    repeated, two occ repeated etc. This function checks 
    whether for given case some deltas try to 
    make two indices that are not supposed to be equal, equal.
    The 10 cases are:
    0. ckdlem
    1. ckdlek
    2. ckdlel
    3. and 4. ckclem
    5. and 6. ckdldm
    7. ckclek
    8. ckclel
    9. ckdldk
    10. ckdkdl
    """

    name = deepcopy(name_org)
    delta = deepcopy(delta_list)

    for i in range(0, len(delta)):
        for j in range(0, len(name)):
            if name[j] == delta[i][1]:
                name = name[0:j] + delta[i][0] + name[j + 1:]
        for j in range(i+1, len(delta)):
            for k in range(0, 2):
                if delta[j][k] == delta[i][1]:
                    delta[j][k] = delta[i][0]

    print(name_org)
    print(name)
    if BRA == 2:
        nv = name[4] + name[6] + name[8]
        no = name[5] + name[7] + name[9]
    elif BRA == 1:
        nv = name[2] + name[4] + name[6]
        no = name[3] + name[5] + name[7]
    if BRA == 3:
        nv = name[0] + name[2] + name[4]
        no = name[1] + name[3] + name[5]



    print(nv, no)
    is_zero = False
    if nr == 0: #ckdlem 
        if nv[0] == nv[1] or nv[0] == nv[2] or nv[1] == nv[2]:
            is_zero = True
            return is_zero
        if no[0] == no[1] or no[0] == no[2] or no[1] == no[2]:
            is_zero = True
            return is_zero
    if nr == 1: #ckdlek
        if nv[0] == nv[1] or nv[0] == nv[2] or nv[1] == nv[2]:
            is_zero = True
            return is_zero
        if no[0] == no[1] or no[1] == no[2]:
            is_zero = True
            return is_zero
    if nr == 2: #ckdlel
        if nv[0] == nv[1] or nv[0] == nv[2] or nv[1] == nv[2]:
            is_zero = True
            return is_zero
        if no[0] == no[1] or no[0] == no[2]:
            is_zero = True
            return is_zero
    if nr == 3: #ckclem 3 and 4
        if nv[0] == nv[2] or nv[1] == nv[2]:
            is_zero = True
            return is_zero
        if no[0] == no[1] or no[0] == no[2] or no[1] == no[2]:
            is_zero = True
            return is_zero
    if nr == 4: #ckdldm 5 and 6
        if nv[0] == nv[1] or nv[0] == nv[2]:
            is_zero = True
            return is_zero
        if no[0] == no[1] or no[0] == no[2] or no[1] == no[2]:
            is_zero = True
            return is_zero
    if nr == 5: #ckclek
        if nv[0] == nv[2] or nv[1] == nv[2]:
            is_zero = True
            return is_zero
        if no[0] == no[1] or no[1] == no[2]:
            is_zero = True
            return is_zero
    if nr == 6: #ckclel
        if nv[0] == nv[2] or nv[1] == nv[2]:
            is_zero = True
            return is_zero
        if no[0] == no[1] or no[0] == no[2]:
            is_zero = True
            return is_zero
    if nr == 7: #ckdldk
        if nv[0] == nv[1] or nv[0] == nv[2]:
            is_zero = True
            return is_zero
        if no[0] == no[1] or no[1] == no[2]:
            is_zero = True
            return is_zero
    if nr == 8: #ckdkdl
        if nv[0] == nv[2] or nv[0] == nv[1]:
            is_zero = True
            return is_zero
        if no[0] == no[2] or no[1] == no[2]:
            is_zero = True
            return is_zero

    return is_zero

def rename_name_fixed(name_org, dl, n0, n1):
    """executes delta in name of                                                                                                                         
    the function i.e. name = 'aibjck'                                                                                                          
    delta_list = ['b', 'a'], """

    name = deepcopy(name_org)
    ban = []
    
    for k in range(0, len(dl)):
        for j in range(n0, n1):
            print('teraz bede zmieniac indeks', name[j])
            if name[j] == dl[k][0] and name[j] not in ban:
                name = name[0:j] + dl[k][1] + name[j + 1:]
                print('zostaje on zmieniony na', dl[k][1])
                ban.append(name[j])
            elif name[j] == dl[k][1] and name[j] not in ban:
                name = name[0:j] + dl[k][0] + name[j + 1:]
                print('zostaje on zmieniony na', dl[k][0])
                ban.append(name[j])

    return name

def rename_name_fixed2(name_org, dl_org, n0, n1):
    """executes delta in name of                                                                                                                          
    the function i.e. name = 'aibjck'                                                                                                                                                   
    delta_list = ['b', 'a'], """

    name = deepcopy(name_org)
    ban = []
    dl = deepcopy(dl_org)
    dl_bra = []
    dl_mix = []

    for i in range(0, len(dl)):
        if dl[i][0] in name_org[n0:n1] and dl[i][1] in name_org[n0:n1]:
            dl_bra.append(dl[i])
        elif dl[i][0] in name_org[n0:n1] and dl[i][1] not in name_org[n0:n1]:
            dl_mix.append(dl[i])
        elif dl[i][0] not in name_org[n0:n1] and dl[i][1] in name_org[n0:n1]:
                dl_mix.append(dl[i])
        else:
            print('DELTA BTW KET ELEMENTS')
            #sys.exit(0)
#    print('dl', dl_org)
#    print('delta bra', dl_bra)
#    print('delta_mix', dl_mix)
    ban = []

    for i in range(n0, n1):
        for j in range(0, len(dl_bra)):
            dl_bra[j].sort()
            if name[i] == dl_bra[j][0] and name[i] not in ban:
 #               print('znalazlem element', i, name[i], 'w delcie')
                idx = name.index(dl_bra[j][1])
  #              print('zmieniam', dl_bra[j][1], 'na', dl_bra[j][0])
                name = name[0:idx] + dl_bra[j][0] + name[idx + 1:]
                ban.append(name[i])
                for k in range(0, len(dl_mix)):
                    if dl_mix[k][0] == dl_bra[j][1]:
   #                   print('w delcie mix zmieniam', dl_mix[k][0], 'na', dl_bra[j][0])
    #                  print('przed', dl_mix)
                      dl_mix[k][0] = dl_bra[j][0]
     #                 print('po', dl_mix)
                    elif dl_mix[k][1] == dl_bra[j][1]:
      #                print('w delcie mix zmieniam', dl_mix[k][1], 'na', dl_bra[j][0])
       #               print('przed', dl_mix, dl_mix[k][1])
                      dl_mix[k][1] = dl_bra[j][0]
          #            print('po', dl_mix, dl_mix[k][1])
                
    ban = []
#    print('wykonuje delte mix')
    for i in range(n0, n1):
        for j in range(0, len(dl_mix)):
            n = 0
            if name[i] == dl_mix[j][0] and name[i] not in ban:
 #               print('zmieniam', name[i], 'na', dl_mix[j][1])
                name = name[0:i] + dl_mix[j][1] + name[i + 1:]
  #              print(name)
                ban.append(dl_mix[j][1])
                n += 1
            if name[i] == dl_mix[j][1] and name[i] not in ban:
   #             print('zmieniam', name[i], 'na', dl_mix[j][0])
                name = name[0:i] + dl_mix[j][0] + name[i + 1:]
    #            print(name)
                ban.append(dl_mix[j][0])
                n += 1
            if n > 1:
                print('DELTA BTW KET ELEMENTS')
                #sys.exit(0)
    return name


def integrate_13(a, i, b, j, c, k):

    w = ugg()
    w.operator_idx.append([a, i])
    w.operator_idx.append([b, j])
    w.operator_idx.append([c, k])
    w.operator_type.append("s")
    w.operator_type.append("s")
    w.operator_type.append("s")
    
    r = evaluate(hamiltonian, w)
    
    rint = r.integrate(bra = ['a', 'i'], braspin = ['s']).scale(0.5)
    
    ars = preprep_for_fortran(rint)
    
    return ars

def integrate_23(a, i, b, j, c, k):

    w = ugg()
    w.operator_idx.append([a, i])
    w.operator_idx.append([b, j])
    w.operator_idx.append([c, k])
    w.operator_type.append("s")
    w.operator_type.append("s")
    w.operator_type.append("s")

    r = evaluate(hamiltoniant, w)

    rint1 = r.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s', 's'])
    rint2 = r.integrate(bra = ['a', 'j', 'b', 'i'], braspin = ['s', 's'])
    rint  = rint1.scale(1./3.) + rint2.scale(1./6.)

    rsimp = simplify(rint)
    rint = deepcopy(rsimp)

    r2 = deepcopy(rint)
    for x in r2:
        x.new_delta("a", "b")
        x.new_delta("i", "j")
    rint += r2.scale(-1./2.)

    ars = preprep_for_fortran(rint)

    return ars


def integrate_3120(r, a, i, b, j, c, k):

    print('integrate_3120 a')
    r01 = r.integrate(bra = [a, i, b, j, c, k], braspin = ['s', 's', 's']).scale(1./3.)
    print('integrate_3120 b')
    r02 = r.integrate(bra = [a, j, b, i, c, k], braspin = ['s', 's', 's']).scale(1./6.)
    r0 = r01 + r02
    ars0 = preprep_for_fortran(r0)

    return ars0

def integrate_3126(r, a, i, b, j, c, k):

    print('integrate_3126 a')
    r61 = r.integrate(bra = [a, i, b, j, c, k], braspin = ['s', 's', 's']).scale(1./3.)
    print('integrate_3126 b')
    r62 = r.integrate(bra = [a, k, b, j, c, i], braspin = ['s', 's', 's']).scale(1./6.)
    r6 = r61 + r62
    ars6 = preprep_for_fortran(r6)

    return ars6

def integrate_31206(r, a, i, b, j, c, k):
    print('integrate here')
    r06 = r.integrate(bra = [a, i, b, j, c, k], braspin = ['s', 's', 's']).scale(1./2.)

    ars06 = preprep_for_fortran(r06)

    return ars06

def integrate_312_prim(r, a, i, b, j, c, k):

    rint = r.integrate(bra = [a, i, b, j, c, k], braspin = ['s', 's', 's'])
    ars = preprep_for_fortran(rint)

    return ars


def integrate_3121(r, a, i, b, j, c, k):

    print('integrate_3121 1')
    rint1 = r.integrate(bra = [a, i, b, j, c, k], braspin = ['s', 's', 's']).scale(1./4.)
    print('integrate_3121 2')
    rint2 = r.integrate(bra = [a, k, b, i, c, j], braspin = ['s', 's', 's']).scale(1./12.)
    print('integrate_3121 3')
    rint3 = r.integrate(bra = [a, k, b, j, c, i], braspin = ['s', 's', 's']).scale(1./6.)
    print('integrate_3121 4')
    rint4 = r.integrate(bra = [a, j, b, i, c, k], braspin = ['s', 's', 's']).scale(1./6)
    print('integrate_3121 5')
    rint5 = r.integrate(bra = [a, j, b, k, c, i], braspin = ['s', 's', 's']).scale(1./12.)
    print('integrate_3121 6')
    r1 = rint1 + rint2 + rint3 + rint4 + rint5
    ars1 = preprep_for_fortran(r1)

    return ars1

def integrate_3122(r, a, i, b, j, c, k):

    print('integrate_3122 1')
    rint1 = r.integrate(bra = [a, j, b, k, c, i], braspin = ['s', 's', 's']).scale(1./12.)
    print('integrate_3122 2')
    rint2 = r.integrate(bra = [a, i, b, j, c, k], braspin = ['s', 's', 's']).scale(1./4.)
    print('integrate_3122 3')
    rint3 = r.integrate(bra = [a, i, b, k, c, j], braspin = ['s', 's', 's']).scale(1./6.)
    print('integrate_3122 4')
    rint4 = r.integrate(bra = [a, k, b, j, c, i], braspin = ['s', 's', 's']).scale(1./6.)
    print('integrate_3122 5')
    rint5 = r.integrate(bra = [a, k, b, i, c, j], braspin = ['s', 's', 's']).scale(1./12.)
    r2 = rint1 + rint2 + rint3 + rint4 + rint5
    ars2 = preprep_for_fortran(r2)

    return ars2

def integrate_3123(r, a, i, b, j, c, k):

    print('integrate_3123 1')
    rint1 = r.integrate(bra = [a, k, b, j, c, i], braspin = ['s', 's', 's']).scale(1./6.)
    print('integrate_3123 2')
    rint2 = r.integrate(bra = [a, i, b, k, c, j], braspin = ['s', 's', 's']).scale(1./6.)
    print('integrate_3123 3')
    rint3 = r.integrate(bra = [a, i, b, j, c, k], braspin = ['s', 's', 's']).scale(1./3.)
    print('integrate_3123 4')
    rint4 = r.integrate(bra = [a, j, b, k, c, i], braspin = ['s', 's', 's']).scale(1./6.)
    print('integrate_3123 5')
    rint5 = r.integrate(bra = [a, j, b, i, c, k], braspin = ['s', 's', 's']).scale(1./6.)
    r3 = rint1 + rint2 + rint3 + rint4 + rint5
    ars3 = preprep_for_fortran(r3)

    return ars3

def integrate_3124(r, a, i, b, j, c, k):

    print('integrate_3124 1')
    rint1 = r.integrate(bra = [a, j, b, i, c, k], braspin = ['s', 's', 's']).scale(1./6.)
    print('integrate_3124 2')
    rint2 = r.integrate(bra = [a, k, b, j, c, i], braspin = ['s', 's', 's']).scale(1./6.)
    print('integrate_3124 3')
    rint3 = r.integrate(bra = [a, k, b, i, c, j], braspin = ['s', 's', 's']).scale(1./6.)
    print('integrate_3124 4')
    rint4 = r.integrate(bra = [a, i, b, j, c, k], braspin = ['s', 's', 's']).scale(1./3.)
    print('integrate_3124 5')
    rint5 = r.integrate(bra = [a, i, b, k, c, j], braspin = ['s', 's', 's']).scale(1./6.)
    r4 = rint1 + rint2 + rint3 + rint4 + rint5
    ars4 = preprep_for_fortran(r4)

    return ars4

def integrate_3125(r, a, i, b, j, c, k):
    
    print('integrate_3125 1')
    rint1 = r.integrate(bra = [a, k, b, i, c, j], braspin = ['s', 's', 's']).scale(1./12.)
    print('integrate_3125 2')
    rint2 = r.integrate(bra = [a, j, b, k, c, i], braspin = ['s', 's', 's']).scale(1./12.)
    print('integrate_3125 3')
    rint3 = r.integrate(bra = [a, j, b, i, c, k], braspin = ['s', 's', 's']).scale(1./6.)
    print('integrate_3125 4')
    rint4 = r.integrate(bra = [a, i, b, k, c, j], braspin = ['s', 's', 's']).scale(1./6.)
    print('integrate_3125 5')
    rint5 = r.integrate(bra = [a, i, b, j, c, k], braspin = ['s', 's', 's']).scale(1./4.)
    r5 = rint1 + rint2 + rint3 + rint4 + rint5
    ars5 = preprep_for_fortran(r5)

    return ars5



def task_eom_mem(BRA, KET,theory, trans, pick, noR, restricted = False):

    if BRA == 3:
        task_eom_triple_mem(BRA, KET, 'cc3', 'trans', pick, noR, restricted)
        return

    print(pick)

    open_flags = ""
    if pick == "dump":
        open_flags = "wb"        
    elif pick == "load":
        open_flags = "rb"


    if BRA == 1:
        if KET == 3:
            out = open('./pickle/ars13_mem.pkl', open_flags)
            outn = open('./pickle/name_13_mem.pkl', open_flags)
            outz = open('./pickle/z_name_13_mem.pkl', open_flags)
        else:
            print('ONLY KET 3 ALLOWED')
            sys.exit(1)

    elif BRA == 2:
        if KET == 3:
            out = open('./pickle/ars23_mem.pkl', open_flags)
            outn = open('./pickle/name_23_mem.pkl', open_flags)
            outz = open('./pickle/z_name_23_mem.pkl', open_flags)
        else:
            print('ONLY KET 3 ALLOWED')
            sys.exit(1)


    if pick == "dump":
        if BRA == 1:
            if KET == 3:
                
                ars_list = []
                z_name = [['z0', 'bjckdl'], ['z1', 'bjckdj'], ['z2', 'bjckdk'], ['z34', 'bjbkdl'], \
                          ['z56', 'bjckcl'], ['z7', 'bjbkdj'], ['z8', 'bjbkdk'], ['z9', 'bjckcj'], ['z10', 'bjcjck']]
                # z0
                ars = integrate_13('b', 'j', 'c', 'k', 'd', 'l')
                ars_list.append(deepcopy(ars))
                # z1
                ars = integrate_13('b', 'j', 'c', 'k', 'd', 'j')
                ars_list.append(deepcopy(ars))
                # z2
                ars = integrate_13('b', 'j', 'c', 'k', 'd', 'k')
                ars_list.append(deepcopy(ars))
                # z3 and z4
                ars = integrate_13('b', 'j', 'b', 'k', 'd', 'l')
                ars_list.append(deepcopy(ars))
                # z5 and z_6
                ars = integrate_13('b', 'j', 'c', 'k', 'c', 'l')
                ars_list.append(deepcopy(ars))
                # z7
                ars = integrate_13('b', 'j', 'b', 'k', 'd', 'j')
                ars_list.append(deepcopy(ars))
                # z8
                ars = integrate_13('b', 'j', 'b', 'k', 'd', 'k')
                ars_list.append(deepcopy(ars))
                # z9                                                                                                                                           
                ars = integrate_13('b', 'j', 'c', 'k', 'c', 'j')
                ars_list.append(deepcopy(ars))
                # z10                                                                                                                                        
                ars = integrate_13('b', 'j', 'c', 'j', 'c', 'k')
                ars_list.append(deepcopy(ars))
                print('z10')
                for x in ars:
                    print(x)


            if noR == True:
                name_list = []
                for y in range(0, len(ars_list)):
                    name = []
                    #print('teraz bedzie, ', y, z_name[y][0])
                    for x in range(0, len(ars_list[y])):
                        print(ars_list[y][x])
                        base_name = 'ai' + z_name[y][1]
                        #print(check_name(base_name, y, ars_list[y][x].delta, BRA))
                        if (check_name(base_name, y, ars_list[y][x].delta, BRA)):
                            ars_list[y][x].num_factor = 0
                        # ars_list[y][x].coefficient.append(EOM_CC_SINGLE_Rl)
                        # ars_list[y][x].coefficient_idx.append(['a', 'i'])
                        # ars_list[y][x].summation.append('a')
                        # ars_list[y][x].summation.append('i')
                        ars_list[y][x].constraints['a'] = set('a')
                        ars_list[y][x].constraints['i'] = set('i')
                        # print(ars_list[y][x])
                        # base_name = 'ai'
                        # sname  = rename_name_fixed2(base_name, ars_list[y][x].delta, 0, 2)
                        # name.append(sname)
    #                    ars_list[y][x].name = deepcopy(sname)
    #                    print('sname', sname)
    #                    print(base_name, sname, ars_list[y][x].delta)
                        # ars_list[y][x].exec_delta()
                        # print('sname', sname, ars_list[y][x])
                        # print('')

                    name_list.append(name)

                    ars_list[y].cleanup()

                    for x in range(0, len(ars_list[y])):
                        print('po clean', ars_list[y][x])
            elif noR == False:
                name_list = []
                for y in range(0, len(ars_list)):
                    name = []
                    for x in range(0, len(ars_list[y])):
                        print(ars_list[y][x])
                        ars_list[y][x].coefficient.append(EOM_CC_SINGLE_Rl)
                        ars_list[y][x].coefficient_idx.append(['a', 'i'])
                        ars_list[y][x].summation.append('a')
                        ars_list[y][x].summation.append('i')
                        ars_list[y][x].constraints['a'] = set('a')
                        ars_list[y][x].constraints['i'] = set('i')
                        base_name = 'ai'
                        sname  = rename_name_fixed2(base_name, ars_list[y][x].delta, 0, 2)
                        name.append(sname)
                        ars_list[y][x].name = deepcopy(sname)
                        ars_list[y][x].exec_delta()

                    name_list.append(name)

                    ars_list[y].cleanup()

                    for x in range(0, len(ars_list[y])):
                        print('po clean', ars_list[y][x])


            pickle.dump(ars_list, out)
            pickle.dump(name_list, outn)
            pickle.dump(z_name, outz)
            sys.exit(0)
        elif BRA == 2:
            if KET == 3:

                ars_list = []
                
                z_name = [['z0', 'ckdlem'], ['z1', 'ckdlek'], ['z2', 'ckdlel'], ['z34', 'ckclem'], \
                          ['z56', 'ckdldm'], ['z7', 'ckclek'], ['z8', 'ckclel'], ['z9', 'ckdldk'], ['z10', 'ckdkdl']]

                print('dump z0')
                ars = integrate_23('c', 'k', 'd', 'l', 'e', 'm')
                ars_list.append(deepcopy(ars))
                # z1            
                print('dump z1')
                ars = integrate_23('c', 'k', 'd', 'l', 'e', 'k')
                ars_list.append(deepcopy(ars))
                # z2          
                print('dump z2')
                ars = integrate_23('c', 'k', 'd', 'l', 'e', 'l')
                ars_list.append(deepcopy(ars))
                # z3 and z4   
                print('dump z34')
                ars = integrate_23('c', 'k', 'c', 'l', 'e', 'm')
                ars_list.append(deepcopy(ars))
                # z5 and z6    
                print('dump z56')
                ars = integrate_23('c', 'k', 'd', 'l', 'd', 'm')
                ars_list.append(deepcopy(ars))
                # z7           
                print('dump z7')
                ars = integrate_23('c', 'k', 'c', 'l', 'e', 'k')
                ars_list.append(deepcopy(ars))
                # z8          
                print('dump z8')
                ars = integrate_23('c', 'k', 'c', 'l', 'e', 'l')
                ars_list.append(deepcopy(ars))
                # z9          
                print('dump z9')
                ars = integrate_23('c', 'k', 'd', 'l', 'd', 'k')
                ars_list.append(deepcopy(ars))
                # z10         
                print('dump z10')
                ars = integrate_23('c', 'k', 'd', 'k', 'd', 'l')
                ars_list.append(deepcopy(ars))

            if noR == True:
                name_list = []
                for y in range(0, len(ars_list)):
                    name = []
                    print('teraz sprawdzam case', y)
                    for x in range(0, len(ars_list[y])):
                        #print('znalazlem element', x, ars_list[y][x])
                        print(ars_list[y][x])
                        base_name = 'aibj' + z_name[y][1] 
                        if (check_name(base_name, y, ars_list[y][x].delta, BRA)):
                            ars_list[y][x].num_factor = 0
                        #print(check_name(base_name, y, ars_list[y][x].delta, BRA))
                        # ars_list[y][x].coefficient.append(EOM_CC_SINGLE_Rl)
                        # ars_list[y][x].coefficient_idx.append(['a', 'i', 'b', 'j'])
                        # ars_list[y][x].summation.append('a')
                        # ars_list[y][x].summation.append('i')
                        # ars_list[y][x].summation.append('b')
                        # ars_list[y][x].summation.append('j')
                        ars_list[y][x].constraints['a'] = set('a')
                        ars_list[y][x].constraints['i'] = set('i')
                        ars_list[y][x].constraints['b'] = set('b')
                        ars_list[y][x].constraints['j'] = set('j')
                        # base_name = 'aibj'
                        # sname  = rename_name_fixed2(base_name, ars_list[y][x].delta, 0, 4)
                        # name.append(sname)
                        # ars_list[y][x].name = deepcopy(sname)
                        # #print(base_name, sname, ars_list[y][x].delta)
                        # print('sname', sname, ars_list[y][x])
                        # ars_list[y][x].exec_delta()
                        # if len(ars_list[y][x].delta) > 0:
                        #     ars_list[y][x].num_factor = 0
                        # print(ars_list[y][x])
                        print('')
                    ars_list[y].cleanup()
                    for z in ars_list[y]:
                        print(z)
                    print('')
                    name_list.append(name)
            elif noR == False:
                name_list = []
                for y in range(0, len(ars_list)):
                    name = []
                    for x in range(0, len(ars_list[y])):
                        print(ars_list[y][x])
                        ars_list[y][x].coefficient.append(EOM_CC_SINGLE_Rl)
                        ars_list[y][x].coefficient_idx.append(['a', 'i', 'b', 'j'])
                        ars_list[y][x].summation.append('a')
                        ars_list[y][x].summation.append('i')
                        ars_list[y][x].summation.append('b')
                        ars_list[y][x].summation.append('j')
                        ars_list[y][x].constraints['a'] = set('a')
                        ars_list[y][x].constraints['i'] = set('i')
                        ars_list[y][x].constraints['b'] = set('b')
                        ars_list[y][x].constraints['j'] = set('j')
                        if not restricted:
                            ars_list[y][x].num_factor = 1./2. * ars_list[y][x].num_factor
                        base_name = 'aibj'
                        sname  = rename_name_fixed2(base_name, ars_list[y][x].delta, 0, 4)
                        name.append(sname)
                        ars_list[y][x].name = deepcopy(sname)
                        # #print(base_name, sname, ars_list[y][x].delta)
                        # print('sname', sname, ars_list[y][x])
                        ars_list[y][x].exec_delta()
                        if len(ars_list[y][x].delta) > 0:
                            ars_list[y][x].num_factor = 0
                        print(ars_list[y][x])
                        print('')
                    ars_list[y].cleanup()
                    for z in ars_list[y]:
                        print(z)
                    print('')
                    name_list.append(name)

            pickle.dump(ars_list, out)
            pickle.dump(name_list, outn)
            pickle.dump(z_name, outz)
            sys.exit(0)
    elif pick == "load":
        ars_list = pickle.load(out)
        name_list = pickle.load(outn)
        z_name = pickle.load(outz)

        k = 1
        for x in ars_list:
            k+= 1
            for y in x:
                if BRA == 1:
                    y.loops = ['a', 'i']
                    y.ibra = ['a', 'i']
                    print(y, y.ibra)
                elif BRA == 2:
                    y.loops = ['a', 'i', 'b', 'j']
                    y.ibra = ['a', 'i', 'b', 'j']
                elif BRA == 3:
                    if KET == 1:
                        y.loops = ['d', 'l']
                        y.iket = ['d', 'l']
                    elif KET == 2:
                        y.loops = ['d', 'l', 'e', 'm']
                        y.iket = ['d', 'l', 'e', 'm']
                print(y, y.constraints)
                y.exec_delta_fixed_mem()
                print(y, y.constraints, y.loops, y.ibra)
                print('---------------------------------------------------------------------------------------------------------------')
            print('')



        # print('sniez sniez')
        # for x in range(0, len(ars_list[7])):
        #     print(x, ars_list[7][x], ars_list[7][x].constraints, ars_list[7][x].loops)


    # ars_list_clean = []
    # name_list_clean = []
    # for x in range(0, len(ars_list)):
    #     ars_mini = arithmetic_string()
    #     name_mini = arithmetic_string()
    #     for y in range(0, len(ars_list[x])):
    #         if ars_list[x][y].num_factor != 0:
    #             ars_mini.append(ars_list[x][y])
    #             name_mini.append(name_list[x][y])
    #     ars_list_clean.append(ars_mini)
    #     name_list_clean.append(name_mini)
    k = 0
    for x in ars_list:
        #x.cleanup()
        #x.as_standarize()
        print('', k)
        for y in x:
            print(y, y.ibra)
        k+=1



    eom_func_mem(ars_list, name_list, z_name, BRA, KET ,theory, trans, triplet = "", lmfold = "", rmfold = "", noR=noR)
    
    # eom_func_mem(ars_list_clean, name_list_clean, z_name, BRA, KET ,theory, trans)
    name_list_clean = []



def task_eom_triple_mem(BRA, KET, theory, trans, pick, noR, restricted=False):
    open_flags = ""
    if pick == "dump":
        open_flags = "wb"        
    elif pick == "load":
        open_flags = "rb"


    if KET == 1:
        o0 = open('./pickle/ars0_mem.pkl', open_flags)
        o6 = open('./pickle/ars6_mem.pkl', open_flags)
        o06 = open('./pickle/ars06_mem.pkl', open_flags)
        oprim = open('./pickle/arsp_mem.pkl', open_flags)

        on0 = open('./pickle/name0_mem.pkl', open_flags)
        on6 = open('./pickle/name6_mem.pkl', open_flags)
        on06 = open('./pickle/name06_mem.pkl', open_flags)
        onprim = open('./pickle/namep_mem.pkl', open_flags)

        on0z = open('./pickle/z_name0_mem.pkl', open_flags)
        on6z = open('./pickle/z_name6_mem.pkl', open_flags)
        on06z = open('./pickle/z_name06_mem.pkl', open_flags)
        onprimz = open('./pickle/z_namep_mem.pkl', open_flags)


    elif KET == 2:
        o0 = open('./pickle/ars0_32_mem.pkl', open_flags)
        o6 = open('./pickle/ars6_32_mem.pkl', open_flags)
        o06 = open('./pickle/ars06_32_mem.pkl', open_flags)
        oprim = open('./pickle/arsp_32_mem.pkl', open_flags)

        on0 = open('./pickle/name0_32_mem.pkl', open_flags)
        on6 = open('./pickle/name6_32_mem.pkl', open_flags)
        on06 = open('./pickle/name06_32_mem.pkl', open_flags)
        onprim = open('./pickle/namep_32_mem.pkl', open_flags)

        on0z = open('./pickle/z_name0_32_mem.pkl', open_flags)
        on6z = open('./pickle/z_name6_32_mem.pkl', open_flags)
        on06z = open('./pickle/z_name06_32_mem.pkl', open_flags)
        onprimz = open('./pickle/z_namep_32_mem.pkl', open_flags)

    if pick == 'dump':
        
        if KET == 1:
            r = evaluate(hamiltoniant, t2, nudl)
        elif KET == 2:
            r = evaluate(hamiltoniant, nudlem)
        ars_list = []

        ars_list_big = []
        name_list_big = []
        z_name_big = []
        ars_list = []

        # z1
        ars = integrate_3120(deepcopy(r), 'a', 'i', 'b', 'j', 'c', 'i')
        ars_list.append(deepcopy(ars))
        # z2                                                                                   
        ars = integrate_3120(deepcopy(r), 'a', 'i', 'b', 'j', 'c', 'j')
        ars_list.append(deepcopy(ars))
        print('pooooooooooooo')
        z_name = [['z1', 'aibjci'], ['z2', 'aibjcj']]

        name_list = []
        for y in range(0, len(ars_list)):

            if noR == False:

                if KET == 1:
                    name = []
                    for x in range(0, len(ars_list[y])):
                        print('przed', ars_list[y][x])
                        ars_list[y][x].coefficient.append(EOM_CC_SINGLE_Rr)
                        ars_list[y][x].coefficient_idx.append(['d', 'l'])
                        ars_list[y][x].summation.append('d')
                        ars_list[y][x].summation.append('l')
                        ars_list[y][x].constraints['d'] = set('d')
                        ars_list[y][x].constraints['l'] = set('l')
                        print('po', ars_list[y][x], ars_list[y][x].constraints)
                        if y == 0:
                            base_name = 'aibjcidl'
                        elif y == 1:
                            base_name = 'aibjcjdl'
                        sname  = rename_name_fixed2(base_name, ars_list[y][x].delta, 6, 8)
                        print('sname', sname)
                        name.append(sname)
                        ars_list[y][x].name = deepcopy(sname)
                        ars_list[y][x].exec_delta()
                        print(ars_list[y][x])
                    name_list.append(name)
                elif KET == 2:
                    name = []
                    for x in range(0, len(ars_list[y])):
                        print('przed', ars_list[y][x])
                        ars_list[y][x].coefficient.append(EOM_CC_SINGLE_Rr)
                        ars_list[y][x].coefficient_idx.append(['d', 'l', 'e', 'm'])
                        ars_list[y][x].summation.append('d')
                        ars_list[y][x].summation.append('l')
                        ars_list[y][x].summation.append('e')
                        ars_list[y][x].summation.append('m')
                        if not restricted:
                            ars_list[y][x].num_factor = 1./2. * ars_list[y][x].num_factor
                        ars_list[y][x].constraints['d'] = set('d')
                        ars_list[y][x].constraints['l'] = set('l')
                        ars_list[y][x].constraints['e'] = set('e')
                        ars_list[y][x].constraints['m'] = set('m')
                        print('po', ars_list[y][x], ars_list[y][x].constraints)

                        if y == 0:
                            base_name = 'aibjcidlem'
                        elif y == 1:
                            base_name = 'aibjcjdlem'
                        sname  = rename_name_fixed2(base_name, ars_list[y][x].delta, 6, 10)
                        print('sname', sname, ars_list[y][x])
                        name.append(sname)
                        ars_list[y][x].name = deepcopy(sname)
                        ars_list[y][x].exec_delta()
                        print(ars_list[y][x])
                        print('')
                    name_list.append(name)
            elif noR == True:
                if KET == 1:
                    name = []
                    for x in range(0, len(ars_list[y])):
                        print('przed', ars_list[y][x])
                        name = []
                        if y == 0:
                            base_name = 'aibjcidl'
                            z = 1
                        elif y == 1:
                            base_name = 'aibjcjdl'
                            z = 2
                        if (check_name(base_name, z, ars_list[y][x].delta, BRA)):
                            ars_list[y][x].num_factor = 0
                        ars_list[y][x].constraints['d'] = set('d')
                        ars_list[y][x].constraints['l'] = set('l')
                        ars_list[y][x].exec_delta()
                        print(ars_list[y][x])
                    name_list.append(name)

                elif KET == 2:
                    name = []
                    for x in range(0, len(ars_list[y])):
                        name = []
                        print('przed', ars_list[y][x])
                        if y == 0:
                            base_name = 'aibjcidlem'
                            z = 1
                        elif y == 1:
                            base_name = 'aibjcjdlem'
                            z = 2
                        if (check_name(base_name, z, ars_list[y][x].delta, BRA)):
                            ars_list[y][x].num_factor = 0
                        ars_list[y][x].constraints['d'] = set('d')
                        ars_list[y][x].constraints['l'] = set('l')
                        ars_list[y][x].constraints['e'] = set('e')
                        ars_list[y][x].constraints['m'] = set('m')
                        ars_list[y][x].exec_delta()
                    name_list.append(name)

        print('dumping 0')
        pickle.dump(deepcopy(ars_list), o0)
        pickle.dump(deepcopy(name_list), on0)
        pickle.dump(deepcopy(z_name), on0z)
        ars_list_big.append(deepcopy(ars_list))
        name_list_big.append(deepcopy(name_list))
        z_name_big.append(deepcopy(z_name))
        ars_list = []

        # z3 and z4                                                                                                                            
        ars = integrate_3126(deepcopy(r), 'a', 'i', 'a', 'j', 'c', 'k')
        ars_list.append(deepcopy(ars))
        # z5 and z6                                                                                                                
        ars = integrate_3126(deepcopy(r), 'a', 'i', 'b', 'j', 'b', 'k')
        ars_list.append(deepcopy(ars))

        z_name = [['z34', 'aiajck'], ['z56', 'aibjbk']]
        name_list = []
        if noR == False:
            for y in range(0, len(ars_list)):
                if KET == 1:
                    name = []
                    for x in range(0, len(ars_list[y])):
                        ars_list[y][x].coefficient.append(EOM_CC_SINGLE_Rr)
                        ars_list[y][x].coefficient_idx.append(['d', 'l'])
                        ars_list[y][x].summation.append('d')
                        ars_list[y][x].summation.append('l')
                        ars_list[y][x].constraints['d'] = set('d')
                        ars_list[y][x].constraints['l'] = set('l')
                        if y == 0:
                            base_name = 'aiajckdl'
                        elif y == 1:
                            base_name = 'aibjbkdl'
                        sname  = rename_name_fixed2(base_name, ars_list[y][x].delta, 6, 8)
                        name.append(sname)
                        ars_list[y][x].name = deepcopy(sname)
                        ars_list[y][x].exec_delta()
                    name_list.append(name)
                elif KET == 2:
                    name = []
                    for x in range(0, len(ars_list[y])):
                        ars_list[y][x].coefficient.append(EOM_CC_SINGLE_Rr)
                        ars_list[y][x].coefficient_idx.append(['d', 'l', 'e', 'm'])
                        ars_list[y][x].summation.append('d')
                        ars_list[y][x].summation.append('l')
                        ars_list[y][x].summation.append('e')
                        ars_list[y][x].summation.append('m')
                        if not restricted:
                            ars_list[y][x].num_factor = 1./2. * ars_list[y][x].num_factor
                        ars_list[y][x].constraints['d'] = set('d')
                        ars_list[y][x].constraints['l'] = set('l')
                        ars_list[y][x].constraints['e'] = set('e')
                        ars_list[y][x].constraints['m'] = set('m')
                        if y == 0:
                            base_name = 'aiajckdlem'
                        elif y == 1:
                            base_name = 'aibjbkdlem'
                        sname  = rename_name_fixed2(base_name, ars_list[y][x].delta, 6, 10)
                        ars_list[y][x].name = deepcopy(sname)
                        print('sname', sname, ars_list[y][x])
                        name.append(sname)
                        ars_list[y][x].exec_delta()
                    name_list.append(name)
        elif noR == True:
            for y in range(0, len(ars_list)):
                if KET == 1:
                    name = []
                    for x in range(0, len(ars_list[y])):
                        if y == 0:
                            base_name = 'aiajckdl'
                            z = 3
                        elif y == 1:
                            base_name = 'aibjbkdl'
                            z = 4
                        if (check_name(base_name, z, ars_list[y][x].delta, BRA)):
                            ars_list[y][x].num_factor = 0

                        ars_list[y][x].constraints['d'] = set('d')
                        ars_list[y][x].constraints['l'] = set('l')
                        ars_list[y][x].exec_delta()
                    name_list.append(name)
                elif KET == 2:
                    name = []
                    for x in range(0, len(ars_list[y])):
                        if y == 0:
                            base_name = 'aiajckdlem'
                            z = 3
                        elif y == 1:
                            base_name = 'aibjbkdlem'
                            z = 4
                        if (check_name(base_name, z, ars_list[y][x].delta, BRA)):
                            ars_list[y][x].num_factor = 0
                        ars_list[y][x].constraints['d'] = set('d')
                        ars_list[y][x].constraints['l'] = set('l')
                        ars_list[y][x].constraints['e'] = set('e')
                        ars_list[y][x].constraints['m'] = set('m')
                        ars_list[y][x].exec_delta()
                    name_list.append(name)


        print('dumping 6')
        pickle.dump(deepcopy(ars_list), o6)
        pickle.dump(deepcopy(name_list), on6)
        pickle.dump(deepcopy(z_name), on6z)
        ars_list_big.append(deepcopy(ars_list))
        name_list_big.append(deepcopy(name_list))
        z_name_big.append(z_name)
        ars_list = []
        name_list = []

        # z7                                                                                                                                  
        ars = integrate_31206(deepcopy(r), 'a', 'i', 'a', 'j', 'c', 'i')
        ars_list.append(deepcopy(ars))
        # z8                                                                                                                                           
        ars = integrate_31206(deepcopy(r), 'a', 'i', 'a', 'j', 'c', 'j')
        ars_list.append(deepcopy(ars))
        ars = integrate_31206(deepcopy(r), 'a', 'i', 'b', 'j', 'b', 'i')
        ars_list.append(deepcopy(ars))
        # z10                                                                                                         
        ars = integrate_31206(deepcopy(r), 'a', 'i', 'b', 'i', 'b', 'j')
        ars_list.append(deepcopy(ars))

        z_name = [['z7', 'aiajci'], ['z8', 'aiajcj'], ['z9', 'aibjbi'], ['z10', 'aibibj']]

        name_list = []
        if noR == False:
            for y in range(0, len(ars_list)):
                if KET == 1:
                    name = []
                    for x in range(0, len(ars_list[y])):
                        ars_list[y][x].coefficient.append(EOM_CC_SINGLE_Rr)
                        ars_list[y][x].coefficient_idx.append(['d', 'l'])
                        ars_list[y][x].summation.append('d')
                        ars_list[y][x].summation.append('l')
                        ars_list[y][x].constraints['d'] = set('d')
                        ars_list[y][x].constraints['l'] = set('l')
                        if y == 0:
                            base_name = 'aiajcidl'
                        elif y == 1:
                            base_name = 'aiajcjdl'
                        elif y == 2:
                            base_name = 'aibjbidl'
                        elif y == 3:
                            base_name = 'aibibjdl'
                        sname  = rename_name_fixed2(base_name, ars_list[y][x].delta, 6, 8)
                        name.append(sname)
                        print('przed', ars_list[y][x], ars_list[y][x].constraints)
                        ars_list[y][x].name = deepcopy(sname)
                        ars_list[y][x].exec_delta()
                        print('po', ars_list[y][x], ars_list[y][x].constraints, sname)
                    name_list.append(name)
                elif KET == 2:
                    name = []
                    for x in range(0, len(ars_list[y])):
                        ars_list[y][x].coefficient.append(EOM_CC_SINGLE_Rr)
                        ars_list[y][x].coefficient_idx.append(['d', 'l', 'e', 'm'])
                        ars_list[y][x].summation.append('d')
                        ars_list[y][x].summation.append('l')
                        ars_list[y][x].summation.append('e')
                        ars_list[y][x].summation.append('m')
                        if not restricted:
                            ars_list[y][x].num_factor = 1./2. * ars_list[y][x].num_factor
                        ars_list[y][x].constraints['d'] = set('d')
                        ars_list[y][x].constraints['l'] = set('l')
                        ars_list[y][x].constraints['e'] = set('e')
                        ars_list[y][x].constraints['m'] = set('m')
                        if y == 0:
                            base_name = 'aiajcidlem'
                        elif y == 1:
                            base_name = 'aiajcjdlem'
                        elif y == 2:
                            base_name = 'aibjbidlem'
                        elif y == 3:
                            base_name = 'aibibjdlem'
                        sname  = rename_name_fixed2(base_name, ars_list[y][x].delta, 6, 10)
                        print('sname', sname, ars_list[y][x])
                        ars_list[y][x].name = deepcopy(sname)
                        name.append(sname)
                        ars_list[y][x].exec_delta()
                    name_list.append(name)
        elif noR == True:
            for y in range(0, len(ars_list)):
                if KET == 1:
                    name = []
                    for x in range(0, len(ars_list[y])):
                        if y == 0:
                            base_name = 'aiajcidl'
                            z = 5
                        elif y == 1:
                            base_name = 'aiajcjdl'
                            z = 6
                        elif y == 2:
                            base_name = 'aibjbidl'
                            z = 7
                        elif y == 3:
                            base_name = 'aibibjdl'
                            z = 8
                        if (check_name(base_name, z, ars_list[y][x].delta, BRA)):
                            ars_list[y][x].num_factor = 0

                        ars_list[y][x].constraints['d'] = set('d')
                        ars_list[y][x].constraints['l'] = set('l')
                        ars_list[y][x].exec_delta()
                    name_list.append(name)
                elif KET == 2:
                    name = []
                    for x in range(0, len(ars_list[y])):
                        if y == 0:
                            base_name = 'aiajcidlem'
                            z = 5
                        elif y == 1:
                            base_name = 'aiajcjdlem'
                            z = 6
                        elif y == 2:
                            base_name = 'aibjbidlem'
                            z = 7
                        elif y == 3:
                            base_name = 'aibibjdlem'
                            z = 8
                        if (check_name(base_name, z, ars_list[y][x].delta, BRA)):
                            ars_list[y][x].num_factor = 0
                        ars_list[y][x].constraints['d'] = set('d')
                        ars_list[y][x].constraints['l'] = set('l')
                        ars_list[y][x].constraints['e'] = set('e')
                        ars_list[y][x].constraints['m'] = set('m')
                        ars_list[y][x].exec_delta()
                    name_list.append(name)

        print('dumping 06')
        pickle.dump(deepcopy(ars_list), o06)
        pickle.dump(deepcopy(name_list), on06)
        pickle.dump(deepcopy(z_name), on06z)
        ars_list_big.append(deepcopy(ars_list))
        name_list_big.append(deepcopy(name_list))
        z_name_big.append(z_name)
        ars_list = []
        name_list = []

        ars = integrate_312_prim(deepcopy(r), 'a', 'i', 'b', 'j', 'c', 'k')
        ars_list.append(deepcopy(ars))
        z_name = [['z0', 'aibjck']]
        name_list = []
        if noR == False:
            for y in range(0, len(ars_list)):
                if KET == 1:
                    name = []
                    for x in range(0, len(ars_list[y])):
                        ars_list[y][x].coefficient.append(EOM_CC_SINGLE_Rr)
                        ars_list[y][x].coefficient_idx.append(['d', 'l'])
                        ars_list[y][x].summation.append('d')
                        ars_list[y][x].summation.append('l')
                        ars_list[y][x].constraints['d'] = set('d')
                        ars_list[y][x].constraints['l'] = set('l')
                        base_name = 'aibjckdl'
                        sname  = rename_name_fixed2(base_name, ars_list[y][x].delta, 6, 8)
                        ars_list[y][x].name = deepcopy(sname)
                        name.append(sname)
                        ars_list[y][x].exec_delta()
                    name_list.append(name)
                elif KET == 2:
                    name = []
                    for x in range(0, len(ars_list[y])):
                        ars_list[y][x].coefficient.append(EOM_CC_SINGLE_Rr)
                        ars_list[y][x].coefficient_idx.append(['d', 'l', 'e', 'm'])
                        ars_list[y][x].summation.append('d')
                        ars_list[y][x].summation.append('l')
                        ars_list[y][x].summation.append('e')
                        ars_list[y][x].summation.append('m')
                        if not restricted:
                            ars_list[y][x].num_factor = 1./2. * ars_list[y][x].num_factor
                        ars_list[y][x].constraints['d'] = set('d')
                        ars_list[y][x].constraints['l'] = set('l')
                        ars_list[y][x].constraints['e'] = set('e')
                        ars_list[y][x].constraints['m'] = set('m')
                        base_name = 'aibjckdlem'
                        sname  = rename_name_fixed2(base_name, ars_list[y][x].delta, 6, 10)
                        ars_list[y][x].name = deepcopy(sname)
                        print('sname', sname, ars_list[y][x])
                        name.append(sname)
                        ars_list[y][x].exec_delta()
                    name_list.append(name)
        if noR == True:
            for y in range(0, len(ars_list)):
                if KET == 1:
                    name = []
                    for x in range(0, len(ars_list[y])):
                        base_name = 'aibjckdl'
                        z = 0
                        if (check_name(base_name, z, ars_list[y][x].delta, BRA)):
                            ars_list[y][x].num_factor = 0
                        ars_list[y][x].constraints['d'] = set('d')
                        ars_list[y][x].constraints['l'] = set('l')
                        ars_list[y][x].exec_delta()
                    name_list.append(name)
                elif KET == 2:
                    name = []
                    for x in range(0, len(ars_list[y])):
                        base_name = 'aibjckdlem'
                        z = 0
                        if (check_name(base_name, z, ars_list[y][x].delta, BRA)):
                            ars_list[y][x].num_factor = 0
                        ars_list[y][x].constraints['d'] = set('d')
                        ars_list[y][x].constraints['l'] = set('l')
                        ars_list[y][x].constraints['e'] = set('e')
                        ars_list[y][x].constraints['m'] = set('m')
                        ars_list[y][x].exec_delta()
                    name_list.append(name)

        print('dumping prim')
        pickle.dump(deepcopy(ars_list), oprim)
        pickle.dump(deepcopy(name_list), onprim)
        pickle.dump(deepcopy(z_name), onprimz)
        ars_list_big.append(deepcopy(ars_list))
        name_list_big.append(deepcopy(name_list))
        z_name_big.append(z_name)
    
        sys.exit(0)

    elif pick == 'load':
        ars0 = pickle.load(o0)
        ars6 = pickle.load(o6)        
        ars06 = pickle.load(o06)
        arsp = pickle.load(oprim)

        # for x in ars0:
        #     for y in x:
        #         print(y)
        # sys.exit(0)

        n0 = pickle.load(on0)
        n6 = pickle.load(on6)
        n06 = pickle.load(on06)
        np = pickle.load(onprim)
        
        z0 = pickle.load(on0z)
        z6 = pickle.load(on6z)
        z06 = pickle.load(on06z)
        zp = pickle.load(onprimz)

        ars_list_big = [ars0, ars6, ars06, arsp]
        name_list_big = [n0, n6, n06, np]
        z_name_big = [z0, z6, z06, zp]


        k = 1
        for x in ars_list_big:
            print('kolejny')
            for y in x:
                for z in y:
                    if BRA == 3:
                        if KET == 1:
                            z.loops = ['d', 'l']
                            z.iket = ['d', 'l']
                        elif KET == 2:
                            z.loops = ['d', 'l', 'e', 'm']
                            z.iket = ['d', 'l', 'e', 'm']
                        print(z, z.constraints, z.iket, z.loops)
                        z.exec_delta_fixed_mem()
                    print(z, z.constraints, z.loops, z.iket)
                print('---------------------------------------------------------------------------------------------------------------')
                print('')

    print('przed standarize')
    for x in range(0, len(ars0)):
        print('la')
        for y in ars0[x]:
            print(y, y.name, y.constraints, y.delta)
        print('')



        #    print('przed standarize', ars06[2][10])

    for x in ars_list_big:
        for y in x:
            y.establish_fixed()
            y.as_standarize()

    for x in ars0:
        for y in x:
            print(y)


    # print('pluszon2')
    # for x in range(0, len(ars06)):
    #     for y in ars06[x]:
    #         print(y, y.name)
    #     print('')
    # print('koniec paldus_eom_mem')

    
    print('po standarize')
    for x in range(0, len(ars06)):
        for y in ars06[x]:
            print(y, y.name, y.constraints, y.delta)
        print('')

    eom_func_triple_mem(ars_list_big, name_list_big, z_name_big, BRA, KET, theory, trans, noR)


   
def integrate_triple_exc_bra(a1, a2, a3, a4, a5, r0):
    
    r = deepcopy(r0)
    
    rint1 = r.integrate(bra = ['a', 'i', 'b', 'k', 'c', 'j'], braspin = ['s', 's', 's']).scale(a1)
    rint2 = r.integrate(bra = ['a', 'j', 'b', 'i', 'c', 'k'], braspin = ['s', 's', 's']).scale(a2)
    rint3 = r.integrate(bra = ['a', 'j', 'b', 'k', 'c', 'i'], braspin = ['s', 's', 's']).scale(a3)
    rint4 = r.integrate(bra = ['a', 'k', 'b', 'i', 'c', 'j'], braspin = ['s', 's', 's']).scale(a4)
    rint5 = r.integrate(bra = ['a', 'k', 'b', 'j', 'c', 'i'], braspin = ['s', 's', 's']).scale(a5)

    rint = rint1 + rint2 + rint3 + rint4 + rint5

    return rint
