#!/usr/bin/env python
import sys
from itertools import *
import math
from copy import deepcopy
from paldus_basic import evaluate
from paldus_classes import * 
from templates import *
from params import *
from paldus_basic import simplify

T_op  = 101
T_opd = 102
S_op  = 103
S_opd = 104



def comb(i, maxexc):

    # Find all possible commutator combinations for given nesting = i
    # and for given maximum excitation = maxexc
    # [...[[X, T_a1],T_a2]...T_a3]_i  where aj \in {1, 2, ..., maxexc}
    #
    # returns i-length list [a1, a2, ..., a3]
    # print('w comb', maxexc)

    lst = list(combinations_with_replacement(range(1, maxexc), i))
#    print('lst', lst)
    for j in range(0, len(lst)):
        lst[j] = list(lst[j])

    lst2 = []
    for j in range(0, len(lst)):
        if (4 not in lst[j]):# and 3 not in lst[j]):
            lst2.append(lst[j])

#    print('lst2', lst2)
    return lst2

def particle(lst):

    # Compute the particle number of the given list
    # if list = [a1, a2, ... a3], 
    # then suma = a1 + a2 + ... + a3

    suma = 0
    for i in range(0, len(lst)):
        if lst[i] == 'tf':
            suma = suma + 2
        elif lst[i] == 'sf':
            suma = suma + 2
        else:
            suma = suma + lst[i]

    return suma

def mbpt_o(amp):

    # returns the lowest mbpt order of given
    # amplitude contriution

    if amp == 1:
        mbpt = 2
    elif amp == 2:
        mbpt = 1
    elif amp == 3:
        mbpt = 2
    else:
        print(amp)
        print('print from paldus commutators - program higher amplitudes and mbpt orders')
        sys.exit(0)
    return mbpt

def mbpt_order(lst_lst):

    # return the mbpt order of the product of all of the
    # operators in the list of lists.
    #
    # lst_lst = [[a1, ...], [b1, ...], ... [c1, ...]]
    # mbpt_order = mbpt_o(a1) + ... +  mbpt_o(b1) + ... +  
    #              mbpt_o(c1) + ...
    suma = 0
    print('s', lst_lst)
    for i in range(0, len(lst_lst)):
        for j in range(0, len(lst_lst[i])):
            suma = suma + mbpt_o(lst_lst[i][j])
    return suma


def check_mw(mw):

    # particle rank: mw, must at least 1

    app = True
    if (mw <= 0):
        app = False
        return app
    return app

def check_swg(sw, cond):

    # check if sw > cond

    app = True
    if( sw <= cond):
        app = False
        return app
    return app

def check_swge(sw, cond):

    # check if sw >= cond

    app = True
    if( sw < cond):
        app = False
        return app
    return app

def check_swl(sw, cond):

    # check if sw < cond

    app = True
    if( sw >= cond):
        app = False
        return app
    return app

def check_swle(sw, cond):

    # check if sw =< cond

    app = True
    if( sw > cond):
        app = False
        return app
    return app

def check_conditions_Wm3(il, ns, n, theory):

    # Check conditions for the expression
    # P(e^(S*) \mu_n e^(-S*))
    # See pdf ....

    # W1 = e^(S*) \mu_n e^(-S*)

    # 1a) 2 * ns - il >= 0
    # 1b) 2 * n - il >= 0
    # 3) mw1 >= 1
    # 4) 3 >= sw1 >=1 CC3 and 2>= sw1 >=1 CCSD

    # 1a)

    if (2 * ns - il) < 0 :
        return False, 0

    # 1b)

    if (2 * n - il) < 0:
        return False, 0

   # 3)

    mw1 = n + ns - il

    if (mw1 < 1):
        return False, 0

    # 4)

    sw1 = n - ns

    if (sw1 < 1):
#        print('d', sw1)
        return False, 0
    # if theory == 'ccsd':
    #     if sw1 > 2:
    #         return False, 0
    # elif theory == 'cc3':
    #     if sw1 > 3:
    #         return False, 0

    return True, sw1

def check_conditions_Wm1(ik, nt, il, ns, n, mx, theory):

    # Check conditions for the expression
    # e^(-S)e^(T*) X e^(-T*) e^(S)
    # See pdf ....

    # W1 = e^(T*) X e^(-T*)
    # W2 = e^(-S) W1 e^(S)

    # 1a) 2 * nt - ik + mx - n >= 0
    # 1b) mx - k + n >= 0
    # 2a) 2 * nt + mx - ik - n - il >=0
    # 2b) 2 * ns - il + mx - ik + n >= 0
    # 4) mw1 >= 1 mw2 >= 1


    # # 1a)
                                                                                                                                                              
    if (2 * nt - ik + mx - n) < 0 :
        print('1a')
        return False, 0

    # # 1b)
    

    if (mx - ik + n) < 0:
        print(mx, ik, n)
        print('1b')
        return False, 0

    # # 2a)                                                                                                                                                             

    if (2 * nt + mx - ik - n - il)  < 0:
        print('2a')
        return False, 0

    # # # 2b)
                                                                                                                                                              
    if (2 * ns - il + mx - ik + n) < 0:
        print('2b')
        return False, 0

    # 4a) 4b)                                                                                                                                                         

    mw1 = mx + nt - ik
    mw2 = mw1 + ns - il
    sw2 = n  - nt + ns
    print('----------', n, nt, ns)
    if (mw1 < 1):
        print('4a')
        return False, 0
    if (mw2 < 1):
        print('4b')
        return False, 0

    return True, sw2


def check_conditions_Wm_ground(k, nt, l, ns, sx, mx, list_t, list_s, theory):

    # Check conditions for the expression
    # e^(S*)e^(-T) X e^(T) e^(-S*)
    # See pdf ....

    # W1 = e^(-T) X e^(T)
    # W2 = e^(S*) W1 e^(-S*)

    #1) k <= 2mx
    #2) l <= 2 mW1 = 2mx + 2nt -2k
    #    2mx + 2nt -2k -l >= 0
    #3) mW1 - nt + k>= 0 bo mx >=0 to wynika z definicji (13.2.51)
    #4) for this [X, Tn1], Tn2]...] check if each inner commutator
    # is an pure excitation commutator. If so, the outer commutator is a commutator
    # of two excitation operators, and is 0. Pure excitation is when mW1 = sW1.
    # sg = excitation rank of [X, Tn1] and so on. sg <=0 to not disappear.
    # for deectiation [X, Sn1*] sg >=0.


    mW1 = mx + nt - k
    mW2 = mW1 - ns -l
    sW2 = sx  + nt - ns

    #0) nt>=k
    if (k > nt):
        # print('zagniezdzenie wieksze od wzbudzenia', k, nt)
        return False, 0
    #1)
    if(k > 2 * mx):
        # print('k jest za duze, k=', k)
        return False, 0
    #2)
    if (l > 2*mW1):
        # print('l jest za duze, l=', l)
        return False, 0

    if (sW2 != 0):
        # print('sW2 jest ! 0', sW2)
        return False, 0

    s_op = sx
    mw_op = mx
    for i in range(0, len(list_t)):
        if mw_op == s_op:
            # print('komutator dwoch operatorow ekscytacji')
            # print(i, mw_op, s_op, list_t, list_s)
            return False, 0
        else:
            s_op += list_t[i]
            mw_op += list_t[i] -1
            
        # if i != len(list_t)-1:
        #     if mg == sg:
        #         print('komutator dwoch operatorow ekscytacji')
        #         print(i, mg, sx, sg, list_t, list_s)
        #         return False, 0
        # else:
        #     if s_op > 0:
        #         print('komutator dwoch operatorow ekscytacji')
        #         print(i, mg, sx, sg, list_t, list_s)
        #         return False, 0

        # s_op = deepcopy(sg)
        # mw_op = deepcopy(mg)
        
    for i in range(0, len(list_s)):
        if mw_op == -s_op:
            # print('komutator dwoch operatorow deekscytacji')
            # print(i, mw_op, s_op, list_t, list_s)
            return False, 0
        else:
            s_op -= list_s[i]
            mw_op += list_s[i] -1


        
    # for i in range(0, len(list_s)):
    #     sg = s_op - list_s[i]
    #     mg = mw_op + list_s[i]-1
    #     if i != len(list_s)-1:
    #         if mg == -sg:
    #             print('komutator dwoch operatorow deekscytacji')
    #             print(i, mg, sx, sg, list_t, list_s)
    #             return False, 0
    #     else:
    #         if s_op < 0:
    #             print('komutator dwoch operatorow deekscytacji')
    #             print(i, mg, sx, sg, list_t, list_s)
    #             return False, 0

    #     s_op = deepcopy(sg)
    #     mw_op = deepcopy(mg)

    

    return True, sW2

def check_conditions_Wm_ground_f12(k, nt, l, ns, sx, mx, list_t, list_s, theory, Gm=0):

    # Check conditions for the expression
    # e^(S*)e^(-T) X e^(T) e^(-S*)
    # See pdf ....

    # W1 = e^(-T) X e^(T)
    # W2 = e^(S*) W1 e^(-S*)

    #1) k <= 2mx
    #2) l <= 2 mW1 = 2mx + 2nt -2k
    #    2mx + 2nt -2k -l >= 0
    #3) mW1 - nt + k>= 0 bo mx >=0 to wynika z definicji (13.2.51)
    #4) for this [X, Tn1], Tn2]...] check if each inner commutator
    # is an pure excitation commutator. If so, the outer commutator is a commutator
    # of two excitation operators, and is 0. Pure excitation is when mW1 = sW1.
    # sg = excitation rank of [X, Tn1] and so on. sg <=0 to not disappear.
    # for deectiation [X, Sn1*] sg >=0.


    mW1 = mx + nt - k
    

    mW2 = mW1 - ns -l

    sW2 = sx  + nt - ns


    #1)
    if(k > 2 * mx):
        #print('k jest za duze, k=', k)
        return False, 0

    #2)
    if (l > 2*mW1):
        #print('l jest za duze, l=', l)
        return False, 0

    if Gm == 0:
        if (sW2 != 0):
            #print('sW2 jest ! 0', sW2)
            return False, 0
    elif Gm == 1:
        if (sW2 > -1):
            #print('sW2 jest > -1', sW2)
            return False, 0
    elif Gm == 2:
        if (sW2 < 1):
            #print('sW2 jest < 1', sW2)
            return False, 0
        # else:
        #     print('sW2 jest wieksze ', sW2)

    s_op = sx
    mw_op = mx
    for i in range(0, len(list_t)):            
        if mw_op == s_op:
            #print('komutator dwoch operatorow ekscytacji')
            #print(i, mw_op, s_op, list_t, list_s)
            return False, 0
        else:
            if list_t[i] == 'tf':
                s_op += 0
                mw_op += 2 -1
            else:
                s_op += list_t[i]
                mw_op += list_t[i] -1

                
        # if i != len(list_t)-1:
        #     if mg == sg:
        #         #print('komutator dwoch operatorow ekscytacji')
        #         print(i, mg, sx, sg, list_t, list_s)
        #         return False, 0
        # else:
        #     if s_op > 0:
        #         print('komutator dwoch operatorow ekscytacji')
        #         print(i, mg, sx, sg, list_t, list_s)
        #         return False, 0

        # s_op = deepcopy(sg)
        # mw_op = deepcopy(mg)
        
    for i in range(0, len(list_s)):
        if mw_op == -s_op:
            #print('komutator dwoch operatorow deekscytacji')
            # print(i, mw_op, s_op, list_t, list_s)
            return False, 0
        else:
            if list_s[i] == 'sf':
                s_op -= 0
                mw_op += 2 -1
            else:
                s_op -= list_s[i]
                mw_op += list_s[i] -1

    return True, sW2


def check_conditions_Wm_ground_ccd(k, nt, l, ns, sx, mx, list_t, list_s, theory, Gm=0):

    # Check conditions for the expression
    # <0|e^(S*)e^(-T) X e^(T) e^(-S*)P(e^(S*)e^(-T) X e^(T) e^(-S*))|0>
    
    # Gm2 = P(e^(S*)e^(-T) X e^(T) e^(-S*))
    # Gm1 = e^(S*)e^(-T) X e^(T) e^(-S*)
    # Gm0 = <0|e^(S*)e^(-T) X e^(T) e^(-S*)|0>
    # See pdf ....

    # W1 = e^(-T) X e^(T)
    # W2 = e^(S*) W1 e^(-S*)

    #1) k <= 2mx
    #2) l <= 2 mW1 = 2mx + 2nt -2k
    #    2mx + 2nt -2k -l >= 0
    #3) mW1 - nt + k>= 0 bo mx >=0 to wynika z definicji (13.2.51)
    #4) for this [X, Tn1], Tn2]...] check if each inner commutator
    # is an pure excitation commutator. If so, the outer commutator is a commutator
    # of two excitation operators, and is 0. Pure excitation is when mW1 = sW1.
    # sg = excitation rank of [X, Tn1] and so on. sg <=0 to not disappear.
    # for deectiation [X, Sn1*] sg >=0.


    mW1 = mx + nt - k
    

    mW2 = mW1 - ns -l

    sW2 = sx  + nt - ns


    #1)
    if(k > 2 * mx):
        # print('k jest za duze, k=', k)
        return False, 0

    #2)
    if (l > 2*mW1):
        # print('l jest za duze, l=', l)
        return False, 0

    if Gm == 0:
        if (sW2 != 0):
            # print('sW2 jest ! 0', sW2)
            return False, 0
    elif Gm == 1:
        if (sW2 > -1):
            # print('sW2 jest > -1', sW2)
            return False, 0
    elif Gm == 2:
        if (sW2 < 1):
            # print('sW2 jest < 1', sW2)
            return False, 0
        else:
            print('sW2 jest wieksze ', sW2)

    s_op = sx
    mw_op = mx
    for i in range(0, len(list_t)):            
        if mw_op == s_op:
            # print('komutator dwoch operatorow ekscytacji')
            # print(i, mw_op, s_op, list_t, list_s)
            return False, 0
        else:
            s_op += list_t[i]
            mw_op += list_t[i] -1
        
    for i in range(0, len(list_s)):
        if mw_op == -s_op:
            # print('komutator dwoch operatorow deekscytacji')
            # print(i, mw_op, s_op, list_t, list_s)
            return False, 0
        else:
            s_op -= list_s[i]
            mw_op += list_s[i] -1

    return True, sW2


def check_conditions_Wm2(ik, nt, il, ns, n, theory):

    # Check conditions for the expression
    # P(e^(-S)e^(T*) \mu_n e^(-T*) e^(S))
    # See pdf ....
    
    # W1 = e^(T*) \mu_n e^(-T*)
    # W2 = e^(-S) W1 e^(S)
    
    # 1a) 2 * nt - ik >= 0
    # 1b) 2 * n - ik >= 0
    # 2a) 2 * nt - il - ik >= 0
    # 2b) 2 * n + 2 * ns - il - ik >= 0
    # 3a)  if nt == 0 and ns != 0 --> return False
    # 3b)  if nt == 0 and ns == 0   3 >= n >= 1 CC3 and 2 >= n >= 1 CCSD
    # 4) mw1 >= 1 mw2 >= 1
    # 5) 3 >= sw2 >=1 CC3 and 2>= sw2 >=1 CCSD


    # 1a)

    if (2 * nt - ik) < 0 :
        return False, 0
    
    # 1b)

    if (2 * n - ik) < 0:
        return False, 0

    # 2a)

    if (2 * nt - il -ik)  < 0:
        return False, 0

    # 2b)

    if (2 * n + 2 * ns - ik - il) < 0:
        return False, 0

    # 3a) 3b)

    if (nt == 0 and ns != 0):
        return False, 0
        
    if (nt == 0 and ns == 0):
        if n < 1:
            return False, 0
        if theory == 'ccsd':
            if n > 2:
                return False, 0
        elif theory == 'cc3':
            if n > 3:
                return False, 0

    # 4a) 4b)
    
    mw1 = n + nt - ik
    mw2 = mw1 + ns - il

    if (mw1 < 1):
        return False, 0
    if (mw2 < 1):
        return False, 0

    # 5a)

    sw2 = n - nt + ns

    if (sw2 < 1):
        return False, 0
    # if theory == 'ccsd':
    #     if sw2 > 2:
    #         return False, 0
    # elif theory == 'cc3':
    #     if sw2 > 3:
    #         return False, 0

    return True, sw2


def check_conditions_gamma1(ik, nt, il, ns, n, mx):

    # Check conditions for the expression
    # e(S*)e(-T) X e(T) e(-S*)
    #
    # The particle rank mw1 for w1 = e^(-T) \mu_n e^(T)
    # must be greater than 0

    # The particle rank mw2 for w2 = e^(S*) w1 e^(-S*)
    # must be greater than 0
    
    # The excitation rank of w2, must be  in (-3, 0)

    # Returns app = True if conditions are satisfied

    # ik and il represent the order of nesting for
    # T and S operators respectively
    
    # mx is the particle rank of X

    # nt and ns represent the total de(excitation)
    # of T and S operators respectively

    # n is given

    # The expressions for sw and mw are from
    # "Molecular Electronic Structure Theory" T. Helgaker
    # vol. 2, pp.664, eq. 13.2.54 and 13.2.55
    
    
    sw = n + nt - ns
    mw1 = mx + nt - ik
    mw2 = mw1 + ns - il

    # if ns = 0, and nt != 0 and n>0 than X is the excitation operator
    # and commutes with T
    if (ns == 0 and nt != 0 and n > 0):
        app = False
        return app, sw

    # if nt = 0, and ns != 0 and n<0 than X is the deexcitation operator
    # and commutes with S*
    if (nt == 0 and ns != 0 and n < 0):
        app = False
        return app, sw

    app1 = check_mw(mw1)
    app2 = check_mw(mw2)

    # check if sW =< 0
    app3 = check_swle(sw, 0)
    # check if sW >= -3
    app4 = check_swge(sw, -3)

    app = True
    if(app1 == False or app2 == False or app3 == False or app4 == False):
        app = False

    return app, sw


# def check_conditions_Wm3(ik, ns, n):

#     # Check conditions for the expression
#     # P(e^(S*) \mu_n e^(-S*))
#     #
#     # The particle mw rank for w = e^(S*) \mu_n e^(-S*)
#     # must be greater than 0

#     # The excitation rank of w, must be greater
#     # than 0, due to the projector P.

#     # Returns app = True if conditions are satisfied

#     # ik represent the order of nesting for
#     # S operators 

#     # ns represent the total deexcitation
#     # of S operators 

#     # n is given

#     # The expressions for sw and mw are from
#     # "Molecular Electronic Structure Theory" T. Helgaker
#     # vol. 2, pp.664, eq. 13.2.54 and 13.2.55

#     sw = n - ns
#     mw = n + ns - ik


#     app1 = check_mw(mw)

#     # check if sW > 0
#     app2 = check_swg(sw, 0)    
#     # check if sW <= 3
#     app3 = check_swle(sw, 3)    

#     app = True
#     if(app1 == False or app2 == False or app3 == False):
#         app = False

#     return app, sw


def check_conditions_Wlr(ik, nt, n, mx):

    # Check conditions for the expression
    # e^(-T) Y e^(T)
    #
    # The particle mw rank for w = e^(-T) \mu_n e^(T)
    # must be greater than 0

    # The excitation rank of w, must be either
    # 1, 2 or 3, ad these are the only possible
    # excitations in ket

    # Returns app = True if conditions are satisfied

    # ik represent the order of nesting for
    # T operators 

    # mx is the particle rank of X
    
    # nt represent the total excitation
    # of T operators 

    # n is given

    # The expressions for sw and mw are from
    # "Molecular Electronic Structure Theory" T. Helgaker
    # vol. 2, pp.664, eq. 13.2.54 and 13.2.55

    sw = n + nt
    mw = mx + nt - ik

    # if ns = 0, and nt != 0 and n>0 than X is an excitation operator 
    # and commutes with T
    if (nt != 0 and n > 0):
        app = False
        return app, sw

    app1 = check_mw(mw)

    # check if sW >= 1
    app2 = check_swge(sw, 1)    
    app3 = check_swle(sw, 3)    

    app = True
    if(app1 == False or app2 == False or app3 == False):
        app = False

    return app, sw

# def check_conditions_Wm1(ik, nt, il, ns, n, mx):


#     # Check conditions for the expression
#     # e^(-S)e^(T*) X e^(-T*) e^(S)
#     #
#     # The particle rank mw1 for w1 = e^(T*) X e^(-T*)
#     # must be greater than 0

#     # The particle rank mw2 for w2 = e^(-S) w1 e^(S)
#     # must be greater than 0
    
#     # Returns app = True if conditions are satisfied

#     # ik and il represent the order of nesting for
#     # T and S operators respectively

#     #  mx is the particle rank of X

#     # nt and ns represent the total de(excitation)
#     # of T and S operators respectively

#     # n is given

#     # The expressions for sw and mw are from
#     # "Molecular Electronic Structure Theory" T. Helgaker
#     # vol. 2, pp.664, eq. 13.2.54 and 13.2.55

#     sw1 = n - nt
#     mw1 = mx + nt - ik
# #    sw1m = (mw1 - sw1) / 2
#     # sw1m (2 * nt -ik) /  2 the division is removed intentionally                                                                                                                      
#     sw1m = (2 * nt - ik)


#     # check the 13.2.59 condition where A is the result of [A, T1*]...T1*]                                                                               
#     #    if (sw1m < il):
#     if (sw1m < 0):
#         app = False
#         return app, 0


#     sw = n - nt + ns
#     mw1 = mx + nt - ik    
#     mw2 = mw1 + ns - il

#     # if ns = 0, and nt != 0 and n<0 than X is the deexcitation operator
#     # and commutes with T
#     if (ns == 0 and nt != 0 and n < 0):
# #        print('false1')
#         app = False
#         return app, sw

#     # if nt = 0, and ns != 0 and n>0 than X is the excitation operator
#     # and commutes with S
#     if (nt == 0 and ns != 0 and n > 0):
# #        print('false2')
#         app = False
#         return app, sw

#     app1 = check_mw(mw1)
#     app2 = check_mw(mw2)

#     app = True
#     if(app1 == False or app2 == False):
# #        print('false3', mw1, mw2)
#         app = False

#     return app, sw



# Gm =  <e(S*)e(-T) E_pq e(T) e(-S*) P(e(S*)e(-T) E_rs e(T) e(-S*))>                                                         
        #--------Gm1----------#  #-------------Gm2-------------#  

        
# def generate_Gm(Gm1, Gm2, maxmbpt):
#     # Find all possible Gm

#     print('GENERATE_GM dupa', maxmbpt)

#     Gm = []

#     jj = 0
#     for a in range(0, len(Gm1)):
#         for b in range(0, len(Gm2)):
#             sw_Gm1 = Gm1[a]['sw'] 
#             sw_Gm2 = Gm2[b]['sw']
#             if (sw_Gm1 == -sw_Gm2):
#                 jj += 1
#                 print('generuje Gm1 i Gm2 i ich rzedy', Gm1[a]['mbpt'] , Gm2[b]['mbpt'] )
#                 mbpt_t = Gm1[a]['mbpt'] + Gm2[b]['mbpt'] 
#                 minidict = {}
#                 minidict['Gm1_S_list'] = Gm1[a]['S_list']
#                 minidict['Gm1_T_list'] = Gm1[a]['T_list']
#                 minidict['Gm2_S_list'] = Gm2[b]['S_list']
#                 minidict['Gm2_T_list'] = Gm2[b]['T_list']
#                 minidict['Gm2_P_list'] = Gm2[b]['sw']
#                 minidict['mbpt_Gm1'] = Gm1[a]['mbpt']
#                 minidict['mbpt_Gm2'] = Gm2[b]['mbpt']
#                 print(minidict)
#                 # if cumulative == True:
#                 #     if (mbpt_t <= maxmbpt):
#                 #         Gm.append(minidict)
#                 # else:
#                 if (mbpt_t == maxmbpt):
#                     print('dodaje', minidict)
#                     Gm.append(minidict)

#     return Gm





# <e(-T)Ye(T)|\mu_m>    <e(-S)e(T*) X e(-T*) e(S) P(e(-S)e(T*) \mu_n e(-T*) e(S))| P(e(S*)\mu_l e(-S*))>    <mu_p|e(-T)Ye(T)>
#-----Wl-----------#   #--------------------------------------Wm-----------------------------------#    #------Wr-------#

# Wm =  <e(-S)e(T*) X e(-T*) e(S) P(e(-S)e(T*) \mu_n e(-T*) e(S))| P(e(S*)\mu_l e(-S*))> 
          #--------Wm1----------#  #-------------Wm2-------------#  #-------Wm3-------#


# Um = <  P(e(S*)\mu_l e(-S*)) | e(-S)e(T*) X e(-T*) e(S) P(e(-S)e(T*) \mu_n e(-T*) e(S))> 
          #--------Wm3---------# #--------Wm1-----------# #----------Wm2---------------#

# Wm = Um*




# Find all possible commutator combinations for Wm2
# for given mbpt order

def generate_Wm2(maxpt, theory):

    # Wm2 = P (e^(-S) e^(T*) \mu_n e^(-T*) e^S)

    Wm2 = []

    if (theory == 'ccsd'):
        maxexc = 2
        maxmu = 2
    elif (theory == 'cc3'):
        maxexc = 3
        maxmu = 3


    # zamienic maxexc na maxpt

    for ik in range(0, maxpt+1):  # loop over nesting of T*
        list_t = comb(ik, maxexc+1)
        lent = len(list_t)
        #print('t', list_t)

        for il in range(0, maxpt+1):  # loop over nesing of S
            if (ik + il) <= maxpt:
                list_s = comb(il, maxexc+1)
                lens = len(list_s)
                #print('s', list_s)

                for jk in range(0, lent):
                    for jl in range(0, lens):

                        # nt and ns are the total excitation rank of
                        # all the operatrors T and S respectively

                        nt = particle(list_t[jk])
                        ns = particle(list_s[jl])


                        for n in range(1, maxmu + 1):
                            #print('sprawdzam warunki dla', list_t[jk], list_s[jl], n)
                            app, sw  = check_conditions_Wm2(ik, nt, il, ns, n, theory)
                            mbpt = mbpt_order([list_t[jk], list_s[jl]])
                            if mbpt > maxpt:
                                app = False

                            if app == True:
                                #print('-----', [list_t[jk], list_s[jl]])
                                mbpt = mbpt_order([list_t[jk], list_s[jl]])
                                minidict = {}
                                minidict['T_list'] = list_t[jk]
                                minidict['S_list'] = list_s[jl]
                                minidict['n'] = n
                                minidict['sw'] = sw
                                minidict['mbpt'] = mbpt
                                #print(minidict)
                                Wm2.append(minidict)
                            # else:
                            #        print('')
                            #        print('nie dodaje bo sw = ', sw)
                            #        print('')
                            #        print('------------------')

#    sys.exit(0)
    # for x in Wm2:
    #     print(x)

    return Wm2


def generate_Wm3(maxpt, theory):


    # Find possible commutator combinations for Wm3

    if (theory == 'ccsd'):
        maxexc = 2
        maxmu = 2
    elif (theory == 'cc3'):
        maxexc = 3
        maxmu = 3

                
    Wm3 = []

    for ik in range(0, 4):

        list_s = comb(ik, maxexc + 1)
        lens = len(list_s)

        for jk in range(0, lens):

            # ns is the total excitation rank of
            # all the operatrors S 

            ns = particle(list_s[jk])
 #           print('lllllllllllllllll', list_s[jk], ns)
            for n in range(1, maxmu  + 1):
                app, sw = check_conditions_Wm3(ik, ns, n, theory)
#                print('---', n, sw, app)
                mbpt = mbpt_order([list_s[jk]])
#                print('mbpt', mbpt)
                if mbpt > maxpt:
                    app = False

                if app == True:
                    mbpt = mbpt_order([list_s[jk]])
                    #                    if (mbpt <= mbpt_max):
                    minidict = {}
                    minidict['S_list'] = list_s[jk]
                    minidict['n'] = n
                    minidict['sw'] = sw
                    minidict['mbpt'] = mbpt
                    Wm3.append(minidict)

    # print('Wm3')                                                                                                                                       
    # for i in Wm3:                                                                                                                                      
    #     print(i)                                                                                                                                       
    # sys.exit(0)  
    return Wm3

def generate_Wm1(maxpt, theory):

    #print('WM1')

    # Find possible commutator combinations for Wm1

    if (theory == 'ccsd'):
        maxexc = 2
        maxmu = 2
    elif (theory == 'cc3'):
        maxexc = 3
        maxmu = 3


    Wm1 = []

    for ik in range(0, maxpt+1):

        list_t = comb(ik, maxexc+1)
        lent = len(list_t)

        for il in range(0, maxpt+1):
            list_s = comb(il, maxexc+1)
            lens = len(list_s)

            if (ik + il) <= maxpt:
                for jk in range(0, lent):
                    for jl in range(0, lens):
                        print('--------------------------------------------------', 't=', list_t[jk], 's=', list_s[jl])
                        # nt and ns are the total excitation rank of
                        # all the operatrors T and S respectively

                        nt = particle(list_t[jk])
                        ns = particle(list_s[jl])


                        # Below n is only in (-1, 0, 1) as Y in Wm1 is
                        # one electron operator
                        for n in range(-1, 2):

                            mbpt = mbpt_order([list_t[jk], list_s[jl]])
                            #if (mbpt <= mbpt_max):                                                                                                       
                            minidict = {}
                            minidict['T_list'] = list_t[jk]
                            minidict['S_list'] = list_s[jl]
                            minidict['n'] = n

                            app, sw  = check_conditions_Wm1(ik, nt, il, ns, n, 2, theory)
                            print('app', app, sw)
                            minidict['sw'] = sw
                            minidict['mbpt'] = mbpt
                            if (mbpt <= maxpt):
                                    print(minidict)


                            if app == True:
                                Wm1.append(minidict)


    print('Wm1')
    for i in Wm1:
        print(i)

    return Wm1


def comp_min_pt(list_s, list_t):

    pt = 0
    for i in list_s:
        if i == 0:
            pt += 0
        elif i == 1:
            pt += 2
        elif i == 2:
            pt += 1
        elif i == 3:
            pt += 2
        elif i == 4:
            pt += 3
        elif i == 'sf':
            pt += 1
        else:
            print('print from paldus_commutators- add higher mbpt orders')
            sys.exit(0)
    for i in list_t:
        if i == 0:
            pt += 0
        elif i == 1:
            pt += 2
        elif i == 2:
            pt += 1
        elif i == 3:
            pt += 2
        elif i == 4:
            pt += 3
        elif i == 'tf':
            pt += 1
        else:
            print('print from paldus_commutators- add higher mbpt orders')
            sys.exit(0)                        
    return pt


def generate_Wm_ground_mbpt(maxpt, theory):
    #Wm_ground < 0 | e(S*)e(-T) X e(T) e(-S*) | 0>
              #--------Wm_ground-----------# 

    #print('WM_ground')

    # Find possible commutator combinations for Wm_ground

    if (theory == 'ccsd'):
        maxexc = 2
    elif (theory == 'cc3'):
        maxexc = 3

    # maxexc = 4
    Wm_ground = []

    for il in range(0, maxpt+1):
        list_s = comb(il, maxexc+1)
        lens = len(list_s)
        # print('list s', il, list_s)
        for ik in range(0, maxpt+1):

            list_t = comb(ik, maxexc+1)
            lent = len(list_t)
            # print('list t', ik, list_t)
                        
            for jk in range(0, lent):
                for jl in range(0, lens):
                    mbpt = comp_min_pt(list_s[jl], list_t[jk])

                    # print('------teraz sprawdzam----------------------', 't=', list_t[jk], 's=', list_s[jl])
                    # print('ich pt =', mbpt, maxpt)

                    if mbpt == maxpt:
                        # print('------teraz sprawdzam----------------------', 't=', list_t[jk], 's=', list_s[jl])
                        # nt and ns are the total excitation rank of
                        # all the operatrors T and S respectively
                        # print('tak mniejsze i licze')
                        nt = particle(list_t[jk])
                        ns = particle(list_s[jl])

                        # Below n is only in (-1, 0, 1) as Y in Wm_ground is
                        # one electron operator
                        for n in range(-1, 2):
                        
                            minidict = {}
                            minidict['T_list'] = list_t[jk]
                            minidict['S_list'] = list_s[jl]
                            minidict['n'] = n
                            # print('')
                            # print('sprawdzam zestaw', n, list_t[jk], nt, list_s[jl], ns)
                            # print('')                                                        
                            app, sw  = check_conditions_Wm_ground(ik, nt, il, ns, n, 1, list_t[jk], list_s[jl], theory)                                                        
                            # print('check conditions, dodaje?', app, sw)
                            # print('')
                            minidict['sw'] = sw
                            minidict['mbpt'] = mbpt

                            if app == True:
                                Wm_ground.append(minidict)


    # print('Wm_ground')
    # for i in Wm_ground:
    #     print(i)

    return Wm_ground

def t2tfsplit(lst, new_list):

    print('')
    print('lst  in new_list na wejsciu', lst, new_list)
    print('')
    idx2 = False
    for x in range(0, len(lst)):
        if lst[x] == 2:
            idx2 = True
            i = deepcopy(x)
            break
    if not idx2:
        new_list.append(lst)
        print('zwracam1', new_list)
        return new_list
    else:
        print('lst jest podwojny', lst )
        lst_t2 = deepcopy(lst)
        lst_t2[i] = 't'
        lst_tf = deepcopy(lst)
        lst_tf[i] = 'f'
        print(lst_t2, lst_tf)
        new_list.append(t2tfsplit(lst_t2, new_list))
        print('po 1 new_list', new_list)
        new_list.append(t2tfsplit(lst_tf, new_list))
        print('zwracam2', new_list)
        return new_list


def generate_Wm_ground_f12(maxpt, theory, Gm1=False, Gm2=False):

    #Wm_ground < 0 | e(-S*)e(-T) X e(T) e(-S*) | 0>
              #--------Wm_ground-----------#
    # T2' = T2 + Tf


    # Find possible commutator combinations for Wm_ground

    if (theory == 'ccsd'):
        maxexc = 2
    elif (theory == 'cc3'):
        maxexc = 3

    # maxexc = 4
    Wm_ground = []

    for il in range(0, maxpt+1):
        list_s = comb(il, maxexc+1)
        lens = len(list_s)
        print('list s', il, list_s)
        for ik in range(0, maxpt+1):

            list_t = comb(ik, maxexc+2)
            lent = len(list_t)
            # print('list t', ik, list_t)
            # number 3 reflects here geminal amplitude, not T3
            for p1 in range(0, len(list_t)):
                for p2 in range(0, len(list_t[p1])):
                    if list_t[p1][p2] == 3:
                        list_t[p1][p2] = 'tf'
            # for p1 in range(0, len(list_s)):
            #     for p2 in range(0, len(list_s[p1])):
            #         if list_s[p1][p2] == 3:
            #             list_s[p1][p2] = 'sf'

            for jk in range(0, lent):
                for jl in range(0, lens):
                    mbpt = comp_min_pt(list_s[jl], list_t[jk])

                    # print('------teraz sprawdzam----------------------', 't=', list_t[jk], 's=', list_s[jl])
                    # print('ich pt =', mbpt, maxpt)

                    if Gm1 == False and Gm2 == False:
                        if mbpt == maxpt:
                            # print('------teraz sprawdzam----------------------', 't=', list_t[jk], 's=', list_s[jl])
                            # nt and ns are the total excitation rank of
                            # all the operatrors T and S respectively
                            # print('tak mniejsze i licze')
                            nt = particle(list_t[jk])
                            ns = particle(list_s[jl])

                            # print(ns, nt, '---')

                            # Below n is only in (-1, 0, 1) as Y in Wm_ground is
                            # one electron operator
                            for n in range(-1, 2):

                                minidict = {}
                                minidict['T_list'] = list_t[jk]
                                minidict['S_list'] = list_s[jl]
                                minidict['n'] = n
                                # print('')
                                # print('sprawdzam zestaw', n, list_t[jk], nt, list_s[jl], ns)

                                app, sw  = check_conditions_Wm_ground_f12(ik, nt, il, ns, n, 1, list_t[jk], list_s[jl], theory, Gm = 0)

                                #print('check conditions, dodaje?', app, sw, n)
                                # print('')
                                minidict['sw'] = sw
                                minidict['mbpt'] = mbpt

                                if app == True:
                                    print('tak')
                                    Wm_ground.append(minidict)
                    else:

                        if mbpt <= maxpt:
                            # print('------teraz sprawdzam----------------------', 't=', list_t[jk], 's=', list_s[jl])
                            # nt and ns are the total excitation rank of
                            # all the operatrors T and S respectively
                            # print('tak mniejsze i licze')
                            nt = particle(list_t[jk])
                            ns = particle(list_s[jl])

                            # print(ns, nt, '---')

                            # Below n is only in (-1, 0, 1) as Y in Wm_ground is
                            # one electron operator
                            for n in range(-1, 2):

                                minidict = {}
                                minidict['T_list'] = list_t[jk]
                                minidict['S_list'] = list_s[jl]
                                minidict['n'] = n
                                # print('')
                                # print('sprawdzam zestaw', n, list_t[jk], nt, list_s[jl], ns)

                                if Gm1 == True:
                                    app, sw  = check_conditions_Wm_ground_f12(ik, nt, il, ns, n, 1, list_t[jk], list_s[jl], theory, Gm = 1)
                                elif Gm2 == True:
                                    app, sw  = check_conditions_Wm_ground_f12(ik, nt, il, ns, n, 1, list_t[jk], list_s[jl], theory, Gm = 2)

                                # print('check conditions, dodaje?', app, sw)
                                # print('')
                                minidict['sw'] = sw
                                minidict['mbpt'] = mbpt

                                if app == True:
                                    Wm_ground.append(minidict)

    print('Wm_ground')
    for i in Wm_ground:
        print(i)

    return Wm_ground


def generate_Wm_ground_ccd(maxpt, theory, Gm1=False, Gm2=False, cumulative = True):

    # Check conditions for the expression
    # <0|e^(S*)e^(-T) X e^(T) e^(-S*)P(e^(S*)e^(-T) X e^(T) e^(-S*))|0>
    
    # Gm2 = P(e^(S*)e^(-T) X e^(T) e^(-S*))
    # Gm1 = e^(S*)e^(-T) X e^(T) e^(-S*)
    # Gm0 = <0|e^(S*)e^(-T) X e^(T) e^(-S*)|0>
    # See pdf ....

    # Find possible commutator combinations for Wm_ground

    if theory == 'ccd':
        maxexc = 2
    elif theory == 'ccsd':
        maxexc = 2
    elif theory == 'cc3':
        maxexc = 3    

    Wm_ground = []

    for il in range(0, maxpt+1):
        list_s = comb(il, maxexc+1)
        lens = len(list_s)
        # print('list s', il, list_s)
        for ik in range(0, maxpt+1):

            list_t = comb(ik, maxexc+1)
            lent = len(list_t)
            # print('list t', ik, list_t)
            # number 3 reflects here geminal amplitude, not T3

            for jk in range(0, lent):
                for jl in range(0, lens):
                    mbpt = comp_min_pt(list_s[jl], list_t[jk])

                    # print('------teraz sprawdzam----------------------', 't=', list_t[jk], 's=', list_s[jl])
                    # print('ich pt =', mbpt, maxpt)

                    cont = False
                    if cumulative == True:
                        # print('jestem w cumulative true')
                        if mbpt <= maxpt:
                            # print('mbpt', mbpt, maxpt)
                            cont = True
                    elif cumulative == False:
                        # print('jestem w cumulative false')
                        if mbpt == maxpt:
                            # print('mbpt', mbpt, maxpt)
                            cont = True

                    if cont == True:
                    # if mbpt <= maxpt:
                    # if mbpt == maxpt:
                        # print('------teraz sprawdzam----------------------', 't=', list_t[jk], 's=', list_s[jl])
                        # nt and ns are the total excitation rank of
                        # all the operatrors T and S respectively
                        # print('tak mniejsze i licze')
                        nt = particle(list_t[jk])
                        ns = particle(list_s[jl])

                        # print(ns, nt, '---')

                        # Below n is only in (-1, 0, 1) as Y in Wm_ground is
                        # one electron operator
                        for n in range(-1, 2):
                        
                            minidict = {}
                            minidict['T_list'] = list_t[jk]
                            minidict['S_list'] = list_s[jl]
                            minidict['n'] = n
                            # print('')
                            # print('sprawdzam zestaw', n, list_t[jk], nt, list_s[jl], ns)

                            if Gm1 == True:
                                app, sw  = check_conditions_Wm_ground_ccd(ik, nt, il, ns, n, 1, list_t[jk], list_s[jl], theory, Gm = 1)
                            elif Gm2 == True:
                                app, sw  = check_conditions_Wm_ground_ccd(ik, nt, il, ns, n, 1, list_t[jk], list_s[jl], theory, Gm = 2)
                            else:
                                app, sw  = check_conditions_Wm_ground_ccd(ik, nt, il, ns, n, 1, list_t[jk], list_s[jl], theory, Gm = 0)
                                
                            # print('check conditions, dodaje?', app, sw)
                            # print('')
                            minidict['sw'] = sw
                            minidict['mbpt'] = mbpt
                            
                            if app == True:
                                Wm_ground.append(minidict)
                                # print('dodaje ten minidict', minidict)

    # print('Wm_ground')
    # for i in Wm_ground:
    #     print(i)

    return Wm_ground


def create_identificator(d):
    
    ls1 = len(d['Wm1_S_list'])
    lt1 = len(d['Wm1_T_list'])
    ls2 = len(d['Wm2_S_list'])
    lt2 = len(d['Wm2_T_list'])
    ls3 = len(d['Wm3_S_list'])
    
    s = ""
    for i in range(0, ls1):
        s = s + "{x}".format(x = d['Wm1_S_list'][i])

    if ls1 == 0:
        s = 0
        
    t = ""
    for i in range(0, lt1):
        t = t + "{x}".format(x = d['Wm1_T_list'][i])

    if lt1 == 0:
        t = 0

    sm1 = "Wm1t{t}s{s}".format(t=t, s=s)

    sm2 = ""
    sm3 = ""

    return sm1, sm2, sm3

# Wm_overlap =  <P(e(-S)e(T*) \mu_n e(-T*) e(S))| P(e(S*)\mu_l e(-S*))> 
                  #--------Wm2----------------#  #-------Wm3-------#

def generate_overlap(Wm2, Wm3, maxmbpt, cumulative):
    
    Wm_overlap = []

    for a in range(0, len(Wm2)):
        for b in range(0, len(Wm3)):
            sw_bra = Wm2[a]['sw']
            sw_ket = Wm3[b]['sw']
            if (sw_bra == sw_ket):
                n_of_t_and_s = len(Wm2[a]['S_list']) + len(Wm2[a]['T_list']) + len(Wm3[b]['S_list'])
                if (n_of_t_and_s <= 20):
                    mbpt_t = Wm2[a]['mbpt'] + Wm3[b]['mbpt']
                    if cumulative == True:
                        if (mbpt_t <= maxmbpt):
                            minidict = {}
                            minidict['Wm2_S_list'] = Wm2[a]['S_list']
                            minidict['Wm2_T_list'] = Wm2[a]['T_list']
                            minidict['Wm3_S_list'] = Wm3[b]['S_list']
                            minidict['Wm2_n'] = Wm2[a]['n']
                            minidict['Wm3_n'] = Wm3[b]['n']
                            minidict['mbpt'] = mbpt_t
                            Wm_overlap.append(minidict)
                    else:
                        if (mbpt_t == maxmbpt):
                            minidict = {}
                            minidict['Wm2_S_list'] = Wm2[a]['S_list']
                            minidict['Wm2_T_list'] = Wm2[a]['T_list']
                            minidict['Wm3_S_list'] = Wm3[b]['S_list']
                            minidict['Wm2_n'] = Wm2[a]['n']
                            minidict['Wm3_n'] = Wm3[b]['n']
                            minidict['mbpt'] = mbpt_t
                            Wm_overlap.append(minidict)


    return Wm_overlap


# Wg =  <e(-S)e(T*) X e(-T*) e(S) | P(e(S*)\mu_l e(-S*))> 
          #--------Wm1----------#  #-------Wm3-------#

def generate_Wg(Wm1, Wm3, maxmbpt, maxcluster, cumulative):
    # Find all possible Wg

    Wg = []
    Wg_nd = []

    print('Wm1')
    for x in Wm1:
        print(x)
    print('Wm3')
    for x in Wm3:
        print(x)

    for a in range(0, len(Wm1)):
        for c in range(0, len(Wm3)):
            sw_bra = Wm1[a]['sw']
            sw_ket = Wm3[c]['sw']
            
            if (sw_bra == sw_ket):
                mbpt_t = Wm1[a]['mbpt'] + Wm3[c]['mbpt']
                cluster_order = sum(Wm1[a]['T_list']) + sum(Wm1[a]['S_list']) + sum(Wm3[c]['S_list'])
                minidict = {}
                minidict['Wm1_S_list'] = Wm1[a]['S_list']
                minidict['Wm1_T_list'] = Wm1[a]['T_list']
                minidict['Wm3_S_list'] = Wm3[c]['S_list']
                minidict['Wm3_n'] = Wm3[c]['n']
                minidict['mbpt'] = mbpt_t
                if cumulative == True:
                    if (mbpt_t <= maxmbpt):
                        if (cluster_order <= maxcluster):
                            Wm.append(minidict)
                        else:
                            print('niedodane przez maxcluster', cluster_order)
                            Wg_nd.append(minidict)
                else:
                    if (mbpt_t == maxmbpt):
                        if (cluster_order <= maxcluster):
                            print('clustclust', cluster_order, maxcluster)
                            Wg.append(minidict)
                        else:
                            print('niedodane przez maxcluster', cluster_order)
                            print(minidict)
                            Wg_nd.append(minidict)
                    else:
                        print('niedodane przez mbpt', mbpt_t)
                        print(minidict)
                                
            
    for x in Wg:
        print(x)

    print('niedodane')
    for x in Wg_nd:
        print(x)
    print('')

    return Wg


# Gm =  <e(S*)e(-T) E_pq e(T) e(-S*) P(e(S*)e(-T) E_rs e(T) e(-S*))>                                                         
         #--------Gm1----------#  #-------------Gm2-------------#  


def generate_Gm(Gm1, Gm2, maxmbpt, cumulative, method):
    # Find all possible Gm

    Gm = []
    Gm_nd = []
    print('generatte')

    for a in range(0, len(Gm1)):
        for c in range(0, len(Gm2)):
            sw_bra = Gm1[a]['sw']
            sw_ket = Gm2[c]['sw']
            print(sw_bra, sw_ket)
            if (sw_bra == -sw_ket):
                mbpt_t = Gm1[a]['mbpt'] + Gm2[c]['mbpt']
                print('tak, zgadzaja sie poziomy wzbudzenia', mbpt_t)
                minidict = {}
                minidict['Gm1_S_list'] = Gm1[a]['S_list']
                minidict['Gm1_T_list'] = Gm1[a]['T_list']
                minidict['Gm2_S_list'] = Gm2[c]['S_list']
                minidict['Gm2_T_list'] = Gm2[c]['T_list']
                minidict['Gm2_P_list'] = Gm2[c]['sw']
                minidict['Gm2_n'] = Gm2[c]['n']
                minidict['Gm1_n'] = Gm1[a]['n']
                minidict['mbpt'] = mbpt_t
                if cumulative == True:
                    if (mbpt_t <= maxmbpt):
                        # if method = CCD, delete all single excitation opedrators
                        if (method == "ccd"):
                            if 1 not in minidict['Gm1_S_list'] and 1 not in minidict['Gm1_T_list'] and 1 not in minidict['Gm2_S_list'] and 1 not in minidict['Gm2_T_list']:
                                print('dodaje, po usunieciu wzbudzen 1', minidict)
                                Gm.append(minidict)
                        else:
                            print('dodaje', minidict)
                            Gm.append(minidict)

                else:
                    if (mbpt_t == maxmbpt):
                        print('zgadza sie rzad mbpt', minidict)
                        # if method = CCD, delete all single excitation opedrators
                        if (method == "ccd"):
                            if 1 not in minidict['Gm1_S_list'] and 1 not in minidict['Gm1_T_list'] and 1 not in minidict['Gm2_S_list'] and 1 not in minidict['Gm2_T_list']:
                                print('')
                                print('dodaje, po usunieciu wzbudzen 1')
                                print(minidict)
                                print('')
                                Gm.append(minidict)
                        else:
                            print('dodaje', minidict)
                            Gm.append(minidict)

    # print('')
    # print('ostatecznie dostaje takie gm_middle')
    # for x in Gm:
    #     print(x)
    # print('')

    return Gm


                       

# Wm =  <e(-S)e(T*) X e(-T*) e(S) P(e(-S)e(T*) \mu_n e(-T*) e(S))| P(e(S*)\mu_l e(-S*))> 
          #--------Wm1----------#  #-------------Wm2-------------#  #-------Wm3-------#

def generate_Wm(Wm1, Wm2, Wm3, maxmbpt, maxcluster, cumulative):
    # Find all possible Wm

    print('GENERATE_WM', maxmbpt)

    Wm = []
    Wm_nd = []

    # print('Wm1')
    # for x in Wm1:
    #     print(x)
    # print('Wm2')
    # for x in Wm2:
    #     print(x)
    # print('Wm3')
    # for x in Wm3:
    #     print(x)

    jj = 0
    for a in range(0, len(Wm1)):

        for b in range(0, len(Wm2)):
            for c in range(0, len(Wm3)):
                sw_bra = Wm1[a]['sw'] + Wm2[b]['sw'] 
                sw_ket = Wm3[c]['sw']
#                print(Wm1[a]['sw'] , Wm2[b]['sw'], Wm3[c]['sw'])
#                print(Wm1[a] , Wm2[b], Wm3[c])
#                sys.exit(0)
                if (sw_bra == sw_ket):
                    jj += 1
                    print('jjjj', jj)
                    mbpt_t = Wm1[a]['mbpt'] + Wm2[b]['mbpt'] + Wm3[c]['mbpt']
                    cluster_order = sum(Wm1[a]['T_list']) + sum(Wm1[a]['S_list']) + sum(Wm2[b]['T_list']) + \
                        sum(Wm2[b]['S_list']) + sum(Wm3[c]['S_list'])
                    minidict = {}
                    minidict['Wm1_S_list'] = Wm1[a]['S_list']
                    minidict['Wm1_T_list'] = Wm1[a]['T_list']
                    minidict['Wm2_S_list'] = Wm2[b]['S_list']
                    minidict['Wm2_T_list'] = Wm2[b]['T_list']
                    minidict['Wm3_S_list'] = Wm3[c]['S_list']
                    minidict['Wm2_n'] = Wm2[b]['n']
                    minidict['Wm3_n'] = Wm3[c]['n']
                    minidict['mbpt'] = mbpt_t
                    print(minidict)
                    if cumulative == True:
                        if (mbpt_t <= maxmbpt):
                            if (cluster_order <= maxcluster):
                                Wm.append(minidict)
#                                print('dodane1', minidict, sw_bra, sw_ket)
                            else:
                               # print('niedodane przez maxcluster', cluster_order)
                                Wm_nd.append(minidict)
                    else:
                        if (mbpt_t == maxmbpt):
                            if (cluster_order <= maxcluster):
                                #print('clustclust', cluster_order, maxcluster)
 #                               print('dodane2', minidict, sw_bra, sw_ket)
  #                              sys.exit(0)
                                Wm.append(minidict)
                            else:
                             #   print('niedodane przez maxcluster', cluster_order)
                              #  print(minidict)
                                Wm_nd.append(minidict)
#                        else:
                            #print('niedodane przez mbpt', mbpt_t)
                            #print(minidict)
                                
            
    # for x in Wm:
    #     print(x)

    # print('niedodane')
    # for x in Wm_nd:
    #     print(x)
   # print('')
    #sys.exit(0)
    return Wm, Wm_nd


# Um = <  P(e(S*)\mu_l e(-S*)) | e(-S)e(T*) X e(-T*) e(S) P(e(-S)e(T*) \mu_n e(-T*) e(S))> 
          #--------Wm3---------# #--------Wm1-----------# #----------Wm2---------------#


def generate_Um(Um1, Um2, Um3, maxmbpt):
    # Find all possible Um

    Um = []

    for a in range(0, len(Um1)):
        for b in range(0, len(Um2)):
            for c in range(0, len(Um3)):
                sw_ket = Um1[a]['sw'] + Um2[b]['sw'] 
                sw_bra = Um3[c]['sw']

                if (sw_bra == sw_ket):
                    n_of_t_and_s = len(Um1[a]['S_list']) + len(Um1[a]['T_list']) + len(Um2[b]['S_list']) + \
                        len(Um2[b]['T_list']) + len(Um3[c]['S_list'])

                    if (n_of_t_and_s <= 2):

                        mbpt_t = Um1[a]['mbpt'] + Um2[b]['mbpt'] + Um3[c]['mbpt']
                        if (mbpt_t <= maxmbpt):
                            minidict = {}
                            minidict['Um1_S_list'] = Um1[a]['S_list']
                            minidict['Um1_T_list'] = Um1[a]['T_list']
                            minidict['Um2_S_list'] = Um2[b]['S_list']
                            minidict['Um2_T_list'] = Um2[b]['T_list']
                            minidict['Um3_S_list'] = Um3[c]['S_list']
                            minidict['Um2_n'] = Um2[b]['n']
                            minidict['Um3_n'] = Um3[c]['n']
                            minidict['mbpt'] = mbpt_t
                            Um.append(minidict)
    return Um



def generate_WlWr(maxpt, theory):
    # Find all possible commutator combinations of Wl (and Wr)


    Wl = []

    if (theory == 'ccsd'):
        maxexc = 2
        maxmu = 2
    elif (theory == 'cc3'):
        maxexc = 3
        maxmu = 3


    for ik in range(0, maxpt+1):
        list_t = comb(ik, maxexc+1)
        lent = len(list_t)

        for jk in range(0, lent):
            nt = particle(list_t[jk])

            # Below n is only in (-1, 0, 1) as Y in Wm1 is
            # one electron operator
            for n in range(-1, 2):
                app, sw = check_conditions_Wlr(ik, nt, n, 1)

                if (sw >= 1 and sw <= maxmu): 
                    if app == True:
                        mbpt = mbpt_order([list_t[jk]])
                        if (mbpt <=maxpt):
                            minidict = {}
                            minidict['T_list'] = list_t[jk]
                            minidict['sw'] = sw
                            minidict['mbpt'] = mbpt
                            Wl.append(minidict)
    Wr = Wl
    return Wl, Wr

def generate_W_gamma(Wg):
    W_gamma = []

    for m in range(0, len(Wg)):
        minidictm = {}
        minidictm['Wm1_S_list'] = Wg[m]['Wm1_S_list']
        minidictm['Wm1_T_list'] = Wg[m]['Wm1_T_list']
        minidictm['Wm3_S_list'] = Wg[m]['Wm3_S_list']
        minidictm['Wm3_n'] = Wg[m]['Wm3_n']
        minidictm['mbpt'] = Wg[m]['mbpt']
        if (minidictm not in W_gamma):
            W_gamma.append(minidictm)

    return W_gamma


def generate_W_middle_for_Wm1(Wm1):
    W_middle = []

    for m in range(0, len(Wm1)):    
        minidictm = {}
        minidictm['Wm1_S_list'] = Wm1[m]['S_list']
        minidictm['Wm1_T_list'] = Wm1[m]['T_list']
        minidictm['mbpt'] = Wm1[m]['mbpt']
        minidictm['n'] = Wm1[m]['n']
        if (minidictm not in W_middle):
            W_middle.append(minidictm)


    return W_middle


def generate_W_middle_with_explicit_S():

    s12 = arithmetic_string(t1)
    for x in s12:
        x.summation = []
        x.operator_idx = []
        x.operator_type = []

    s13 = evaluate(t1c, t2)
    rint_s13 = s13.integrate(bra = ['a', 'i'], braspin =['s'])
    rsimp_s13 = simplify(rint_s13)
    rsimp_s13.cleanup()

    s21 = arithmetic_string(t2)
    for x in s21:
        x.summation = []
        x.operator_idx = []
        x.operator_type = []
    s23 = evaluate(t2c, t2, t2).scale(0.5)

    rint_s23a = s23.integrate(bra = ['a', 'i', 'b', 'j'], braspin =['s', 's'])
    rint_s23b = s23.integrate(bra = ['a', 'j', 'b', 'i'], braspin =['s', 's'])
    sint_s23  = rint_s23a.scale(1./3.) + rint_s23b.scale(1./6.)
    sint_s23.exec_delta()
    
    ssimp_s23 = simplify(sint_s23)
    ssimp_s23.cleanup()

    for x in ssimp_s23:
        print(x)


    return s21, ssimp_s23, s12, rsimp_s13

def generate_W_middle(Wm):
    W_middle = []

    for m in range(0, len(Wm)):    
        minidictm = {}
        minidictm['Wm1_S_list'] = Wm[m]['Wm1_S_list']
        minidictm['Wm1_T_list'] = Wm[m]['Wm1_T_list']
        minidictm['Wm2_S_list'] = Wm[m]['Wm2_S_list']
        minidictm['Wm2_T_list'] = Wm[m]['Wm2_T_list']
        minidictm['Wm3_S_list'] = Wm[m]['Wm3_S_list']
        minidictm['Wm2_n'] = Wm[m]['Wm2_n']
        minidictm['Wm3_n'] = Wm[m]['Wm3_n']
        minidictm['mbpt'] = Wm[m]['mbpt']
        if (minidictm not in W_middle):
            W_middle.append(minidictm)

    return W_middle

def generate_W_left_W_middle_W_right(Wl, Wm, Wr, maxmbpt, cumulative):
    W_left = []
    W_middle = []
    W_right = []


    for l in range(0, len(Wl)):
        for m in range(0, len(Wm)):
            for r in range(0, len(Wr)):
                mbpt_total = Wl[l]['mbpt'] + Wm[m]['mbpt'] + Wr[r]['mbpt']
                #                print('mbpt_total',  Wm[m]['mbpt'], mbpt_total)

                #if (mbpt_total == maxmbpt):
                if cumulative == True:
                    if  (mbpt_total <= maxmbpt):
                        minidictl = {}
                        minidictm = {}
                        minidictr = {}

                        minidictl['T_list'] = Wl[l]['T_list']
                        minidictl['n'] = Wl[l]['sw']
                        minidictl['mbpt'] = Wl[l]['mbpt']

                        minidictm['Wm1_S_list'] = Wm[m]['Wm1_S_list']
                        minidictm['Wm1_T_list'] = Wm[m]['Wm1_T_list']
                        minidictm['Wm2_S_list'] = Wm[m]['Wm2_S_list']
                        minidictm['Wm2_T_list'] = Wm[m]['Wm2_T_list']
                        minidictm['Wm3_S_list'] = Wm[m]['Wm3_S_list']
                        minidictm['Wm2_n'] = Wm[m]['Wm2_n']
                        minidictm['Wm3_n'] = Wm[m]['Wm3_n']
                        minidictm['mbpt'] = Wm[m]['mbpt']

                        minidictr['T_list'] = Wr[r]['T_list']
                        minidictr['n'] = Wr[r]['sw']
                        minidictr['mbpt'] = Wr[r]['mbpt']

                        if (minidictl not in W_left):
                            W_left.append(minidictl)
                        if (minidictm not in W_middle):
                            W_middle.append(minidictm)
                        if (minidictr not in W_right):
                            W_right.append(minidictr)
                else:
                    if (mbpt_total == maxmbpt):
                        minidictl = {}
                        minidictm = {}
                        minidictr = {}

                        minidictl['T_list'] = Wl[l]['T_list']
                        minidictl['n'] = Wl[l]['sw']
                        minidictl['mbpt'] = Wl[l]['mbpt']

                        minidictm['Wm1_S_list'] = Wm[m]['Wm1_S_list']
                        minidictm['Wm1_T_list'] = Wm[m]['Wm1_T_list']
                        minidictm['Wm2_S_list'] = Wm[m]['Wm2_S_list']
                        minidictm['Wm2_T_list'] = Wm[m]['Wm2_T_list']
                        minidictm['Wm3_S_list'] = Wm[m]['Wm3_S_list']
                        minidictm['Wm2_n'] = Wm[m]['Wm2_n']
                        minidictm['Wm3_n'] = Wm[m]['Wm3_n']
                        minidictm['mbpt'] = Wm[m]['mbpt']

                        minidictr['T_list'] = Wr[r]['T_list']
                        minidictr['n'] = Wr[r]['sw']
                        minidictr['mbpt'] = Wr[r]['mbpt']

                        if (minidictl not in W_left):
                            W_left.append(minidictl)
                        if (minidictm not in W_middle):
                            W_middle.append(minidictm)
                        if (minidictr not in W_right):
                            W_right.append(minidictr)

                        


    return W_left, W_middle, W_right


def generate_U_middle(Wl, Um, Wr, maxmbpt):
    U_middle = []

    for l in range(0, len(Wl)):
        for m in range(0, len(Um)):
            for r in range(0, len(Wr)):
                mbpt_total = Wl[l]['mbpt'] + Um[m]['mbpt'] + Wr[r]['mbpt']

                if (mbpt_total <= maxmbpt):

                    minidictm = {}

                    minidictm['Um1_S_list'] = Um[m]['Um1_S_list']
                    minidictm['Um1_T_list'] = Um[m]['Um1_T_list']
                    minidictm['Um2_S_list'] = Um[m]['Um2_S_list']
                    minidictm['Um2_T_list'] = Um[m]['Um2_T_list']
                    minidictm['Um3_S_list'] = Um[m]['Um3_S_list']
                    minidictm['Um2_n'] = Um[m]['Um2_n']
                    minidictm['Um3_n'] = Um[m]['Um3_n']
                    minidictm['mbpt'] = Um[m]['mbpt']

                    if (minidictm not in U_middle):
                        U_middle.append(minidictm)
    for x in U_middle:
        print(x)
    return U_middle

def print_super_P(n):

    pp = """\\hat{{\\mathcal{{P}}}}_{n}
    """.format(n=n)

    return pp

def print_amp(amp):
    
    if (amp == T_op):
        T = 'T'
    elif (amp == T_opd):
        T = 'T^{\\dagger}'
    elif (amp == S_op):
        T = 'S'
    elif (amp == S_opd):
        T = 'S^{\\dagger}'

    return T

def print_comm(T, lst, op, lr):

    s = op
    if lr == 'l':
        for i in range(0, len(lst)):
            ss = "{T}_{x}".format(T=T, x=lst[i])
            s = print_comm_basic(ss, s)
    elif lr == 'r':
        for i in range(0, len(lst)):
            if lst[i] == 'tf':
                ss = f"{T}^{{'}}_2"
                s = print_comm_basic(s, ss)
            else:
                ss = "{T}_{x}".format(T=T, x=lst[i])
                s = print_comm_basic(s, ss)

    return s
        

def print_comm_basic(op1, op2):
    
    s = "[{op1}, {op2}]".format(op1=op1, op2=op2)

    return s


def compute_num_factor(lst):

    l = len(lst)

    nf = 1
    for i in range(1, len(lst)+1):
        nf = nf * i
    return nf

def latex_W_left(W_left):

    for x in range(0, len(W_left)):

        T = print_amp(T_op)
        op = 'Y'
        s = print_comm(T, W_left[x]['T_list'], op, 'r')

        nf = compute_num_factor(W_left[x]['T_list'])
        n = W_left[x]['n']

        if (nf > 1):
            strg = "\\frac1{nf}\\langle {s} | \\mu_{n}\\rangle + ".format(nf=nf, s=s, n=n)
            print(strg)
        else:
            strg = "\\langle {s} | \\mu_{n}\\rangle + ".format(nf=nf, s=s, n=n)
            print(strg)

def latex_W_right(W_right):

    for x in range(0, len(W_right)):

        T = print_amp(T_op)
        op = 'Y'
        s = print_comm(T, W_right[x]['T_list'], op, 'r')

        nf = compute_num_factor(W_right[x]['T_list'])
        n = W_right[x]['n']

        if (nf > 1):
            strg = "\\frac1{nf}\\langle\\mu_{n}| {s}  \\rangle + ".format(nf=nf, s=s, n=n)
            print(strg)
        else:
            strg = "\\langle  \\mu_{n}|{s} \\rangle + ".format(s=s, n=n)
            print(strg)

def computational_cost_W_middle(W_middle):

    # compute particle of Wm1
    wm1_p1 = particle(W_middle['Wm1_T_list'])
    wm1_p2 = particle(W_middle['Wm1_S_list'])
    wm1_p1l = len(W_middle['Wm1_T_list'])
    wm1_p2l = len(W_middle['Wm1_S_list'])
    wm1_m1 = 1 + wm1_p1 + wm1_p2 - wm1_p1l - wm1_p2l
    
    
    # compute particle of Wm2
    wm2_p1 = particle(W_middle['Wm2_T_list'])
    wm2_p2 = particle(W_middle['Wm2_S_list'])
    wm2_p1l = len(W_middle['Wm2_T_list'])
    wm2_p2l = len(W_middle['Wm2_S_list'])
    wm2_m2 = W_middle['Wm2_n'] + wm2_p1 + wm2_p2 - wm2_p1l - wm2_p2l
    
    # compute particle of Wm3
    wm3_p2 = particle(W_middle['Wm3_S_list'])
    wm3_p2l = len(W_middle['Wm3_S_list'])
    wm3_m3 = W_middle['Wm3_n'] + wm3_p2 - wm3_p2l
    
    # compute excitation of Wm3
    wm3_s3 = W_middle['Wm3_n'] - wm3_p2

    # compute computational cost
    
    sum_idx = (wm1_m1 + wm2_m2 + wm3_m3 - wm3_s3 - 1)*2
    return sum_idx
        

def latex_Wm_overlap(Wm_overlap, cumulative, maxmbpt):

    z = 1
    for x in range(0, len(Wm_overlap)):

        # W_middle =  <P(e(-S)e(T*) \mu_n e(-T*) e(S))| P(e(S*)\mu_l e(S))> 
                       #-------------Wm2-------------#  #-------Wm3-------#

        # Create string for Wm1

        Td = print_amp(T_opd)
        S = print_amp(S_op)
        Sd = print_amp(S_opd)

        op = "\\mu_{n}".format(n = Wm_overlap[x]['Wm2_n'])
        s2_inner = print_comm(Td, Wm_overlap[x]['Wm2_T_list'], op, 'r')
        s2 = print_comm(S,  Wm_overlap[x]['Wm2_S_list'], s2_inner, 'l')

        op = "\\mu_{n}".format(n = Wm_overlap[x]['Wm3_n'])
        s3 = print_comm(Sd,  Wm_overlap[x]['Wm3_S_list'], op, 'r')

        nf2t = compute_num_factor(Wm_overlap[x]['Wm2_T_list'])
        nf2s = compute_num_factor(Wm_overlap[x]['Wm2_S_list'])

        nf3s = compute_num_factor(Wm_overlap[x]['Wm3_S_list'])

        nf = nf2t * nf2s * nf3s

        mbpt = Wm_overlap[x]['mbpt']

        if cumulative == True:
            if (mbpt <= maxmbpt):


                if (nf > 1):
                    strg = "{z}&\\frac1{nf}\\langle {s2} | {s3}\\rangle + \\\\"\
                    .format(z=z, nf=nf, s2=s2, s3=s3)
                    print(strg)
                    z += 1
                else:
                    strg = "{z}& \\langle {s2} | {s3} \\rangle + \\\\".\
                           format(z=z, s2=s2, s3=s3)
                    print(strg)
                    z += 1


def latex_U_middle(U_middle):

    for x in range(0, len(U_middle)):

        # Um = <  P(e(S*)\mu_l e(-S*)) | e(-S)e(T*) X e(-T*) e(S) P(e(-S)e(T*) \mu_n e(-T*) e(S))> 
               #--------Wm3---------# #--------Wm1-----------# #----------Wm2---------------#


        # Create string for Um1

        Td = print_amp(T_opd)
        S = print_amp(S_op)
        Sd = print_amp(S_opd)

        op = 'X'

        s1_inner = print_comm(Td, U_middle[x]['Um1_T_list'], op, 'r')
        s1 = print_comm(S,  U_middle[x]['Um1_S_list'], s1_inner, 'l')


        op = "\\mu_{n}".format(n = U_middle[x]['Um2_n'])
        s2_inner = print_comm(Td, U_middle[x]['Um2_T_list'], op, 'r')
        s2 = print_comm(S,  U_middle[x]['Um2_S_list'], s2_inner, 'l')

        op = "\\mu_{n}".format(n = U_middle[x]['Um3_n'])
        s3 = print_comm(Sd,  U_middle[x]['Um3_S_list'], op, 'r')


        nf1t = compute_num_factor(U_middle[x]['Um1_T_list'])
        nf1s = compute_num_factor(U_middle[x]['Um1_S_list'])

        nf2t = compute_num_factor(U_middle[x]['Um2_T_list'])
        nf2s = compute_num_factor(U_middle[x]['Um2_S_list'])

        nf3s = compute_num_factor(U_middle[x]['Um3_S_list'])

        nf = nf1t * nf1s * nf2t * nf2s * nf3s
        mbpt = U_middle[x]['mbpt']
     
        if (nf > 1):
            strg = "&{x}\\qquad \\frac1{nf}\\langle {s3} |{s1}{s2}\\rangle + \\\\"\
                .format(x=x, nf=nf, s1=s1, s2=s2, s3=s3, mbpt=mbpt)
            print(strg)
        else:
            strg = "&{x}\\qquad \\langle {s3} |{s1}{s2} \\rangle + \\\\".\
                format(x=x, s1=s1, s2=s2, s3=s3, mbpt=mbpt)
            print(strg)

def latex_W_middle(W_middle, n_int, n_op, maxmbpt, cumulative):

    print('LATEX W MIDDLE')

    for x in W_middle:
        print(x)
    
    cost_list = []
    op_list = []
    for x in range(0, len(W_middle)):

        # W_middle =  <e(-S)e(T*) X e(-T*) e(S) P(e(-S)e(T*) \mu_n e(-T*) e(S))| P(e(S*)\mu_l e(S))> 
                       #--------Wm1----------#  #-------------Wm2-------------#  #-------Wm3-------#

        # Create string for Wm1

        Td = print_amp(T_opd)
        S = print_amp(S_op)
        Sd = print_amp(S_opd)

        op = 'X'

        s1_inner = print_comm(Td, W_middle[x]['Wm1_T_list'], op, 'r')
        s1 = print_comm(S,  W_middle[x]['Wm1_S_list'], s1_inner, 'l')


        op = "\\mu_{n}".format(n = W_middle[x]['Wm2_n'])
        s2_inner = print_comm(Td, W_middle[x]['Wm2_T_list'], op, 'r')
        s2 = print_comm(S,  W_middle[x]['Wm2_S_list'], s2_inner, 'l')

        op = "\\mu_{n}".format(n = W_middle[x]['Wm3_n'])
        s3 = print_comm(Sd,  W_middle[x]['Wm3_S_list'], op, 'r')


        nf1t = compute_num_factor(W_middle[x]['Wm1_T_list'])
        nf1s = compute_num_factor(W_middle[x]['Wm1_S_list'])

        nf2t = compute_num_factor(W_middle[x]['Wm2_T_list'])
        nf2s = compute_num_factor(W_middle[x]['Wm2_S_list'])

        nf3s = compute_num_factor(W_middle[x]['Wm3_S_list'])

        nf = nf1t * nf1s * nf2t * nf2s * nf3s
        mbpt = W_middle[x]['mbpt']
     
        #sum_idx = computational_cost_W_middle(W_middle[x])
        #cost_list.append(sum_idx)
        #n_of_int = n_int[x]
        #n_of_op  = n_op[x]
        #op_list.append(n_of_op)
        #        print('nfff', nf, nf1t , nf1s , nf2t , nf2s , nf3s)
        cluster_op = len(W_middle[x]['Wm1_T_list']) + len(W_middle[x]['Wm1_S_list']) + len(W_middle[x]['Wm2_T_list']) + \
                     len(W_middle[x]['Wm2_S_list']) + len(W_middle[x]['Wm3_S_list'])
        if cumulative == True:
            #if (n_of_op <= maxop and mbpt <= maxmbpt):
            if (mbpt <= maxmbpt):
                if (nf > 1):
                    #strg = "{x}& $\\frac1{nf}\\langle {s1}{s2} | {s3}\\rangle $& ({mbpt})& $N^{sum_idx} $& {n_of_int} calek & {n_of_op} \\\\".format(x=x, nf=nf, s1=s1, s2=s2, s3=s3, mbpt=mbpt, sum_idx=sum_idx, n_of_int=n_of_int, n_of_op=n_of_op)
                    if (x%3 == 0):
                        strg = "&+ \\frac1{nf}\\langle {s1}{s2} | {s3}\\rangle ".format(x=x, nf=nf, s1=s1, s2=s2, s3=s3)
                    elif(x%3 == 2):
                        strg = "+\\frac1{nf}\\langle {s1}{s2} | {s3}\\rangle \\\\".format(x=x, nf=nf, s1=s1, s2=s2, s3=s3)
                    else:
                        strg = "+\\frac1{nf}\\langle {s1}{s2} | {s3}\\rangle".format(x=x, nf=nf, s1=s1, s2=s2, s3=s3)
                    print(strg)
                else:
                    #strg = "{x}& $\\langle {s1}{s2} | {s3} \\rangle$ &({mbpt})& $N^{sum_idx}$ & {n_of_int} calek& {n_of_op}  \\\\".format(x=x, s1=s1, s2=s2, s3=s3, mbpt=mbpt, sum_idx=sum_idx, n_of_int=n_of_int, n_of_op=n_of_op)
                    if (x%3 == 0):
                        strg = "&+ \\langle {s1}{s2} | {s3} \\rangle ".format(x=x, s1=s1, s2=s2, s3=s3)
                    elif (x%3 ==2):
                        strg = "+\\langle {s1}{s2} | {s3} \\rangle \\\\".format(x=x, s1=s1, s2=s2, s3=s3)
                    else:
                        strg = "+\\langle {s1}{s2} | {s3} \\rangle".format(x=x, s1=s1, s2=s2, s3=s3)
                    print(strg)
        else:
            #if (n_of_op <= maxop and mbpt == maxmbpt):
            if (mbpt == maxmbpt):
                if (nf > 1):
                    #strg = "{x}& $\\frac1{nf}\\langle {s1}{s2} | {s3}\\rangle $& ({mbpt})& $N^{sum_idx} $& {n_of_int} calek & {n_of_op} \\\\".format(x=x, nf=nf, s1=s1, s2=s2, s3=s3, mbpt=mbpt, sum_idx=sum_idx, n_of_int=n_of_int, n_of_op=n_of_op)
                    if (x%3 == 0):
                        strg = "{x}&+ \\frac1{nf}\\langle {s1}{s2} | {s3}\\rangle ".format(x=x, nf=nf, s1=s1, s2=s2, s3=s3)
                    elif(x%3 == 2):
                        strg = "{x}+\\frac1{nf}\\langle {s1}{s2} | {s3}\\rangle \\\\".format(x=x, nf=nf, s1=s1, s2=s2, s3=s3)
                    else:
                        strg = "{x}+\\frac1{nf}\\langle {s1}{s2} | {s3}\\rangle ".format(x=x, nf=nf, s1=s1, s2=s2, s3=s3)
                    print(strg)
                else:
                    #strg = "{x}& $\\langle {s1}{s2} | {s3} \\rangle$ &({mbpt})& $N^{sum_idx}$ & {n_of_int} calek& {n_of_op}  \\\\".format(x=x, s1=s1, s2=s2, s3=s3, mbpt=mbpt, sum_idx=sum_idx, n_of_int=n_of_int, n_of_op=n_of_op)
                    if (x%3 == 0):
                        strg = "{x}&+ \\langle {s1}{s2} | {s3} \\rangle ".format(x=x, s1=s1, s2=s2, s3=s3)
                    elif (x%3 ==2):
                        strg = "{x}+\\langle {s1}{s2} | {s3} \\rangle \\\\".format(x=x, s1=s1, s2=s2, s3=s3)
                    else:
                        strg = "{x}+\\langle {s1}{s2} | {s3} \\rangle ".format(x=x, s1=s1, s2=s2, s3=s3)
                    print(strg)



def latex_Gm(Gm1, Gm2, maxmbpt, f):

    print('LATEX Gm ')
    
    cost_list = []
    op_list = []
    latex_comm = ""
    kk = -1
    for x1 in range(0, len(Gm1)):
        for x2 in range(0, len(Gm2)):
            # Gm =  <e(S*)e(-T) E_pq e(T) e(-S*) P(e(S*)e(-T) E_rs e(T) e(-S*))>                                                         
            #--------Gm1----------#  #-------------Gm2-------------#  

            # Create string for Wm1

            Td = print_amp(T_opd)
            T = print_amp(T_op)
            S = print_amp(S_op)
            Sd = print_amp(S_opd)
            p = print_super_P(Gm2[x2]['sw'])

            opx = 'E_{pq}'
            opy = 'E_{rs}'

            s2_inner = print_comm(T, Gm2[x2]['T_list'], opy, 'r')
            s2 = print_comm(Sd, Gm2[x2]['S_list'], s2_inner, 'l')

            s2_p = f"{p}\\left({s2}\\right)"


            s1_inner = print_comm(T, Gm1[x1]['T_list'], opx, 'r')
            s1 = print_comm(Sd, Gm1[x1]['S_list'], s1_inner, 'l')

            nf1t = compute_num_factor(Gm1[x1]['T_list'])
            nf1s = compute_num_factor(Gm1[x1]['S_list'])

            nf2t = compute_num_factor(Gm2[x2]['T_list'])
            nf2s = compute_num_factor(Gm2[x2]['S_list'])


            nf = nf1t * nf1s * nf2t * nf2s 
            mbpt = Gm1[x1]['mbpt'] + Gm2[x2]['mbpt']
            if Gm1[x1]['sw'] == - Gm2[x2]['sw']:
                if (mbpt == maxmbpt):
                    kk += 1
                    if kk == 0:
                        latex_comm = "\equa{\n"
                    if (nf > 1):
                        if (kk%2 == 0):
                            strg = f"&+ \\frac1{nf}\\langle {s1}{s2_p} \\rangle"
                        elif(kk%2 == 1):
                            strg = f"+\\frac1{nf}\\langle {s1}{s2_p} \\rangle \\\\ "
                        else:
                            strg = f"+\\frac1{nf}\\langle {s1}{s2_p} \\rangle "
                        print(strg)
                        latex_comm += strg
                    else:
                        if (kk%2 == 0):
                            strg = f"&+ \\langle {s1}{s2_p}  \\rangle "
                        elif (kk%2 ==1):
                            strg = f"+\\langle {s1}{s2_p}  \\rangle \\\\ "
                        else:
                            strg = f"+\\langle {s1}{s2_p}  \\rangle "
                        print(strg)
                        latex_comm += strg
                if (kk%30) == 0:
                    if kk != 0:
                        latex_comm += """
                        }
                        \\newpage
                        \equag
                        {
                        """
    if latex_comm != "":
        latex_comm += "\n}"
        f.write(latex_comm)

def latex_WmWm(Gm1, Gm2, idx1, idx2, maxmbpt, f):

    print('LATEX Wm WM ')
    
    cost_list = []
    op_list = []
    latex_comm = ""
    kk = -1
    for x1 in range(0, len(Gm1)):
        for x2 in range(0, len(Gm2)):
            #    <e(S*)e(-T) E_pq e(T) e(-S*)>    <e(S*)e(-T) E_rs e(T) e(-S*)>                                                         
            #--------Gm1----------#               #-------------Gm2-------------#  

            # Create string for Wm1

            Td = print_amp(T_opd)
            T = print_amp(T_op)
            S = print_amp(S_op)
            Sd = print_amp(S_opd)

            op1 = 'E_{{{p}{q}}}'.format(p=idx1[0], q=idx1[1])
            op2 = 'E_{{{p}{q}}}'.format(p=idx2[0], q=idx2[1])

            s2_inner = print_comm(T, Gm2[x2]['Wm1_T_list'], op2, 'r')
            s2 = print_comm(Sd, Gm2[x2]['Wm1_S_list'], s2_inner, 'l')


            s1_inner = print_comm(T, Gm1[x1]['Wm1_T_list'], op1, 'r')
            s1 = print_comm(Sd, Gm1[x1]['Wm1_S_list'], s1_inner, 'l')

            nf1t = compute_num_factor(Gm1[x1]['Wm1_T_list'])
            nf1s = compute_num_factor(Gm1[x1]['Wm1_S_list'])

            nf2t = compute_num_factor(Gm2[x2]['Wm1_T_list'])
            nf2s = compute_num_factor(Gm2[x2]['Wm1_S_list'])

            
            nf_com = 0.5 # bo wyraz ma być mnożony przez 1/2
            
            nf = 1./nf1t * nf1s * nf2t * nf2s *nf_com
            mbpt = Gm1[x1]['mbpt'] + Gm2[x2]['mbpt']
            kk += 1
            if kk == 0:
                latex_comm = "\equa{\n"
            if (1== 1):
                if (kk%2 == 0):
                    strg = f"&+ {nf}\\langle {s1}\\rangle\\langle{s2} \\rangle"
                elif(kk%2 == 1):
                    strg = f"+{nf}\\langle {s1}\\rangle\\langle{s2} \\rangle \\\\ "
                else:
                    strg = f"+{nf}\\langle {s1}\\rangle\\langle{s2} \\rangle "
                print(strg)
                latex_comm += strg
            # else:
            #     if (kk%2 == 0):
            #         strg = f"&+ \\langle {s1}\\rangle\\langle {s2}  \\rangle "
            #     elif (kk%2 ==1):
            #         strg = f"+\\langle {s1}\\rangle\\langle {s2}  \\rangle \\\\ "
            #     else:
            #         strg = f"+\\langle {s1}\\rangle\\langle {s2}  \\rangle "
            #     print(strg)
            #     latex_comm += strg
        if (kk%30) == 0:
            if kk != 0:
                latex_comm += """
                }
                \\newpage
                \equag
                {
                """
    if latex_comm != "":
        latex_comm += "\n}"
        f.write(latex_comm)

        
def latex_WmWm_delta(Gm2, idx1, idx2, maxmbpt, f):

    print('LATEX WmWmdelta ')
    
    cost_list = []
    op_list = []
    latex_comm = ""
    kk = -1

    for x2 in range(0, len(Gm2)):
        print('gowno123', x2)
        # 1/2   delta_pq                          <e(S*)e(-T) E_rs e(T) e(-S*)>                                                         
        #--------Gm1----------#               #-------------Gm2-------------#  

        # Create string for Wm1

        Td = print_amp(T_opd)
        T = print_amp(T_op)
        S = print_amp(S_op)
        Sd = print_amp(S_opd)

        s1 = '\\delta_{{{p}{q}}}'.format(p=idx1[0], q=idx1[1])
        op2 = 'E_{{{p}{q}}}'.format(p=idx2[0], q=idx2[1])

        s2_inner = print_comm(T, Gm2[x2]['Wm1_T_list'], op2, 'r')
        s2 = print_comm(Sd, Gm2[x2]['Wm1_S_list'], s2_inner, 'l')

        nf1s = -1.0

        nf2t = compute_num_factor(Gm2[x2]['Wm1_T_list'])
        nf2s = compute_num_factor(Gm2[x2]['Wm1_S_list'])


        nf = 1./nf1s * nf2t * nf2s
        mbpt = Gm2[x2]['mbpt']
        kk += 1
        if kk == 0:
            latex_comm = "\equa{\n"
        if (1==1):
            if (kk%2 == 0):
                strg = f"&+ {nf}{s1}\\langle {s2} \\rangle"
            elif(kk%2 == 1):
                strg = f"+{nf}{s1}\\langle {s2} \\rangle \\\\ "
            else:
                strg = f"+{nf}{s1}\\langle {s2} \\rangle "
            print(strg)
            latex_comm += strg
        # else:
        #     if (kk%2 == 0):
        #         strg = f"&+ {s1}\\langle{s2}  \\rangle "
        #     elif (kk%2 ==1):
        #         strg = f"+{s1}\\langle{s2}  \\rangle \\\\ "
        #     else:
        #         strg = f"+{s1}\\langle{s2}  \\rangle "
        #     print(strg)
        #     latex_comm += strg
    if (kk%30) == 0:
        if kk != 0:
            latex_comm += """
            }
            \\newpage
            \equag
            {
            """
    if latex_comm != "":
        latex_comm += "\n}"
        f.write(latex_comm)

def find_explicit_op_form(n1, n2, withCABS = False):

    idx_gm1 = []
    idx_gm2 = []
    banned = []

    if withCABS:
        vblock = completev
    else:
        vblock = virtual

    if n1 == -1:
        a1 = free_idx(occupied, [])
        a2 = free_idx(vblock, [])                    
        banned.append(a1)
        banned.append(a2)
        idxl = [a1, a2]
        idx_gm1.append(idxl)
        # if withCABS:
        #     a3 = free_idx(CABS, [])
        #     banned.append(a3)
        #     idx_gm1.append([a1, a3])

    elif n1 == 1:
        a1 = free_idx(vblock, [])
        a2 = free_idx(occupied, [])
        banned.append(a1)
        banned.append(a2)
        idxl = [a1, a2]
        idx_gm1.append(idxl)
        # if withCABS:
        #     a3 = free_idx(CABS, [])
        #     banned.append(a3)
        #     idx_gm1.append([a3, a2])

    elif n1 == 0:

        a1 = free_idx(vblock, [])
        a2 = free_idx(vblock, [a1])
        banned.append(a1)
        banned.append(a2)
        idxl = [a1, a2]
        # if withCABS:
        #     a3 = free_idx(CABS, [])
        #     banned.append(a3)
        #     a4 = free_idx(CABS, banned)
        #     banned.append(a3)
        #     idx_gm1.append([a1, a3])
        #     idx_gm1.append([a3, a1])
        #     idx_gm1.append([a3, a4])
        

        idx_gm1.append(idxl)
        a3 = free_idx(occupied, [])
        a4 = free_idx(occupied, [a3])
        banned.append(a3)
        banned.append(a4)
        idxl = [a3, a4]

        idx_gm1.append(idxl)

    if n2 == -1:

        a1 = free_idx(occupied, banned)
        a2 = free_idx(vblock, banned)
        idxl = [a1, a2]    
        idx_gm2.append(idxl)
        # if withCABS:
        #     a3 = free_idx(CABS, banned)
        #     banned.append(a3)
        #     idx_gm2.append([a1, a3])

    elif n2 == 1:
        a1 = free_idx(vblock, banned)
        a2 = free_idx(occupied, banned)
        idxl = [a1, a2]
        idx_gm2.append(idxl)
        # if withCABS:
        #     a3 = free_idx(CABS, banned)
        #     banned.append(a3)
        #     idx_gm2.append([a3, a1])

    elif n2 == 0:
        a1 = free_idx(vblock, banned)
        banned.append(a1)
        a2 = free_idx(vblock, banned)
        idxl = [a1, a2]
        idx_gm2.append(idxl)

        # if withCABS:
        #     a3 = free_idx(CABS, banned)            
        #     banned.append(a3)
        #     a4 = free_idx(CABS, banned)
        #     idx_gm2.append([a1, a3])
        #     idx_gm2.append([a3, a1])
        #     idx_gm2.append([a3, a4])
        
        a3 = free_idx(occupied, banned)
        banned.append(a3)
        a4 = free_idx(occupied, banned)
        idxl = [a3, a4]
        idx_gm2.append(idxl)
    # print(idx_gm2)
    # print('~~~~~~~~~~~~~~~~~~~~~')
    return idx_gm1, idx_gm2

def latex_Gmiddle(Gm, maxmbpt, f, withCABS=False):

    print('LATEX Gm ')
    
    cost_list = []
    op_list = []
    latex_comm = ""
    kk = -1
    for x1 in range(0, len(Gm)):
        # Gm =  <e(S*)e(-T) E_pq e(T) e(-S*) P(e(S*)e(-T) E_rs e(T) e(-S*))>                                                         
                #--------Gm1----------#  #-------------Gm2-------------#  

        # Create string for Wm1

        Td = print_amp(T_opd)
        T = print_amp(T_op)
        S = print_amp(S_op)
        Sd = print_amp(S_opd)
        p = print_super_P(Gm[x1]['Gm2_P_list'])
#        print('lalala')
        idx_gm1, idx_gm2 = find_explicit_op_form(Gm[x1]['Gm1_n'], Gm[x1]['Gm2_n'], withCABS)
 #       print(idx_gm1)
  #      print(idx_gm2)
#        idx_gm1 = [['κ', 'λ']] - to bylo tylko do ladnego printowania do artykulu
 #       idx_gm2 = [['ζ', 'τ']]
        # print('gg', Gm[x1])
        # print('n1n2', Gm[x1]['Gm1_n'], Gm[x1]['Gm2_n'], idx_gm1, idx_gm2)
        # print(Gm[x1])

        for i in range(0, len(idx_gm1)):
            for j in range(0, len(idx_gm2)):
                
                op1 = 'E_{{{p}{q}}}'.format(p=idx_gm1[i][0], q = idx_gm1[i][1])
                op2 = 'E_{{{p}{q}}}'.format(p=idx_gm2[j][0], q = idx_gm2[j][1])

                s2_inner = print_comm(T, Gm[x1]['Gm2_T_list'], op2, 'r')
                s2 = print_comm(Sd, Gm[x1]['Gm2_S_list'], s2_inner, 'l')

                s2_p = f"{p}\\left({s2}\\right)"


                s1_inner = print_comm(T, Gm[x1]['Gm1_T_list'], op1, 'r')
                s1 = print_comm(Sd, Gm[x1]['Gm1_S_list'], s1_inner, 'l')

                nf1t = compute_num_factor(Gm[x1]['Gm1_T_list'])
                nf1s = compute_num_factor(Gm[x1]['Gm1_S_list'])

                nf2t = compute_num_factor(Gm[x1]['Gm2_T_list'])
                nf2s = compute_num_factor(Gm[x1]['Gm2_S_list'])


                nf = nf1t * nf1s * nf2t * nf2s 
                # mbpt = Gm[x1]['Gm1_mbpt'] + Gm[x]['mbpt']
                kk += 1
                if kk == 0:
                    latex_comm = "\equa{\n"
                if (nf > 1):
                    if (kk%2 == 0):
                        strg = f"&+ \\frac1{nf}\\langle {s1}{s2_p} \\rangle"
                    elif(kk%2 == 1):
                        strg = f"+\\frac1{nf}\\langle {s1}{s2_p} \\rangle \\\\ "
                    else:
                        strg = f"+\\frac1{nf}\\langle {s1}{s2_p} \\rangle "
                    print(strg)
                    latex_comm += strg
                else:
                    if (kk%2 == 0):
                        strg = f"&+ \\langle {s1}{s2_p}  \\rangle "
                    elif (kk%2 ==1):
                        strg = f"+\\langle {s1}{s2_p}  \\rangle \\\\ "
                    else:
                        strg = f"+\\langle {s1}{s2_p}  \\rangle "
                    print(strg)
                    latex_comm += strg
                if (kk%30) == 0:
                    if kk != 0:
                        latex_comm += """
                        }
                        \\newpage
                        \equag
                        {
                        """
    if latex_comm != "":
        latex_comm += "\n}"
        f.write(latex_comm)

    
    
def latex_W_middle_Wm1(W_middle2, mbpt_in, cumulative, f, start=0):

    print('LATEX W MIDDLE')
    print('this is input')
    for x in W_middle2:
        print(x)
    print('')
    maxmbpt = 10
    W_middle = []
    # for x in W_middle2:
    #     # print(x['Wm1_T_list'])
    #     if 1 not in x['Wm1_S_list']:
    #         if 1 not in x['Wm1_T_list']:
    #             W_middle.append(x)
    #     print(x)
    W_middle = W_middle2

    cost_list = []
    op_list = []
    latex_comm = "\equa{\n"

    for x in range(0, len(W_middle)):

        # W_middle_Wm1 =  <e(S*)e(-T) X e(T) e(-S*)> 
                       #--------Wm1----------#  

        # Create string for Wm1

        T = print_amp(T_op)
        S = print_amp(S_op)
        Sd = print_amp(S_opd)

        op = 'X'
        # print('rykrykryk', W_middle[x]['Wm1_T_list'])
        s1_inner = print_comm(T, W_middle[x]['Wm1_T_list'], op, 'r')
        s1 = print_comm(Sd,  W_middle[x]['Wm1_S_list'], s1_inner, 'l')

        nf1t = compute_num_factor(W_middle[x]['Wm1_T_list'])
        nf1s = compute_num_factor(W_middle[x]['Wm1_S_list'])

        nf = nf1t * nf1s
        mbpt = W_middle[x]['mbpt']
     

        cluster_op = len(W_middle[x]['Wm1_T_list']) + len(W_middle[x]['Wm1_S_list'])

        # print('cumulative', cumulative)
        if cumulative == True:
            if (mbpt <= maxmbpt):
                if x == len(W_middle)-1:
                    strg = "&+ \\langle 0|{s1} | 0\\rangle ".format(x=x, s1=s1)
                else:
                    strg = "&+ \\langle 0|{s1} | 0\\rangle\\\\\n ".format(x=x, s1=s1)
                latex_comm += strg
                print(strg)
        else:
            # print('tak1', mbpt_in)
            if (mbpt == mbpt_in):
                # print('tak2')
                if x == len(W_middle)-1:
                    # print('tak3')
                    strg = "{x}&+ \\langle 0|{s1} | 0\\rangle ".format(x=x, s1=s1)
                else:
                    # print('tak4')
                    strg = "{x}&+ \\langle 0|{s1} | 0\\rangle\\\\\n ".format(x=x, s1=s1)
                latex_comm += strg
                print(strg)

    latex_comm += "}\n"
    if start == 0:
        f.write(latex_comm)



#<e(S*)e(-T) Y e(T) e(-S*) P(e(S*)\mu_l e(S))>    <mu_p|e(-T)Xe(T)>
#----------------gamma-----------------------#    #------xi-------#

# gamma = <e(S*)e(-T) Y e(T) e(-S*) P(e(S*)\mu_l e(-S*))>
          #--------gamma1---------#  #-----gamma2-------#


def generate_gamma1():

    # Find all possible commutator combinations for gamma1
    # for given mbpt order


    k = 3         # - nesting of T
    l = 3         # - nesting of S
    maxexc = 3    # - maximum excitation of T or S
    mbpt_max = 4  # - maximum mbpt order
    gamma1 = []

    for ik in range(0, k+1):

        list_t = comb(ik, maxexc+1)
        lent = len(list_t)

        for il in range(0, l+1):

            list_s = comb(il, maxexc+1)
            lens = len(list_s)

            for jk in range(0, lent):
                for jl in range(0, lens):

                    # nt and ns are the total excitation rank of
                    # all the operatrors T and S respectively

                    nt = particle(list_t[jk])
                    ns = particle(list_s[jl])

                    for n in range(-1, 2):
                        app, sw  = check_conditions_gamma1(ik, nt, il, ns, n, 1)

                        if app == True:
                            mbpt = mbpt_order([list_t[jk], list_s[jl]])
                            if (mbpt <= mbpt_max):
                                minidict = {}
                                minidict['T_list'] = list_t[jk]
                                minidict['S_list'] = list_s[jl]
                                minidict['n'] = n
                                minidict['sw'] = sw
                                minidict['mbpt'] = mbpt
                                gamma1.append(minidict)
    return gamma1

def generate_gamma2():
    # Find possible commutator combinations for gamma2
    gamma2 = generate_Wm3()
    return gamma2

def generate_xi():

    # Find possible commutator combinations for xi
    xil, xir = generate_WlWr()

    return xir

def generate_gamma():
    # Find all possible gamma

    gamma = []
    for a in range(0, len(gamma1)):
        for b in range(0, len(gamma2)):
                sw_g1 = gamma1[a]['sw'] 
                sw_g2 = gamma2[b]['sw']

                if (sw_g1 + sw_g2 == 0):
                    # print(gamma1[a]['T_list'], gamma1[a]['S_list'], sw_g1)
                    # print(gamma2[b]['S_list'], sw_g2)
                    # print()
                    mbpt_t = gamma1[a]['mbpt'] + gamma2[b]['mbpt']
                    if (mbpt_t <= 3):
                        minidict = {}
                        minidict['g1_S_list'] = gamma1[a]['S_list']
                        minidict['g1_T_list'] = gamma1[a]['T_list']
                        minidict['g2_S_list'] = gamma2[b]['S_list']
                        minidict['g2_n'] = gamma2[b]['n']
                        minidict['mbpt'] = mbpt_t
                        gamma.append(minidict)
    return gamma
def generate_W_xi_W_gamma(gamma, xi):
    W_xi = []
    W_gamma = []

    for i in range(0, len(gamma)):
        for j in range(0, len(xi)):

            mbpt_total = gamma[i]['mbpt'] + xi[j]['mbpt']

            if (mbpt_total <= 3):

                minidictxi = {}
                minidictga = {}

                minidictxi['T_list'] = xi[j]['T_list']
                minidictxi['n'] = xi[j]['sw']
                minidictxi['mbpt'] = xi[j]['mbpt']

                minidictga['g1_S_list'] = gamma[i]['g1_S_list']
                minidictga['g1_T_list'] = gamma[i]['g1_T_list']
                minidictga['g2_S_list'] = gamma[i]['g2_S_list']
                minidictga['g2_n'] = gamma[i]['g2_n']
                minidictga['mbpt'] = gamma[i]['mbpt']

                if (minidictxi not in W_xi):
                    W_xi.append(minidictxi)
                if (minidictga not in W_gamma):
                    W_gamma.append(minidictga)

    return W_xi, W_gamma

def latex_gamma(W_gamma):

    for x in range(0, len(W_gamma)):

        # W_gamma =  <e(S*)e(-T) X e(T) e(-S*) P(e(S*)\mu_l e(-S*))>
                     #--------gamma1--------#  #-----gamma2-------#
        # Create string for gamma1

        T = print_amp(T_op)
        Sd = print_amp(S_opd)

        op = 'Y'
        # print(W_middle[x]['Wm1_T_list'])

        s1_inner = print_comm(T, W_gamma[x]['g1_T_list'], op, 'r')
        s1 = print_comm(Sd,  W_gamma[x]['g1_S_list'], s1_inner, 'l')


        op = "\\mu_{n}".format(n = W_gamma[x]['g2_n'])
        s2 = print_comm(Sd,  W_gamma[x]['g2_S_list'], op, 'l')


        nf1t = compute_num_factor(W_gamma[x]['g1_T_list'])
        nf1s = compute_num_factor(W_gamma[x]['g1_S_list'])

        nf2s = compute_num_factor(W_gamma[x]['g2_S_list'])

        nf = nf1t * nf1s * nf2s
        mbpt = W_gamma[x]['mbpt']



        if (nf > 1):
            strg = "&\\frac1{nf}\\langle {s1}{s2}\\rangle + \\qquad({mbpt})\\\\".format(nf=nf, s1=s1, s2=s2, mbpt=mbpt)
            print(strg)
        else:
            strg = "&\\langle {s1}{s2} \\rangle + \\qquad({mbpt})\\\\".format(s1=s1, s2=s2, mbpt=mbpt)
            print(strg)

def latex_W_xi(W_xi):    

    for x in range(0, len(W_xi)):

        T = print_amp(T_op)
        op = 'Y'
        s = print_comm(T, W_xi[x]['T_list'], op, 'r')

        nf = compute_num_factor(W_xi[x]['T_list'])
        n = W_xi[x]['n']

        if (nf > 1):
            strg = "\\frac1{nf}\\langle\\mu_{n}| {s}  \\rangle + ".format(nf=nf, s=s, n=n)
            print(strg)
        else:
            strg = "\\langle  \\mu_{n}|{s} \\rangle + ".format(s=s, n=n)
            print(strg)




