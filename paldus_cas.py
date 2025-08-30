# ------------------------------------------------------
# Authors: Aleksandra Tucholska, University of Warsaw
# ------------------------------------------------------
from params import * 
from copy import deepcopy
from itertools import permutations
from itertools import product
from itertools import combinations
from paldus_classes import ugg
from paldus_classes import arithmetic_string
from paldus_classes import nswaps_in_two_lists
import collections
import math
import sys
from wick import integrate_wick_basic
from fractions import Fraction
# import numba
# from numba import jit
import re

class cas:
    # class of creators and anihilators strings

    def __init__(self):
        self.summation = []
        self.coefficient = []
        self.coefficient_idx = []
        self.coefficient_spin = []
        self.operator_idx = []
        self.operator_type = []
        self.num_factor = 1.0
        self.delta = []
        self.delta_spin = []
        self.coefficient_rep = []

    def __str__(self):

        if self is None:
            return str(0)

        summ = ""
        opidx = ""
        cidx = ""
        temp = ""
        tempo = ""
        tempc = ""
        delidx = ""
        permidx = ""
        sumidx = ""

        if self.summation != []:
            for x in self.summation:
                summ = summ + str(x) + str("")
            sumidx = "\sum_{" + summ + "}"
        upperi = ""
        loweri = ""
        if len(self.operator_idx) != 0:
            for x in range(0, len(self.operator_idx)):
                temp = self.operator_idx[x]
                if self.operator_type[x] == CRE:
                    opidx = opidx + temp+"^+"
                elif self.operator_type[x] == ANI:
                    opidx = opidx + temp

        if len(self.coefficient) != 0:
            for x in range(0, len(self.coefficient)):
                temp = ""
                for y in range(0, len(self.coefficient_idx[x])):
                    temp = temp + str(self.coefficient_idx[x][y]) + str("")

                if len(self.coefficient_spin) != 0:
                    temp2 = ""
                    for y in range(0, len(self.coefficient_spin[x])):
                        temp2 = temp2 + str(self.coefficient_spin[x][y]) + str("")

                coef = self.coefficient[x]
                if self.coefficient[x] == TWOEL_INT_DIRAC:
                    if len(self.coefficient_spin) != 0:
                        minitemp = ""
                        for y in range(0, len(self.coefficient_idx[x])):
                            minitemp = minitemp + str(self.coefficient_idx[x][y]) + str("")
                            minitemp = minitemp + str(self.coefficient_spin[x][y]) + str("")
                        cidx = cidx + str(coef)+ "<"+minitemp[0:4]+"|"+minitemp[4:8]+">"
                    else:
                            
                        cidx = cidx + str(coef)+ "<"+temp[0:2]+"|"+temp[2:4]+">"
                # elif self.coefficient[x] == DENS3 or self.coefficient[x] == DENS4:
                #     strup = ""
                #     strdown=""
                #     for elem in range(0, len(temp),2):
                #         strup = strup + temp[elem]
                #         strdown = strdown + temp[elem+1]
                #     cidx = cidx + str(coef)+"_{"+strdown+"}^{"+strup+"}"
                        
                else:
                    if len(self.coefficient_spin) != 0:
                        cidx = cidx + str(coef)+ "_{"+temp+"}"+ "^{"+temp2+"}"
                    else:
                        cidx = cidx + str(coef)+ "_{"+temp+"}"
            
        for x in range(0, len(self.delta)):
            
            temp = ""
            for y in self.delta[x]:
                temp = temp + str(y) + str("")
            if len(self.coefficient_spin) != 0:
                
                temp2 = ""
                for y in self.delta_spin[x]:
                    temp2 = temp2 + str(y) + str("")                
            if len(self.coefficient_spin) != 0:
                delidx = delidx +"\delta_{"  + temp+"}"+ "^{"+temp2+"}"
            else:
                delidx = delidx + "\delta_{" + temp + "}"

        if len(self.coefficient) == 0:
            if len(self.operator_idx) == 0:
                if len(self.delta) == 0:
                    if len(self.summation) != 0:
                        if self.summation[0] in occupied:
                            sumidx = 'nocc'
                        elif self.summation[0] in virtual:
                            sumidx = 'nvirt'
                        for i in range(1, len(self.summation)):
                            if self.summation[i] in occupied:
                                sumidx = sumidx + '*nocc'
                            elif self.summation[i] in virtual:
                                sumidx = sumidx + '*nvirt'


        if self.num_factor != 0.0:
            if self.num_factor > 0 :
                if self.num_factor == 1.0:
                    return "+"+permidx + sumidx + cidx + opidx + delidx 
                else:    
                    f = Fraction(self.num_factor).limit_denominator()
                    num = f.numerator
                    den = f.denominator
                    if den == 1:
                        nf = "{num}".format(num=num)
                    else:
                        nf = "\\frac{{{num}}}{{{den}}}".format(num=num, den=den)
                    return "+"+ nf +permidx + sumidx + cidx + opidx + delidx 

            elif self.num_factor < 0:
                if self.num_factor == -1.0:
                    return "-"+permidx+ sumidx + cidx + opidx + delidx 
                else:    
                    f =Fraction(self.num_factor).limit_denominator()
                    num= f.numerator
                    num *= (-1)
                    den= f.denominator
                    if den == 1:
                        nf = "-{num}".format(num=num)
                    else:
                        nf = "-\\frac{{{num}}}{{{den}}}".format(num=num, den=den)
                    return nf +permidx+ sumidx + cidx + opidx + delidx 
        else:
            return str(0)

    def exec_delta_fixed(self, lst):
        changed = lst[1]
        unchanged = lst[0]
        self.substitute(changed, unchanged)
        
    def exec_delta(self):
        """Perform substitution of indices / reduction of 
        summation indices resulting from Kronecker delta,
        \delta_{pq}. Numerical factor may be changed to zero 
        if p and q indices belong to disjoint sets.
        """
        delta = deepcopy(self.delta)
        ch = []
        for p, q in delta:
            for x in range(0, len(ch)):
                if p == ch[x][0]:
                    p = ch[x][1]
                if q == ch[x][0]:
                    q = ch[x][1]

            if p not in self.summation and \
                    q not in self.summation:
                #
                # Both indices fixed
                #
                pq = set([p,q])
                # if not (set(virtualall) and pq == pq or set(occupied) and pq == pq \
                #         or set(general) and pq == pq):
                #     self.num_factor = 0.0
                if (p in virtualall and q in occupied) or  \
                   (p in occupied and q in virtual):
                    self.num_factor = 0.0
            else:
                if p in self.summation and \
                        q not in self.summation:
                    #
                    # Summation over p; q is fixed
                    #
                    changed = p
                    unchanged = q
                elif q in self.summation and \
                        p not in self.summation:
                    #
                    # Summation over q; p is fixed
                    #

                    changed = q
                    unchanged = p
                else:
                    #
                    # There is summation over
                    # both indices
                    #
                    # print('ZOBACZ CZY TO DOBRZE DZIALA, WYCHODZE')
                    # sys.exit(0)
                    # print('zobacz tu1', self)
                    l1 = [p,q]
                    # print('zobacz tu2', l1)
                    l1.sort()
                    # print('zobacz tu3', l1)
                    changed = l1[0]
                    unchanged = l1[1]


                if changed in self.summation:
                    self.summation.remove(changed)
                #
                # Remove delta symbol
                #
                del self.delta[self.delta.index([p, q])]
                #
                # Substitute CHANGED index for UNCHANGED
                #
                self.substitute(changed, unchanged)
                # print('zobacz tu4', self)
                # print('')

                ch.append([changed, unchanged])


    def indices(self):
        """Return set of indices occuring in self. Both
        summation and fixed indices are returned.
        """
       
        idx = set(self.summation)
        for x in self.operator_idx:
            idx = idx | set(x)

        for x in self.coefficient_idx:
            idx = idx | set(x)

        for x in self.delta:
            idx = idx | set(x)

        return idx

    def set_fixed(self):
        a = set(self.summation)
        b = set(self.summation)
        for w in self.coefficient_idx:
            b = b | set(w)
        for w in self.operator_idx:
            b = b | set(w)
        for w in self.delta:
            b = b | set(w)
        #
        # Set c contains fixed indices
        #
        c = b - a

        return c


    def substitute(self, p, q):
        """Substitute any instance of index p for index q.
        Substitution is performed in both summation and fixed
        indices subsets.
        """

        for x in range(0, len(self.coefficient)):
            for y in self.coefficient_idx[x]:
                if p in self.coefficient_idx[x]:
                    cidx = self.coefficient_idx[x].index(p)
                    self.coefficient_idx[x][cidx] = q

        if p in self.operator_idx:
            opidx = self.operator_idx.index(p)
            self.operator_idx[opidx] = q

        for delta in self.delta:
            if delta[0] == p:
                delta[0] = q
            if delta[1] == p:
                delta[1] = q
        #
        # Search for redundant \delta_{qq}
        #

        for i in range(0, len(self.delta)):
            if self.delta[i] == [q, q]:
                del self.delta[i]
                break
        
        if p in self.summation:
            idx = self.summation.index(p)
            self.summation[idx] = q


    def new_delta(self, p, q, spin1='', spin2=''):
        """
        Append Kronecker delta, \delta_{pq} to the cas
        instance. If appending \delta_{pq} is pointless
        regarding items that are already in self.delta list,
        do nothing.
        """

        if p == q:
            return

        if len(self.delta) == 0:
            self.delta.append([p, q])
            if spin1 != '':
                self.delta_spin.append([spin1, spin2])

        else:
            pointless = False
            for d in self.delta:
                if d == [p, q] or d == [q, p]:
                    pointless = True
                    break

            if not pointless:
                self.delta.append([p, q])
                if spin1 != '':
                    self.delta_spin.append([spin1, spin2])

    def left_split(self):
        """LEFT >> RIGHT = SELF, LEFT = p^+q """
        left = cas()
        right = deepcopy(self)

        left.operator_idx = [deepcopy(right.operator_idx[0])]
        left.operator_idx.append(deepcopy(right.operator_idx[1]))
        left.operator_type = [deepcopy(right.operator_type[0])]
        left.operator_type.append(deepcopy(right.operator_type[1]))

        right.operator_idx.remove(right.operator_idx[0])
        right.operator_idx.remove(right.operator_idx[0])

        right.operator_type.remove(right.operator_type[0])
        right.operator_type.remove(right.operator_type[0])

        return left, right


    def right_split(self):
        """ LEFT << RIGHT = SELF, RIGHT = p^+q """

        left = deepcopy(self)
        right = cas()

        s = len(left.operator_idx) - 1
        right.operator_idx = [deepcopy(left.operator_idx[s-1])]
        right.operator_idx.append(deepcopy(left.operator_idx[s]))

        right.operator_type = [deepcopy(left.operator_type[s-1])]
        right.operator_type.append(deepcopy(left.operator_type[s]))

        left.operator_idx.remove(left.operator_idx[-1])
        left.operator_idx.remove(left.operator_idx[-1])

        left.operator_type.remove(left.operator_type[-1])
        left.operator_type.remove(left.operator_type[-1])

        return left, right



    def fromleft(self, e):
        """Concatenate from left: E >> SELF"""

        res = cas()
        res.summation       = deepcopy(e.summation) + deepcopy(self.summation)
        res.coefficient     = deepcopy(e.coefficient) + deepcopy(self.coefficient)
        res.coefficient_idx = deepcopy(e.coefficient_idx) + deepcopy(self.coefficient_idx)
        res.operator_idx    = deepcopy(e.operator_idx) + deepcopy(self.operator_idx)
        res.coefficient_spin    = deepcopy(e.coefficient_spin) + deepcopy(self.coefficient_spin)
        res.delta_spin    = deepcopy(self.delta_spin) + deepcopy(e.delta_spin)
        res.operator_type  = deepcopy(e.operator_type) + deepcopy(self.operator_type)
        res.num_factor      = e.num_factor * self.num_factor

        if self.delta != []:
            for delta in self.delta:
                res.new_delta(delta[0], delta[1])

        if e.delta != []:
            for delta in e.delta:
                res.new_delta(delta[0], delta[1])
        return res 
    
    
    def fromright(self, e):
        """Concatenate from right: SELF << E"""

        res = cas()
        res.summation       = deepcopy(self.summation) + deepcopy(e.summation)
        res.coefficient     = deepcopy(self.coefficient) + deepcopy(e.coefficient)
        res.coefficient_idx = deepcopy(self.coefficient_idx) + deepcopy(e.coefficient_idx)
        res.operator_idx    = deepcopy(self.operator_idx) + deepcopy(e.operator_idx)
        res.coefficient_spin    = deepcopy(self.coefficient_spin) + deepcopy(e.coefficient_spin)
        res.delta_spin    = deepcopy(self.delta_spin) + deepcopy(e.delta_spin)
        res.operator_type  = deepcopy(self.operator_type) + deepcopy(e.operator_type)
        res.num_factor      = self.num_factor * e.num_factor
        
        if self.delta != []:
            for delta in self.delta:
                res.new_delta(delta[0], delta[1])

        if e.delta != []:
            for delta in e.delta:
                res.new_delta(delta[0], delta[1])
        return res       

    
    def rename_as_density(self):
        # print('')
        # print('rendens', self)
        if len(self.operator_idx) == 2:
            self.coefficient.append(DENS1)
            idx_list = [self.operator_idx[1], self.operator_idx[0]]            
            self.coefficient_idx.append(idx_list)
            
            self.operator_idx = []
        if len(self.operator_idx) == 4:
            self.coefficient.append(DENS2)
            idx_list = [self.operator_idx[3], self.operator_idx[2], self.operator_idx[0], self.operator_idx[1]]            
            self.coefficient_idx.append(idx_list)
            self.operator_idx = []

        if len(self.operator_idx) == 6:
            self.coefficient.append(DENS3)
            # wersja I
            # p q t
            # s r v
            # jest rownowazne p*q*t* vrs
            #                 0 1 2  345
            #                 [0,5] [1,4],[2,3] 
            # indeksy idą parami [p,s],[q,r],[t,v]
            # idx_list = [self.operator_idx[0], self.operator_idx[5], \
            #             self.operator_idx[1], self.operator_idx[4], \
            #             self.operator_idx[2], self.operator_idx[3]]
            # wersja II, zeby spiny sie zmiescily
            idx_list = [self.operator_idx[5], self.operator_idx[4], \
                        self.operator_idx[3], self.operator_idx[0], \
                        self.operator_idx[1], self.operator_idx[2]]

            
            self.coefficient_idx.append(idx_list)
            self.operator_idx = []

            
        if len(self.operator_idx) == 8:
            self.coefficient.append(DENS4)
            # wersja I
            # p q t w
            # s r v x
            # jest rownowazne p*q*t*w* vrsx
            #                 0 1 2 3  4567
            #                 [0,7] [1,6],[2,5],[3,4] 
            # indeksy idą parami [p,s],[q,r],[t,v]
            # idx_list = [self.operator_idx[0], self.operator_idx[7], \
            #             self.operator_idx[1], self.operator_idx[6], \
            #             self.operator_idx[2], self.operator_idx[5], \
            #             self.operator_idx[3], self.operator_idx[4]]
            # wersja II, zeby spiny sie zmiescily
            idx_list = [self.operator_idx[7], self.operator_idx[6], \
                        self.operator_idx[5], self.operator_idx[4], \
                        self.operator_idx[0], self.operator_idx[1], \
                        self.operator_idx[2], self.operator_idx[3]]

            self.coefficient_idx.append(idx_list)
            self.operator_idx = []


        # print('rendens2', self)
        # print('')


    def contraction(self, list_of_c):
        # --> tak wyglada list kontrakcji[[0,1],[2,3]..]

        new_operator_idx = []
        new_operator_type = []
        deleted_list = []
        # print(self.operator_type, list_of_c[0][0])
        # print('l-------------------------------l', list_of_c)
        deleted = []
        for i in range(0, len(list_of_c)):
#            print(i, 'list of c', self, list_of_c[i])
            if self.operator_type[list_of_c[i][0]] == ANI and self.operator_type[list_of_c[i][1]]== CRE:
#                print('tak')
                idx1 = self.operator_idx[list_of_c[i][0]]
                idx2 = self.operator_idx[list_of_c[i][1]]
                self.delta.append([idx1, idx2])
                num = (list_of_c[i][1]-list_of_c[i][0]) - 1
                for j in range(list_of_c[i][0]+1, list_of_c[i][1]):
                    if j in deleted:
                        num = num - 1
#                print('numnum', num)
                if (num%2!= 0):
                    self.num_factor *= -1.0
 #                   print('zmieniam znak')
                deleted.append(list_of_c[i][0])
                deleted.append(list_of_c[i][1])
#                print('sra', self)


            deleted_list.append(list_of_c[i][0])
            deleted_list.append(list_of_c[i][1])

        # print('del', deleted_list)
        for i in range(0, len(self.operator_idx)):
            # print(self.operator_idx[i])
            if i not in deleted_list:
                new_operator_idx.append(self.operator_idx[i])
                new_operator_type.append(self.operator_type[i])

        # print(new_operator_idx)
        self.operator_idx = new_operator_idx
        self.operator_type = new_operator_type


    def wick_ca(self):
        # print('wynik WIKA', self)
        contr = self.all_contractions()
        nonzero_contr = self.remove_null_contractions(contr)
#        print('contr', nonzero_contr)
        n = len(nonzero_contr)
        # print('')
        # print(self, '----->', nonzero_contr)
        
        for x in nonzero_contr:
            x.sort()
          #  print(x)
        # print('')
#        sys.exit(0)
        result_list = []
        for i in range(0, n+1):
            result_list.append(deepcopy(self))

        
        result_list[0].normal_order()
        # print('res0', result_list[0])
        for i in range(1, n+1):
            # print('')
            # print('przed',i,  result_list[i], nonzero_contr[i-1])            
            result_list[i].contraction(nonzero_contr[i-1])
            result_list[i].normal_order()
            # print('po', result_list[i])
            # print('')
            
        return result_list
    

    def remove_null_contractions(self, contr):

        nonzero_single_contr = [ANI, CRE]
        nonzero_contr = []
        for i in range(0, len(contr)):
            iszero = False
            for j in range(0, len(contr[i])):
                op0 = self.operator_type[contr[i][j][0]]
                op1 = self.operator_type[contr[i][j][1]]                
                single_contr = [op0, op1]
                if single_contr != nonzero_single_contr:
                    iszero = True
                    # print('usuwam', self, contr[i], self.operator_idx[contr[i][j][0]], self.operator_type[contr[i][j][0]], \
                    #       self.operator_idx[contr[i][j][1]], self.operator_type[contr[i][j][1]])
                    break
            if not iszero:
                nonzero_contr.append(contr[i])

        for x in nonzero_contr:
            x.sort()

        contr_without_rep = []
        for x in nonzero_contr:
            if x not in contr_without_rep:
                contr_without_rep.append(x)
                

        # print('niezerowe')
        # for x in nonzero_contr:
        #     print(x)
        # return(nonzero_contr)
        return(contr_without_rep)
                           
            
    def all_contractions(self):
        # print('kontrakcje', self)

        if len(self.operator_idx) == 2:
            return []
        idx_list = list(range(len(self.operator_idx)))

        comb_set = list(combinations(idx_list, 2))
        for i in range(0, len(comb_set)):
            comb_set[i] = [list(comb_set[i])]
        
        comb_set_all = comb_set
        comb_set2=[]

        # for i in range(1, len(self.operator_idx)//2):
        if len(self.operator_idx) > 2:
            # print('comb_set', comb_set)
            
            idx_list2 = list(range(len(comb_set)))
            comb_set2_temp = list(combinations(idx_list2, 2))
            # print(idx_list2)

            comb_set2=[]
            for i in range(0, len(comb_set2_temp)):
                comb_set2_temp[i] = list(comb_set2_temp[i])

                bigl1 = comb_set[comb_set2_temp[i][0]]
                bigl2 = comb_set[comb_set2_temp[i][1]]
                # print('bb', bigl1)
                # print('bb', bigl2)
                l1 = []
                l2 = []
                for k1 in bigl1:
                    l1 = l1 + k1
                for k2 in bigl2:
                    l2 = l2 + k2


                # print('l1l2', l1,l2)
                llist = l1+l2
                llset = set(l1+l2)
                contain_duplicates = len(llist)!=len(llset)
                ap_list = []
                if not contain_duplicates:
                    for k1 in bigl1:                        
                        ap_list.append(k1)
                    for k2 in bigl2:                        
                        ap_list.append(k2)    

                    comb_set2.append(ap_list)
            # print('bede dodawac to', comb_set2)
            comb_set_all = comb_set_all + comb_set2
            # print('po dodaniu')
            # print('')

        # for x in comb_set_all:
        #     print(x)

        if len(comb_set2) > 0:
            comb_set_last = comb_set2
        else:
            comb_set_last = []
        comb_set_new = []
        for k in range(1, len(self.operator_idx)//2-1):
            # print('ROBIE',k, 'RAZ')
            # print('comb_set', comb_set)

            for i in range(0, len(comb_set)):
                for j in range(0, len(comb_set_last)):

                    bigl1 = comb_set[i]
                    bigl2 = comb_set_last[j]

                    # print('bb', bigl1)
                    # print('bb', bigl2)
                    l1 = []
                    l2 = []
                    for k1 in bigl1:
                        l1 = l1 + k1
                    for k2 in bigl2:
                        l2 = l2 + k2


                    # print('l1l2', l1,l2)
                    llist = l1+l2
                    llset = set(l1+l2)
                    contain_duplicates = len(llist)!=len(llset)
                    ap_list = []
                    if not contain_duplicates:
                        for k1 in bigl1:                        
                            ap_list.append(k1)
                        for k2 in bigl2:                        
                            ap_list.append(k2)    

                    if len(ap_list)>0:
                        comb_set_new.append(ap_list)
            comb_set_all = comb_set_all + deepcopy(comb_set_new)
            comb_set_last = deepcopy(comb_set_new)


        # print('po dodaniu')
        # for x in comb_set_all:
        #     print(x)
        # print('')

        return comb_set_all

    def copy(self):

        elem2 = cas()
        elem2.summation = deepcopy(self.summation)
        elem2.coefficient = deepcopy(self.coefficient)
        elem2.coefficient_idx = deepcopy(self.coefficient_idx)
        elem2.coefficient_spin = deepcopy(self.coefficient_spin)
        elem2.num_factor = deepcopy(self.num_factor)
        elem2.delta = deepcopy(self.delta)
        elem2.delta_spin = deepcopy(self.delta_spin)
            
        return(elem2)
            
    def normal_order(self):

        creator_list = []
        anlator_list = []
        original_list = deepcopy(self.operator_idx)

        for i in range(0, len(self.operator_idx)):
            if self.operator_type[i] == CRE:

                creator_list.append(self.operator_idx[i])
            elif self.operator_type[i] == ANI:
                anlator_list.append(self.operator_idx[i])


        # print(creator_list)
        # print(anlator_list)        
        c_trans = swap_count(creator_list)
        a_trans = swap_count(anlator_list)
        # print(c_trans)
        # print(a_trans)

        c_num = 1.0
        a_num = 1.0
        if (c_trans%2!=0):
            c_num = -1.0
        if (a_trans%2!=0):
            a_num = -1.0
        num = c_num*a_num
        # print(c_num, a_num, num)

        num_anlator_pass = count_anlator_pass(self.operator_idx, self.operator_type)
        # print(num_anlator_pass)
        if (num_anlator_pass%2!=0):
            num *= -1.0
        # print('num', num)

        creator_list.sort()
        anlator_list.sort()
        ca_list = creator_list + anlator_list
        original_op_type_list = deepcopy(self.operator_type)
        
        # print('-----------')
        # print(original_list)
        # print(ca_list)
        self.operator_type = []
        for i in range(0, len(ca_list)):
            # print(i, ca_list[i], original_list.index(ca_list[i]))
            k = original_list.index(ca_list[i])
            self.operator_type.append(original_op_type_list[k])

        # print('-----------')

        
        self.operator_idx = creator_list + anlator_list
        # print(creator_list, anlator_list)
        self.num_factor *= num

def count_anlator_pass(crean_list, crean_type):

    l = len(crean_list)-1

    pass_count = 0
    # print('oto')
    # print(crean_list)
    # print(crean_type)
    for i in range(l, -1, -1):
        if crean_type[i] == CRE:            
            for j in range(i-1, -1, -1):
                if crean_type[j] == ANI:
                    pass_count += 1

    return pass_count

def swap_count(input_arr):

    pos = sorted(list(enumerate(input_arr)), key=lambda x: x[1])
    cnt = 0

    for index in range(len(input_arr)):
        while True:
            if (pos[index][0] == index):
                break
            else:
                cnt += 1
                swap_index = pos[index][0]
                pos[index], pos[swap_index] = pos[swap_index], pos[index]

    return cnt


def cas_to_ugg(res):

    res2 = arithmetic_string()
    for x in range(0, len(res)):

        if (len(res[x].operator_idx) > 0):
            print('OPERATORA STILL PRESENT, CANNOT TRANSFORM CAS TO UGG')
            sys.exit(1)
        y = ugg()
        y.summation = res[x].summation
        y.coefficient = res[x].coefficient
        y.coefficient_idx = res[x].coefficient_idx
        y.num_factor = res[x].num_factor
        y.delta = res[x].delta
        res2.append(deepcopy(y))
        # res2 = res2 + arithmetic_string(deepcopy(y))
    return res2

def ugg_to_cas(res):

    res2 = arithmetic_string()
    for x in range(0, len(res)):
        
        if (len(res[x].operator_idx) > 0):
            print('OPERATORA STILL PRESENT, CANNOT TRANSFORM CAS TO UGG')
            sys.exit(1)
        # print('o  taki', res[x])
        y = cas()
        y.summation = res[x].summation
        y.coefficient = res[x].coefficient
        y.coefficient_idx = res[x].coefficient_idx
        y.num_factor = res[x].num_factor
        y.delta = res[x].delta
        # print('na taki', y)
        res2.append(deepcopy(y))
        
        # res2 = res2 + arithmetic_string(deepcopy(y))
    return res2



def add_spin_from_list(res, spin_dict, numf):


    for i in range(0, len(res)):
        elem = res[i]
        # print('')
        # print('sspin przed', elem, spin_dict)
        elem.num_factor *= numf
        for a in range(0, len(elem.coefficient_idx)):
            elem.coefficient_spin.append([])
        for a in range(0, len(elem.delta)):
            elem.delta_spin.append([])


        for a in range(0, len(elem.coefficient_idx)):
            for b in range(0, len(elem.coefficient_idx[a])):
                if elem.coefficient_idx[a][b] in spin_dict.keys():
                    elem.coefficient_spin[a].append(spin_dict[elem.coefficient_idx[a][b]])
                else:
                    elem.coefficient_spin[a].append('0')
        for a in range(0, len(elem.delta)):
            for b in range(0, len(elem.delta[a])):
                if elem.delta[a][b] in spin_dict.keys():
                    elem.delta_spin[a].append(spin_dict[elem.delta[a][b]])
                else:
                    elem.delta_spin[a].append('0')
        # print('sspin poo', elem)
        # print('')


    return res

def add_spin_summ(res):

    res2 = arithmetic_string()
    for i in range(0, len(res)):
        elem = res[i]
        # print('elem000', elem)
        if len(elem.summation)==0:
            res2.append(elem)
            # res2 = res2 + arithmetic_string(elem)            
        else:
            # if len(elem.summation)==3:
                # print('')

            result = add_spin_summ_basic(elem, elem.summation)
            for x in result:
                res2.append(x)
                # res2 = res2 + add_spin_summ_basic(elem, elem.summation)
                # print('')
                # for f in res3:
                #     print('ppp', f)
    return res2

def clean_cas(res):

    res2 = []
    for x in res:
        if x.num_factor != 0.0:
            res2.append(x)
    return res2

def res_1gam_only(res):

    res_1g = []
    rest = []

    for x in res:
        gam2 = False
        for y in x.coefficient:
            if y == DENS2:
                gam2 = True

        if gam2:
            rest.append(x)
        else:
            res_1g.append(x)

    return res_1g, rest

def res_divied_gammas(res):


    # this works only if three is one gamma2 in res[i]
    res_pppp = []
    res_pmpm = []
    res_mpmp = []

    for x in res:
        for y in range(0, len(x.coefficient)):
            if x.coefficient[y] == DENS2:
                s = gam2_spin_type(x.coefficient_spin[y])
                if s == '4p' or s == '4m':
                    res_pppp.append(x)
                    break
                # elif s == '4m':
                #     res_mmmm.append(x)
                #     break
                elif s == 'pm':
                    res_pmpm.append(x)
                    break
                elif s == 'mp':
                    res_mpmp.append(x)
                    break
                
    return res_pppp, res_pmpm, res_mpmp

def gam2_spin_type(lst):

    p = 0
    m = 0
    
    for x in lst:
        if x == '+':
            p+=1
        elif x == '-':
            m+=1
        else:
            print('something weird in spin list')
            print(lst)
            sys.exit(1)


    if p == 4 and m == 0:
        return '4p'
    elif m == 4 and p == 0:
        return '4m'
    elif p == 2 and m == 2:
        if lst[0] == '+':
            return 'pm'
        elif lst[0] == '-':
            return 'mp'

    else:
        print('something weird in spin list')
        print('a', lst)
        sys.exit(1)
    

def add_spin_summ_basic(elem, sumi):

    ll = len(sumi)

    elem2 = elem.copy()
    for a in range(0, len(elem.coefficient_idx)):
        for b in range(0, len(elem.coefficient_idx[a])):
            if sumi[-1] == elem.coefficient_idx[a][b]:
                elem.coefficient_spin[a][b] = '+'
                print('elem', elem)
                elem2.coefficient_spin[a][b] = '-'
                print('elem', elem)
                print('elem2', elem2)
                print('------------------------------------------------------------------------------')

    if len(sumi[0:ll-1])==0:
        result = arithmetic_string(elem, elem2)
        return result
    else:
        return add_spin_summ_basic(elem, sumi[0:ll-1]) + add_spin_summ_basic(elem2, sumi[0:ll-1])

def simplify_spin_initially(res):

    res2 = arithmetic_string()
    for x in range(0, len(res)):

        
        # print(res[x].coefficient_spin)
        # print(res[x].delta_spin)
        elem = res[x]

        if ['+', '-'] in elem.delta_spin or ['-', '+'] in elem.delta_spin:
            elem.num_factor = 0
            continue
        else:
            if ['+', '-'] in elem.coefficient_spin or ['-', '+'] in elem.coefficient_spin:
                elem.num_factor = 0
                continue
        res2.append(elem)

                
            
        # for minidel in elem.delta_spin:
        #     if minidel[0] != minidel[1]:
        #         if minidel[0]!=0 and minidel[1]!=0:
        #             elem.num_factor = 0

        # for x in range(0, len(elem.coefficient)):
        #     coef = deepcopy(elem.coefficient[x])
        #     coef_spin = deepcopy(elem.coefficient_spin[x])
            
        #     if coef == BARENUCL_HAM:
        #         if coef_spin[0] != coef_spin[1]:
        #             if coef_spin[0] != '0' and coef_spin[1] != '0':
        #                 elem.num_factor = 0
        #     if coef == DENS1:
        #         if coef_spin[0] != coef_spin[1]:
        #             if coef_spin[0] != '0' and coef_spin[1] != '0':
        #                 elem.num_factor = 0

        #     if coef == TWOEL_INT_DIRAC:
        #         if coef_spin[0] != coef_spin[2]:
        #             if coef_spin[0] != '0' and coef_spin[2] != '0':
        #                 elem.num_factor = 0
        #         if coef_spin[1] != coef_spin[3]:
        #             if coef_spin[1] != '0' and coef_spin[3] != '0':
        #                 elem.num_factor = 0

    #     res2.append(elem)
    #     # res2 = res + arithmetic_string(elem)                                

    # resout = arithmetic_string()
    # for x in res:
    #     if (x.num_factor) != 0.0:
    #         resout
    #         resout = resout + arithmetic_string(x)
            
    return res2



def simplify_spin_onedet(res):

    res2 = arithmetic_string()
    for x in range(0, len(res)):

        elem = res[x]
        if ['+', '-'] in elem.delta_spin or ['-', '+'] in elem.delta_spin:
            elem.num_factor = 0
            continue
        else:
            if ['+', '-'] in elem.coefficient_spin or ['-', '+'] in elem.coefficient_spin:
                elem.num_factor = 0
                continue

        for x in range(0, len(elem.coefficient)):
            if coef == TWOEL_INT_DIRAC:
                if coef_spin[0] != coef_spin[2]:
                    elem.num_factor = 0
                    continue
                if coef_spin[1] != coef_spin[3]:
                    elem.num_factor = 0
                    continue

        
        res2.append(elem)

    return res2

def simplify_spin(res):

    res2 = arithmetic_string()
    for x in range(0, len(res)):

        elem = res[x]
        #print('elem przed', elem)

#        print('sspin przed', elem)
        # find delta = 0
        for minidel in elem.delta_spin:
            if minidel[0] != minidel[1]:
                elem.num_factor = 0

        for x in range(0, len(elem.coefficient)):
            coef = deepcopy(elem.coefficient[x])
            coef_spin = deepcopy(elem.coefficient_spin[x])
            
            if coef == BARENUCL_HAM:
                if coef_spin[0] != coef_spin[1]:
                    elem.num_factor = 0
            if coef == DENS1:
                if coef_spin[0] != coef_spin[1]:
                    if coef_spin[0] != '0' and coef_spin[1] != '0':
                        # print('z powodu DENS1')
                        elem.num_factor = 0
            if coef == TWOEL_INT_DIRAC:
                if coef_spin[0] != coef_spin[2]:
                    elem.num_factor = 0
                if coef_spin[1] != coef_spin[3]:
                    elem.num_factor = 0
            if coef == DENS2:
                # print('sss', coef_spin)
                # +--+
                # -++-
                # -+-+                
                if coef_spin[2]=='-' and coef_spin[3]=='+':
                    sp = deepcopy(elem.coefficient_spin[x][3])
                    id = deepcopy(elem.coefficient_idx[x][3])
                    elem.coefficient_idx[x][3]=deepcopy(elem.coefficient_idx[x][2])
                    elem.coefficient_spin[x][3]=deepcopy(elem.coefficient_spin[x][2])
                    elem.coefficient_idx[x][2]=deepcopy(id)
                    elem.coefficient_spin[x][2]=deepcopy(sp)
                    elem.num_factor *= -1.0
                if coef_spin[0]=='-' and coef_spin[1]=='+':
                    sp = deepcopy(elem.coefficient_spin[x][1])
                    id = deepcopy(elem.coefficient_idx[x][1])
                    elem.coefficient_idx[x][1]=deepcopy(elem.coefficient_idx[x][0])
                    elem.coefficient_spin[x][1]=deepcopy(elem.coefficient_spin[x][0])
                    elem.coefficient_idx[x][0]=deepcopy(id)
                    elem.coefficient_spin[x][0]=deepcopy(sp)
                    elem.num_factor *= -1.0
                # if coef_spin[0]==coef_spin[3] and \
                #    coef_spin[0]!=coef_spin[1] and \
                #    coef_spin[1]==coef_spin[2]:
                #     print('przedelem', elem)

                #     sp = deepcopy(elem.coefficient_spin[x][3])
                #     id = deepcopy(elem.coefficient_idx[x][3])
                #     elem.coefficient_idx[x][3]=deepcopy(elem.coefficient_idx[x][2])
                #     elem.coefficient_spin[x][3]=deepcopy(elem.coefficient_spin[x][2])
                #     elem.coefficient_idx[x][2]=deepcopy(id)
                #     elem.coefficient_spin[x][2]=deepcopy(sp)

                #     elem.num_factor *= -1.0
        #print('spin po    ', elem)
        res2.append(elem)                                

    # resout = arithmetic_string()
    # for x in res:
    #     if (x.num_factor) != 0.0:
    #         resout = resout + arithmetic_string(x)

    return res2

def simplify_2rdm_occ(res):

    res2 = arithmetic_string()
    for x in range(0, len(res)):
        elem = res[x]

        if DENS2 not in elem.coefficient and DENS1 not in elem.coefficient:
            res2 = res2 + arithmetic_string(elem)
        else:
            for x in range(0, len(elem.coefficient)):
                coef = deepcopy(elem.coefficient[x])
                coef_idx = deepcopy(elem.coefficient_idx[x])
                coef_spin = deepcopy(elem.coefficient_spin[x])
                if coef == DENS1:
                    if coef_idx[0]!=coef_idx[1]:
                        coef_sort = deepcopy(coef_idx)
                        coef_sort.sort()
                        elem.delta.append([coef_idx[0], coef_idx[1]])
                        elem.delta_spin.append([coef_spin[0], coef_spin[1]])
                        res2 = res2 + arithmetic_string(elem)
                        
                        #                print(coef, coef_spin)
                elif coef == DENS2:                    
                    if coef_spin[0]==coef_spin[1] and \
                       coef_spin[1]==coef_spin[2] and \
                       coef_spin[2]==coef_spin[3]:
                        # print('gamma plusowa', elem)
                        # print('elemstart', elem)
                        elem.coefficient_idx.pop(x)
                        elem.coefficient.pop(x)
                        elem.coefficient_spin.pop(x)                        
                        elem2 = elem.copy()
                        # print(coef, coef_spin)
                        elem.coefficient.append(DENS1)                        
                        elem.coefficient_idx.append([coef_idx[0], coef_idx[2]])
                        elem.coefficient_spin.append(['+', '+'])
                        elem.delta.append([coef_idx[0], coef_idx[2]])
                        elem.delta_spin.append(['+', '+'])
                        # print(elem)
                        elem.coefficient.append(DENS1)                        
                        elem.coefficient_idx.append([coef_idx[1], coef_idx[3]])
                        elem.coefficient_spin.append(['+', '+'])
                        elem.delta.append([coef_idx[1], coef_idx[3]])
                        elem.delta_spin.append(['+', '+'])

                        elem2.coefficient.append(DENS1)
                        elem2.coefficient_idx.append([coef_idx[0], coef_idx[3]])
                        elem2.coefficient_spin.append(['+', '+'])
                        elem2.delta.append([coef_idx[0], coef_idx[3]])
                        elem2.delta_spin.append(['+', '+'])

                        elem2.coefficient.append(DENS1)
                        elem2.coefficient_idx.append([coef_idx[1], coef_idx[2]])
                        elem2.coefficient_spin.append(['+', '+'])
                        elem2.delta.append([coef_idx[1], coef_idx[2]])
                        elem2.delta_spin.append(['+', '+'])

                        elem2.num_factor *= -1.0
                        res2 = res2 + arithmetic_string(elem)
                        res2 = res2 + arithmetic_string(elem2)
                        # print('')
                        # print(elem)
                        # print(elem2)
                    if coef_spin[0]==coef_spin[2] and \
                       coef_spin[0]!=coef_spin[1] and \
                       coef_spin[1]==coef_spin[3]:
                        # print('gamma minsowa', elem)
                        # print('elemstart', elem)
                        elem.coefficient_idx.pop(x)
                        elem.coefficient.pop(x)
                        elem.coefficient_spin.pop(x)                        

                        # print(coef, coef_spin)
                        elem.coefficient.append(DENS1)                        
                        elem.coefficient_idx.append([coef_idx[0], coef_idx[2]])
                        elem.coefficient_spin.append(['+', '+'])
                        elem.delta.append([coef_idx[0], coef_idx[2]])
                        elem.delta_spin.append(['+', '+'])

                        # print(elem)
                        elem.coefficient.append(DENS1)                        
                        elem.coefficient_idx.append([coef_idx[1], coef_idx[3]])
                        elem.coefficient_spin.append(['-', '-'])
                        elem.delta.append([coef_idx[1], coef_idx[3]])
                        elem.delta_spin.append(['-', '-'])

                        res2 = res2 + arithmetic_string(elem)
                        # print('')
                        # print(elem)




    return res2

# def simplify_34rdm(res):

def simplify_1rdm_mult(res):

    res2 = arithmetic_string()
    for x in range(0, len(res)):
        elem = res[x]

        if DENS1 not in elem.coefficient:
            res2 = res2 + arithmetic_string(elem)
        else:
            for x in range(0, len(elem.coefficient)):
                coef = deepcopy(elem.coefficient[x])
                coef_idx = deepcopy(elem.coefficient_idx[x])
                coef_spin = deepcopy(elem.coefficient_spin[x])
                if coef == DENS1:
                    if coef_idx[0]!=coef_idx[1]:
                        elem.delta.append([coef_idx[0], coef_idx[1]])
                        elem.delta_spin.append([coef_spin[0], coef_spin[1]])
                        if coef_spin[0] != coef_spin[1]:
                            elem.num_factor = 0.0
                        else:
                            if coef_spin[0] == '+':
                                elem.coefficient[x] = DENS1P
                            elif coef_spin[0] == '-':
                                elem.coefficient[x] = DENS1M
                                
                        res2 = res2 + arithmetic_string(elem)

    return res2

def simplify_2rdm_mult_occ(res):

    res2 = arithmetic_string()
    for x in range(0, len(res)):
        elem = res[x]

        if DENS2 not in elem.coefficient:
            res2 = res2 + arithmetic_string(elem)
        else:
            for x in range(0, len(elem.coefficient)):
                coef = deepcopy(elem.coefficient[x])
                coef_idx = deepcopy(elem.coefficient_idx[x])
                coef_spin = deepcopy(elem.coefficient_spin[x])
                        
                if coef == DENS2:                    
                    if coef_spin[0]==coef_spin[1] and \
                       coef_spin[1]==coef_spin[2] and \
                       coef_spin[2]==coef_spin[3]:
#                        print('gamma plusowa', elem)
                        # print('elemstart', elem)
                        elem.coefficient_idx.pop(x)
                        elem.coefficient.pop(x)
                        elem.coefficient_spin.pop(x)                        
                        elem2 = elem.copy()
                        # print(coef, coef_spin)
                        elem.coefficient.append(DENS1)                        
                        elem.coefficient_idx.append([coef_idx[0], coef_idx[2]])
                        elem.coefficient_spin.append(['+', '+'])
                        elem.delta.append([coef_idx[0], coef_idx[2]])
                        elem.delta_spin.append(['+', '+'])
                        # print(elem)
                        elem.coefficient.append(DENS1)                        
                        elem.coefficient_idx.append([coef_idx[1], coef_idx[3]])
                        elem.coefficient_spin.append(['+', '+'])
                        elem.delta.append([coef_idx[1], coef_idx[3]])
                        elem.delta_spin.append(['+', '+'])

                        elem2.coefficient.append(DENS1)
                        elem2.coefficient_idx.append([coef_idx[0], coef_idx[3]])
                        elem2.coefficient_spin.append(['+', '+'])
                        elem2.delta.append([coef_idx[0], coef_idx[3]])
                        elem2.delta_spin.append(['+', '+'])

                        elem2.coefficient.append(DENS1)
                        elem2.coefficient_idx.append([coef_idx[1], coef_idx[2]])
                        elem2.coefficient_spin.append(['+', '+'])
                        elem2.delta.append([coef_idx[1], coef_idx[2]])
                        elem2.delta_spin.append(['+', '+'])

                        elem2.num_factor *= -1.0
                        res2 = res2 + arithmetic_string(elem)
                        res2 = res2 + arithmetic_string(elem2)
                        # print('')
                        # print(elem)
                        # print(elem2)
                    if coef_spin[0]==coef_spin[2] and \
                       coef_spin[0]!=coef_spin[1] and \
                       coef_spin[1]==coef_spin[3]:
                        # print('gamma minsowa', elem)
                        # print('elemstart', elem)
                        elem.coefficient_idx.pop(x)
                        elem.coefficient.pop(x)
                        elem.coefficient_spin.pop(x)                        

                        # print(coef, coef_spin)
                        elem.coefficient.append(DENS1)                        
                        elem.coefficient_idx.append([coef_idx[0], coef_idx[2]])
                        elem.coefficient_spin.append(['+', '+'])
                        elem.delta.append([coef_idx[0], coef_idx[2]])
                        elem.delta_spin.append(['+', '+'])

                        # print(elem)
                        elem.coefficient.append(DENS1)                        
                        elem.coefficient_idx.append([coef_idx[1], coef_idx[3]])
                        elem.coefficient_spin.append(['-', '-'])
                        elem.delta.append([coef_idx[1], coef_idx[3]])
                        elem.delta_spin.append(['-', '-'])

                        res2 = res2 + arithmetic_string(elem)
                        # print('')
                        # print(elem)
                        
    return res2

def simplify_2rdm_mult_act(res):

    
    res2 = arithmetic_string()
    for x in range(0, len(res)):

        elem = res[x]
        
        if DENS2 not in elem.coefficient:
            res2 = res2 + arithmetic_string(elem)
        else:
            for x in range(0, len(elem.coefficient)):
                coef = deepcopy(elem.coefficient[x])
                coef_idx = deepcopy(elem.coefficient_idx[x])
                coef_spin = deepcopy(elem.coefficient_spin[x])
                        
                if coef == DENS2:
  #                  print('elem', elem)
                    if coef_spin[0]==coef_spin[1] and \
                       coef_spin[1]==coef_spin[2] and \
                       coef_spin[2]==coef_spin[3]:
                        if coef_spin[0] == '+':
                            elem.coefficient[x] = DENS2+"p"
                            res2 = res2 + arithmetic_string(elem)
                        elif coef_spin[0] == '-':
                            elem.coefficient[x] = DENS2+"m"
                            res2 = res2 + arithmetic_string(elem)

                    if coef_spin[0]==coef_spin[2] and \
                       coef_spin[0]!=coef_spin[1] and \
                       coef_spin[1]==coef_spin[3]:
                        if coef_spin[0] == '+':
                            elem.coefficient[x] = DENS2+"pm"
                            res2 = res2 + arithmetic_string(elem)
                        if coef_spin[0] == '-':
                            elem.coefficient[x] = DENS2+"mp"
                            res2 = res2 + arithmetic_string(elem)
#                    print('elempo', elem)
 #                   print()

    return res2

def simplify_3rdm_mult_act(res):

    res2 = arithmetic_string()
    for y in range(0, len(res)):
        elem = res[y]

        if DENS3 not in elem.coefficient:
            res2 = res2 + arithmetic_string(elem)
        else:
            for x in range(0, len(elem.coefficient)):
                coef = deepcopy(elem.coefficient[x])
                coef_idx = deepcopy(elem.coefficient_idx[x])
                coef_spin = deepcopy(elem.coefficient_spin[x])
                if coef == DENS3:
                    pl_cre = 0
                    mn_cre = 0
                    for i in range(0, 3):
                        if coef_spin[i] == "+":
                            pl_cre +=1
                        else:
                            mn_cre +=1
                    pl_ani = 0
                    mn_ani = 0
                    for i in range(3, 6):
                        if coef_spin[i] == "+":
                            pl_ani +=1
                        else:
                            mn_ani +=1
#                    print('DDDDens3', pl_cre, mn_cre, pl_ani, mn_ani, elem)

                if coef == DENS3:
                    if elem.num_factor != 0:
                        if pl_cre+pl_ani==4 and mn_cre+mn_ani ==2:
                            #print('Gm4p2m', elem, coef_idx, coef_spin)
                            sl1, sl2, s1 = perm_drag(coef_spin[0:3], coef_idx[0:3])
                            sll1, sll2, s2 = perm_drag(coef_spin[3:6], coef_idx[3:6])
                            new_coef = sl2 + sll2
                            new_spin = sl1 + sll1
                            elem.coefficient_idx[x] = new_coef
                            elem.coefficient_spin[x] = new_spin
                            elem.num_factor *= s1*s2
                            elem.coefficient[x] = DENS3+"pm"
                            #print('Gm4p2m', elem)
                            #print()
                        elif pl_cre+pl_ani==2 and mn_cre+mn_ani ==4:
                            #print('Gm4p2m zupelnie przed', elem, coef_idx, coef_spin)
                            for mini_spin_list in elem.coefficient_spin:
                                #print(elem, mini_spin_list, 'mini_spin_list')
                                for y in range(0, len(mini_spin_list)):
                                    if mini_spin_list[y] == '+':
                                        mini_spin_list[y] = '-'
                                    elif mini_spin_list[y] == '-':
                                        mini_spin_list[y] = '+'


                            for delspin in elem.delta_spin:
                                #print(delspin, 'delspin')
                                for y in range(0, len(delspin)):
                                    if delspin[y] == '+':
                                        delspin[y] = '-'
                                    elif delspin[y] == '-':
                                        delspin[y] = '+'

                            #print('Gm4p2m przed', elem, coef_idx, coef_spin)
                            coef = deepcopy(elem.coefficient[x])
                            coef_idx = deepcopy(elem.coefficient_idx[x])
                            coef_spin = deepcopy(elem.coefficient_spin[x])
                            
                            sl1, sl2, s1 = perm_drag(coef_spin[0:3], coef_idx[0:3])
                            sll1, sll2, s2 = perm_drag(coef_spin[3:6], coef_idx[3:6])
                            new_coef = sl2 + sll2
                            new_spin = sl1 + sll1
                            elem.coefficient_idx[x] = new_coef
                            elem.coefficient_spin[x] = new_spin
                            elem.num_factor *= s1*s2
                            elem.coefficient[x] = DENS3+"pm"
                            #print('Gm4p2m ppppo', elem)
                            #print()                                        

                        elif pl_cre+pl_ani==6 and mn_cre+mn_ani ==0:
                            elem.coefficient = [coef if coef != DENS3 else DENS3P for coef in elem.coefficient]
                        elif pl_cre+pl_ani==0 and mn_cre+mn_ani ==6:
                            elem.coefficient = [coef if coef != DENS3 else DENS3M for coef in elem.coefficient]

                    else:
                        if elem.num_factor != 0:
                            print('DENS3', elem, pl_cre+pl_ani, mn_cre+mn_ani)
            if (pl_cre == pl_ani and mn_cre == mn_ani):

                res2 = res2 + arithmetic_string(elem)

    return res2

def perm_sign(lst, lst2):

    sorted_lst = sorted(lst)
    perm = [sorted_lst.index(x) for x in lst]
    inversions = sum(1 for i in range(len(perm)) 
                    for j in range(i + 1, len(perm)) 
                    if perm[i] > perm[j])
    
    sign = 1 if inversions % 2 == 0 else -1
    return sorted_lst, sign


def perm_drag(lst1, lst2):

    perm_indices = sorted(range(len(lst1)), key=lambda k: lst1[k])
    inversions = sum(1 for i in range(len(perm_indices)) 
                    for j in range(i + 1, len(perm_indices)) 
                    if perm_indices[i] > perm_indices[j])
    sign = 1 if inversions % 2 == 0 else -1

    sorted_lst1 = sorted(lst1)
    sorted_lst2 = [lst2[i] for i in perm_indices]
    
    return sorted_lst1, sorted_lst2, sign


def simplify_4rdm_mult_act(res):

    res2 = arithmetic_string()
    for y in range(0, len(res)):
        elem = res[y]

        if DENS4 not in elem.coefficient:
            res2 = res2 + arithmetic_string(elem)
        else:
            for x in range(0, len(elem.coefficient)):
                coef = deepcopy(elem.coefficient[x])
                coef_idx = deepcopy(elem.coefficient_idx[x])
                coef_spin = deepcopy(elem.coefficient_spin[x])
                if coef == DENS4:
                    pl_cre = 0
                    mn_cre = 0
                    for i in range(0, 4):
                        if coef_spin[i] == "+":
                            pl_cre +=1
                        else:
                            mn_cre +=1
                    pl_ani = 0
                    mn_ani = 0
                    for i in range(4, 8):
                        if coef_spin[i] == "+":
                            pl_ani +=1
                        else:
                            mn_ani +=1
#                    print('DDDDens4', pl_cre, mn_cre, pl_ani, mn_ani, elem)
                if coef == DENS4:
                    if elem.num_factor != 0:
                        if pl_cre+pl_ani==6 and mn_cre+mn_ani ==2:
                            sl1, sl2, s1 = perm_drag(coef_spin[0:4], coef_idx[0:4])
                            sll1, sll2, s2 = perm_drag(coef_spin[4:8], coef_idx[4:8])
                            new_coef = sl2 + sll2
                            new_spin = sl1 + sll1
                            elem.coefficient_idx[x] = new_coef
                            elem.coefficient_spin[x] = new_spin
                            elem.num_factor *= s1*s2
                            elem.coefficient[x] = DENS4+"ppm"
                                                        
                            #print('Gm4p2m', elem)
                            #print()
                        elif pl_cre+pl_ani==2 and mn_cre+mn_ani ==6:
                            #print('Gm4p2m zupelnie przed', elem, coef_idx, coef_spin)
                            for mini_spin_list in elem.coefficient_spin:
                                #print(elem, mini_spin_list, 'mini_spin_list')
                                for y in range(0, len(mini_spin_list)):
                                    if mini_spin_list[y] == '+':
                                        mini_spin_list[y] = '-'
                                    elif mini_spin_list[y] == '-':
                                        mini_spin_list[y] = '+'

                            for delspin in elem.delta_spin:
                                #print(delspin, 'delspin')
                                for y in range(0, len(delspin)):
                                    if delspin[y] == '+':
                                        delspin[y] = '-'
                                    elif delspin[y] == '-':
                                        delspin[y] = '+'

                            #print('Gm4p2m przed', elem, coef_idx, coef_spin)
                            coef = deepcopy(elem.coefficient[x])
                            coef_idx = deepcopy(elem.coefficient_idx[x])
                            coef_spin = deepcopy(elem.coefficient_spin[x])
                            
                            sl1, sl2, s1 = perm_drag(coef_spin[0:4], coef_idx[0:4])
                            sll1, sll2, s2 = perm_drag(coef_spin[4:8], coef_idx[4:8])
                            new_coef = sl2 + sll2
                            new_spin = sl1 + sll1
                            elem.coefficient_idx[x] = new_coef
                            elem.coefficient_spin[x] = new_spin
                            elem.num_factor *= s1*s2
                            elem.coefficient[x] = DENS4+"ppm"
                            #print('Gm4p2m ppppo', elem)
                            #print()
                        elif pl_cre+pl_ani==4 and mn_cre+mn_ani ==4:
                            sl1, sl2, s1 = perm_drag(coef_spin[0:4], coef_idx[0:4])
                            sll1, sll2, s2 = perm_drag(coef_spin[4:8], coef_idx[4:8])
                            new_coef = sl2 + sll2
                            new_spin = sl1 + sll1
                            elem.coefficient_idx[x] = new_coef
                            elem.coefficient_spin[x] = new_spin
                            elem.num_factor *= s1*s2
                            elem.coefficient[x] = DENS4+"pm"

                        elif pl_cre+pl_ani==8 and mn_cre+mn_ani ==0:
                            elem.coefficient = [coef if coef != DENS4 else DENS4P for coef in elem.coefficient]
                        elif pl_cre+pl_ani==0 and mn_cre+mn_ani ==8:
                            elem.coefficient = [coef if coef != DENS4 else DENS4M for coef in elem.coefficient]


                    else:
                        if elem.num_factor != 0:
                            print('DENS4-whooo', elem, pl_cre+pl_ani, mn_cre+mn_ani)
#            if elem.num_factor != 0:
#                print('DENS4', elem, pl_cre+pl_ani, mn_cre+mn_ani)                            
            if (pl_cre == pl_ani and mn_cre == mn_ani):
                if elem.num_factor != 0:
                    res2 = res2 + arithmetic_string(elem)


    return res2

# def simplify_4rdm_mult_act(res):

#     res2 = arithmetic_string()
#     for x in range(0, len(res)):
#         elem = res[x]

#         if DENS4 not in elem.coefficient:
#             res2 = res2 + arithmetic_string(elem)
#         else:
#             for x in range(0, len(elem.coefficient)):
#                 coef = deepcopy(elem.coefficient[x])
#                 coef_idx = deepcopy(elem.coefficient_idx[x])
#                 coef_spin = deepcopy(elem.coefficient_spin[x])
#                 if coef == DENS4:
#                     pl_cre = 0
#                     mn_cre = 0
#                     for i in range(0, 4):
#                         if coef_spin[i] == "+":
#                             pl_cre +=1
#                         else:
#                             mn_cre +=1
#                     pl_ani = 0
#                     mn_ani = 0
#                     for i in range(4, 8):
#                         if coef_spin[i] == "+":
#                             pl_ani +=1
#                         else:
#                             mn_ani +=1

#             print('DENS4', elem, pl_cre+pl_ani, mn_cre+mn_ani)                            
#             if (pl_cre == pl_ani and mn_cre == mn_ani):
#                 res2 = res2 + arithmetic_string(elem)

#     return res2

def dirac_to_coulomb(res):

    res2 = arithmetic_string()
    for x in range(0, len(res)):
        elem = res[x]

        for x in range(0, len(elem.coefficient)):
            # print('elem', elem)
            coef = deepcopy(elem.coefficient[x])
            coef_idx = deepcopy(elem.coefficient_idx[x])            
            if coef == TWOEL_INT_DIRAC: 
                new_coef_idx = []
                new_coef_idx.append(coef_idx[0])
                new_coef_idx.append(coef_idx[2])
                new_coef_idx.append(coef_idx[1])
                new_coef_idx.append(coef_idx[3])
                elem.coefficient_idx[x] = new_coef_idx
                elem.coefficient[x] = TWOEL_INT
                # print('elempo', elem)
        res2 = res2 + arithmetic_string(elem)

    return res2



def remove_spin(res):

    for x in range(0, len(res)):
        res[x].coefficient_spin = []
        res[x].delta_spin = []

    return res

def rename_gm1(res):

    for x in res:
        indices_dens1 = [i for i, x in enumerate(x.coefficient) if x == DENS1]
        indices_dens1P = [i for i, x in enumerate(x.coefficient) if x == DENS1P]
        indices_dens1M = [i for i, x in enumerate(x.coefficient) if x == DENS1M]
        indices_densn = [i for i, x in enumerate(x.coefficient) if x == DENSN]
        for i in indices_dens1:
            if len(x.coefficient_idx[i]) >1:
                print(x, x.coefficient_idx[i])
                if x.coefficient_idx[i][0] ==  x.coefficient_idx[i][1]:
                    x.coefficient_idx[i] = [x.coefficient_idx[i][0]]
                else:
                    xcoef = sorted(x.coefficient_idx[i])
                    for j in x.delta:
                        js = sorted(j)
                        if js == xcoef:
                            x.coefficient_idx[i] = [xcoef[0]]
            x.coefficient[i] = DENSN1

        for i in indices_dens1P:
            if len(x.coefficient_idx[i]) >1:
                print(x, x.coefficient_idx[i])
                if x.coefficient_idx[i][0] ==  x.coefficient_idx[i][1]:
                    x.coefficient_idx[i] = [x.coefficient_idx[i][0]]
                else:
                    xcoef = sorted(x.coefficient_idx[i])
                    for j in x.delta:
                        js = sorted(j)
                        if js == xcoef:
                            x.coefficient_idx[i] = [xcoef[0]]
            x.coefficient[i] = DENSN1P

        for i in indices_dens1M:            
            if len(x.coefficient_idx[i]) >1:
                print(x, x.coefficient_idx[i])
                if x.coefficient_idx[i][0] ==  x.coefficient_idx[i][1]:
                    x.coefficient_idx[i] = [x.coefficient_idx[i][0]]
                else:
                    xcoef = sorted(x.coefficient_idx[i])
                    for j in x.delta:
                        js = sorted(j)
                        if js == xcoef:
                            x.coefficient_idx[i] = [xcoef[0]]
            x.coefficient[i] = DENSN1M


        for i in indices_densn:
            if len(x.coefficient_idx[i]) >1:
                print('la')
                if x.coefficient_idx[i][0] ==  x.coefficient_idx[i][1]:
                    x.coefficient_idx[i] = [x.coefficient_idx[i][0]]
                else:
                    xcoef = sorted(x.coefficient_idx[i])
                    for j in x.delta:
                        js = sorted(j)
                        if js == xcoef:
                            x.coefficient_idx[i] = [xcoef[0]]

            x.coefficient[i] = DENSN1


    return res
        
        
def rename_and_rescale_gamma(res):

    res2 = arithmetic_string()
    for x in range(0, len(res)):
        elem = res[x]
        # print('baddam', elem)
        elem2 = deepcopy(elem)
        elem2.coefficient_idx = []
        elem2.coefficient = []
        for x in range(0, len(elem.coefficient)):
            coef = deepcopy(elem.coefficient[x])
            coef_idx = deepcopy(elem.coefficient_idx[x])
            # print('ten', coef)
            if coef == DENS1:
                
                # 1/2 pochodz z zamiany gamma++ na gamma
                # 2 pochodzi z zamiany gamma na n
                elem.num_factor *= 0.5
                elem.num_factor *= 2.0

                elem2.coefficient.append(DENSN)
                elem2.coefficient_idx.append([coef_idx[0]])
                # print('')
                # print('mh', elem)
                # print('')
            else:
                elem2.coefficient.append(coef)
                elem2.coefficient_idx.append(coef_idx)

        # print('dodaje', elem2)
        res2 = res2 + arithmetic_string(elem2)


    return res2


#---------------------NEW APPROXIMATE GAMMA-------------
def approximate_rdm3(res, onedet = False):

    # Approximate all  3-RDM
    # in accordance with
    # Eq (40 a), Normal order by Mukherjee
    # Gamma3 = Lambda3  + sum(-1)^p Gm1 * Lambda2
    #          + det{Gm1Gm1Gm1}

    # Part A = Lambda3 ---> 0 due to approximation
    # Part B = sum(-1)^p Gm1 * Lambda2
    # Part C = det{Gm1Gm1Gm1}
    # if onedet = True, approximate by 1rdm only

    # res is arithmetic tring
    #pluszeron
    rdm_order = 3
    res_mega_big = []                        
    for i in range(0, len(res)):        
        
        place = []
        order = []        
        for k in range(0, len(res[i].coefficient)):
            if res[i].coefficient[k] == DENS3:                
                place.append(k)
                order.append(rdm_order)

        res_big = []                    
  #      print('placeorder', place, order)
        if len(place) > 0:
            lenpl = len(place)
            resxa = arithmetic_string(res[i])

            for k in range(0, len(place)):
                result_big = arithmetic_string()

                for kk in range(0, len(resxa)):

                    if not onedet:
 #                       print('part B')
                        split_first = False
                        partB = generate_partition(deepcopy(resxa[kk]), place[k], rdm_order, 1, 2, split_first)
                    else:
                        partB = arithmetic_string()
#                    print('part C')
                    partC = generate_1rdm_prod(deepcopy(resxa[kk]), place[k], rdm_order)
                for x in partB:
#                    print('dodajeB', x)
                    res_big.append(deepcopy(x))                           
                for x in partC:
#                    print('dodajeC', x)
                    res_big.append(deepcopy(x))
            resxa = deepcopy(res_big)
        else:
            res_big.append(res[i])

        for x in res_big:
            res_mega_big.append(x)


    return res_mega_big


def approximate_rdm4(res, onedet = False):

    # Approximate all  4-RDM
    # in accordance with
    # Eq (40 b), Normal order by Mukherjee
    # Gamma4 = Lambda4  + sum(-1)^p Gm1 * Lambda3
    #          + 0.5 sum(-1)^p Lambda2 * Lambda2
    #          + sum (-1)^p Gm1 Gm1 Lambda2
    #          + det{Gm1Gm1Gm1}

    # Part A = Lambda4 ---> 0 due to approximation
    # Part B = sum(-1)^p Gm1 * Lambda3 ---> 0 due to approximation
    # Part C = 0.5 sum(-1)^p Lambda2 * Lambda2
    # Part D = sum (-1)^p Gm1 Gm1 Lambda2
    # Part e = det{Gm1Gm1Gm1}
    # if onedet = True, approximate by 1rdm only

    # res is arithmetic tring
    #pluszeron
    rdm_order = 4
    res_mega_big = []                        
    for i in range(0, len(res)):        
        
        place = []
        order = []        
        for k in range(0, len(res[i].coefficient)):
            if res[i].coefficient[k] == DENS4:                
                place.append(k)
                order.append(rdm_order)

        res_big = []                    
  #      print('placeorder', place, order)
        if len(place) > 0:
            lenpl = len(place)
            resxa = arithmetic_string(res[i])

            for k in range(0, len(place)):
                result_big = arithmetic_string()

                for kk in range(0, len(resxa)):

                    if not onedet:
                        split_first = False
                        partC = generate_partition(deepcopy(resxa[kk]), place[k], rdm_order, 2, 2, split_first, 0.5)
                    else:
                        partC = arithmetic_string()
                    if not onedet:
                        split_first = True
                        partD = generate_partition(deepcopy(resxa[kk]), place[k], rdm_order, 2, 2, split_first)
                    else:
                        partD = arithmetic_string()
                        
                    partE = generate_1rdm_prod(deepcopy(resxa[kk]), place[k], rdm_order)
                for x in partC:
                    res_big.append(deepcopy(x))                           
                for x in partD:
                    res_big.append(deepcopy(x))
                for x in partE:
                    res_big.append(deepcopy(x))

            resxa = deepcopy(res_big)
        else:
            res_big.append(res[i])

        for x in res_big:
            res_mega_big.append(x)


    return res_mega_big

#-------------------------------------------------------

def approximate_rdm3_by_rdm12(res):

    # Approximate all 3-RDM 
    # by all possible 2-RDM and 1-RDM products

    res_big = []
    zz = 1
    for i in range(0, len(res)):
        print('lala', res[i])
        place = []
        order = []
        for k in range(0, len(res[i].coefficient)):
            print('-->', res[i].coefficient[k])
            if res[i].coefficient[k] == DENS3:
                place.append(k)
                order.append(3)
        # print(place, order)
        if len(place) > 0:
#            print('z takiego', i, res[i])
            newres = generate_12rdm_prod(res[i], place, order)
#            print('mam taki')
            for x in newres:
 #               print('qq', zz, res[i])
                zz +=1
                res_big.append(deepcopy(x))
                # print(x)
#            print('')
#            sys.exit(0)
        else:
#            print('qq', zz, res[i])
            zz +=1
            res_big.append(res[i])

    return res_big

def approximate_rdm34_by_rdm1(res, onlyFour = False):

    # Approximate all 3-RDM and 4-RDM
    # by all possible 1-RDM products
    # if onlyFour = True, approximate only 4-RDM.
    # 3-RDM are approximated in another way

    res_big = []
    for i in range(0, len(res)):
        # print('lala', res[i])
        place = []
        order = []
        for k in range(0, len(res[i].coefficient)):
            if not onlyFour:
                if res[i].coefficient[k] == DENS3:
                    place.append(k)
                    order.append(3)
            if res[i].coefficient[k] == DENS4:
                place.append(k)
                order.append(4)
        # print(place, order)
        if len(place) > 0:
            newres = generate_1rdm_prod(res[i], place, order)
            for x in newres:
                res_big.append(deepcopy(x))
        else:
            res_big.append(res[i])

    return res_big

def approximate_rdm2_by_rdm1_when_occ(res, nospin=False, onedet=False):

    # Approximate all 2-RDMs that does not have
    # all indices active
    # if onedet = True - approximate everything
    #             there are no active orbitals
    # if onedet = False - aproximate only those
    #             that have not active orbitals

    res_big = []
    zz = 1
    for i in range(0, len(res)):
        place = []
        order = []
        for k in range(0, len(res[i].coefficient)):
            if res[i].coefficient[k] == DENS2:
                coef_idx = res[i].coefficient_idx[k]
                if onedet:
                    place.append(k)
                    order.append(3)
                else:
                    do_approx = check_if_all_active(coef_idx)
                    if do_approx:
                        place.append(k)
                        order.append(3)
        # print(place, order)
        if len(place) > 0:
#            print('z takiego2', i, res[i])
            if nospin:
                newres = generate_1rdm_prod_for_2RDM_nospin(res[i], place, order)
            else:
                newres = generate_1rdm_prod_for_2RDM(res[i], place, order)
#            print('mam taki2')
            for x in newres:
#                print(x)
 #               print('zz', zz, deepcopy(x))
                zz+=1
                res_big.append(deepcopy(x))
#            print('')
        else:
  #          print('zz', zz, res[i])
            zz+=1
            res_big.append(res[i])

    return res_big


def cum2_to_rdm2(res, nospin=False, onedet=False):

    # express cum2 in rdm2
    # if onedet = True cum2 = 0

    res_big = []

    for i in range(0, len(res)):
        place = []
        order = []
        # print(i, 'robie', res[i])
        for k in range(0, len(res[i].coefficient)):
            if res[i].coefficient[k] == CUM2:
                coef_idx = res[i].coefficient_idx[k]
                if onedet:
                    res[i].num_factor = 0
                else:

                    place.append(k)

        if len(place) > 0:
            if nospin:
                newres = rdm1_rdm2_prod_nospin(res[i], place)
            else:
                print('bardzo mi przykro ale musisz doprogramowac te wersje')
            # print('')
            for x in newres:
                # print('dodaje', x)
                res_big.append(deepcopy(x))
            # print('')
                            
        else:
            res_big.append(res[i])

        # print('teraz mma')
        # print('')
        # for x in res_big:
        #     print(x)
        # print('')

    
        
    return res_big

def check_if_all_active(coef_idx):

    for x in coef_idx:
        if x in occupied:
            return False
        if x in virtual:
            return False

    return True

def generate_1rdm_prod(resx, place, order):

    # resx - one ugg which contain RDM of
    # order = order as the place-th coefficient
    # return arithmetic string of terms
    # with all possible 1RDM permutations
    # place - list of indices with higher RDMs
    # order - list of orders of RDMs in places
    
    # print('bede robic', resx, place, order)

    result = []
    
    lst = []
    for i in range(0, order):
        lst.append(i)

    tup = permutations(lst)

    plist = []
    # print('halo')
    for x in tup:
        plist.append(list(x))
        # print(list(x))

    result = arithmetic_string()
        
    n  = len(plist)

    for i in range(0, n):
        result.append(deepcopy(resx))

    perm_idx = 0
    for i in range(0, n):

        idx_ani = deepcopy(result[i].coefficient_idx[place][0:order])
        idx_cre = deepcopy(result[i].coefficient_idx[place][order:])

        new_idx_ani = []
        for anidx in range(0, len(idx_ani)):
            # print('anidx', anidx)
            whichno = plist[perm_idx][anidx]
            # print(whichno, plist[perm_idx], anidx)
            new_idx_ani.append(idx_ani[whichno])
        # print('i jestem')
        result[i].coefficient.pop(place)
        result[i].coefficient_idx.pop(place)

        for t in range(0, order):
            result[i].coefficient.append(DENS1)
            result[i].coefficient_idx.append([new_idx_ani[t], idx_cre[t]])

        # print('***', plist[0], plist[perm_idx])
        nswaps = nswaps_in_two_lists(deepcopy(plist[0]), plist[perm_idx])
        # print('nswaps', nswaps)
        if (nswaps%2) != 0:
            result[i].num_factor *= -1.0
#            print('minus')
        # else:
        #     print('plus')

        perm_idx += 1

    return result

def generate_partition(resx, place, order, kl, kr, split_first=False, numf = 1.0):
    
    
    # resx - one ugg which contain RDM of
    # order = order as the place-th coefficient
    # return arithmetic string of terms
    # with all possible 1RDM permutations
    # place - list of indices with higher RDMs
    # order - list of orders of RDMs in places

    
    # print('bede robic', resx, place, order)

    if kl == 1:
        namel = DENS1
    else:
        namel = 'Cm' + str(kl)
    namer = 'Cm' + str(kr)

    if split_first:
        namel = DENS1
        namem = DENS1

    idx_ani = deepcopy(resx.coefficient_idx[place][0:order])
    idx_cre = deepcopy(resx.coefficient_idx[place][order:])

    idx_ani_tup = combinations(idx_ani, kl)
    idx_cre_tup = combinations(idx_cre, kl)
    idx_ani_comb = []
    idx_cre_comb = []
    
    for x in idx_ani_tup:
        temp = list(x)
        idx_ani_comb.append(temp[0:kl])

    for x in idx_cre_tup:
        temp = list(x)
        idx_cre_comb.append(temp[0:kl])


    lstl_cre = []
    lstl_ani = []
    lstp_cre = []
    lstp_ani = []
    cre_perm = []
    ani_perm = []
    
    pc_idx = 0
    for jc_i in range(0, len(idx_cre_comb)):
        for ja_i in range(0, len(idx_ani_comb)):
            pc_idx += 1
            jc = idx_cre_comb[jc_i]
            ja = idx_ani_comb[ja_i]
            lstl_cre.append(jc)
            lstl_ani.append(ja)
            
            lstp_cre_mini = []
            lstp_ani_mini = []
            cre_perm_mini = []
            ani_perm_mini = []
            for kk in jc:
                cre_perm_mini.append(idx_cre.index(kk))
            for kk in ja:
                ani_perm_mini.append(idx_ani.index(kk))

            for jcr_i in range(0, len(idx_cre)):
                if idx_cre[jcr_i] not in jc:
                    jcr = idx_cre[jcr_i]
                    lstp_cre_mini.append(jcr)
                    cre_perm_mini.append(jcr_i)

            for jar_i in range(0, len(idx_ani)):
                if idx_ani[jar_i] not in ja:
                    jar = idx_ani[jar_i]
                    lstp_ani_mini.append(jar)
                    ani_perm_mini.append(jar_i)
            lstp_cre.append(lstp_cre_mini)
            lstp_ani.append(lstp_ani_mini)
            cre_perm.append(cre_perm_mini)
            ani_perm.append(ani_perm_mini)
    

    sign_list = []
    result = arithmetic_string()

    # print('resx', resx)
    for i in range(0, pc_idx):
        result.append(deepcopy(resx))

        nswaps_cre = nswaps_in_two_lists(deepcopy(cre_perm[0]), deepcopy(cre_perm[i]))
        nswaps_ani = nswaps_in_two_lists(deepcopy(ani_perm[0]), deepcopy(ani_perm[i]))
        cres = 1.0
        anis = 1.0
        if (nswaps_cre%2) != 0:
            cres =  -1.0
        if (nswaps_ani%2) != 0:
            anis =  -1.0
            
        sign_list.append(cres * anis)

        result[i].coefficient.pop(place)
        result[i].coefficient_idx.pop(place)
        if split_first:
 #           print('taksplit')
            result[i].coefficient.append(namel)
            result[i].coefficient_idx.append([lstl_ani[i][0], lstl_cre[i][0]])
#            print(result[i].coefficient_idx)
            result[i].coefficient.append(namem)
            result[i].coefficient_idx.append([lstl_ani[i][1], lstl_cre[i][1]])
        else:
            result[i].coefficient.append(namel)
            result[i].coefficient_idx.append(lstl_ani[i]+ lstl_cre[i])
        result[i].coefficient.append(namer)
        result[i].coefficient_idx.append(lstp_ani[i]+ lstp_cre[i])
        result[i].num_factor *= sign_list[i]


    # print('a to wynik')
    i = 1
    for x in range(0, len(result)):
        result[x].num_factor *= numf
        # print(i, result[x])
        i+= 1
    # print('')

                    

    return result




def generate_1rdm_prod_for_2RDM(resx, place, order):

    # resx - one ugg which contain RDM of
    # order = order as the place-th coefficient
    # return arithmetic string of terms
    # with all possible 1RDM permutations
    # place - list of indices with higher RDMs
    # order - list of orders of RDMs in places
    
    # print('bede robic', resx, place, order)

    result = []
    start = 0
    lenpl = len(place)
    start = 0

    for j in range(0, lenpl):

        result = arithmetic_string()

        this_2RDM_idx = place[j]
        this_2RDM_idx_list = deepcopy(resx.coefficient_idx[this_2RDM_idx] )
        this_RDM_spin = deepcopy(resx.coefficient_spin[this_2RDM_idx])
        spin_type = gam2_spin_type(this_RDM_spin)


        if spin_type == '4p' or spin_type == '4m':
            for i in range(0, 2):
                result.append(deepcopy(resx))

            
            this_fragment = cas()
            this_fragment.coefficient.append(DENS1)
            this_fragment.coefficient.append(DENS1)
            l1 = [this_2RDM_idx_list[0], this_2RDM_idx_list[2]]
            l2 = [this_2RDM_idx_list[1], this_2RDM_idx_list[3]]
            this_fragment.coefficient_idx.append(l1)
            this_fragment.coefficient_idx.append(l2)
            if spin_type == '4p':
                this_fragment.coefficient_spin.append(['+', '+'])
                this_fragment.coefficient_spin.append(['+', '+'])
                this_fragment.delta_spin.append(['+', '+'])
                this_fragment.delta_spin.append(['+', '+'])
            elif spin_type == '4m':
                this_fragment.coefficient_spin.append(['-', '-'])
                this_fragment.coefficient_spin.append(['-', '-'])
                this_fragment.delta_spin.append(['-', '-'])
                this_fragment.delta_spin.append(['-', '-'])


            result[0].coefficient.pop(this_2RDM_idx)
            result[0].coefficient_idx.pop(this_2RDM_idx)
            result[0].coefficient_spin.pop(this_2RDM_idx)

            #print('this-a', this_fragment)
        
            result[0] = result[0].fromright(this_fragment)
            #print('res[0]', result[0])


            this_fragment = cas()
            this_fragment.coefficient.append(DENS1)
            this_fragment.coefficient.append(DENS1)
            l1 = [this_2RDM_idx_list[0], this_2RDM_idx_list[3]]
            l2 = [this_2RDM_idx_list[1], this_2RDM_idx_list[2]]
            this_fragment.coefficient_idx.append(l1)
            this_fragment.coefficient_idx.append(l2)
            if spin_type == '4p':
                this_fragment.coefficient_spin.append(['+', '+'])
                this_fragment.coefficient_spin.append(['+', '+'])
                this_fragment.delta_spin.append(['+', '+'])
                this_fragment.delta_spin.append(['+', '+'])                            
            elif spin_type == '4m':
                this_fragment.coefficient_spin.append(['-', '-'])
                this_fragment.coefficient_spin.append(['-', '-'])
                this_fragment.delta_spin.append(['-', '-'])
                this_fragment.delta_spin.append(['-', '-'])

            result[1].num_factor *= -1
            result[1].coefficient.pop(this_2RDM_idx)
            result[1].coefficient_idx.pop(this_2RDM_idx)
            result[1].coefficient_spin.pop(this_2RDM_idx)


            #print('this-b', this_fragment)
            
            result[1] = result[1].fromright(this_fragment)
            #print('res[0]', result[1])

        elif spin_type == 'mp' or spin_type == 'pm':
            for i in range(0, 1):
                result.append(deepcopy(resx))

            
            this_fragment = cas()
            this_fragment.coefficient.append(DENS1)
            this_fragment.coefficient.append(DENS1)
            l1 = [this_2RDM_idx_list[0], this_2RDM_idx_list[2]]
            l2 = [this_2RDM_idx_list[1], this_2RDM_idx_list[3]]
            this_fragment.coefficient_idx.append(l1)
            this_fragment.coefficient_idx.append(l2)
            if spin_type == 'pm':
                this_fragment.coefficient_spin.append(['+', '+'])
                this_fragment.coefficient_spin.append(['-', '-'])
                this_fragment.delta_spin.append(['+', '+'])
                this_fragment.delta_spin.append(['-', '-'])
            elif spin_type == 'mp':
                this_fragment.coefficient_spin.append(['-', '-'])
                this_fragment.coefficient_spin.append(['+', '+'])
                this_fragment.delta_spin.append(['-', '-'])
                this_fragment.delta_spin.append(['+', '+'])

            result[0].coefficient.pop(this_2RDM_idx)
            result[0].coefficient_idx.pop(this_2RDM_idx)
            result[0].coefficient_spin.pop(this_2RDM_idx)
        
            #print('this-c', this_fragment)

            result[0] = result[0].fromright(this_fragment)
            #print('res[0]', result[0])

    return result


def generate_1rdm_prod_for_2RDM_nospin(resx, place, order):

    # resx - one ugg which contain RDM of
    # order = order as the place-th coefficient
    # return arithmetic string of terms
    # with all possible 1RDM permutations
    # place - list of indices with higher RDMs
    # order - list of orders of RDMs in places
    
    # print('bede robic', resx, place, order)

    result = []
    start = 0
    lenpl = len(place)
    start = 0

    for j in range(0, lenpl):

        result = arithmetic_string()

        this_2RDM_idx = place[j]
        this_2RDM_idx_list = deepcopy(resx.coefficient_idx[this_2RDM_idx] )

        for i in range(0, 2):
            result.append(deepcopy(resx))
            
        this_fragment = cas()
        this_fragment.coefficient.append(DENS1)
        this_fragment.coefficient.append(DENS1)
        l1 = [this_2RDM_idx_list[0], this_2RDM_idx_list[2]]
        l2 = [this_2RDM_idx_list[1], this_2RDM_idx_list[3]]
        this_fragment.coefficient_idx.append(l1)
        this_fragment.coefficient_idx.append(l2)

        result[0].coefficient.pop(this_2RDM_idx)
        result[0].coefficient_idx.pop(this_2RDM_idx)
        
        result[0] = result[0].fromright(this_fragment)

        this_fragment = cas()
        this_fragment.coefficient.append(DENS1)
        this_fragment.coefficient.append(DENS1)
        l1 = [this_2RDM_idx_list[0], this_2RDM_idx_list[3]]
        l2 = [this_2RDM_idx_list[1], this_2RDM_idx_list[2]]
        this_fragment.coefficient_idx.append(l1)
        this_fragment.coefficient_idx.append(l2)

        result[1].num_factor *= -1
        result[1].coefficient.pop(this_2RDM_idx)
        result[1].coefficient_idx.pop(this_2RDM_idx)
        
        result[1] = result[1].fromright(this_fragment)


    return result


def rdm1_rdm2_prod_nospin(resx, place):

    # resx - one ugg which contain RDM of
    # order = order as the place-th coefficient
    # return arithmetic string of terms
    # with all possible 1RDM permutations
    # place - list of indices with higher RDMs
    # order - list of orders of RDMs in places
    
    # print('bede robic', resx, place, order)

    lenpl = len(place)
    resxa = arithmetic_string(resx)

    for j in range(0, lenpl):

        result_big = arithmetic_string()

        this_cum2_idx = place[j]

        for k in range(0, len(resxa)):
            this_cum2_idx_list = deepcopy(resxa[k].coefficient_idx[this_cum2_idx] )

            result = arithmetic_string()

            for i in range(0, 3):
                result.append(deepcopy(resxa[k]))
        
            this_fragment = cas()
            this_fragment.coefficient.append(DENS1)
            this_fragment.coefficient.append(DENS1)
            l1 = [this_cum2_idx_list[0], this_cum2_idx_list[2]]
            l2 = [this_cum2_idx_list[1], this_cum2_idx_list[3]]
            this_fragment.coefficient_idx.append(l1)
            this_fragment.coefficient_idx.append(l2)
            
       #     print('')
        #    print(this_fragment)
            # result[0].coefficient.pop(this_cum2_idx)
            # result[0].coefficient_idx.pop(this_cum2_idx)

         #   print('')            
          #  print('res0', result[0])
            result[0] = result[0].fromright(this_fragment)
           # print('res0', result[0])
            #print('')
            #        sys.exit(0)

            this_fragment = cas()
            this_fragment.coefficient.append(DENS1)
            this_fragment.coefficient.append(DENS1)
            l1 = [this_cum2_idx_list[0], this_cum2_idx_list[3]]
            l2 = [this_cum2_idx_list[1], this_cum2_idx_list[2]]
            this_fragment.coefficient_idx.append(l1)
            this_fragment.coefficient_idx.append(l2)
            
            result[1].num_factor *= -1
            # result[1].coefficient.pop(this_cum2_idx)
            # result[1].coefficient_idx.pop(this_cum2_idx)
            
            result[1] = result[1].fromright(this_fragment)
    #        print('res1', result[1])
            #        sys.exit(0)


            this_fragment = cas()
            this_fragment.coefficient.append(DENS2)
            this_fragment.coefficient_idx.append(this_cum2_idx_list)
            result[2] = result[2].fromright(this_fragment)
                       
#            result[2].coefficient[this_cum2_idx] = DENS2

            # result[2] = result[0].fromright(this_fragment)
     #       print('res2', result[2])
      #      print('')

            for i in range(0, 3):
             #   print('do big dodaje', result[i])
                result_big.append(deepcopy(result[i]))
        resxa = deepcopy(result_big)



    for k in range(0, lenpl):
        this_cum2_idx = place[k]
        #print('kkk?', k)
        for j in range(0, len(result_big)):
            result_big[j].coefficient.pop(this_cum2_idx)
            result_big[j].coefficient_idx.pop(this_cum2_idx)
        for j in range(k+1, len(place)):
            place[j] -=1

    # print('wynika')
    # for x in result_big:
    #     print(x)


    return result_big



def generate_12rdm_prod(resx, place, order):

    # resx - one ugg which contain RDM of
    # order = order at the place-th coefficient
    # return arithmetic string of terms
    # with the approximation of 3-RDM by 1 and 2-RDM
    # from Kutzelnig paper Conditions on the cumulants
    # eq. 3.3 with 3-body cumulant = 0
    # place - list of indices with higher RDMs
    # order - list of orders of RDMs in places
    
    # print('bede robic rdm33', resx, place, order)

    result = []
    start = 0
    lenpl = len(place)
    start = 0
    for j in range(0, lenpl):


        lst_12 = generate_lst_12()
        lst_1 = generate_lst_1()
    
    
        result = arithmetic_string()

        n_12  = len(lst_12)
        n_1  = len(lst_1)

        for i in range(0, n_12):
            result.append(deepcopy(resx))
        for i in range(0, n_1):
            result.append(deepcopy(resx))

        # print('zaczynam od skopiowania')
        # for x in result:
        #     print(x)


        for i in range(0, n_12):

            this_3RDM_idx = place[j]
            this_3RDM_idx_list = result[i].coefficient_idx[this_3RDM_idx]

            this_fragment = cas()
            this_fragment.coefficient.append(DENS1)
            this_fragment.coefficient.append(CUM2)
            i1 = lst_12[i][0]
            i2 = lst_12[i][1]
            # print(i1, i2, this_3RDM_idx_list)
            l1 = [this_3RDM_idx_list[i1], this_3RDM_idx_list[i2]]
            this_fragment.coefficient_idx.append(l1)
            i1 = lst_12[i][2]
            i2 = lst_12[i][3]
            i3 = lst_12[i][4]                
            i4 = lst_12[i][5]
            l2 = [this_3RDM_idx_list[i1], this_3RDM_idx_list[i2], this_3RDM_idx_list[i3], this_3RDM_idx_list[i4]]
            this_fragment.coefficient_idx.append(l2)
            i1 = lst_12[i][6]
#            print('iii', int(i1))
            this_fragment.num_factor = int(i1)
            # print('dys', this_fragment)

            # usun Gamma_3
            result[i].coefficient.pop(this_3RDM_idx)
            result[i].coefficient_idx.pop(this_3RDM_idx)

            result[i] = result[i].fromright(this_fragment)

        k = 0
        for i in range(n_12, n_12+n_1):


            this_3RDM_idx = place[j]
            this_3RDM_idx_list = result[i].coefficient_idx[this_3RDM_idx]

            this_fragment = cas()
            this_fragment.coefficient.append(DENS1)
            this_fragment.coefficient.append(DENS1)
            this_fragment.coefficient.append(DENS1)
            i1 = lst_1[k][0]
            i2 = lst_1[k][1]
            l1 = [this_3RDM_idx_list[i1], this_3RDM_idx_list[i2]]
            this_fragment.coefficient_idx.append(l1)
            i1 = lst_1[k][2]
            i2 = lst_1[k][3]
            l2 = [this_3RDM_idx_list[i1], this_3RDM_idx_list[i2]]
            this_fragment.coefficient_idx.append(l2)
            i1 = lst_1[k][4]
            i2 = lst_1[k][5]
            l3 = [this_3RDM_idx_list[i1], this_3RDM_idx_list[i2]]
            this_fragment.coefficient_idx.append(l3)
            i1 = lst_1[k][6]

            this_fragment.num_factor = int(i1)
            print('this_frag')
            print(this_fragment)


            # usun Gamma_3
            result[i].coefficient.pop(this_3RDM_idx)
            result[i].coefficient_idx.pop(this_3RDM_idx)

            result[i] = result[i].fromright(this_fragment)
            k+=1


    # print('a to wynik')
    # for x in result:
    #     print(x)
    # print('')

    return result

def generate_lst_12():

    lst = []
    lst.append([0,3,1,2,4,5, 1])
    lst.append([1,4,0,2,3,5, 1])
    lst.append([2,5,0,1,3,4, 1])
    lst.append([1,3,0,2,4,5, -1])
    lst.append([2,4,0,1,3,5, -1])
    lst.append([0,5,2,1,3,4, -1])
    lst.append([2,3,1,0,4,5,-1])
    lst.append([0,4,1,2,3,5, -1])
    lst.append([1,5,0,2,3,4, -1])

    return lst

def generate_lst_1():

    lst = []
    # lst.append([0,3,1,4,2,5, -2])
    # lst.append([0,3,2,4,1,5, 2])
    # lst.append([1,3,0,4,2,5, 2])
    # lst.append([1,3,2,4,0,5, -2])
    # lst.append([2,3,0,4,1,5, -2])
    # lst.append([2,3,1,4,0,5, 2])
    lst.append([0,3,1,4,2,5, 1])
    lst.append([0,3,2,4,1,5, -1])
    lst.append([1,3,0,4,2,5, -1])
    lst.append([1,3,2,4,0,5, 1])
    lst.append([2,3,0,4,1,5, 1])
    lst.append([2,3,1,4,0,5, -1])

    return lst

def nswaps_in_two_lists(nums, lst2):

    #compute_n_of_swaps_in_sort

    table = {}
    # print('nums', nums)
    # print('lst', lst2)
    for i in range(0, len(nums)):
        table[nums[i]] = i
    # print('table', table)
    swaps = 0
    for i in range(len(nums)):
        n = nums[i]
        s_n = lst2[i]
        s_i = table[s_n]

        if s_n != n:
            swaps += 1
            nums[s_i] = n
            nums[i] = s_n
            table[n] = s_i
            table[s_n] = i
    #     print(i, nums, swaps)
    # print('=========================')
    return swaps

def are_spins_equal(lst):

    p = 0
    m = 0

    for x in lst:
        if x == '+':
            p+=1
        elif x == '-':
            m +=1
#    print(p, m)
    if p!=m:
        if p==0 or m==0:
            return True
        else:
            return False
    else:
        return True



def read_cas_from_str(line, spin=True):

    strike1 = False
    rs = cas()
    print(line)
    z =  re.split('}|{|_|\\\\|\^|<|>', line)
    s = []
    for x in z:
        if x != '':
            s.append(x)
    print(s)

    if "sum" not in s:
        if s[0][0] == '+':

            rs.num_factor = 1.0
            s[0] = s[0][1:]
            
        elif s[0][0] == '-':
            rs.num_factor = -1.0
            s[0] = s[0][1:]
        else:
            rs.num_factor = 1.0
        

        
    
    if s[0] == '+':
        rs.num_factor = 1.0
    elif s[0] == '-':
        rs.num_factor = -1.0
    else:

        if 'sum' in s:
            if s[0][0] == '+' and len(s[0])==1:
                rs.num_factor = 1.0
                s[0] = s[0][1:]
            elif s[0][0] == '+' and len(s[0])!=1:
                t = int(s[0][1:])
                rs.num_factor = t
                
            elif s[0][0] == '-' and len(s[0])==1:
                rs.num_factor = -1.0
                s[0] = s[0][1:]
            elif s[0][0] == '-' and len(s[0])!=1:
                t = int(s[0][1:])
                rs.num_factor = -t                        
            else:        
                # print(s)
                # print('no num factor')
                strike1 = True

    if s[0] == 'frac':
        nom = float(s[1])
        denom = float(s[2])
        rs.num_factor = rs.num_factor * nom/denom
    else:
        if strike1 == True:
            print(s)
            print('no num factor')
            sys.exit(0)


    if s[1] == 'frac':
        nom = float(s[2])
        denom = float(s[3])
        rs.num_factor = rs.num_factor * nom/denom

    for k in range(0, len(s)):
        if s[k] == 'sum':
            rs.summation = list(s[k+1])
        elif s[k] == 'delta':
            rs.delta.append(list(s[k+1]))
            if spin:
                rs.delta_spin.append(list(s[k+2]))
        elif s[k] == DENS1:
            rs.coefficient.append(DENS1)
            rs.coefficient_idx.append(list(s[k+1]))
            if spin:
                rs.coefficient_spin.append(list(s[k+2]))
        elif s[k] == DENS2:
            rs.coefficient.append(DENS2)
            rs.coefficient_idx.append(list(s[k+1]))
            if spin:
                rs.coefficient_spin.append(list(s[k+2]))
        elif s[k] == DENS4:
            rs.coefficient.append(DENS4)
            rs.coefficient_idx.append(list(s[k+1]))
            if spin:
                rs.coefficient_spin.append(list(s[k+2]))
        elif s[k] == DENS3:
            rs.coefficient.append(DENS3)
            rs.coefficient_idx.append(list(s[k+1]))
            if spin:
                rs.coefficient_spin.append(list(s[k+2]))                
        elif s[k] == BARENUCL_HAM:
            rs.coefficient.append(BARENUCL_HAM)
            rs.coefficient_idx.append(list(s[k+1]))
            if spin:
                rs.coefficient_spin.append(list(s[k+2]))
        elif s[k] == TWOEL_INT_DIRAC:
            rs.coefficient.append(TWOEL_INT_DIRAC)
            idx = list(s[k+1])
            # ['i', '+', 'j', '-', '|', 't', '+', 'u', '-']
            #   0    1    2    3    4    5    6    7    8
            if spin:
                idx_lst = []
                idx_lst.append(idx[0])
                idx_lst.append(idx[2])
                idx_lst.append(idx[5])
                idx_lst.append(idx[7])

                spin_lst = []
                spin_lst.append(idx[1])
                spin_lst.append(idx[3])
                spin_lst.append(idx[6])
                spin_lst.append(idx[8])

                rs.coefficient_idx.append(idx_lst)
                rs.coefficient_spin.append(spin_lst)
            else:
                # ['i', 'j', '|', 't', 'u']
                #   0    1    2    3    4 

                idx_lst = []
                idx_lst.append(idx[0])
                idx_lst.append(idx[1])
                idx_lst.append(idx[3])
                idx_lst.append(idx[4])
                rs.coefficient_idx.append(idx_lst)



#    print('rsssss', rs)
    return rs

        
