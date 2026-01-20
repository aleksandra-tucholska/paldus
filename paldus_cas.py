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
        self.operator_stat = []

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
        if len(right.operator_stat) > 0:
            left.operator_stat = [deepcopy(right.operator_stat[0])]
            left.operator_stat.append(deepcopy(right.operator_stat[1]))
            right.operator_stat.remove(right.operator_stat[0])
            right.operator_stat.remove(right.operator_stat[0])

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

        if len(left.operator_stat) > 0:
            right.operator_stat = [deepcopy(left.operator_stat[s-1])]
            right.operator_stat.append(deepcopy(left.operator_stat[s]))
            left.operator_stat.remove(left.operator_stat[-1])
            left.operator_stat.remove(left.operator_stat[-1])

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
        res.operator_stat = deepcopy(e.operator_stat) + deepcopy(self.operator_stat)

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
        res.coefficient     = deepcopy(e.coefficient) + deepcopy(e.coefficient)
        res.coefficient_idx = deepcopy(self.coefficient_idx) + deepcopy(e.coefficient_idx)
        res.operator_idx    = deepcopy(self.operator_idx) + deepcopy(e.operator_idx)
        res.coefficient_spin    = deepcopy(self.coefficient_spin) + deepcopy(e.coefficient_spin)
        res.delta_spin    = deepcopy(self.delta_spin) + deepcopy(e.delta_spin)
        res.operator_type  = deepcopy(self.operator_type) + deepcopy(e.operator_type)
        res.num_factor      = self.num_factor * e.num_factor
        res.operator_stat = deepcopy(self.operator_stat) + deepcopy(e.operator_stat)
        
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
        new_operator_stat = []
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
                    elif self.operator_stat and self.operator_stat[j] == boso:
                        num = num - 1 # Bosons don't contribute to sign change
#                print('numnum', num)

                # Check if current contraction pair involves bosons.
                # Usually contraction is only non-zero for same statistics.
                # If either is boson, sign change logic might be different if we crossed them?
                # Actually, crossing a boson never changes sign. Crossing a fermion changes sign if crossing another fermion.
                # The loop above counts how many active operators are crossed.
                # If we cross a boson, we decremented num so it doesn't flip sign.
                # BUT, we also need to check if the operators *being contracted* affect sign?
                # No, contraction itself is a scalar replacement. The sign comes from moving them adjacent.
                # If both contracted ops are bosons, moving them doesn't generate sign.
                # If both are fermions, moving them generates sign based on how many fermions they cross.
                # My logic above: `num` counts crossed operators. If I subtract 1 for every boson crossed,
                # then `num` counts only fermions crossed.
                # Note: if operator_stat is empty, we assume fermions, so condition is false, num stays same.

                is_fermion_pair = True
                if self.operator_stat:
                     if self.operator_stat[list_of_c[i][0]] == boso:
                         is_fermion_pair = False

                if is_fermion_pair and (num%2!= 0):
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
                if len(self.operator_stat) > 0:
                    new_operator_stat.append(self.operator_stat[i])

        # print(new_operator_idx)
        self.operator_idx = new_operator_idx
        self.operator_type = new_operator_type
        if len(self.operator_stat) > 0:
            self.operator_stat = new_operator_stat


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

                # Check statistics if available
                if self.operator_stat:
                    stat0 = self.operator_stat[contr[i][j][0]]
                    stat1 = self.operator_stat[contr[i][j][1]]
                    if stat0 != stat1:
                        iszero = True
                        break

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

        # Determine stats for elements in lists
        creator_stats = []
        anlator_stats = []
        if self.operator_stat:
            for i in range(0, len(self.operator_idx)):
                if self.operator_type[i] == CRE:
                    creator_stats.append(self.operator_stat[i])
                elif self.operator_type[i] == ANI:
                    anlator_stats.append(self.operator_stat[i])

        c_trans = swap_count_stats(creator_list, creator_stats)
        a_trans = swap_count_stats(anlator_list, anlator_stats)
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

        num_anlator_pass = count_anlator_pass(self.operator_idx, self.operator_type, self.operator_stat)
        # print(num_anlator_pass)
        if (num_anlator_pass%2!=0):
            num *= -1.0
        # print('num', num)

        creator_list.sort()
        anlator_list.sort()
        ca_list = creator_list + anlator_list
        original_op_type_list = deepcopy(self.operator_type)
        original_op_stat_list = deepcopy(self.operator_stat) if self.operator_stat else []
        
        # print('-----------')
        # print(original_list)
        # print(ca_list)
        self.operator_type = []
        if original_op_stat_list:
            self.operator_stat = []

        for i in range(0, len(ca_list)):
            # print(i, ca_list[i], original_list.index(ca_list[i]))
            k = original_list.index(ca_list[i])
            self.operator_type.append(original_op_type_list[k])
            if original_op_stat_list:
                self.operator_stat.append(original_op_stat_list[k])

        # print('-----------')

        
        self.operator_idx = creator_list + anlator_list
        # print(creator_list, anlator_list)
        self.num_factor *= num

def count_anlator_pass(crean_list, crean_type, crean_stat=[]):

    l = len(crean_list)-1

    pass_count = 0
    # print('oto')
    # print(crean_list)
    # print(crean_type)
    for i in range(l, -1, -1):
        if crean_type[i] == CRE:            
            for j in range(i-1, -1, -1):
                if crean_type[j] == ANI:
                    # Only increment pass_count if both are fermions (or if stats not provided, assume fermion)
                    if not crean_stat or (crean_stat[i] == ferm and crean_stat[j] == ferm):
                        pass_count += 1

    return pass_count

def swap_count_stats(input_arr, input_stats=[]):
    if not input_stats:
        return swap_count(input_arr)
        
    # Standard swap count (Bubble sort based) to count sign changes
    # Only swaps of two fermions count as sign change.
    
    # We need to sort input_arr based on value, but we also need to carry stats with it
    # to know what we are swapping.
    
    arr = list(zip(input_arr, input_stats))
    n = len(arr)
    swaps = 0
    
    # Simple bubble sort to count inversions with stats awareness
    # Note: swap_count used sorting to determine target positions.
    # Bubble sort is appropriate for counting inversions if we want to reach sorted state.
    
    # However, standard swap_count logic was:
    # pos = sorted(enumerate(arr), key=x:x[1])
    # ...
    # This logic calculates minimum swaps to reach sorted state.
    # But for mixed statistics, we can't just count total swaps.
    # We need to simulate the swaps.
    
    # Let's perform a bubble sort and count specific swaps.
    
    temp_arr = deepcopy(arr)
    
    for i in range(n):
        for j in range(0, n-i-1):
            val1, stat1 = temp_arr[j]
            val2, stat2 = temp_arr[j+1]
            
            if val1 > val2:
                # Swap needed
                temp_arr[j], temp_arr[j+1] = temp_arr[j+1], temp_arr[j]

                # Check sign change
                if stat1 == ferm and stat2 == ferm:
                    swaps += 1
                    
    return swaps

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
