# ------------------------------------------------------
# Authors: Marcin Modrzejewski , University of Warsaw
#          Aleksandra Tucholska, University of Warsaw
# ------------------------------------------------------
from params import *
from copy import deepcopy
from itertools import permutations
from itertools import product
import collections
import math
import sys
from wick import integrate_wick_basic
from fractions import Fraction
#import numba
import re
import time
import itertools
#from numba import jit


class ugg:
    
    def __init__(self):
        self.summation = []
        self.coefficient = []
        self.coefficient_idx =[]
        self.coefficient_spin =[]
        self.operator_idx = []
        self.operator_type = []
        self.num_factor = 1.0
        self.delta = []
        self.long_delta = []
        self.hash_tuple = ()
        self.fixed_virtual = set()
        self.fixed_occupied = set()
        self.fixed_complete = set()
        self.fixed_completev = set()
        self.fixed_cabs = set()
        self.permutation_op = []
        self.number = ""
        self.name = ""
        self.constraints = {}
        self.loops = []
        self.ibra = []
        self.iket = []
        self.binary_hash = ''
        self.descriptor = 0
        self.df = []
        self.orbital_type = {}
        self.spin_dict = {}


    def remove_duplicate_summation(self):
        s1 = set(self.summtaion)
        # print(s1, list(s1))
        self.summation = list(s1)

        
    def set_permutation(self, p):
        self.permutation_op = deepcopy(p)

    def remove_coefficient(self, name):
        other = deepcopy(self)
        self.coefficient = []
        for x in other.coefficient:
            if(x != name):
                self.coefficient.append(x)

    def remove_coefficient_idx(self, i):
        other = deepcopy(self)
        self.coefficient_idx = []
        for k in range(0, len(other.coefficient_idx)):
            if k!= i:
                 self.coefficient_idx.append(other.coefficient_idx[k])
        
    def transpose(self):
        """
        Operator transposition:
        E_{pq}E_{rs} -> E_{sr}E_{qp}
        """
        
        if len(self.operator_idx) > 0:
            a = deepcopy(self.operator_idx)
            b = []
            for p, q in self.operator_idx:
                b.insert(0, [q, p])
            c = []
            for p in self.operator_type:
                c.insert(0, p)
                
            self.operator_idx = b
            self.operator_type = c

    def scale(self, s):
        self.num_factor *= s
            
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

    def is_repeated(self):
        """Return true if occupied index is
        repeated in two electron integral"""

        for x in range(0, len(self.coefficient)):
            if self.coefficient[x] == TWOEL_INT_TRANS:
                for k in range(0, len(self.coefficient_idx[x])):
                    for l in range(k+1, len(self.coefficient_idx[x])):
                        if self.coefficient_idx[x][k] == self.coefficient_idx[x][l]:
                            if self.coefficient_idx[x][k] in occupied:
                                return True
        return False

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

                            
    def ovsplit(self):
        """Split sum over general indices into sums
        over occupied and virtual indices, separately:
        \sum_p E_{ip} -> \sum_j E_{ij} + \sum_a E_{ia}
        """
#        print(self)
        genidx = False
        for k in range(len(self.summation)):
            p = self.summation[k]
            if p in general or p in complete:
                genidx = True
                break
        if not genidx:
            return arithmetic_string(self)
        else:
            eo = deepcopy(self)
            ev = deepcopy(self)
            excluded = self.indices()
            i = free_idx(occupied, excluded)
            if p in general:
                a = free_idx(virtual, excluded)
            elif p in complete:
                a = free_idx(completev, excluded)
            eo.substitute(p, i)
            ev.substitute(p, a)
            return eo.ovsplit() + ev.ovsplit()


    def cabssplit(self):

        # find pair αβ and transform to aB + Ab + AB
        # find pair αb and transform to Ab
        # print('rozpatruje self', self)
        k = 0
        for x in self.coefficient_idx:
            # print('taki wspolczynnik', x)
            k += 1
            cvlst = []
            vlst = []
            olst = []                        
            cabslst = []
            cvno = 0
            vno = 0
            ono = 0
            cabsno = 0
            for idx in x:
                # print('idx', idx)
                if idx in completev:
                    # print('completev')
                    if idx not in cvlst:
                        cvlst.append(idx)
                        cvno += 1
                elif idx in virtual:
                    # print('virt')
                    if idx not in vlst:
                        vlst.append(idx)
                        vno += 1
                elif idx in CABS:
                    if idx not in cabslst:
                        cabslst.append(idx)
                        cabsno += 1
                elif idx in occupied:
                    if idx not in olst:
                        olst.append(idx)
                        ono += 1

            # print('------------->', cvno, vno, cabsno)
            if vno  == 0 and cvno == 2:
                e1 = deepcopy(self)
                e2 = deepcopy(self)
                e3 = deepcopy(self)
                e4 = deepcopy(self)
                excluded = self.indices()
                a1 = free_idx(CABS, excluded)
                excluded = excluded | set(a1)
                a2 = free_idx(CABS, excluded)
                a3 = free_idx(virtual, excluded)
                excluded = excluded | set(a3)
                a4 = free_idx(virtual, excluded)
                alpha = cvlst[0]
                beta = cvlst[1]
                
                e1.substitute(alpha, a1)
                e1.substitute(beta, a2)
                
                e2.substitute(alpha, a1)
                e2.substitute(beta, a3)
                
                e3.substitute(alpha, a3)
                e3.substitute(beta, a1)

                e4.substitute(alpha, a3)
                e4.substitute(beta, a4)
                # print('')
                # print('utworzony 3 wyrazy')
                # print(e1)
                # print(e2)
                # print(e3)
                # print('')
                return e1.cabssplit() + e2.cabssplit() + e3.cabssplit() + e4.cabssplit()

            elif cvno == 1:
                # print('cvno1')
                e1 = deepcopy(self)
                e2 = deepcopy(self)
                excluded = self.indices()
                a1 = free_idx(CABS, excluded)
                a2 = free_idx(virtual, excluded)
                alpha =cvlst[0]
                
                e1.substitute(alpha, a1)
                e2.substitute(alpha, a2)
                # print('')
                # print('uworzone jeden wyraz male a')
                # print(e1)
                # print('')
                return e1.cabssplit() + e2.cabssplit()
            # elif ono ==1 and cvno == 1:
            #     e1 = deepcopy(self)
            #     excluded = self.indices()
            #     a1 = free_idx(CABS, excluded)
            #     alpha =cvlst[0]
                
            #     e1.substitute(alpha, a1)
            #     # print('')
            #     # print('uworzone jeden wyraz male a')
            #     # print(e1)
            #     # print('')
            #     return e1.cabssplit()

            # elif cvno ==1 and cabsno == 1:
            #     e1 = deepcopy(self)
            #     excluded = self.indices()
            #     a1 = free_idx(CABS, excluded)
            #     alpha =cvlst[0]

            #     e1.substitute(alpha, a1)
            #     # print('')
            #     # print('uworzone jeden wyraz duze A')
            #     # print(e1)
            #     # print('')
            #     return e1.cabssplit()



        return arithmetic_string(self)
        # else:
        #     if k == len(self.coefficient_idx):
        #         return self
        #     else:
        #         return self.cabssplit()
            # print("ERROR - paldus classes, cabssplit. Ugg coefficient have more than 2 virtual indices")



        # complidx = False
        # for k in range(len(self.summation)):
        #     p = self.summation[k]
        #     if p in completev:
        #         complidx = True
        #         break
        # if not complidx:
        #     return arithmetic_string(self)
        # else:
        #     ea = deepcopy(self)
        #     eAa = deepcopy(self)
        #     excluded = self.indices()
        #     i = free_idx(virtual, excluded)
        #     a = free_idx(CABS, excluded)
        #     ea.substitute(p, i)
        #     eAa.substitute(p, a)
        #     return ea.cabssplit() + eAa.cabssplit()



    def excitation_rank(self):
        """
        Return excitation rank.
        """
        # print('exc self', self)
        rank = 0
        for x in self.operator_idx:
            if x[0] == x[1]:
                rank += 0
            else:
                if x[0] in virtualall and x[1] in virtualall:
                    rank += 0
                if x[0] in occupied and x[1] in occupied:
                    rank += 0
                elif x[0] in virtualall and x[1] in occupied:
                    rank += 1
                elif x[1] in virtualall and x[0] in occupied:
                    rank += -1
                elif x[0] in general or x[1] in general:
                    print("ERROR. EXCITATION RANK CANNOT BE COMPUTED")
                    print("WHEN GENERAL INDICES ARE PRESENT")
                    sys.exit(1)
                    rank = None
#        print('rank', rank)
        return rank        


    def new_delta(self, p, q):
        """
        Append Kronecker delta, \delta_{pq} to the ugg
        instance. If appending \delta_{pq} is pointless
        regarding items that are already in self.delta list,
        do nothing.
        """

        if p == q:
            return

        if len(self.delta) == 0:
            self.delta.append([p, q])
        else:
            pointless = False
            for d in self.delta:
                if d == [p, q] or d == [q, p]:
                    pointless = True
                    break

            if not pointless:
                self.delta.append([p, q])


    def copy_without(self, lst, summation_idx=[]):
        # makes a copy of self
        # without coefficients in the lst

        e = ugg()

        e.num_factor = self.num_factor
        e.summation = []
        e.delta = self.delta
        for i in range(0, len(self.coefficient)):
            if i not in lst:
#                print('czy tu w ogole jest', i, lst)
                e.coefficient.append(self.coefficient[i])
                e.coefficient_idx.append(self.coefficient_idx[i])
                #               print(e)
        # print('teraz gole e', e)
        # print('tymczasem self.sum', self.summation)
        # for s in self.summation:
        #     for i in e.coefficient_idx:
        #         if s in i:
        #             if s not in e.summation:
        #                 print('dodaje do sum e')
        #                 e.summation.append(s)
        # print('nowy e', e)
        return e

    # def int_type(self, i):
    #     if self.coefficient_idx[i][0] in virt_all:
    #         if self.coefficient_idx[i][1] in virt_all:
    #             if self.coefficient_idx[i][2] in virt_all:
    #                 if self.coefficient_idx[i][3] in virt_all:
    #                     g_type = "vvvv"
    #                 else:
    #                     g_type = "vvvo"
    #             else:
    #                 g_type = "vvoo"
    #         else:
    #             if self.coefficient_idx[i][2] in virt_all:
    #                 g_type = "vovo"
    #             else:
    #                 g_type = "vooo"
    #     else:
    #         g_type = "oooo"

    #     return g_type


    # def int_type_f12(self, i):


        
    #     if self.coefficient_idx[i][0] in CABS_all:
    #         if self.coefficient_idx[i][2] in virt_all:
    #             g_type = "covo"
    #         elif self.coefficient_idx[i][2] in CABS_all:
    #             g_type = "coco"
    #         else:
    #             print('unknown f12 integral type. Please check paldus_classes int_type_f12')
    #     if self.coefficient_idx[i][3] in CABS_all:
    #         if self.coefficient_idx[i][0] in virt_all:
    #             g_type = "vooc"
    #         else:
    #             print('unknown f12 integral type. Please check paldus_classes int_type_f12')

    #     if self.coefficient_idx[i][0] in virt_all:
    #         if self.coefficient_idx[i][1] in CABS_all:
    #             g_type = "vcoo"
    #         else:
    #             print('unknown f12 integral type. Please check paldus_classes int_type_f12')

    #     if self.coefficient_idx[i][1] in virt_all:
    #         if self.coefficient_idx[i][3] in CABS_all:
    #             g_type = "ovoc"
    #         else:
    #             print('unknown f12 integral type. Please check paldus_classes int_type_f12')
    #     if self.coefficient_idx[i][1] in CABS_all:
    #         if self.coefficient_idx[i][3] in virt_all:
    #             g_type = "ocov"
    #         else:
    #             print('unknown f12 integral type. Please check paldus_classes int_type_f12')
    #     if self.coefficient_idx[i][1] in CABS_all:
    #         if self.coefficient_idx[i][3] in CABS_all:
    #             g_type = "ococ"
    #         else:
    #             print('unknown f12 integral type. Please check paldus_classes int_type_f12')

                
    #     print(self)

    #     return g_type

    
    def nvirtocc(self,i):

        nvirt = 0
        nocc = 0
        ncabs = 0
        
        if self.coefficient_idx[i][0] in occ_all:
            nocc += 1
        else:
            if self.coefficient_idx[i][0] in CABS_all:
                ncabs += 1
            else:
                nvirt += 1
        if self.coefficient_idx[i][1] in occ_all:
            nocc += 1
        else:
            if self.coefficient_idx[i][1] in CABS_all:
                ncabs += 1
            else:
                nvirt += 1
        if self.coefficient_idx[i][2] in occ_all:
            nocc += 1
        else:
            if self.coefficient_idx[i][2] in CABS_all:
                ncabs += 1
            else:
                nvirt += 1
        if self.coefficient_idx[i][3] in occ_all:
            nocc += 1
        else:
            if self.coefficient_idx[i][3] in CABS_all:
                ncabs += 1
            else:
                nvirt += 1
            
        return nvirt, nocc, ncabs

    def int_type(self, prefix, i):
        nvirt, nocc, ncabs = self.nvirtocc(i)
        v_type = ""
        
        if nvirt == 4:
            v_type = "read_ftvvvv"
            return v_type
        if nocc == 4:
            v_type = prefix+"oooo"
            return v_type


        if nocc == 3:
            if nvirt == 1:
                if self.coefficient_idx[i][0] in virt_all:                    
                    v_type = "vooo"
                    v_type = prefix+v_type
                    return v_type

                elif self.coefficient_idx[i][1] in virt_all:                    
                    v_type = "ovoo"
                    v_type = prefix+v_type
                    return v_type

                else:
                    pairs = [self.coefficient_idx[i][0:2],\
                             self.coefficient_idx[i][2:4]]
                    self.coefficient_idx[i] = pairs[1] + pairs[0]
                    v_type = self.int_type(prefix, i)
 
            if ncabs == 1:
                if self.coefficient_idx[i][0] in CABS_all:                    
                    v_type = "cooo"
                    v_type = prefix+v_type
                    return v_type

                elif self.coefficient_idx[i][1] in CABS_all:                    
                    v_type = "ocoo"
                    v_type = prefix+v_type
                    return v_type

                else:
                    pairs = [self.coefficient_idx[i][0:2],\
                             self.coefficient_idx[i][2:4]]
                    self.coefficient_idx[i] = pairs[1] + pairs[0]
                    v_type = self.int_type(prefix, i)
 
        if nocc == 2:
            if ncabs == 0:
                if self.coefficient_idx[i][0] in virt_all and self.coefficient_idx[i][1] in virt_all:
                    v_type = "vvoo"
                    v_type = prefix+v_type
                    return v_type

                elif self.coefficient_idx[i][0] in virt_all and self.coefficient_idx[i][2] in virt_all:
                    v_type = "vovo"
                    v_type = prefix+v_type
                    return v_type

                elif self.coefficient_idx[i][0] in virt_all and self.coefficient_idx[i][3] in virt_all:
                    v_type = "voov"
                    v_type = prefix+v_type
                    return v_type

                elif self.coefficient_idx[i][1] in virt_all and self.coefficient_idx[i][3] in virt_all:
                    v_type = "ovov"
                    v_type = prefix+v_type
                    return v_type

                else:
                    pairs = [self.coefficient_idx[i][0:2],\
                             self.coefficient_idx[i][2:4]]
                    self.coefficient_idx[i] = pairs[1] + pairs[0]
                    v_type = self.int_type(prefix, i)

            if nvirt == 0:
                if self.coefficient_idx[i][0] in CABS_all and self.coefficient_idx[i][1] in CABS_all:
                    v_type = "ccoo"
                    v_type = prefix+v_type
                    return v_type

                elif self.coefficient_idx[i][0] in CABS_all and self.coefficient_idx[i][2] in CABS_all:
                    v_type = "coco"
                    v_type = prefix+v_type
                    return v_type

                elif self.coefficient_idx[i][0] in CABS_all and self.coefficient_idx[i][3] in CABS_all:
                    v_type = "cooc"
                    v_type = prefix+v_type
                    return v_type

                elif self.coefficient_idx[i][1] in CABS_all and self.coefficient_idx[i][3] in CABS_all:
                    v_type = "ococ"
                    v_type = prefix+v_type
                    return v_type

                else:
                    pairs = [self.coefficient_idx[i][0:2],\
                             self.coefficient_idx[i][2:4]]
                    self.coefficient_idx[i] = pairs[1] + pairs[0]
                    v_type = self.int_type(prefix, i)

            if nvirt == 1 and ncabs == 1:
                # print('tak tutaj dwa dwa', self)
                if self.coefficient_idx[i][0] in virt_all and self.coefficient_idx[i][1] in CABS_all:
                    # print('h1')
                    v_type = "vcoo"
                    v_type = prefix+v_type
                    return v_type

                elif self.coefficient_idx[i][0] in CABS_all and self.coefficient_idx[i][1] in virt_all:
                                        
                    v_type = "cvoo"
                    v_type = prefix+v_type
                    return v_type

                elif self.coefficient_idx[i][0] in CABS_all and self.coefficient_idx[i][2] in virt_all:

                    v_type = "covo"
                    v_type = prefix+v_type
                    return v_type

                elif self.coefficient_idx[i][0] in CABS_all and self.coefficient_idx[i][3] in virt_all:
                    v_type = "coov"
                    v_type = prefix+v_type
                    return v_type

                elif self.coefficient_idx[i][1] in CABS_all and self.coefficient_idx[i][3] in virt_all:

                    v_type = "ocov"
                    v_type = prefix+v_type
                    return v_type

                elif self.coefficient_idx[i][1] in CABS_all and self.coefficient_idx[i][2] in virt_all:

                    v_type = "ocvo"
                    v_type = prefix+v_type
                    return v_type

                else:
                                        
                    pairs = [self.coefficient_idx[i][0:2],\
                             self.coefficient_idx[i][2:4]]
                    self.coefficient_idx[i] = pairs[1] + pairs[0]
                    v_type = self.int_type(prefix, i)


        if nocc == 1:
            if nvirt==3:
                if self.coefficient_idx[i][3] in occ_all:
                    v_type = "vvvo"
                    v_type = prefix+v_type
                    return v_type

                elif self.coefficient_idx[i][2] in occ_all:
                    v_type = "vvov"
                    v_type = prefix+v_type
                    return v_type

                else:
                    pairs = [self.coefficient_idx[i][0:2],\
                             self.coefficient_idx[i][2:4]]
                    self.coefficient_idx[i] = pairs[1] + pairs[0]
                    v_type = self.int_type(prefix, i)

            if ncabs == 3:
                if self.coefficient_idx[i][3] in occ_all:
                    v_type = "ccco"
                    v_type = prefix+v_type
                    return v_type

                elif self.coefficient_idx[i][2] in occ_all:
                    v_type = "ccoc"
                    v_type = prefix+v_type
                    return v_type

                else:
                    pairs = [self.coefficient_idx[i][0:2],\
                             self.coefficient_idx[i][2:4]]
                    self.coefficient_idx[i] = pairs[1] + pairs[0]
                    v_type = self.int_type(prefix, i)

            if nvirt == 2 and ncabs == 1:
                if self.coefficient_idx[i][1] in occ_all:
                    if self.coefficient_idx[i][0] in CABS_all:
                        v_type = "covv"
                        v_type = prefix+v_type
                        return v_type

                elif self.coefficient_idx[i][0] in occ_all:
                    if self.coefficient_idx[i][1] in CABS_all:
                        v_type = "ocvv"
                        v_type = prefix+v_type
                        return v_type

                else:
                    pairs = [self.coefficient_idx[i][0:2],\
                             self.coefficient_idx[i][2:4]]
                    self.coefficient_idx[i] = pairs[1] + pairs[0]
                    v_type = self.int_type(prefix, i)



            if nvirt == 1 and ncabs == 2:
                if self.coefficient_idx[i][2] in occ_all:
                    if self.coefficient_idx[i][3] in virt_all:
                        v_type = "ccov"
                        v_type = prefix+v_type
                        return v_type

                    else:
                        pairs = [self.coefficient_idx[i][0:2],\
                                 self.coefficient_idx[i][2:4]]
                        self.coefficient_idx[i] = pairs[1] + pairs[0]
                        v_type = self.int_type(prefix, i)
                else:
                    pairs = [self.coefficient_idx[i][0:2],\
                             self.coefficient_idx[i][2:4]]
                    self.coefficient_idx[i] = pairs[1] + pairs[0]
                    v_type = self.int_type(prefix, i)


                if self.coefficient_idx[i][3] in occ_all:
                    if self.coefficient_idx[i][2] in virt_all:
                        v_type = "ccvo"
                        v_type = prefix+v_type
                        return v_type

                    else:
                        pairs = [self.coefficient_idx[i][0:2],\
                                 self.coefficient_idx[i][2:4]]
                        self.coefficient_idx[i] = pairs[1] + pairs[0]
                        v_type = self.int_type(prefix, i)
                else:
                    pairs = [self.coefficient_idx[i][0:2],\
                             self.coefficient_idx[i][2:4]]
                    self.coefficient_idx[i] = pairs[1] + pairs[0]
                    v_type = self.int_type(prefix, i)

                
        if v_type == "":
            sys.exit(0)
        return v_type

    def int1_type_trans(self, i):
        if self.coefficient_idx[i][0] in virt_all:
            if self.coefficient_idx[i][1] in virt_all:
                k_type = "tvv"
            else:
                k_type = "tvo"
        else:
            if self.coefficient_idx[i][1] in virt_all:
                k_type = "tov"
            else:
                k_type = "too"

        return k_type
    
    def establish_fixed(self, fixed_fx = []):
        """ establish fixed inidces
        """

        # fixed  = ['a', 'b', 'i', 'j', 'A', 'B']
            
        for i in range(0, len(self.coefficient)):
            for j in self.coefficient_idx[i]:
                if j not in self.summation and j not in fixed:
                    fixed.append(j)

        for i in self.operator_idx:
            for j in i:
                if j not in self.summation and j not in fixed:
                    fixed.append(j)
                    
        for i in self.delta:
            for j in i:
                if j not in self.summation and j not in fixed:
                    fixed.append(j)

        if fixed_fx != []:
            self.clear_fixed()
            for x in fixed_fx:
                fixed.append(x)



    def establish_fixed_and_return(self):
        """ establish fixed inidces
        """


        fixed = []
        
        for i in range(0, len(self.coefficient)):
            for j in self.coefficient_idx[i]:
                if j not in self.summation and j not in fixed:
                    fixed.append(j)

        for i in self.operator_idx:
            for j in i:
                if j not in self.summation and j not in fixed:
                    fixed.append(j)
                    
        for i in self.delta:
            for j in i:
                if j not in self.summation and j not in fixed:
                    fixed.append(j)

        return fixed

            
    def clear_fixed(self):
        """ clear fixed inidces
        """
        fixed.clear()
        # fx = deepcopy(fixed)
        # for x in fx:
        #     fixed.remove(x)
    
    def clear_deltas(self):

        if len(self.delta) > 0:
            for i in range(0, len(self.delta)):
                if self.delta[i][0] != self.delta[i][1]:
                    self.num_factor = 0.0
                    
    def amplitude_type(self, i):
        ln = str(len(self.coefficient_idx[i])//2)
        if self.coefficient[i] == CC_AMPLITUDE:
            at = CC_AMPLITUDE + ln
        if self.coefficient[i] == EOM_CC_AMPLITUDE_R:
            at = EOM_CC_AMPLITUDE_R + ln
        if self.coefficient[i] == EOM_CC_SINGLE_Rl:
            at = EOM_CC_AMPLITUDE_R + ln
        if self.coefficient[i] == EOM_CC_SINGLE_Rr:
            at = EOM_CC_AMPLITUDE_R + ln
        if self.coefficient[i] == EOM_CC_SINGLE_Rl_plus:
            at = EOM_CC_AMPLITUDE_R + ln
        if self.coefficient[i] == EOM_CC_SINGLE_Rr_plus:
            at = EOM_CC_AMPLITUDE_R + ln
        if self.coefficient[i] == EOM_CC_SINGLE_Rl_minus:
            at = EOM_CC_AMPLITUDE_R + ln
        if self.coefficient[i] == EOM_CC_SINGLE_Rr_minus:
            at = EOM_CC_AMPLITUDE_R + ln
        if self.coefficient[i] == EOM_CC_SINGLE_Ll:
            at = EOM_CC_AMPLITUDE_L + ln
        if self.coefficient[i] == EOM_CC_SINGLE_Lr:
            at = EOM_CC_AMPLITUDE_L + ln
        if self.coefficient[i] == EOM_CC_AMPLITUDE_L:
            at = EOM_CC_AMPLITUDE_L + ln
        elif self.coefficient[i] == S_AMPLITUDE:
            at = S_AMPLITUDE + ln
        return ln, at


    def integrate(self, bra = [], ket = [], braspin = [], ketspin = [], nodelta=False):
        """Compute <0|self|0> integral. """
        return integrate(self, bra, ket, braspin, ketspin, nodelta)

    def cabstransform(self):
        return cabstransform(self)

    def integrate_wick(self, bra = [], ket = [], braspin = [], ketspin = []):
        """Compute <0|self|0> integral. """
        return integrate_wick(self, bra, ket, braspin, ketspin)

    def binary_hash_gen(self):
        """Generates binary code for string
        generated from standarized ugg"""

        if self is None:
            return str(0)
        if self.num_factor == 0.0 or self.num_factor == 0:
            return str(0)

        summ = ""
        opidx = ""
        cidx = ""
        temp = ""
        delidx = ""
        sumidx = ""

        if self.summation != []:
            for x in self.summation:
                summ = summ + str(x) + str("")
            sumidx = "sum" + summ + "|"
            upperi = ""
            loweri = ""
            if len(self.operator_idx) != 0:
                    for x in range(0, len(self.operator_idx)):
                        temp = ""
                        for y in self.operator_idx[x]:
                            temp = temp + str(y) + str("")
                    if self.operator_type[x] == "t0":
                        opidx = opidx + "T"+temp+"|"
                    elif self.operator_type[x] == "t1":
                        opidx = opidx + "T1"+temp+"|"
                    elif self.operator_type[x] == "tm1":
                        opidx = opidx + "Tm1"+temp+"|"
                    elif self.operator_type[x] == "s":
                        opidx = opidx + "E"+temp+"|"


        if len(self.coefficient) != 0:
            for x in range(0, len(self.coefficient)):
                temp = ""

                if self.coefficient[x] == CC_AMPLITUDE or self.coefficient[x] == S_AMPLITUDE\
                        or self.coefficient[x] == TEMP1 or self.coefficient[x] == TEMP2 \
                        or self.coefficient[x] == EOM_TRIPLET_R3 \
                        or self.coefficient[x] == EOM_TRIPLET_R2p \
                        or self.coefficient[x] == EOM_TRIPLET_R2m \
                        or self.coefficient[x] == EOM_TRIPLET_R1 \
                        or self.coefficient[x] == EOM_CC_SINGLE_Rr \
                        or self.coefficient[x] == EOM_CC_SINGLE_Rl \
                        or self.coefficient[x] == EOM_CC_SINGLE_Rr_plus \
                        or self.coefficient[x] == EOM_CC_SINGLE_Rl_plus \
                        or self.coefficient[x] == EOM_CC_SINGLE_Rr_minus \
                        or self.coefficient[x] == EOM_CC_SINGLE_Rl_minus \
                        or self.coefficient[x] == F12_AMPLITUDE \
                        or self.coefficient[x] == F12_TWOEL \
                        or self.coefficient[x] == CI_AMPLITUDE \
                        or self.coefficient[x] == TWOEL_INT_AS \
                        or self.coefficient[x] == SLAT2_SYM1 \
                        or self.coefficient[x] == SLAT2_SYM2 \
                        or self.coefficient[x] == SLAT3 \
                        or self.coefficient[x] == SLAT4:
                    for y in self.coefficient_idx[x] :
                        temp += str(y) + str("")
                else:
                    for y in self.coefficient_idx[x]:
                        temp = temp + str(y) + str("")
                if (self.coefficient[x] == EOM_CC_SINGLE_Rr or self.coefficient[x] == EOM_CC_SINGLE_Rl):
                    coef = 'R'
                elif (self.coefficient[x] == EOM_CC_SINGLE_Rr_plus or self.coefficient[x] == EOM_CC_SINGLE_Rl_plus):
                    coef = '+R'
                elif (self.coefficient[x] == EOM_CC_SINGLE_Rr_minus or self.coefficient[x] == EOM_CC_SINGLE_Rl_minus):
                    coef = '-R'
                elif (self.coefficient[x] == SLAT2_SYM1 or self.coefficient[x] == SLAT2_SYM2 or self.coefficient[x] == SLAT3 or \
                          self.coefficient[x] == SLAT4):
                    coef = "o"
                else:
                    coef = self.coefficient[x]

                cidx = cidx + str(coef)+ ""+temp+"|"

        for x in range(0, len(self.delta)):
            temp = ""
            for y in self.delta[x]:
                temp = temp + str(y) + str("")
            delidx = delidx + "D" + temp + "|"

        orbtype = ""
        if self.orbital_type != {}:
            # print('tak robie orbital type')
            for c_list in self.coefficient_idx:
                for idx in c_list:
                    orbtype = orbtype + str(idx)+ str(self.orbital_type[idx])
                    
            

        ascii = sumidx + cidx + opidx + delidx + orbtype
        self.binary_hash = bin(int.from_bytes(ascii.encode(), 'big'))



    def hash(self):
        """Contents of hash n-tuple:
        1. Number of summation indices,
        2. Number of fixed occupied indices,
        3. Number of fixed virtual indices,
        4. Ugg order (number of Epq's),
        5. Number of coefficients,
        6. Number of Kronecker deltas,
        7. Set of permuted pairs

        A subset of ugg instances yielding the same hash n-tuple
        belong to the same cluster of the arhtmetic_string.cluster
        method.
        """

        l = []
        l.append(len(self.summation))
#        l.append(self.summation)
        
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
        
        self.fixed_occupied = set()
        self.fixed_virtual = set()

        for p in c:
            if p in occupied:
                self.fixed_occupied.add(p)
            if p in virtual:
                self.fixed_virtual.add(p)
            if p in complete:
                self.fixed_complete.add(p)
            if p in completev:
                self.fixed_completev.add(p)
            if p in cabs:
                self.fixed_cabs.add(p)

        l.append(len(self.fixed_occupied))
        l.append(len(self.fixed_virtual))
        l.append(len(self.fixed_complete))
        l.append(len(self.fixed_completev))
        l.append(len(self.fixed_cabs))
        l.append(len(self.operator_idx))
        l.append(len(self.coefficient))
        l.append(len(self.delta))
        l.append(len(self.permutation_op))
#        l.append(self.coefficient)

        self.hash_tuple = tuple(l)

    def recreate_fixed(self, fixed_list, typ):
        # print('recfix1', self)
        self.multisubst(fixed_list, ['y', 'z'])        

        if typ == 'oo':
            if 'i' in self.summation:
                for j in range(2, len(occupied)):
                    if occupied[j] not in self.summation:
                        self.multisubst(['i'], [occupied[j]])
                        break
            if 'j' in self.summation:
                for j in range(2, len(occupied)):
                    if occupied[j] not in self.summation:
                        self.multisubst(['j'], [occupied[j]])
                        break
            self.multisubst(['y', 'z'], ['i', 'j'])

        if typ == 'ov':
            if 'i' in self.summation:
                for j in range(1, len(occupied)):
                    if occupied[j] not in self.summation:
                        self.multisubst(['i'], [occupied[j]])
                        break
            if 'a' in self.summation:
                for j in range(1, len(virtual)):
                    if virtual[j] not in self.summation:
                        self.multisubst(['a'], [virtual[j]])
                        break
            self.multisubst(['y', 'z'], ['i', 'a'])

        if typ == 'vo':
            if 'a' in self.summation:
                for j in range(1, len(virtual)):
                    if virtual[j] not in self.summation:
                        self.multisubst(['a'], [virtual[j]])
                        break
            if 'i' in self.summation:
                for j in range(1, len(occupied)):
                    if occupied[j] not in self.summation:
                        self.multisubst(['i'], [occupied[j]])
                        break
            self.multisubst(['y', 'z'], ['a', 'i'])

        if typ == 'vv':
            if 'a' in self.summation:
                for j in range(1, len(virtual)):
                    if virtual[j] not in self.summation:
                        self.multisubst(['a'], [virtual[j]])
                        break
            if 'b' in self.summation:
                for j in range(2, len(virtual)):
                    if virtual[j] not in self.summation:
                        self.multisubst(['b'], [virtual[j]])
                        break
            self.multisubst(['y', 'z'], ['a', 'b'])

    def general_split(self, fixed_list, typ):
        
        print("THIS IS NOT IMPLEMENTED FOR COMPLETE AND CABS INDICES - FIX")
        print("exiting...")
        sys.exit(1)
        
        rsimp = arithmetic_string()
        
        if typ == 'po':
            for i in range(0, 2):
                rsimp.append(self)
            for i in occupied:
                if i not in self.summation:
                    rsimp[0].multisubst([fixed_list[0]], [i])
                    break
            for i in virtual:
                if i not in self.summation:
                    rsimp[1].multisubst([fixed_list[0]], [i])
                    break
        elif typ == 'op':
            for i in range(0, 2):
                rsimp.append(self)
            for i in occupied:
                if i not in self.summation:
                    rsimp[0].multisubst([fixed_list[1]], [i])
                    break
            for i in virtual:
                if i not in self.summation:
                    rsimp[1].multisubst([fixed_list[1]], [i])
                    break
        elif typ == 'pp':
            for i in range(0, 4):
                rsimp.append(self)
            for i in range(0, len(occupied)):
                if occupied[i] not in self.summation:
                    rsimp[0].multisubst([fixed_list[0]], [i])
                    for j in range(i+1, len(occupied)):
                        if occupied[j] not in self.summation:
                            rsimp[0].multisubst([fixed_list[1]], [j])
                            break
            for i in range(0, len(occupied)):
                if occupied[i] not in self.summation:
                    rsimp[1].multisubst([fixed_list[0]], [i])
                    for i in range(0, len(virtual)):
                        if virtual[i] not in self.summation:
                            rsimp[1].multisubst([fixed_list[1]], [i])
                            break
            for i in range(0, len(occupied)):
                if occupied[i] not in self.summation:
                    rsimp[2].multisubst([fixed_list[1]], [i])
                    for i in range(0, len(virtual)):
                        if virtual[i] not in self.summation:
                            rsimp[2].multisubst([fixed_list[0]], [i])
                            break
            for i in range(0, len(virtual)):
                if virtual[i] not in self.summation:
                    rsimp[3].multisubst([fixed_list[0]], [i])
                    for j in range(i+1, len(virtual)):
                        if virtual[j] not in self.summation:
                            rsimp[3].multisubst([fixed_list[1]], [j])
                            break


        return rsimp

    def multisubst(self, refer, subst):
        """
        Changes all indices in refer list to corresponding indices
        in subst. 
        """
        # print('')
        # print('multisubst')
        # print(self, refer, subst, len(refer), len(subst))
        if len(refer) != len(subst):
            print("Error. Nonequal lenghts of reference and permuted lists.")
            return

        refer = deepcopy(refer)
        changed_list = []
        j = 1
        for i in range(0, len(refer)):
            if subst[i] != refer[i]:
                if subst[i] in changed_list:
                    self.substitute(refer[i], subst[i])
                elif subst[i] not in refer:
                    self.substitute(refer[i], subst[i])
                else:
                    self.substitute(subst[i], j)
                    m = refer.index(subst[i])
                    refer[m] = j
                    j += 1
                    self.substitute(refer[i], subst[i])

                changed_list.append(refer[i])


    def fixed_free_pairs(self, fixed_indices, free_indices):
        #
        # Generate PAIRS dictionary. PAIRS[P], where P is
        # a fixed index and Counter is a class instance, which contains
        # info how many times each summation (free) index appears
        # in pair with P.
        #
#        print('')
#        print(self)
        pairs = {}
#        # print('fixed ind', fixed_indices)
        for p in fixed_indices:
            pairs[p] = collections.Counter()
        # print('FXX', pairs)
        if len(self.coefficient) > 0:
            for i in range(0, len(self.coefficient)):
                name = self.coefficient[i]
                # print(name)
                indices = self.coefficient_idx[i]
#                # print(indices)
                nidx = len(indices)
                nfixed = 0
                for p in indices:
                    if p in fixed_indices:
                        nfixed += 1
                # print(nfixed, nidx)
                #
                # Test if there can be any fixed-free pair
                #
                if nfixed == 0 or nfixed == nidx:
                    continue
#                # print('indices', indices, HAS_PAIR_INDICES)
                if name in HAS_PAIR_INDICES:
                    for k in range(0, nidx, 2):
                        # print(name, indices[k:k+2], k, nidx, fixed_indices)
                        f = filterout(indices[k:k+2], fixed_indices)
                        # print('po filterout', f)
                        if len(f) == 1:
                            # print('tak, mam jeden sumacyjny i jeden fixed')
                            s = filterout(indices[k:k+2], free_indices)
                            # print('i to jest wyfiltrowany ten sumacyjny', s)
                            if len(s) == 1:
                                # print('updejtuje indeks', f[0], 'indeksem', s[0])
                                pairs[f[0]].update(s[0])
#                    # print('laaa', pairs)
                else:
#                    # print('')
#                    # print(name)
                    #
                    # Coefficients of unknown permutational symmetry
                    #
                    f = filterout(indices, fixed_indices)
#                    # print('fwww', f)
                    if len(f) > 0:
                        s = filterout(indices, free_indices)
#                        # print('sw', s)
                        if len(s) > 0:
                            for p in f:
                                for q in s:
#                                    # print(p, q)
                                    pairs[p].update(q)
#        # print('FXX', pairs)
        return pairs

    def clean_pairs_from_equivalent(self, pairs):

        # If two summation indices in pairs, 
        # have the same counting, do not change their names
    
        # 1. Find all sets of occ, virt, comp, compv 
        # present in fixed pairs
    
        occ_pairs = []
        virt_pairs = []
        comp_pairs = []
        compv_pairs = []
        gen_pairs = []
    
        for x in pairs:
            for y in pairs[x]:
                # print(x, pairs[x], y)
                if y in occupied:
                    occ_pairs.append([y, pairs[x][y]])
                elif y in virtual:
                    virt_pairs.append([y, pairs[x][y]])
                elif y in complete:
                    comp_pairs.append([y, pairs[x][y]])
                elif y in completev:
                    compv_pairs.append([y, pairs[x][y]])
                elif y in general:
                    gen_pairs.append([y, pairs[x][y]])


        # print('dupa1', occ_pairs)
        # print('dupa2', virt_pairs)

        del_list = []
        for x in range(0, len(occ_pairs)):
            for y in range(x+1, len(occ_pairs)):
                if occ_pairs[x][1] == occ_pairs[y][1]:
                    if occ_pairs[x][0] not in del_list:
                        del_list.append(occ_pairs[x][0])
                    if occ_pairs[y][0] not in del_list:
                        del_list.append(occ_pairs[y][0])
        for x in range(0, len(virt_pairs)):
            for y in range(x+1, len(virt_pairs)):
                if virt_pairs[x][1] == virt_pairs[y][1]:
                    if virt_pairs[x][0] not in del_list:
                        del_list.append(virt_pairs[x][0])
                    if virt_pairs[y][0] not in del_list:
                        del_list.append(virt_pairs[y][0])
        for x in range(0, len(comp_pairs)):
            for y in range(x+1, len(comp_pairs)):
                if comp_pairs[x][1] == comp_pairs[y][1]:
                    if comp_pairs[x][0] not in del_list:
                        del_list.append(comp_pairs[x][0])
                    if comp_pairs[y][0] not in del_list:
                        del_list.append(comp_pairs[y][0])
        for x in range(0, len(compv_pairs)):
            for y in range(x+1, len(compv_pairs)):
                if compv_pairs[x][1] == compv_pairs[y][1]:
                    if compv_pairs[x][0] not in del_list:
                        del_list.append(compv_pairs[x][0])
                    if compv_pairs[y][0] not in del_list:
                        del_list.append(compv_pairs[y][0])
        for x in range(0, len(gen_pairs)):
            for y in range(x+1, len(gen_pairs)):
                if gen_pairs[x][1] == gen_pairs[y][1]:
                    if gen_pairs[x][0] not in del_list:
                        del_list.append(gen_pairs[x][0])
                    if gen_pairs[y][0] not in del_list:
                        del_list.append(gen_pairs[y][0])

        # print('occ_list', occ_pairs)
        # print('ddddd', del_list)
        for x in del_list:
            for y in pairs:
                if x in pairs[y]:
                    del pairs[y][x]

        # print('pairs po  zmianie', pairs)
        return pairs


    def generate_descriptor(self, fixed_indices):
        """
        Generate descriptor of coefficients and indices in
        a form of the following list:
        [ ..., [X, N, NOCC, [[A, I], [B, C], ...], ... ],
        where X denotes name of a coefficient, N is the number of all indices of
        the coefficient, and NOCC is nuber of occupied indices. A descriptor of
        a single coefficient is generated as follows:
        
        (a, b - fixed) g_{bakl}    -> ["g", 4, 2, [["a", "b"], [   ]]]
        (a, i - fixed) t^{ab}_{ji} -> ["t", 4, 2, [["a"     ], ["i"]]]

        -------------------------------------------------------------------------
        Any two UGG instances that are eqivalent generate identical coefficient
        descriptors. It is possible that nonequivalent UGGs generate identical 
        descriptors. The sole purpose of using coefficient descriptors is to
        filter out the UGG instances that are nonequvalent in an obvious way.
        -------------------------------------------------------------------------
        """

        descriptor = []
        if len(self.coefficient) > 0:
            for i in range(0, len(self.coefficient)):
                name = self.coefficient[i]
                indices = self.coefficient_idx[i]
                nidx = len(indices)
                nocc = 0
                for p in indices:
                    if p in occupied:
                        nocc += 1
                idxclasses = []

                if name in HAS_PAIR_INDICES:
                    for k in range(0, nidx, 2):
                        f = filterout(indices[k:k+2], fixed_indices)
                        f.sort()
                        idxclasses.append(f)
                else:
                    #
                    # Coefficients of unknown permutational symmetry
                    #
                    f = filterout(indices, fixed_indices)
                    f.sort()
                    idxclasses.append(f)

                idxclasses.sort()
                descriptor.append([name, nidx, nocc, idxclasses])
        #
        # Empty list is returned if there are
        # no symbolic coefficients
        #
        descriptor.sort()
        return descriptor


    def __eq__(self, other):
        """
        Note: Commutation relations are not utilized during comparison.
        Equivalent strings of operators can return 0.

        ------------------------------------------------------------------
        SELF and OTHER are assumed to be in standard form (*standarized*).
        ------------------------------------------------------------------

        Equivalence *modulo* multiplicative factor. Phase is 
        returned because of the presence of antisymmetric coefficients.
        If C is to be a sum of A and B equivalent UGG instances,
        C can be generated as follows:

        PHASE = A == B
        if PHASE != 0:
              C = deepcopy(A)
              C.num_factor += PHASE * B.num_factor
              #
              # Now C is equivalent to A + B
              #
        """
        #
        # If permutation operator is present, extract
        # all terms permuted from the base one and compare.
        # Permutation operator can by only by dependent on
        # fixed indices.
        #
        if self.permutation_op != []:
            if self.permutation_op.sort() == other.permutation_op.sort():
                a = deepcopy(self)
                a.permutation_op = []
                b = deepcopy(other)
                b.permutation_op = []
                #
                # Generate all permutations
                #
                l = a.permute(self.permutation_op)
                l1 = 0
                for x in l:
                    phase = b == x
                    if phase != 0:
                        l1 = phase
                        break

                return l1
            else:
                return 0
        #
        # Code below is executed if there are no
        # permutation operators
        #
        indices_a = self.indices()
        indices_b = other.indices()
        summation_a = set(self.summation)
        summation_b = set(other.summation)
        fixed_a = indices_a - summation_a
        fixed_b = indices_b - summation_b

        if not (indices_a == indices_b and fixed_a == fixed_b):
            return 0
        coeffs_a = deepcopy(self.coefficient)
        coeffs_b = deepcopy(other.coefficient)
        coeffs_a.sort()
        coeffs_b.sort()
        if not coeffs_a == coeffs_b:
            return 0
        #
        # Generate descriptors of coefficients.
        # Descriptors provide an approximate way of comparing UGGs, and
        # perform well in filtering out UGGs that are for sure not equal.
        #
        descriptor_self = self.generate_descriptor(fixed_a)
        descriptor_other = other.generate_descriptor(fixed_b)
        if descriptor_self != descriptor_other:
            return 0
        #
        # Compare coefficients
        #
        phase = 1
        for k in range(len(self.coefficient)):
            name_a = self.coefficient[k]
            name_b = other.coefficient[k]
            idx_a = self.coefficient_idx[k]
            idx_b = other.coefficient_idx[k]
            phase *= compare_coeffs(name_a, idx_a, name_b, idx_b)
            if phase == 0:
                return 0

        #
        # Compare E_{pq} operators and Kronecker deltas
        #
        l1 = self.operator_idx == other.operator_idx
        l3 = self.operator_type == other.operator_type
        l2 = self.delta == other.delta

        if l1 and l2 and l3:
            return phase
        else:
            return 0
        
        
    def rename_fixed(self, ar):
        """ Renames all fixed indices form "i" to "0i"
        or, from "0i" to "i" if ar = 'remove'
        """
        changed = []

        for i in range(0, len(self.coefficient_idx)):
            for j in range(0, len(self.coefficient_idx[i])):
                old = self.coefficient_idx[i][j]

                if old not in self.summation:
                    if ar == "add":
                        if old[1:] not in changed:
                            new = str("0") + self.coefficient_idx[i][j]
                            self.substitute(old, new)
                            changed.append(old)
                    elif ar == "remove":
                        if old not in changed:
                            new = self.coefficient_idx[i][j][1:]
                            self.substitute(old, new)
                            changed.append(new)
        for i in range(0, len(self.operator_idx)):
            for j in range(0, len(self.operator_idx[i])):
                old = self.operator_idx[i][j]
                if old not in self.summation:
                    if old[1:] not in changed:
                        if ar == "add":
                            new = str("0") + self.operator_idx[i][j]
                        elif ar == "remove":
                            new = self.operator_idx[i][j][1:] 
                        self.substitute(old, new)
                        changed.append(old)
        for i in range(0, len(self.delta)):
            for j in range(0, len(self.delta[i])):
                old = self.delta[i][j]
                if old not in self.summation and '0' in old:
                    if old[1:] not in changed:
                        if ar == "add":
                            new = str("0") + self.delta[i][j]
                            changed.append(old)
                        elif ar == "remove":
                            new = self.delta[i][j][1:] 
                            changed.append(new)
                        self.substitute(old, new)
                        changed.append(old)


    def rename_dummy(self, coef_dict,  ar):
        """ 
        Write appropriate priority number in front of summation indices" 
        """
        summation_dict = dict()
        for p in self.summation:
            summation_dict[p] = 0

        for p in self.summation:
            for i in range(0, len(self.coefficient_idx)):
                for j in self.coefficient_idx[i]:
                    if j == p:
#                        print('j', j, self.coefficient[i], coef_dict[self.coefficient[i]])
                        summation_dict[p] += coef_dict[self.coefficient[i]]
        for p in self.summation:
            new_p = str(summation_dict[p]) + str(p)
            self.substitute(p, new_p)

            
    def rename_dummy_back(self, ar):

        """ Removes redundant numbers arosen from rename_dummy, or
        from rename_dummy_opt if ar = 'opt'
        """
        
        numbers = set(['1','2','3','4','5','6','7','8','9', '0'])
        for p in range(0, len(self.summation)):
            i = 0
            for s in self.summation[p]:
                if s in numbers:
                    i += 1
                else:
                    break
            new_p = self.summation[p][i:]
            self.substitute(self.summation[p], new_p)
        if ar == 'opt' :
            for q in range(0, len(self.coefficient)):
                if self.coefficient[q] == TWOEL_INT:
                    for p in range(0, len(self.coefficient_idx[q])):
                        i = 0
                        for s in self.coefficient_idx[q][p]:
                            if s in numbers:
                                i += 1
                            else:
                                break
                        self.coefficient_idx[q][p] = self.coefficient_idx[q][p][i:]


    def sort_indices_fixed_priority(self, fixed_fx):
        """Sorts all indices alphabetically, but with fixed priority   
          according to permutation symmetry                                                                                               
           of each type of coefficient. Sort indices only for those coeffs
           that have fixed indices. 
        """

        
        indices = self.indices()
        fixed_set = indices - set(self.summation)
        if fixed_fx != []:
            fixed = fixed_fx
        else:
            fixed= list(fixed_set)
        for i in range(0, len(fixed)):
            fixed[i] = '0'+fixed[i]

        self.rename_fixed("add")

        for x in range(0, len(self.coefficient)):
            have_fixed = False
            for i in self.coefficient_idx[x]:
                if i in fixed:
                    have_fixed = True

            if have_fixed == True:
                self.standarize_selected(self.coefficient[x])


        self.rename_fixed("remove")
        
    def standarize_selected(self, coef):

        if coef == TWOEL_INT_TRANS:
            self.standarize_v()
        elif coef == EOM_TRIPLET_R3:
            self.standarize_triplet_R3()
        elif coef == EOM_TRIPLET_R2p:
            self.standarize_triplet_R2p()
        elif coef == EOM_TRIPLET_R2m:
            self.standarize_triplet_R2m()
        elif coef == TWOEL_INT: 
            self.standarize_g()
        elif coef == FOCK_MATRIX:
            self.standarize_f()
        elif coef == FOCK_MATRIX_TRANS:
            self.standarize_f()
        elif coef == TWOEL_INT_TRANS:
            self.standarize_w()
        elif coef == ANY_OP: 
            self.standarize_b()
        elif coef == CC_AMPLITUDE:
            self.standarize_t()
        elif coef == F12_AMPLITUDE:
            self.standarize_tf()
        elif coef == S_AMPLITUDE:
            self.standarize_s()
        elif coef == OBSERVABLE_X:
            self.standarize_sm()
        elif coef == OBSERVABLE_Y:
            self.standarize_sm()
        elif coef == OBSERVABLE_Z:
            self.standarize_sm()
        elif coef == BARENUCL_HAM:
            self.standarize_sm()

    def sort_indices(self):
        """Sorts all indices alphabetically
           according to permutation symmetry
           of each type of coefficient.
        """
        #
        # Sorting coefficient v - t1_transformed two electron integrals
        #

        self.standarize_v()

        #
        self.standarize_triplet_R3()
        self.standarize_triplet_R2p()
        self.standarize_triplet_R2m()
        #
        # Sorting coefficient g - two electron integrals
        #
#        print('au1')
        self.standarize_g()
        self.standarize_f()
        self.standarize_w()
        self.standarize_b()
        #
        # Sorting  coefficient t - for more than 1 pair of indices
        #
        self.standarize_t()
        self.standarize_tf()
        #
        # Sorting  coefficient s - for more than 1 pair of indices
        #
        self.standarize_s()
        #
        # Sort indices of coefficients belonging to one-electron
        # symmetric matrices
        #
        self.standarize_sm()
        self.standarize_delta()

    def standarize_delta(self):

        for x in self.delta:
            x.sort()
        self.delta.sort()


    def standarize_g(self):

#        print('self', self)
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == TWOEL_INT :
                pairs = [self.coefficient_idx[i][0:2],\
                             self.coefficient_idx[i][2:4]]
                #print(pairs)
                pairs[0].sort()
                pairs[1].sort()
                pairs.sort()
                self.coefficient_idx[i] = pairs[0] + pairs[1]

    def standarize_gm2(self):
        # print('')
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == DENS2 or self.coefficient[i] == DENS2P:
                # print('stand dens2', self)
                pairs_original = [[self.coefficient_idx[i][0], self.coefficient_idx[i][1]], [self.coefficient_idx[i][2], self.coefficient_idx[i][3]]]
                pairs = deepcopy(pairs_original)
                pairs[0].sort()
                pairs[1].sort()
                if pairs[0] != pairs_original[0]:
                    self.num_factor *= -1.0
                if pairs[1] != pairs_original[1]:
                    self.num_factor *= -1.0

                pairs.sort()
                self.coefficient_idx[i] = [pairs[0][0], pairs[0][1], pairs[1][0], pairs[1][1]]
#                print('plusio', self)
        #         print('stand dens2', self)
        # print('')

    def standarize_cm2(self):
        # print('')
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == CUM2:
                # print('stand dens2', self)
                pairs_original = [[self.coefficient_idx[i][0], self.coefficient_idx[i][1]], [self.coefficient_idx[i][2], self.coefficient_idx[i][3]]]
                pairs = deepcopy(pairs_original)
                pairs[0].sort()
                pairs[1].sort()
                if pairs[0] != pairs_original[0]:
                    self.num_factor *= -1.0
                if pairs[1] != pairs_original[1]:
                    self.num_factor *= -1.0
                self.coefficient_idx[i] = [pairs[0][0], pairs[0][1], pairs[1][0], pairs[1][1]]
        #         print('stand dens2', self)
        # print('')


    def standarize_gm234spin_new(self):
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == DENS2PM or self.coefficient[i] == DENS3PM or self.coefficient[i] == DENS4PM or self.coefficient[i] == DENS4PPM:
#                print('przed newspin', self)
                spl = len(self.coefficient_idx[i])//2
                coef_idx = self.coefficient_idx[i]
                spin_idx = self.coefficient_spin[i]
                new_coef_idx, num = canonicalize_tensor(coef_idx, spin_idx, self.coefficient[i], split=spl)

                self.coefficient_idx[i] = new_coef_idx

                self.num_factor *= num
 #               print('po newspin', self)


    def standarize_gm34spin(self):
        # print('przedspluszspin', self)
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == DENS3PM:
                # ++- ++-
                # 012 345

                ani_plus_list = self.coefficient_idx[i][0:2]
                cre_plus_list = self.coefficient_idx[i][3:5]
                ani_minus_list = [self.coefficient_idx[i][2]]
                cre_minus_list = [self.coefficient_idx[i][5]]

                
                
            elif self.coefficient[i] == DENS4PM:
                # ++-- ++--
                # 0123 4567

                ani_plus_list = self.coefficient_idx[i][0:2]
                ani_minus_list = self.coefficient_idx[i][2:4]
                cre_plus_list = self.coefficient_idx[i][4:6]
                cre_minus_list = self.coefficient_idx[i][6:8]

            elif self.coefficient[i] == DENS4PPM:
                # +++- +++-
                # 0123 4567

                ani_plus_list = self.coefficient_idx[i][0:3]
                cre_plus_list = self.coefficient_idx[i][4:7]
                ani_minus_list = [self.coefficient_idx[i][3]]
                cre_minus_list = [self.coefficient_idx[i][7]]


            if self.coefficient[i] == DENS3PM or self.coefficient[i] == DENS4PM or self.coefficient[i] == DENS4PPM:
                ani_plus_list_sorted = deepcopy(ani_plus_list)
                cre_plus_list_sorted = deepcopy(cre_plus_list)

                ani_minus_list_sorted = deepcopy(ani_minus_list)
                cre_minus_list_sorted = deepcopy(cre_minus_list)

                ani_plus_list_sorted.sort()
                cre_plus_list_sorted.sort()

                ani_minus_list_sorted.sort()
                cre_minus_list_sorted.sort()

            
                nswaps_plus_cre = nswaps_in_two_lists(deepcopy(cre_plus_list), deepcopy(cre_plus_list_sorted))
                nswaps_plus_ani = nswaps_in_two_lists(deepcopy(ani_plus_list), deepcopy(ani_plus_list_sorted))

                nswaps_minus_cre = nswaps_in_two_lists(deepcopy(cre_minus_list), deepcopy(cre_minus_list_sorted))
                nswaps_minus_ani = nswaps_in_two_lists(deepcopy(ani_minus_list), deepcopy(ani_minus_list_sorted))


                nswaps = nswaps_plus_cre +  nswaps_plus_ani + nswaps_minus_cre + nswaps_minus_ani
                num = (-1.0)**nswaps

                self.num_factor *= num
                self.coefficient_idx[i] = ani_plus_list_sorted + ani_minus_list_sorted + cre_plus_list_sorted+ cre_minus_list_sorted
                # print('po', self)
                # print()

                
                
    def standarize_gm34(self):
        # print('')
        # print('stand dens34', self)
        for i in range(0, len(self.coefficient)):
            # print(self.coefficient[i])
            if self.coefficient[i] == DENS3 or self.coefficient[i] == DENS4:
                # print('self.coefficient_idx[i]', self.coefficient_idx[i])
                pairs_original = []
                # wersja I
                # for j in range(0, len(self.coefficient_idx[i]), 2):
                #     pairs_original.append([self.coefficient_idx[i][j], self.coefficient_idx[i][j+1]])
                # wersja II zeby spiny sie zmiescily

                lencoef = len(self.coefficient_idx[i])
                ani_list = self.coefficient_idx[i][0:int(lencoef/2)]
                cre_list = self.coefficient_idx[i][int(lencoef/2):lencoef]

                ani_list_sorted = deepcopy(ani_list)
                cre_list_sorted = deepcopy(cre_list)

                ani_list_sorted.sort()
                cre_list_sorted.sort()
#                print('anicre', self)

                nswaps_ani = nswaps_in_two_lists(deepcopy(ani_list), deepcopy(ani_list_sorted))
                nswaps_cre = nswaps_in_two_lists(deepcopy(cre_list), deepcopy(cre_list_sorted))

                num_ani = 1.0
                num_cre = 1.0
                if (nswaps_ani%2) == 1:
                    num_ani = -1.0
                if (nswaps_cre%2) ==1:
                    num_cre = -1.0

                self.num_factor *= num_ani*num_cre
                self.coefficient_idx[i] = ani_list_sorted + cre_list_sorted
#                print('anicrepo', self)
                #----------------------------------------------------------------------------------------
                # This is incorrect part of the code were the wrong symmetries are tested
                # when we swap pairs of creators and anihilator pairs together
                # the sign is always plus, like in CC amplitudes
                # I leave this part for now only for testing
                # and to remind myself how stupid I am

                # lenself = len(self.coefficient_idx[i])
                # for j in range(0, lenself//2):
                #     pairs_original.append([self.coefficient_idx[i][j], self.coefficient_idx[i][lenself-1-j]])


                # # compute n of swaps
                # nswaps = nswaps_in_sort(deepcopy(pairs_original))
                # print('pairs original', pairs_original)
                # pairs_sorted = deepcopy(pairs_original)
                # pairs_sorted.sort()
                # print('pairs_sorted', pairs_sorted)

                # print(self, 'nswaps', nswaps)
                # # if (nswaps%2) == 1:
                # #     self.num_factor *= -1.0
                # # print('pairs_sorted', pairs_sorted)
                # new_idx_list = [None] * lenself

                # for j in range(0, lenself//2):
                #     new_idx_list[j] = pairs_sorted[j][0]
                #     new_idx_list[lenself-1-j] = pairs_sorted[j][1]
                # # for pair in pairs_sorted:
                # #     for elem in pair:
                # #         # print('dodaje', elem)
                # #         new_idx_list.append(elem)
                # self.coefficient_idx[i] = new_idx_list
                # print(self, 'nswaps', nswaps)
                # print('')
                #----------------------------------------------------------------------------------------
                # print('222self.coefficient_idx[i]', self.coefficient_idx[i])
                                        
        # print('poand dens34', self)
        # print('')



    def standarize_gs(self):

        # print('')
        # print('self w stand gs', self)
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == TWOEL_INT_DIRAC:
                pairs = [[self.coefficient_idx[i][0], self.coefficient_idx[i][2]], [self.coefficient_idx[i][1], self.coefficient_idx[i][3]]]
                # print(pairs)
                pairs[0].sort()
                pairs[1].sort()
                # print(pairs)
                pairs.sort()
                # print(pairs)
                self.coefficient_idx[i] = [pairs[0][0], pairs[1][0], pairs[0][1], pairs[1][1]]

                
    def standarize_f(self):        
#        print('standf', self)
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == F12_TWOEL :
                pairs = [self.coefficient_idx[i][0:2],\
                             self.coefficient_idx[i][2:4]]
#                print('pairs', pairs)
                pairs[0].sort()
                pairs[1].sort()
                pairs.sort()
                self.coefficient_idx[i] = pairs[0] + pairs[1]
#        print('standf2', self)

    # def standarize_f(self):
    #     #                                                                                                                                                            
    #     # Sort indices of F12 AMPLITUDE pairwise in alphabetical order                                                                                             
    #     #                                                                                                                                                             
    #     for i in range(0, len(self.coefficient)):
    #         print('standf', self)
    #         if self.coefficient[i] == F12_AMPLITUDE:
    #             print(self.coefficient_idx[i])
    #             pairs = []
    #             pairs.append([self.coefficient_idx[i][0], self.coefficient_idx[i][2]])
    #             pairs.append([self.coefficient_idx[i][1], self.coefficient_idx[i][3]])
    #             pairs_sort = deepcopy(pairs)
    #             pairs_sort[0].sort()
    #             pairs_sort[1].sort()
    #             print(pairs)
    #             print(pairs_sort)
    #             nf1 = 1
    #             nf2 = 1
    #             if pairs_sort[0] != pairs[0]:
    #                 nf1 = -1
    #             if  pairs_sort[1] != pairs[1]:
    #                 nf2 = -1
    #             nf = nf1 * nf2
    #             self.num_factor *= nf
    #             self.coefficient_idx[i] = []
    #             self.coefficient_idx[i].append(pairs_sort[0][0])
    #             self.coefficient_idx[i].append(pairs_sort[1][0])
    #             self.coefficient_idx[i].append(pairs_sort[0][1])
    #             self.coefficient_idx[i].append(pairs_sort[1][1])
    #         print('standf2', self)



    def standarize_b(self):

        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == SLAT2_SYM1:
                pairs = []
                for y in range(0, len(self.coefficient_idx[i]), 2):
                    pairs.append([self.coefficient_idx[i][y], self.coefficient_idx[i][y+1]])
                for y in range(0, len(pairs)):
                    pairs[y].sort()
                pairs.sort()
                
                self.coefficient_idx[i] = []
                for y in range(0, len(pairs)):
                    self.coefficient_idx[i].append(pairs[y][0])
                    self.coefficient_idx[i].append(pairs[y][1])

        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == SLAT2_SYM2 or self.coefficient[i] == SLAT3:
                pairs = []
                for y in range(0, len(self.coefficient_idx[i]), 2):
                    pairs.append([self.coefficient_idx[i][y], self.coefficient_idx[i][y+1]])
                for y in range(0, len(pairs)):
                    pairs[y].sort()

                self.coefficient_idx[i] = []
                for y in range(0, len(pairs)):
                    self.coefficient_idx[i].append(pairs[y][0])
                    self.coefficient_idx[i].append(pairs[y][1])

        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == SLAT4:
                pairs = []
                for y in range(0, len(self.coefficient_idx[i]), 2):
                    pairs.append([self.coefficient_idx[i][y], self.coefficient_idx[i][y+1]])
                for y in range(0, len(pairs)):
                    pairs[y].sort()
                pairs1 = pairs[0:2]
                pairs2 = pairs[2:4]
                pairs1.sort()
                pairs2.sort()
                pairs0 = [pairs1[0], pairs1[1], pairs2[0], pairs2[1]]
 
                self.coefficient_idx[i] = []
                for y in range(0, len(pairs0)):
                    self.coefficient_idx[i].append(pairs[y][0])
                    self.coefficient_idx[i].append(pairs[y][1])


    def standarize_w(self):
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == TWOEL_INT_AS:
                pairs = [self.coefficient_idx[i][0:2],\
                             self.coefficient_idx[i][2:4]]
                pairs_start = deepcopy(pairs)
                pairs[0].sort()
                pairs[1].sort()
                if pairs[0] != pairs_start[0] and pairs[1] == pairs_start[1]:
                    self.num_factor *= -1.0
                elif pairs[0] == pairs_start[0] and pairs[1] != pairs_start[1]:
                    self.num_factor *= -1.0
                pairs.sort()
                self.coefficient_idx[i] = pairs[0] + pairs[1]

    def standarize_v(self):
        # print(self)
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == TWOEL_INT_TRANS or self.coefficient[i] == TWOEL_INT_COMBO:
                pairs = [self.coefficient_idx[i][0:2],\
                             self.coefficient_idx[i][2:4]]
                pairs.sort()
                self.coefficient_idx[i] = pairs[0] + pairs[1]


    def standarize_g_for_fortran(self):
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == TWOEL_INT :
                pairs = [self.coefficient_idx[i][0:2],\
                             self.coefficient_idx[i][2:4]]
                pairs[0].sort()
                pairs[1].sort()
                ds01 = 0
                ds02 = 0
                for j in pairs[0]:
                    for k in j:
                        if k == '0':
                            ds01 += 1
                        elif k in virtual:
                            ds02 += 1
                ds11 = 0
                ds12 = 0
                for j in pairs[1]:
                    for k in j:
                        if k == '0':
                            ds11 += 1
                        elif k in virtual:
                            ds12 += 1
                # print(pairs)

                if (ds01 + ds02) >= (ds11 + ds12) and ds02 > ds12:
                    self.coefficient_idx[i] = pairs[0] + pairs[1]
                else:
                    self.coefficient_idx[i] = pairs[1] + pairs[0]


    def standarize_v_for_fortran(self):
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == TWOEL_INT_TRANS :
                pairs = [self.coefficient_idx[i][0:2],\
                             self.coefficient_idx[i][2:4]]

                ds01 = 0
                ds02 = 0
                for j in pairs[0]:
                    for k in j:
                        if k == '0':
                            ds01 += 1
                        elif k in virtual:
                            ds02 += 1
                ds11 = 0
                ds12 = 0
                for j in pairs[1]:
                    for k in j:
                        if k == '0':
                            ds11 += 1
                        elif k in virtual:
                            ds12 += 1
                            
                if (ds01 + ds02) >= (ds11 + ds12) and ds02 > ds12:
                    self.coefficient_idx[i] = pairs[0] + pairs[1]
                elif (ds11 + ds12) >= (ds01 + ds02) and ds02 > ds12:
                    self.coefficient_idx[i] = pairs[0] + pairs[1]
                else:
                    self.coefficient_idx[i] = pairs[1] + pairs[0]        
                    
    def standarize_TT(self):
        #
        # Sort indices fronm Torun code
        # only upper lower indices change possible
        #
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i][0:2] == 'TT':
                off = int(len(self.coefficient_idx[i])/2)
                upperi = self.coefficient_idx[i][0:off]
                loweri = self.coefficient_idx[i][off:len(self.coefficient_idx[i])]
                sort_idx = [upperi, loweri]
                sort_idx.sort()
                self.coefficient_idx[i] = []
                for k in sort_idx:
                    for j in k:
                        self.coefficient_idx[i].append(j)


    def standarize_t(self):
        #
        # Sort indices of CC amplitudes pairwise in alphabetical order
        #

        indices = self.indices()
        fixed_set = indices - set(self.summation)
        fixed= list(fixed_set)

        for i in range(0, len(self.coefficient)):
            last_idx_base = 'zzz'
            if self.coefficient[i] == CC_AMPLITUDE:
                if len(self.coefficient_idx[i]) > 2:
                    pairs = []
                    pairscp = []
                    for j in range(0, len(self.coefficient_idx[i]) // 2):
                        p = self.coefficient_idx[i][2 * j]
                        q = self.coefficient_idx[i][2 * j + 1]
                        last_idx = last_idx_base + '{j}'.format(j=j)
                        pairs.append([p, q])
                        if p in fixed:
                            p0 = '0' + p
                        else:
                            p0 = p
                        if q in fixed:
                            q0 = '0' + q
                        else:
                            q0 = q
                        pairscp.append([p0, q0, last_idx])
                    for j in range(0, len(pairscp)):
                        pairscp[j].sort()
                    pairscp.sort()
                    sorted_list = []
#                    print('prprprpr', pairscp)
                    for j in pairscp:
                        li = j[2][3:] #'zzz1' or 'zzz2', or 'zzz0' or...
                        li = int(li)
                        sorted_list.append(deepcopy(li))
#                    print('sorteeeeeed', sorted_list)
                    self.coefficient_idx[i] = []
                    for m in sorted_list:
                        for k in range(0, 2):
                            self.coefficient_idx[i].append(pairs[m][k])
#                    print('palace', pairs)
                    # pairs.sort()
                    # print(pairs)
                    # self.coefficient_idx[i] = []
                    # for p, q in pairs:
                    #     self.coefficient_idx[i].append(p)
                    #     self.coefficient_idx[i].append(q)

    def standarize_tf(self):
        indices = self.indices()
        fixed_set = indices - set(self.summation)
        fixed= list(fixed_set)

        for i in range(0, len(self.coefficient)):
            last_idx_base = 'zzz'
            if self.coefficient[i] == F12_AMPLITUDE:
                if len(self.coefficient_idx[i]) > 2:
                    pairs = []
                    pairscp = []
                    for j in range(0, len(self.coefficient_idx[i]) // 2):
                        p = self.coefficient_idx[i][2 * j]
                        q = self.coefficient_idx[i][2 * j + 1]
                        last_idx = last_idx_base + '{j}'.format(j=j)
                        pairs.append([p, q])
                        if p in fixed:
                            p0 = '0' + p
                        else:
                            p0 = p
                        if q in fixed:
                            q0 = '0' + q
                        else:
                            q0 = q
                        pairscp.append([p0, q0, last_idx])
                    for j in range(0, len(pairscp)):
                        pairscp[j].sort()
                    pairscp.sort()
                    sorted_list = []
                    # print('prprprpr', pairscp)
                    for j in pairscp:
                        li = j[2][3:] #'zzz1' or 'zzz2', or 'zzz0' or...                                                                                                                                                             
                        li = int(li)
                        sorted_list.append(deepcopy(li))
                    # print('sorteeeeeed', sorted_list)
                    self.coefficient_idx[i] = []
                    for m in sorted_list:
                        for k in range(0, 2):
                            self.coefficient_idx[i].append(pairs[m][k])
                    # print('palace', pairs)



#     def standarize_tf(self):
#         #                                                                                                                                                          
#         # Sort indices of F12 amplitudes pairwise in alphabetical order                                                                               
#         #          
#         # print('foto', self)
#         for i in range(0, len(self.coefficient)):
#             if self.coefficient[i] == F12_AMPLITUDE:
# #                print(self.coefficient_idx[i])
#                 pairs = []
#                 pairs.append([self.coefficient_idx[i][0], self.coefficient_idx[i][2]])
#                 pairs.append([self.coefficient_idx[i][1], self.coefficient_idx[i][3]])
#                 pairs_sort = deepcopy(pairs)
#                 pairs_sort[0].sort()
#                 pairs_sort[1].sort()
#                 # print(pairs)
#                 # print(pairs_sort)
#                 nf1 = 1
#                 nf2 = 1
#                 if pairs_sort[0] != pairs[0]:
#                     nf1 = -1
#                 if  pairs_sort[1] != pairs[1]:
#                     nf2 = -1
#                 nf = nf1 * nf2
#                 self.num_factor *= nf
#                 self.coefficient_idx[i] = []
#                 self.coefficient_idx[i].append(pairs_sort[0][0])
#                 self.coefficient_idx[i].append(pairs_sort[1][0])
#                 self.coefficient_idx[i].append(pairs_sort[0][1])
#                 self.coefficient_idx[i].append(pairs_sort[1][1])

                
#        print('tofo', self)

    def standarize_s(self):
        #
        # Sort indices of CC amplitudes pairwise in alphabetical order
        #
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == S_AMPLITUDE:
                if len(self.coefficient_idx[i]) > 2:
                    pairs = []
                    for j in range(0, len(self.coefficient_idx[i]) // 2):
                        p = self.coefficient_idx[i][2 * j]
                        q = self.coefficient_idx[i][2 * j + 1]
                        pairs.append((p, q))
                    pairs.sort()
                    self.coefficient_idx[i] = []
                    for p, q in pairs:
                        self.coefficient_idx[i].append(p)
                        self.coefficient_idx[i].append(q)

    def standarize_rl(self):
        #
        # Sort indices of CC amplitudes pairwise in alphabetical order
        #
        # print('przed rl', self)
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == EOM_CC_AMPLITUDE_R or self.coefficient[i] == EOM_CC_AMPLITUDE_L \
               or self.coefficient[i] == EOM_CC_SINGLE_Rr or self.coefficient[i] == EOM_CC_SINGLE_Rl \
               or self.coefficient[i] == EOM_CC_SINGLE_Rr_minus or self.coefficient[i] == EOM_CC_SINGLE_Rl_minus:
                # print('tatttttttk')
                if len(self.coefficient_idx[i]) > 2:
                    pairs = []
                    for j in range(0, len(self.coefficient_idx[i]) // 2):
                        p = self.coefficient_idx[i][2 * j]
                        q = self.coefficient_idx[i][2 * j + 1]
                        pairs.append((p, q))
                    pairs.sort()
                    self.coefficient_idx[i] = []
                    for p, q in pairs:
                        self.coefficient_idx[i].append(p)
                        self.coefficient_idx[i].append(q)
        # print('po rl', self)


    def standarize_triplet_R3(self):
        #
        # Sort indices of EOM_TRIPLET
        # according to the following symmetry rules
        # in alphabetical order
        # R_(ai)bjck = R_(ai)ckbj = - R_(ai)bkcj = - R_(ai)cjbk
        #
        #


       for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == EOM_TRIPLET_R3:


                pairs_virt = [self.coefficient_idx[i][2], self.coefficient_idx[i][4]]
                pairs_occ = [self.coefficient_idx[i][3], self.coefficient_idx[i][5]]

                pairs_virt_sort = deepcopy(pairs_virt)
                pairs_virt_sort.sort()
                pairs_occ_sort = deepcopy(pairs_occ)
                pairs_occ_sort.sort()

                if pairs_virt_sort != pairs_virt:
                    self.num_factor = - self.num_factor
                    self.coefficient_idx[i][2] = pairs_virt_sort[0]
                    self.coefficient_idx[i][4] = pairs_virt_sort[1]


                if pairs_occ_sort != pairs_occ:
                    self.num_factor = - self.num_factor
                    self.coefficient_idx[i][3] = pairs_occ_sort[0]
                    self.coefficient_idx[i][5] = pairs_occ_sort[1]



        # for i in range(0, len(self.coefficient)):
        #     if self.coefficient[i] == EOM_TRIPLET_R3:
        #         pairs = []

        #         for j in range(1, len(self.coefficient_idx[i]) // 2):
        #             p = self.coefficient_idx[i][2 * j]
        #             q = self.coefficient_idx[i][2 * j + 1]
        #             pairs.append([p, q])
        #         pairs.sort()

        #         l1 = [pairs[0][1], pairs[1][1]]
        #         l2 = deepcopy(l1)
        #         l2.sort()

        #         if l2 != l1:
        #             pairs_copy = deepcopy(pairs)
        #             pairs_copy[0][1] = pairs[1][1]
        #             pairs_copy[1][1] = pairs[0][1]
        #             self.num_factor = - self.num_factor
        #             pairs = deepcopy(pairs_copy)
                
                
        #         self.coefficient_idx[i][2] = pairs[0][0]
        #         self.coefficient_idx[i][3] = pairs[0][1]
        #         self.coefficient_idx[i][4] = pairs[1][0]
        #         self.coefficient_idx[i][5] = pairs[1][1]



    def standarize_triplet_R2p(self):
        #
        # Sort indices of EOM_TRIPLET
        # according to the following symmetry rules
        # in alphabetical order
        # R_bjck = R_ckbj = - R_bkcj = - R_cjbk
        #
        #

        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == EOM_TRIPLET_R2p:

                pairs_virt = [self.coefficient_idx[i][0], self.coefficient_idx[i][2]]
                pairs_occ = [self.coefficient_idx[i][1], self.coefficient_idx[i][3]]

                pairs_virt_sort = deepcopy(pairs_virt)
                pairs_virt_sort.sort()
                pairs_occ_sort = deepcopy(pairs_occ)
                pairs_occ_sort.sort()

                if pairs_virt_sort != pairs_virt:
                    self.num_factor = - self.num_factor
                    self.coefficient_idx[i][0] = pairs_virt_sort[0]
                    self.coefficient_idx[i][2] = pairs_virt_sort[1]


                if pairs_occ_sort != pairs_occ:
                    self.num_factor = - self.num_factor
                    self.coefficient_idx[i][1] = pairs_occ_sort[0]
                    self.coefficient_idx[i][3] = pairs_occ_sort[1]


    def standarize_triplet_R2m(self):
        #
        # Sort indices of EOM_TRIPLET
        # according to the following symmetry rules
        # in alphabetical order
        # R_bjck = R_ckbj = - R_bkcj = - R_cjbk
        #
        #
        
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == EOM_TRIPLET_R2m:


                pairs_virt = [self.coefficient_idx[i][0], self.coefficient_idx[i][2]]
                pairs_occ = [self.coefficient_idx[i][1], self.coefficient_idx[i][3]]

                # [d, b]  --> [b, d]    tak (zamien oba)  [
                # [j, k]  --> [j, k]

                # [d, b]  --> [b, d]    tak (zamien oba)
                # [k, j]  --> [j, k]


                # [b, d]  --> [b, d]    nie zamieniaj
                # [j, k]  --> [j, k]


                # [b, d]  --> [b, d]    nie zamieniaj
                # [k, j]  --> [j, k]


                pairs_virt_sort = deepcopy(pairs_virt)
                pairs_virt_sort.sort()
                pairs_occ_sort = deepcopy(pairs_occ)
                pairs_occ_sort.sort()


                if pairs_virt_sort != pairs_virt and pairs_occ_sort != pairs_occ:
                    self.num_factor = - self.num_factor
                    self.coefficient_idx[i][0] = pairs_virt_sort[0]
                    self.coefficient_idx[i][2] = pairs_virt_sort[1]
                    self.coefficient_idx[i][1] = pairs_occ_sort[0]
                    self.coefficient_idx[i][3] = pairs_occ_sort[1]

                if pairs_virt_sort != pairs_virt and pairs_occ_sort == pairs_occ:
                    self.num_factor = - self.num_factor
                    self.coefficient_idx[i][0] = pairs_virt_sort[0]
                    self.coefficient_idx[i][2] = pairs_virt_sort[1]
                    self.coefficient_idx[i][1] = pairs_occ[1]
                    self.coefficient_idx[i][3] = pairs_occ[0]


                # pairs = []

                # for j in range(0, len(self.coefficient_idx[i]) // 2):
                #     p = self.coefficient_idx[i][2 * j]
                #     q = self.coefficient_idx[i][2 * j + 1]
                #     pairs.append([p, q])
                # print('pairs', pairs)
                # pairs.sort()
                # print('pairs', pairs)
                # l1 = [pairs[0][1], pairs[1][1]]
                # l2 = deepcopy(l1)
                # l2.sort()

                # print('l1', l1)
                # print('l2', l2)

                # if l2 != l1:
                #     print('l1!=l2')
                #     pairs_copy = deepcopy(pairs)
                #     pairs_copy[0][1] = pairs[1][1]
                #     pairs_copy[1][1] = pairs[0][1]
                #     self.num_factor = - self.num_factor
                #     pairs = deepcopy(pairs_copy)
                #     print(self.num_factor)
                #     print(pairs)
                
                
                # self.coefficient_idx[i][0] = pairs[0][0]
                # self.coefficient_idx[i][1] = pairs[0][1]
                # self.coefficient_idx[i][2] = pairs[1][0]
                # self.coefficient_idx[i][3] = pairs[1][1]


    def standarize_sm(self):
        #
        # Sort indices of coefficients belonging to SYMMETRIC_MATRIX
        # subset
        #
        for i in range(0, len(self.coefficient)):
            name = self.coefficient[i]
            if name in SYMMETRIC_MATRIX:
                if name!=OBSERVABLE_X:
                    self.coefficient_idx[i].sort()

    def sort_nosort_coeffs(self, last_idx = 'dup', cas=False):
        """
        Sort_nosort coefficients according to:
        [COEFF. NAME, NO. OF INDICES, FIXED INDICES, ALL INDICES]
        Always put two electron integrals of any kind last.
        """
        # print('')
        # print('')
        # print('-------------------------------', cas)
        #        print('lastlast', last_idx)
        twoelset = [TWOEL_INT, TWOEL_INT_COMBO, TWOEL_INT_TRANS, F12_TWOEL, TWOEL_INT_DIRAC]

        if len(self.coefficient) > 0:
            priority_list = []
            indices = self.indices()
            fixed = indices - set(self.summation)

            for i in range(0, len(self.coefficient)):
                name_true = deepcopy(self.coefficient[i])

                #name = deepcopy(self.coefficient[i])

                name = self.coefficient[i]
                if name_true in twoelset:
                    name = 'z'+ name

                if name_true == last_idx:
                    name = 'zz'+name
                nidx = len(self.coefficient_idx[i])
                fidx = 0
                flst = []
                for p in self.coefficient_idx[i]:
                    if p in fixed:
                        fidx += 1
                        flst.append(p)

                idx = self.coefficient_idx[i]
                sidx = self.coefficient_spin[i]
                
                if name_true == last_idx:
                    priority_list.append([100, flst, name, nidx, idx, name_true, sidx])
                else:
                    if cas:
                        priority_list.append([ name, 100-fidx, flst, nidx, idx, name_true, sidx])

                    else:
                        priority_list.append([100-fidx, flst, name, nidx, idx, name_true, sidx])

            priority_list.sort()

            #self.coefficient = []
            coefficient_idx = []
            #self.coefficient_spin = []
            for c in priority_list:
                coefficient_idx.append(c[4])
        refine_idx = []
        if len(self.coefficient) > 0:
        
            for coef_idx_lst in coefficient_idx:
                for idx in coef_idx_lst:
                    if idx in self.summation:
                        if idx not in refine_idx:
                            refine_idx.append(idx)

        return refine_idx

    def sort_coeffs(self, last_idx = 'dup', cas=False):
        """
        Sort coefficients according to:
        [COEFF. NAME, NO. OF INDICES, FIXED INDICES, ALL INDICES]
        Always put two electron integrals of any kind last.
        """
        # print('')
        # print('')
        # print('-------------------------------', cas)
        #        print('lastlast', last_idx)
        twoelset = [TWOEL_INT, TWOEL_INT_COMBO, TWOEL_INT_TRANS, F12_TWOEL, TWOEL_INT_DIRAC]

        if len(self.coefficient) > 0:
            priority_list = []
            indices = self.indices()
            fixed = indices - set(self.summation)

            for i in range(0, len(self.coefficient)):
                name_true = deepcopy(self.coefficient[i])

                name = deepcopy(self.coefficient[i])
                if name_true in twoelset:
                    name = 'z'+ name

                if name_true == last_idx:
                    name = 'zz'+name
                nidx = len(self.coefficient_idx[i])
                fidx = 0
                flst = []
                for p in self.coefficient_idx[i]:
                    if p in fixed:
                        fidx += 1
                        flst.append(p)

                idx = self.coefficient_idx[i]
                #print(self)
                #print(self.coefficient_spin)
                if len(self.coefficient_spin) > 0:
                    sidx = self.coefficient_spin[i]
                else:
                    sidx = []
                
                if name_true == last_idx:
                    priority_list.append([100, flst, name, nidx, idx, name_true, sidx])
                else:
                    if cas:
                        priority_list.append([ name, 100-fidx, flst, nidx, idx, name_true, sidx])

                    else:
                        priority_list.append([100-fidx, flst, name, nidx, idx, name_true, sidx])

            priority_list.sort()

            self.coefficient = []
            self.coefficient_idx = []
            for c in priority_list:
                self.coefficient.append(c[5])
                self.coefficient_idx.append(c[4])

            if len(self.coefficient_spin) >0:
                self.coefficient_spin = []
            for c in priority_list:
                self.coefficient_spin.append(c[6])
                
        self.delta.sort()
        self.summation.sort()

        

    def reallocate_dummy(self, cas=False):
        """
        Reallocate summation indices to canonical form. The set
        of summation indices remains unchanged after reallocation.
        """
        #print('')
        #print('')
        #print('reallo1', self)
        all_indices = self.indices()
        free_indices = set(self.summation)
        fixed_indices = all_indices - free_indices

        #print('--------------------------------------------------------')
        #print(all_indices)
        #print(free_indices)
        #print(fixed_indices)
        
        opools = set(occupied) & free_indices
        vpools = set(virtual) & free_indices
        gpools = set(general) & free_indices
        cvpools = set(completev) & free_indices
        cpools = set(complete) & free_indices
        cabspool = set(CABS) & free_indices

        
        opool = list(opools)
        vpool = list(vpools)
        gpool = list(gpools)
        cpool = list(cpools)
        cvpool = list(cvpools)
        cabspool = list(cabspool)


        opool.sort()
        vpool.sort()
        gpool.sort()
        cpool.sort()
        cvpool.sort()
        cabspool.sort()

        pairs = self.fixed_free_pairs(fixed_indices, free_indices)
        #print('pairs', pairs)
        pairs = self.clean_pairs_from_equivalent(pairs)
        ref_o = []
        ref_v = []
        ref_g = []
        ref_c = []
        ref_cv = []
        new_o = []
        new_v = []
        new_g = []
        new_c = []
        new_cv = []
        #print('reallo2', self)
        #print('NEW', new_v)
        for p in sorted(pairs.keys()):
            most_common = pairs[p].most_common()
            #print('mostcommon', most_common)
            for k in range(0, len(most_common)):
                #print('kkkk', k)
                isunique = True
                if k > 0:
                    if most_common[k - 1][1] == most_common[k][1]:
                        isunique = False
                if k < len(most_common) - 1:
                    if most_common[k + 1][1] == most_common[k][1]:
                        isunique = False
                if isunique:
                    #print('TR', isunique)
                    #print(most_common)
                    refidx = most_common[k][0]               
                    #print('refidx', refidx)
                    if refidx in occupied:
                        if refidx not in ref_o:
                            new_o.append(free_idx(opool, new_o))
                            ref_o.append(refidx)
                            #print('zamieniam to na najczestszy', ref_o, new_o)
                    elif refidx in virtual:
                        if refidx not in ref_v:
                            #print('new_v przed', new_v)
                            new_v.append(free_idx(vpool, new_v))
                            ref_v.append(refidx)
                            #print('new_v po', new_v)
                            #print('zamieniam to na najczestszy', ref_v, new_v)
                    elif refidx in completev:
                        if refidx not in ref_cv:
                            new_cv.append(free_idx(cvpool, new_cv))
                            ref_cv.append(refidx)
                            #print('zamieniam to na najczestszy', ref_cv, new_cv)
                    elif refidx in complete:
                        if refidx not in ref_c:
                            new_c.append(free_idx(cpool, new_c))
                            ref_c.append(refidx)
                    elif refidx in CABS:
                        if refidx not in ref_c:
                            new_c.append(free_idx(cabspool, new_c))
                            ref_c.append(refidx)
                    else:
                        if refidx not in ref_g:
                            new_g.append(free_idx(gpool, new_g))
                            ref_g.append(refidx)
        # #print('reallo3', self)
       # #print('NEW2', new_v)

        idx_list = self.coefficient_idx + self.operator_idx + self.delta
        #print('reallo4', self)
        #print('idx_list', idx_list)

        for j in idx_list:
            for i in j:
                #print(i, free_indices)
                if i in free_indices:
                    if i in virtual:
                        #print('h1', i)
                        if i not in ref_v:
                            new_v.append(free_idx(vpool, new_v))
                            ref_v.append(i)
                    if i in occupied:
                        #print('h2')
                        if i not in ref_o:
                            new_o.append(free_idx(opool, new_o))
                            ref_o.append(i)        
                    if i in completev:
                        #print('h3', ref_cv)
                        if i not in ref_cv:
                            #print('ano')
                            new_cv.append(free_idx(cvpool, new_cv))
                            ref_cv.append(i)
                    if i in complete:
#                        # #print('h4')
                        if i not in ref_c:
                            new_c.append(free_idx(cpool, new_c))
                            ref_c.append(i)
                    if i in CABS:
                        if i not in ref_c:
                            new_c.append(free_idx(cabspool, new_c))
                            ref_c.append(i)

                    if i in general:
                        if i not in ref_g:
                            new_g.append(free_idx(gpool, new_g))
                            ref_g.append(i)
                            # #print(new_v, new_o, new_g, new_c, new_cv)
        #print(ref_o, new_o)
        #print('reallo5', self, ref_v, new_v)
        self.multisubst(ref_v, new_v)
        #print('po_multisubst', self)
        # to sie okazalo niepotrzebne bo wystarczy dwa razy zrobic symplifikacje
        self.standarize_t()
        #print('po stand t', self)
        self.sort_coeffs('', cas)
        #print('po_ sort coef', self)                
        #print('po_v', self)
        self.multisubst(ref_o, new_o)
        #print('po_o', self)
        self.multisubst(ref_g, new_g)
        #print(ref_c, new_c)
        self.multisubst(ref_c, new_c)
        #print('reallo6', self)
        self.multisubst(ref_cv, new_cv)
        #print('reallo7', self)
        self.summation.sort()



    def move_summation_to_the_left(self):
        
        for i in range(0, len(self.coefficient)):
            for j in range(0, len(self.coefficient_idx[i])):
                if self.coefficient_idx[i][j] in self.summation:
                    self.coefficient_idx[i][j] = "0{strg}".format(strg = self.coefficient_idx[i][j])

        self.sort_indices()
#        # #print('au2')
        self.standarize_g()
        self.standarize_t()
        self.standarize_tf()
        self.standarize_s()
        self.standarize_v()
        self.standarize_delta()
        self.standarize_rl()
#        self.standarize_r_plus_minus()
        
        for i in range(0, len(self.coefficient)):
            for j in range(0, len(self.coefficient_idx[i])):
                self.coefficient_idx[i][j] = self.coefficient_idx[i][j].replace('0', '')


    def standarize(self, last_idx = '', fixed_fx = [], cas=False, drag_fixed=[]):
        """Arrange indices in self.summation, self.coefficient, self.coefficient_idx,
        self.delta into standard form.
        """
        #
        # Create dictionary of differend values for different
        # coefficients
        #
        # print('self przed standarize', cas, self)
        # t0 = time.time()
        self.clear_fixed()
        self.establish_fixed()
        # t1 = time.time()
        coef_dict = dict()
        k = 1
        for i in range(0, len(self.coefficient)):
            if i not in list(coef_dict.keys()):
                coef_dict[self.coefficient[i]] = 10**k
                k += 2

        # t2 = time.time()
        # print()
        # print('two time standarize,', t2-t1)
        # print()

        # t1 = time.time()
        self.sort_coeffs(last_idx, cas)
        # t2 = time.time()
        # print()
        # print('sor_coefs time standarize,', t2-t1)
        # print()
        # t1 = time.time()

        self.sort_indices_fixed_priority(fixed_fx)
        # t2 = time.time()
        # print()
        # print('sort_fixed time standarize,', t2-t1)
        # print()
        # t1 = time.time()
        self.sort_indices()
        # t2 = time.time()
        # print()
        # print('sort indices time standarize,', t2-t1)
        # print()
#        print('przed refine', self, self.orbital_type)
        # t1 = time.time()
        self.refine_summation(last_idx, cas, drag_fixed)
#        print('po refine', self, self.orbital_type)
        # t2 = time.time()
        # print()
        # print('refine summ time standarize,', t2-t1)
        # print()
        # t1 = time.time()
        self.rename_fixed("add")
        # t2 = time.time()
        # print()
        # print('rename fixed time standarize,', t2-t1)
        # print()
        # t1 = time.time()
        self.rename_dummy(coef_dict, "add")
        # t2 = time.time()
        # print()
        # print('rename dummy time standarize,', t2-t1)
        # print()
        # t1 = time.time()
        # t2 = time.time()
        # print()
        # print('three time standarize,', t2-t1)
        # print()

        self.sort_indices()
        # t2 = time.time()
        # print()
        # print('three1 time standarize,', t2-t1)
        # print()
        # t1 = time.time()
        self.rename_dummy_back('usual')
        # t2 = time.time()
        # print()
        # print('three2 time standarize,', t2-t1)
        # print()
        # t1 = time.time()
        self.rename_fixed("remove")
        # t2 = time.time()
        # print()
        # print('three3 time standarize,', t2-t1)
        # print()
        # t1 = time.time()
        self.summation.sort()
        # t2 = time.time()
        # print()
        # print('half time standarize,', t2-t1)
        # print()
        # t1 = time.time()

        self.standarize_triplet_R3
        self.standarize_triplet_R2p
        self.standarize_triplet_R2m
        self.standarize_rl()
        self.standarize_g()
        self.standarize_gs()
 #       print('przed gm2', self)
        self.standarize_gm2()
  #      print('po gm2', self)
        self.standarize_gm34()
    #    print('po gm34', self)
#        self.standarize_gm34spin()
#        self.standarize_gm234spin_new()
   #     print('po gm34 spin new', self)
        self.standarize_f()
        self.standarize_t()
        self.standarize_tf()
        self.standarize_s()
        self.standarize_v()
        self.standarize_delta()
        self.standarize_TT()
        # t2 = time.time()
        # print()
        #print('WHOLE time standarize,', t2-t0)
     #   print('po standarize', self)
      #  print()
        # print()



    def standarize_names(self):
        self.standarize_triplet_R3()
        self.standarize_triplet_R2p()
        self.standarize_triplet_R2m()
        self.standarize_rl()
        self.standarize_g()
        self.standarize_f()
        self.standarize_t()
        self.standarize_tf()
        self.standarize_s()
        self.standarize_v()
        self.standarize_delta()

    def best_rearrange(self):
        #print('standarro', self)
        check = False
        sflist = [self]
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == F12_AMPLITUDE:
                for elem in sflist:
                    check = True
                    sf1 = deepcopy(elem)
                    sf2 = deepcopy(elem)
                    sf3 = deepcopy(elem)
                    sf4 = deepcopy(elem)

                    pair1 = [self.coefficient_idx[i][0], self.coefficient_idx[i][2]]
                    pair2 = [self.coefficient_idx[i][1], self.coefficient_idx[i][3]]
                    # print(pair1, pair2)
                    if pair1[0] in self.summation and pair1[1] in self.summation:
                        # print('zamieniam pare', pair1)
                        sf2.multisubst([pair1[1]], ['z'])
                        sf2.multisubst([pair1[0]], [pair1[1]])
                        sf2.multisubst(['z'], [pair1[0]])
                    if pair2[0] in self.summation and pair2[1] in self.summation:
                        # print('zamieniam pare', pair2)
                        sf3.multisubst([pair2[1]], ['z'])
                        sf3.multisubst([pair2[0]], [pair2[1]])
                        sf3.multisubst(['z'], [pair2[0]])
                        
                    if pair1[0] in self.summation and pair1[1] in self.summation and\
                            pair2[0] in self.summation and pair2[1] in self.summation:
                        # print('zamieniam pare', pair1, pair2)
                        sf4.multisubst([pair1[1]], ['z'])
                        sf4.multisubst([pair1[0]], [pair1[1]])
                        sf4.multisubst(['z'], [pair1[0]])
                        sf4.multisubst([pair2[1]], ['z'])
                        sf4.multisubst([pair2[0]], [pair2[1]])
                        sf4.multisubst(['z'], [pair2[0]])
                    
                    sf2.standarize_names()
                    sf3.standarize_names()
                    sf4.standarize_names()

                    # #print('----------------------------------------------------------------')
                    # print(i, len(sflist), len(self.coefficient), self.coefficient_idx[i], pair1, pair2)
                    # for elem2 in sflist:
                    #     print(elem2)
                    # print(self)
                    # print(sf2)
                    # print(sf3)
                    # print(sf4)
                    #print('--------------------------------------------------------------------')
                sflist.append(sf2)
                sflist.append(sf3)
                sflist.append(sf4)
                # print('----------------------------------------------------------------')
                # print(i, len(sflist), len(self.coefficient), self.coefficient_idx[i], pair1, pair2)
                # for elem2 in sflist:
                #     print(elem2)
                # print('--------------------------------------------------------------------')

        #print('przed check', self)
        if check:
            # print('')
            # print('------1', self)
            # print('------2', sf2)
            # print('------3', sf3)
            # print('------4', sf4)
            # print('')

            sflist = [self, sf2, sf3, sf4]
            pointlist = []
            for x in range(0, len(sflist)):
                pointlist.append(100)
            # pointlist = [100, 100, 100, 100]
            for j in range(0, len(sflist)):
                for i in range(0, len(sflist[j].coefficient)):
                    if self.coefficient[i] != F12_AMPLITUDE:
                        mini = deepcopy(sflist[j].coefficient_idx[i])
                        minisort = deepcopy(mini)
                        minisort.sort()
                        # print('minisort', mini, minisort)
                        for k in range(0, len(minisort)):
                            # print('vegeta', k, mini.index(minisort[k]))
                            pointlist[j] = pointlist[j] - (abs(k - mini.index(minisort[k])))

            # print('ppppp', pointlist)
            p = pointlist.index(max(pointlist))
            # print('OTO PP', sflist[p])
            # print('po check1', self)
            self.coefficient = sflist[p].coefficient
            self.coefficient_idx = sflist[p].coefficient_idx
            self.num_factor = sflist[p].num_factor
            #print('po check2', self)

        
    def standarize_drag_list(self, lst, xfx):
        coef_dict = dict()
        k = 1
        for i in range(0, len(self.coefficient)):
            if i not in list(coef_dict.keys()):
                coef_dict[self.coefficient[i]] = 10**k
                k += 2

                
        for x in lst:
            self.summation.append(x)
        self.refine_summation_drag(lst, xfx)

        self.rename_dummy(coef_dict, "add")

        self.sort_indices()

        self.sort_coeffs()

        self.rename_dummy_back('usual')
#        self.rename_fixed("remove")

        self.reallocate_dummy()

        self.summation.sort()

        self.standarize_triplet_R3

        self.standarize_triplet_R2p
        self.standarize_triplet_R2m
        self.standarize_rl()
        #print('au4')
        self.standarize_g()
        self.standarize_f()
        self.standarize_t()
        self.standarize_t()

        self.standarize_s()
        self.standarize_v()
        self.standarize_delta()
#        #print('e', self)
        



    def standarize_interm(self):
        """Arrange indices in self.summation, self.coefficient, self.coefficient_idx,
        self.delta into standard form.
        """
        #
        # Create dictionary of differend values for different
        # coefficients
        #
        coef_dict = dict()
        k = 1
        for i in range(0, len(self.coefficient)):
            if i not in list(coef_dict.keys()):
                coef_dict[self.coefficient[i]] = 10**k
                k += 2
#        print('self', self)
        self.refine_fixed_interm_0()
#        print('po refine fixed0        ', self)
        self.refine_summation()
#        print('po refine sum         ', self)
        self.refine_fixed_interm()
#        print('po refine fixed         ', self)
        self.rename_fixed("add")
#        print('po rename fixed', self)
        self.rename_dummy(coef_dict, "add")
#        print('po rename dum  ', self)
        self.sort_indices()
#        print('po sord ix     ', self)
        self.sort_coeffs()
#        print('po sort coeff  ', self)
        self.rename_dummy_back('usual')
#        print('po r back      ', self)
        self.rename_fixed("remove")
#        print('po r fixed     ', self)
        self.reallocate_dummy()
#        print('po reallo      ', self)
        self.summation.sort()
        self.standarize_triplet_R3
        self.standarize_triplet_R2p
        self.standarize_triplet_R2m
        self.standarize_rl()
#        print('au5')
        self.standarize_g()
        self.standarize_f()
#        print('stand g        ', self)
        self.standarize_t()
        self.standarize_tf()

#        print('stand t        ', self)
        self.standarize_s()
#        print('stand s        ', self)
        self.standarize_v()
        self.standarize_delta()
#        print('stand v        ', self)


    def standarize_t_for_fortran(self):

        for j in range(0, len(self.coefficient)):
            if self.coefficient[j] == CC_AMPLITUDE and len(self.coefficient_idx[j]) > 2:
                tempv = []
                tempo = []
                for k in range(0, len(self.coefficient_idx[j]), 2):
                    tempv.append(self.coefficient_idx[j][k])
                    tempo.append(self.coefficient_idx[j][k+1])
                temp = tempv + tempo
                self.coefficient_idx[j] = temp

    def standarize_s_for_fortran(self):
        for j in range(0, len(self.coefficient)):
            if self.coefficient[j] == S_AMPLITUDE and len(self.coefficient_idx[j]) > 2:
                tempv = []
                tempo = []
                for k in range(0, len(self.coefficient_idx[j]), 2):
                    tempv.append(self.coefficient_idx[j][k])
                    tempo.append(self.coefficient_idx[j][k+1])
                temp = tempv + tempo
                self.coefficient_idx[j] = temp


    def sum_priority(self):
        """ Sorts summation indices so they are optimal
        for fortran code.
        """
        slow = dict()
        for i in self.summation:
            val = 0
            if i in occupied:
                ov = 'o'
            if i in virtual:
                ov = 'v'
            if i in complete:
                ov = 'c'
            if i in completev:
                ov = 'cv'
            if i in CABS or i in fortran_CABS:
                ov = 'cbs'        
            for j in range(0, len(self.coefficient)):
                #
                # Creating temporary list with Fortran amplitude indices convention
                # i.e. taking into account that Fortran arrays take the following
                # index ordering: t2(a, b, i, j)
                #
                if self.coefficient[j] == CC_AMPLITUDE and len(self.coefficient_idx[j]) > 2:
                    tempv = []
                    tempo = []
                    for k in range(0, len(self.coefficient_idx[j]), 2):
                        tempv.append(self.coefficient_idx[j][k])
                        tempo.append(self.coefficient_idx[j][k+1])
                    temp = tempv + tempo
                    for k in range(0, len(temp)):
                        if i == temp[k]:
                            if k > val:
                                val = k

                elif self.coefficient[j] == S_AMPLITUDE and len(self.coefficient_idx[j]) > 2:
                    tempv = []
                    tempo = []
                    for k in range(0, len(self.coefficient_idx[j]), 2):
                        tempv.append(self.coefficient_idx[j][k])
                        tempo.append(self.coefficient_idx[j][k+1])
                    temp = tempv + tempo
                    for k in range(0, len(temp)):
                        if i == temp[k]:
                            if k > val:
                                val = k
    
                else:
                    for k in range(0, len(self.coefficient_idx[j])):
                        if i == self.coefficient_idx[j][k]:
                            if k > val:
                                val = k
            slow[i] = str(val) + ov + i

        self.summation = sorted(slow, key = slow.__getitem__, reverse = True)


    def rename_dummy_opt(self):
#        print('self', self)
        for i in range(0, len(self.summation)):
            old =  str(self.summation[i])
            new = "0" + str(self.summation[i])
            self.substitute(old, new)
            
        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == TWOEL_INT:
                for j in range(0, len(self.coefficient_idx[i])):
                    if self.coefficient_idx[i][j] in virtualall:
                        self.coefficient_idx[i][j] = "0" +  str(self.coefficient_idx[i][j])
                    elif (len(self.coefficient_idx[i][j])) > 1:
                        if  self.coefficient_idx[i][j][1:] in virtualall:
                            self.coefficient_idx[i][j] = "0" +  str(self.coefficient_idx[i][j])


    def optimize(self):

        """ Optimize self for fortran coding
        """
        self.rename_dummy_opt()
        """ Standarize_g_for fortran is made in order to obtain
        two electron integrals in FORTRAN convenction i.e.
        all virtual orbitals are shifted to the left as much
        as it is possible in accordance with integral symmetry.
        Only available integrals are: vvvv, vvvo, vvoo, vooo, vovo, oooo
        """
        # print('au6', self)
        self.standarize_g_for_fortran()
        self.standarize_v_for_fortran()
        self.standarize_t_for_fortran()
        self.standarize_s_for_fortran()
        # print(self)
        self.rename_dummy_back('opt')
        # print(self)
        self.sum_priority()
        # print('au66', self)

    def optimize_f12(self):

        """ Optimize self for fortran coding
        """
#        self.rename_dummy_opt()
        """ Standarize_g_for fortran is made in order to obtain
        two electron integrals in FORTRAN convenction i.e.
        all virtual orbitals are shifted to the left as much
        as it is possible in accordance with integral symmetry.
        Only available integrals are: vvvv, vvvo, vvoo, vooo, vovo, oooo
        """
        # print('au6', self)
        self.standarize_v_for_fortran()
        self.standarize_t_for_fortran()
        # print(self)
 #       self.rename_dummy_back('opt')
        # print(self)
  #      self.sum_priority()
        # print('au66', self)

    def arrange_number(self, i):
        ref_oc = []
        ref_vt = []
        ref_gn = []
        for j in range(i+1, len(self.coefficient)):
            for k in self.coefficient_idx[j]:
                if k in self.summation:
                    if k in occupied:
                        ref_oc.append(k)
                    if k in virtual:
                        ref_vt.append(k)
                    if k in general:
                        ref_gn.append(k)
                    
        n = bubblesort(ref_oc)
        r = bubblesort(ref_vt)
        s = bubblesort(ref_gn)
        return n+r+s

    def number_of_g_idx_in_sum(self, i):
        result = 0
        for k in self.coefficient_idx[i]:
            if k in self.summation:
                result += 1
        return result

    def standarize_all(self, tpl):
        """ first part: standarize, second part: 
        standarization of operators according to given 
        alphaber (tpl). Gives back self, and result in form 
        of arithmetic string with terms arosen from 
        anticomutation relation beetwen operators
        """
        self.standarize()
        
        lista = list(tpl)
        result = arithmetic_string()
        counter = 1
        while counter != 0:
            counter = 0
            for i in range(len(lista)-1, 0, -1):
                if lista[i] < i+1 and lista[i-1] >= (i):
                    idx_tmp = deepcopy(lista[i])
                    lista[i] = deepcopy(lista[i-1])
                    lista[i-1] = deepcopy(idx_tmp)
                    res, selfnew = self.swap(i)
                    res.exec_delta()
                    res.as_standarize()
                    result = result + res
                    self = selfnew
                    counter += 1
        return result


    # def swap(self, i):
    #     """Interchange operator i and (i-1):
    #     ...E_{pq} E_{rs}... = ...E_{rs} E_{pq}... + \delta_{rq}... E_{ps}...
    #     - \delta_{ps} ...E_{qr}...,
    #     where dots denote other unitary group generators. Right hand side of
    #     the above identity is returned as a result of this function. Sum of
    #     two terms containing Kronecker deltas is returned as arithmetic_string,
    #     and the remaining term is returned as ugg instance.
    #     """
    #     res0 = deepcopy(self)
    #     res1 = deepcopy(self)
    #     res2 = deepcopy(self)

    #     idx_list = [self.operator_idx[i-1], self.operator_idx[i]]
    #     spin_list = [self.operator_type[i-1], self.operator_type[i]]
        
    #     if spin_list[0] == "s" and spin_list[1] == "s":
    #         spin = 0
    #     elif spin_list[0] == "t" and spin_list[1] == "s":
    #         spin = 1
    #     elif spin_list[0] == "s" and spin_list[1] == "t":
    #         spin = 1
    #     elif spin_list[0] == "t" and spin_list[1] == "t":
    #         spin = 0
    #     else:
    #         print("Spin operator not yet implemented")
    #         sys.exit(1)

        
        # #
        # # res0 <- E_{rs}E_{pq}
        # #
        # temp = deepcopy(res0.operator_idx[i-1])
        # res0.operator_idx[i-1] = deepcopy(res0.operator_idx[i])
        # res0.operator_idx[i] = deepcopy(temp)

        # temp = deepcopy(res0.operator_type[i-1])
        # res0.operator_type[i-1] = deepcopy(res0.operator_type[i])
        # res0.operator_type[i] = deepcopy(temp)

        # del res1.operator_idx[i]
        # del res2.operator_idx[i]

        # del res1.operator_type[i]
        # del res2.operator_type[i]

        # contracted = arithmetic_string()
        # #
        # # \delta_{rq} E_{ps}
        # #
        # res1.operator_idx[i-1] = [idx_list[0][0], idx_list[1][1]]
        # res1.new_delta(idx_list[1][0], idx_list[0][1])
        # res1_delta = [idx_list[1][0], idx_list[0][1]]
        
        # if(spin) == 1:
        #     res1.operator_type[i-1] = "t"
        # if(spin) == 0:
        #     res1.operator_type[i-1] = "s"

        # if not is_delta_zero(res1_delta):
        #     contracted.append(res1)
        # #
        # # -\delta{ps} E_{qr}
        # #
        # res2.operator_idx[i-1] = [idx_list[1][0], idx_list[0][1]]
        # res2.new_delta(idx_list[0][0], idx_list[1][1])
        # res2_delta = [idx_list[0][0], idx_list[1][1]]
        # res2.num_factor *= -1
        # if(spin) == 1:
        #     res2.operator_type[i-1] = "t"
        # if(spin) == 0:
        #     res2.operator_type[i-1] = "s"

        # if not is_delta_zero(res2_delta):
        #     contracted.append(res2)

        # return contracted, res0


    def swap(self, i):
        """Interchange operator i and (i-1):
        ...E_{pq} E_{rs}... = ...E_{rs} E_{pq}... + \delta_{rq}... E_{ps}...
        - \delta_{ps} ...E_{qr}...,
        where dots denote other unitary group generators. Right hand side of
        the above identity is returned as a result of this function. Sum of
        two terms containing Kronecker deltas is returned as arithmetic_string,
        and the remaining term is returned as ugg instance.

        E_mn E_pq =     E_pq E_mn    + E_mq   \delta_pn   - E_pn \delta_mq

        T0_mn T0_pq =   T0_pq T0_mn  + E_mq   \delta_pn   - E_pn \delta_mq

        T0_mn E_pq =    E_pq T0_mn   + T0_mq   \delta_pn  - T0_pn \delta_mq
        E_mn T0_pq =    T0_pq E_mn   + T0_mq   \delta_pn  - T0_pn \delta_mq
        
        E_mn Tm1_pq =   Tm1_pq E_mn  + Tm1_mq  \delta_pn  - Tm1_pn \delta_mq
        Tm1_mn E_pq =   E_pq Tm1_mn  + Tm1_mq  \delta_pn  - Tm1_pn \delta_mq

        E_mn T1_pq =    T1_pq E_mn   + T1_mq   \delta_pn  - T1_pn \delta_mq
        T1_mn E_pq =    E_pq T1_mn   + T1_mq   \delta_pn  - T1_pn \delta_mq

        T0_mn T1_pq =   T1_pq T0_mn  + T1_mq   \delta_pn  + T1_pn \delta_mq
        T1_mn T0_pq =   T0_pq T1_mn  - T1_mq   \delta_pn  - T1_pn \delta_mq
        
        T0_mn Tm1_pq =  Tm1_pq T0_mn - Tm1_mq  \delta_pn  - Tm1_pn \delta_mq
        Tm1_mn T0_pq =  T0_pq Tm1_mn + Tm1_mq  \delta_pn  + Tm1_pn \delta_mq

        T1_mn Tm1_pq =  Tm1_pq T1_mn + 0.5 E_pn \delta_mq - 0.5 T0_pn \delta_mq
        Tm1_mn T1_pq =  T1_pq Tm1_mn - 0.5 E_mq \delta_pn + 0.5 T0_mq \delta_pn
        
        T1_mn T1_pq =   T1_pq T1_mn
        
        Tm1_mn Tm1_pq = Tm1_pq Tm1_pq
        
        """
        res0 = deepcopy(self)
        res1 = deepcopy(self)
        res2 = deepcopy(self)
        # print('rex', res0)
        idx_list = [self.operator_idx[i-1], self.operator_idx[i]]
        spin_list = [self.operator_type[i-1], self.operator_type[i]]

        #
        # res0 <- E_{rs}E_{pq}
        #
        temp = deepcopy(res0.operator_idx[i-1])
        res0.operator_idx[i-1] = deepcopy(res0.operator_idx[i])
        # print('rex1', res0)
        res0.operator_idx[i] = deepcopy(temp)
        # print('rex2', res0)

        temp = deepcopy(res0.operator_type[i-1])
        res0.operator_type[i-1] = deepcopy(res0.operator_type[i])
        # print('rex3', res0)
        res0.operator_type[i] = deepcopy(temp)
        # print('rex4', res0)
        # print(temp, res0)
        del res1.operator_idx[i]
        del res2.operator_idx[i]

        del res1.operator_type[i]
        del res2.operator_type[i]
#        print('sra')
        contracted = arithmetic_string()
        #
        # \delta_{rq} E_{ps}
        #
        if spin_list[0] == "t1" and spin_list[1] == "tm1":
            res1.operator_idx[i-1] = [idx_list[1][0], idx_list[0][1]]
            res1.new_delta(idx_list[0][0], idx_list[1][1])
            res1_delta = [idx_list[0][0], idx_list[1][1]]
        else:
            res1.operator_idx[i-1] = [idx_list[0][0], idx_list[1][1]]
            res1.new_delta(idx_list[1][0], idx_list[0][1])
            res1_delta = [idx_list[1][0], idx_list[0][1]]
        

        if spin_list[0] == "s" and spin_list[1] == "s":
            res1.operator_type[i-1] = "s"
        elif spin_list[0] == "t0" and spin_list[1] == "t0":
            res1.operator_type[i-1] = "s"
        elif spin_list[0] == "t0" and spin_list[1] == "s":
            res1.operator_type[i-1] = "t0"
        elif spin_list[0] == "s" and spin_list[1] == "t0":
            res1.operator_type[i-1] = "t0"
        elif spin_list[0] == "s" and spin_list[1] == "tm1":
            res1.operator_type[i-1] = "tm1"
        elif spin_list[0] == "tm1" and spin_list[1] == "s":
            res1.operator_type[i-1] = "tm1"
        elif spin_list[0] == "s" and spin_list[1] == "t1":
            res1.operator_type[i-1] = "t1"
        elif spin_list[0] == "t1" and spin_list[1] == "s":
            res1.operator_type[i-1] = "t1"
        elif spin_list[0] == "t0" and spin_list[1] == "t1":
            res1.operator_type[i-1] = "t1"
        elif spin_list[0] == "t1" and spin_list[1] == "t0":
            res1.operator_type[i-1] = "t1"
            res1.num_factor *= -1
        elif spin_list[0] == "t0" and spin_list[1] == "tm1":
            res1.operator_type[i-1] = "tm1"
            res1.num_factor *= -1
        elif spin_list[0] == "tm1" and spin_list[1] == "t0":
            res1.operator_type[i-1] = "tm1"
        elif spin_list[0] == "t1" and spin_list[1] == "tm1":
            res1.operator_type[i-1] = "s"
            res1.num_factor *= 0.5
        elif spin_list[0] == "tm1" and spin_list[1] == "t1":
            res1.operator_type[i-1] = "s"
            res1.num_factor *= -0.5
        # print('res1', res1)
        if not is_delta_zero(res1_delta):
            if not (spin_list[0] == "t1" and spin_list[1] == "t1") \
               and not (spin_list[0] == "tm1" and spin_list[1] == "tm1"):
                contracted.append(res1)
        #
        # -\delta{ps} E_{qr}
        #

        if spin_list[0] == "tm1" and spin_list[1] == "t1":
            res2.operator_idx[i-1] = [idx_list[1][1], idx_list[0][0]]
            res2.new_delta(idx_list[0][1], idx_list[1][0])
            res2_delta = [idx_list[0][1], idx_list[1][0]]
        else:
            res2.operator_idx[i-1] = [idx_list[1][0], idx_list[0][1]]
            res2.new_delta(idx_list[0][0], idx_list[1][1])
            res2_delta = [idx_list[0][0], idx_list[1][1]]
            


        if spin_list[0] == "s" and spin_list[1] == "s":
            res2.operator_type[i-1] = "s"
            res2.num_factor *= -1
        elif spin_list[0] == "t0" and spin_list[1] == "t0":
            res2.operator_type[i-1] = "s"
            res2.num_factor *= -1
        elif spin_list[0] == "t0" and spin_list[1] == "s":
            res2.operator_type[i-1] = "t0"
            res2.num_factor *= -1
        elif spin_list[0] == "s" and spin_list[1] == "t0":
            res2.operator_type[i-1] = "t0"
            res2.num_factor *= -1
        elif spin_list[0] == "s" and spin_list[1] == "tm1":
            res2.operator_type[i-1] = "tm1"
            res2.num_factor *= -1
        elif spin_list[0] == "tm1" and spin_list[1] == "s":
            res2.operator_type[i-1] = "tm1"
            res2.num_factor *= -1
        elif spin_list[0] == "s" and spin_list[1] == "t1":
            res2.operator_type[i-1] = "t1"
            res2.num_factor *= -1
        elif spin_list[0] == "t1" and spin_list[1] == "s":
            res2.operator_type[i-1] = "t1"
            res2.num_factor *= -1
        elif spin_list[0] == "t0" and spin_list[1] == "t1":
            res2.operator_type[i-1] = "t1"
        elif spin_list[0] == "t1" and spin_list[1] == "t0":
            res2.operator_type[i-1] = "t1"
            res2.num_factor *= -1
        elif spin_list[0] == "t0" and spin_list[1] == "tm1":
            res2.operator_type[i-1] = "tm1"
            res2.num_factor *= -1
        elif spin_list[0] == "tm1" and spin_list[1] == "t0":
            res2.operator_type[i-1] = "tm1"
        elif spin_list[0] == "t1" and spin_list[1] == "tm1":
            res2.operator_type[i-1] = "t0"
            res2.num_factor *= -0.5
        elif spin_list[0] == "tm1" and spin_list[1] == "t1":
            res2.operator_type[i-1] = "t0"
            res2.num_factor *= 0.5


        # print('res2', res2)
        if not is_delta_zero(res2_delta):
            if not (spin_list[0] == "t1" and spin_list[1] == "t1") \
               and not (spin_list[0] == "tm1" and spin_list[1] == "tm1"):
                contracted.append(res2)
        # print('rex5', res0)
        return contracted, res0


    def __add__(self, other):
        a = arithmetic_string(self, other)
        return a

    # def weigh_idx(self, i):        
    #     wi = 0
    #     for k in range(0, len(self.coefficient)):
    #         coef_w = AMP_NAME_DICT[self.coefficient[k]]
    #         for idx in self.coefficient_idx[k]:
    #             if i == idx:
    #                 wi += coef_w
    #     return wi
    

    def refine_summation(self, last_idx = '', cas=False, drag_fixed=[]):
        """ \sum_{rp}E_{rp}E_{ai} --> \sum_{pq}E_{pq}E_{ai}
        Changes summation indices so they are first letters
        in each alphabets (virtual, occupied, general).
        """
        #
        # Append fixed indices to the 
        # list of indices excluded from
        # substitutions
        #

        indices = self.indices()
        banned = indices - set(self.summation) 
        banned = banned | set(fixed)

        summation_old = self.summation
        summation_new = []

        sum_virt = []
        sum_occ = []
        sum_gen = []
        sum_comp = []
        sum_compv = []
        

        # t1 = time.time()
        refine_idx = []
        # self_sort = deepcopy(self)
        # self_sort.sort_coeffs(last_idx, cas)
        refine_idx = self.sort_nosort_coeffs(last_idx, cas)
        # t2 = time.time()
        # print('deepcop', t2-t1)



        summation_old = refine_idx
        if drag_fixed !=[]:
            self.df = []
            self.df.append(drag_fixed)
            self.df.append(drag_fixed)
        for i in summation_old:
            if i in virtual:
                j = free_idx(virtual, banned)
            elif i in occupied:
                j = free_idx(occupied, banned)
            elif i in complete:
                j = free_idx(complete, banned)
            elif i in completev:
                j = free_idx(completev, banned)
            elif i in general:
                j = free_idx(general, banned)
            elif i in CABS:
                j = free_idx(CABS, banned)
            else:
                pool = idx_type(i)
                j = free_idx(pool, banned)
            summation_new.append(j)
            if drag_fixed != []:
                if i in drag_fixed:
                    self.df[1][drag_fixed.index(i)] = j
            banned = banned | set([j])

        s_old = deepcopy(summation_old)
        self.multisubst(summation_old, summation_new)
        print('summation_old', summation_old)
        print('summation_new', summation_new)


        temp_values = {}
        for i in range(len(summation_old)):
            if summation_old[i] in self.orbital_type:
                temp_values[summation_new[i]] = self.orbital_type[summation_old[i]]
                del self.orbital_type[summation_old[i]]
                
        self.orbital_type.update(temp_values)

        
        if self.name != "":
            name_old_lst = []
            for i in self.name:
                name_old_lst.append(i)
            self.refine_name(name_old_lst, s_old, summation_new)


    def refine_summation_drag(self, lst, xfx):
        """ \sum_{rp}E_{rp}E_{ai} --> \sum_{pq}E_{pq}E_{ai}                                                              
        Changes summation indices so they are first letters                                                                 
        in each alphabets (virtual, occupied, general).                                                   
        """
        indices = self.indices()
        #print('indices', indices)
        #print('ssum', self.summation)
        banned = indices - set(self.summation)
        #print('banned1', banned)
        #print('fixed', fixed)
        banned = banned | set(fixed)
        #print('banned2', banned)
        summation_old = self.summation
        summation_new = []

        #print('sefl przed', self)
        #print('xfx przed', xfx)
        #print('banned', banned)
        for i in summation_old:
            if i in virtual:
                j = free_idx(virtual, banned)
            elif i in occupied:
                j = free_idx(occupied, banned)
            elif i in complete:
                j = free_idx(complete, banned)
            elif i in completev:
                j = free_idx(completev, banned)
            elif i in general:
                j = free_idx(general, banned)
            elif i in CABS:
                j = free_idx(CABS, banned)
            else:
                pool = idx_type(i)
                # print(pool)
                j = free_idx(pool, banned)
            if i in lst:
                lstidx = lst.index(i)
                lst[lstidx] = j
            if i in xfx:
                xfxidx = xfx.index(i)
                xfx[xfxidx] = j



            summation_new.append(j)
            banned = banned | set([j])
        #print('soldnew', summation_old, summation_new)
        #print('xfx po ', xfx)
        s_old = deepcopy(summation_old)
        self.multisubst(summation_old, summation_new)
        #print('pomsub', self)
        if self.name != "":
            name_old_lst = []
            for i in self.name:
                name_old_lst.append(i)
            self.refine_name(name_old_lst, s_old, summation_new)

        #print(lst)
        sumold = deepcopy(self.summation)
        for x in sumold:
            # print(sumold, x)
            if x in lst:
                # print('jest')
                del self.summation[self.summation.index(x)]

                
    def refine_name(self, name, s_old, summation_new):
        
        for i in range(0, len(name)):
            for j in range(0, len(s_old)):
                if name[i] == s_old[j]:
                    name[i] = summation_new[j]

        self.name = ''.join(name)



    def refine_fixed_interm_0(self):

        len_virt = 0 
        len_occ = 0
        for x in self.summation:
            if x in virtual:
                len_virt += 1
            elif x in occupied:
                len_occ += 1
            else:
                print('General index in summation: aborting process')
                sys.exit(1)

        banned = set()
        for i in range(0, len_virt):
            banned = banned | set(virtual[i])
        for i in range(0, len_occ):
            banned = banned | set(occupied[i])

        banned = banned | set(self.summation)



        indices_old = []
        
        for i in range(0, len(self.coefficient_idx)):
            for j in self.coefficient_idx[i]:
                if j not in self.summation and j not in indices_old:
                    indices_old.append(j)
        
        indices_new = []

        for i in indices_old:
            if i in virtual:
                j = free_idx(virtual, banned)
            elif i in occupied:
                j = free_idx(occupied, banned)

            indices_new.append(j)
            banned = banned | set([j])

        #
        # Substitute all changed indices
        #
        self.multisubst(indices_old, indices_new)


    def refine_fixed_interm(self):

        banned = set(self.summation)

        indices_old = []
        
        for i in range(0, len(self.coefficient_idx)):
            for j in self.coefficient_idx[i]:
                if j not in banned and j not in indices_old:
                    indices_old.append(j)
        
        indices_new = []

        for i in indices_old:
            if i in virtual:
                j = free_idx(virtual, banned)
            elif i in occupied:
                j = free_idx(occupied, banned)

            indices_new.append(j)
            banned = banned | set([j])
        #
        # Substitute all changed indices
        #
        self.multisubst(indices_old, indices_new)


    def fromleft(self, e):
        """Concatenate from left: E >> SELF"""

        res = ugg()
        res.summation       = deepcopy(e.summation) + deepcopy(self.summation)
        res.coefficient     = deepcopy(e.coefficient) + deepcopy(self.coefficient)
        res.coefficient_idx = deepcopy(e.coefficient_idx) + deepcopy(self.coefficient_idx)
        res.operator_idx    = deepcopy(e.operator_idx) + deepcopy(self.operator_idx)
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
        res = ugg()
        res.summation       = deepcopy(self.summation) + deepcopy(e.summation)
        res.coefficient     = deepcopy(self.coefficient) + deepcopy(e.coefficient)
        res.coefficient_idx = deepcopy(self.coefficient_idx) + deepcopy(e.coefficient_idx)
        res.operator_idx    = deepcopy(self.operator_idx) + deepcopy(e.operator_idx)
        res.operator_type  = deepcopy(self.operator_type) + deepcopy(e.operator_type)
        res.num_factor      = self.num_factor * e.num_factor
        
        if self.delta != []:
            for delta in self.delta:
                res.new_delta(delta[0], delta[1])

        if e.delta != []:
            for delta in e.delta:
                res.new_delta(delta[0], delta[1])
        return res       


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
                if x in monab:
                    summ = summ + Tdict[x]+ str("")
                else:
                    summ = summ + str(x) + str("")                
            sumidx = "\sum_{" + summ + "}"
        upperi = ""
        loweri = ""
        if self.permutation_op != []:
            for y in self.permutation_op:
                upperi += str(y[0]) + str("")
                loweri += str(y[1]) + str("")
            permidx += "P_{"+loweri+"}^{"+upperi+"}"
        if len(self.operator_idx) != 0:
            for x in range(0, len(self.operator_idx)):
                temp = ""
                for y in self.operator_idx[x]:
                    temp = temp + str(y) + str("")
                if self.operator_type[x] == "t0":
                    opidx = opidx + "T_{"+temp+"}"
                elif self.operator_type[x] == "t1":
                    opidx = opidx + "T1_{"+temp+"}"
                elif self.operator_type[x] == "tm1":
                    opidx = opidx + "Tm1_{"+temp+"}"
                elif self.operator_type[x] == "s":
                    opidx = opidx + "E_{"+temp+"}"

        if len(self.coefficient) != 0:
            for x in range(0, len(self.coefficient)):
                temp = ""

                if self.coefficient[x] == CC_AMPLITUDE or self.coefficient[x] == S_AMPLITUDE\
                        or self.coefficient[x] == TEMP1 or self.coefficient[x] == TEMP2 \
                        or self.coefficient[x] == EOM_TRIPLET_R3 \
                        or self.coefficient[x] == EOM_TRIPLET_R2p \
                        or self.coefficient[x] == EOM_TRIPLET_R2m \
                        or self.coefficient[x] == EOM_TRIPLET_R1 \
                        or self.coefficient[x] == EOM_CC_SINGLE_Rr \
                        or self.coefficient[x] == EOM_CC_SINGLE_Rl \
                        or self.coefficient[x] == EOM_CC_SINGLE_Rr_plus \
                        or self.coefficient[x] == EOM_CC_SINGLE_Rl_plus \
                        or self.coefficient[x] == EOM_CC_SINGLE_Rr_minus \
                        or self.coefficient[x] == EOM_CC_SINGLE_Rl_minus \
                        or self.coefficient[x] == F12_AMPLITUDE\
                        or self.coefficient[x] == CI_AMPLITUDE\
                        or self.coefficient[x] == CC_TAU:
                        # or self.coefficient[x] == DENS3\
                        # or self.coefficient[x] == DENS4:
                    # or self.coefficient[x] == F12_TWOEL\
                                                
                    
                    
                    upperi = ""
                    loweri = ""
                    for y in range(0, len(self.coefficient_idx[x]), 2):
                        upperi += str(str(self.coefficient_idx[x][y]) + str(""))
                        loweri += str(str(self.coefficient_idx[x][y+1]) + str(""))
                    temp = loweri + "}^{" + upperi
                    #                    print('selfk', self.coefficient[x], self.coefficient_idx[x], upperi, loweri)

                elif self.coefficient[x] == DENS2PM:
                    upperi = ""
                    loweri = ""
                    for y in range(0, 2):
                        loweri += str(str(self.coefficient_idx[x][y]) + str(""))
                    for y in range(2, 4):
                        upperi += str(str(self.coefficient_idx[x][y]) + str(""))
                    #temp = loweri + "}^{" + upperi
                    temp = loweri + upperi
                elif self.coefficient[x] == DENS2P:
                    upperi = ""
                    loweri = ""
                    for y in range(0, 2):
                        loweri += str(str(self.coefficient_idx[x][y]) + str(""))
                    for y in range(2, 4):
                        upperi += str(str(self.coefficient_idx[x][y]) + str(""))
                    #temp = loweri + "}^{" + upperi
                    temp = loweri + upperi
                elif self.coefficient[x] == DENS3PM or self.coefficient[x] == DENS3P:
                    upperi = ""
                    loweri = ""
                    for y in range(0, 3):
                        loweri += str(str(self.coefficient_idx[x][y]) + str(""))
                    for y in range(3, 6):
                        upperi += str(str(self.coefficient_idx[x][y]) + str(""))
                    #temp = loweri + "}^{" + upperi
                    temp = loweri + upperi
                elif self.coefficient[x] == DENS4PM or self.coefficient[x] == DENS4PPM:
                    upperi = ""
                    loweri = ""
                    for y in range(0, 4):
                        loweri += str(str(self.coefficient_idx[x][y]) + str(""))
                    for y in range(4, 8):
                        upperi += str(str(self.coefficient_idx[x][y]) + str(""))
                    temp = loweri + "}^{" + upperi
                    #temp = loweri + upperi
                elif self.coefficient[x] == F12_TWOEL:
                    upperi = ""
                    loweri = ""
                    for y in self.coefficient_idx[x]:
                        if y in virtualall:
                            if y in completev:
                                yy = completev_latex_dict[y]
                            else:
                                yy = y
                            upperi += str(yy + str(""))
                        else:
                            if y in completev:
                                yy = completev_latex_dict[y]
                            else:
                                yy = y
                            loweri += str(yy + str(""))
                    temp = loweri + "}^{" + upperi
                elif self.coefficient[x] == TWOEL_INT_AS:
                    upperi = ""
                    loweri = ""
                    for y in range(0, 2):
                        upperi += str(self.coefficient_idx[x][y] + str(""))
                        loweri += str(self.coefficient_idx[x][y+2] + str(""))
                    temp = loweri + "}^{" + upperi
                elif self.coefficient[x] == INTERM_Ft_F12 or self.coefficient[x] == INTERM_Ftt_F12:
                    upperi = ""
                    loweri = ""
                    for y in range(0, 2):
                        upperi += str(self.coefficient_idx[x][y] + str(""))
                    for y in range(2, 4):
                        loweri += str(self.coefficient_idx[x][y] + str(""))

                    temp = loweri + "}^{" + upperi

                # elif self.coefficient[x] == F12_TWOEL:
                #     upperi = ""
                #     loweri = ""
                #     for y in range(0, 2):
                #         upperi += str(self.coefficient_idx[x][y] + str(""))
                #     for y in range(2, 4):
                #         loweri += str(self.coefficient_idx[x][y] + str(""))

                #     temp = loweri + "}^{" + upperi

                elif (self.coefficient[x] == SLAT2_SYM1 or self.coefficient[x] == SLAT2_SYM2 or self.coefficient[x] == SLAT3 or \
                          self.coefficient[x] == SLAT4):
                    upperi = ""
                    loweri = ""
                    for y in range(0, len(self.coefficient_idx[x]), 2):
                        upperi += str(self.coefficient_idx[x][y] + str(""))
                        loweri += str(self.coefficient_idx[x][y+1] + str(""))
                    temp = loweri + "}^{" + upperi

                else:
                    if self.coefficient[x][0:2] == 'TT':
                        if 'interm' in self.coefficient[x]:
                            for y in range(0, len(self.coefficient_idx[x])):
                                temp = temp + str(Tdict[self.coefficient_idx[x][y]]) + str("")
                        else:
                            upperi = ""
                            loweri = ""
                            off = int(len(self.coefficient_idx[x])/2)
                            for y in range(0, off):
                                upperi += str(Tdict[self.coefficient_idx[x][y]] + str(""))
                                loweri += str(Tdict[self.coefficient_idx[x][off+y]] + str(""))
                                temp = loweri + "}^{" + upperi                       
                    else:
                        for y in range(0, len(self.coefficient_idx[x])):
                            temp = temp + str(self.coefficient_idx[x][y]) + str("")


                if (self.coefficient[x] == EOM_CC_SINGLE_Rr or self.coefficient[x] == EOM_CC_SINGLE_Rl):
                    coef = 'R'
                elif (self.coefficient[x] == EOM_CC_SINGLE_Rr_plus or self.coefficient[x] == EOM_CC_SINGLE_Rl_plus):
                    coef = '{}^{(+)}\!R'
                elif (self.coefficient[x] == EOM_CC_SINGLE_Rr_minus or self.coefficient[x] == EOM_CC_SINGLE_Rl_minus):
                    coef = ' {}^{(-)}\!R'
                elif (self.coefficient[x] == SLAT2_SYM1 or self.coefficient[x] == SLAT2_SYM2 or self.coefficient[x] == SLAT3 or \
                          self.coefficient[x] == SLAT4):
                    coef = "o"
                elif (self.coefficient[x] == F12_AMPLITUDE):
                    coef = "t"
                elif (self.coefficient[x] == CC_TAU):
                    coef = "\tau"
                elif self.coefficient[x] == TWOEL_INT:
                    coef = ""
                elif self.coefficient[x] == TWOEL_INT_DIRAC or self.coefficient[x] == TWOEL_INT_DIRAC_A:
                    coef = ""
                elif self.coefficient[x] == DENS1:
                    coef = "\gamma"
                elif self.coefficient[x] == DENS2:
                    coef = "\Gamma"
                elif self.coefficient[x] == DENS3:
                    coef = "\Gamma"
                elif self.coefficient[x] == DENS4:
                    coef = "\Gamma"
                elif self.coefficient[x] == DENS3PM:
#                    coef = "\Gamma^{++-++-}"
                    coef = "\Gamma^{\\sixplusmin}"
#                    coef = "\Gamma^{\\alpha\\alpha\\beta\\alpha\\alpha\\beta}"
                    #coef = "\Gamma"
                elif self.coefficient[x] == DENS3P:
                    #coef = "\Gamma^{++++++}"
                    coef = "\Gamma^{\\sixplus}"
#                    coef = "\Gamma^{\\alpha\\alpha\\alpha\\alpha\\alpha\\alpha}"
                    #coef = "\\tilde{\Gamma}"
                elif self.coefficient[x] == DENS4PPM:
                    #coef = "\Gamma^{+++-+++-}"
                    #coef = "\Gamma^{\\alpha\\alpha\\alpha\\beta\\alpha\\alpha\\alpha\\beta}"
                    #coef = "\widetilde{\Gamma}"
                    coef = "\Gamma"
                elif self.coefficient[x] == DENS4PM:
                    #coef = "\Gamma^{++--++--}"
                    #coef = "\Gamma^{\\alpha\\alpha\\beta\\beta\\alpha\\alpha\\beta\\beta}"
                    coef = "\\tilde{\Gamma}"
                    #coef = "\Gamma"
                elif self.coefficient[x] == DENS2P:
                    #coef = "\Gamma^{++++}"
                    coef = "\Gamma^{\\fourplus}"
                    #coef = "\Gamma^{\\alpha\\alpha\\alpha\\alpha}"
                    #coef = "\\tilde{\Gamma}"
                elif self.coefficient[x] == DENS2M:
                    #coef = "\Gamma^{\\beta\\beta\\beta\\beta}"
                    #coef = "\Gamma^{----}"
                    coef = "\Gamma^{\\fourmin}"
                elif self.coefficient[x] == DENS2PM:
                    coef = "\Gamma^{\\plusmin}"
                    #coef = "\Gamma^{+-+-}"
                    #coef = "\Gamma"
                    #coef = "\Gamma^{\\alpha\\beta\\alpha\\beta}"
                elif self.coefficient[x] == DENS2MP:
                    coef = "\Gamma^{-+-+}"
                    coef = "\Gamma^{\\minplus}"
                    #coef = "\Gamma^{\\beta\\alpha\\beta\\alpha}"
                elif 'l' in self.coefficient[x]:
                    lenl = 0
                    for s in self.coefficient[x]:
                        if s == "l":
                            lenl += 1
                    no = self.coefficient[x][lenl:len(self.coefficient[x])]
                    if lenl == 1:
                        coef = f"(\Theta_{{{no}}})" 
                    elif lenl == 2:
                        coef = f"(\Xi_{{{no}}})" 
                    elif lenl == 3:
                        coef = f"(\\vartheta_{{{no}}})"
                    elif lenl == 4:
                        coef = f"(\\upsilon_{{{no}}})" 
                    else:
                        # print(self, lenl)
                        print('implement name for higher level intermediates')
                        sys.exit(0)
                elif 'interm' in self.coefficient[x]:
                    pattern = r".*interm(\d+)$"
                    num =  re.match(pattern, self.coefficient[x]).group(1)
                    #                    leni = len('interm')
                    #                    coef = "I^{"+ self.coefficient[x][leni:len(self.coefficient[x])]+"}"
                    if 'TTs' in self.coefficient[x]:
                        coef = "K^{" + num + "}"
                    elif 'TTz' in self.coefficient[x]:
                        coef = "Y^{" + num + "}"
                    else:
                        coef = "I^{" + num + "}"
                elif (self.coefficient[x] == F12_TWOEL):
                    coef = 'F'
                elif (self.coefficient[x] == INTERM_Ft_F12):
                    coef = '\\tilde{t}'
                elif (self.coefficient[x] == INTERM_Ftt_F12):
                    coef = '\\tilde{\\bar{t}}'
                else:
                    if self.coefficient[x][0:2] == 'TT':
                        coef = deepcopy(self.coefficient[x])
                        coef = coef.replace('TT', '')
                    else:
                        coef = self.coefficient[x]

                if self.coefficient[x] == TWOEL_INT:
                    cidx = cidx + str(coef)+ "("+temp[0:2]+"|"+temp[2:4]+")"
                elif self.coefficient[x] == TWOEL_INT_DIRAC:
                    cidx = cidx + str(coef)+ "<"+temp[0:2]+"|"+temp[2:4]+">"
                elif self.coefficient[x] == TWOEL_INT_DIRAC_A:
                    cidx = cidx + str(coef)+ "<"+temp[0:2]+"||"+temp[2:4]+">"      
                else:
                    cidx = cidx + str(coef)+ "_{"+temp+"}"
            
        for x in range(0, len(self.delta)):
            temp = ""
            for y in self.delta[x]:
                temp = temp + str(y) + str("")
            delidx = delidx + "\delta_{" + temp + "}"

        for x in range(0, len(self.long_delta)):
            temp = "\delta_"
            for y in self.long_delta[x]:
                temp = temp + str(y) + str("")
            delidx = delidx +  temp + ","

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


    def exec_delta_fixed(self, delta_subst):
        """ nonzero - list of nonzero delta. 
        All indices in deltas are fixed. When
        in self.delta list is at least one delta which is
        not on nonzero list, term is equal to zero.
        Deltas from nonzero list are executed.      
        """
        for p in self.delta:
            p.sort()
            if p[0] in delta_subst.keys() and p[1] in delta_subst.keys():
                if delta_subst[p[0]] != delta_subst[p[1]]:
                    self.num_factor = 0.0
            else:
                self.num_factor = 0.0
                break

        if self.num_factor != 0.0:
            for i in delta_subst.keys():
                self.substitute(i, delta_subst[i])
            self.delta = []

        for i in range(0, len(self.coefficient)):
            if self.coefficient[i] == FOCK_MATRIX:
                if self.coefficient_idx[i][0] != self.coefficient_idx[i][1]:
                    self.num_factor = 0.0


    def exec_delta_fixed_mem(self):
        # print(self)
        while len(self.delta) != 0:
            p = self.delta[0]
            # print(p, self.delta, len(self.delta))
            if p[0] in self.loops and p[1] in self.loops:
                self.loops.remove(p[1])
                if len(self.ibra) > 0:
                    indices = [i for i, x in enumerate(self.ibra) if x == p[1]]
                    for z in indices:
                        self.ibra[z] = p[0]
                else:
                    indices = [i for i, x in enumerate(self.iket) if x == p[1]]
                    for z in indices:
                        self.iket[z] = p[0]
                self.substitute(p[1], p[0])
            elif p[0] in self.loops and p[1] not in self.loops:
                self.loops.remove(p[0])
                if len(self.ibra) > 0:
                    indices = [i for i, x in enumerate(self.ibra) if x == p[0]]
                    for z in indices:
                        self.ibra[z] = p[1]
                else:
                    indices = [i for i, x in enumerate(self.iket) if x == p[0]]
                    for z in indices:
                        self.iket[z] = p[1]
                self.substitute(p[0], p[1])
            elif p[0] not in self.loops and p[1] in self.loops:
                self.loops.remove(p[1])
                if len(self.ibra) > 0:
                    indices = [i for i, x in enumerate(self.ibra) if x == p[1]]
                    for z in indices:
                        self.ibra[z] = p[0]
                else:
                    indices = [i for i, x in enumerate(self.iket) if x == p[1]]
                    for z in indices:
                        self.iket[z] = p[0]
                self.substitute(p[1], p[0])
            elif p[0] not in self.loops and p[1] not in self.loops:
                if p[0] == p[1]:
                    del self.delta[0]
                else:
                    self.num_factor = 0
                    #            print('jest delta miedzy dwoma indeksami ket, usuwam wyraz')
            else:
                # print('aha', self, p, self.loops, p[0], p[1])
                sys.exit(0)

        

                
    # def exec_fock(self):

    #     for p in range(0, len(self.summation)):
    #         for q in range(0, len(self.coefficient)):
    #             if self.coefficient[q] == FOCK_MATRIX:
    #                 if self.coefficient_idx[q][0] == self.summation[p]:
    #                     self.substitute(self.summation[p], self.coefficient_idx[q][1])
    #                 if self.coefficient_idx[q][1] == self.summation[p]:
    #                     self.substitute(self.summation[p], self.coefficient_idx[q][0])
    #     self.summation = []

    def exec_delta_cas_fixed(self):

        for i, l  in enumerate(self.delta):
            # print('il', i, l)
            
            p = self.delta[i][0]
            q = self.delta[i][1]

            p_type = self.orbital_type[p]
            q_type = self.orbital_type[q]
            # print('pqpq', p, q)


            if p_type != q_type:
                self.num_factor = 0.0
                return

            if p  not in self.summation and q not in self.summation:
                if p_type == q_type:
                    self.coefficient_idx = [[(p if x == q else x) for x in sublist] for sublist in self.coefficient_idx]

                else:
                    self.num_factor = 0.0
                    return
        
    def exec_delta_cas(self):

        # print('delta cas for this')
        # print(self)
        dellist = []
        for i, l  in enumerate(self.delta):
            # print('il', i, l)
            
            p = self.delta[i][0]
            q = self.delta[i][1]

            p_type = self.orbital_type[p]
            q_type = self.orbital_type[q]
            # print('pqpq', p, q)


            if p_type != q_type:
                self.num_factor = 0.0
                return

            if p in self.summation and q in self.summation:
#                print('a1', self.orbital_type)
                if p_type == q_type:
                    self.summation.remove(q)
                    self.coefficient_idx = [[(p if x == q else x) for x in sublist] for sublist in self.coefficient_idx]
                    dellist.append(i)
                else:
                    self.num_factor = 0.0
                    return

            elif p in self.summation and q not in self.summation:
#                print('a2', self.orbital_type)
                if p_type == q_type:
                    self.summation.remove(p)
                    self.coefficient_idx = [[(q if x == p else x) for x in sublist] for sublist in self.coefficient_idx]
                    dellist.append(i)
                else:
                    self.num_factor = 0.0
                    return
                    
            elif p not in self.summation and q in self.summation:
#                print('a3', self.orbital_type)
                if p_type == q_type:
                    self.summation.remove(q)
                    self.coefficient_idx = [[(p if x == q else x) for x in sublist] for sublist in self.coefficient_idx]
                    dellist.append(i)
                else:
                    self.num_factor = 0.0
                    return
            else:
#                print('a4', self.orbital_type)
                if p_type != q_type:
                    
                    self.num_factor = 0.0
#                    print(self)
                    return

        new_delta = []
        for i, l  in enumerate(self.delta):
            if i not in dellist:
                new_delta.append(l)
#        print('self przed delta', self, new_delta)
        self.delta = new_delta
 #       print('self po delta', self, new_delta)
            
            
    def exec_delta(self, general_cond=False):
        """Perform substitution of indices / reduction of 
        summation indices resulting from Kronecker delta,
        \delta_{pq}. Numerical factor may be changed to zero 
        if p and q indices belong to disjoint sets.
        """

#        print('self przed', self)
        change_both = False
        delta = deepcopy(self.delta)
        ch = []
        for p, q in delta:
            leave = False
            # print('pq', p, q)
            for x in range(0, len(ch)):
                if p == ch[x][0]:
                    p = ch[x][1]
                if q == ch[x][0]:
                    q = ch[x][1]
            if not general_cond:
                if p in general and p not in self.summation:
                    print('GENERAL INDEX CANNOT BE FIXED - b')
                    sys.exit(1)
                if q in general and q not in self.summation:
                    print('COMPLETE INDEX CANNOT BE FIXED')
                    sys.exit(1)
            if p in virtualall and q in occupied:
                self.num_factor = 0.0
            elif p in occupied and q in virtualall:
                self.num_factor = 0.0
            elif p in virtual and q in CABS:
                self.num_factor = 0.0
            elif p in CABS and q in virtual:
                self.num_factor = 0.0
            elif p not in self.summation and \
                    q not in self.summation:
                #
                # Both indices fixed
                #
                pq = set([p,q])
                if not (set(virtualall) and pq == pq or set(occupied) and pq == pq \
                        or set(general) and pq == pq):
                    self.num_factor = 0.0
                else:
                    if p in completev and q in completev:
                        leave = True
            else:
                if p in self.summation and \
                        q not in self.summation:
                    #
                    # Summation over p; q is fixed
                    #
                    if q in completev and p in virtual:
                        leave = True
                    else:
                        changed = p
                        unchanged = q
                elif q in self.summation and \
                        p not in self.summation:
                    #
                    # Summation over q; p is fixed
                    #
                    if p in completev and q in virtual:
                        leave = True
                    else:
                        changed = q
                        unchanged = p
                else:
                    #
                    # There is summation over
                    # both indices
                    #

                    if p in generalall and q not in generalall:
                        changed = p
                        unchanged = q
                    elif p not in generalall and q in generalall:
                        changed = q
                        unchanged = p
                    if p not in generalall and q not in generalall:
                        changed = p
                        unchanged = q
                    # print(p, q, complete, completev)
                    # print(p in completev, q in complete)
                    if p in generalall and q in generalall:
                        if p in general and q in complete:
                            changed = q
                            unchanged = p
                        elif p in complete and q in general:
                            changed = p
                            unchanged = q
                        elif p in general and q in completev:
                            change_both = True
                        elif p in completev and q in general:
                            change_both = True
                        elif p in complete and q in completev:
                            changed = p
                            unchanged = q
                        elif p in completev and q in complete:
                            changed = q
                            unchanged = p
                        else:
                            changed = p
                            unchanged = q

                    # if p not in general:
                    #     changed = q
                    #     unchanged = p
                    # elif q not in general:
                    #     changed = p
                    #     unchanged = q
                    # else:
                    #     changed = p
                    #     unchanged = q

                if change_both == True:
                    if p in self.summation:
                        self.summation.remove(p)
                    if q in self.summation:
                        self.summation.remove(q)
                        
                    indices = self.indices()
                    banned = indices - set(self.summation)
                    banned = banned | set(fixed)

                    j = free_idx(virtual, banned)
                    self.summation.append(j)
                    if p in self.constraints and q in self.constraints:
                        self.constraints[j] = self.constraints[p] | self.constraints[q]
                        del self.constraints[p]
                        del self.constraints[q]
                    elif p in self.constraints and q not in self.constraints:
                        self.constraints[j] = self.constraints[p] 
                        del self.constraints[p]
                    elif p not in self.constraints and q in self.constraints:
                        self.constraints[j] = self.constraints[q]
                        del self.constraints[q]
                    del self.delta[self.delta.index([p, q])]
                    self.substitute(p, j)
                    self.substitute(q, j)
                    ch.append([p, j])
                    ch.append([q, j])

                else:
                    if leave == True:
                        continue
                    else:
                        if changed in self.summation:
                            self.summation.remove(changed)

                        if changed in self.constraints and unchanged in self.constraints:

                            #                    print('tutaj gowno', self, self.constraints)
                            #                    sys.exit(1)
                            self.constraints[unchanged]  = self.constraints[unchanged] | self.constraints[changed]
                            del self.constraints[changed]
                        #
                        # Remove delta symbol
                        #
                        del self.delta[self.delta.index([p, q])]
                        #
                        # Substitute CHANGED index for UNCHANGED
                        #
                        self.substitute(changed, unchanged)
                        ch.append([changed, unchanged])
        # print('self po   ', self)
        # print('')

    def substitute(self, p, q):
        """Substitute any instance of index p for index q.
        Substitution is performed in both summation and fixed
        indices subsets.
        """
        if p != q:
            if p in self.constraints and q not in self.constraints:
                self.constraints[q] = self.constraints.pop(p)
            elif p in self.constraints and q in self.constraints:
                self.constraints[q] = self.constraints[q] | self.constraints[p]
#                self.constraints[q] = list(self.constraints[q])
#                print(self.constraints[p])
#                self.constraints[q].append(list(self.constraints[p]))
 #               del self.constraints[p]
                #self.constraints[q] = self.constraints.pop(p)
        #         print('paldus_classes.substitute: cannot substitute two constrained indices')
        #         sys.exit(1)
        # print('a1', self, 'zamien', p, 'na', q)
        for x in range(0, len(self.coefficient)):
            for y in self.coefficient_idx[x]:
                if p in self.coefficient_idx[x]:
                    cidx = self.coefficient_idx[x].index(p)
                    self.coefficient_idx[x][cidx] = q
        # print('a2', self)
        for x in range(0, len(self.operator_idx)):
            for y in self.operator_idx[x]:
                if p in self.operator_idx[x]:
                    opidx = self.operator_idx[x].index(p)
                    self.operator_idx[x][opidx] = q
        # print('a3', self)
        for delta in self.delta:
            if delta[0] == p:
                delta[0] = q
            if delta[1] == p:
                delta[1] = q
        #
        # Search for redundant \delta_{qq}
        #
        # print('a4', self)
        for i in range(0, len(self.delta)):
            if self.delta[i] == [q, q]:
                del self.delta[i]
                break
        
        if p in self.summation:
            idx = self.summation.index(p)
            self.summation[idx] = q


    def rename_summation_with_restriction(self, restr1, restr2):

        excluded = restr1 + restr2
        refer = self.summation
        subst = []
        # print('w summation', self.summation)
        for x in self.summation:
            if x in virtual:
                fx = free_idx(virtual, excluded)
            elif x in occupied:
                fx = free_idx(occupied, excluded)
            elif x in general:
                fx = free_idx(general, excluded)
            else:
                print('Unknown idx - paldus_classes')
                sys.exit(1)
            subst.append(fx)
            # self.substitute(x, fx)
            excluded.append(fx)

        self.multisubst(refer, subst)
            

            
    def operator_string(self):
        """From given ugg creates arithmetic string, containg
        only operators i.e.
        \sum_aibi g_aibj E_ai E_bj --> (E_ai, E_bj)
        """
        temp = ugg()
        os = arithmetic_string()
        for i in range(0, len(self.operator_idx)):
            temp.operator_idx = []
            temp.operator_idx.append(self.operator_idx[i])
            temp.operator_type = []
            temp.operator_type.append(self.operator_type[i])
            os.append(temp)

        return(os)

    def left_split(self):
        """LEFT >> RIGHT = SELF, LEFT = E_{pq} """
        left = ugg()
        right = deepcopy(self)

        left.operator_idx.append(deepcopy(right.operator_idx[0]))
        left.operator_type.append(deepcopy(right.operator_type[0]))

        right.operator_idx.remove(right.operator_idx[0])
        right.operator_type.remove(right.operator_type[0])

        return left, right


    def right_split(self):
        """ LEFT << RIGHT = SELF, RIGHT = E_{pq} """

        left = deepcopy(self)
        right = ugg()

        s = len(left.operator_idx) - 1
        right.operator_idx.append(deepcopy(left.operator_idx[s]))
        right.operator_type.append(deepcopy(left.operator_type[s]))
        left.operator_idx.remove(left.operator_idx[s])
        left.operator_type.remove(left.operator_type[s])

        return left, right


    def permute(self, fixed_pairs):
        """Construct arithmetic string containing
        ugg objects differing from self ugg by a permutation
        of fixed indices. Permutations are restricted to a
        subset that does not disjoin pairs defined in fixed_pairs.
        """

        a = arithmetic_string()
        n = len(fixed_pairs)
        for p in permutations(fixed_pairs):
            x = deepcopy(self)
            substitutions = {}
            for i in range(0, n):
                substitutions[fixed_pairs[i][0]] = p[i][0]
                substitutions[fixed_pairs[i][1]] = p[i][1]

            for l in x.coefficient_idx:
                    for m in range(0, len(l)):
                        if l[m] in substitutions:
                            l[m] = substitutions[l[m]]

            for l in x.operator_idx:
                    for m in range(0, len(l)):
                        if l[m] in substitutions:
                            l[m] = substitutions[l[m]]

            for l in x.delta:
                    for m in range(0, len(l)):
                        if l[m] in substitutions:
                            l[m] = substitutions[l[m]]

            a.append(x)

        return a

def nswaps_lists(lst1, lst2):
    index_map = {v: i for i, v in enumerate(lst1)}
    
    lst2_indices = [index_map[v] for v in lst2]

    n_swaps = sum(1 for i in range(len(lst2_indices)) for j in range(i + 1, len(lst2_indices)) if lst2_indices[i] > lst2_indices[j])
    
    return n_swaps

def canonical_group(indices, signs):

#    print(indices, signs)
    n = len(indices)
#    print(n)
    perm = itertools.permutations(range(n))
#    print('perm', perm)
    cand_ind  = []
    cand_sign  = []
    csw = []
    original = list(range(0, n))
    for x in perm:
        nn = nswaps_lists(original, x)
        ml1 = []
        ml2 = []
        for i in x:
            ml1.append(indices[i])
            ml2.append(signs[i])
        cand_ind.append(ml1)
        cand_sign.append(ml2)

        csw.append(nn)
#        print('oto x', x, nn)
#    print('ci', cand_ind)
#    print('cs', cand_sign)

    return cand_ind, cand_sign, csw

def canonicalize_tensor(b, a, g, split):

    if g == DENS2PM:
        pattern =     ['+', '-', '+', '-']
        neg_pattern = ['-', '+', '-', '+']
    elif g == DENS2P:
        pattern =     ['+', '+', '+', '+']
        neg_pattern = ['-', '-', '-', '-']
    elif g == DENS3PM:
        pattern =     ['+', '+', '-', '+', '+', '-']
        neg_pattern = ['-', '-', '+', '-', '-', '+']
    elif g == DENS3P:
        pattern =     ['+', '+', '+', '+', '+', '+']
        neg_pattern = ['-', '-', '-', '-', '-', '-']
    elif g == DENS4PPM:
        pattern =     ['+', '+', '+', '-', '+', '+', '+', '-']
        neg_pattern = ['-', '-', '-', '+', '-', '-', '-', '+']
    elif g == DENS4PM:
        pattern =     ['+', '+', '-', '-', '+', '+', '-', '-']
        neg_pattern = ['-', '-', '+', '+', '-', '-', '+', '+']

#    print(b, split)
    b1, b2 = b[:split], b[split:]
    a1, a2 = pattern[:split], pattern[split:]

#    print(b1, b2)
#    print(a1, a2)

    can_b1, can_a1, fac1 = canonical_group(b1, a1)
    can_b2, can_a2, fac2 = canonical_group(b2, a2)

    candidates_idx = []
    candidates_ns = []
    n = len(can_b1)
    for i in range(0, n):
        for j in range(0, n):
            joins = can_a1[i]+ can_a2[j]
            if joins == pattern or joins == neg_pattern:

                ns = (fac1[i] + fac2[j])
                candidates_idx.append(can_b1[i]+can_b2[j])
                candidates_ns.append(ns)
            joins = can_a2[j]+ can_a1[i]
            if joins == pattern or joins == neg_pattern:
                ns = (fac1[i] + fac2[j])
                candidates_idx.append(can_b2[j]+can_b1[i])
                candidates_ns.append(ns)
#    print(candidates_idx)
#    print(candidates_ns)
#    print(min(candidates_idx))
    best = min(candidates_idx)
    min_index = candidates_idx.index(min(candidates_idx))
#    print(candidates_ns[min_index])

    num = (-1)**candidates_ns[min_index]
    
    return best, num

    
class arithmetic_string:
    
    def __init__(self, *args):
        self.string = []
        self.clusters = {}
        self.einstein = False
        for x in range(0, len(self)):
            self[x].einstein = self.einstein
        s = len(args)
        for i in range(0, s):
            self.string.append(deepcopy(args[i]))

    def __str__(self):
        s = ""
        if len(self) == 0:
            print("0")
        else:
            for x in self:
                s += x.__str__()
            return s


    # def __add__(self, other):
    #     a = deepcopy(self.string)
    #     b = deepcopy(other.string)
    #     c = arithmetic_string()
    #     c.string = a + b
    #     return c

    def __add__(self, other):
        c = arithmetic_string()
        c.string = self.string + other.string
        return c


    def __mul__(self, other):
        #
        # C <- A * B. A, B need not be disambiguated.
        #
        c = arithmetic_string()
        for xa in self:
            for xb in other:
                a = deepcopy(xa)
                b = deepcopy(xb)
                disambiguate(a, b)
                c.append(a.fromright(b))
        return c

    
    def __eq__(self, other):
        """Determine if given arithmetic_string 
        instances are *identical*.
        """
        if len(self) == len(other):
            excluded = []
            equal = True

            for i in range(len(self)):
                a = self[i]
                l = False
                for j in range(len(other)):
                    b = other[j]
                    na = a.num_factor
                    nb = b.num_factor
                    if a == b and na == nb and \
                            j not in excluded:
                        l = True
                        excluded.append(j)
                        break
                if not l:
                    equal = False
                    break
        else:
            equal = False

        return equal
                


    def __contains__(self, item):
        """Check if arithmetic_string contains ugg object item.
        Ugg objects which differ merely by numerical factor
        are considered *different*.
        """
        
        contains = False
        for x in self:
            if x == item and x.num_factor == item.num_factor:
                contains = True
                break

        return contains

        
    def __next__(self):
        for x in self.string:
            yield x


    def __iter__(self):
        return self.__next__()


    def __len__(self):
        return len(self.string)


    def true_len(self):
        tr_len = 0
        for x in self:
            if len(x.permutation_op) > 0:
                tr_len += math.factorial(len(x.permutation_op))
            else:
                tr_len += 1

        return tr_len    

    
    def __getitem__(self, key):
        return self.string[key]

    
    def __setitem__(self, key, value):
        self.string[key] = value


    def __delitem__(self, key):
        del self.string[key]


    def scale(self, s):
        c = deepcopy(self)
        for x in c:
            x.num_factor *= s
        return c

    

    def establish_fixed(self, fixed_fx=[]):
        """ establish fixed indices for whole
        arithmetic string
        """
        fixed = []
        for x in self:
            x.establish_fixed(fixed_fx)

    def clear_fixed(self):
        """ establish fixed indices for whole
        arithmetic string
        """
        for x in self:
            x.clear_fixed()

    def disambiguate(self):
        #print('taktu')
        disambiguate(*self.string)     


    def fromleft(self, e):
        """Concatenate from left: E >> SELF"""
        a = deepcopy(self)
        for x in range(0, len(a.string)):
            a.string[x] = a.string[x].fromleft(e)
        return a


    def fromright(self, e):
        """Concatenate from right: SELF << E"""
        a = deepcopy(self)        
        for x in range(0, len(a.string)):
            a.string[x] = a.string[x].fromright(e)
        return a

    def append(self, item):
        self.string.append(deepcopy(item))


    def as_deepcopy(self):
        result = deepcopy(self)
        return result
    

    def cleanup(self):
        """Remove items with zero numerical factor."""
        remove = []
        for k in range(0, len(self)):
            if abs(self[k].num_factor) < EPSILON:
                remove.append(k)
            
        if len(remove) > 0:
            n = 0
            for k in remove:
                del self[k - n]
                n = n + 1


    def transpose(self):
        """
        Transpose ugg operators.
        """
        for x in self:
            x.transpose()


    def exec_delta_cas(self):

        for e in self:
            e.exec_delta_cas()
        
        
    def exec_delta(self, general_cond=False):
        """
        Execute Kronecker deltas in every
        ugg in self.
        """
        
        i = 0
        for e in self:
            i += 1
            if general_cond:
                e.exec_delta(True)
            else:
                e.exec_delta()
        self.cleanup()

    def optimize(self):
        
        for x in self:
            x.optimize()
            
    def optimize_f12(self):
                    
        for x in self:
            x.optimize_f12()

    def exec_delta_fixed(self, nonzero):
        """
        Execute Kronecker deltas in every
        ugg in self.
        """
        for e in self:
            e.exec_delta_fixed(nonzero)
        self.cleanup()

    def clear_deltas(self):
        for i in self:
            i.clear_deltas()
        self.cleanup()


    def as_standarize(self, last_idx = '', fixed_fx=[], cas = False, drag_fixed=[]):
        # t1 = time.time()
        for x in range(0, len(self)):

            if self[x].num_factor != 0.0:
                if not cas:
                    self[x].exec_delta()
                # t1a = time.time()
                self[x].standarize(last_idx, cas=cas, drag_fixed=drag_fixed)
                # t2a = time.time()
                # print('mid-as_stand', t2a-t1a)


        # t2 = time.time()
        # print('as_stand', t2-t1, len(self))



    def as_standarize_all(self):
        result = arithmetic_string()
        q = False
        j = 0
        while q!= True:
            temp = deepcopy(self)
            i = 0
            for x in range(0, len(self)):
                if self[x].num_factor != 0.0:
                    i +=1
                    result  = self[x].standarize_all()
                    self = self + result
            s2 = self.exec_delta()
            self = s2
            j += 1
            if self == temp:
                q = True
        return self, result


    def excitation_rank(self, er):
        result = arithmetic_string()

        for x in self:
            rank = x.excitation_rank()
            if rank == er:
                result.append(x)
        return result


    def cluster(self):
        self.clusters = {}

        k = 0
        for x in self:
            k += 1
            x.hash()
            h = x.hash_tuple
 #           if h in self.clusters:
            if h in list(self.clusters.keys()):
                self.clusters[h].append(x)
            else:
                self.clusters[h] = arithmetic_string(x)


    def cluster_for_fortran(self):
        """ from given result (self) gives creates
        list of clusters, which have the same summation indices
        in the same order. 
        """
        banned = []
        self.clusters = {}
        for x in range(0, len(self)):
            if x not in banned:
                self.clusters[x] = arithmetic_string(self[x])
                for y in range(x +1, len(self)):
                    if y not in banned:
                        if self[y].summation == self[x].summation:
                            self.clusters[x].append(self[y])
                            banned.append(y)

    def cluster_for_fortran_f12(self, all_interm_dict):
        """ from given result (self) gives creates
        list of clusters, which have the same summation indices
        in the same order. The clusters are divided in blocks, depending
        on how much intermediates are used in a cluster. If one need to load
        intermediates witch memory bigger than the thershold, cluster is divided.
        """
        banned = []
        self.clusters = {}
        k = -1
        for x in range(0, len(self)):
            if x not in banned:
                k += 1
                z = 1
                self.clusters[k] = arithmetic_string(self[k])
                for y in range(x +1, len(self)):
                    if y not in banned:
                        if self[y].summation == self[k].summation:
                            if z==3:
                                k += 1                                                                
                                self.clusters[k] = arithmetic_string(self[y])
                                z = 1
                            else:
                                self.clusters[k].append(self[y])
                                banned.append(y)
                                z += 1



    def cluster_for_fortran_f12(self, all_interm_dict):
        """ from given result (self) gives creates
        list of clusters, which have the same summation indices
        in the same order. The clusters are divided in blocks, depending
        on how many intermediates are used in a cluster. If one need to load
        intermediates witch memory bigger than the thershold, cluster is divided 
        into subgroups.
        """

        print('clusters for fortran f12')
        banned = []
        self.clusters = {}
        k = -1
        self.clusters[0] = arithmetic_string(self[0])
        print('poczatkowy cluster to', self.clusters[0])
        for x in range(0, len(self)):
            if x not in banned:
                print('i jade od nowa')
                k += 1
                z = 1
                self.clusters[k] = arithmetic_string(self[x])
                print(k, self[x], 'taki jest nowy poczatkowy')
                disk_coef_name = []
                for coef in self[x].coefficient:
                    if 'l' in coef:
                        coef_name = extract_coef(coef)
                        if all_interm_dict[coef_name]['mem'] == 'disk':
                            disk_coef_name.append(deepcopy(coef_name))                            
                if len(disk_coef_name) != 0:
                    print('jest czytanie z dysku?', disk_coef_name)
                    find_all_disk_interm = True
                    n = 0
                    coef_name = disk_coef_name[0]
                    while find_all_disk_interm:
                    # for coef_name in disk_coef_name:
                        print('teraz szukam innego wyrazu z tym wspolczynnikiem', coef_name)
                        for y in range(x +1, len(self)):
                            if y not in banned:
                                if self[y].summation == self[x].summation:
                                    print(y, 'mam kandydata y', self[y])
                                    for coefy in self[y].coefficient:
                                        if 'l' in coefy:
                                            coef_namey = extract_coef(coefy)
                                            if coef_name == coef_namey:
                                                print('dodaje wyraz ktory tez ma wspolczynnik')
                                                # print('dodaje, k i zwiekszam z', k, y)
                                                self.clusters[k].append(self[y])
                                                banned.append(y)
                                                z += 1
                                                for coefd in self[y].coefficient:
                                                    if 'l' in coefd:
                                                        coef_named = extract_coef(coefd)
                                                    if all_interm_dict[coef_named]['mem'] == 'disk':
                                                        if coef_named not in disk_coef_name:
                                                            disk_coef_name.append(deepcopy(coef_named))
                                                            print('dodaje na liste disk_coef')
                                        else:
                                            print('nie jest zapisany na dysku')
                        n += 1 # n numeruje petle while
                        if n == len(disk_coef_name):
                            find_all_disk_interm = False
                print('banned sa', banned)
                for y in range(x +1, len(self)):
                    print('y', y)
                    if y not in banned:
                        print('not in')
                        print('tak summ', z, self[y], self[x])
                        if self[y].summation == self[x].summation:
                            if z >=4:
                                print('maja rowne summation, ale to juz nowy interm bedzie')
                                # k += 1                                                                
                                # self.clusters[k] = arithmetic_string(self[y])
                                # print('dodaje go jako nowy poczatek', k, self[y])
                                break
                                z = 1
                            else:
                                print('dodaje do danego klastera', self[y])
                                self.clusters[k].append(self[y])
                                banned.append(y)

                                disk_coef_name = []
                                for coef in self[y].coefficient:
                                    if 'l' in coef:
                                        coef_name = extract_coef(coef)
                                        if all_interm_dict[coef_name]['mem'] == 'disk':
                                            disk_coef_name.append(deepcopy(coef_name))                            
                                if len(disk_coef_name) != 0:
                                    print('jest czytanie z dysku?', disk_coef_name)
                                    z += 1
                                    find_all_disk_interm = True
                                    n = 0
                                    coef_name = disk_coef_name[0]
                                    while find_all_disk_interm:
                                        print('teraz szukam innego wyrazu z tym wspolczynnikiem',n, len(disk_coef_name), coef_name)
                                        for p in range(y + 1, len(self)):
                                            if p not in banned:
                                                if self[p].summation == self[y].summation:
                                                    print(p, 'mam kandydata p', self[p])
                                                    for coefp in self[p].coefficient:
                                                        if 'l' in coefp:
                                                            coef_namep = extract_coef(coefp)
                                                            print('nowy coef z l', coef_namep)
                                                            
                                                            if coef_name == coef_namep:
                                                                print('dodaje wyraz ktory tez ma wspolczynnik')
                                                                # print('dodaje, k i zwiekszam z', k, y)
                                                                self.clusters[k].append(self[p])
                                                                banned.append(p)
                                                                z += 1
                                                                for coefd in self[p].coefficient:
                                                                    if 'l' in coefd:
                                                                        coef_named = extract_coef(coefd)
                                                                        if all_interm_dict[coef_named]['mem'] == 'disk':
                                                                            if coef_named not in disk_coef_name:
                                                                                disk_coef_name.append(deepcopy(coef_named))
                                                                                print('dodaje na liste disk_coef')
                                                        else:
                                                            print('nie jest zapisany na dysku')
                                        n += 1 # n numeruje petle while
                                        print(disk_coef_name)
                                        if n == len(disk_coef_name):
                                            find_all_disk_interm = False





                                    
                                # # z += 1
                                # for coef in self[y].coefficient:
                                #     if 'l' in coef:
                                #         coef_name = extract_coef(coef)
                                #         disk_coef_name = []
                                #         if all_interm_dict[coef_name]['mem'] == 'disk':
                                #             print('w nowododanym jest interm z disk', coef)
                                #             z += 1
                                #             disk_coef_name.append(deepcopy(coef_name))
                                #             if len(disk_coef_name) != 0:
                                #                 print('nowododany ma interm z disk')
                                #                 for coef_name in disk_coef_name:
                                #                     for m in range(x +1, len(self)):
                                #                         if m not in banned:
                                #                             if self[m].summation == self[x].summation:
                                #                                 for coefm in self[m].coefficient:
                                #                                     if 'l' in coefm:                                                                    
                                #                                         coef_namem = extract_coef(coefm)
                                #                                         if coef_name == coef_namem:
                                #                                             print('tak, dodaje')
                                #                                             self.clusters[k].append(self[m])
                                #                                             z += 1
                                #                                             banned.append(m)
                                   
                # for y in range(x +1, len(self)):
                #     if y not in banned:
                #         if self[y].summation == self[k].summation:
                #             if z==3:
                #                 k += 1                                                                
                #                 self.clusters[k] = arithmetic_string(self[y])
                #                 z = 1
                #             else:
                #                 self.clusters[k].append(self[y])
                #                 banned.append(y)
                #                 z += 1



    def cluster_for_fortran_mem(self, name_list, cluster_name):
        """ from given result (self) gives creates                                                                                             
        list of clusters, which have the same summation indices                                                                                
        in the same order.                                                                                                                 
        """
        banned = []        
        self.clusters = {}
        for x in range(0, len(self)):
            if x not in banned:
                self.clusters[x] = arithmetic_string(self[x])
                cluster_name.append(name_list[x])
                for y in range(x +1, len(self)):
                    if y not in banned:
                        if self[y].summation == self[x].summation:
                            if name_list[y] == name_list[x]:
                                self.clusters[x].append(self[y])
                                banned.append(y)


    def compress(self, selector):
        s = len(self)
        for x in range(s - 1, 0, -1):
            if selector[x] == 0:
                del self[x]


    def expand(self):
        """ input: self with permutation operators, output:
        all permutation operators expanded. Resulted lenght: longer"""
        res_new = arithmetic_string()
        for x in self:
            if len(x.permutation_op) != 0:
                a = x.permute(x.permutation_op)
                for y in a:
                    y.permutation_op = []
                for k in a:
                    res_new.append(k)
            else:
                a = x
                res_new.append(a)

        return res_new


    def integrate(self, bra = [], ket = [], braspin = [], ketspin = [], nodelta=False):
        s = arithmetic_string()
        for x in self:
            d =  x.integrate(bra, ket, braspin, ketspin, nodelta=False)
            if len(d) > 0:
                print('sss', d)
            else:
                print('sss', 0)
            s = s + d
        return s


    def cabstransform(self):
        s = arithmetic_string()
        for x in self:
            s = s + x.cabstransform()

        # because regular amplitudes t2 can only have virtual indices ab
        # and cannot have (by definition) aA, Aa or AB, I'm removing
        # all terms that have t^aA, t^Aa, t^AB. situation is different
        # for the S amplitudes, where S = T + T', therefore S^aA
        # is in fact T'aA
        
        for x in s:
            for coef_no in range(0, len(x.coefficient)):
                if x.coefficient[coef_no] == CC_AMPLITUDE:                    
                    no_cabs = 0
                    for k in x.coefficient_idx[coef_no]:
                        if k in CABS:
                            x.num_factor = 0.0
                            break

        # I'm removing all instances of ff{akbl}, where a and b are both wirtual indices
        # According to definition of F in Shiozaki article, these instances are equal to zero
        for x in s:
            for j in range(0, len(x.coefficient)):
                y = x.coefficient[j]
                if y == F12_TWOEL or y ==F12_TWOEL_COMBO:
                    n_of_virt = 0
                    for z in x.coefficient_idx[j]:
                        if z in virtual:
                            n_of_virt +=1
                    if n_of_virt >=2:
                        x.num_factor =0

        return s



    def integrate_wick(self, bra = [], ket = [], braspin = [], ketspin = []):
        s = arithmetic_string()
        for x in self:
            s = s + x.integrate_wick(bra, ket, braspin, ketspin)
        return s


    def multisubst(self, refer, subst):
        """
        Changes all indices in refer list to corresponding indices
        in subst. 
        """
        # print(refer, subst)
        if len(refer) != len(subst):
            print("Error. Nonequal lenghts of reference and permuted lists.")
            return
        for i in range(0, len(self)):
            self[i].multisubst(refer, subst)

def free_idx(pool, excluded):
    # print('ejstem w free')
    for y in pool:
        #print(y, excluded)
        if y not in excluded:
            return y

def idx_type(idx):

    if idx in virtual:
        return virtual
    elif idx in occupied:
        return occupied
    elif idx in CABS:
        return CABS
    elif idx in complete:
        return complete
    elif idx in completev:
        return completev
    elif idx in mona_occup_original:
        return mona_occup#_original
    elif idx in monb_occup_original:
        return monb_occup#_original
    elif idx in mona_virt_original:
        return mona_virt#_original
    elif idx in monb_virt_original:
        return  monb_virt#_original
    elif idx in mona_occup:
        return mona_occup
    elif idx in monb_occup:
        return monb_occup
    elif idx in mona_virt:
        return mona_virt
    elif idx in monb_virt:
        return  monb_virt
    elif idx in general:
        return  general

    else:
        return None



def basic_disambiguate(e1, e2, all_fixed):
    """Perform index substitutions in e1 so that there are
    no common summation indices in e1 and e2, and summation
    indices in e2 are different from fixed indices in e1 and
    e2.
    """ 
#    print('e', e1, e2)
    for x in e2.summation:
        se1 = e1.indices()
        se2 = e2.indices()
        #print('sese', x)
#        print(se1, se2, a)
        #print(se1 | se2|all_fixed)
        if x in se1:
            #print('sprawdzam', x)
            #
            # Substitute x for y
            #
            if x in virtual:
                y = free_idx(virtual, se1 | se2|all_fixed)
                e2.substitute(x, y)
            elif x in occupied:
                # print(x, occupied, se1 | se2|all_fixed)
                y = free_idx(occupied, se1 | se2|all_fixed)
                # print('oto', y)
                e2.substitute(x, y)
            elif x in completev:
                y = free_idx(completev, se1 | se2|all_fixed)
                e2.substitute(x, y)
            elif x in complete:
                y = free_idx(complete, se1 | se2|all_fixed)
                e2.substitute(x, y)
            elif x in general:
                y = free_idx(general, se1 | se2|all_fixed)
                e2.substitute(x, y)
            elif x in CABS:
                y = free_idx(CABS, se1 | se2 )
                e2.substitute(x, y)

    #print(e1, e2)
    #print('wychodze')
def disambiguate_with_list(e1, idx_list):
    """Perform index substitutions in e1 so that 
    in e1 there are no indices from idx_list
    """ 

    for x in e1.summation:
        se1 = e1.indices()
        se2 = set(idx_list)

        if x in se2:
            #
            # Substitute x for y
            #
            if x in virtual:
                y = free_idx(virtual, se1 | se2)
                e1.substitute(x, y)
            elif x in occupied:
                y = free_idx(occupied, se1 | se2)
                e1.substitute(x, y)
            elif x in completev:
                y = free_idx(completev, se1 | se2)
                e1.substitute(x, y)
            elif x in complete:
                y = free_idx(complete, se1 | se2)
                e1.substitute(x, y)
            elif x in general:
                y = free_idx(general, se1 | se2)
                e1.substitute(x, y)

def disambiguate(*args):
    s = len(args)
#    print('na wejsciu do disambiguate', s)#, args)

    all_fixed = set()
    for f in args:
        c = f.set_fixed()
        if len(c.intersection(all_fixed)) != 0:
            print('please change the names of the fixed indices in', f)
#            sys.exit(0)
        else:
            all_fixed = all_fixed | deepcopy(c)

        #print(f)
    #print('all fixed', all_fixed)
    if s < 2:
        #print('s < 2')
        return
    elif s == 2:
        #print('s =2')
        basic_disambiguate(args[0], args[1], all_fixed)
        basic_disambiguate(args[1], args[0], all_fixed)
    else:
        #print('s > 2')
        #print('kle')
        a = deepcopy(args[0])
        #print('kra', a)
        for k in range(0, s - 1):
            #print('kkk-tu', k, 'z', s-1)
            #print(a, args[k + 1])
            basic_disambiguate(a, args[k + 1], all_fixed)
            a = args[k + 1].fromright(a)
            #print('po fromriggh', a, args[0], args[1], args[2])
        a = deepcopy(args[s - 1])
        #print('')
        #print('i zaczynam od tylu', a)
        #print('')
        for k in range(s - 1, 0, -1):
            #print('kkk', k, a, args[k - 1])
            basic_disambiguate(a, args[k - 1], all_fixed)
            a = args[k - 1].fromleft(a)
            #print(a, 'po fromleft')

            
def disambiguate_two_list(e1, e2, idx_list):

    basic_disambiguate(e1, e2)
    basic_disambiguate(e2, e1)


def pair_permutations(list1, list2):
    """ makes all possible permutations of list2, and then appends all of the permutations with list1"""
    list1 = list(list1)
    list2 = list(list2)
    if len(list1) != len(list2):
        print("ERROR: Number of virtual indices does not equal number of occupied indices.")

    s1 = list(permutations(list2, len(list2)))

    s = deque()
    for x in range(0, len(s1)):
        s.append([])
    for x in range(0, len(s1)):
        for y in range(0, len(s1[x])):
            s[x].append([list1[y],s1[x][y]])
    return s


def bubblesort(l):
    """Sort list l in ascending order. Return
    number of elementary interchanges of indices.
    """
    if len(l) < 2:
        return 0

    done = False
    n = 0
    while not done:
        done = True
        for k in range(len(l) - 1):
            if l[k] > l[k + 1]:
                t = l[k + 1]
                l[k + 1] = l[k]
                l[k] = t
                n += 1
                done = False

    return n


def same_subset(p, q):
    if (p in virtual) and (q in virtual):
        l = True
    elif (p in occupied) and (q in occupied):
        l = True
    elif (p in general) and (q in general):
        l = True
    else:
        l = False

    return l


def compare_coeffs(name_a, idx_a, name_b, idx_b):
    """Compare values of coefficients. Permutational symmetry
    of electron-repulsion integrals is taken into account. Antisymmetry
    of excitation amplitudes is taken into account.

    Returned values:
     1 -- Coefficients are identical
    -1 -- Coefficients differ by sign
     0 -- Coefficients are not equal under any permutation
    """

    idx_a0 = deepcopy(idx_a)
    idx_b0 = deepcopy(idx_b)
    idx_a0.sort()
    idx_b0.sort()

    if (name_a == name_b) and (idx_a0 == idx_b0):
        if name_a == TWOEL_INT:
            #
            # Electron repulsion integral
            #
            sorted_a = [idx_a[0:2], idx_a[2:4]]
            sorted_a[0].sort()
            sorted_a[1].sort()
            sorted_a.sort()

            sorted_b = [idx_b[0:2], idx_b[2:4]]
            sorted_b[0].sort()
            sorted_b[1].sort()
            sorted_b.sort()

            if sorted_a == sorted_b:
                return 1
            else:
                return 0

        elif name_a == CC_AMPLITUDE:
            #
            # Orbital excitation amplitude
            # Assumption:
            # t^{abc}_{ijk} = t^{bac}_{jik} = t^{bca}_{jki} = ...
            # (indices can be interchanged provided that pairs
            # belonging to the same E_{ia} operator are not
            # disjoined)
            #
            pairs_a = []
            pairs_b = []
            #
            # t^{ab}_{ij} => ["a", "i", "b", "j"]
            #
            for k in range(0, len(idx_a) // 2):
                virt_a = idx_a[2 * k]
                occ_a  = idx_a[2 * k + 1]
                pairs_a.append((virt_a, occ_a))

                virt_b = idx_b[2 * k]
                occ_b  = idx_b[2 * k + 1]
                pairs_b.append((virt_b, occ_b))

            pairs_a.sort()
            pairs_b.sort()

            if pairs_a == pairs_b:
                return 1
            else:
                return 0

        elif name_a == S_AMPLITUDE:
            #
            # Orbital excitation amplitude
            # Assumption:
            # t^{abc}_{ijk} = t^{bac}_{jik} = t^{bca}_{jki} = ...
            # (indices can be interchanged provided that pairs
            # belonging to the same E_{ia} operator are not
            # disjoined)
            #
            pairs_a = []
            pairs_b = []
            #
            # t^{ab}_{ij} => ["a", "i", "b", "j"]
            #
            for k in range(0, len(idx_a) // 2):
                virt_a = idx_a[2 * k]
                occ_a  = idx_a[2 * k + 1]
                pairs_a.append((virt_a, occ_a))

                virt_b = idx_b[2 * k]
                occ_b  = idx_b[2 * k + 1]
                pairs_b.append((virt_b, occ_b))

            pairs_a.sort()
            pairs_b.sort()

            if pairs_a == pairs_b:
                return 1
            else:
                return 0
        
        elif name_a in SYMMETRIC_MATRIX:
            sorted_a = deepcopy(idx_a)
            sorted_b = deepcopy(idx_b)
            if sorted_a == sorted_b:
                return 1
            else:
                return 0

        else:
            if idx_a == idx_b:
                return 1
            else:
                return 0
    else:
        return 0

#@jit
def locate_epa(a):
    """
    Return index of the first E_{ia} or E_{ba} operator
    to the right. If the operator is not found, return -1.
    """
    r = -1
    for i in range(len(a.operator_idx)):
        if a.operator_idx[i][1] in virtualall:
            r = i    

    return r

#@jit
def moveright(a, nodelta = False):
    """
    epa = E_{pa} or T_{pa} - excitation operator
    with first index general, and second virtual.

    Check if epa is acting directly on the
    reference determinant.

    If epa is in the middle, move it to the
    first place on the right,
    using swap method, and return all contracted
    terms produced by the E (T) operator
    comutation rules.
    
    """
    ugg_order = len(a.operator_idx)
    pos = locate_epa(a)
    if pos == ugg_order - 1:
        #
        # Virtual annihilator acting
        # directly on reference determinant
        #
        empty = arithmetic_string()
        # print('Virtual annihilator acting directly on reference determinant')
        return empty
    elif pos == -1:
        #
        # No instances of E_{ba} or
        # E_{ia} operators. Do not
        # change anything
        # #
        # print('No instances of E_{ba} or E_{ia} operators. Do not change anything')
        return arithmetic_string(a)
    else:
        # print('MOVE')
        contracted = arithmetic_string()
        transposed = deepcopy(a)
        # print('transp', transposed)
        for k in range(pos + 1, ugg_order):
            c, transposed = transposed.swap(k)
            # if len(c) > 0:
            #     print('ctra', c)
            # print('ctrb',transposed)
            if not nodelta:
                c.exec_delta()
            # print('dupa')
            contracted = contracted + c

        return contracted

#@jit
def projectout_epa(s, nodelta = False):


    old = deepcopy(s)
#    print('')
#    print('old-beg', old)
#    print('')
    l = True
    new = arithmetic_string()
    while l:
        new = arithmetic_string()
 #       kkk = 0
        for x in old:
#            print('x', x)
  #          kkk += 1
            move = moveright(x, nodelta)
            if len(move) !=0:
                new = new + move
#            print('lnln', len(move))
#        print(kkk)
        old = deepcopy(new)
#        print('new', new)
        pos = -1
        for x in new:
            pos = max(pos, locate_epa(x))


        if pos == -1:
            l = False
        else:
            l = True
 #   print('old', old)
    return old

#@jit
def projectout(e, nodelta=False):
    """
    Return arithmetic_string containing ugg instances
    that do not contain E_{ia} or E_{ba} unitary group
    generators. Ugg instances with E_{ia} or E_{ba} 
    operators acting on ket vectors are projected out.
    Transpose and project out E_{ai} operators acting on
    bra vector.
    """
    s0 = arithmetic_string(e)
#    print('projectout, s0', s0)
    s1 = projectout_epa(s0, nodelta)
    # if len(s1)>0:
    #     print('s1', s1)
    s1.transpose()
    s2 = projectout_epa(s1, nodelta)
    s2.transpose()
    # if len(s2)>0:
    #     print('s2', s2)
    return s2

#@jit
def intprojected(e):
    """
    Calculate <0| e |0> integral, where
    e does not contain E_{ia} or E_{ba}
    unitary group generators.
    """
    #
    # Search for index of virtual orbital
    #
    for p, q in e.operator_idx:
        if p in virtualall or q in virtualall:
            print("ERROR: INVALID ARGUMENT.")
            print("UGG INSTANCE CONTAINING VIRTUAL ORBITAL INDEX.")
            sys.exit(1)


    f = deepcopy(e)
    if len(f.operator_idx) == 0:
        #
        # No operators acting on reference HF state
        #
        return arithmetic_string(f)
    else:
        #
        # String of E_{ij} operators:
        # E_{ij} E_{kl} ... E_{mn}.
        # Every unitary group generator contributes
        # 2 \delta_{ij} to the <0| . |0> integral
        #
        if "t0" in f.operator_type:
            empty = arithmetic_string()
            return empty
        else:
            for ijpair in f.operator_idx:
                i = ijpair[0]
                j = ijpair[1]
                f.new_delta(i, j)
                f.num_factor *= 2

            f.operator_idx = []
            f.operator_type = []

            return arithmetic_string(f)

def integrate_wick(e0, bra = [], ket = [], braspin = [], ketspin = []):
    """Calculate <0| e |0> integral, where e is ugg object,
    and |0> denotes HF reference state. Result of the
    integration is represented as an arithmetic string.
    By default, integral is computed with nonexcited
    determinants. Nondiagonal integral with excited
    determinants:
    
    integrate(e, bra = ["a", "i", "b", "j", "c", "k"], \
    ket = ["a", "i", "b", "j"]
    => <\Psi^{abc}_{ijk} | e | \Psi^{ab}_{ij}

    """
    #
    # Represent general summation indices as
    # occupied + virtual. This algorithm does
    # not work with general indices
    #
    genidx = False
    for p in e0.summation:
        if p in general:
            genidx = True
            break
    if genidx:
        s = e0.ovsplit()
        ss = 1
        for x in s:
            ss += 1
        return s.integrate_wick(bra, ket, braspin, ketspin)


    e = deepcopy(e0)
    er = e.excitation_rank()
    bra_ugg = ugg()
    ket_ugg = ugg()
    if bra != []:
        for i in range(0, len(bra)-1, 2):
            bra_ugg.operator_idx.append([bra[i], bra[i+1]])
            bra_ugg.operator_type.append(braspin[i//2])

    if ket != []:
        for i in range(0, len(ket)-1, 2):
            ket_ugg.operator_idx.append([ket[i], ket[i+1]])
            ket_ugg.operator_type.append(ketspin[i//2])

    er_bra = bra_ugg.excitation_rank()
    er_ket = ket_ugg.excitation_rank()

    if er - er_bra + er_ket != 0:
        empty = arithmetic_string()
        return empty
    #
    # Integral with excited determinants
    #

    if bra != []:
        x = ugg()
        for k in range(0, len(bra) // 2):
            a = bra[2 * k]
            i = bra[2 * k + 1]
            #
            # Transpose and revert order
            #
            x.operator_idx.insert(0, [i, a])
            x.operator_type.insert(0, braspin[k])
        disambiguate(e, x)
        e = e.fromleft(x)
    if ket != []:
        y = ugg()
        for k in range(0, len(ket) // 2):
            a = ket[2 * k]
            i = ket[2 * k + 1]
            y.operator_idx.append([a, i])
            y.operator_type.append(ketspin[k])
        
        disambiguate(e, y)
        e = e.fromright(y)

    nonzero_contractions = integrate_wick_basic(e.operator_idx, e.operator_type)

    integral = arithmetic_string()
    for i in range(0, len(nonzero_contractions)):
        e_temp = deepcopy(e)
        e_temp.operator_idx = []
        e_temp.operator_type = []
        e_temp.delta = nonzero_contractions[i]['delta']
        e_temp.num_factor *= float(nonzero_contractions[i]['n_factor'])
        integral = integral + arithmetic_string(e_temp)

    return integral
                    
def cabstransform(res):

    complidx = False
#    print('kgowno')
    #    result = arithmetic_string()
    # print('reserre', res)

    for p in res.coefficient_idx:
        for q in p:
            if q in completev:
                complidx = True
                break
    
    for p in res.summation:
        # print('p', p)
        if p in completev:
            complidx = True
            break
    if complidx:
        print('tak, complix', complidx)
        s  = res.cabssplit()
        return s.cabstransform()
    else:
        # print('zwracam res', res)
        return arithmetic_string(res)

#@jit
def integrate(e0, bra = [], ket = [], braspin = [], ketspin = [], nodelta=False):
    """Calculate <0| e |0> integral, where e is ugg object,
    and |0> denotes HF reference state. Result of the
    integration is represented as an arithmetic string.
    By default, integral is computed with nonexcited
    determinants. Nondiagonal integral with excited
    determinants:
    
    integrate(e, bra = ["a", "i", "b", "j", "c", "k"], \
    ket = ["a", "i", "b", "j"]
    => <\Psi^{abc}_{ijk} | e | \Psi^{ab}_{ij}

    """
    print('e0', e0, ket)
    #
    # Represent general summation indices as
    # occupied + virtual. This algorithm does
    # not work with general indices
    #
    genidx = False
    for p in e0.summation:
        if p in general or p in complete:
            genidx = True
            break
    if genidx:
        s = e0.ovsplit()
        ss = 1
        for x in s:
            ss += 1
        return s.integrate(bra, ket, braspin, ketspin)


    e = deepcopy(e0)
#    print('po split', e)

    #
    # Check excitation rank
    #
    er = e.excitation_rank()
#    print('er', er)
    bra_ugg = ugg()
    ket_ugg = ugg()
    if bra != []:
        for i in range(0, len(bra)-1, 2):
            bra_ugg.operator_idx.append([bra[i], bra[i+1]])
            bra_ugg.operator_type.append(braspin[i//2])

    if ket != []:
        for i in range(0, len(ket)-1, 2):
            ket_ugg.operator_idx.append([ket[i], ket[i+1]])
            ket_ugg.operator_type.append(ketspin[i//2])
    er_bra = bra_ugg.excitation_rank()
    er_ket = ket_ugg.excitation_rank()
    # print(bra, ket)
    # print(bra_ugg)
    # print(ket_ugg)
    # print('e2', e, er_bra, er_ket, er, er - er_bra + er_ket)
    if er - er_bra + er_ket != 0:
        #    if er - len(bra) // 2 + len(ket) // 2 != 0:
        #        print(e0, bra, ket, bra_ugg, ket_ugg)
        #        print('tutaj1', er, len(bra)//2, len(ket)//2, er - len(bra) // 2 + len(ket) // 2)
        #        print('tutaj2', er, er_bra, er_ket, er - er_bra + er_ket) 

        empty = arithmetic_string()
        return empty
#    print('e1', e)
    #
    # Integral with excited determinants
    #
    if bra != []:
        x = ugg()
        for k in range(0, len(bra) // 2):
            a = bra[2 * k]
            i = bra[2 * k + 1]
            #
            # Transpose and revert order
            #
            x.operator_idx.insert(0, [i, a])
            x.operator_type.insert(0, braspin[k])
        disambiguate(e, x)
        e = e.fromleft(x)
    if ket != []:
        y = ugg()
        for k in range(0, len(ket) // 2):
            a = ket[2 * k]
            i = ket[2 * k + 1]
            y.operator_idx.append([a, i])
            y.operator_type.append(ketspin[k])
        
        disambiguate(e, y)
        e = e.fromright(y)
#    print('eeeeee', e)
    #
    # Project out E_{ia}, E_{ai}, and E_{ba}
    # operators
    #
    projected = projectout(e, nodelta)
    # print('proj?')
    # if len(projected)>0:
    #     print('proj', projected)
    k = 1
    for x in projected:
        k += 1

    integral = arithmetic_string()

    k = 1
    for x in projected:
#        print('x', x)
        k+= 1
        integral2 = intprojected(x)
#        print('int2', integral2)
        integral = integral + intprojected(x)
#        print('integral', integral)
        if len(integral2) > 0:
            for y in integral2:
                y.exec_delta()
#        print('integralw', integral)
        # if len(intprojected) > 0:
        #     print('intprojected(x)', intprojected(x))

    return integral  
 

def is_delta_zero(delta):
    """
    Returns True if delta is 0, and False 
    if delta is != 0
    """
    if (delta[0] in occupied) and (delta[1] in virtualall):
        return True
    elif (delta[0] in virtualall) and (delta[1] in occupied):
        return True
    else:
        return False


def filterout(l, retained):
    f = []
    for p in l:
        if p in retained:
            f.append(p)

    return f

def free_fixed():
    fixed = []


# def join(this, other):

#     res = arithmetic_string()

#     for x in this:
#         for y in other:
#             result = ugg()
#             result.coefficient = x.coefficient + y.coefficient
#             result.summation = x.summation
#             for z in y.summation:
#                 if z not in x.summation:
#                     result.summation.append(z)
#             result.num_factor = x.num_factor * y.num_factor
#             result.coefficient_idx = x.coefficient_idx + y.coefficient_idx
# #            result.delta = x.delta + y.delta
#             for z in y.delta:
#                 if z not in x.delta:
#                     result.delta.append(z)
#             result.operator_idx = x.operator_idx + y.operator_idx
#             res.append(result)
#     return res

def add_coef(ars0, coef):

    ars = deepcopy(ars0)

    for x in range(0, len(ars)):
        ars[x].coefficient.append(coef.coefficient[0])
        ars[x].coefficient_idx.append(coef.coefficient_idx[0])



def merge_temp_into_w(rsimp):
    rsimp2 = deepcopy(rsimp)
    for x in range(0,len(rsimp)):

        rsimp2[x].coefficient = []
        rsimp2[x].coefficient_idx = []
        for i in range (0, len(rsimp[x].coefficient)):
            if rsimp[x].coefficient[i] == TEMP1:
                a = rsimp[x].coefficient_idx[i][0]
                c = rsimp[x].coefficient_idx[i][1]
            elif rsimp[x].coefficient[i] == TEMP2:
                b = rsimp[x].coefficient_idx[i][0]
                d = rsimp[x].coefficient_idx[i][1]
            else:
                rsimp2[x].coefficient.append(rsimp[x].coefficient[i])
                rsimp2[x].coefficient_idx.append(rsimp[x].coefficient_idx[i])
        rsimp2[x].coefficient.append('w')
        rsimp2[x].coefficient_idx.append([a, b, c, d])

    return rsimp2

def extract_coef(coef):
                                    
    coef_name =""
    no = "0123456789"
    for s in coef:
        if s == "_":
            break
        else:
            coef_name += s

    return coef_name

def nswaps_in_sort(nums2):

    #compute_n_of_swaps_in_sort

    nums = []
    for x in nums2:
        nums.append(str(x))
    
#    print('nums', nums)
    sort_seq = sorted(nums)
#    print('sort_seq', sort_seq)
    table = {}

    for i in range(0, len(nums)):
        table[nums[i]] = i
        
#    print('table', table)
    swaps = 0
    for i in range(len(nums)):
        n = nums[i]
        s_n = sort_seq[i]
        s_i = table[s_n]
#        print(n, s_n, s_i)

        if s_n != n:
#            print('plusz 1')
            swaps += 1
            nums[s_i] = n
            nums[i] = s_n
            table[n] = s_i
            table[s_n] = i
#            print('nums', nums)
#            print('table', table)

#        print('')
    return swaps


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

def vspace(size):

    if size == 0:
        print('')
        print('')
    elif size ==1:
        print('')
        print('')
        print('')
    elif size == 2:
        print('')
        print('')
        print('')
        print('')
        print('')
