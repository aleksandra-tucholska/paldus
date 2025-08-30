#
# Automated coding of loops over elements of
# the EOM-CC Jacobian matrix
# ------------------------------------------------------
# Authors: Marcin Modrzejewski, University of Warsaw
#          Aleksandra Tucholska, University of Warsaw
# ------------------------------------------------------
#
# Occupied index changes faster than the virtual one.
# Whenever permutation symmetry can be utilized, lower
# triangle matrix elements are computed.
#
# Single electron index: ai
# --------------------------
# NOCC0, NVIRT0 - First occupied (virtual) index
# NOCC1, NVIRT1 - Last occupied (virtual) index
# NOCC = NOCC1 - NOCC0 + 1
# NVIRT = NVIRT1 - NVIRT0 + 1
#
# ai = NOCC * (a - NVIRT0) + (i - NOCC0) + 1
#
# Two electron index: aibj
# --------------------------
# In case of <aibj| X> matrix element, only those
# elements for which ai >= bj holds are computed.
#
# NPAIR - Number of distinct one electron indices
# NPAIR = NOCC * NVIRT
# 
# aibj = [(2 * NPAIR - bj + 2) * (bj - 1)] / 2 + ai - bj + 1
#
# Three electron index: aibjck
# --------------------------
# In case of <aibjck| X> matrix element, only those
# elements for which ai >= bj >= ck holds are computed.
#
# NPAIR - Number of distinct one electron indices
# NPAIR = NOCC * NVIRT
# 
# aibjck = (6 * ai - 3 * bj**2 + ck**3 - 3 * ck**2 (1 + NPAIR) &
#          -3 * NPAIR * (3 + NPAIR) + bj * (3 + 6 * NPAIR) &
#          + ck * (2 + 3 * NPAIR * (2 + NPAIR))) / 6
# 
# 
from copy import deepcopy
from params import occupied
from params import virtual
import datetime
import sys
from collections import Counter
#
# Constants defining types of pending operations
#
PEND_COMPIDX = 1
PEND_LOBOUND = 2
PEND_UPBOUND = 3
PEND_INTDECL = 4
PEND_CYCLTRI = 5
PEND_EXITTRI = 6
#
# Maximum number of loops per single Fortran file
#
MAX_LOOPS = 50

def listsep(l, separator):
    if len(l) == 0:
        s = ""
    else:
        len_l = len(l)
        s = ""
        for k in range(0, len_l):
            if len(s) > 0:
                s += (separator + l[k])
            else:
                s += l[k]
    return s


def commasep(l):
    """
    Convert list of indices ["i", "j", "k"] into
    comma separated string: "i, j, k".
    """
    return listsep(l, ", ")


def newlinesep(l):
    """
    Convert list of strings ["AA", "BB", "CC"] into
    "AA\nBB\nCC".
    """
    return listsep(l, "\n")


def pendop_add(key, code, pending_dict):
    """
    Add new pending operations.
    """
    if key in pending_dict:
        pending_dict[key].append(code)
    else:
        l = [code]
        pending_dict[key] = l


def pendop_count(key, pending_dict):
    """
    Count pending operations of specified type.
    """
    n = 0
    if key in pending_dict:
        n = len(pending_dict[key])

    return n


def pendop_addskiptriplecond(i, j, k, specname, pending_dict):
    """
    Add code that skips triples if occupied indices are 
    ordered in decreasing order: E^{abc}_{ijk} (i > j > k).
    It is assumed that i != j != k. I, J, K indices passed
    as arguments to this subroutine should be substituted
    symbols (free indices).
    """
    l = [i, j, k]
    l.sort()
    #
    # IDX0 is the index in the most deeply 
    # placed loop. 
    #
    idx0 = l[0]
    loopname = "{0}_{1}".format(idx0, specname)
    d = {"i":i, "j":j, "k":k, "loopname":loopname}
    code = -1
    if i == l[0]:
        d["cmd"] = "exit"
        code = PEND_EXITTRI
    else:
        d["cmd"] = "cycle"
        code = PEND_CYCLTRI

    s = "if ({i} > {j} .and. {j} > {k}) {cmd} {loopname}".format(**d)
    pendop_add((idx0, code), s, pending_dict)


def pendop_addinequality(ineq, pending_dict):
    """
    Process inequality to get loop bounds. INEQ represents
    inequality as a list of indices:
    INEQ = ["a", "b"] # a >= b
    A and B are free indices, or indices after substitution
    """
    a = ineq[0]
    b = ineq[1]

    if a != b:
        #
        # Decide when the inequality will be applied.
        # The condition is that the values of both indices
        # must be known
        #
        if a < b:
            #
            # If A < B, then the loop over A will be deeper
            # in the nested loop structure than loop over B.
            # The value of B will be known when the loop over
            # A will be executed.
            # ------------------------------------------------
            # Apply the inequality as a lower bound of the
            # loop over A. Inequality is added only if it
            # is nonredundant
            #
            n = pendop_count((a, PEND_LOBOUND), pending_dict)
            if n == 0:
                pendop_add((a, PEND_LOBOUND), b, pending_dict)
            else:
                if b not in pending_dict[(a, PEND_LOBOUND)]:
                    pendop_add((a, PEND_LOBOUND), b, pending_dict)
        else:
            #
            # Apply the inequality as a upper bound of the
            # loop over B. Inequality is added only if it
            # is nonredundant
            #
            n = pendop_count((b, PEND_UPBOUND), pending_dict)
            if n == 0:
                pendop_add((b, PEND_UPBOUND), a, pending_dict)
            else:
                if a not in pending_dict[(b, PEND_UPBOUND)]:
                    pendop_add((b, PEND_UPBOUND), a, pending_dict)
    

def maxval_code(l, a):
    args = commasep(l)
    d = {"a":a, "args":args}
    s = "{a} = max({args})".format(**d)
    return s


def minval_code(l, a):
    args = commasep(l)
    d = {"a":a, "args":args}
    s = "{a} = min({args})".format(**d)
    return s


def isfree(a, subst):
    if subst[a] == a:
        return True
    else:
        return False


def nondiag_cond(a, skip_nondiag_list, looplabel, braket_type, subst):
    n = braket_type[0] + braket_type[1]
    l = []
    if a in occupied:
        pos = occupied.index(a)
        for k in range(pos + 1, n):
            b = occupied[k]
            if b not in skip_nondiag_list:
                if isfree(b, subst):
                    l.append(a + " == " + b)
    elif a in virtual:
        pos = virtual.index(a)
        for k in range(pos + 1, n):
            b = virtual[k]
            if b not in skip_nondiag_list:
                if isfree(b, subst):
                    l.append(a + " == " + b)

    if len(l) > 0:
        equalities = listsep(l, " .or. ")
        f = {"equalities":equalities, "looplabel":looplabel}
        s = "if ({equalities}) cycle {looplabel}".format(**f)
    else:
        s = ""
    
    return s


def braidx(braket_type):
    bratype = braket_type[0]
    if bratype == 1:
        s = "ibra = braoffset + ai"
    elif bratype == 2:
        s = """ibra = braoffset + &
((2 * npair - bj + 2) * (bj - 1)) / 2 + ai - bj + 1"""
    elif bratype == 3:
        s = "ibra = braoffset + mu3(ai, bj, ck)"
    return s


def ketidx(braket_type):
    kettype = braket_type[1]
    ketidx1 = braket_type[0]
    ketidx2 = braket_type[0] + 1
    ketidx3 = braket_type[0] + 2
    if kettype == 1:
        ck = virtual[ketidx1] + occupied[ketidx1]
        f = {"ck":ck}
        s = "iket = ketoffset + {ck}".format(**f)
    elif kettype == 2:
        ck = virtual[ketidx1] + occupied[ketidx1]
        dl = virtual[ketidx2] + occupied[ketidx2]
        f = {"ck":ck, "dl":dl}
        s = """iket = ketoffset + &
((2 * npair - {dl} + 2) * ({dl} - 1)) / 2 + {ck} - {dl} + 1""".format(**f)
    elif kettype == 3:
        ck = virtual[ketidx1] + occupied[ketidx1]
        dl = virtual[ketidx2] + occupied[ketidx2]
        em = virtual[ketidx3] + occupied[ketidx3]
        f = {"ck":ck, "dl":dl, "em":em}
        s = "iket = ketoffset + mu3({ck}, {dl}, {em})".format(**f)
    return s


# def permute_arglist(arglist, ijk):
#     arglist_permuted = []

#     for x in arglist:
#         if x == occupied[0]:
#             arglist_permuted.append(ijk[0])
#         elif x == occupied[1]:
#             arglist_permuted.append(ijk[1])
#         elif x == occupied[2]:
#             arglist_permuted.append(ijk[2])
#         else:
#             arglist_permuted.append(x)

#     return arglist_permuted


def whichvk(braket_type, kronecker):
    #
    # WHICHVK value
    # --------------
    # -1 Redundant bra or ket vector
    #    OR
    #    Single index repeated three times in bra or ket
    #  0 Bra vector: v_0
    #  1 Bra vectors: v_1 ... v_5
    #  6 Bra vector: v6
    #
    substitutions = {}
    occ_summ_idx = []
    virt_summ_idx = []
    substitutions, occ_summ_idx, virt_summ_idx = comp_substitutions_dict(braket_type, kronecker)
    
    if braket_type[0] == 3 or braket_type[1] == 3:
        isthreefold = is_threefold(braket_type, substitutions)
        if isthreefold:
            return -1

        isredundant = isredundant_iik(braket_type, substitutions)
        if isredundant:
            return -1
    
    if braket_type[0] == 3:
        return tripletype(substitutions)


def tripletype(substitutions):
    bra_occupied = occupied[0:3]
    bra_virtual = virtual[0:3]

    s_bra_occupied = []
    s_bra_virtual = []
    for x in bra_occupied:
        s_bra_occupied.append(substitutions[x])
    for x in bra_virtual:
        s_bra_virtual.append(substitutions[x])

    cv = Counter(s_bra_virtual)
    co = Counter(s_bra_occupied)
    
    vk = -2

    if (cv.most_common(1)[0][1] == 1) and (co.most_common(1)[0][1] == 1):
        #
        # CASE 1
        #
        vk = 1
    elif (cv.most_common(1)[0][1] == 1) and (co.most_common(1)[0][1] == 2):
        #
        # CASE 2
        #
        vk = 0
    elif (cv.most_common(1)[0][1] == 2) and (co.most_common(1)[0][1] == 1):
        #
        # CASE 3
        # or
        # CASE 4
        #
        vk = 6
    elif (cv.most_common(1)[0][1] == 2) and (co.most_common(1)[0][1] == 2):
            if (s_bra_virtual[0] == s_bra_virtual[1]) and (s_bra_occupied[0] == s_bra_occupied[2]):
                #
                # CASE 5
                #
                vk = 6
            elif (s_bra_virtual[0] == s_bra_virtual[1]) and (s_bra_occupied[1] == s_bra_occupied[2]):
                #
                # CASE 6
                #
                vk = 0
            elif (s_bra_virtual[1] == s_bra_virtual[2]) and (s_bra_occupied[0] == s_bra_occupied[2]):
                #
                # CASE 7
                #
                vk = 6
            elif (s_bra_virtual[1] == s_bra_virtual[2]) and (s_bra_occupied[0] == s_bra_occupied[1]):
                #
                # CASE 8
                #
                vk = 0

    return vk


def matrix_element_code(braket_type, substitutions, d, vanish_case1):
    #
    # Determine arguments (free indices)
    #
    arglist = []
    for k in range(braket_type[0]+braket_type[1]):
        a = virtual[k]
        i = occupied[k]

        if isfree(a, substitutions):
            arglist.append(a)
        if isfree(i, substitutions):
            arglist.append(i)

    s_arglist = commasep(arglist)
    f = deepcopy(d)
    f["args"] = s_arglist

    if braket_type[0] != 3:
        #
        # Single or double excitation in bra
        #
        s = "jac(ibra, iket) = {genname}_{specname}(eorb, t2, t1, nocc, nactive, {args})".format(**f)
    else:
        #
        # Triples manifold
        # ---
        # Determine class of triple excitation
        #
        vk = tripletype(substitutions)

        if vk != 1:
            f["vk"] = vk
            s = "jac(ibra, iket) = v{vk}_{genname}_{specname}(eorb, t2, t1, nocc, nactive, {args})".format(**f)
        else:
            # #
            # # Bra vector belonging to <E^{abc}_{ijk}| class (no repeated indices)
            # # ---
            # # Permuted arguments (see v1...v5 vectors)
            # #
            # arglist_kij = permute_arglist(arglist, "kij")
            # arglist_ikj = permute_arglist(arglist, "ikj")
            # arglist_jik = permute_arglist(arglist, "jik")
            # arglist_kji = permute_arglist(arglist, "kji")
            # arglist_jki = permute_arglist(arglist, "jki")

            # s_arglist_kij = commasep(arglist_kij)
            # s_arglist_ikj = commasep(arglist_ikj)
            # s_arglist_jik = commasep(arglist_jik)
            # s_arglist_kji = commasep(arglist_kji)
            # s_arglist_jki = commasep(arglist_jki)

            # f["args_kij"] = s_arglist_kij
            # f["args_ikj"] = s_arglist_ikj
            # f["args_jik"] = s_arglist_jik
            # f["args_kji"] = s_arglist_kji
            # f["args_jki"] = s_arglist_jki

            a = virtual[0]
            b = virtual[1]
            c = virtual[2]
            i = occupied[0]
            j = occupied[1]
            k = occupied[2]

            f.update({"a":a, "b":b, "c":c, "i":i, "j":j, "k":k})
            v1_funcname = "v1_{genname}_{specname}".format(**f)
            v2_funcname = "v2_{genname}_{specname}".format(**f)
            v3_funcname = "v3_{genname}_{specname}".format(**f)
            v4_funcname = "v4_{genname}_{specname}".format(**f)
            v5_funcname = "v5_{genname}_{specname}".format(**f)
            #
            # ----------------------------------------------------------------
            # Indices in decreasing order       Label of bra vector V_k label
            # according to numerical value      (see PDF)           (see PDF)
            # ----------------------------------------------------------------
            # 1. ijk                            Not used
            # 2. kij                            jki                 v_3
            # 3. ikj                            ikj                 v_1
            # 4. jik                            jik                 v_2
            # 5. kji                            kji                 v_5
            # 6. jki                            kij                 v_4
            # ----------------------------------------------------------------
            #
            if v1_funcname not in vanish_case1:
                v1_call = "v1_{genname}_{specname}(eorb, t2, t1, nocc, nactive, {args})".format(**f)
            else:
                v1_call = "0.d+0"

            if v2_funcname not in vanish_case1:
                v2_call = "v2_{genname}_{specname}(eorb, t2, t1, nocc, nactive, {args})".format(**f)
            else:
                v2_call = "0.d+0"

            if v3_funcname not in vanish_case1:
                v3_call = "v3_{genname}_{specname}(eorb, t2, t1, nocc, nactive, {args})".format(**f)
            else:
                v3_call = "0.d+0"

            if v4_funcname not in vanish_case1:
                v4_call = "v4_{genname}_{specname}(eorb, t2, t1, nocc, nactive, {args})".format(**f)
            else:
                v4_call = "0.d+0"

            if v5_funcname not in vanish_case1:
                v5_call = "v5_{genname}_{specname}(eorb, t2, t1, nocc, nactive, {args})".format(**f)
            else:
                v5_call = "0.d+0"

            f.update({"v1_call":v1_call, "v2_call":v2_call, "v3_call":v3_call, "v4_call":v4_call, "v5_call":v5_call})

            s = """
if ({i} > {j}) then
      if ({j} > {k}) then
            exit {i}_{specname}
      else if ({k} > {i}) then
            jac(ibra, iket) = {v3_call}
      else
            jac(ibra, iket) = {v1_call}
      end if
else if ({i} > {k}) then
      jac(ibra, iket) = {v2_call}
else if ({k} > {j}) then
      jac(ibra, iket) = {v5_call}
else
      jac(ibra, iket) = {v4_call}
end if
""".format(**f)

    return s


def jacobian_loops(genname, specname, level, occ_idx, 
                   virt_idx, occ_summ_idx, virt_summ_idx, substitutions, 
                   braket_type, pending, vanish_case1):
    #
    # Test if current occ-virt pair of indices belongs
    # to bra or ket subset
    #
    nbra = braket_type[0]
    nket = braket_type[1]
    isbra = level > nket
    isket = level <= nket

    maxlevel = sum(braket_type)
    occ_this = occ_idx[level - 1] 
    virt_this = virt_idx[level - 1]
    sum_i = occ_this in occ_summ_idx
    sum_a = virt_this in virt_summ_idx
    if sum_i:
        occ_this_subst = occ_this
    else:
        occ_this_subst = substitutions[occ_this]

    if sum_a:
        virt_this_subst = virt_this
    else:
        virt_this_subst = substitutions[virt_this]
    #
    # Generalized indices of the Jacobian
    #
    ibra = braidx(braket_type)
    iket = ketidx(braket_type)
    d = {"i0":"nocc0", "i1":"nocc1", "a0":"nvirt0", "a1":"nvirt1", \
             "occ_this":occ_this, "virt_this":virt_this, \
             "occ_this_subst":occ_this_subst, "virt_this_subst":virt_this_subst, \
             "specname":specname, "genname":genname, \
             "ibra":ibra, "iket":iket}
    #        
    # Determine lower and upper bounds for the loop
    # over occupied index
    #
    occ_pending_lobound = ""
    i0 = ""
    nlo = pendop_count((occ_this, PEND_LOBOUND), pending)
    #
    # OCC_SKIP_NONDIAG_LIST is a list of indices that will
    # be excluded from diagonality test. In those cases
    # non-diagonality is already ensured by the form of
    # upper/lower bounds
    #
    occ_skip_nondiag_list = []
    if nlo == 0:
        #
        # No constraints on the lower bound found.
        # Use generic bound
        #
        i0 = "nocc0"
    elif nlo == 1:
        #
        # A single inequality found. The inequality I >= X
        # together with non-diagonality condition I != X
        # is transformed into I >= X + 1 to decrease
        # the number of non-diagonality tests.
        #
        i0 = pending[(occ_this, PEND_LOBOUND)][0] + " + 1"
        #
        # Add index excluded from diagonality test
        #
        occ_skip_nondiag_list += pending[(occ_this, PEND_LOBOUND)]
    else:
        #
        # Multiple inequalities found. MAXVAL_CODE will be used
        # to find the highest lower bound. This bound will be
        # denoted as I0. "I0 + 1" lower bound is used due to
        # non-diagonality conditions
        #
        i0 = occ_this + "0" + " + 1"
        #
        # Compute the highest lower bound. This code is executed
        # before entering the loop
        #
        occ_pending_lobound = maxval_code(pending[(occ_this, PEND_LOBOUND)], occ_this+"0")
        #
        # Declare new integer variable if necessary
        #
        pendop_add(PEND_INTDECL, occ_this+"0", pending)
        #
        # Add multiple indieces to excluded list
        #
        occ_skip_nondiag_list += pending[(occ_this, PEND_LOBOUND)]

    d["i0"] = i0
    d["occ_pending_lobound"] = occ_pending_lobound

    occ_pending_upbound = ""
    i1 = ""
    nup = pendop_count((occ_this, PEND_UPBOUND), pending)
    if nup == 0:
        i1 = "nocc1"
    elif nup == 1:
        i1 = pending[(occ_this, PEND_UPBOUND)][0] + " - 1"
        occ_skip_nondiag_list += pending[(occ_this, PEND_UPBOUND)]
    else:
        i1 = occ_this + "1" + " - 1"
        occ_pending_upbound = minval_code(pending[(occ_this, PEND_UPBOUND)], occ_this+"1")
        pendop_add(PEND_INTDECL, occ_this+"1", pending)
        occ_skip_nondiag_list += pending[(occ_this, PEND_UPBOUND)]

    d["i1"] = i1
    d["occ_pending_upbound"] = occ_pending_upbound
    #
    # Determine lower and upper bounds for the loop over
    # virtual index
    #
    virt_pending_lobound = ""
    a0 = ""
    nlo = pendop_count((virt_this, PEND_LOBOUND), pending)
    virt_skip_nondiag_list = []
    if nlo == 0:
        a0 = "nvirt0"
    elif nlo == 1:
        a0 = pending[(virt_this, PEND_LOBOUND)][0] + " + 1"
        virt_skip_nondiag_list += pending[(virt_this, PEND_LOBOUND)]
    else:
        a0 = virt_this + "0" + " + 1"
        virt_pending_lobound = maxval_code(pending[(virt_this, PEND_LOBOUND)], virt_this+"0")
        pendop_add(PEND_INTDECL, virt_this+"0", pending)
        virt_skip_nondiag_list += pending[(virt_this, PEND_LOBOUND)]

    d["a0"] = a0
    d["virt_pending_lobound"] = virt_pending_lobound

    virt_pending_upbound = ""
    a1 = ""
    nup = pendop_count((virt_this, PEND_UPBOUND), pending)
    if nup == 0:
        a1 = "nvirt1"
    elif nup == 1:
        a1 = pending[(virt_this, PEND_UPBOUND)][0] + " - 1"
        virt_skip_nondiag_list += pending[(virt_this, PEND_UPBOUND)]
    else:
        a1 = virt_this + "1" + " - 1"
        virt_pending_upbound = minval_code(pending[(virt_this, PEND_UPBOUND)], virt_this+"1")
        pendop_add(PEND_INTDECL, virt_this+"1", pending)
        virt_skip_nondiag_list += pending[(virt_this, PEND_UPBOUND)]

    d["a1"] = a1
    d["virt_pending_upbound"] = virt_pending_upbound
    #
    # Determine pending operations for the current pair of indices
    #
    if pendop_count((occ_this, PEND_COMPIDX), pending) > 0:
        occ_pending_compidx = newlinesep(pending[(occ_this, PEND_COMPIDX)])
    else:
        occ_pending_compidx = ""

    if pendop_count((virt_this, PEND_COMPIDX), pending) > 0:
        virt_pending_compidx = newlinesep(pending[(virt_this, PEND_COMPIDX)])
    else:
        virt_pending_compidx = ""
    d["occ_pending_compidx"] = occ_pending_compidx
    d["virt_pending_compidx"] = virt_pending_compidx
    #
    # Determine whether i > j > k condition should
    # be applied. This is relevant if triples are present,
    # and linearly independent basis is requested.
    #
    if pendop_count((occ_this, PEND_EXITTRI), pending) > 0:
        occ_pending_skiptriple_exit = newlinesep(pending[(occ_this, PEND_EXITTRI)])
    else:
        occ_pending_skiptriple_exit = ""
    if pendop_count((occ_this, PEND_CYCLTRI), pending) > 0:
        occ_pending_skiptriple_cycle = newlinesep(pending[(occ_this, PEND_CYCLTRI)])
    else:
        occ_pending_skiptriple_cycle = ""
    d["occ_pending_skiptriple_exit"] = occ_pending_skiptriple_exit
    d["occ_pending_skiptriple_cycle"] = occ_pending_skiptriple_cycle
    #
    # Loop labels
    #
    occ_looplabel = "{occ_this}_{specname}".format(**d)
    virt_looplabel = "{virt_this}_{specname}".format(**d)
    d["occ_looplabel"] = occ_looplabel
    d["virt_looplabel"] = virt_looplabel
    #
    # Non-diagonality conditions
    #
    occ_nondiag = nondiag_cond(occ_this, occ_skip_nondiag_list, occ_looplabel, braket_type, substitutions)
    virt_nondiag = nondiag_cond(virt_this, virt_skip_nondiag_list, virt_looplabel, braket_type, substitutions)
    d["occ_nondiag"] = occ_nondiag
    d["virt_nondiag"] = virt_nondiag

    if level < maxlevel:
        innercode = jacobian_loops(genname, specname, level + 1, occ_idx, 
                                   virt_idx, occ_summ_idx, virt_summ_idx, substitutions, 
                                   braket_type, pending, vanish_case1)
    else:
        comp_matrix_element =  matrix_element_code(braket_type, substitutions, d, vanish_case1)
        d["comp_matrix_element"] = comp_matrix_element

        innercode = """
{ibra}
{iket}
{comp_matrix_element}
""".format(**d)
    
    d["innercode"] = innercode

    if sum_i and sum_a:
        s = """
{virt_pending_lobound}
{virt_pending_upbound}
{virt_looplabel}: do {virt_this} = {a0}, {a1}
{virt_nondiag}
{virt_pending_compidx}
{occ_pending_lobound}
{occ_pending_upbound}
{occ_looplabel}: do {occ_this} = {i0}, {i1}
{occ_nondiag}
{occ_pending_skiptriple_exit}
{occ_pending_skiptriple_cycle}
{occ_pending_compidx}
{innercode}
end do {occ_looplabel}
end do {virt_looplabel}
""".format(**d)

    elif sum_i and (not sum_a):
        s = """
{occ_pending_lobound}
{occ_pending_upbound}
{occ_looplabel}: do {occ_this} = {i0}, {i1}
{occ_nondiag}
{occ_pending_skiptriple_exit}
{occ_pending_skiptriple_cycle}
{occ_pending_compidx}
{innercode}
end do {occ_looplabel}
""".format(**d)

    elif (not sum_i) and sum_a:
        s = """
{virt_pending_lobound}
{virt_pending_upbound}
{virt_looplabel}: do {virt_this} = {a0}, {a1}
{virt_nondiag}
{virt_pending_compidx}
{innercode}
end do {virt_looplabel}
""".format(**d)

    elif (not sum_i) and (not sum_a):
        s = """
{innercode}
""".format(**d)
        
    return s


def skip_triple_inequality(virtidx, occidx, substitutions, specname, pending_dict):
    s_virt = []
    s_occ = []
        
    for x in virtidx:
        s_virt.append(substitutions[x])
        
    for x in occidx:
        s_occ.append(substitutions[x])

    c_virt = Counter(s_virt)
    c_occ = Counter(s_occ)

    if c_occ.most_common(1)[0][1] == 1:
        if c_virt.most_common(1)[0][1] == 1:
            #
            # CASE 1.
            # ---
            # No two identical virtual indices.
            # No two identical occupied indices:
            # ---
            # If no additional constraints were applied,
            # there would be six vectors of type
            # E^{abc}_{ijk} (permute i, j, k).
            # The vector with occupied indices
            # in decreasing order is skipped so that
            # only linearly independent triples are present.
            #
            iprime = s_occ[0]
            jprime = s_occ[1]
            kprime = s_occ[2]
            #
            # Code that exits loop if i' > j' > k'
            #
            pendop_addskiptriplecond(iprime, jprime, kprime, specname, pending_dict)
                
        elif c_virt.most_common(1)[0][1] == 2:
            #
            # Two identical virtual indices: E^{aac}_{ijk} or E^{abb}_{ijk} 
            # No two identical occupied indices
            # ---
            # The following triples would be present if no
            # additional constraints would be applied:
            # CASE 3. E^{aac}_{ijk}, E^{aac}_{ikj}, E^{aac}_{jki}
            # CASE 4. E^{abb}_{ijk}, E^{abb}_{jik}, E^{abb}_{kij}
            # (i > j > k assumed above)
            # In each group, the vector with occupied indices
            # in descending order (E^{aac}_{ijk} and E^abb}_{ijk})
            # is skipped so that only linearly independent triples
            # are present. Let i', j', and k' denote first, second,
            # and third occupied index from the left. If we enforce
            # 1. k' > j' (Group 1)
            # 2. j' > i' (Group 2)
            # then we exclude the E^{aac}_{ijk} and E^{baa}_{ijk}
            # vectors.
            #
            if s_virt[0] == s_virt[1]:
                #
                # E^{aac}_{ijk}
                #
                kprime = s_occ[2]
                jprime = s_occ[1]
                #
                # Apply k' > j' inequality
                #
                pendop_addinequality([kprime, jprime], pending_dict)
            elif s_virt[1] == s_virt[2]:
                #
                # E^{abb}_{ijk} 
                #
                jprime = s_occ[1]
                iprime = s_occ[0]
                #
                # Apply j' > i'
                # 
                pendop_addinequality([jprime, iprime], pending_dict)   


def nonredundant_ijk(braket_type, substitutions, specname, pending_dict):
    #
    # Add inequalities to eliminate redundant triples
    # (both bra and ket vectors) in case when three
    # different occupied indices are present:
    # CASE 3. E^{aac}_{ijk} + ijk permutations
    # CASE 4. E^{abb}_{ijk} + ijk permutations
    # CASE 1. E^{abc}_{ijk} + ijk permutations
    #
    bra_virtual = virtual[0:braket_type[0]]
    bra_occupied = occupied[0:braket_type[0]]
    ket_virtual = virtual[braket_type[0]:braket_type[0]+braket_type[1]]
    ket_occupied = occupied[braket_type[0]:braket_type[0]+braket_type[1]]

    if braket_type[1] == 3:
        skip_triple_inequality(ket_virtual, ket_occupied, substitutions, specname, pending_dict)

    if braket_type[0] == 3:
        skip_triple_inequality(bra_virtual, bra_occupied, substitutions, specname, pending_dict)



def loopbounds(braket_type, substitutions, pending_dict):
    """
    Determine upper and lower bounds in every loop over
    a free one-orbital index. This subroutine works for
    arbitrary braket type, including triples,
    quadruples, etc.
    """
    #
    # First index of ket vector
    #
    firstket = braket_type[0]
    #
    # Number of compound (occ-virt pair) indices
    #
    n = braket_type[0] + braket_type[1]
    for k in range(0, n - 1):
        if (k + 1) == firstket:
            #
            # Cycle the loop because there is no
            # inequality between bra and ket indices
            #
            continue

        a = virtual[k]
        b = virtual[k + 1]
        subst_a = substitutions[a]
        subst_b = substitutions[b]
        if subst_a != subst_b:
            #
            # Add to pending operations code that enforces
            # a >= b
            #
            pendop_addinequality([subst_a, subst_b], pending_dict)
        else:
            #
            # If virtual orbital indices are identical,
            # apply i >= j. Note that in ai pair, the occupied
            # index changes faster
            #
            i = occupied[k]
            j = occupied[k + 1]
            subst_i = substitutions[i]
            subst_j = substitutions[j]
            if subst_i != subst_j:
                pendop_addinequality([subst_i, subst_j], pending_dict)


def isredundant_iik(braket_type, substitutions):
    #
    # Check if loop over a given class of matrix
    # elements should be skipped due to the fact
    # that bra/ket vectors are not present in the 
    # the independent triples' basis. The following
    # cases are considered:
    #
    # CASE 2. EaiEbiEck, EaiEbkEci, EakEbiEci (Exclude EaiEbiEck)
    # CASE 5. EaiEaiEck, EaiEakEci (Exclude EaiEaiEck)
    # CASE 6. EaiEajEcj, EajEajEci (Exclude EajEajEci)
    # CASE 7. EaiEbjEbj, EajEbiEbj (Exclude EaiEbjEbj)
    # CASE 8. EaiEbiEbj, EajEbiEbi (Exclude EajEbiEbi)
    # 
    # For more details see description of the nonredundant
    # basis in a separete PDF file.
    #
    #
    isredundant = False

    bra_occupied = occupied[0:braket_type[0]]
    bra_virtual = virtual[0:braket_type[0]]
    ket_occupied = occupied[braket_type[0]:braket_type[0]+braket_type[1]]
    ket_virtual = virtual[braket_type[0]:braket_type[0]+braket_type[1]]

    if braket_type[0] == 3:
        s_bra_occupied = []
        s_bra_virtual = []
        for x in bra_occupied:
            s_bra_occupied.append(substitutions[x])
        for x in bra_virtual:
            s_bra_virtual.append(substitutions[x])

        cv = Counter(s_bra_virtual)
        co = Counter(s_bra_occupied)
        if (cv.most_common(1)[0][1] == 1) and (co.most_common(1)[0][1] == 2):
            if s_bra_occupied[0] == s_bra_occupied[1]:
                #
                # CASE 2. E^{abc}_{iik}
                #
                isredundant = True
        elif (cv.most_common(1)[0][1] == 2) and (co.most_common(1)[0][1] == 2):
            if (s_bra_virtual[0] == s_bra_virtual[1]) and (s_bra_occupied[0] == s_bra_occupied[1]):
                #
                # CASE 5. E^{aac}_{iik}
                # or
                # CASE 6. E^{aac}_{jji}
                #
                isredundant = True
            elif (s_bra_virtual[1] == s_bra_virtual[2]) and (s_bra_occupied[1] == s_bra_occupied[2]):
                #
                # CASE 7. E^{abb}_{ijj}
                # or
                # CASE 8. E^{abb}_{jii}
                #
                isredundant = True

    if braket_type[1] == 3:
        s_ket_occupied = []
        s_ket_virtual = []
        for x in ket_occupied:
            s_ket_occupied.append(substitutions[x])
        for x in ket_virtual:
            s_ket_virtual.append(substitutions[x])

        cv = Counter(s_ket_virtual)
        co = Counter(s_ket_occupied)
        if (cv.most_common(1)[0][1] == 1) and (co.most_common(1)[0][1] == 2):
            if s_ket_occupied[0] == s_ket_occupied[1]:
                #
                # CASE 2. E^{abc}_{iik}
                #
                isredundant = True
        elif (cv.most_common(1)[0][1] == 2) and (co.most_common(1)[0][1] == 2):
            if (s_ket_virtual[0] == s_ket_virtual[1]) and (s_ket_occupied[0] == s_ket_occupied[1]):
                #
                # CASE 5. E^{aac}_{iik}
                # or
                # CASE 6. E^{aac}_{jji}
                #
                isredundant = True
            elif (s_ket_virtual[1] == s_ket_virtual[2]) and (s_ket_occupied[1] == s_ket_occupied[2]):
                #
                # CASE 7. E^{abb}_{ijj}
                # or
                # CASE 8. E^{abb}_{jii}
                #
                isredundant = True

    return isredundant


def pairindices(braket_type, subst, pending_dict):
    n = braket_type[0] + braket_type[1]
    for k in range(0, n):
        a = virtual[k]
        i = occupied[k]
        ai = a + i
        subst_a = subst[a]
        subst_i = subst[i]
        idxa = virtual.index(subst_a)
        idxi = occupied.index(subst_i)
        f = {"ai":ai, "a":subst_a, "i":subst_i}
        s = "{ai} = ({a} - nvirt0) * nocc + ({i} - nocc0) + 1".format(**f)
        if idxi <= idxa:
            pendop_add((subst_i, PEND_COMPIDX), s, pending_dict)
        else:
            pendop_add((subst_a, PEND_COMPIDX), s, pending_dict)


def verify_deltalist(deltalist):
    #
    # Search for diagonal deltas, \delta_{pp}
    #
    for l in deltalist:
        for pair in l:
            if pair[0] == pair[1]:
                print("EOMCCJAC ERROR: LIST OF KRONECKER DELTAS CONTAINS DIAGONAL TERMS")
                return False
    
    #
    # Search for deltas that are not sorted in alphabetical order
    #
    for l in deltalist:
        for pair in l:
            if pair[0] > pair[1]:
                print("EOMCCJAC ERROR: LIST OF KRONECKER DELTAS CONTAINS UNSORTED INDICES")
                return False
    #
    # Search for duplicate indices
    #
    for l in deltalist:
        for k in range(0, len(l)):
            q = l[k][1]
            for m in range(0, len(l)):
                if m != k:
                    if q in l[m]:
                        print("EOMCCJAC ERROR: SUBSTITUTED INDEX PRESENT IN MORE THAN ONE KRONECKER DELTA")
                        return False

    return True


def is_threefold(braket_type, substitutions):
    #
    # Search for threefold occurence of the same orbital
    # in bra/ket vector. This class of integrals is skipped
    # because <aibici| X> or <aiajak|X> is zero for any X.
    #
    bra_virtual = virtual[0:braket_type[0]]
    bra_occupied = occupied[0:braket_type[0]]
    ket_virtual = virtual[braket_type[0]:braket_type[0]+braket_type[1]]
    ket_occupied = occupied[braket_type[0]:braket_type[0]+braket_type[1]]

    threefold = False

    if braket_type[0] >= 3:
        s_bra = []
        for x in bra_virtual:
            s_bra.append(substitutions[x])
        for x in bra_occupied:
            s_bra.append(substitutions[x])

        c = Counter(s_bra)
        #
        # Check if any symbol (virtual or occupied index) occures
        # three times or more
        #
        if c.most_common(1)[0][1] >= 3:
            threefold = True
        

    if braket_type[1] >= 3:
        s_ket = []
        for x in ket_virtual:
            s_ket.append(substitutions[x])
        for x in ket_occupied:
            s_ket.append(substitutions[x])

        c = Counter(s_ket)
        #
        # Check if any symbol (virtual or occupied index) occures
        # three times or more
        #
        if c.most_common(1)[0][1] >= 3:
            threefold = True

    return threefold


def comp_substitutions_dict(braket_type, deltas):
    n = braket_type[0] + braket_type[1]
    deltas_indices = set([])
    substitutions = {}
    occ_summ_idx = []
    virt_summ_idx = []
    #
    # Compute dictionary of index substitutions.
    # Summation indices (free indices) are substituted
    # for themselves by definition, so every index has
    # its key in SUBSTITUTIONS dictionary.
    #
    for deltapair in deltas:
        p = deltapair[0]
        q = deltapair[1]
        deltas_indices.add(p)
        deltas_indices.add(q)
        if (p in occupied) and (p not in occ_summ_idx):
            occ_summ_idx.append(p)
        elif (p in virtual) and (p not in virt_summ_idx):
            virt_summ_idx.append(p)
        substitutions[q] = p
        substitutions[p] = p
        
    for o in occupied[0:n]:
        if o not in deltas_indices:
            occ_summ_idx.append(o)
            substitutions[o] = o

    for v in virtual[0:n]:
        if v not in deltas_indices:
            virt_summ_idx.append(v)
            substitutions[v] = v 

    return substitutions, occ_summ_idx, virt_summ_idx


def writefile(s, filename):
    f = open(filename, "w")
    #
    # Skip blank lines
    #
    lines = s.splitlines()
    for line in lines:
        if line != "":
            f.write(line + "\n")
    f.close()


def jacobian_loop(braket_type, nonzero_kronecker, genname, specnames,  vanish_case1=[], independent_triples=True):
    """
    BRAKET_TYPE       - 2-tuple determining type of matrix element:
                        (1, 1) -> <ia| X | jb>
                        (2, 1) -> <iajb| X |kc>
                        (2, 2) -> <iajb| X |kcld>
                        (3, 1) -> <iajbkc| X |ld>
                        (3, 2) -> <iajbkc| X |ldme>
                        (3, 3) -> <iajbkc| X |ldmenf>
                        ...
    NONZERO_KRONECKER - List of lists of nonzero kronecker deltas
                        in canonical form. Example of canonical 
                        list:
                        
                        [["a", "b"], ["a", "c"], ["a", "d"]]

                        instead of

                        [["a", "b"], ["b", "c"], ["c", "d"]].

                        Every possible delta that is not explicitly
                        or implicitly included in
                        NONZERO_KRONECKER[K] is assumed to be ZERO.
                        Determines which elementary loops should be
                        coded. Simplistic example:
                        NONZERO_KRONECKER == [[], [["a", "b"]]]
                        (assume a and b are the only sulmmation
                        indices, ["a", "b"] is a representation 
                        of Kronecker delta). Two elementary loops
                        will be coded:
                        !
                        ! Machine generated Fortran code
                        !
                        do a = 1, N
                              do b = 1, N
                                   if (b == a) cycle
                                   !
                                   ! Execute code dependent on the assumption
                                   ! b /= a: GENNAME_SPECNAMES[0]
                                   !
                              end do
                        end do
    
                        do a = 1, N
                              !
                              ! Execute code dependent on the assumption
                              ! b == a: GENNAME_SPECNAMES[1]
                              !
                        end do
    
    GENNAME           - Generic name of subroutine to call inside each
                        elementary loop. It is assumed that subroutines
                        to call have names of the form GENNAME_SPECNAME
    
    SPECNAMES         - SPECNAME chunk of subroutine name that is called in
                        each elementary loop

    VANISH_CASE1      - (Relevant if INDEPENDENT_TRIPLES==True and BRAKET_TYPE==(3,n))
                        List of functions that correspond to vanishing matrix elements
                        with bra vector belonging to the set: {v_1, v_2, ..., v_5}.
                        If function's name is present in VANISH_CASE1 list,
                        it will not be called. See the PDF with details on the
                        triples basis.
    INDEPENDENT_
    TRIPLES           - Optional argument. If set to True (default value), 
                        linearly independent basis of triples is used. 
                        """
    if (braket_type[0] == 3) and (vanish_case1 == []):
        print("ERROR: VANISH_CASE1 LIST NOT PROVIDED")
        sys.exit()
    #
    # Verify whether NONZERO_KRONECKER list is in canonical form
    #
    iscanonical = verify_deltalist(nonzero_kronecker)
    if not iscanonical:
        print("ERROR: NONSTANDARD LIST OF KRONECKER DELTAS.")
        sys.exit()

    kronecker_list = []
    specnames_list = []
    kronecker_part = []
    specnames_part = []
    
    nloops = 0
    
    for k in range(0, len(nonzero_kronecker)):
        deltas = nonzero_kronecker[k]
        substitutions = {}
        occ_summ_idx = []
        virt_summ_idx = []
        substitutions, occ_summ_idx, virt_summ_idx = comp_substitutions_dict(braket_type, deltas)
        #
        # Check if any index occures three times
        # or more in bra/ket vector. Such matrix
        # elements are identically equal to zero.
        #
        if is_threefold(braket_type, substitutions):
            print("NOTE: {0} CLASS OF MATRIX ELEMENTS IS SKIPPED. THESE ELEMENTS VANISH IDENTICALLY.".format(specnames[k]))
            #
            # Go to the next class of matrix elements
            #
            continue
        #
        # Skip the loop over current class of matrix elements if bra/ket
        # vector belongs to the triple excitations not present in the 
        # independent triples' basis set (occupied index on repeated
        # twice on the left)
        #
        if independent_triples:
            if isredundant_iik(braket_type, substitutions):
                print("NOTE: {0} CLASS OF MATRIX ELEMENTS IS SKIPPED. BRA/KET VECTOR IS NOT PRESENT IN THE INDEPENDENT TRIPLES' SET.".format(specnames[k]))
                #
                # Go to the next class of matrix elements
                #
                continue
        
        nloops += 1
        kronecker_part.append(nonzero_kronecker[k])
        specnames_part.append(specnames[k])
        
        if (nloops == MAX_LOOPS):
            kronecker_list.append(kronecker_part)
            specnames_list.append(specnames_part)
            nloops = 0
            kronecker_part = []
            specnames_part = []
    if (nloops > 0):
        kronecker_list.append(kronecker_part)
        specnames_list.append(specnames_part)

    print("DIVIDING LOOPS INTO {0} FILES.".format(len(kronecker_list)))
    #
    # Number of parts the loops over Jacobian are divided into
    #
    nparts = len(kronecker_list)
    modulelist = []
    calllist = []

    for k in range(0, len(kronecker_list)):
        if nparts > 1:
            module_postfix = "_part{0}".format(k+1)
        else:
            module_postfix = ""
            
        f = {"bra_type":braket_type[0], "ket_type":braket_type[1], "module_postfix":module_postfix}
        sm = "use ccjac_block_{bra_type}{ket_type}{module_postfix}".format(**f)
        ss = "call ccjac_{bra_type}{ket_type}{module_postfix}(jac, eorb, t1, t2, nocc0, nocc1, nvirt0, nvirt1, bra0, ket0)".format(**f)
        modulelist.append(sm)
        calllist.append(ss)
        #
        # Generate Fortran file
        #
        kronecker_part = kronecker_list[k]
        specnames_part = specnames_list[k]
        jacobian_loop_part(braket_type, kronecker_part, genname,
                           module_postfix, specnames_part, independent_triples, vanish_case1)

    if nparts > 1:
        f = {"bra_type":braket_type[0], "ket_type":braket_type[1]}
        f["modulelist"] = newlinesep(modulelist)
        f["calllist"] = newlinesep(calllist)
        utc_datetime = datetime.datetime.utcnow()
        datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")
        f["datestr"] = datestr
        
        s = """
module ccjac_block_{bra_type}{ket_type}
{modulelist}
implicit none
!
! File generated automatically on {datestr} UTC.
!
contains
 
subroutine ccjac_{bra_type}{ket_type}(jac, eorb, t1, t2, nocc0, nocc1, nvirt0, nvirt1, bra0, ket0)
double precision, dimension(:,:), intent(out)       :: jac
double precision, dimension(:), intent(in)          :: eorb
double precision, dimension(:, :), intent(in)       :: t1
double precision, dimension(:, :, :, :), intent(in) :: t2
integer, intent(in)                                 :: nocc0, nocc1
integer, intent(in)                                 :: nvirt0, nvirt1
integer, intent(in)                                 :: bra0, ket0
!
! Generating full CC Jacobian
!
{calllist}
end subroutine ccjac_{bra_type}{ket_type}
end module ccjac_block_{bra_type}{ket_type}
""".format(**f)
        filename = "./fortran_codes/ccjac_bez_dav/ccjac_block_{bra_type}{ket_type}.f90".format(**f)
        writefile(s, filename)
        
    

def jacobian_loop_part(braket_type, nonzero_kronecker, genname,
    module_postfix, specnames, independent_triples, vanish_case1):
    
    bra_type = braket_type[0]
    ket_type = braket_type[1]
    n = sum(braket_type)
    
    occ_idx = occupied[0:n]
    virt_idx = virtual[0:n]
    ovpair_decl = ""
    occ_idx_decl = ""
    virt_idx_decl = ""
    #
    # Strings used in declaration of compound indices
    #
    for k in range(0, n):
        if k > 0:
            ovpair_decl += (", " + virt_idx[k] + occ_idx[k])
        else:
            ovpair_decl += (virt_idx[k] + occ_idx[k])
    #
    # List of indices in inverted alphabetical order
    # (required by subroutine geneargin loops: loops
    # over last ket pair comes first, etc.)
    #
    occ_idx_inv = deepcopy(occ_idx)
    virt_idx_inv = deepcopy(virt_idx)
    occ_idx_inv.sort(reverse=True)
    virt_idx_inv.sort(reverse=True)
    occ_summ_idx_set = set([])
    virt_summ_idx_set = set([])
    #
    # Generate loops covering whole CC Jacobian matrix
    #
    loops = ""

    idxbound_decl_set = set()

    print("EOMCCJAC MODULE. NUMBER OF LOOPS TO PROCESS: {0}".format(len(nonzero_kronecker)))
    loop_number = 1
    for k in range(0, len(nonzero_kronecker)):
        print("LOOP {0}".format(specnames[k]))
        deltas = nonzero_kronecker[k]
        specname = specnames[k]
        #
        # Determine equivalence classes of indices
        #
        deltas_indices = set([])
        substitutions = {}
        occ_summ_idx = []
        virt_summ_idx = []
        #
        # Compute dictionary of index substitutions.
        # Summation indices (free indices) are substituted
        # for themselves by definition, so every index has
        # its key in SUBSTITUTIONS dictionary.
        #
        substitutions, occ_summ_idx, virt_summ_idx = comp_substitutions_dict(braket_type, deltas)
        #
        # Check if any index occures three times
        # or more in bra/ket vector. Such matrix
        # elements are identically equal to zero.
        #
        if is_threefold(braket_type, substitutions):
            print("NOTE: {0} CLASS OF MATRIX ELEMENTS IS SKIPPED. THESE ELEMENTS VANISH IDENTICALLY.".format(specnames[k]))
            #
            # Go to the next class of matrix elements
            #
            continue
        #
        # Skip the loop over current class of matrix elements if bra/ket
        # vector belongs to the triple excitations not present in the 
        # independent triples' basis set (occupied index on repeated
        # twice on the left)
        #
        if independent_triples:
            if isredundant_iik(braket_type, substitutions):
                print("NOTE: {0} CLASS OF MATRIX ELEMENTS IS SKIPPED. BRA/KET VECTOR IS NOT PRESENT IN THE INDEPENDENT TRIPLES' SET.".format(specnames[k]))
                #
                # Go to the next class of matrix elements
                #
                continue
        #
        # Determine variable loop bounds that
        # result from lower-triangularity condition:
        # (ai) >= (bj) in element <aibj| X |ck>,
        # (ai) >= (bj) >= (ck), (dl) >= (ek)
        # in <aibjck| X | dlek>, etc.
        #
        pending = {}
        loopbounds(braket_type, substitutions, pending)
        #
        # Add proper inequalities to skip redundant triples
        # in cases where three different (i, j, k) occupied
        # indices are present
        #
        if independent_triples:
            nonredundant_ijk(braket_type, substitutions, specnames[k], pending)
        #
        # Generate pending code for computing pair indices
        # (compound indices)
        #
        pairindices(braket_type, substitutions, pending)
        #
        # Update global list of variables
        # to declare at the beginning
        # of the subprogram
        #
        occ_summ_idx_set |= set(occ_summ_idx)
        virt_summ_idx_set |= set(virt_summ_idx)
        #
        # Generate comment string in from of 
        # elementary loop
        #
        eqs = []
        for deltapair in deltas:
            f = {"p": deltapair[1], "q": deltapair[0]}
            eqs.append("{p} == {q}".format(**f))

        c = {"virtual": commasep(virt_summ_idx), \
                 "occupied": commasep(occ_summ_idx), \
                 "equalities": commasep(eqs), \
                 "loop_number":loop_number}
        loop_number += 1

        if c["equalities"] == "":
            c["equalities"] = "none"
        

        comment = """
!
! Elementary loop {loop_number}
! --------------------
! Free virtual indices: {virtual}
! Free occupied indices: {occupied}
! Equalities: {equalities}
! No equalities independent of the above can hold.
!
""".format(**c)

        loops += comment
        loops += jacobian_loops(genname, specname, 1, occ_idx_inv, 
                               virt_idx_inv, occ_summ_idx, virt_summ_idx, substitutions, 
                               braket_type, pending, vanish_case1)
        if pendop_count(PEND_INTDECL, pending) > 0:
            idxbound_decl_set |= set(pending[PEND_INTDECL])

    idxbound_decl_list = list(idxbound_decl_set)
    idxbound_decl_list.sort()
    if len(idxbound_decl_list) > 0:
        idxbound_decl = "integer :: " + commasep(idxbound_decl_list)
    else:
        idxbound_decl = ""
    virt_summ_idx_list = list(virt_summ_idx_set)
    occ_summ_idx_list = list(occ_summ_idx_set)
    virt_summ_idx_list.sort()
    occ_summ_idx_list.sort()
    virt_idx_decl = commasep(virt_summ_idx_list)
    occ_idx_decl = commasep(occ_summ_idx_list)

    utc_datetime = datetime.datetime.utcnow()
    datestr = utc_datetime.strftime("%Y-%m-%d %H:%M:%S")
    #
    # Additional variable/function declarations are needed if
    # there are loops over triples
    #
    tidx_decl = ""
    tidx_init = ""
    tidx_func = ""
    if (braket_type[0] == 3) or (braket_type[1] == 3):
        tidx_decl = """
integer :: qbj, qbj2
integer :: qck, qck2
integer :: q00"""
        tidx_init = """
        qbj  = 3 + 6 * npair
        qbj2 = -3
        qck  = 2 + 3 * npair * (2 + npair)
        qck2 = -3 * (1 + npair)
        q00  = -3 * npair * (3 + npair)
"""
        tidx_func = """
contains
!
! Locally visible functions
!
function mu3(ai, bj, ck)
      !
      ! Compute compound three-electron index
      ! (assumed that ai >= bj >= ck)
      !
      integer :: mu3
      integer, intent(in) :: ai, bj, ck
      integer :: mu31, mu32
      !
      ! Compound index is compouted relative
      ! to the first element of matrix block.
      ! Braoffset/ketoffset should be added
      ! to get the absolute position
      !
      mu31 = (qbj + qbj2 * bj) * bj 
      mu32 = (qck + (qck2 + ck) * ck) * ck
      mu3  = ai + (mu31 + mu32 + q00) / 6
end function mu3
"""

    d = {"virt_idx_decl":virt_idx_decl, "occ_idx_decl":occ_idx_decl, "ovpair_decl":ovpair_decl, \
             "loops":loops, "bra_type":bra_type, "ket_type":ket_type, "datestr":datestr, \
             "idxbound_decl":idxbound_decl, "tidx_decl":tidx_decl, "tidx_init":tidx_init, \
             "tidx_func":tidx_func, "module_postfix":module_postfix}

    s = """
module ccjac_block_{bra_type}{ket_type}{module_postfix}
use eom_cc3_{bra_type}{ket_type}_trans
use davidson_main

 
implicit none
!
! File generated automatically on {datestr} UTC.
!
contains
 
subroutine ccjac_{bra_type}{ket_type}{module_postfix}(jac, eorb, t1, t2, nocc0, nocc1, nvirt0, nvirt1, bra0, ket0)
double precision, dimension(:,:), intent(out)       :: jac
double precision, dimension(:), intent(in)          :: eorb
double precision, dimension(:, :), intent(in)       :: t1
double precision, dimension(:, :, :, :), intent(in) :: t2
integer, intent(in)                                 :: nocc0, nocc1
integer, intent(in)                                 :: nvirt0, nvirt1
integer, intent(in)                                 :: bra0, ket0
!
! Local variables
!
integer :: {virt_idx_decl}
integer :: {occ_idx_decl}
integer :: {ovpair_decl}
{idxbound_decl}
integer :: nocc, nvirt
integer :: npair, nactive
integer :: ibra, iket
integer :: braoffset, ketoffset
{tidx_decl}
!
! Offset of the jacobian blocks
!
braoffset = bra0 - 1
ketoffset = ket0 - 1
!
! Number of occupied and virtual orbitals
! present in calculations
!
nocc = nocc1 - nocc0 + 1
nvirt = nvirt1 - nvirt0 + 1
npair = nocc * nvirt
nactive = nocc + nvirt
{tidx_init}
{loops}
{tidx_func}
end subroutine ccjac_{bra_type}{ket_type}{module_postfix}
end module ccjac_block_{bra_type}{ket_type}{module_postfix}
""".format(**d)
    #
    # Write textfiles
    #
    filename = "./fortran_codes/ccjac_bez_dav/ccjac_block_{bra_type}{ket_type}{module_postfix}.f90".format(**d)
    writefile(s, filename)


def demo():
    #
    # Demonstration of CC jacobian loops
    # Compute <aijb| X |ck> element
    #
    braket_type = (3, 3)
    nonzero_kronecker = [[]]
    genname = "eomccjac"
    specnames = ["aibjckdlemfn"]
    
    jacobian_loop(braket_type, nonzero_kronecker, genname, specnames)
