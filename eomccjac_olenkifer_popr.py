#
# Automated coding of loops over elements of
# the EOM-CC Jacobian matrix
# ------------------------------------------------------
# Author: Marcin Modrzejewski, University of Warsaw
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
# elements for which ai >= bj holds are commputed.
#
# NPAIR - Number of distinct one electron indices
# NPAIR = NOCC * NVIRT
# 
# aibj = [(2 * NPAIR - bj + 2) * (bj - 1)] / 2 + ai - bj + 1
#
from copy import deepcopy
from paldus_classes import occupied
from paldus_classes import virtual
import datetime
import sys
#
# Constants defining types of pending operations
#
PEND_COMPIDX = 1
PEND_LOBOUND = 2
PEND_UPBOUND = 3
PEND_INTDECL = 4

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


def pendop_addinequality(ineq, pending_dict):
    """
    Process inequality to get loop bounds. INEQ represents
    inequality as a list of indices:
    INEQ = ["a", "b"] # a >= b
    """
    a = ineq[0]
    b = ineq[1]

    if a < b:
        pendop_add((a, PEND_LOBOUND), b, pending_dict)
    else:
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
    return s


def ketidx(braket_type):
    kettype = braket_type[1]
    ketidx1 = braket_type[0]
    ketidx2 = braket_type[0] + 1
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
    return s
        


def jacobian_loops(genname, specname, level, occ_idx, 
                   virt_idx, occ_summ_idx, virt_summ_idx, substitutions, 
                   braket_type, pending):
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
    occ_skip_nondiag_list = []
    if nlo == 0:
        i0 = "nocc0"
    elif nlo == 1:
        i0 = pending[(occ_this, PEND_LOBOUND)][0] + " + 1"
        occ_skip_nondiag_list += pending[(occ_this, PEND_LOBOUND)]
    else:
        i0 = occ_this + "0"
        occ_pending_lobound = maxval_code(pending[(occ_this, PEND_LOBOUND)], i0)
        pendop_add(PEND_INTDECL, i0, pending)
        occ_skip_nondiag_list += pending[(occ_this, PEND_LOBOUND)]
    d["i0"] = i0
    d["occ_pending_lobound"] = occ_pending_lobound

    occ_pending_upbound = ""
    i1 = ""
    nup = pendop_count((occ_this, PEND_UPBOUND), pending)
    if nup == 0:
        i1 = "nocc1"
    elif nup == 1:
        i1 = pending[(occ_this, PEND_UPBOUND)][0]
    else:
        i1 = occ_this + "1"
        occ_pending_upbound = minval_code(pending[(occ_this, PEND_UPBOUND)], i1)
        pendop_add(PEND_INTDECL, i1, pending)
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
        a0 = virt_this + "0"
        virt_pending_lobound = maxval_code(pending[(virt_this, PEND_LOBOUND)], a0)
        pendop_add(PEND_INTDECL, a0, pending)
        virt_skip_nondiag_list += pending[(virt_this, PEND_LOBOUND)]
    d["a0"] = a0
    d["virt_pending_lobound"] = virt_pending_lobound

    virt_pending_upbound = ""
    a1 = ""
    nup = pendop_count((virt_this, PEND_UPBOUND), pending)
    if nup == 0:
        a1 = "nvirt1"
    elif nup == 1:
        a1 = pending[(virt_this, PEND_UPBOUND)][0]
    else:
        a1 = virt_this + "1"
        virt_pending_upbound = minval_code(pending[(virt_this, PEND_UPBOUND)], a1)
        pendop_add(PEND_INTDECL, a1, pending)
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
                                   braket_type, pending)
    else:
        #
        # Determine list of arguments for called function
        #
        args = ""
        for k in range(maxlevel - 1, -1, -1):
            if virt_idx[k] in virt_summ_idx:
                if len(args) > 0:
                    args += ", "
                args += virt_idx[k]
            if occ_idx[k] in occ_summ_idx:
                if len(args) > 0:
                    args += ", "
                args += occ_idx[k]

        d["args"] = args
        innercode = """
{ibra}
{iket}
jacbig = {genname}_{specname}(eorb, t2, t1, {args})

if(ibra.eq.iket)diag(ibra) = jacbig

do p = 1, m
abm(ibra, p)  = abm(ibra, p)  +  jacbig * rt(iket, p)
abmb(iket, p) = abmb(iket, p) +  jacbig * lf(ibra, p)
end do
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


def loopbounds(braket_type, substitutions, pending_dict):
    firstket = braket_type[0]
    n = braket_type[0] + braket_type[1]
    for k in range(0, n - 1):
        if (k + 1) == firstket:
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
            i = occupied[k]
            j = occupied[k + 1]
            subst_i = substitutions[i]
            subst_j = substitutions[j]
            if subst_i != subst_j:
                pendop_addinequality([subst_i, subst_j], pending_dict)


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


def jacobian_loop(braket_type, nonzero_kronecker, genname, specnames):
    """
    BRAKET_TYPE       - 2-tuple determining type of matrix element:
                        (1, 1) -> <ia| X | jb>
                        (2, 1) -> <iajb| X |kc>
                        (2, 2) -> <iajb| X |kcld>
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
                        """

    #
    # Verify whether NONZERO_KRONECKER list is in canonical form
    #
    iscanonical = verify_deltalist(nonzero_kronecker)
    if not iscanonical:
        sys.exit()
    
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

    for k in range(0, len(nonzero_kronecker)):
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
        #
        # Determine variable loop bounds that
        # result from lower-triangularity condition:
        # (ai) >= (bj) in element <aibj| X |ck> etc.
        #
        pending = {}
        loopbounds(braket_type, substitutions, pending)
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
                 "k": k + 1}

        if c["equalities"] == "":
            c["equalities"] = "none"
        

        comment = """
!
! Elementary loop {k}
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
                               braket_type, pending)
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
    
    d = {"virt_idx_decl":virt_idx_decl, "occ_idx_decl":occ_idx_decl, "ovpair_decl":ovpair_decl, \
             "loops":loops, "bra_type":bra_type, "ket_type":ket_type, "datestr":datestr, \
             "idxbound_decl":idxbound_decl}

    s = """
module ccjac_block_{bra_type}{ket_type}
use eom_{bra_type}{ket_type}
 
implicit none
!
! File generated automatically on {datestr} UTC.
!
contains

subroutine ccjac_{bra_type}{ket_type}(eorb, t1, t2, rt, lf, nocc0, nocc1, nvirt0, nvirt1, bra0, ket0, m)
double precision, dimension(:), intent(in)          :: eorb
double precision, dimension(:, :), intent(in)       :: t1
double precision, dimension(:, :, :, :), intent(in) :: t2
double precision, dimension(:, :), intent(in)       :: rt
double precision, dimension(:, :), intent(in)       :: lf
integer, intent(in)                                 :: nocc0, nocc1
integer, intent(in)                                 :: nvirt0, nvirt1
integer, intent(in)                                 :: bra0, ket0
integer, intent(in)                                 :: m
!
! Local variables
!
integer :: {virt_idx_decl}
integer :: {occ_idx_decl}
integer :: {ovpair_decl}
{idxbound_decl}
integer :: nocc, nvirt
integer :: npair
integer :: ibra, iket
integer :: braoffset, ketoffset
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

{loops}

end subroutine ccjac_{bra_type}{ket_type}
end module ccjac_block_{bra_type}{ket_type}
""".format(**d)

    filename = "./ccjac_block_{bra_type}{ket_type}.f90".format(**d)
    f = open(filename, "w")
    #
    # Skip blank lines
    #
    lines = s.splitlines()
    for line in lines:
        if line != "":
            f.write(line + "\n")
    f.close()


def demo():
    #
    # Demonstration of CC jacobian loops
    # Compute <aijb| X |ck> element
    #
    braket_type = (2, 2)
    nonzero_kronecker = [[]]
    genname = "eomccjac"
    specnames = ["aibjckdl"]
    
    jacobian_loop(braket_type, nonzero_kronecker, genname, specnames)
