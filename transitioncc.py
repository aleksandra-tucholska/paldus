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
from eomccjac import *

def transition_loop(braket_type, nonzero_kronecker, genname, specnames,  vanish_case1=[], independent_triples=True):
    """
    BRAKET_TYPE       - 2-tuple determining type of matrix element:
                        (0, 3) -> <| X | aibjck>
                        (3, 0) -> <aibjck| X |>
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
        sm = "use ccjac_block_{bra_type}{ket_type}_dav{module_postfix}".format(**f)
        ss = "call ccjac_{bra_type}{ket_type}_dav{module_postfix}(eorb, t1, t2, rt, lf, nocc0, nocc1, nvirt0, nvirt1, bra0, ket0, nvec)".format(**f)
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
module ccjac_block_{bra_type}{ket_type}_dav
{modulelist}
implicit none
!
! File generated automatically on {datestr} UTC.
!
contains
 
subroutine ccjac_{bra_type}{ket_type}_dav(eorb, t1, t2, rt, lf, nocc0, nocc1, nvirt0, nvirt1, bra0, ket0, nvec)
double precision, dimension(:), intent(in)          :: eorb
double precision, dimension(:, :), intent(in)       :: t1
double precision, dimension(:, :, :, :), intent(in) :: t2
double precision, dimension(:, :), intent(in)       :: rt
double precision, dimension(:, :), intent(in)       :: lf
integer, intent(in)                                 :: nocc0, nocc1
integer, intent(in)                                 :: nvirt0, nvirt1
integer, intent(in)                                 :: bra0, ket0
integer, intent(in)                                 :: nvec
!
! Generating full CC Jacobian
!
{calllist}
end subroutine ccjac_{bra_type}{ket_type}_dav
end module ccjac_block_{bra_type}{ket_type}_dav
""".format(**f)
        filename = "./fortran_codes/transition/{trans_name}.f90".format(**f)
        writefile(s, filename)
