#
# Disjoint sets of indices.
# --
# IT IS ASSUMED THROUGHOUT THE CODE THAT FOR ANY
# A IN VIRTUAL, I IN OCCUPIED, AND P IN GENERAL
# THE FOLLOWING RELATION HOLDS:
# A < I < P,
# I.E. A ALWAYS PRECEDES I AND P IN ALPHABETICAL 
# ASCENDING ORDERING.
# --
# GENERAL INDICES SHOULD NOT BE USED AS FIXED 
# INDICES. IT IS ASSUMED THAT FIXED INDICES
# ARE EITHER VIRTUAL OR OCCUPIED.
#

""" the name of the indices can't be combination of numbers
"""

virtual  = ["a", "b", "c", "d", "e", "f", "g", \
                "a>", "b>", "c>", "d>", "e>", "f>", \
                "a>>", "b>>", "c>>", "d>>", "e>>", "f>>", \
                "a>>>", "b>>>", "c>>>", "d>>>", "e>>>", "f>>>"]

fortran_virtual  = ["a1", "a2", "a3", "a4", "a5", "a6", "a7", \
                    "a8", "a9", "a10", "a11", "a12", "a13", "a14"]

fortran_occupied = ["i1", "i2", "i3","i4","i5","i6","i7","i8","i9","i10","i11","i12"]

fortran_CABS = ["X1", "X2", "X3","X4","X5","X6","X7","X8","X9","X10","X11","X12"]

latex_virtual  = ["a_1", "a_2", "a_3", "a_4", "a_5", "a_6", "a_7", \
                    "a_8", "a_9", "a_{10}", "a_{11}", "a_{12}", "a_{13}", "a_{14}"]

latex_occupied = ["i_1", "i_2", "i_3", "i_4", "i_5", "i_6", "i_7", \
                    "i_8", "i_9", "i_{10}", "i_{11}", "i_{12}", "i_{13}", "i_{14}"]


latex_CABS = ["X_1", "X_2", "X_3", "X_4", "X_5", "X_6", "X_7", \
                    "X_8", "X_9", "X_{10}", "X_{11}", "X_{12}", "X_{13}", "X_{14}"]


occupied = ["i", "j", "k", "l", "m", "n", \
                "i>", "j>", "k>", "l>", "m>", "n>",\
                "i>>", "j>>", "k>>", "l>>", "m>>", "n>>", \
                "i>>>", "j>>>", "k>>>", "l>>>", "m>>>", "n>>>"]

general  = ["p", "q", "r", "s", "t", "u", "v", "w","x","y","z",\
                "p>", "q>", "r>", "s>", "t>", "u>",\
                "p>>", "q>>", "r>>", "s>>", "t>>", "u>>"]

complete = ["á", "č", "ě", "í", "ň", "ó", "ř" ]

completev = ["α", "β", "γ", "δ", "ε", "ζ", "η", \
                 "α>", "β>", "γ>", "δ>", "ε>", "ζ>", "η>"]


completev_latex_dict = {"α":"\\alpha ", "β":"\\beta ", "γ":"\\gamma ", "δ":"\\delta ", "ε": "\\epsilon ", "ζ": "\\zeta ", "η": "\\eta", \
                        "α>":"\\alpha' ", "β>":"\\beta' ", "γ>":"\\gamma' ", "δ>":"\\delta' ", \
                        "ε>": "\\epsilon' ", "ζ>": "\\zeta' ", "η>": "\\eta'"}


CABS = ["A", "B", "C", "D", "E", "F", "G", \
            "A>", "B>", "C>", "D>", "E>", "F>", "G>"\
            "A>>", "B>>", "C>>", "D>>", "E>>", "F>>", "G>>"\
            "A>>>", "B>>>", "C>>>", "D>>>", "E>>>", "F>>>", "G>>>"]


#""a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10", "a11", "a12", "a13", "a14"]

#"p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10", "p11", "p12"]
cabs = ["b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10", "b11", "b12", "b13", "b14"]

nsvd = ['X', 'Y', 'Z', 'X>', 'Y>', 'Z>']

aux = ['Q', 'R', 'S', 'T']



Ablock_dict_pp = {
    'ooao': ['o', 'o', 'a', 'o'],
    'aooo': ['a', 'o', 'o', 'o'],
}


Ablock_dict = {

#    'vaao': ['v', 'a', 'a', 'o'],
     'aaaa': ['a', 'a', 'a', 'a']
#    'aova': ['a', 'o', 'v', 'a'],

#    'aoao': ['a', 'o', 'a', 'o'],
#    'aooa': ['a', 'o', 'o', 'a'],
#    'oaao': ['o', 'a', 'a', 'o'],
#    'oaoa': ['o', 'a', 'o', 'a'],

#    'aoaa': ['a', 'o', 'a', 'a'],
#    'oaaa': ['o', 'a', 'a', 'a'],
#    'aaao': ['a', 'a', 'a', 'o'],
#    'aaoa': ['a', 'a', 'o', 'a'],

#    'oaaa': ['o', 'a', 'a', 'a'],
#    'aaao': ['a', 'a', 'a', 'o'],
#    'aaoa': ['a', 'a', 'o', 'a'],
#

#
#    'aova': ['a', 'o', 'v', 'a'],
#    'vaao': ['v', 'a', 'a', 'o'],
#    'vaoo': ['v', 'a', 'o', 'a'],
#    'avao': ['a', 'v', 'a', 'o'],
#    'avoa': ['a', 'v', 'o', 'a'],
    

#    'aoav': ['a', 'o', 'a', 'v'],
#    'oava': ['o', 'a', 'v', 'a'],
#    'oaav': ['o', 'a', 'a', 'v'],
#    
#
    #
#    'aovo': ['a', 'o', 'v', 'o'],
#    'voao': ['v', 'o', 'a', 'o'],
#    'vooa': ['v', 'o', 'o', 'a'],
#    'ovao': ['o', 'v', 'a', 'o'],
#    'ovoa': ['o', 'v', 'o', 'a'],




#    'aoov': ['a', 'o', 'o', 'v'],
#    'oavo': ['o', 'a', 'v', 'o'],
#    'oaov': ['o', 'a', 'o', 'v'],
    #'aovo': ['a', 'o', 'v', 'o'],
#
#
#    'aaaa': ['a', 'a', 'a', 'a'],
#    'aava': ['a', 'a', 'v', 'a'],
#    'vaaa': ['v', 'a', 'a', 'a'],
#    'avaa': ['a', 'v', 'a', 'a'],

#    'aaav': ['a', 'a', 'a', 'v'],

#    'aavo': ['a', 'a', 'v', 'o'],
#    'voaa': ['v', 'o', 'a', 'a'],
#    'ovaa': ['o', 'v', 'a', 'a'],

#    'aaov': ['a', 'a', 'o', 'v'],

#    'vava': ['v', 'a', 'v', 'a'],
#    'vaav': ['v', 'a', 'a', 'v'],
#    'avva': ['a', 'v', 'v', 'a'],
#    'avav': ['a', 'v', 'a', 'v'],

    #    
#
#    'vavo': ['v', 'a', 'v', 'o'],
#    'vova': ['v', 'o', 'v', 'a'],
#    'voav': ['v', 'o', 'a', 'v'],
#    'ovva': ['o', 'v', 'v', 'a'],
#   'ovav': ['o', 'v', 'a', 'v'],



#    'vaov': ['v', 'a', 'o', 'v'],
#    'avvo': ['a', 'v', 'v', 'o'],
#    'avov': ['a', 'v', 'o', 'v'],
    
    
#    'ovov': ['o', 'v', 'o', 'v'],    
#    'vovo': ['v', 'o', 'v', 'o'],
#    'voov': ['v', 'o', 'o', 'v'],
#    'ovvo': ['o', 'v', 'v', 'o'],
#    'ovov': ['o', 'v', 'o', 'v'],   
#
#

}



#-----------------------------------------TORUN INDICES--------------------------------
mona_occup_original = []
monb_occup_original = []
mona_virt_original = []
monb_virt_original = []

mona_occup = ['ia', 'ib', 'ic', 'id', 'ie', 'if', 'ig', 'ih', 'ii', 'ij', 'ik']
monb_occup = ['ja', 'jb', 'jc', 'jd', 'je', 'jf', 'jg', 'jh', 'ji', 'jj', 'jk']
mona_virt =  ['aa', 'ab', 'ac', 'ad', 'ae', 'af', 'ag', 'ah', 'ai', 'aj', 'ak']
monb_virt =  ['ba', 'bb', 'bc', 'bd', 'be', 'bf', 'bg', 'bh', 'bi', 'bj', 'bk']

monab = mona_occup+monb_occup+mona_virt+monb_virt
Tdict = {}

for i in range(0, 11):
    mona_occup_original.append('i'+str(i+1))
    monb_occup_original.append('j'+str(i+1))
    mona_virt_original.append('a'+str(i+1))
    monb_virt_original.append('b'+str(i+1))

    Tdict[mona_occup[i]] = 'i_'+str(i+1)
    Tdict[monb_occup[i]] = 'j_'+str(i+1)
    Tdict[mona_virt[i]]  = 'a_'+str(i+1)
    Tdict[monb_virt[i]]  = 'b_'+str(i+1)
for x in occupied:
    Tdict[x]=x
for x in virtual:
    Tdict[x]=x
for x in general:
    Tdict[x]=x
#------------------------------------------------------------------------------------------------------

fixed  = []

occ_all = occupied + fortran_occupied
virt_all = virtual + fortran_virtual
CABS_all = CABS + fortran_CABS

virtualall = virtual + completev + cabs + CABS

generalall = general + complete + completev

#
# Available matrices/tensors:
# ----------------------------
# OBSERVABLE_X,  - Symmetric real matrices of one-electron operators
# OBSERVABLE_Y,    representing measurable quantities
# OBSERVABLE_Z
# 
# NONSYM_MTRX
#
# FOCK_MATRIX    - Symmetric real matrix of the Fock operator
#
# CC_AMPLITUDE   - Spin restricted coupled cluster amplitude. Number
#                  of declared indices determines whether it is
#                  T1, T2, or higher-order amplitude. Invariant under
#                  interchanging occupied-virtual pairs of indices.
# EOM-CC_AMPLITUDE   - Spin restricted coupled cluster amplitude. Number
#                  of declared indices determines whether it is
#                  r1, r2, or higher-order amplitude. Invariant under
#                  interchanging occupied-virtual pairs of indices.
#
# S_AMPLITUDE    - Spin restricted coupled cluster S amplitude. Number
#                  of declared indices determines whether it is
#                  S1, S2, or higher-order amplitude. Invariant under
#                  interchanging occupied-virtual pairs of indices.
#
# TWOEL_INT      - Real two-electron electron repulsion integral.
#                  Invariant under interchanging indices inside bra
#                  or ket, or interchanging bra and ket pairs.
# TWOEL_INT_DIRAC- Real two-electron electron repulsion integral.
#                  Invariant under interchanging indices according to
#                  Dirac notation <pq|rs> = <rq|ps>
#
# TWOEL_INT_TRANS- Real t1_transformed two-electron electron repulsion integral.
#                  Invariant under interchanging bra and ket pairs.
#
# BARENUCL_HAM   - Real symmetric matrix of bare nuclei Hamiltonian,
#                  that is operator representing the potential of nuclei-
#                  electrons interaction.
#

EINSTEIN = False
CRE = '+'
ANI = '0'
OCC_CRE = 1
OCC_ANI = 2
VIR_CRE = 3
VIR_CRE = 4
ANY_OP                = "b"
SLAT2_SYM1   = "b1"
SLAT2_SYM2   = "b2"
SLAT3        = "b3"
SLAT4        = "b4"
CI_AMPLITUDE = 'c'
COREL                 = 'fcc'
OBSERVABLE_X          = "x"
OBSERVABLE_AX         = "Ax"
OBSERVABLE_AY         = "Ay"
OBSERVABLE_Y          = "y"
OBSERVABLE_Z          = "z"
NONSYM_MTRX           = "q"
OBSERVABLE_X_ASYM     = "u"
FOCK_MATRIX           = "f"
FOCK_MATRIX_TRANS     = "fc"
CC_AMPLITUDE          = "t"
CC_TAU                = "tu"
EOM_CC_AMPLITUDE_R    = "r"
EOM_CC_SINGLE_Rl      = "rl"
EOM_CC_SINGLE_Rr      = "rr"
EOM_CC_SINGLE_Rl_plus   = "rlp"
EOM_CC_SINGLE_Rl_minus  = "rlm"
EOM_CC_SINGLE_Rr_plus   = "rrp"
EOM_CC_SINGLE_Rr_minus  = "rrm"
EOM_CC_SINGLE_Ll      = "Ll"
EOM_CC_SINGLE_Lr      = "Ll"
EOM_CC_AMPLITUDE_L    = "l"
S_AMPLITUDE           = "s"
TWOEL_INT             = "vv"
TWOEL_INT_DIRAC       = "gs"
TWOEL_INT_DIRAC_SPINRES       = "gssr"
TWOEL_INT_DIRAC_SPINRES_AAAA       = "gspp"
TWOEL_INT_DIRAC_SPINRES_BBBB       = "gsmm"
TWOEL_INT_DIRAC_SPINRES_ABAB       = "gspm"
TWOEL_INT_SPINRES                  = "vvsr"
TWOEL_INT_SPINRES_AAAA             = "vvpp"
TWOEL_INT_SPINRES_BBBB             = "vvmm"
TWOEL_INT_SPINRES_ABAB             = "vvpm"
TWOEL_INT_DIRAC_A       = "gsa"
TWOEL_INT_COMBO       = "L"
TWOEL_INT_COMBO2      = "LL"
F12_TWOEL_COMBO       = "F"
F12_TWOEL_COMBO2      = "FF"
TEMP1                 = "D"
TEMP2                 = "J"
TWOEL_INT_TRANS       = "v"
BARENUCL_HAM          = "h"
BARENUCL_HAM_TRANS    = "k"
EOM_TRIPLET_R3        = "R3"
EOM_TRIPLET_R1        = "R1"
EOM_TRIPLET_R2p       = "R2p"
EOM_TRIPLET_R2m       = "R2m"

DENS1                 = "gm1"
DENS1P                 = "gm1p"
DENS1M                 = "gm1m"
DENS2                 = "Gm2"
CUM2                  = "Cm2"
CUM3                  = "Cm3"
DENS2P                = "Gm2p"
DENS2M                = "Gm2m"
DENS2PM               = "Gm2pm"
DENS2MP               = "Gm2mp"
DENS3                 = "Gm3"
DENS4                 = "Gm4"

DENSN                   = "nn"
DENSN1                  = "n"
DENSN1P                 = "np"
DENSN1M                 = "nm"

CL2PM                   = "CL2pm"
CL2PP                   = "CL2pp"

CL2                     = "CL2"
CL3                     = "CL3"

CL3PM                   = "CL2pm"
CL3PP                   = "CL2pp"
ferm = 'F'
boso = 'B'


DENS3PM                 = "Gm3pm"
DENS3P                 = "Gm3p"
DENS3M                 = "Gm3m"
DENS4PPM                 = "Gm4ppm"
DENS4PM                 = "Gm4pm"
DENS4P                 = "Gm4p"
DENS4M                 = "Gm4m"

F12_AMPLITUDE         = "tf"
F12_SAMPLITUDE         = "sf"
F12_TWOEL             = "ff"
INTERM_V_1_F12          = "V_1"
INTERM_V_F12          = "V"
INTERM_X_F12          = "X"
INTERM_X_1_F12          = "X_1"
INTERM_Z_F12          = "Z"
INTERM_B_F12          = "B"
INTERM_B_1_F12          = "B_1"
INTERM_P_F12          = "P"
INTERM_P_1_F12          = "P_1"
INTERM_Ft_F12          = "Ft"
INTERM_Ftt_F12          = "Ftt"

SVD_B = "B"
SVD_C = "C"
SVD_Y = "Y"
SVD_X = "X"
SVD_U = "U"
SVD_A = "A"
SVD_Z = "Z"
SVD_t = "ttt"
SVD_chi = "CHI"

AMP_NAME_DICT = {FOCK_MATRIX_TRANS:0, FOCK_MATRIX:1, F12_TWOEL:2, F12_TWOEL_COMBO:3, \
BARENUCL_HAM:4, BARENUCL_HAM_TRANS:5, TWOEL_INT:6, TWOEL_INT_COMBO:7, \
TWOEL_INT_TRANS :8,CC_AMPLITUDE: 9, F12_AMPLITUDE:10, S_AMPLITUDE:11, \
TWOEL_INT_TRANS:12, 'Pp':13, 'Pm':14, NONSYM_MTRX:15, \
                 TWOEL_INT_COMBO2 :16, F12_TWOEL_COMBO2: 17, INTERM_Ft_F12:18, INTERM_V_F12:19, INTERM_X_F12:20, INTERM_B_F12:21, INTERM_P_F12:22,
                 INTERM_V_1_F12:23, INTERM_B_1_F12:24, INTERM_X_1_F12:25, INTERM_P_1_F12:26, OBSERVABLE_X:27}

#--------tg-------------
TWOEL_INT_AS          = "w"
#
# Subset of matrices/tensors that have indices grouped
# into consecutive electron pairs
#
DENS2_SYM = set([DENS2, DENS2P, DENS2M, DENS2MP, DENS2PM])
TWOEL = set([TWOEL_INT, TWOEL_INT_TRANS])
HAS_PAIR_INDICES = set([OBSERVABLE_X, OBSERVABLE_Y, OBSERVABLE_Z, \
                            FOCK_MATRIX, CC_AMPLITUDE, TWOEL_INT, TWOEL_INT_TRANS, F12_TWOEL, BARENUCL_HAM, \
                            EOM_CC_AMPLITUDE_R, EOM_CC_AMPLITUDE_L, EOM_CC_SINGLE_Rl, EOM_CC_SINGLE_Rr, EOM_CC_SINGLE_Lr, EOM_CC_SINGLE_Ll])
SYMMETRIC_MATRIX = set([OBSERVABLE_X, OBSERVABLE_Y, OBSERVABLE_Z, FOCK_MATRIX, BARENUCL_HAM])

NONSYMMETRIC_MATRIX = set([NONSYM_MTRX])

ASYMMETRIC_MATRIX = set([OBSERVABLE_X_ASYM])

OBSERVABLE_MATRIX = set([OBSERVABLE_X, OBSERVABLE_Y, OBSERVABLE_Z, OBSERVABLE_X_ASYM])
#
# Terms with numerical prefactor less than EPSILON are
# regarded negligible and eliminated from calculations / display
#
EPSILON = 10**(-14)

tex_preamble = """
\\documentclass[aip,jcp,preprint,amsmath,amssymb, groupedaddress]{revtex4-1}
\\usepackage[T1]{fontenc}
\\usepackage[utf8]{inputenc}
\\usepackage{amsfonts}
\\usepackage{multirow,dcolumn}
\\usepackage{booktabs}
\\usepackage[version=3]{mhchem}
\\usepackage{setspace} 
\\usepackage{lmodern}
\\usepackage{leftidx}
\\usepackage{mathtools}
\\usepackage{caption}
\\usepackage{tabularx}
\\usepackage{threeparttable}
\\usepackage[dvipsnames]{xcolor}
\\usepackage{blindtext}
\\usepackage{booktabs}
\\usepackage[version=3]{mhchem}
\\usepackage{lmodern}
\\usepackage{natbib}
\\usepackage[pdftex]{graphicx}
\\usepackage[amssymb]{SIunits} 

\\usepackage{paralist} 
\\usepackage[normalem]{ulem}
\\usepackage{xcccom}
\\usepackage{braket}

"""
