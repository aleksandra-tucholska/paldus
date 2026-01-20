from params import *
from paldus_classes import ugg
from paldus_cas import cas
from paldus_classes import arithmetic_string
from math import sqrt
from copy import deepcopy
#
# Real symmetric matrix of a one-electron observable
#
obs = ugg()
obs.summation = ['p','q']
obs.coefficient = [OBSERVABLE_X]
obs.coefficient_idx.append(['p','q'])
obs.operator_idx.append(['p','q'])
obs.operator_type.append("s")

obsy = ugg()
obsy.summation = ['p','q']
obsy.coefficient = [OBSERVABLE_Y]
obsy.coefficient_idx.append(['p','q'])
obsy.operator_idx.append(['p','q'])
obsy.operator_type.append("s")


obst = ugg()
obst.summation = ['p','q']
obst.coefficient = [OBSERVABLE_X]
obst.coefficient_idx.append(['p','q'])
obst.operator_idx.append(['p','q'])
obst.operator_type.append("t0")

aobst = ugg()
aobst.summation = ['p','q']
aobst.coefficient = [OBSERVABLE_X_ASYM]
aobst.coefficient_idx.append(['p','q'])
aobst.operator_idx.append(['p','q'])
aobst.operator_type.append("t0")


obsij = ugg()
obsij.summation = ['i','j']
obsij.coefficient = [OBSERVABLE_X]
obsij.coefficient_idx.append(['i','j'])
obsij.operator_idx.append(['i','j'])
obsij.operator_type.append("s")

obsia = ugg()
obsia.summation = ['i','a']
obsia.coefficient = [OBSERVABLE_X]
obsia.coefficient_idx.append(['i','a'])
obsia.operator_idx.append(['i','a'])
obsia.operator_type.append("s")


obsai = ugg()
obsai.summation = ['a','i']
obsai.coefficient = [OBSERVABLE_X]
obsai.coefficient_idx.append(['a','i'])
obsai.operator_idx.append(['a','i'])
obsai.operator_type.append("s")


obsab = ugg()
obsab.summation = ['a','b']
obsab.coefficient = [OBSERVABLE_X]
obsab.coefficient_idx.append(['a','b'])
obsab.operator_idx.append(['a','b'])
obsab.operator_type.append("s")



nobs = ugg()
nobs.summation = ['p','q']
nobs.coefficient = [NONSYM_MTRX]
nobs.coefficient_idx.append(['p','q'])
nobs.operator_idx.append(['p','q'])
nobs.operator_type.append("s")

nobst = ugg()
nobst.summation = ['p','q']
nobst.coefficient = [NONSYM_MTRX]
nobst.coefficient_idx.append(['p','q'])
nobst.operator_idx.append(['p','q'])
nobst.operator_type.append("t0")

obij = ugg()
obij.operator_idx.append(['i','j'])
obij.operator_type.append("s")

obia = ugg()
obia.operator_idx.append(['i','a'])
obia.operator_type.append("s")

obai = ugg()
obai.operator_idx.append(['a','i'])
obai.operator_type.append("s")

obab = ugg()
obab.operator_idx.append(['a','b'])
obab.operator_type.append("s")

obs2 = ugg()
obs2.summation = ['p','q']
obs2.coefficient = [OBSERVABLE_X]
obs2.coefficient_idx.append(['p','q'])
obs2.operator_idx.append(['p','q'])
obs2.operator_type.append("s")

eomr1 = ugg()
eomr1.summation = ['a','i']
eomr1.coefficient = [EOM_CC_AMPLITUDE_R]
eomr1.coefficient_idx.append(['a','i'])
eomr1.operator_idx.append(['a','i'])
eomr1.operator_type.append("s")

eomrl1 = deepcopy(eomr1)
eomrl1.coefficient = [EOM_CC_SINGLE_Rl]

eomrr1 = deepcopy(eomr1)
eomrr1.coefficient = [EOM_CC_SINGLE_Rr]

eomrr1_triplet = ugg()
eomrr1_triplet.summation = ['a','i']
eomrr1_triplet.coefficient = [EOM_CC_SINGLE_Rr]
eomrr1_triplet.coefficient_idx.append(['a','i'])
eomrr1_triplet.operator_idx.append(['a','i'])
eomrr1_triplet.operator_type.append("t0")



eomrl1_triplet = ugg()
eomrl1_triplet.summation = ['a','i']
eomrl1_triplet.coefficient = [EOM_CC_SINGLE_Rl]
eomrl1_triplet.coefficient_idx.append(['a','i'])
eomrl1_triplet.operator_idx.append(['a','i'])
eomrl1_triplet.operator_type.append("t0")

eoml1 = ugg()
eoml1.summation = ['a','i']
eoml1.coefficient = [EOM_CC_AMPLITUDE_L]
eoml1.coefficient_idx.append(['a','i'])
eoml1.operator_idx.append(['a','i'])
eoml1.operator_type.append("s")

eomll1 = deepcopy(eoml1)
eomll1.coefficient = [EOM_CC_SINGLE_Ll]

eomlr1 = deepcopy(eoml1)
eomlr1.coefficient = [EOM_CC_SINGLE_Lr]

eomr2 = ugg()
eomr2.summation = ['a','i', 'b', 'j']
eomr2.coefficient = [EOM_CC_AMPLITUDE_R]
eomr2.coefficient_idx.append(['a','i', 'b', 'j'])
eomr2.operator_idx.append(['a','i'])
eomr2.operator_idx.append(['b','j'])
eomr2.operator_type.append("s")
eomr2.operator_type.append("s")
eomr2.num_factor = 1./2.

eomrl2 = deepcopy(eomr2)
eomrl2.coefficient = [EOM_CC_SINGLE_Rl]

eomrl2_triplet_plus = ugg()
eomrl2_triplet_plus.summation = ['a','i', 'b', 'j']
eomrl2_triplet_plus.coefficient = [EOM_CC_SINGLE_Rl_plus]
eomrl2_triplet_plus.coefficient_idx.append(['a','i', 'b', 'j'])
eomrl2_triplet_plus.operator_idx.append(['a','i'])
eomrl2_triplet_plus.operator_idx.append(['b','j'])
eomrl2_triplet_plus.operator_type.append("t0")
eomrl2_triplet_plus.operator_type.append("s")
eomrl2_triplet_plus.num_factor = 1./2.

eomrl2_triplet_minus = ugg()
eomrl2_triplet_minus.summation = ['a','i', 'b', 'j']
eomrl2_triplet_minus.coefficient = [EOM_CC_SINGLE_Rl_minus]
eomrl2_triplet_minus.coefficient_idx.append(['a','i', 'b', 'j'])
eomrl2_triplet_minus.operator_idx.append(['a','i'])
eomrl2_triplet_minus.operator_idx.append(['b','j'])
eomrl2_triplet_minus.operator_type.append("t0")
eomrl2_triplet_minus.operator_type.append("s")

eomrr2 = deepcopy(eomr2)
eomrr2.coefficient = [EOM_CC_SINGLE_Rr]


eomrr2_triplet_plus = ugg()
eomrr2_triplet_plus.summation = ['a','i', 'b', 'j']
eomrr2_triplet_plus.coefficient = [EOM_CC_SINGLE_Rr_plus]
eomrr2_triplet_plus.coefficient_idx.append(['a','i', 'b', 'j'])
eomrr2_triplet_plus.operator_idx.append(['a','i'])
eomrr2_triplet_plus.operator_idx.append(['b','j'])
eomrr2_triplet_plus.operator_type.append("t0")
eomrr2_triplet_plus.operator_type.append("s")
eomrr2_triplet_plus.num_factor = 1./2.

eomrr2_triplet_minus = ugg()
eomrr2_triplet_minus.summation = ['a','i', 'b', 'j']
eomrr2_triplet_minus.coefficient = [EOM_CC_SINGLE_Rr_minus]
eomrr2_triplet_minus.coefficient_idx.append(['a','i', 'b', 'j'])
eomrr2_triplet_minus.operator_idx.append(['a','i'])
eomrr2_triplet_minus.operator_idx.append(['b','j'])
eomrr2_triplet_minus.operator_type.append("t0")
eomrr2_triplet_minus.operator_type.append("s")


eoml2 = ugg()
eoml2.summation = ['a','i', 'b', 'j']
eoml2.coefficient = [EOM_CC_AMPLITUDE_L]
eoml2.coefficient_idx.append(['a','i', 'b', 'j'])
eoml2.operator_idx.append(['a','i'])
eoml2.operator_idx.append(['b','j'])
eoml2.operator_type.append("s")
eoml2.operator_type.append("s")
eoml2.num_factor = 1./2.

eomll2 = deepcopy(eoml2)
eomll2.coefficient = [EOM_CC_SINGLE_Ll]

eomlr2 = deepcopy(eoml2)
eomlr2.coefficient = [EOM_CC_SINGLE_Lr]

eomr3 = ugg()
eomr3.summation = ['a','i', 'b', 'j','c','k']
eomr3.coefficient = [EOM_CC_AMPLITUDE_R]
eomr3.coefficient_idx.append(['a','i', 'b', 'j','c','k'])
eomr3.operator_idx.append(['a','i'])
eomr3.operator_idx.append(['b','j'])
eomr3.operator_idx.append(['c','k'])
eomr3.operator_type.append("s")
eomr3.operator_type.append("s")
eomr3.operator_type.append("s")
eomr3.num_factor = 1./6.

eomrl3 = deepcopy(eomr3)
eomrl3.coefficient = [EOM_CC_SINGLE_Rl]

eomrr3 = deepcopy(eomr3)
eomrr3.coefficient = [EOM_CC_SINGLE_Rr]


eomrr3_triplet = ugg()
eomrr3_triplet.summation = ['a','i', 'b', 'j','c','k']
eomrr3_triplet.coefficient = [EOM_CC_SINGLE_Rr]
eomrr3_triplet.coefficient_idx.append(['a','i', 'b', 'j','c','k'])
eomrr3_triplet.operator_idx.append(['a','i'])
eomrr3_triplet.operator_idx.append(['b','j'])
eomrr3_triplet.operator_idx.append(['c','k'])
eomrr3_triplet.operator_type.append("s")
eomrr3_triplet.operator_type.append("t0")
eomrr3_triplet.operator_type.append("s")
eomrr3_triplet.num_factor = 1./2.


eoml3 = ugg()
eoml3.summation = ['a','i', 'b', 'j','c','k']
eoml3.coefficient = [EOM_CC_AMPLITUDE_L]
eoml3.coefficient_idx.append(['a','i', 'b', 'j','c','k'])
eoml3.operator_idx.append(['a','i'])
eoml3.operator_idx.append(['b','j'])
eoml3.operator_idx.append(['c','k'])
eoml3.operator_type.append("s")
eoml3.operator_type.append("s")
eoml3.operator_type.append("s")
eoml3.num_factor = 1./6.

eomll3 = deepcopy(eoml3)
eomll3.coefficient = [EOM_CC_SINGLE_Ll]

eomlr3 = deepcopy(eoml3)
eomlr3.coefficient = [EOM_CC_SINGLE_Lr]

eomr1_amp = ugg()
eomr1_amp.summation = ['a','i']
eomr1_amp.coefficient = [EOM_CC_AMPLITUDE_R]
eomr1_amp.coefficient_idx.append(['a','i'])

eomr2_amp = ugg()
eomr2_amp.summation = ['a','i', 'b', 'j']
eomr2_amp.coefficient = [EOM_CC_AMPLITUDE_R]
eomr2_amp.coefficient_idx.append(['a','i', 'b', 'j'])
eomr2_amp.num_factor = 1./2.

eomr3_amp = ugg()
eomr3_amp.summation = ['a','i', 'b', 'j','c','k']
eomr3_amp.coefficient = [EOM_CC_AMPLITUDE_R]
eomr3_amp.coefficient_idx.append(['a','i', 'b', 'j','c','k'])
eomr3_amp.num_factor = 1./6.


a1 = ugg()
a1.summation = ['a','i', 'b', 'j','c','k']
a1.coefficient = [EOM_CC_AMPLITUDE_L]
a1.coefficient_idx.append(['a','i', 'b', 'j','c','k'])


a2 = ugg()
a2.summation = ['a','i', 'b', 'j','c','k']
a2.coefficient = [EOM_CC_AMPLITUDE_L]
a2.coefficient_idx.append(['b','j', 'a', 'i','c','k'])

obs3 = ugg()
obs3.operator_idx.append(['b','j'])
obs3.operator_type.append("s")
#
# Bare-nuclei hamiltonian
#

obsc = ugg()
obsc.summation = ["á","č"]
obsc.coefficient = [OBSERVABLE_X]
obsc.coefficient_idx.append(["á","č"])
obsc.operator_idx.append(["á","č"])
obsc.operator_type.append("s")


obscy = ugg()
obscy.summation = ["ě","í"]
obscy.coefficient = [OBSERVABLE_Y]
obscy.coefficient_idx.append(["ě","í"])
obscy.operator_idx.append(["ě","í"])
obscy.operator_type.append("s")


hc = ugg()
hc.summation = ["á","č"]
hc.coefficient = [BARENUCL_HAM]
hc.coefficient_idx.append(["á","č"])
hc.operator_idx.append(["á","č"])
hc.operator_type.append("s")

h = ugg()
h.summation = ["p","q"]
h.coefficient = [BARENUCL_HAM]
h.coefficient_idx.append(["p","q"])
h.operator_idx.append(["p","q"])
h.operator_type.append("s")

# hc = ugg()
# hc.summation = ["á","p2"]
# hc.coefficient = [BARENUCL_HAM]
# hc.coefficient_idx.append(["á","p2"])
# hc.operator_idx.append(["á","p2"])
# hc.operator_type.append("s")


ht = ugg()
ht.summation = ["p","q"]
ht.coefficient = [BARENUCL_HAM_TRANS]
ht.coefficient_idx.append(["p","q"])
ht.operator_idx.append(["p","q"])
ht.operator_type.append("s")

htc = ugg()
htc.summation = ["á","č"]
htc.coefficient = [BARENUCL_HAM_TRANS]
htc.coefficient_idx.append(["á","č"])
htc.operator_idx.append(["á","č"])
htc.operator_type.append("s")

#
# Two-electron part of the electronic Hamiltonian
#
g_1 = ugg()
g_1.summation = ["p","q","r","s"]
g_1.coefficient = [TWOEL_INT]
g_1.coefficient_idx.append(["p","q","r","s"])
g_1.operator_idx.append(["p", "q"])
g_1.operator_idx.append(["r", "s"])
g_1.operator_type.append("s")
g_1.operator_type.append("s")
g_1.num_factor = 1./2.

g_2 = ugg()
g_2.summation = ["p","q","r","s"]
g_2.coefficient = [TWOEL_INT]
g_2.coefficient_idx.append(["p","q","r","s"])
g_2.operator_idx.append(["p", "s"])
g_2.operator_type.append("s")
g_2.delta.append(["q","r"])
g_2.num_factor = -1./2.

g_1c = ugg()
g_1c.summation = ["á","č","ě","í"]
g_1c.coefficient = [TWOEL_INT]
g_1c.coefficient_idx.append(["á","č","ě","í"])
g_1c.operator_idx.append(["á", "č"])
g_1c.operator_idx.append(["ě", "í"])
g_1c.operator_type.append("s")
g_1c.operator_type.append("s")
g_1c.num_factor = 1./2.

g_2c = ugg()
g_2c.summation = ["á","č","ě","í"]
g_2c.coefficient = [TWOEL_INT]
g_2c.coefficient_idx.append(["á","č","ě","í"])
g_2c.operator_idx.append(["á", "í"])
g_2c.operator_type.append("s")
g_2c.delta.append(["č","ě"])
g_2c.num_factor = -1./2.

g_1t = ugg()
g_1t.summation = ["p","q","r","s"]
g_1t.coefficient = [TWOEL_INT_TRANS]
g_1t.coefficient_idx.append(["p","q","r","s"])
g_1t.operator_idx.append(["p", "q"])
g_1t.operator_idx.append(["r", "s"])
g_1t.operator_type.append("s")
g_1t.operator_type.append("s")
g_1t.num_factor = 1./2.

g_2t = ugg()
g_2t.summation = ["p","q","r","s"]
g_2t.coefficient = [TWOEL_INT_TRANS]
g_2t.coefficient_idx.append(["p","q","r","s"])
g_2t.operator_idx.append(["p", "s"])
g_2t.operator_type.append("s")
g_2t.delta.append(["q","r"])
g_2t.num_factor = -1./2.

g_1tc = ugg()
g_1tc.summation = ["á","č","ě","í"]
g_1tc.coefficient = [TWOEL_INT_TRANS]
g_1tc.coefficient_idx.append(["á","č","ě","í"])
g_1tc.operator_idx.append(["á", "č"])
g_1tc.operator_idx.append(["ě", "í"])
g_1tc.operator_type.append("s")
g_1tc.operator_type.append("s")
g_1tc.num_factor = 1./2.

g_2tc = ugg()
g_2tc.summation = ["á","č","ě","í"]
g_2tc.coefficient = [TWOEL_INT_TRANS]
g_2tc.coefficient_idx.append(["á","č","ě","í"])
g_2tc.operator_idx.append(["á", "í"])
g_2tc.operator_type.append("s")
g_2tc.delta.append(["č","ě"])
g_2tc.num_factor = -1./2.

#
# f_implicit = f_pq
#    \sum_{pq} f_implicit E_{pq}
# f_explicit = h_{pq} + \sum_i(2 * g_{pqii} - g_{piiq})
#    \sum_{pq} f_explicit E_{pq}
#
fock_imp = ugg()
fock_imp.summation = ["p","q"]
fock_imp.coefficient = [FOCK_MATRIX]
fock_imp.coefficient_idx.append(["p","q"])
fock_imp.operator_idx.append(["p","q"])
fock_imp.operator_type.append("s")

fock_impc = ugg()
fock_impc.summation = ["á","č"]
fock_impc.coefficient = [FOCK_MATRIX]
fock_impc.coefficient_idx.append(["á","č"])
fock_impc.operator_idx.append(["á","č"])
fock_impc.operator_type.append("s")

fock_impt = ugg()
fock_impt.summation = ["p","q"]
fock_impt.coefficient = [FOCK_MATRIX_TRANS]
fock_impt.coefficient_idx.append(["p","q"])
fock_impt.operator_idx.append(["p","q"])
fock_impt.operator_type.append("s")

fock_imptc = ugg()
fock_imptc.summation = ["á","č"]
fock_imptc.coefficient = [FOCK_MATRIX_TRANS]
fock_imptc.coefficient_idx.append(["á","č"])
fock_imptc.operator_idx.append(["á","č"])
fock_imptc.operator_type.append("s")


fock_exp1 = ugg()
fock_exp1.summation = ["p","q"]
fock_exp1.coefficient = [BARENUCL_HAM]
fock_exp1.coefficient_idx.append(["p", "q"])
fock_exp1.operator_idx.append(["p", "q"])
fock_exp1.operator_type.append("s")
fock_exp1.num_factor = -1.0

fock_exp1c = ugg()
fock_exp1c.summation = ["á","č"]
fock_exp1c.coefficient = [BARENUCL_HAM]
fock_exp1c.coefficient_idx.append(["á", "č"])
fock_exp1c.operator_idx.append(["á", "č"])
fock_exp1c.operator_type.append("s")
fock_exp1c.num_factor = -1.0

fock_exp1tc = ugg()
fock_exp1tc.summation = ["á","č"]
fock_exp1tc.coefficient = [BARENUCL_HAM_TRANS]
fock_exp1tc.coefficient_idx.append(["á", "č"])
fock_exp1tc.operator_idx.append(["á", "č"])
fock_exp1tc.operator_type.append("s")
fock_exp1tc.num_factor = -1.0


fock_exp1t = ugg()
fock_exp1t.summation = ["p","q"]
fock_exp1t.coefficient = [BARENUCL_HAM_TRANS]
fock_exp1t.coefficient_idx.append(["p", "q"])
fock_exp1t.operator_idx.append(["p", "q"])
fock_exp1t.operator_type.append("s")
fock_exp1t.num_factor = -1.0

fock_exp2 = ugg()
fock_exp2.summation = ["p","q", "i"]
fock_exp2.coefficient = [TWOEL_INT]
fock_exp2.coefficient_idx.append(["p", "q", "i", "i"])
fock_exp2.operator_idx.append(["p", "q"])
fock_exp2.operator_type.append("s")
fock_exp2.num_factor = -2.0

fock_exp3 = ugg()
fock_exp3.summation = ["p","q", "i"]
fock_exp3.coefficient = [TWOEL_INT]
fock_exp3.coefficient_idx.append(["p", "i", "i","q"])
fock_exp3.operator_idx.append(["p", "q"])
fock_exp3.operator_type.append("s")
fock_exp3.num_factor = 1.0

fock_exp2c = ugg()
fock_exp2c.summation = ["á","č", "i"]
fock_exp2c.coefficient = [TWOEL_INT]
fock_exp2c.coefficient_idx.append(["á", "č", "i", "i"])
fock_exp2c.operator_idx.append(["á", "č"])
fock_exp2c.operator_type.append("s")
fock_exp2c.num_factor = -2.0

fock_exp3c = ugg()
fock_exp3c.summation = ["á","č", "i"]
fock_exp3c.coefficient = [TWOEL_INT]
fock_exp3c.coefficient_idx.append(["á", "i", "i","č"])
fock_exp3c.operator_idx.append(["á", "č"])
fock_exp3c.operator_type.append("s")
fock_exp3c.num_factor = 1.0

fock_exp2tc = ugg()
fock_exp2tc.summation = ["á","č", "i"]
fock_exp2tc.coefficient = [TWOEL_INT_TRANS]
fock_exp2tc.coefficient_idx.append(["á", "č", "i", "i"])
fock_exp2tc.operator_idx.append(["á", "č"])
fock_exp2tc.operator_type.append("s")
fock_exp2tc.num_factor = -2.0

fock_exp3tc = ugg()
fock_exp3tc.summation = ["á","č", "i"]
fock_exp3tc.coefficient = [TWOEL_INT_TRANS]
fock_exp3tc.coefficient_idx.append(["á", "i", "i","č"])
fock_exp3tc.operator_idx.append(["á", "č"])
fock_exp3tc.operator_type.append("s")
fock_exp3tc.num_factor = 1.0




fock_exp2t = ugg()
fock_exp2t.summation = ["p","q", "i"]
fock_exp2t.coefficient = [TWOEL_INT_TRANS]
fock_exp2t.coefficient_idx.append(["p", "q", "i", "i"])
fock_exp2t.operator_idx.append(["p", "q"])
fock_exp2t.operator_type.append("s")
fock_exp2t.num_factor = -2.0

fock_exp3t = ugg()
fock_exp3t.summation = ["p","q", "i"]
fock_exp3t.coefficient = [TWOEL_INT_TRANS]
fock_exp3t.coefficient_idx.append(["p", "i", "i","q"])
fock_exp3t.operator_idx.append(["p", "q"])
fock_exp3t.operator_type.append("s")
fock_exp3t.num_factor = 1.0

fc_1 = ugg()
fc_1.summation = ["p","q","r","s"]
fc_1.coefficient = [COREL]
fc_1.coefficient_idx.append(["p","q","r","s"])
fc_1.operator_idx.append(["p", "q"])
fc_1.operator_idx.append(["r", "s"])
fc_1.operator_type.append("s")
fc_1.operator_type.append("s")
fc_1.num_factor = 1./2.

fc_2 = ugg()
fc_2.summation = ["p","q","r","s"]
fc_2.coefficient = [COREL]
fc_2.coefficient_idx.append(["p","q","r","s"])
fc_2.operator_idx.append(["p", "s"])
fc_2.operator_type.append("s")
fc_2.delta.append(["q","r"])
fc_2.num_factor = -1./2.


ham = arithmetic_string(h, g_1, g_2)
hamiltonian = arithmetic_string(h, g_1, g_2, fock_imp, fock_exp1, fock_exp2, fock_exp3)
hamiltonian_comp = arithmetic_string(hc, g_1c, g_2c, fock_impc, fock_exp1c, fock_exp2c, fock_exp3c)
hamiltoniant_comp = arithmetic_string(htc, g_1tc, g_2tc, fock_imptc, fock_exp1tc, fock_exp2tc, fock_exp3tc)
hamiltonian_comp_v = arithmetic_string(g_1c, g_2c,)
hamiltonian_comp_fock = arithmetic_string(fock_impc)
#
# T1-transformed full Hamiltonian: Exp(-T1) H Exp(T1) where H = sum_pq h_pq + sum_pqrs g_pqrs(E_pqE_rs - delta_qrE_ps)
# Keep in mind that in a lot of nordic publications operator F = h_pq + sum_k(2g_kkpq - g_kqpk)
#
hamiltoniant = arithmetic_string(ht, g_1t, g_2t)#, fock_impt, fock_exp1t, fock_exp2t, fock_exp3t)
hamiltoniantc = arithmetic_string(htc, g_1tc, g_2tc)
#hamiltoniant_comp = arithmetic_string(htc, g_1tc, g_2tc, fock_imptc, fock_excp1tc, fock_exp2tc, fock_exp3tc)

flukt_potential = arithmetic_string(g_1, g_2, fock_exp2, fock_exp3)
flukt_potentialc = arithmetic_string(g_1c, g_2c, fock_exp2c, fock_exp3c) 
flukt_potentialt = arithmetic_string(g_1t, g_2t)#, fock_exp2t, fock_exp3t) 
fock = arithmetic_string(fock_imp)

g12 = arithmetic_string(g_1, g_2)
#
# Coupled cluster T1 amplitude
#
t1 = ugg()
t1.summation = ["a","i"]
t1.coefficient = [CC_AMPLITUDE]
t1.coefficient_idx.append(["a","i"])
t1.operator_idx.append(["a", "i"])
t1.operator_type.append("s")

t1c = ugg()
t1c.summation = ["c","k"]
t1c.coefficient = [CC_AMPLITUDE]
t1c.coefficient_idx.append(["c","k"])
t1c.operator_idx.append(["k", "c"])
t1c.operator_type.append("s")

#
# Coupled cluster T2 amplitude
#
t2 = ugg()
t2.summation = ["a","i","b","j"]
t2.coefficient = [CC_AMPLITUDE]
t2.coefficient_idx.append(["a","i", "b", "j"])
t2.operator_idx.append(["a", "i"])
t2.operator_idx.append(["b", "j"])
t2.operator_type.append("s")
t2.operator_type.append("s")
t2.num_factor = 1./2.


t2ef = ugg()
#t2ef.summation = ["e","m","f","n"]
t2ef.coefficient = [CC_AMPLITUDE]
t2ef.coefficient_idx.append(["e","m", "f", "n"])
t2ef.operator_idx.append(["e", "m"])
t2ef.operator_idx.append(["f", "n"])
t2ef.operator_type.append("s")
t2ef.operator_type.append("s")
t2ef.num_factor = 1./2.



# t2tt = ugg()
# t2tt.summation = ["α","i","β","j", "k", "l"]
# t2tt.coefficient = [TWOEL_INT]
# t2tt.coefficient_idx.append(["α","k", "β", "l"])
# t2tt.operator_idx.append(["α", "i"])
# t2tt.operator_idx.append(["β", "j"])
# t2tt.operator_type.append("s")
# t2tt.operator_type.append("s")
# t2tt.num_factor = 1./2.




t2c = ugg()
t2c.summation = ["a","i","b","j"]
t2c.coefficient = [CC_AMPLITUDE]
t2c.coefficient_idx.append(["a","i", "b", "j"])
t2c.operator_idx.append(["i", "a"])
t2c.operator_idx.append(["j", "b"])
t2c.operator_type.append("s")
t2c.operator_type.append("s")
t2c.num_factor = 1./2.

#
# Coupled cluster T3 amplitude
#
t3 = ugg()
t3.summation = ["a","i","b","j", "c", "k"]
t3.coefficient = ["t"]
t3.coefficient_idx.append(["a","i", "b", "j", "c", "k"])
t3.operator_idx.append(["a", "i"])
t3.operator_idx.append(["b", "j"])
t3.operator_idx.append(["c", "k"])
t3.operator_type.append("s")
t3.operator_type.append("s")
t3.operator_type.append("s")
t3.num_factor = 1.0/6.

t3c = ugg()
t3c.summation = ["a","i","b","j", "c", "k"]
t3c.coefficient = ["t"]
t3c.coefficient_idx.append(["a","i", "b", "j", "c", "k"])
t3c.operator_idx.append(["i", "a"])
t3c.operator_idx.append(["j", "b"])
t3c.operator_idx.append(["k", "c"])
t3c.operator_type.append("s")
t3c.operator_type.append("s")
t3c.operator_type.append("s")
t3c.num_factor = 1.0/6.



#                                                                                                            
# Coupled cluster T2-F12 amplitude                                                                  
#                                      
t2fa = ugg()
t2fa.summation = ["k","i","l","j", "α", "β"]
t2fa.coefficient = [F12_AMPLITUDE]
t2fa.coefficient.append(F12_TWOEL)
t2fa.coefficient_idx.append(["k","i", "l", "j"])
t2fa.coefficient_idx.append(["α","k", "β", "l"])
t2fa.operator_idx.append(["α", "i"])
t2fa.operator_idx.append(["β", "j"])
t2fa.operator_type.append("s")
t2fa.operator_type.append("s")
t2fa.num_factor = 1./2.

t2fac = ugg()
t2fac.summation = ["k","i","l","j", "α", "β"]
t2fac.coefficient = [F12_AMPLITUDE]
t2fac.coefficient.append(F12_TWOEL)
t2fac.coefficient_idx.append(["k","i", "l", "j"])
t2fac.coefficient_idx.append(["α","k", "β", "l"])
t2fac.operator_idx.append(["i", "α"])
t2fac.operator_idx.append(["j","β"])
t2fac.operator_type.append("s")
t2fac.operator_type.append("s")
t2fac.num_factor = 1./2.

s1fa = ugg()
s1fa.summation = ["α", "i"]
s1fa.coefficient = [S_AMPLITUDE]
s1fa.coefficient_idx.append(["α","i"])
s1fa.operator_idx.append(["α", "i"])
s1fa.operator_type.append("s")

s1fac = ugg()
s1fac.summation = ["α", "i"]
s1fac.coefficient = [S_AMPLITUDE]
s1fac.coefficient_idx.append(["α","i"])
s1fac.operator_idx.append(["i","α"])
s1fac.operator_type.append("s")


s_v1 = ugg()
s_v1.summation = ['a', 'b', 'i', 'j']
s_v1.coefficient = [CC_AMPLITUDE]
s_v1.coefficient.append(CC_AMPLITUDE)
s_v1.coefficient_idx.append(["a","j", "b", "i"])
s_v1.coefficient_idx.append(["b", "j"])
s_v1.operator_idx.append(["i","a"])
s_v1.operator_type.append("s")
s_v1.num_factor = -1.0


s_v2 = ugg()
s_v2.summation = ['a', 'b', 'i', 'j']
s_v2.coefficient = [CC_AMPLITUDE]
s_v2.coefficient.append(CC_AMPLITUDE)
s_v2.coefficient_idx.append(["a","i", "b", "j"])
s_v2.coefficient_idx.append(["b", "j"])
s_v2.operator_idx.append(["i","a"])
s_v2.operator_type.append("s")
s_v2.num_factor = 2.0

s_c1 = ugg()
s_c1.summation = ["α", 'a', 'k', 'i', 'j', 'l']
s_c1.coefficient = [F12_AMPLITUDE]
s_c1.coefficient.append(F12_TWOEL)
s_c1.coefficient.append(CC_AMPLITUDE)
s_c1.coefficient_idx.append(["j","i", "k", "l"])
s_c1.coefficient_idx.append(["a","j", "α", "k"])
s_c1.coefficient_idx.append(["a", "l"])
s_c1.operator_idx.append(["i","α"])
s_c1.operator_type.append("s")
s_c1.num_factor = -1.0

s_c2 = ugg()
s_c2.summation = ["α", 'a', 'k', 'i', 'j', 'l']
s_c2.coefficient = [F12_AMPLITUDE]
s_c2.coefficient.append(F12_TWOEL)
s_c2.coefficient.append(CC_AMPLITUDE)
s_c2.coefficient_idx.append(["j","i", "k", "l"])
s_c2.coefficient_idx.append(["a","k", "α", "j"])
s_c2.coefficient_idx.append(["a", "l"])
s_c2.operator_idx.append(["i","α"])
s_c2.operator_type.append("s")
s_c2.num_factor = 2.0


s_A1 = ugg()
s_A1.summation = ["α", 'a', 'k', 'i', 'j', 'l']
s_A1.coefficient = [F12_AMPLITUDE]
s_A1.coefficient.append(F12_TWOEL)
s_A1.coefficient.append(CC_AMPLITUDE)
s_A1.coefficient_idx.append(["j","i", "k", "l"])
s_A1.coefficient_idx.append(["a","j", "A", "k"])
s_A1.coefficient_idx.append(["a", "l"])
s_A1.operator_idx.append(["i","A"])
s_A1.operator_type.append("s")
s_A1.num_factor = -1.0

s_A2 = ugg()
s_A2.summation = ["α", 'a', 'k', 'i', 'j', 'l']
s_A2.coefficient = [F12_AMPLITUDE]
s_A2.coefficient.append(F12_TWOEL)
s_A2.coefficient.append(CC_AMPLITUDE)
s_A2.coefficient_idx.append(["j","i", "k", "l"])
s_A2.coefficient_idx.append(["a","k", "A", "j"])
s_A2.coefficient_idx.append(["a", "l"])
s_A2.operator_idx.append(["i","A"])
s_A2.operator_type.append("s")
s_A2.num_factor = 2.0





s2fa = ugg()
s2fa.summation = ["k","i","l","j", "α", "β"]
s2fa.coefficient = [F12_SAMPLITUDE]
# s2fa.coefficient.append(F12_TWOEL)
# s2fa.coefficient_idx.append(["k","i", "l", "j"])
s2fa.coefficient_idx.append(["α","i", "β", "j"])
s2fa.operator_idx.append(["α", "i"])
s2fa.operator_idx.append(["β", "j"])
s2fa.operator_type.append("s")
s2fa.operator_type.append("s")
s2fa.num_factor = 1./2.


gemi = ugg()
gemi.summation = ["α", "β"]
gemi.coefficient.append(F12_TWOEL)
gemi.coefficient_idx.append(["α","k", "β", "l"])
gemi.operator_idx.append(["j", "β"])
gemi.operator_idx.append(["i", "α"])
gemi.operator_type.append("s")
gemi.operator_type.append("s")
gemi.num_factor = 1./2.


         
# t2fa = ugg()
# t2fa.summation = ["k","i","l","j", "a1", "a2"]
# t2fa.coefficient = [CC_AMPLITUDE]
# t2fa.coefficient.append(TWOEL_INT_TRANS)
# t2fa.coefficient_idx.append(["k","i", "l", "j"])
# t2fa.coefficient_idx.append(["a1","a2", "k", "l"])
# t2fa.operator_idx.append(["a1", "i"])
# t2fa.operator_idx.append(["a2", "j"])
# t2fa.operator_type.append("s")
# t2fa.operator_type.append("s")
# t2fa.num_factor = 1./2.


t2fb = ugg()
t2fb.summation = ["k","i","l","j", "a", "b"]
t2fb.coefficient = [F12_AMPLITUDE]
t2fb.coefficient.append(F12_TWOEL)
t2fb.coefficient_idx.append(["k","i", "l", "j"])
t2fb.coefficient_idx.append(["a","b", "k", "l"])
t2fb.operator_idx.append(["a", "i"])
t2fb.operator_idx.append(["b", "j"])
t2fb.operator_type.append("s")
t2fb.operator_type.append("s")
t2fb.num_factor = -1./2.



t2c = ugg()
t2c.summation = ["a","i","b","j"]
t2c.coefficient = [CC_AMPLITUDE]
t2c.coefficient_idx.append(["a","i", "b", "j"])
t2c.operator_idx.append(["i", "a"])
t2c.operator_idx.append(["j", "b"])
t2c.operator_type.append("s")
t2c.operator_type.append("s")
t2c.num_factor = 1./2.


#omega
omega1 = ugg()
omega1.summation = ["a","i"]
omega1.coefficient = ["Og"]
omega1.coefficient_idx.append(["a","i"])
omega1.operator_idx.append(["a", "i"])
omega1.operator_type.append("s")

omega2 = ugg()
omega2.summation = ["a","i","b","j"]
omega2.coefficient = ["Og"]
omega2.coefficient_idx.append(["a","i", "b", "j"])
omega2.operator_idx.append(["a", "i"])
omega2.operator_idx.append(["b", "j"])
omega2.operator_type.append("s")
omega2.operator_type.append("s")
omega2.num_factor = 1./2.

omega3 = ugg()
omega3.summation = ["a", "i", "b", "j", "c", "k"]
omega3.coefficient = ["Og"]
omega3.coefficient_idx.append(["a","i", "b", "j", "c", "k"])
omega3.operator_idx.append(["a", "i"])
omega3.operator_idx.append(["b", "j"])
omega3.operator_idx.append(["c", "k"])
omega3.operator_type.append("s")
omega3.operator_type.append("s")
omega3.operator_type.append("s")
omega3.num_factor = 1.0/6.


s1 = ugg()
s1.summation = ["a","i"]
s1.coefficient = [S_AMPLITUDE]
s1.coefficient_idx.append(["a","i"])
s1.operator_idx.append(["a", "i"])
s1.operator_type.append("s")

s1A = ugg()
s1A.summation = ["A","i"]
s1A.coefficient = [S_AMPLITUDE]
s1A.coefficient_idx.append(["A","i"])
s1A.operator_idx.append(["A", "i"])
s1A.operator_type.append("s")

s1c = ugg()
s1c.summation = ["c","k"]
s1c.coefficient = [S_AMPLITUDE]
s1c.coefficient_idx.append(["c","k"])
s1c.operator_idx.append(["k", "c"])
s1c.operator_type.append("s")

s1Ac = ugg()
s1Ac.summation = ["A","i"]
s1Ac.coefficient = [S_AMPLITUDE]
s1Ac.coefficient_idx.append(["A","i"])
s1Ac.operator_idx.append(["i", "A"])
s1Ac.operator_type.append("s")


s1cc = ugg()
s1cc.summation = ["a","i"]
s1cc.coefficient = [S_AMPLITUDE]
s1cc.coefficient_idx.append(["a","i"])
s1cc.operator_idx.append(["i", "a"])
s1cc.operator_type.append("s")



#
# Coupled cluster T2 amplitude
#
s2 = ugg()
s2.summation = ["a","i","b","j"]
s2.coefficient = [S_AMPLITUDE]
s2.coefficient_idx.append(["a","i", "b", "j"])
s2.operator_idx.append(["a", "i"])
s2.operator_idx.append(["b", "j"])
s2.operator_type.append("s")
s2.operator_type.append("s")
s2.num_factor = 1./2.


s2Aa = ugg()
s2Aa.summation = ["A","i","a","j"]
s2Aa.coefficient = [S_AMPLITUDE]
s2Aa.coefficient_idx.append(["A","i", "a", "j"])
s2Aa.operator_idx.append(["A", "i"])
s2Aa.operator_idx.append(["a", "j"])
s2Aa.operator_type.append("s")
s2Aa.operator_type.append("s")
s2Aa.num_factor = 1./2.

s2AB = ugg()
s2AB.summation = ["A","i","B","j"]
s2AB.coefficient = [S_AMPLITUDE]
s2AB.coefficient_idx.append(["A","i", "B", "j"])
s2AB.operator_idx.append(["A", "i"])
s2AB.operator_idx.append(["B", "j"])
s2AB.operator_type.append("s")
s2AB.operator_type.append("s")
s2AB.num_factor = 1./2.

s2bB = ugg()
s2bB.summation = ["a","i","A","j"]
s2bB.coefficient = [S_AMPLITUDE]
s2bB.coefficient_idx.append(["a","i", "A", "j"])
s2bB.operator_idx.append(["a", "i"])
s2bB.operator_idx.append(["A", "j"])
s2bB.operator_type.append("s")
s2bB.operator_type.append("s")
s2bB.num_factor = 1./2.


s2c = ugg()
s2c.summation = ["a","i","b","j"]
s2c.coefficient = [S_AMPLITUDE]
s2c.coefficient_idx.append(["a","i", "b", "j"])
s2c.operator_idx.append(["i", "a"])
s2c.operator_idx.append(["j", "b"])
s2c.operator_type.append("s")
s2c.operator_type.append("s")
s2c.num_factor = 1./2.


s2Aac = ugg()
s2Aac.summation = ["A","i","a","j"]
s2Aac.coefficient = [S_AMPLITUDE]
s2Aac.coefficient_idx.append(["A","i", "a", "j"])
s2Aac.operator_idx.append(["i", "A"])
s2Aac.operator_idx.append(["j", "a"])
# s2Aac.operator_idx.append(["A", "i"])
# s2Aac.operator_idx.append(["a", "j"])
s2Aac.operator_type.append("s")
s2Aac.operator_type.append("s")
s2Aac.num_factor = 1./2.

s2ABc = ugg()
s2ABc.summation = ["A","i","B","j"]
s2ABc.coefficient = [S_AMPLITUDE]
s2ABc.coefficient_idx.append(["A","i", "B", "j"])
s2ABc.operator_idx.append(["i", "A"])
s2ABc.operator_idx.append(["j", "B"])

# s2ABc.operator_idx.append(["A", "i"])
# s2ABc.operator_idx.append(["B", "j"])
s2ABc.operator_type.append("s")
s2ABc.operator_type.append("s")
s2ABc.num_factor = 1./2.

s2bBc = ugg()
s2bBc.summation = ["a","i","A","j"]
s2bBc.coefficient = [S_AMPLITUDE]
s2bBc.coefficient_idx.append(["a","i", "A", "j"])
s2bBc.operator_idx.append(["i", "a"])
s2bBc.operator_idx.append(["j", "A"])
# s2bBc.operator_idx.append(["a", "i"])
# s2bBc.operator_idx.append(["A", "j"])
s2bBc.operator_type.append("s")
s2bBc.operator_type.append("s")
s2bBc.num_factor = 1./2.


#
# Coupled cluster T3 amplitude
#
s3 = ugg()
s3.summation = ["a","i","b","j", "c", "k"]
s3.coefficient = [S_AMPLITUDE]
s3.coefficient_idx.append(["a","i", "b", "j", "c", "k"])
s3.operator_idx.append(["a", "i"])
s3.operator_idx.append(["b", "j"])
s3.operator_idx.append(["c", "k"])
s3.operator_type.append("s")
s3.operator_type.append("s")
s3.operator_type.append("s")
s3.num_factor = 1.0/6.

s3c = ugg()
s3c.summation = ["a","i","b","j", "c", "k"]
s3c.coefficient = [S_AMPLITUDE]
s3c.coefficient_idx.append(["a","i", "b", "j", "c", "k"])
s3c.operator_idx.append(["i", "a"])
s3c.operator_idx.append(["j", "b"])
s3c.operator_idx.append(["k", "c"])
s3c.operator_type.append("s")
s3c.operator_type.append("s")
s3c.operator_type.append("s")
s3c.num_factor = 1.0/6.



e1cs = ugg()
e1cs.summation = ["a","i"]
e1cs.operator_idx.append(["i", "a"])
e1cs.operator_type.append("s")

e2csa = ugg()
e2csa.summation = ["a","i", "b", "j"]
e2csa.operator_idx.append(["i", "a"])
e2csa.operator_idx.append(["j", "b"])
e2csa.operator_type.append("s")
e2csa.operator_type.append("s")

e2csb = ugg()
e2csb.summation = ["a","i", "b", "j"]
e2csb.operator_idx.append(["j", "a"])
e2csb.operator_idx.append(["i", "b"])
e2csb.operator_type.append("s")
e2csb.operator_type.append("s")

e3cs = ugg()
e3cs.summation = ["a","i", "b", "j", "c", "k"]
e3cs.operator_idx.append(["i", "a"])
e3cs.operator_idx.append(["j", "b"])
e3cs.operator_idx.append(["k", "c"])
e3cs.operator_type.append("s")
e3cs.operator_type.append("s")
e3cs.operator_type.append("s")



#
# Coupled cluster T3 amplitude - linear combination
#
t3p1 = ugg()
t3p1.summation = ["a","i","b","j", "c", "k"]
t3p1.coefficient = ["t"]
t3p1.coefficient_idx.append(["a","i", "b", "j", "c", "k"])
t3p1.operator_idx.append(["a", "i"])
t3p1.operator_idx.append(["b", "j"])
t3p1.operator_idx.append(["c", "k"])
t3p1.operator_type.append("s")
t3p1.operator_type.append("s")
t3p1.operator_type.append("s")
t3p1.num_factor = -1.0/6.


t3p2 = ugg()
t3p2.summation = ["a","i","b","j", "c", "k"]
t3p2.coefficient = ["t"]
t3p2.coefficient_idx.append(["a","i", "b", "j", "c", "k"])
t3p2.operator_idx.append(["a", "j"])
t3p2.operator_idx.append(["b", "i"])
t3p2.operator_idx.append(["c", "k"])
t3p2.operator_type.append("s")
t3p2.operator_type.append("s")
t3p2.operator_type.append("s")
t3p2.num_factor = 1.0/6.

#
# Hermitian-conjutage of the T1 amplitude operator
#
s1c = ugg()
s1c.summation = ["c","k"]
s1c.coefficient = [S_AMPLITUDE]
s1c.coefficient_idx.append(["c","k"])
s1c.operator_idx.append(["k", "c"])
s1c.operator_type.append("s")
#
# Hermitian-conjugate of the T2 amplitude operator
#
s2c = ugg()
s2c.summation = ["a","i","b","j"]
s2c.coefficient = [S_AMPLITUDE]
s2c.coefficient_idx.append(["a","i", "b", "j"])
s2c.operator_idx.append(["i", "a"])
s2c.operator_idx.append(["j", "b"])
s2c.operator_type.append("s")
s2c.operator_type.append("s")
s2c.num_factor = 1./2.

eab = ugg()
eab.operator_idx.append(["a", "b"])
eab.operator_type.append("s")

mu1 = ugg()
mu1.operator_idx.append(['a', 'i'])
mu1.operator_type.append("s")

mu1_triplet = ugg()
mu1_triplet.operator_idx.append(['a', 'i'])
mu1_triplet.operator_type.append("t0")

mu1p = ugg()
mu1p.operator_idx.append(['c', 'k'])
mu1p.operator_type.append("s")

mu2 = ugg()
mu2.operator_idx.append(['a', 'i'])
mu2.operator_idx.append(['b', 'j'])
mu2.operator_type.append("s")
mu2.operator_type.append("s")

mu2_triplet = ugg()
mu2_triplet.operator_idx.append(['a', 'i'])
mu2_triplet.operator_idx.append(['b', 'j'])
mu2_triplet.operator_type.append("t0")
mu2_triplet.operator_type.append("s")

mu3_triplet = ugg()
mu3_triplet.operator_idx.append(['a', 'i'])
mu3_triplet.operator_idx.append(['b', 'j'])
mu3_triplet.operator_idx.append(['c', 'k'])
mu3_triplet.operator_type.append("s")
mu3_triplet.operator_type.append("t0")
mu3_triplet.operator_type.append("s")


mu2p = ugg()
mu2p.operator_idx.append(['c', 'k'])
mu2p.operator_idx.append(['d', 'l'])
mu2p.operator_type.append("s")
mu2p.operator_type.append("s")


mu3 = ugg()
mu3.operator_idx.append(['a', 'i'])
mu3.operator_idx.append(['b', 'j'])
mu3.operator_idx.append(['c', 'k'])
mu3.operator_type.append("s")
mu3.operator_type.append("s")
mu3.operator_type.append("s")


nukl = ugg()
nukl.operator_idx.append(["k", "l"])
nukl.operator_type.append("s")


nubj = ugg()
nubj.operator_idx.append(["b", "j"])
nubj.operator_type.append("s")

nubj_triplet = ugg()
nubj_triplet.operator_idx.append(["b", "j"])
nubj_triplet.operator_type.append("t0")

nuck = ugg()
nuck.operator_idx.append(["c", "k"])
nuck.operator_type.append("s")

nudl = ugg()
nudl.operator_idx.append(["d", "l"])
nudl.operator_type.append("s")

nuai = ugg()
nuai.operator_idx.append(["a", "i"])
nuai.operator_type.append("s")

nuaiai = ugg()
nuaiai.operator_idx.append(["a", "i"])
nuaiai.operator_idx.append(["a", "i"])
nuaiai.operator_type.append("s")
nuaiai.operator_type.append("s")

nubjck = ugg()
nubjck.operator_idx.append(["b", "j"])
nubjck.operator_idx.append(["c", "k"])
nubjck.operator_type.append("s")
nubjck.operator_type.append("s")

nuckdl = ugg()
nuckdl.operator_idx.append(["c", "k"])
nuckdl.operator_idx.append(["d", "l"])
nuckdl.operator_type.append("s")
nuckdl.operator_type.append("s")

nudlem = ugg()
nudlem.operator_idx.append(["d", "l"])
nudlem.operator_idx.append(["e", "m"])
nudlem.operator_type.append("s")
nudlem.operator_type.append("s")

nuaiam = ugg()
nuaiam.operator_idx.append(["a", "l"])
nuaiam.operator_idx.append(["b", "i"])
nuaiam.operator_type.append("s")
nuaiam.operator_type.append("s")

nubjckdl = ugg()
nubjckdl.operator_idx.append(["b", "j"])
nubjckdl.operator_idx.append(["c", "k"])
nubjckdl.operator_idx.append(["d", "l"])
nubjckdl.operator_type.append("s")
nubjckdl.operator_type.append("s")
nubjckdl.operator_type.append("s")

nubjbkdj = ugg()
nubjbkdj.operator_idx.append(["b", "j"])
nubjbkdj.operator_idx.append(["b", "k"])
nubjbkdj.operator_idx.append(["d", "j"])
nubjbkdj.operator_type.append("s")
nubjbkdj.operator_type.append("s")
nubjbkdj.operator_type.append("s")

nubjckdj = ugg()
nubjckdj.operator_idx.append(["b", "j"])
nubjckdj.operator_idx.append(["c", "k"])
nubjckdj.operator_idx.append(["d", "j"])
nubjckdj.operator_type.append("s")
nubjckdj.operator_type.append("s")
nubjckdj.operator_type.append("s")

nubjckdk = ugg()
nubjckdk.operator_idx.append(["b", "j"])
nubjckdk.operator_idx.append(["c", "k"])
nubjckdk.operator_idx.append(["d", "k"])
nubjckdk.operator_type.append("s")
nubjckdk.operator_type.append("s")
nubjckdk.operator_type.append("s")




nuaibjck = ugg()
nuaibjck.operator_idx.append(["a", "i"])
nuaibjck.operator_idx.append(["b", "j"])
nuaibjck.operator_idx.append(["c", "k"])
nuaibjck.operator_type.append("s")
nuaibjck.operator_type.append("s")
nuaibjck.operator_type.append("s")


nuckdlem = ugg()
nuckdlem.operator_idx.append(["c", "k"])
nuckdlem.operator_idx.append(["d", "l"])
nuckdlem.operator_idx.append(["e", "m"])
nuckdlem.operator_type.append("s")
nuckdlem.operator_type.append("s")
nuckdlem.operator_type.append("s")

nuckclek = ugg()
nuckclek.operator_idx.append(["c", "k"])
nuckclek.operator_idx.append(["c", "l"])
nuckclek.operator_idx.append(["e", "k"])
nuckclek.operator_type.append("s")
nuckclek.operator_type.append("s")
nuckclek.operator_type.append("s")


nuckckck = ugg()
nuckckck.operator_idx.append(["c", "k"])
nuckckck.operator_idx.append(["c", "k"])
nuckckck.operator_idx.append(["c", "k"])
nuckckck.operator_type.append("s")
nuckckck.operator_type.append("s")
nuckckck.operator_type.append("s")



nudlemfn = ugg()
nudlemfn.operator_idx.append(["d", "l"])
nudlemfn.operator_idx.append(["e", "m"])
nudlemfn.operator_idx.append(["f", "n"])
nudlemfn.operator_type.append("s")
nudlemfn.operator_type.append("s")
nudlemfn.operator_type.append("s")


nudlemfn1 = ugg()
nudlemfn1.operator_idx.append(["a", "i"])
nudlemfn1.operator_idx.append(["b", "j"])
nudlemfn1.operator_idx.append(["c", "k"])
nudlemfn1.operator_type.append("s")
nudlemfn1.operator_type.append("s")
nudlemfn1.operator_type.append("s")
#nudlemfn1.num_factor = (-1.0)

nudlemfn2 = ugg()
nudlemfn2.operator_idx.append(["a", "k"])
nudlemfn2.operator_idx.append(["b", "i"])
nudlemfn2.operator_idx.append(["c", "j"])
nudlemfn2.operator_type.append("s")
nudlemfn2.operator_type.append("s")
nudlemfn2.operator_type.append("s")
nudlemfn2.num_factor = (-1.0)


es = ugg()
es.summation = ['a','i']
es.coefficient = [OBSERVABLE_X]
es.coefficient_idx.append(['i','a'])
es.operator_idx.append(['i','a'])
es.operator_type.append("s")

eai_comp = ugg()
eai_comp.summation = ["α","i","β","j"]
eai_comp.operator_idx.append(["i", "α"])
eai_comp.operator_idx.append(["j", "β" ])
eai_comp.operator_type.append("s")
eai_comp.operator_type.append("s")
eai_comp.num_factor = 1./2.

ejbias = ugg()
ejbias.summation = ['a', 'b', 'i', 'j']
ejbias.operator_idx.append(["j","b"])
ejbias.operator_idx.append(["i","a"])
ejbias.operator_type.append("s")
ejbias.operator_type.append("s")

eaibjs = ugg()
eaibjs.summation = ['a', 'b', 'i', 'j']
eaibjs.operator_idx.append(["a","i"])
eaibjs.operator_idx.append(["b","j"])
eaibjs.operator_type.append("s")
eaibjs.operator_type.append("s")


eaibjcks = ugg()
eaibjcks.summation = ['d', 'e', 'f', 'l', 'm', 'n']
eaibjcks.operator_idx.append(["d","l"])
eaibjcks.operator_idx.append(["e","m"])
eaibjcks.operator_idx.append(["f","n"])
eaibjcks.operator_type.append("s")
eaibjcks.operator_type.append("s")
eaibjcks.operator_type.append("s")


eai = ugg()
eai.operator_idx.append(["a","i"])
eai.operator_type.append("s")

eia = ugg()
eia.operator_idx.append(["i","a"])
eia.operator_type.append("s")

eiabj = ugg()                                                                                                               
eiabj.operator_idx.append(["i","a"])                                                                                           
eiabj.operator_idx.append(["b","j"])                                                                                            
eiabj.operator_type.append("s")                                                                                                            
eiabj.operator_type.append("s")


eiajb = ugg()
eiajb.operator_idx.append(["i","a"])
eiajb.operator_idx.append(["j","b"])
eiajb.operator_type.append("s")
eiajb.operator_type.append("s")

ejaib = ugg()
ejaib.operator_idx.append(["j","a"])
ejaib.operator_idx.append(["i","b"])
ejaib.operator_type.append("s")
ejaib.operator_type.append("s")


eaibj = ugg()
eaibj.operator_idx.append(["a","i"])
eaibj.operator_idx.append(["b","j"])
eaibj.operator_type.append("s")
eaibj.operator_type.append("s")

ejbia = ugg()
ejbia.operator_idx.append(["j","b"])
ejbia.operator_idx.append(["i","a"])
ejbia.operator_type.append("s")
ejbia.operator_type.append("s")

eibja = ugg()
eibja.operator_idx.append(["i","b"])
eibja.operator_idx.append(["j","a"])
eibja.operator_type.append("s")
eibja.operator_type.append("s")


ebjai = ugg()
ebjai.operator_idx.append(["b","j"])
ebjai.operator_idx.append(["a","i"])
ebjai.operator_type.append("s")
ebjai.operator_type.append("s")

ebiaj = ugg()
ebiaj.operator_idx.append(["b","i"])
ebiaj.operator_idx.append(["a","j"])
ebiaj.operator_type.append("s")
ebiaj.operator_type.append("s")


eckdls = ugg()
eckdls.summation = ['c','k','d','l']
eckdls.operator_idx.append(["k","c"])
eckdls.operator_idx.append(["l","d"])
eckdls.operator_type.append("s")
eckdls.operator_type.append("s")
eckdls.num_factor = 1.0/4.0

ekclds = ugg()
#ekclds.summation = ['c','k','d','l']
ekclds.operator_idx.append(["c","k"])
ekclds.operator_idx.append(["d","l"])
ekclds.operator_type.append("s")
ekclds.operator_type.append("s")
ekclds.num_factor = 1.0/2.0

e3s = ugg()
e3s.operator_idx.append(["a","i"])
e3s.operator_idx.append(["b","j"])
e3s.operator_idx.append(["c","k"])
e3s.operator_type.append("s")
e3s.operator_type.append("s")
e3s.operator_type.append("s")
e3s.num_factor = 1.0/6.0


sw1 = ugg()
sw1.operator_idx.append(["l","d"])
sw1.operator_idx.append(["k","c"])
sw1.operator_type.append("s")
sw1.operator_type.append("s")
sw1.num_factor = 1./3.

sw2 = ugg()
sw2.operator_idx.append(["k","d"])
sw2.operator_idx.append(["l","c"])
sw2.operator_type.append("s")
sw2.operator_type.append("s")
sw2.num_factor = 1./6.



eiajbs = ugg()
#eiajbs.summation = ['a', 'i', 'b', 'j']
eiajbs.coefficient.append('Y')
eiajbs.coefficient_idx.append(['a', 'b', 'i', 'j'])
eiajbs.operator_idx.append(["i","a"])
eiajbs.operator_idx.append(["j","b"])
eiajbs.operator_type.append("s")
eiajbs.operator_type.append("s")

eiajbkcs = ugg()
eiajbkcs.coefficient.append('Y')
eiajbkcs.coefficient_idx.append(['a>', 'b>', 'c>', 'i>', 'j>', 'k>'])
eiajbkcs.operator_idx.append(["i>","a>"])
eiajbkcs.operator_idx.append(["j>","b>"])
eiajbkcs.operator_idx.append(["k>","c>"])
eiajbkcs.operator_type.append("s")
eiajbkcs.operator_type.append("s")
eiajbkcs.operator_type.append("s")


eckdls2 = ugg()
eckdls2.summation = ['c','k','d','l']
eckdls2.operator_idx.append(["l","c"])
eckdls2.operator_idx.append(["k","d"])
eckdls2.operator_type.append("s")
eckdls2.operator_type.append("s")


eckdl = ugg()
eckdl.operator_idx.append(["c","k"])
eckdl.operator_idx.append(["d","l"])
eckdl.operator_type.append("s")
eckdl.operator_type.append("s")

ekcld = ugg()
ekcld.operator_idx.append(["k","c"])
ekcld.operator_idx.append(["l","d"])
ekcld.operator_type.append("s")
ekcld.operator_type.append("s")


ejb = ugg()
ejb.operator_idx.append(["j","b"])
ejb.operator_type.append("s")


ebj = ugg()
ebj.operator_idx.append(["b","j"])
ebj.operator_type.append("s")

ecdkl_1 = ugg()
ecdkl_1.operator_idx.append(["c", "k"])
ecdkl_1.operator_idx.append(["d", "l"])
ecdkl_1.operator_type.append("s")
ecdkl_1.operator_type.append("s")

tg_1 = ugg()
tg_1.summation = ["p","q","r","s"]
tg_1.coefficient = [TWOEL_INT]
tg_1.coefficient_idx.append(["r","s","p","q"])
tg_1.operator_idx.append(["p", "r"])
tg_1.operator_idx.append(["q", "s"])
tg_1.operator_type.append("s")
tg_1.operator_type.append("s")
tg_1.num_factor = 1./4.

tg_2 = ugg()
tg_2.summation = ["p","q","r","s"]
tg_2.coefficient = [TWOEL_INT]
tg_2.coefficient_idx.append(["r","s","p","q"])
tg_2.operator_idx.append(["p", "s"])
tg_2.operator_type.append("s")
tg_2.delta.append(["q","r"])
tg_2.num_factor = -1./4.

t1ai = ugg()
t1ai.coefficient = ["t"]
t1ai.coefficient_idx.append(["a","i"])
t1ai.operator_idx.append(["a", "i"])
t1ai.operator_type.append("s")


t2ai = ugg()
t2ai.coefficient = ["t"]
t2ai.coefficient_idx.append(["a","i", "a", "i"])
t2ai.operator_idx.append(["a", "i"])
t2ai.operator_idx.append(["a", "i"])
t2ai.operator_type.append("s")
t2ai.operator_type.append("s")
t2ai.num_factor = 1./2.

eaa = ugg()
eaa.operator_idx.append(["a","i"])
eaa.operator_idx.append(["a","i"])
eaa.operator_idx.append(["a","i"])
eaa.operator_type.append("s")
eaa.operator_type.append("s")
eaa.operator_type.append("s")


eck = ugg()
eck.operator_idx.append(["c","k"])
eck.operator_type.append("s")

eckr = ugg()
eckr.operator_idx.append(["c","k"])
eckr.operator_type.append("s")

ekcr = ugg()
ekcr.operator_idx.append(["k","c"])
ekcr.operator_type.append("s")

edl = ugg()
edl.operator_idx.append(["d","l"])
edl.operator_type.append("s")

eii = ugg()
eii.operator_idx.append(["i","i"])
eii.operator_type.append("s")

eij = ugg()
eij.operator_idx.append(["i","j"])
eij.operator_type.append("s")

eijkl = ugg()
eijkl.operator_idx.append(["i","j"])
eijkl.operator_type.append("s")
eijkl.operator_idx.append(["k","l"])
eijkl.operator_type.append("s")


kk = ugg()
kk.operator_idx.append(["d", "l"])
kk.operator_idx.append(["e", "m"])
kk.operator_idx.append(["f", "n"])
kk.operator_type.append("s")
kk.operator_type.append("s")
kk.operator_type.append("s")

epq3 = ugg()
epq3.summation = ["b", "j", "k"]
epq3.coefficient = ["t", "g"]
epq3.coefficient_idx.append(["a","j", "b", "k"])
epq3.coefficient_idx.append(["i","k", "b", "j"])

epq4 = ugg()
epq4.summation = ["b", "j", "k"]
epq4.coefficient = ["t", "g"]
epq4.coefficient_idx.append(["a","j", "b", "k"])
epq4.coefficient_idx.append(["i","j", "b", "k"])

eps = ugg()
#eps.summation = ["a", "b", "j"]
eps.operator_idx.append(["c","k"])
eps.operator_idx.append(["d","l"])
eps.operator_type.append("s")
eps.operator_type.append("s")
eps.num_factor = 2.08166817117e-17

w2 = ugg()
w2.summation = ["b", "c", 'j']
w2.coefficient = ["g", "t", "t"]
w2.coefficient_idx.append(["a","c", "b", "j"])
w2.coefficient_idx.append(["c","i"])
w2.coefficient_idx.append(["b","j"])
#w2.operator_idx.append(["p","k"])
#w2.operator_idx.append(["b", "i"])

u2 = ugg()
u2.summation = ["p"]
u2.coefficient = ["g"]
u2.coefficient_idx.append(["a","p", "b", "k"])
u2.operator_idx.append(["p","i"])
u2.operator_idx.append(["c", "j"])
u2.operator_type.append("s")
u2.operator_type.append("s")
u2.ugg_order = 2
u2.set_permutation([["a","i"],["b","j"]])


o1 = ugg()
o1.summation = ['p','q']
o1.operator_idx.append(['p','q'])
o1.operator_type.append("s")

o2 = ugg()
o2.summation = ['r','s']
o2.operator_idx.append(['r','s'])
o2.operator_type.append("s")

o3 = ugg()
o3.summation = ['t','u']
o3.operator_idx.append(['t','u'])
o3.operator_type.append("s")


b3d = ugg()
b3d.operator_idx.append(["b","l"])
b3d.operator_idx.append(["c","k"])
b3d.operator_idx.append(["d","j"])
b3d.operator_type.append("s")
b3d.operator_type.append("s")
b3d.operator_type.append("s")
b3d.num_factor = 1.0

b4d = ugg()
b4d.operator_idx.append(["b","k"])
b4d.operator_idx.append(["c","j"])
b4d.operator_idx.append(["d","l"])
b4d.operator_type.append("s")
b4d.operator_type.append("s")
b4d.operator_type.append("s")
b4d.num_factor = 1.0

p1ad = arithmetic_string(b3d)
p1bd = arithmetic_string(b4d)
p1d = p1ad.scale(1.0/sqrt(12)) + p1bd.scale(-1.0/sqrt(12))

b1 = ugg()
b1.operator_idx.append(["a","i"])
b1.operator_idx.append(["b","j"])
b1.operator_idx.append(["c","k"])
b1.operator_type.append("s")
b1.operator_type.append("s")
b1.operator_type.append("s")
b1.num_factor = 1.0

b2 = ugg()
b2.operator_idx.append(["a","j"])
b2.operator_idx.append(["b","k"])
b2.operator_idx.append(["c","i"])
b2.operator_type.append("s")
b2.operator_type.append("s")
b2.operator_type.append("s")
b2.num_factor = 1.0

b3 = ugg()
b3.operator_idx.append(["a","k"])
b3.operator_idx.append(["b","j"])
b3.operator_idx.append(["c","i"])
b3.operator_type.append("s")
b3.operator_type.append("s")
b3.operator_type.append("s")
b3.num_factor = 1.0

b4 = ugg()
b4.operator_idx.append(["a","j"])
b4.operator_idx.append(["b","i"])
b4.operator_idx.append(["c","k"])
b4.operator_type.append("s")
b4.operator_type.append("s")
b4.operator_type.append("s")
b4.num_factor = 1.0

b5 = ugg()
b5.operator_idx.append(["a","i"])
b5.operator_idx.append(["b","k"])
b5.operator_idx.append(["c","j"])
b5.operator_type.append("s")
b5.operator_type.append("s")
b5.operator_type.append("s")
b5.num_factor = 1.0

b6 = ugg()
b6.operator_idx.append(["a","k"])
b6.operator_idx.append(["b","i"])
b6.operator_idx.append(["c","j"])
b6.operator_type.append("s")
b6.operator_type.append("s")
b6.operator_type.append("s")
b6.num_factor = 1.0

p1 = ugg()
p1.operator_idx.append(["a","i"])
p1.operator_idx.append(["b","k"])
p1.operator_idx.append(["c","j"])
p1.operator_type.append("s")
p1.operator_type.append("s")
p1.operator_type.append("s")

p2 = ugg()
p2.operator_idx.append(["a","j"])
p2.operator_idx.append(["b","i"])
p2.operator_idx.append(["c","k"])
p2.operator_type.append("s")
p2.operator_type.append("s")
p2.operator_type.append("s")

p3 = ugg()
p3.operator_idx.append(["a","j"])
p3.operator_idx.append(["b","k"])
p3.operator_idx.append(["c","i"])
p3.operator_type.append("s")
p3.operator_type.append("s")
p3.operator_type.append("s")

p4 = ugg()
p4.operator_idx.append(["a","k"])
p4.operator_idx.append(["b","i"])
p4.operator_idx.append(["c","j"])
p4.operator_type.append("s")
p4.operator_type.append("s")
p4.operator_type.append("s")

p5 = ugg()
p5.operator_idx.append(["a","k"])
p5.operator_idx.append(["b","j"])
p5.operator_idx.append(["c","i"])
p5.operator_type.append("s")
p5.operator_type.append("s")
p5.operator_type.append("s")



n211a = arithmetic_string(b5, b6)
n211b = arithmetic_string(b1, b2, b3, b4)
n211 = n211a.scale(1.0/3) + n211b.scale(-1.0/6)

n212a = arithmetic_string(b1, b3)
n212b = arithmetic_string(b2, b4)
n212 = n212a.scale(sqrt(3.0)/6) + n212b.scale(-sqrt(3.0)/6)

n221a = arithmetic_string(b1, b4)
n221b = arithmetic_string(b2, b3) 
n221 = n221a.scale(sqrt(3.0)/6) + n221b.scale(-sqrt(3.0)/6)

n222a = arithmetic_string(b6)
n222b = arithmetic_string(b5)
n222c = arithmetic_string(b4, b3)
n222d = arithmetic_string(b1, b2)
n222 = n222a.scale(-1.0/3.) + n222b.scale(1.0/3.) + n222c.scale(-1.0/6.) + n222d.scale(1.0/6.)

n4a = arithmetic_string(b5, b4, b3)
n4b = arithmetic_string(b6, b2, b1)
n4 = n4a.scale(1.0/12) + n4b.scale(-1.0/12)


#------------------------------------------tg---------------

tga = ugg()
tga.summation = ["a", "i", "b", "j"]
tga.coefficient = [TWOEL_INT_AS]
tga.coefficient_idx.append(["i","j", "a", "b"])
tga.operator_idx.append(["a","j"])
tga.operator_idx.append(["b","i"])
tga.operator_type.append("s")
tga.operator_type.append("s")

tga.num_factor = -1.0

tgc = ugg()
tgc.summation = ["c", "k", "d", "l"]
tgc.coefficient = [TWOEL_INT_AS]
tgc.coefficient_idx.append(["k","l", "c", "d"])
tgc.operator_idx.append(["c","l"])
tgc.operator_idx.append(["d","k"])
tgc.operator_type.append("s")
tgc.operator_type.append("s")
tgc.num_factor = -1.0

Wg1 = ugg()
Wg1.summation = ["p", "q", "r", "s"]
Wg1.coefficient = [TWOEL_INT_AS]
Wg1.coefficient_idx.append(["p","q", "r", "s"])
Wg1.operator_idx.append(["r","p"])
Wg1.operator_type.append("s")
Wg1.delta.append(["q","s"])
Wg1.num_factor = 0.25

Wg2 = ugg()
Wg2.summation = ["p", "q", "r", "s"]
Wg2.coefficient = [TWOEL_INT_AS]
Wg2.coefficient_idx.append(["p","q", "r", "s"])
Wg2.operator_idx.append(["r","q"])
Wg2.operator_idx.append(["s","p"])
Wg2.operator_type.append("s")
Wg2.operator_type.append("s")
Wg2.num_factor = -0.25

Wg3 = ugg()
Wg3.summation = ["p", "q", "i"]
Wg3.coefficient = [TWOEL_INT_AS]
Wg3.coefficient_idx.append(["p","i", "q", "i"])
Wg3.operator_idx.append(["p","q"])
Wg3.operator_type.append("s")
Wg3.num_factor = -1.0

epq = ugg()
epq.summation = ["p", "q"]
epq.coefficient = [TEMP1]
epq.coefficient_idx.append(["p","q"])
epq.operator_idx.append(["p","q"])
epq.operator_type.append("s")

ekc = ugg()
ekc.summation = ["k", "c"]
ekc.coefficient = [TEMP1]
ekc.coefficient_idx.append(["c","k"])
ekc.operator_idx.append(["c","k"])
ekc.operator_type.append("s")

eld = ugg()
eld.summation = ["l", "d"]
eld.coefficient = [TEMP2]
eld.coefficient_idx.append(["d","l"])
eld.operator_idx.append(["d","l"])
eld.operator_type.append("s")

ers = ugg()
ers.summation = ["r", "s"]
ers.coefficient = [TEMP2]
ers.coefficient_idx.append(["s","r"])
ers.operator_idx.append(["r","s"])
ers.operator_type.append("s")

esr = ugg()
esr.summation = ["r", "s"]
esr.coefficient = [TEMP2]
esr.coefficient_idx.append(["r","s"])
esr.operator_idx.append(["s","r"])
esr.operator_type.append("s")

ers2 = ugg()
ers2.summation = ["r", "s"]
ers2.coefficient = [TEMP2]
ers2.coefficient_idx.append(["r","s"])
ers2.operator_idx.append(["r","s"])
ers2.operator_type.append("s")

eajbi = ugg()
eajbi.summation = ["a", "i", "b", "j"]
eajbi.coefficient = [TWOEL_INT_AS]
eajbi.coefficient_idx.append(["a","b", "i", "j"])
eajbi.operator_idx.append(["a","j"])
eajbi.operator_idx.append(["b","i"])
eajbi.operator_type.append("s")
eajbi.operator_type.append("s")
eajbi.num_factor = -1.0

eajbir = ugg()
eajbir.operator_idx.append(["a","j"])
eajbir.operator_idx.append(["b","i"])
eajbir.operator_type.append("s")
eajbir.operator_type.append("s")

eenfm = ugg()
eenfm.summation = ["e", "f", "m", "n"]
eenfm.coefficient = [TWOEL_INT_AS]
eenfm.coefficient_idx.append(["m","n", "e", "f"])
eenfm.operator_idx.append(["e","n"])
eenfm.operator_idx.append(["f","m"])
eenfm.operator_type.append("s")
eenfm.operator_type.append("s")
eenfm.num_factor = -1.0

ecldk = ugg()
ecldk.operator_idx.append(["c","l"])
ecldk.operator_idx.append(["d","k"])
ecldk.operator_type.append("s")
ecldk.operator_type.append("s")
ecldk.num_factor = -1.0

ekdlc = ugg()
ekdlc.operator_idx.append(["k","d"])
ekdlc.operator_idx.append(["l","c"])
ekdlc.operator_type.append("s")
ekdlc.operator_type.append("s")
ekdlc.num_factor = -1.0

wrpqs = ugg()
wrpqs.coefficient = [TWOEL_INT]
wrpqs.coefficient_idx.append(["r","p", 'q', 'q'])

wrsqp = ugg()
wrsqp.coefficient = [TWOEL_INT]
wrsqp.coefficient_idx.append(["r","s", 'p', 'q'])

plusz = ugg()
plusz.summation = ["a", "i", "b", "j", "c", "k"]
plusz.coefficient = [EOM_TRIPLET_R3]
plusz.coefficient_idx.append(["a", "i", "b", "j", "c", "k"])
plusz.operator_idx.append(["a", "i"])
plusz.operator_idx.append(["b", "j"])
plusz.operator_idx.append(["c", "k"])
plusz.operator_type.append("s")
plusz.operator_type.append("t0")
plusz.operator_type.append("s")
plusz.num_factor = 1.0/2.0


pluszR1 = ugg()
pluszR1.summation = ["a", "i"]
pluszR1.coefficient = [EOM_TRIPLET_R1]
pluszR1.coefficient_idx.append(["a", "i"])
pluszR1.operator_idx.append(["a", "i"])
pluszR1.operator_type.append("t0")
#pluszRj.num_factor = 1.0/2.0


pluszR2p = ugg()
pluszR2p.summation = ["a", "i", "b", "j"]
pluszR2p.coefficient = [EOM_TRIPLET_R2p]
pluszR2p.coefficient_idx.append(["a", "i", "b", "j"])
pluszR2p.operator_idx.append(["a", "i"])
pluszR2p.operator_idx.append(["b", "j"])
pluszR2p.operator_type.append("t0")
pluszR2p.operator_type.append("s")
pluszR2p.num_factor = 1.0/2.0

pluszR2m = ugg()
pluszR2m.summation = ["a", "i", "b", "j"]
pluszR2m.coefficient = [EOM_TRIPLET_R2m]
pluszR2m.coefficient_idx.append(["a", "i", "b", "j"])
pluszR2m.operator_idx.append(["a", "i"])
pluszR2m.operator_idx.append(["b", "j"])
pluszR2m.operator_type.append("t0")
pluszR2m.operator_type.append("s")

pluszR3 = ugg()
pluszR3.summation = ["a", "i", "b", "j", "c", "k"]
pluszR3.coefficient = [EOM_TRIPLET_R3]
pluszR3.coefficient_idx.append(["a", "i", "b", "j", "c", "k"])
pluszR3.operator_idx.append(["a", "i"])
pluszR3.operator_idx.append(["b", "j"])
pluszR3.operator_idx.append(["c", "k"])
pluszR3.operator_type.append("s")
pluszR3.operator_type.append("t0")
pluszR3.operator_type.append("s")
pluszR3.num_factor = 1.0/2.0

polnR3 = ugg()
polnR3.summation = ["a", "i", "b", "j", "c", "k"]
polnR3.coefficient = [EOM_CC_AMPLITUDE_R]
polnR3.coefficient_idx.append(["a", "i", "b", "j", "c", "k"])
polnR3.operator_idx.append(["a", "i"])
polnR3.operator_idx.append(["b", "j"])
polnR3.operator_idx.append(["c", "k"])
polnR3.operator_type.append("s")
polnR3.operator_type.append("s")
polnR3.operator_type.append("s")
polnR3.num_factor = 1.0/2.0

polnR1 = ugg()
polnR1.summation = ["d", "l"]
polnR1.coefficient = [EOM_CC_AMPLITUDE_R]
polnR1.coefficient_idx.append(["d", "l"])
polnR1.operator_idx.append(["d", "l"])
polnR1.operator_type.append("s")



def operat1(op_idx, spin):

    op = ugg()
    op.operator_idx.append(op_idx)
    op.operator_type.append(spin)

    return op

def operat1x(op_idx, spin):

    op = ugg()
    op.operator_idx.append(op_idx)
    op.operator_type.append(spin)
    op.coefficient.append(OBSERVABLE_X)
    op.coefficient_idx.append(op_idx)
    return op

def operat1xs(op_idx, spin):

    op = ugg()
    op.summation = op_idx
    op.operator_idx.append(op_idx)
    op.operator_type.append(spin)
    op.coefficient.append(OBSERVABLE_AX)
    op.coefficient_idx.append(op_idx)
    return op


def operat1ys(op_idx, spin):

    op = ugg()
    op.summation = op_idx
    op.operator_idx.append(op_idx)
    op.operator_type.append(spin)
    op.coefficient.append(OBSERVABLE_AY)
    op.coefficient_idx.append(op_idx)
    return op



def operat2(op_idx, spin):

    op = ugg()
    op.operator_idx.append(op_idx[0])
    op.operator_idx.append(op_idx[1])
    op.operator_type.append(spin[0])
    op.operator_type.append(spin[1])

    return op

def operat3(op_idx, spin):

    op = ugg()
    op.operator_idx.append(op_idx[0])
    op.operator_idx.append(op_idx[1])
    op.operator_idx.append(op_idx[2])
    op.operator_type.append(spin[0])
    op.operator_type.append(spin[1])
    op.operator_type.append(spin[2])

    return op


#-------------------------------------------------------------CLASS CAS--------------------


AArp = cas()
AArp.operator_idx = ['r', 'p']
AArp.operator_type = ['+', '0']

AAqs = cas()
AAqs.operator_idx = ['q', 's']
AAqs.operator_type = ['+', '0']



AArs = cas()
AArs.operator_idx = ['r', 's']
AArs.operator_type = ['+', '0']


AAsr = cas()
AAsr.operator_idx = ['s', 'r']
AAsr.operator_type = ['+', '0']


AArsg = cas()
AArsg.operator_idx = ['s', 'r']
AArsg.operator_type = ['+', '0']

AApq = cas()
AApq.operator_idx = ['p', 'q']
AApq.operator_type = ['+', '0']


AAps = cas()
AAps.operator_idx = ['p', 's']
AAps.operator_type = ['+', '0']


XXp = cas()
XXp.operator_idx = ['p']
XXp.operator_type = ['+']

YYq = cas()
YYq.operator_idx = ['q']
YYq.operator_type = ['0']


AAqp = cas()
AAqp.operator_idx = ['q', 'p']
AAqp.operator_type = ['+', '0']


# AAqpqp = cas()
# AAqpqp.operator_idx = ['q', 'p', 'r', 's']
# AAqpqp.operator_type = ['+', '0', '+', '0']

# AArsrs = cas()
# AArsrs.operator_idx = ['t', 'u', 'v', 'w']
# AArsrs.operator_type = ['+', '0', '+', '0']

AASqpqp = cas()
AASqpqp.operator_idx = ['p', 'q', 'p', 'q']
AASqpqp.operator_type = ['+', '0', '+', '0']

AATqpqp = cas()
AATqpqp.operator_idx = ['q', 'p', 'u', 't']
AATqpqp.operator_type = ['+', '0', '+', '0']


AApqpq = cas()
AApqpq.operator_idx = ['p', 'q', 't', 'u']
AApqpq.operator_type = ['+', '0', '+', '0']

AAqpqp = cas()
AAqpqp.operator_idx = ['p', 'q', 't', 'u']
AAqpqp.operator_type = ['+', '0', '+', '0']

AArsrs = cas()
AArsrs.operator_idx = ['r', 's', 'v', 'w']
AArsrs.operator_type = ['+', '0', '+', '0']

AAsrsr = cas()
AAsrsr.operator_idx = ['s', 'r', 'w', 'v']
AAsrsr.operator_type = ['+', '0', '+', '0']



AArsc = cas()
AArsc.operator_idx = ['r', 's']
AArsc.operator_type = ['0', '+']

AApqc = cas()
AApqc.operator_idx = ['p', 'q']
AApqc.operator_type = ['0', '+']


AAtuvw = cas()
AAtuvw.operator_idx = ['i', 'j', 'r', 's']
AAtuvw.operator_type = ['+', '+', '0', '0']

AAqprs = cas()
AAqprs.operator_idx = ['q', 'p', 'k', 'l']
AAqprs.operator_type = ['+', '+', '0', '0']


h1cas = cas()
h1cas.summation = ['t', 'u']
h1cas.coefficient = ['h']
h1cas.coefficient_idx.append(['t', 'u'])
h1cas.operator_idx = ['t', 'u']
h1cas.operator_type = ['+', '0']


h2cas = cas()
h2cas.summation = ['t', 'u', 'v', 'w']
h2cas.coefficient = ['gs']
h2cas.coefficient_idx.append(['t', 'u', 'w', 'v'])
h2cas.operator_idx = ['t', 'u', 'v', 'w']
h2cas.operator_type = ['+', '+', '0', '0']
h2cas.num_factor = 0.5

# spin respolved
h2cassr = cas()
h2cassr.summation = ['t', 'u', 'v', 'w']
h2cassr.coefficient = [TWOEL_INT_DIRAC_SPINRES]
h2cassr.coefficient_idx.append(['t', 'u', 'w', 'v'])
h2cassr.operator_idx = ['t', 'u', 'v', 'w']
h2cassr.operator_type = ['+', '+', '0', '0']
h2cassr.num_factor = 0.5


AApq_pp = cas()
AApq_pp.operator_idx = ['p', 'q']
AApq_pp.operator_type = ['+', '+']

AApq_aa = cas()
AApq_aa.operator_idx = ['p', 'q']
AApq_aa.operator_type = ['0', '0']

AAqp_pp = cas()
AAqp_pp.operator_idx = ['q', 'p']
AAqp_pp.operator_type = ['+', '+']

AAqp_aa = cas()
AAqp_aa.operator_idx = ['q', 'p']
AAqp_aa.operator_type = ['0', '0']


AArs_aa = cas()
AArs_aa.operator_idx = ['r', 's']
AArs_aa.operator_type = ['0', '0']

AArs_pp = cas()
AArs_pp.operator_idx = ['r', 's']
AArs_pp.operator_type = ['+', '+']

AAsr_aa = cas()
AAsr_aa.operator_idx = ['s', 'r']
AAsr_aa.operator_type = ['0', '0']

AAsr_pp = cas()
AAsr_pp.operator_idx = ['s', 'r']
AAsr_pp.operator_type = ['+', '+']


AArs_cc = cas()
AArs_cc.operator_idx = ['r', 's']
AArs_cc.operator_type = ['+', '+']

AAr_cc = cas()
AAr_cc.operator_idx = ['r']
AAr_cc.operator_type = ['+']


AAqp_hh = cas()
AAqp_hh.operator_idx = ['q', 'p']
AAqp_hh.operator_type = ['0', '0']

AAp_hh = cas()
AAp_hh.operator_idx = ['p']
AAp_hh.operator_type = ['0']

ttia = ugg()
ttia.operator_idx.append(["i","a"])
ttia.operator_type.append("t0")

ttai = ugg()
ttai.operator_idx.append(["a","i"])
ttai.operator_type.append("t0")

ttkc = ugg()
ttkc.operator_idx.append(["k","c"])
ttkc.operator_type.append("t0")

ttck = ugg()
ttck.operator_idx.append(["c","k"])
ttck.operator_type.append("t0")


etjbia = ugg()
etjbia.operator_idx.append(["j","b"])
etjbia.operator_idx.append(["i","a"])
etjbia.operator_type.append("s")
etjbia.operator_type.append("t0")

etiajb = ugg()
etiajb.operator_idx.append(["i","a"])
etiajb.operator_idx.append(["j","b"])
etiajb.operator_type.append("s")
etiajb.operator_type.append("t0")


etbjai = ugg()
etbjai.operator_idx.append(["b","j"])
etbjai.operator_idx.append(["a","i"])
etbjai.operator_type.append("s")
etbjai.operator_type.append("t0")


etaibj = ugg()
etaibj.operator_idx.append(["a","i"])
etaibj.operator_idx.append(["b","j"])
etaibj.operator_type.append("s")
etaibj.operator_type.append("t0")


etckdl = ugg()
etckdl.operator_idx.append(["c","k"])
etckdl.operator_idx.append(["d","l"])
etckdl.operator_type.append("t0")
etckdl.operator_type.append("s")


etdlck = ugg()
etdlck.operator_idx.append(["d","l"])
etdlck.operator_idx.append(["c","k"])
etdlck.operator_type.append("t0")
etdlck.operator_type.append("s")

etkcld = ugg()
etkcld.operator_idx.append(["k","c"])
etkcld.operator_idx.append(["l","d"])
etkcld.operator_type.append("t0")
etkcld.operator_type.append("s")

etldkc = ugg()
etldkc.operator_idx.append(["l","d"])
etldkc.operator_idx.append(["k","c"])
etldkc.operator_type.append("t0")
etldkc.operator_type.append("s")


letckdl = ugg()
letckdl.operator_idx.append(["c","k"])
letckdl.operator_idx.append(["d","l"])
letckdl.operator_type.append("s")
letckdl.operator_type.append("t0")


letdlck = ugg()
letdlck.operator_idx.append(["d","l"])
letdlck.operator_idx.append(["c","k"])
letdlck.operator_type.append("s")
letdlck.operator_type.append("t0")

letkcld = ugg()
letkcld.operator_idx.append(["k","c"])
letkcld.operator_idx.append(["l","d"])
letkcld.operator_type.append("s")
letkcld.operator_type.append("t0")

letldkc = ugg()
letldkc.operator_idx.append(["l","d"])
letldkc.operator_idx.append(["k","c"])
letldkc.operator_type.append("s")
letldkc.operator_type.append("t0")
