import sys
from templates import *
from paldus_main import evaluate
from paldus_main import execute

global_dict = {'nu' : nu, 'nu2' : nu2, 'hamiltonian' : hamiltonian}

def calculate_max_nt(X, Y, Bra, Rop):
    """ If Y = T, then this procedure
    is called to calculate max excitation
    caused by T
    """

    if Bra == '1':
        Bra_exc = 1
    
    if X == 'hamiltonian':
        min_s_x = -2
    
    if Rop != '':
        k = 0
        for i in global_dict[Rop].operator_idx:
            k +=1
    
    print(k,'k')
    nt = Bra_exc - min_s_x - k
    print('nt', nt)
    
    return nt

def nested_number(X, Y):
    """ Returns number of how many nested
    commutators are possible computing
    e^(-Y) X e^(Y) and using Baker-Campbell-Hausdorf
    expansion.
    """

    if X == 'hamiltonian' and Y == 'T':
        k = 4

    return k

name = sys.argv[1]

f = open(name, 'r')

for line in f:
    if '@BCH' in line:
        BCH = line[len('@BCH='):len(line)-1]
        print(BCH)
    if '@X-operator' in line:
        X = line[len('@X-operator='):len(line)-1]
        print(X)
                       
    if '@Y-operator' in line:
        Y = line[len('@Y-operator='):len(line)-1]
        print(Y)
    if '@Right_commutator_operator' in line:
        Rop = line[len('@Right_commutator_operator='):len(line)-1]
        R = True
        print(Rop)
    if '@Bra' in line:
        Bra = line[len('@Bra='):len(line)-1]
        print(Bra)
    if '@Ket' in line:
        Ket = line[len('@Ket='):len(line)-1]
        print(Ket)

if R == False:
    Rop = ''

if BCH == True:
    N_nest = nested_number(X, Y)
    
if R == True:
    if Y == 'T':
        NT = calculate_max_nt(X, Y, Bra, Rop)

list_to_eval = list_to_evaluate(




#r = evaluate(hamiltonian, nu2) + evaluate(hamiltonian, t1, nu2)
#execute(r, 1, 0, 2)
