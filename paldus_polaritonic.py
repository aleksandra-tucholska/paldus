from params import *
import paldus_classes
from paldus_cas import cas
from paldus_cas import swap_count
from paldus_cas import cas_to_ugg
from paldus_cas import ugg_to_cas
from paldus_cas import add_spin_from_list
from paldus_cas import *
from paldus_basic import *
from paldus_classes import ugg
from copy import deepcopy
from paldus_classes import arithmetic_string
from paldus_classes import disambiguate
from collections import deque
from collections import Counter
from paldus_classes import pair_permutations
import math
from itertools import product
from itertools import permutations
from multiprocessing import Pool
from paldus_classes import integrate
from paldus_classes import virtual, occupied, general
from fortran_code import *
from templates import *
import sys
import time
import io
from eomccjac import jacobian_loop
from paldus_cas import are_spins_equal
import pickle
import time




def polarit_test():

    A = cas()
    A.operator_idx = ['p']
    A.operator_type = ['0']
    A.operator_stat = [ferm]
    print('A', A)
    B = cas()
    B.operator_idx = ['p']
    B.operator_type = ['+']
    B.operator_stat = [ferm]

#    K = evaluate(A,B)
    K = evaluate(XXp, YYq)
#    print('aars', AArs)
    print('result')
    for elem in K:
        print(elem)
