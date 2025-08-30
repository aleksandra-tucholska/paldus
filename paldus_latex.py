from params import *
import paldus_classes
from paldus_classes import ugg
from copy import deepcopy
from paldus_classes import arithmetic_string
from paldus_classes import disambiguate
from paldus_classes import pair_permutations
import math
from paldus_classes import integrate
from paldus_classes import virtual, occupied, general
from fortran_code import *
from fortran_code_f12 import *
from templates import *
from random import shuffle
from factor import *
from paldus_basic import *


def latex_S_amplitudes(n, interms_lev, all_hash, outer):

    interms = arithmetic_string()
    for i in range(0, n):
        for j in range(0, len(interms_lev[i])):
            interms.append(interms_lev[i][j])


    print('')
    print('all interms')
    print('len przed', len(interms))
    interms = simplify(interms, fixed_fx=['a', 'i', 'b', 'j'])
    print('len po', len(interms))

    print('\equa{')
    for x in range(0, len(interms)):
        if x > 0:
            if x%15 == 0:
                print(interms[x])
                print('}')
                print('\equa{')
            else:
                if x%3 == 0:
                    print(interms[x], '\\\\')
                elif(x-1)%3==0:
                    print('&', interms[x])
                else:
                    print(interms[x])
        else:
            print('&', interms[x], '\\\\')
                            
    print('}')
            







