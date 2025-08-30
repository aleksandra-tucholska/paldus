from params import *
import paldus_classes
from paldus_basic import *
from paldus_classes import ugg
from copy import deepcopy
from paldus_classes import arithmetic_string
from paldus_classes import disambiguate
from collections import deque
from paldus_classes import pair_permutations
import math
from itertools import product
from itertools import permutations
from paldus_classes import integrate
from paldus_classes import virtual, occupied, general
from fortran_code import *
from templates import *
import sys
import time
import pickle
from eomccjac import jacobian_loop

def execute_ccsd_f12():

    t2_f12 = arithmetic_string(t2fa, t2fb)
    print(t2_f12)

    r1 = hamiltoniant + evaluate(hamiltoniant, t2_f12)
    r2 = hamiltoniant + evaluate(hamiltoniant, t2_f12) +evaluate(hamiltoniant, t2_f12, t2_f12).scale(0.5)

    
    r1list = []
    r2list = []
    nu1c = ugg()
    nu1c.operator_idx.append(['i', 'a'])
    nu1c.operator_type.append('s')
    

    nu2c = ugg()
    nu2c.operator_idx.append(['i', 'a'])
    nu2c.operator_idx.append(['j', 'b'])
    nu2c.operator_type.append('s')
    nu2c.operator_type.append('s')



    for x in r1:
        disambiguate(x, nu1c)
        x = x.fromleft(nu1c)
        r1list.append(x)
        print(x)


    for x in r2:
        disambiguate(x, nu2c)
        x = x.fromleft(nu2c)
        r2list.append(x)
        print(x)

    r2listcp = deepcopy(r2list)



    print('do calkowania r1 bedzie', len(r1))
    print('do calkowania r2 bedzie', len(r2))
    

    z = Parallel(n_jobs=30,verbose=100)(delayed(integrate)(r1list[i]) for i in range(0, len(r1list)))
    print('r111')
    for i in range(0, len(z)):
        rint = deepcopy(z[i])

        rint.exec_delta()
        rsimp = simplify(rint)
        for x in rsimp:
            print(x)

    z2 = Parallel(n_jobs=30,verbose=100)(delayed(integrate)(r2list[i]) for i in range(0, len(r2list)))

    
    print('r222')
    for i in range(0, len(z2)):
        rint = deepcopy(z2[i])

        rint.exec_delta()
        rsimp = simplify(rint)
        for x in rsimp:
            print(x)


    sys.exit(0)
    print('rint1')
    rint1 = r1.integrate(bra = ['a', 'i'], braspin = ['s']).scale(0.5)
    print('rint21')
    rint21  = r2.integrate(bra = ['a', 'i', 'b', 'j'], braspin = ['s', 's']).scale(1.0/3.0)
    print('rint22')
    rint22 = r2.integrate(bra = ['a', 'j', 'b', 'i'], braspin = ['s', 's']).scale(1.0/6.0)
    print('sum')
    rint2 = (rint21 + rint22)

    rint1.exec_delta()
    rint2.exec_delta()


    print('simplify')
    rsimp1 = simplify(rint1)
    rsimp2 = simplify(rint2)

    for x in rsimp1:
        x.optimize()
    for x in rsimp2:
        x.optimize()

    rsimp1.cleanup()
    rsimp2.cleanup()


    print('rsimp2')
    for x in rsimp1:
        print(x)
    print('rsimp2')
    for x in rsimp2:
        print(x)

    sys.exit(0)

    function_template_ccsd(rsimp1, 1)
    function_template_ccsd(rsimp2, 2)


