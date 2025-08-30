from params import *
from paldus_classes import ugg
from paldus_classes import arithmetic_string
from math import sqrt
from copy import deepcopy


bb1 = ugg()
bb1.coefficient = [SVD_U]
bb1.coefficient_idx.append(['i','a', 'X'])

bb2 = ugg()
bb2.coefficient = [SVD_U]
bb2.coefficient_idx.append(['j','b', 'Y'])


bb3 = ugg()
bb3.coefficient = [SVD_X]
bb3.coefficient_idx.append(['k','Q'])

bb4 = ugg()
bb4.coefficient = [SVD_Y]
bb4.coefficient_idx.append(['c','Q'])

bb5 = ugg()
bb5.coefficient = [SVD_C]
bb5.coefficient_idx.append(['l', 'd', 'Q'])

bb6 = ugg()
bb6.coefficient = [SVD_U]
bb6.coefficient_idx.append(['i','c', 'Z'])

bb7 = ugg()
bb7.coefficient = [SVD_U]
bb7.coefficient_idx.append(['j', 'd','W'])

bb8 = ugg()
bb8.coefficient = [SVD_t]
bb8.coefficient_idx.append(['Z','W'])

bb9 = ugg()
bb9.coefficient = [SVD_U]
bb9.coefficient_idx.append(['k','a', 'V'])

bb10 = ugg()
bb10.coefficient = [SVD_U]
bb10.coefficient_idx.append(['l','b', 'U'])

bb11 = ugg()
bb11.coefficient = [SVD_t]  
bb11.coefficient_idx.append(['U','V'])

bb12 = ugg()                                                                                                                                                   
bb12.coefficient = [SVD_U]                                                                              
bb12.coefficient_idx.append(['j','b', 'W'])

bb13 = ugg()                                                                                                                                           
bb13.coefficient = [SVD_t]                                                                                                                                
bb13.coefficient_idx.append(['Z','W'])

bb5a = ugg()
bb5a.coefficient = [SVD_t]
bb5a.coefficient_idx.append(['X', 'Y', 'Y>'])




ww1 = ugg()
ww1.coefficient = [SVD_U]
ww1.coefficient_idx.append(['c','k', 'Z'])

ww2 = ugg()
ww2.coefficient = [CC_AMPLITUDE]
ww2.coefficient_idx.append(['b','b', 'l', 'm'])


ww3 = ugg()
ww3.coefficient = [SVD_B]
ww3.coefficient_idx.append(['l','k', 'Q'])

ww4 = ugg()
ww4.coefficient = [SVD_B]
ww4.coefficient_idx.append(['m','d', 'Q'])

b1 = ugg()
b1.summation = ['a', 'i', 'b', 'j']
b1.coefficient = [EOM_CC_SINGLE_Rl]
b1.coefficient_idx.append(['a','i', 'b', 'j'])


b2 = ugg()
b2.summation = ['k', 'c', 'd', 'l']
b2.coefficient = [EOM_CC_SINGLE_Rl]
b2.coefficient_idx.append(['a','k', 'c', 'i', 'd', 'l'])

b3 = ugg()
b3.summation = ['e', 'm']
b3.coefficient = [S_AMPLITUDE]
b3.coefficient_idx.append(['d','m', 'e', 'k'])


b4 = ugg()
b4.coefficient = [CC_AMPLITUDE]
b4.coefficient_idx.append(['b','j', 'e', 'm'])

b5 = ugg()
b5.coefficient = [OBSERVABLE_X]
b5.coefficient_idx.append(['c','l'])
