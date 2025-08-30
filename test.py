from paldus_classes import *




sys.exit(0)

a = ugg()
a.summation = ["b", "j", "c", "k"]
a.coefficient = ["g", "t", "t"]
a.coefficient_idx = [["b", "j", "c", "k"], ["a", "j"], ["b", "k", "c", "i"]]
a.num_factor = 1.

b = ugg()
b.summation = ["b", "c", "j", "k"]
b.coefficient = ["g", "t", "t"]
b.coefficient_idx = [["b", "j", "c", "k"], ["a", "k"], ["b", "i", "c", "j"]]
b.num_factor = 1.

a.standarize()
b.standarize()

print(a)
print(b)
l = a == b
print(l)
