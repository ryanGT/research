from scipy import *
import sympy
reload(sympy)
from sympy import Wild, Basic

import copy

import sympy_utils
reload(sympy_utils)

import sympy_TMM# import find_row

s = sympy.Symbol("s")

N = 3

def Us(k):
    U = sympy_utils.eye(N)
    U[0,1] = 1/k
    return U


def Um(m):
    U = sympy_utils.eye(N)
    U[1,0] = m*s**2
    return U


def Uf(F):
    U = sympy_utils.eye(N)
    U[1,2] = -F
    return U


def z(x, F):
    z_mat = sympy.Matrix([[x],[F],[1]])
    return z_mat


def find_row(sym_in, var_list, N=3):
    return sympy_TMM.find_row(sym_in, var_list, N=N)

    
def Symbols(eq):
    if type(eq) is str:
        import re
        return re.findall('([A-Za-z_][A-Za-z_0-9]*)',eq)
    else:
        return [x for x in eq.atoms() if x.is_Symbol]
    

def as_whole_frac(eq):
    n,d = eq.as_numer_denom()
    q,r = sympy.div(n,d,Symbols(eq))
    return q+r/d


def as_fraction(eq):
    n,d = eq.as_numer_denom()
    q,r = sympy.div(n,d,Symbols(eq))
    return (q*d+r)/d


def my_compare(a, b):
   main_var = s
   p1, p2, p3 = Wild("p1"), Wild("p2"), Wild("p3")
   r_a = a.match(p1 * s**p3)
   r_b = b.match(p1 * s**p3)
   if r_a is not None and r_b is not None:
       c = Basic.compare(r_a[p3], r_b[p3])
       if c!=0:
           return c

   return Basic._compare_pretty(a,b)


my_list = ['m1','m2','F','k','k1','k2','Gc','v','x1','x2','F1','F2','xd']

sympy_utils.declare_many_sims(my_list, globals())
