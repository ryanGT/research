from scipy import *
import sympy

#from IPython.Debugger import Pdb
#my_syms = ['Gact','v','c','s','k','th_d', 'tha','xb','thb','Mb','Vb']


## for sym in my_syms:
##     cmd = '%s = sympy.Symbol("%s")' % (sym, sym)
##     exec(cmd)

#sympy.var(my_syms)

s = sympy.Symbol('s')

def diag(i, j):
    if i == j:
        return 1
    else:
        return 0


def eye(N):
    return sympy.Matrix(N, N, diag)


def find_row(sym_in, var_list, N=5):
    row_out = []
    for var in var_list:
        temp = sym_in.subs(var, 1)
        for other in var_list:
            if other != var:
                temp = temp.subs(other, 0)
        row_out.append(temp)
    if len(row_out) < N:
        row_out.append(1)
    return row_out


def U_hyd_act(u):
    """Return a 5x5 hydraulic acutator transfer matrix with input u."""
    U = eye(5)
    U[1,4] = u
    return U


def U_tsd(k, c, n=5):
    """Return an nxn transfer matrix for a torsional spring and damper."""
    U = eye(n)
    U[1,2] = 1.0/(c*s+k)
    return U

def z_sub(substr, n=5):
    """Generate a z vector with the subscript substr."""
    states = ['x','th','M','V']
    varnames = [item+substr for item in states]
    mylist = []
    for item in varnames:
        mylist.append([item])
    if n == 5:
        mylist.append(['1'])
    mat_str = str(mylist)
    code = 'z_vect = sympy.Matrix(%s)' % mat_str
    exec(code)
    return z_vect


def U_rigid(m, L, I, r, n=5):
    MRx = -m*s**2*(L-r)
    MRth = I*s**2-m*r*s**2*(L-r)
    U = eye(n)
    U[0,1] = L
    U[2,0] = MRx
    U[2,1] = MRth
    U[2,3] = -L
    U[3,0] = m*s**2
    U[3,1] = m*r*s**2
    return U
