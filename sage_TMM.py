from sage.all import *
import numpy

from IPython.Debugger import Pdb

import txt_mixin
reload(txt_mixin)

import time

def sym_eye(N):
    U = identity_matrix(SR, N)
    return U


def aug_wrap(matin, N=4):
    """Take a symbolic matrix and augment it with a row and column of
    zeros, with a 1 in the bottom right corner, so that a TMM matrix
    is compatible with a forcing input."""
    #last_col = zeros((N,1))
    last_col = zero_matrix(SR, nrows=N, ncols=1)
    #Utemp = matin.col_insert(N, last_col)
    mat_out = identity_matrix(SR,N+1)
    mat_out[0:4,0:4] = matin
    ## bottom_row = zeros((1,N+1))
    ## bottom_row[N] = 1
    ## mat_out = Utemp.row_insert(N, bottom_row)
    return mat_out

def diag(i, j):
    if i == j:
        return 1
    else:
        return 0

def eye(N):
    return Matrix(N, N, diag)

s = var('s')


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
    U = sym_eye(5)
    U[1,4] = u
    return U


def U_tsd(k, c, n=5):
    """Return an nxn transfer matrix for a torsional spring and damper."""
    U = sym_eye(n)
    U[1,2] = 1.0/(c*s+k)
    return U

def U_rigid(m, L, I, r, n=5):
    MRx = -m*s**2*(L-r)
    MRth = I*s**2-m*r*s**2*(L-r)
    U = sym_eye(n)
    U[0,1] = L
    U[2,0] = MRx
    U[2,1] = MRth
    U[2,3] = -L
    U[3,0] = m*s**2
    U[3,1] = m*r*s**2
    return U

def z_sub(substr, n=5):
    """Generate a z vector with the subscript substr."""
    states = ['w','th','M','V']
    varnames = [item+substr for item in states]
    mylist = []
    for item in varnames:
        mylist.append([item])
    if n == 5:
        mylist.append(['1'])
    mat_str = str(mylist)
    code = 'z_vect = Matrix(%s)' % mat_str
    exec(code)
    return z_vect


class sage_TMM_Element(object):
    def __init__(self, params, label='', N=4):
        self.params = params
        self.label = label
        self.N = N

    def Get_Mat(self, s):
        raise NotImplementedError
        
    def Get_Aug_Mat(self, s):
        U = self.Get_Mat(s)
        augU = aug_wrap(U, self.N)
        self.augU = augU
        return augU
        

class sage_Beam_Element(sage_TMM_Element):
    def Get_Mat(self, s):
        params = self.params
        label = self.label
##     beta = Symbol('beta')
##     a = Symbol('a')
        d1 = Symbol('d1'+label)
        d2 = Symbol('d2'+label)
        d3 = Symbol('d3'+label)
        d4 = Symbol('d4'+label)
        mu = params['mu']
        EI = params['EI']
        L = params['L']
        #beta = (-s**2*L**4*mu/EI)**0.25
        beta = params['beta']
##     d1 = 0.5*(cos(beta)+cosh(beta))
##     d2 = 0.5*(sinh(beta)-sin(beta))
##     d3 = 0.5*(cosh(beta)-cos(beta))
##     d4 = 0.5*(sin(beta)+sinh(beta))
        a = L**2/EI
        B = Matrix([[d1, L*d4/beta, a*d3/beta**2, -L*a*d2/beta**3], \
                    [beta*d2/L, d1, a*d4/(L*beta), -a*d3/beta**2], \
                    [d3*beta**2/a, L*beta*d2/a, d1, -L*d4/beta], \
                    [-d4*beta**3/(L*a), -d3*beta**2/a, -beta*d2/L, d1]])
        self.U = B
        return B


class sage_TSD_Element(sage_TMM_Element):
    def Get_Mat(self, s):
        k = self.params['k']
        c = self.params['c']
        S = Matrix([[1.0, 0, 0, 0],\
                    [0, 1.0, 1.0/(c*s+k), 0],\
                    [0, 0, 1.0, 0],\
                    [0, 0, 0, 1.0]])
        self.U = S
        return S


class sage_TSD_Generic_Element(sage_TMM_Element):
    def Get_Mat(self, s):
        D_s = self.params['D_s']
        S = Matrix([[1.0, 0, 0, 0],\
                    [0, 1.0, 1.0/D_s, 0],\
                    [0, 0, 1.0, 0],\
                    [0, 0, 0, 1.0]])
        self.U = S
        return S

class sage_Rigid_Mass_Element(sage_TMM_Element):
    def Get_Mat(self, s):
        L = self.params['L']
        m = self.params['m']
        r = self.params['r']
        Iz = self.params['I']
        R = Matrix([[1.,L,0,0],\
                    [0,1.,0,0],\
                    [-m*s**2*(L-r),s**2*Iz-m*s**2*r*(L-r),1.,-L],\
                    [m*s**2,m*s**2*r,0,1.]])
        self.U = R
        return R


class sage_AVS_Element(sage_TMM_Element):
    def Get_Mat(self, s):
        U = sym_eye(self.N)
        self.U = U
        return U

    def Get_Aug_Mat(self, s):
        K_act = self.params['K_act']
        tau = self.params['tau']
        U = sym_eye(self.N+1)
        U[1,self.N] = K_act*tau/(s*(s+tau))
        self.augU = U
        return U

class sage_AVS2_Element(sage_AVS_Element):
    def Get_Aug_Mat(self, s):
        K_act = self.params['K_act']
        p_act1 = self.params['p_act1']
        p_act2 = self.params['p_act2']
        U = sym_eye(self.N+1)
        s1 = 1.0*2.0j*pi#magnitude of s at 1 Hz - fix this point for
                        #changes in p's
        m1 = abs(s1+p_act1)
        m2 = abs(s1+p_act2)
        num = K_act*m1*m2
        U[1,self.N] = num/(s*(s+p_act1)*(s+p_act2))
        self.augU = U
        return U

class sage_AVS1_Element(sage_AVS_Element):
    def Get_Aug_Mat(self, s):
        K_act = self.params['K_act']
        p_act1 = self.params['p_act1']
        #p_act2 = self.params['p_act2']
        U = sym_eye(self.N+1)
        s1 = 1.0*2.0j*pi#magnitude of s at 1 Hz - fix this point for
                        #changes in p's
        m1 = abs(s1+p_act1)
        #m2 = abs(s1+p_act2)
        num = K_act*m1#*m2
        U[1,self.N] = num/(s*(s+p_act1))
        self.augU = U
        return U


class sage_AVS1N_Element(sage_AVS_Element):
    def Get_Aug_Mat(self, s):
        #K_act = self.params['K_act']
        N = self.params['num_act']
        p_act1 = self.params['p_act1']
        #p_act2 = self.params['p_act2']
        U = sym_eye(self.N+1)
##         s1 = 1.0*2.0j*pi#magnitude of s at 1 Hz - fix this point for
##                         #changes in p's
##         m1 = abs(s1+p_act1)
##         #m2 = abs(s1+p_act2)
##         num = K_act*m1#*m2
##         U[1,self.N] = num/(s*(s+p_act1))
        U[1,self.N] = N/(s*(s+p_act1))
        self.augU = U
        return U

class sage_AVS_Generic_Element(sage_AVS_Element):
    def Get_Aug_Mat(self, s):
        G_act = self.params['G_act']
        U = sym_eye(self.N+1)
        U[1,self.N] = G_act
        self.augU = U
        return U

class sage_AVS3_Element(sage_AVS_Element):
    def Get_Aug_Mat(self, s):
        K_act = self.params['K_act']
        p_act1 = self.params['p_act1']
        p_act2 = self.params['p_act2']
        z_act = self.params['z_act']
        U = sym_eye(self.N+1)
        s1 = 1.0*2.0j*pi#magnitude of s at 1 Hz - fix this point for
                        #changes in p's
        m1 = abs(s1+p_act1)
        m2 = abs(s1+p_act2)
        mz = abs(s1+z_act)
        num = K_act*m1*m2/mz
        U[1,self.N] = num*(s+z_act)/(s*(s+p_act1)*(s+p_act2))
        self.augU = U
        return U


class sage_AVS_ThetaFB_Element(sage_AVS_Element):
    def __init__(self, params, Gth=None, label='', N=4):
        sage_AVS_Element.__init__(self, params, label=label, N=N)
        if Gth is None:
            Gth = Symbol('Gth')
        self.Gth = Gth


    def Get_Aug_Mat(self, s):
        ##--------------------------
        ## From numeric TMM code:
        ##--------------------------
        ## Gact = self.Gact_func(s, self.params)
        ## Gth = self.Gth(s)
        ## k_spring = self.params['k_spring']
        ## c_spring = self.params['c_spring']
        ## H = self.params['H']
        ## term1 = 1.0/((1.0 + Gact*Gth*H)*(k_spring + c_spring*s))
        ## term2 = Gact*Gth/(1.0 + Gact*Gth*H)
        ## matout[myrow,2] = term1
        ## matout[myrow,N] = term2
        ##--------------------------
        K_act = self.params['K_act']
        p_act1 = self.params['p_act1']
        H = self.params['H']
        k = self.params['k']
        c = self.params['c']
        #p_act2 = self.params['p_act2']
        U = sym_eye(self.N+1)
        s1 = 1.0*2.0j*pi#magnitude of s at 1 Hz - fix this point for
            #changes in p's
        m1 = abs(s1+p_act1)
        #m2 = abs(s1+p_act2)
        num = K_act*m1#*m2
        Gact = num/(s*(s+p_act1))
        Gth = self.Gth
        term1 = 1.0/((1.0 + Gact*Gth*H)*(k + c*s))
        term2 = Gact*Gth/(1.0 + Gact*Gth*H)
        U[1,2] = term1
        U[1,self.N] = term2
        self.augU = U
        return U


class sage_Forcing_Element(sage_TMM_Element):
    def Get_Mat(self, s):
        U = sym_eye(self.N)
        self.U = U
        return U

    def Get_Aug_Mat(self, s):
        fv = self.params['fv']
        U = sym_eye(self.N+1)
        U[0:4,4] = fv
        self.augU = U
        return U


def find_submat(Uin):
    submat = Uin[2:4, 2:4]
    return submat

def find_submat_inv(Uin):
    submat = find_submat(Uin)
    submati = submat.inv()
    return submati

def find_base_vector(Uin):
    submati = find_submat_inv(Uin)
    Uc4 = Uin[2:4,4]
    MbVb = -1.0*(submati*Uc4)
    z_b = zeros((5,1))
    z_b[2] = MbVb[0]
    z_b[3] = MbVb[1]
    z_b[-1] = 1.0
    return z_b

    
def cse_tuples_to_txtlist(tuplelist, ws=" "*4):
    """Take a list of tuples returned as the first output from
    sympy.cse and convert it to a txtlist of valid Python code."""
    mylist = None
    for var, expr in tuplelist:
        curline = ws + '%s = %s' % (var, expr)
        if mylist is None:
            mylist = [curline]
        else:
            mylist.append(curline)
    return mylist


def list_to_array_str(listin, outlabel, ws=" "*4):
    """Assume that the list contains elements of a square matrix and
    return valid Python code to convert the list back to an array.  Do
    this using numpy.reshape."""
    T = numpy.array(listin)
    N = T.shape[0]
    n = numpy.sqrt(N)
    mat = numpy.reshape(T, (n,n))
    str_mat = None
    for row in mat:
        row_list = None
        for ent in row:
            ent_str = str(ent)
            if row_list is None:
                row_list = [ent_str]
            else:
                row_list.append(ent_str)
        row_str = '[' + ', '.join(row_list)+']'
        if str_mat is None:
            str_mat = [row_str]
        else:
            str_mat.append(row_str)
    n1 = len(outlabel)
    n2 = len(' = array([')
    ws2 = ws + " "*(n1+n2)
    mat_str = 'array([' + (', \\\n'+ws2).join(str_mat)+'])'
    return mat_str

def create_output_lines(outlist, outlabels, ws=" "*4):
    last_line = ''
    for expr, label in zip(outlist, outlabels):
        curline = ws + '%s = %s\n' % (label, expr)
        last_line += curline
    return last_line

def cse_to_txtlist(expr_list, outlabels, ws=" "*4):
    if type(outlabels) == str:
        outlabels = [outlabels]
    t1 = time.time()
    tuplist, out = cse(expr_list)
    t2 = time.time()
    print('cse time='+str(t2-t1))
    mylist = cse_tuples_to_txtlist(tuplist)
    if len(out) > len(outlabels):
        out_str = list_to_array_str(out, outlabel, ws=ws)
        last_line = ws + outlabels[0] + ' = ' + out_str
    else:
        last_line = create_output_lines(out, outlabels, ws=ws)
    return_line = ws + 'return ' + ', '.join(outlabels)
    if mylist is None:#no tuples were returned by cse
        mylist = [last_line]
    else:
        mylist.append(last_line)
    mylist.append(return_line)
    return mylist

def cse_to_file(expr_list, filename, outlabels, funcname, \
                inputs=[], ws=' '*4, headerfile=None, \
                replace_dict={}):
    line0 = 'from __future__ import division'
    line1 = 'from scipy import *'
    line2 = 'def '+funcname +'(' + ', '.join(inputs) + '):'
    preamble = [line0, line1, '', line2]
    mylist = []
    if headerfile:
        headerlist = txt_mixin.read(headerfile)
        mylist.extend(headerlist)
    mylist.extend(cse_to_txtlist(expr_list, outlabels, ws=ws))
    if replace_dict:
        mylist = txt_mixin.txt_list(mylist)
        for key, value in replace_dict.iteritems():
            mylist.replaceall(key,value)
    mylist = preamble + mylist#don't do the search and replace in the
                              #preamble
    txt_mixin.dump(filename, mylist)
    
                  

if __name__ == '__main__':

    from optparse import OptionParser

    usage = 'usage: %prog [options]'
    parser = OptionParser(usage)


    parser.add_option("-n", "--name", dest="name", \
                      help="module name for the numeric Bode module", \
                      default="sympy_bodes_base_mass.py", type=str)

    (options, args) = parser.parse_args()

    mod_name = options.name

    t_start = time.time()
    s = Symbol('s')
    ##################################################
    #
    # Create the Elements
    #
    ##################################################
    #---------------------
    mu = Symbol('mu')
    EI = Symbol('EI')
    L1 = Symbol('L1')
    beta1 = Symbol('beta1')
    params1 = {'mu':mu, 'EI':EI, 'L':L1, 'beta':beta1}
    L2 = Symbol('L2')
    beta2 = Symbol('beta2')
    params2 = {'mu':mu, 'EI':EI, 'L':L2, 'beta':beta2}
    beam1 = sage_Beam_Element(params1, label='_1')
    beam2 = sage_Beam_Element(params2, label='_2')
    #---------------------
    k_clamp = Symbol('k_clamp')
    c_clamp = Symbol('c_clamp')
    TSDparams = {'k':k_clamp, 'c':c_clamp}
    TSD_clamp = sage_TSD_Element(TSDparams)
    #---------------------
    a_m = Symbol('a_m')
    a_L = Symbol('a_L')
    a_r = Symbol('a_r')
    a_I = Symbol('a_I')
    a_gain = Symbol('a_gain')
    am_params = {'m':a_m, 'L':a_L, 'r':a_r, 'I':a_I}
    Accel_Mass = sage_Rigid_Mass_Element(am_params)
    #---------------------
    b_m = Symbol('b_m')
    b_L = Symbol('b_L')
    b_r = Symbol('b_r')
    b_I = Symbol('b_I')
    b_gain = Symbol('b_gain')
    bm_params = {'m':b_m, 'L':b_L, 'r':b_r, 'I':b_I}
    Base_Mass = sage_Rigid_Mass_Element(bm_params)
    #----------------------
    k_spring = Symbol('k_spring')
    c_spring = Symbol('c_spring')
    TSDparams = {'k':k_spring, 'c':c_spring}
    TSD_spring = sage_TSD_Element(TSDparams)
    #---------------------
    K_act = Symbol('K_act')
    p_act1 = Symbol('p_act1')
    p_act2 = Symbol('p_act2')
    tau = Symbol('tau')
    z_act = Symbol('z_act')
    AVS_params = {'K_act':K_act, 'p_act1':p_act1, 'p_act2':p_act2, \
                  'z_act':z_act, 'tau':tau}
    #AVS = sage_AVS3_Element(AVS_params)
    AVS = sage_AVS1_Element(AVS_params)
    #AVS = sage_AVS_Element(AVS_params)
    ##################################################
    #
    # Forced Response
    #
    ##################################################
    U0 = AVS.Get_Aug_Mat(s)
    U1 = TSD_spring.Get_Aug_Mat(s)
    U2 = Base_Mass.Get_Aug_Mat(s)
    U3 = TSD_clamp.Get_Aug_Mat(s)
    U4 = beam1.Get_Aug_Mat(s)
    U5 = Accel_Mass.Get_Aug_Mat(s)
    U6 = beam2.Get_Aug_Mat(s)
    ta = time.time()
    Uaug = U6*(U5*(U4*(U3*(U2*(U1*U0)))))
    tb = time.time()
    U_LR = Uaug[2:4, 2:4]
    U_LRi = U_LR.inv()
    tc = time.time()
    Uc4 = Uaug[2:4,4]
    MbVb = -1.0*(U_LRi*Uc4)
    z_b = zeros((5,1))
    z_b[2] = MbVb[0]
    z_b[3] = MbVb[1]
    z_b[-1] = 1.0
    z_enc = U2*(U1*(U0*z_b))
    #z_enc = U0*z_b
    z_accel = U5*(U4*(U3*z_enc))
    a_out = s**2*z_accel[0]*a_gain
    enc_gain = 180.0/pi*1024.0/360.0
    th_out = z_enc[1]*enc_gain

    tcse_start = time.time()
##     cse_to_file([th_out, a_out], 'sympy_bodes_debug.py',\
##                 ['th_out','a_out'],'Bodes',\
##                 inputs=['s','params'], headerfile='header.py')
    cse_to_file([th_out, a_out], mod_name,\
                ['th_out','a_out'],'Bodes',\
                inputs=['s','params'], headerfile='header.py')


##     Ulist = [U0, U1, U2, U3, U4, Uaug]
##     Unames = ['U0','U1','U2','U3','U4','Uaug']
##     for U, name in zip(Ulist, Unames):
##         cse_to_file(U, name+'.py', 'U', name, inputs=['s','params'], \
##                     headerfile='header.py')


    tend = time.time()
    print('total time='+str(tend-t_start))
    #unforced subdet analysis
##     Bz = beam1.Get_Mat(s)
##     Bz2 = beam2.Get_Mat(s)
##     Sclamp = TSD.Get_Mat(s)
##     R = RigidMass.Get_Mat(s)
##     U2 = Bz2*(R*(Bz*Sclamp))
##     U = R*(Bz*Sclamp)
##     U22 = U[2,2]
##     U33 = U[3,3]
##     U23 = U[2,3]
##     U32 = U[3,2]
    
##     det = U22*U33-U23*U32
##     cse_to_file(U, 'Usympy.py', 'U', 'U_sympy', inputs=['s','params'])
##     cse_to_file(Bz, 'Bzsympy.py','B','Bz_sympy', inputs=['s','params'])
##     cse_to_file(U2, 'Usympy_two_piece.py', 'U', 'U_sympy_two_piece', inputs=['s','params'])
##     det2 = U2[2,2]*U2[3,3]-U2[2,3]*U2[3,2]
##     cse_to_file(det2, 'det_sympy_two_piece.py','det','det_two_piece',\
##                 inputs=['s','params'])
