from pylab import *
from scipy import *
from scipy import linalg, signal, integrate, optimize
from scipy.linalg import eigvals
import numpy

import controls

import rwkmisc, rwkbode, txt_mixin
from rwkmisc import colwise, rowwise
from IPython.core.debugger import Pdb

import copy, sys, os

import SFLR_TF_models
#reload(SFLR_TF_models)

import plotting_mixin


class SS_model(plotting_mixin.item_that_plots):
    def _init2(self):
        nr,nc = self.A.shape
        self.N = nr
        nrC, ncC = self.C.shape
        self.M = nrC
        self.I = eye(self.N)


    def __init__(self, A, B, C, D=0.0, **kwargs):
        """Note that kwargs is just used to ignore extra keyword
        arguements passed in by inherited classes."""
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.use_dig = False
        self.mysat = 200.0
        self.accel = False
        if self.A is not None:
            self._init2()


    def calc_Ghat(self, K):
        """Find the Ghat matrix that would correspond to using K for
        feedback: Ghat = self.G - dot(self.B, K)"""
        self.Ghat = self.G - dot(self.H, K)
        return self.Ghat


    def calc_Ahat(self, K=None):
        """Find the Ahat matrix that would correspond to using K for
        feedback: Ahat = self.A - dot(self.B, K)

        If K is None, use self.K"""
        if K is None:
            K = self.K
        self.Ahat = self.A - dot(self.B, K)
        return self.Ahat


    def build_K(self, C, inds=[0]):
        """Create a state feedback vector K by assigning the elements
        of C to the elements of K that correspond to inds."""
        K = zeros((1,self.N))
        for ci, ind in zip(C, inds):
            K[0,ind] = ci
        return K


    def find_Ghat_poles(self, K):
        Ghat = self.calc_Ghat(K)
        poles = eigvals(Ghat)
        return poles


    def pole_cost(self, C, inds=[0]):
        K = self.build_K(C, inds=inds)
        Ghat = self.calc_Ghat(K)
        poles = eigvals(Ghat)
        real_ind = abs(imag(poles)) < 1e-8
        real_poles = poles[real_ind]
        imag_ind = abs(imag(poles)) > 1e-8
        imag_poles = poles[imag_ind]
        if len(imag_poles) > 2:
            imag_poles.sort()
            imag_poles = imag_poles[-2:]
        #cost = abs_poles[-1] + abs_poles[-2]
        cost = sum(abs(imag_poles))
        if real_poles.any():
            if abs(real_poles).max() > 0.95:
                cost += 1e4
        if any(abs(poles) > 1.0):
            cost += 1e6
        return cost



    def plot_states_many_figs(self, fi=3, clear=True, \
                              attr='X_dig'):
        mat = getattr(self, attr)
        mat = colwise(mat)
        nr, nc = mat.shape
        for n in range(nc):
            cur_col = mat[:,n]
            label = attr + '%i' % n
            ax = self._prep_ax(fi=fi, clear=clear)
            ax.plot(self.t, cur_col, label=label)
            ax.set_title(label)
            fi += 1

    def _ones_on_super(self):
        mat = zeros((self.N,self.N))
        for i in range(self.N-1):
            mat[i, i+1] = 1.0#ones on super-diagonal
        return mat


    def calc_freq_resp_one_s(self, s_i):
        mat = s_i*self.I - self.A
        mati = linalg.inv(mat)
        M = dot(mati, self.B)
        comp = dot(self.C,M)
        return comp


    def calc_freq_resp(self, f):
        if not hasattr(self, 'I'):
            self.I = eye(self.N)
        svect = 2.0j*pi*f
        comp_mat = zeros((self.M, len(svect)), dtype='complex128')

        for i, s_i in enumerate(svect):
            comp_i = self.calc_freq_resp_one_s(s_i)
            comp_mat[:,i] = squeeze(comp_i)
        self.comp_mat = comp_mat
        return comp_mat


    def get_params(self, attrlist):
        mydict = {}
        for key in attrlist:
            val = getattr(self, key)
            mydict[key] = val
        return mydict


    def save_params(self, pklpath, attrlist):
        mydict = self.get_params(attrlist)
        rwkmisc.SavePickle(mydict, pklpath)


    def discretize(self, T):
        """This function follows the approach of Ogata's
        'Discrete-Time Control Systems 2nd Edition' pages 314-315.

        x((k+1)T) = G(T)x(kT) + H(T)u(kT)
        y(kT) = C(x(kT) + D u(kT)

        G(T) = exp(A*T) (matrix exponential)
        H(T) = (integral from 0 to T of exp(A*lambda) d lambda)*B"""
        self.T = T
        self.G = linalg.expm(self.A*T)

        def Hprime(t, i):
            M = linalg.expm(self.A*t)
            H = dot(M, self.B)
            return H[i]

        self.H = zeros((self.N,1))
        for i in range(self.N):
            out = integrate.quad(Hprime, 0, T, args=(i,))
            self.H[i] = out[0]
        return self.G, self.H


    def _get_A_and_B(self):
        if self.use_dig:
            A = self.G
            B = self.H
        else:
            A = self.A
            B = self.B
        return A, B


    def create_ax(self, fi=1, clear=True):
        fig = figure(fi)
        if clear:
            fig.clf()
        self.ax = fig.add_subplot(1,1,1)
        return self.ax


    def get_or_create_ax(self, fi=1, clear=True):
        fig = figure(fi)
        if clear:
            fig.clf()
        if len(fig.axes) == 1:
            self.ax = fig.axes[0]
        else:
            self.ax = fig.add_subplot(1,1,1)
        return self.ax


    def check_pole_locations(self, K):
        """Find the pole locations that would result from using
        feedback vector K."""
        A, B = self._get_A_and_B()
        Atemp = A + dot(rwkmisc.colwise(B), K)
        poles = linalg.eigvals(Atemp)
        return poles


    def find_poles(self):
        A, B = self._get_A_and_B()
        self.poles = linalg.eigvals(A)
        return self.poles


    def _plot_poles(self, poles, lt='bo', fi=1, clear=True, label=None):
        #Pdb().set_trace()
        ## if not hasattr(self, 'ax'):
        ##     self.create_ax(fi=fi, clear=clear)
        ## if self.ax.figure.number != fi:
        ##    self.create_ax(fi=fi, clear=clear)#<-- I think we should do this no matter what
        self.get_or_create_ax(fi=fi, clear=clear)
        if not hasattr(self, 'poles'):
            self.find_poles()
        self.ax.plot(real(poles), imag(poles), \
                     lt, label=label)


    def plot_poles(self, *args, **kwargs):
        self._plot_poles(self.poles, *args, **kwargs)


    def plot_CL_poles(self, *args, **kwargs):
        self._plot_poles(self.CL_poles, *args, **kwargs)


    def plot_OL_poles(self, *args, **kwargs):
        if not hasattr(self, 'OL_poles'):
            self.find_OL_poles()
        self._plot_poles(self.OL_poles, *args, **kwargs)


    def state_space_digital_root_locus(self, fi=2, clear=False, \
                                       K_mul=None, pole_lt='k^'):
        orig_K = copy.copy(self.K)
        if K_mul is None:
            steps = 10.0
            dK = 1.0/steps
            K_mul = arange(0,1+dK/2.0, dK)
        if clear:
            figure(fi)
            clf()

        for mul in K_mul:
            self.K = mul*orig_K
            self.find_CL_poles()
            self.plot_CL_poles(lt=pole_lt, fi=fi, \
                               clear=False)
        self.K = orig_K
        self.find_CL_poles()


    def try_K(self, K, u, pole_lt='k^', \
              pole_fi=2, step_fi=31, \
              clear_step=True, clear_poles=False, \
              label=None, E=None):
        self.K = K
        self.find_CL_poles()
        self.plot_CL_poles(lt=pole_lt, fi=pole_fi, \
                           clear=clear_poles, \
                           label=label)

        if E is None:
            self.E = 1.0
            self.digital_lsim(u, closed_loop=True)
            self.E = 200/self.Y_dig[0,-1]
        else:
            self.E = E

        self.digital_lsim_with_obs(u)
        self.plot_digital_lsim_results(fi=step_fi, clear=clear_step)

        #self.digital_lsim(u, closed_loop=True)
        #self.plot_digital_lsim_results(fi=step_fi, clear=clear_step)




    def phi_des_of_A(self, descoeffs):
        """Evalute the charcteristic polynominal using
        descoeffs[0]*A**n+descoeffs[1]*A**(n-1)"""
        A, B = self._get_A_and_B()
        N = len(descoeffs)-1
        assert N == self.N, "bad length of descoeffs"
        for n in range(N+1):
            a_n = descoeffs[N-n]
            if n == 0:
                phi = eye(N)*a_n
            else:
                if n == 1:
                    An = A
                else:
                    An = dot(An, A)
                phi += An*a_n
        return phi


    def controlability(self):
        A, B = self._get_A_and_B()
        N = self.N
        for n in range(N):
            if n == 0:
                collist = [B]
            else:
                if n == 1:
                    An = A
                else:
                    An = dot(An, A)
                curcol = dot(An, B)
                collist.append(curcol)
        self.P = column_stack(collist)
        return self.P


    def _get_C(self):
        """For now, I am handling multiple output systems by using
        just the first output.  So, this method just return self.C for
        most systems, but grabs the first row of C for
        systems deliberately ignoring the accel or other outputs."""
        return self.C


    def observability(self, out_ind=0, tryMO=False):
        A, B = self._get_A_and_B()
        C = self._get_C()
        #force using only first output for now
        nrC, ncC = C.shape
        if not tryMO:
            if nrC > 1:
                C = atleast_2d(C[out_ind,:])
        N = self.N
        for n in range(N):
            if n == 0:
                rowlist = [C]
            else:
                if n == 1:
                    CAn = dot(C, A)
                else:
                    CAn = dot(CAn, A)
                currow = CAn
                rowlist.append(currow)
        self.Q = row_stack(rowlist)
        return self.Q


    def acker(self, despoles, discretize=True, \
              obs=False, out_ind=0, tryMO=False, debug=1):
        """obs True means find the observer gain vector L and is used
        in observer design.  obs False means find the state-feedback
        gain vector used in pole-placement control design."""
        if self.use_dig:
            if discretize:
                z_poles = [exp(s*self.T) for s in despoles]
            else:
                z_poles = despoles
            descoeffs = poly(z_poles)
        else:
            descoeffs = poly(despoles)
        phi = self.phi_des_of_A(descoeffs)
        if obs:
            Q = self.observability(out_ind=out_ind, tryMO=tryMO)
            u,s,v = linalg.svd(Q)
            r = numpy.sum(s > 1e-9)
            mytup = (self.N, r, str(s))
            assert r == self.N, 'Problem with the rank of Q: self.N = %s, r = %s, singular values = %s' % mytup
            if debug:
                msg = 'self.N = %s, r = %s, singular values = %s' % mytup
                print(msg)

            Qinv = linalg.inv(Q)
            en = zeros((self.N,1))
            en[-1,0] = 1.0
            M = dot(Qinv, en)
            K = dot(phi, M)
            self.Ke = K
        else:
            P = self.controlability()
            Pinv = linalg.inv(P)
            en = zeros((1, self.N))
            en[0,-1] = 1.0
            M = dot(Pinv, phi)
            K = dot(en, M)
            self.K = K
        return K


    def digital_lsim(self, u, X0=None, closed_loop=False):
        self.r = u
        if not hasattr(self, 'E'):
            self.E = 1.0
        if X0 is None:
            X0 = zeros((self.N,1))
        N2 = len(u)
        X = zeros((self.N, N2))
        X[:,0] = squeeze(X0)

        if closed_loop:
            self.calc_Ghat(K=self.K)
            G = self.Ghat
        else:
            G = self.G
        H = self.H
        C = self.C

        Ny, Nx = self.C.shape
        Y = zeros((Ny, N2))
        Y[:,0] = squeeze(dot(C, X0))

        prev_x = X0
        for k in range(1,N2):
            next_x = dot(G, prev_x) + H*u[k-1]*self.E
            Y[:,k] = squeeze(dot(C, next_x))
            X[:,k] = squeeze(next_x)
            prev_x = next_x

        self.X_dig = X
        self.Y_dig = Y
        self.v = squeeze(self.E*u - dot(self.K, self.X_dig))
        return self.Y_dig


    def calc_feeback_matrices(self, K, E=1.0):
        """Calculate feedback matrices for A and B or G and H
        depending on whether or not self.use_dig is True."""
        if self.use_dig:
            A = self.G
            B = self.H
            if not hasattr(self, 'G_ol'):
                self.G_ol = copy.copy(A)
            if not hasattr(self, 'H_ol'):
                B_temp = copy.copy(B)
                self.H_ol = rwkmisc.colwise(B_temp)
            A_ol = self.G_ol
            B_ol = self.H_ol
        else:
            A = self.A
            B = self.B
            if not hasattr(self, 'A_ol'):
                self.A_ol = copy.copy(A)
            if not hasattr(self, 'B_ol'):
                B_temp = copy.copy(B)
                self.B_ol = rwkmisc.colwise(B_temp)
            A_ol = self.A_ol
            B_ol = self.B_ol


        self.K = rwkmisc.rowwise(K)
        self.E = E

        A = A_ol - dot(B_ol, self.K)
        if isscalar(E):
            B = B_ol*E
        else:
            B = dot(B_ol, E)

        if self.use_dig:
            self.G = A
            self.H = B
        else:
            self.A = A
            self.B = B
            self.create_LTI()

        return A, B


    def sat(self, vin):
        if not hasattr(self, 'mysat'):
            self.mysat = 200
        vmax = self.mysat
        vmin = -1*self.mysat
        if vin > vmax:
            return vmax
        elif vin < vmin:
            return vmin
        else:
            return vin


    def calc_v(self, k):
        vtemp = self.E*self.r[k] - squeeze(dot(self.K, self.X_tilde[:,k-1]))
        return self.sat(vtemp)


    def make_t(self):
        N2 = len(self.r)
        self.nvect = arange(N2)
        self.t = self.nvect*self.T


    def digital_lsim_with_obs(self, r, X0=None, X0_tilde=None, \
                              exp_Y=None):
        self.r = r
        if X0 is None:
            X0 = zeros((self.N,1))
        if X0_tilde is None:
            X0_tilde = zeros((self.N,1))
        N2 = len(r)
        self.v = zeros_like(r)
        self.X_dig = zeros((self.N, N2))
        self.X_tilde = zeros((self.N, N2))
        self.X_dig[:,0] = squeeze(X0)
        self.X_tilde[:,0] = squeeze(X0_tilde)

        if hasattr(self, 'G_ol'):
            G = self.G_ol
        else:
            G = self.G

        if hasattr(self, 'H_ol'):
            H = self.H_ol
        else:
            H = self.H
        #C = self.C
        C = self._get_C()
        #Ke = self.Ke
        #K = self.K

        Ny, Nx = self.C.shape
        if exp_Y is not None:
            use_exp_Y = True
            self.Y_dig = exp_Y
        else:
            use_exp_Y = False
            self.Y_dig = zeros((Ny, N2))
            self.Y_dig[:,0] = squeeze(dot(C, X0))

        FO = G - dot(self.Ke,C)
        prev_x = X0
        prev_x_tilde = X0_tilde
        for k in range(1,N2):
            self.v[k] = self.calc_v(k)
            term1 = dot(FO, prev_x_tilde)
            term2 = self.H*self.v[k]
            term3 = colwise(dot(self.Ke, self.Y_dig[:,k-1]))
            ## if term1.any() or term2.any() or term3.any():
            ##     Pdb().set_trace()
            next_x_tilde =  term1 + term2 + term3
            next_x = dot(G, prev_x) + self.H*self.v[k]
            if not use_exp_Y:
                self.Y_dig[:,k] = squeeze(dot(C, next_x))
            self.X_dig[:,k] = squeeze(next_x)
            self.X_tilde[:,k] = squeeze(next_x_tilde)
            prev_x = next_x
            prev_x_tilde = next_x_tilde

        #self.X_dig = X
        #self.X_tilde = X_tilde
        #self.Y_dig = Y
        #self.v = squeeze(self.E*u + dot(self.K, self.X_dig))
        #self.v = V
        return self.Y_dig


    def bodes_from_digital_lsim(self, attr='X_tilde'):
        mat = rowwise(getattr(self, attr))
        denom = fft(self.r)#self.r should be defined from running
                           #digital_lsim_with_obs
        self.lsim_bodes = []
        for row in mat:
            num = fft(row)
            comp = num/denom
            curbode = rwkbode.rwkbode(compin=comp)
            self.lsim_bodes.append(curbode)

        y_fft = fft(squeeze(self.Y_dig))
        y_comp = y_fft/denom
        self.lsim_Y_bode = rwkbode.rwkbode(compin=y_comp)
        return self.lsim_bodes


    def plot_digital_lsim_results(self, fi=1, clear=True, \
                                  legloc=5):
        ax = self.create_ax(fi=fi, clear=clear)
        self.make_t()
        ax.plot(self.t, self.r, label='$r$')
        for n, Yi in enumerate(self.Y_dig):
            if self.M > 1:
                label = '$y_{%i}$' % n
            else:
                label = '$y$'
            ax.plot(self.t, Yi, label=label)
        ax.plot(self.t, self.v, label='$v$')
        self.label_axis()
        self.ax.legend(loc=legloc)


    def lsim_from_exp_file_w_obs(self, filepath, fi=1, plot=True, \
                                 clear=True):
        self.load_exp_time_file(filepath)
        u = self.data_file.u
        t = self.data_file.t
        self.digital_lsim_with_obs(u)


    def _calc_next_P(self, P, Q, R):
        mat1 = dot(self.G.T, dot(P,self.G))
        mat2 = dot(self.G.T, dot(P,self.H))
        mat3 = R + dot(self.H.T, dot(P,self.H))
        mat3i = linalg.inv(mat3)
        mat4 = dot(self.H.T, dot(P,self.G))
        next_P = Q + mat1 - dot(mat2, dot(mat3i,mat4))
        return next_P


    def digital_Ricatti_P(self, Q, R, P0=None, \
                          eps=1e-10, max_N=1e5):
        if P0 is None:
            P = zeros((self.N, self.N))
        else:
            P = P0
        n = 0

        while n < max_N:
            next_P = self._calc_next_P(P, Q, R)
            P_diff = next_P - P
            P = next_P
            if abs(P_diff).max() < eps:
                break
            else:
                n += 1
        if n > max_N - 1:
            print('_calc_next_P failed to converge after %i iterations' % n)
            print('current max of P_diff = ' + str(abs(P_diff).max()))
        return P


    def find_digital_optimal_K(self, Q=None, R=1.0, P0=None):
        """Following the approach of Ogata 'Discrete-Time Control
        Systems', Second Edition MATLAB Program 8-2 on pages 590-591,
        find the steady-state optimal gain vector K based on the
        solution P of the Ricatti equation."""
        if Q is None:
            Q = eye(self.N)
        P = self.digital_Ricatti_P(Q, R, P0)
        mat1 = R + dot(self.H.T, dot(P,self.H))
        mat2 = dot(self.H.T, dot(P,self.G))
        mat1i = linalg.inv(mat1)
        K = dot(mat1i, mat2)
        self.Ricatti_P = P
        return K


class SFLR_Motor_Only_Model(SS_model, \
                            SFLR_TF_models.SFLR_Time_File_Mixin):
    def __init__(self, *args, **kwargs):
        SS_model.__init__(self, *args, **kwargs)
        self.accel = False
        self.E = 1.0


    def load_exp_time_file(self, filepath, \
                           col_map={0:'t', 1:'n', 2:'u', \
                                    3:'v', 4:'theta'}):
        SFLR_TF_models.SFLR_Time_File_Mixin.load_exp_time_file(self, \
                                                               filepath, \
                                                               col_map=col_map)

class SS_P_control_mixin(object):
    def calc_v(self, k):
        vtemp = self.r[k] - self.Y_dig[0,k-1]
        return self.sat(vtemp)


class SS_model_P_control(SS_P_control_mixin, SS_model):
    pass

class SFLR_Motor_Only_P_Control_Model(SS_model_P_control, SFLR_Motor_Only_Model):
    pass




def model_from_pickle(pklpath, model_class=SS_model):
    mydict = rwkmisc.LoadPickle(pklpath)
    ss = model_class(**mydict)
    return ss


def SS_model_from_pickle(pklpath):
    return model_from_pickle(pklpath, model_class=SS_model)


def SS_model_P_control_from_pickle(pklpath):
    return model_from_pickle(pklpath, model_class=SS_model_P_control)


def Motor_Only_P_control_model_from_pickle(pklpath):
    return model_from_pickle(pklpath, model_class=SFLR_Motor_Only_P_Control_Model)


def Motor_Only_model_from_pickle(pklpath):
    return model_from_pickle(pklpath, model_class=SFLR_Motor_Only_Model)


class CCF_SS_Model_from_poles_and_zeros(SS_model):
    def _vector_to_code_list(self, myvect, pat, \
                             eps=1.0e-10):
        mylist = []
        for n in range(self.N):
            curval = myvect[n]
            if abs(curval) < eps:
                curval = 0.0
            curline = pat % (n, curval)
            mylist.append(curline)
        return mylist

    def dump_params_to_file(self, pathout, a_prefix='a', \
                            b_prefix_list=['bth','ba'], \
                            fmt = '%0.12e', eps=1.0e-10):
        a_pat = a_prefix + '%i = ' + fmt
        avect = -1.0*self.A[self.N-1,:]
        alist = self._vector_to_code_list(avect, a_pat, \
                                          eps=eps)
        big_list = alist + ['']
        for r, label in enumerate(b_prefix_list):
            b_pat = label + '%i = ' + fmt
            bvect = self.C[r,:]
            blist = self._vector_to_code_list(bvect, b_pat, \
                                              eps=eps)
            big_list.extend(blist)
            big_list.append('')
        txt_mixin.dump(pathout, big_list)


    def _build_A(self):
        self.A = zeros((self.N,self.N))
        for i in range(self.N-1):
            self.A[i, i+1] = 1.0#ones on super-diagonal

        for i in range(self.N):
            self.A[-1,i] = -1.0*self.coeffs[self.N-i]


    def _build_C(self):
        myzeros = self.zeros_in
        if myzeros is None:
            self.C_coeffs = [array([1.0])]
            nr = 1
        else:
            if isscalar(myzeros[0]):#there is only one output
                self.C_coeffs = [poly(myzeros)]
                nr = 1
            else:
                #myzeros = array(myzeros)
                #nr, nc = myzeros.shape
                nr = len(myzeros)#assume a list of lists
                #self.C_coeffs = numpy.zeros((nr, self.N+1))
                clist = []
                for row in myzeros:
                    coeffs = poly(row)
                    clist.append(coeffs)
                self.C_coeffs = clist
        self.C = numpy.zeros((nr, self.N))
        for r, row in enumerate(self.C_coeffs):
            N = len(row)-1
            for i, val in enumerate(row):
                self.C[r,N-i] = val


    def __init__(self, poles, zeros=None):
        self.poles_in = poles
        self.zeros_in = zeros
        coeffs = poly(poles)
        an = coeffs[0]
        self.coeffs = coeffs/an#force an to equal 1
        self.N = len(coeffs) - 1
        self._build_A()
        self.B = numpy.zeros((self.N,1))
        self.B[-1] = 1.0
        self._build_C()
        self.C_hat = copy.copy(self.C)#save values before scaling self.C
        nrC, junk = self.C.shape
        self.M = nrC
        self.I = eye(self.N)



class SFLR_model_w_bodes:
    def find_opt(self, output, input):
        found = 0
        for opt in self.bode_opts:
            if (opt.output_label == output) and \
               (opt.input_label == input):
                found = 1
                return opt
        #if we got this far, we didn't find a match
        assert found, "Did not find a bode with output %s and input %s." % \
               (output, input)


    def find_bode(self, bodeopt):
        output = bodeopt.output_label
        input = bodeopt.input_label
        found = 0
        for bode in self.bodes:
            if (bode.output == output) and \
                   (bode.input == input):
                found = 1
                return bode
        #if we got this far, we didn't find a match
        assert found, "Did not find a bode with output %s and input %s." % \
               (output, input)


    def calc_bodes(self, f):
        if hasattr(self, 'theta_C_ind') and (self.theta_C_ind is not None):
            theta_C_ind = self.theta_C_ind
        else:
            theta_C_ind = 0
        if hasattr(self, 'accel_C_ind') and (self.accel_C_ind is not None):
            accel_C_ind = self.accel_C_ind
        else:
            accel_C_ind = 1

        comp_mat = self.calc_freq_resp(f)
        th_v_comp = comp_mat[theta_C_ind,:]
        a_v_comp = comp_mat[accel_C_ind,:]
        if not hasattr(self, 'bode_input'):
            self.bode_input = 'v'
        in_str = self.bode_input
        th_v_opts = self.find_opt('theta',in_str)
        self.th_v_bode = rwkbode.rwkbode(output='theta', \
                                         input=in_str, \
                                         compin=th_v_comp, \
                                         seedfreq=th_v_opts.seedfreq, \
                                         seedphase=th_v_opts.seedphase)
        self.th_v_bode.PhaseMassage(f)

        a_v_opts = self.find_opt('a',in_str)
        self.a_v_bode = rwkbode.rwkbode(output='a', \
                                        input=in_str, \
                                        compin=a_v_comp, \
                                        seedfreq=a_v_opts.seedfreq, \
                                        seedphase=a_v_opts.seedphase)
        self.a_v_bode.PhaseMassage(f)
        self.bodes = [self.th_v_bode, self.a_v_bode]

        return self.bodes



class CCF_SFLR_model_from_TF_coeffs(SS_P_control_mixin, \
                                    CCF_SS_Model_from_poles_and_zeros, \
                                    SFLR_model_w_bodes):
    def __init__(self, params_mod_name, mod_folder, \
                 b_prefixes=['bth','ba'], N=4, \
                 bode_opts=None):
        self.N = N
        self.load_coeffs(mod_folder, params_mod_name, b_prefixes)
        self.A = zeros((self.N,self.N))
        self._build_A()
        self.B = zeros((self.N, 1))
        self.B[-1,0] = 1.0
        self.M = len(b_prefixes)
        self.C = zeros((self.M, self.N))
        for r, prefix in enumerate(b_prefixes):
            attr = prefix + '_coeffs'
            vect = getattr(self, attr)
            self.C[r,:] = vect
        self.C_all = copy.copy(self.C)
        self.C = rowwise(self.C[0,:])
        self.use_dig = False
        self.bode_opts = bode_opts


    def load_coeffs(self, mod_folder, \
                    params_mod_name, \
                    b_prefixes):
        if mod_folder not in sys.path:
            sys.path.insert(1, mod_folder)
        TF_params = rwkmisc.my_import(params_mod_name)
        #import params_mod_name as TF_params
        reload(TF_params)
        self.TF_params = TF_params
        self.a_coeffs = self.build_coeff_list('a')
        for r, prefix in enumerate(b_prefixes):
            attr = prefix + '_coeffs'
            vect = self.build_coeff_list(prefix)
            setattr(self, attr, vect)


    def build_coeff_list(self, prefix):
        vect = zeros((self.N,))
        for n in range(self.N):
            attr = prefix + str(n)
            val = getattr(self.TF_params, attr)
            vect[n] = val
        return vect


    def _build_A(self):
        self.A = zeros((self.N,self.N))
        for i in range(self.N-1):
            self.A[i, i+1] = 1.0#ones on super-diagonal

        for i in range(self.N):
            self.A[-1,i] = -1.0*self.a_coeffs[i]


    def calc_bodes(self, f):
        comp_mat = self.calc_freq_resp(f)
        th_v_comp = comp_mat[0,:]
        #a_v_comp = comp_mat[1,:]
        th_v_opts = self.find_opt('theta','v')
        self.th_v_bode = rwkbode.rwkbode(output='theta', \
                                         input='v', \
                                         compin=th_v_comp, \
                                         seedfreq=th_v_opts.seedfreq, \
                                         seedphase=th_v_opts.seedphase)
        self.th_v_bode.PhaseMassage(f)

        ## a_v_opts = self.find_opt('a','v')
        ## self.a_v_bode = rwkbode.rwkbode(output='a', \
        ##                                 input='v', \
        ##                                 compin=a_v_comp, \
        ##                                 seedfreq=a_v_opts.seedfreq, \
        ##                                 seedphase=a_v_opts.seedphase)
        ## self.a_v_bode.PhaseMassage(f)
        ## self.bodes = [self.th_v_bode, self.a_v_bode]

        self.bodes = [self.th_v_bode]

        return self.bodes


class OCF_SFLR_model_from_TF_coeffs(CCF_SFLR_model_from_TF_coeffs):
    def __init__(self, params_mod_name, mod_folder, \
                 b_prefixes=['bth'], N=4, \
                 bode_opts=None):
        self.bode_opts = bode_opts
        self.N = N
        self.M = len(b_prefixes)
        self.load_coeffs(mod_folder, params_mod_name, b_prefixes)
        self._build_A()
        self._build_B()
        self._build_C()


    def _build_A(self):
        A = self._ones_on_super()
        for i in range(self.N):
            A[i,0] = -1.0*self.a_coeffs[self.N-1-i]
        self.A = A
        return self.A


    def _build_B(self, prefix='bth'):
        attr = prefix+'_coeffs'
        vect = getattr(self, attr)
        B = zeros((self.N,1))
        for i in range(self.N):
            B[i,0] = vect[self.N-1-i]
        self.B = B
        return self.B


    def _build_C(self):
        self.C = zeros((self.M, self.N))
        self.C[0,0] = 1.0
        return self.C


class non_P_control_mixin:
    def calc_v(self, k):
        vtemp = self.E*self.r[k] - squeeze(dot(self.K, self.X_tilde[:,k-1]))
        return self.sat(vtemp)


class CCF_from_TF_no_P_control(non_P_control_mixin, \
                               CCF_SFLR_model_from_TF_coeffs):
    pass


class OCF_from_TF_no_P_control(non_P_control_mixin, \
                               OCF_SFLR_model_from_TF_coeffs):
    pass


class SFLR_model_that_saves(rwkmisc.object_that_saves):
    def build_dict(self):
        mydict = {}
        for attr in self.saveattrs:
            if hasattr(self, attr):
                val = getattr(self, attr)
            else:
                val = None
            mydict[attr] = val
        return mydict


    def save(self, filepath, protocol=2):
        if not hasattr(self, 'saveattrs'):
            self.saveattrs = ['A','B','C','D',\
                              'bode_opts', \
                              'accel_C_ind', \
                              'theta_C_ind']

        rwkmisc.object_that_saves.save(self, \
                                       filepath, \
                                       protocol=protocol)


    def load(self, filepath):
        rwkmisc.object_that_saves.load(self, filepath)
        self._init2()


class SFLR_model_with_dig_obs(SFLR_model_that_saves):
    def save(self, filepath, protocol=2):
        if not hasattr(self, 'saveattrs'):
            self.saveattrs = ['A','B','C','D',\
                              'F','g','L','L_dig']

        rwkmisc.object_that_saves.save(self, \
                                       filepath, \
                                       protocol=protocol)


    def load(self, filepath):
        rwkmisc.object_that_saves.load(self, filepath)
        self._init2()


class SFLR_CCF_model(CCF_SS_Model_from_poles_and_zeros, \
                     SFLR_model_w_bodes, \
                     SFLR_model_that_saves):
    def __init__(self, poles, zeros=None, bode_opts=None):
        CCF_SS_Model_from_poles_and_zeros.__init__(self, poles, \
                                                   zeros=zeros)
        self.bode_opts = bode_opts


    def find_fit_bode(self, bode, fitbodes):
        for curfit in fitbodes:
            if curfit.input == bode.input and \
                   curfit.output == bode.output:
                return curfit


    def set_fit_parames(self, fitbodes, ffit):
        self.fitbodes = fitbodes
        self.ffit = ffit


    def get_fit_params(self, fitbodes=None, ffit=None):
        if fitbodes is None:
            fitbodes = self.fitbodes
        if ffit is None:
            ffit = self.ffit
        return fitbodes, ffit


    def find_C_gains(self, fitbodes=None, ffit=None):
        """Find the best gain to use for each output (i.e. a scaling
        factor for each row of C), by comparing the model bode
        magnitudes with the corresponding ones from fitbodes.

        This method assumes that bode[i] corresponds to row i of
        self.C"""
        fitbodes, ffit = self.get_fit_params(fitbodes, ffit)

        self.calc_bodes(ffit)
        for i, bode in enumerate(self.bodes):
            fit_bode = self.find_fit_bode(bode, fitbodes)
            ig = numpy.average(fit_bode.dBmag()/bode.dBmag())

            def my_cost(C):
                X = C[0]
                evect = fit_bode.dBmag()-(bode.dBmag()+20.0*log10(X))
                e = sum(evect**2)
                return e

            C_opt = optimize.fmin(my_cost, [ig])
            C_i = C_opt[0]
            self.C[i,:] *= C_i


    def phase_error_i(self, i, fitbodes=None, ffit=None, recalc=False):
        """Calculate the phase error for output i.  This method
        assumes that bode[i] corresponds to row i of self.C."""
        fitbodes, ffit = self.get_fit_params(fitbodes, ffit)

        if recalc:
            self.calc_bodes(ffit)
        bode = self.bodes[i]
        fit_bode = self.find_fit_bode(bode, fitbodes)
        phe_vect = fit_bode.phase - bode.phase
        e = sum(phe_vect**2)
        return e


    def phase_error(self, fitbodes=None, ffit=None, recalc=True):
        """Add up the phase error for each output.  This method
        assumes that bode[i] corresponds to row i of self.C"""
        fitbodes, ffit = self.get_fit_params(fitbodes, ffit)

        if recalc:
            self.calc_bodes(ffit)
        e = 0
        for i, bode in enumerate(self.bodes):
            e_i = self.phase_error_i(i, fitbodes, ffit, recalc=0)
            e += e_i
        return e



    def dBmag_error_i(self, i, fitbodes=None, ffit=None, recalc=False):
        """Calculate the phase error for output i.  This method
        assumes that bode[i] corresponds to row i of self.C."""
        fitbodes, ffit = self.get_fit_params(fitbodes, ffit)

        if recalc:
            self.calc_bodes(ffit)
        bode = self.bodes[i]
        fit_bode = self.find_fit_bode(bode, fitbodes)
        e_vect = fit_bode.dBmag() - bode.dBmag()
        e = sum(e_vect**2)
        return e


    def dBmag_error(self, fitbodes=None, ffit=None, recalc=True):
        """Add up the phase error for each output.  This method
        assumes that bode[i] corresponds to row i of self.C"""
        fitbodes, ffit = self.get_fit_params(fitbodes, ffit)

        if recalc:
            self.calc_bodes(ffit)
        e = 0
        for i, bode in enumerate(self.bodes):
            e_i = self.dBmag_error_i(i, fitbodes, ffit, recalc=0)
            e += e_i
        return e


    def total_error(self, fitbodes=None, ffit=None, recalc=True, \
                    phaseweight=0.02):
        mag_e = self.dBmag_error(fitbodes, ffit, recalc=recalc)
        phase_e = self.phase_error(fitbodes, ffit, recalc=0)
        total_e = mag_e + phase_e*phaseweight
        return total_e


    def check_C_signs(self, fitbodes=None, ffit=None):
        """Determine whether or not the phase error could be
        multipying each row of C by -1."""
        fitbodes, ffit = self.get_fit_params(fitbodes, ffit)

        for i, bode in enumerate(self.bodes):
            e1 = self.phase_error_i(i, fitbodes, ffit, recalc=1)
            self.C[i,:] *= -1.0
            e2 = self.phase_error_i(i, fitbodes, ffit, recalc=1)
            if e1 < e2:
                self.C[i,:] *= -1.0
        self.calc_bodes(ffit)


class SFLR_CCF_model_closed_loop(SFLR_CCF_model):
    def __init__(self, *args, **kwargs):
        SFLR_CCF_model.__init__(self, *args, **kwargs)
        self.bode_input = 'u'


class SFLR_CCF_model_theta_only(SFLR_CCF_model):
    def __init__(self, poles, \
                 theta_zeros=None, \
                 bode_opts=None):
        SFLR_CCF_model.__init__(self, poles, zeros=theta_zeros, \
                                bode_opts=bode_opts)
        self.C_theta = copy.copy(self.C)
        self.C_theta_hat = copy.copy(self.C)


class SFLR_CCF_model_w_accel_semi_Jordan(SFLR_CCF_model_theta_only):
    def __init__(self, poles, \
                 flex_pole1, \
                 flex_pole2, \
                 theta_zeros=None, \
                 B1=1.0, B2=-1.0, \
                 bode_opts=None):
        SFLR_CCF_model_theta_only.__init__(self, poles, \
                                           theta_zeros=theta_zeros, \
                                           bode_opts=bode_opts)
        self.flex_pole1 = flex_pole1
        self.flex_pole2 = flex_pole2
        self.B1 = B1
        self.B2 = B2
        self._build_accel_C()
        self.M = 2


    def find_theta_C_gain(self, fitbodes, ffit):
        """This method assumes that the theva/v bode is self.bodes[0]."""
        self.calc_bodes(ffit)
        i = 0
        bode = self.bodes[i]
        fit_bode = self.find_fit_bode(bode, fitbodes)
        ig = fit_bode.mag/bode.mag

        def my_cost(C):
            X = C[0]
            evect = fit_bode.mag-bode.mag*X
            return sum(evect**2)

        C_opt = optimize.fmin(my_cost, [ig])
        C_i = C_opt[0]
        self.C[i,:] *= C_i
        self.C_theta *= C_i
        return C_i


    def find_accel_C_row(self, fitbodes, ffit):
        """This method assumes that the accel/v bode is self.bodes[1]."""
        self.calc_bodes(ffit)
        i = 1
        bode = self.bodes[i]
        fit_bode = self.find_fit_bode(bode, fitbodes)
        ig = [self.B1, self.B2]

        def my_cost(C):
            self.B1 = C[0]
            self.B2 = C[1]
            self._build_accel_C()
            self.calc_bodes(ffit)
            evect = fit_bode.mag - self.bodes[i].mag
            return sum(evect**2)

        C_opt = optimize.fmin(my_cost, [ig])
        self.B1 = C_opt[0]
        self.B2 = C_opt[1]
        self._build_accel_C()
        return self.C[1,:]



    def find_C_gains(self, fitbodes, ffit):
        """Find the best gain to use for each output (i.e. a scaling
        factor for each row of C), by comparing the model bode
        magnitudes with the corresponding ones from fitbodes.

        This method assumes that bode[i] corresponds to row i of
        self.C"""
        self.calc_bodes(ffit)
        for i, bode in enumerate(self.bodes):
            fit_bode = self.find_fit_bode(bode, fitbodes)
            ig = fit_bode.mag/bode.mag

            def my_cost(C):
                X = C[0]
                evect = fit_bode.mag-bode.mag*X
                return sum(evect**2)

            C_opt = optimize.fmin(my_cost, [ig])
            C_i = C_opt[0]
            self.C[i,:] *= C_i


    def _build_accel_C(self):
        term1 = self.B1*self.flex_pole2
        term2 = self.B2*self.flex_pole1
        bvect = term1+term2
        self.a_term = copy.copy(bvect)
        self.a_zeros = roots(self.a_term)
        bvect = numpy.append(bvect, [0,0])
        C_accel = zeros((1, self.N))
        N = len(bvect)-1
        for i, val in enumerate(bvect):
            C_accel[0,N-i] = val
        self.C_accel = C_accel
        self.C = row_stack([self.C_theta, self.C_accel])
        return self.C


    def find_accel_C_row(self, fitbodes, ffit):
        """This method assumes that the accel/v bode is self.bodes[1]."""
        self.calc_bodes(ffit)
        i = 1
        bode = self.bodes[i]
        fit_bode = self.find_fit_bode(bode, fitbodes)
        ig = [self.B1, self.B2]

        def my_cost(C):
            self.B1 = C[0]
            self.B2 = C[1]
            self._build_accel_C()
            self.calc_bodes(ffit)
            #evect = fit_bode.mag - self.bodes[i].mag
            evect = fit_bode.dBmag() - self.bodes[i].dBmag()
            return sum(evect**2)

        C_opt = optimize.fmin(my_cost, [ig])
        self.B1 = C_opt[0]
        self.B2 = C_opt[1]
        self._build_accel_C()
        return self.C[1,:]


class SFLR_CCF_model_arb_accel(SFLR_CCF_model_w_accel_semi_Jordan):
    def __init__(self, poles, \
                 theta_zeros=None, \
                 C_ig = [], \
                 C_inds = [], \
                 bode_opts=None, \
                 phaseweight=0.0):
        SFLR_CCF_model_theta_only.__init__(self, poles, \
                                           theta_zeros=theta_zeros, \
                                           bode_opts=bode_opts)
        self.C_ig = C_ig
        self.C_accel = copy.copy(C_ig)
        self.C_inds = C_inds
        self._build_accel_C()
        self.M = 2
        self.phaseweight = phaseweight


    def _build_accel_C(self):
        C_accel = zeros((1, self.N))
        for ind, C_i in zip(self.C_inds, self.C_accel):
            C_accel[0,ind] = C_i
        self.C_accel = C_accel
        self.C = row_stack([self.C_theta, self.C_accel])
        return self.C


    def find_accel_C_row(self, fitbodes, ffit):
        """This method assumes that the accel/v bode is self.bodes[1]."""
        self.calc_bodes(ffit)
        i = 1
        bode = self.bodes[i]
        fit_bode = self.find_fit_bode(bode, fitbodes)
        ig = self.C_ig

        def my_cost(C):
            self.C_accel = C
            self._build_accel_C()
            self.calc_bodes(ffit)
            evect = fit_bode.dBmag() - self.bodes[i].dBmag()
            e = sum(evect**2)
            if self.phaseweight > 0.0:
                ephase = fit_bode.phase - self.bodes[i].phase
                e += sum(ephase**2)
            return e

        C_opt = optimize.fmin(my_cost, ig)
        self.C_accel = C_opt
        self._build_accel_C()
        return self.C[1,:]


class digital_SS_model(SS_model):
    def __init__(self, G, H, C, D=0.0, T=1.0/500):
        self.G = G
        self.H = H
        self.C = C
        self.D = D
        self.T = T
        nr,nc = G.shape
        self.N = nr
        self.use_dig = True


def digital_SS_model_from_pickle(pklpath):
    return model_from_pickle(pklpath, model_class=digital_SS_model)


class digital_SFLR_SS_model(digital_SS_model, \
                            SFLR_TF_models.SFLR_Time_File_Mixin_w_accel):
    def __init__(self, G, H, C, D=0.0, K=None, Ke=None, E=1.0, T=1.0/500):
        self.G = G
        self.H = H
        self.C = C
        self.D = D
        self.K = K
        self.Ke = Ke
        self.E = E
        self.T = T
        nr,nc = G.shape
        self.N = nr
        self.use_dig = True


    def lsim_from_exp_file(self, filepath, fi=1, plot=True, \
                           clear=True):
        self.load_exp_time_file(filepath)
        u = self.data_file.u
        t = self.data_file.t
        self.digital_lsim(u)


    def plot_time_domain_exp_vs_model(self, fi=1, legloc=5):
        self.create_ax(fi=fi)
        self.plot_exp_time_data()
        self.plot_model_data()
        self.label_time_domain_plot()
        self.ax.legend(loc=legloc)


    def digital_lsim(self, u, X0=None):
        SS_model.digital_lsim(self, u, X0=X0)
        self.theta = squeeze(self.Y_dig[0,:])
        self.accel = squeeze(self.Y_dig[1,:])


class digital_SFLR_SS_model_ignoring_accel(digital_SFLR_SS_model):
    def _get_C(self):
        """For now, I am handling multiple output systems by using
        just the first output.  So, this method just return self.C for
        systems with one output, but grabs the first row of C for
        multiple output systems."""
        temp = atleast_2d(self.C)[0,:]
        C = atleast_2d(temp)
        return C


    def calc_v(self, k):
        vtemp = self.E*self.r[k] - squeeze(dot(self.K, self.X_tilde[:,k]))
        return vtemp


    def plot_obs_terms_from_verification(self, fi=300, clear=True, \
                                         ind=0, \
                                         label_suffix=''):
        keep_going = True
        n = 1
        while keep_going:
            attr = 'term%i' % n
            label_attr = attr + label_suffix
            if hasattr(self, attr):
                ax = self._prep_ax(fi=fi, clear=clear)
                vect = getattr(self, attr)
                ax.plot(self.t, vect[ind,:], label=label_attr)
                ax.set_title(attr)
                self.label_axis()
                n += 1
                fi += 1
            else:
                keep_going = False
                break
        ax = self._prep_ax(fi=fi, clear=clear)
        #vect = self.term1[ind, :] + self.term2[ind,:] + self.term3[ind,:]
        vect = self.term1[ind, :] + self.term3[ind,:]
        ax.plot(self.t, vect, label='all terms '+label_suffix)
        ax.plot(self.t, self.X_tilde[ind, :], label='X_tile_' + str(ind) + ' ' + label_suffix)


        self.label_axis()


    def digital_lsim_with_obs(self, r, X0=None, X0_tilde=None):
        self.r = r
        if X0 is None:
            X0 = zeros((self.N,1))
        if X0_tilde is None:
            X0_tilde = zeros((self.N,1))
        N2 = len(r)
        V = zeros_like(r)
        X = zeros((self.N, N2))
        self.X_tilde = zeros((self.N, N2))
        X[:,0] = squeeze(X0)
        self.X_tilde[:,0] = squeeze(X0_tilde)
        self.term1 = zeros((self.N, N2))
        self.term2 = zeros((self.N, N2))
        self.term3 = zeros((self.N, N2))


        if hasattr(self, 'G_ol'):
            G = self.G_ol
        else:
            G = self.G

        if hasattr(self, 'H_ol'):
            H = self.H_ol
        else:
            H = self.H
        #C = self.C
        C = self._get_C()
        Ke = self.Ke
        K = self.K

        Ny, Nx = self.C.shape
        self.Y_dig = zeros((Ny, N2))
        self.Y_dig[:,0] = squeeze(dot(C, X0))

        FO = G - dot(Ke,C)
        prev_x = X0
        prev_x_tilde = X0_tilde
        debug_ind = 0
        for k in range(1,N2):
            #V[k] = self.E*r[k] - squeeze(dot(K, prev_x_tilde))
            V[k] = self.calc_v(k)
            self.term1[:,k] = squeeze(dot(FO, prev_x_tilde))
            self.term2[:,k] = squeeze(H*V[k-1])
            self.term3[:,k] = squeeze(colwise(dot(Ke, self.Y_dig[0,k-3].astype(int))))
            next_x_tilde =  self.term1[:,k] + self.term2[:,k] + self.term3[:,k]
            ## if term1.any() or term2.any() or term3.any():
            ##     if debug_ind < 2:
            ##         print('k = '+str(k))
            ##         print('term1 = ' + str(term1))
            ##         print('term2 = ' + str(term2))
            ##         print('term3 = ' + str(term3))
            ##         print('next_x_tilde = ' + str(next_x_tilde))
            ##         print('-'*20)
            ##         debug_ind += 1


            ## next_x_tilde = dot(FO, prev_x_tilde) + H*V[k] + \
            ##                colwise(dot(Ke, self.Y_dig[0,k-1]))
            next_x = dot(G, prev_x) + H*V[k]
            self.Y_dig[:,k] = squeeze(dot(C, next_x))
            X[:,k] = squeeze(next_x)
            self.X_tilde[:,k] = squeeze(next_x_tilde)
            prev_x = next_x
            prev_x_tilde = next_x_tilde

        self.X_dig = X
        #self.Y_dig = Y
        #self.v = squeeze(self.E*u + dot(self.K, self.X_dig))
        self.v = V
        self.theta = squeeze(self.Y_dig[0,:])
        return self.Y_dig


    def digital_lsim(self, u, X0=None):
        SS_model.digital_lsim(self, u, X0=X0)
        self.theta = squeeze(self.Y_dig[0,:])


    def plot_model_data(self):
        SFLR_TF_models.SFLR_Time_File_Mixin.plot_model_data(self, \
                                                            accel=False)
        ax = self.ax
        t = self.t
        ax.plot(t, self.v, label='$v_{model}$')

    ## def plot_model_data(self):
    ##     SFLR_TF_models.SFLR_Time_File_Mixin.plot_model_data(self, accel=False)


class digital_SFLR_SS_model_ignoring_accel_P_control(digital_SFLR_SS_model_ignoring_accel):
    def calc_v(self, k):
        vtemp = self.r[k] - self.Y_dig[0,k-1]
        return vtemp


def digital_SFLR_model_from_pickle(pklpath):
    return model_from_pickle(pklpath, model_class=digital_SFLR_SS_model)


def digital_SFLR_model_from_pickle_ignore_accel(pklpath):
    return model_from_pickle(pklpath, digital_SFLR_SS_model_ignoring_accel)


def digital_SFLR_model_from_pickle_ignore_accel_P_control(pklpath):
    return model_from_pickle(pklpath, digital_SFLR_SS_model_ignoring_accel_P_control)


class SFLR_SS_model(SFLR_model_w_bodes, \
                    SFLR_TF_models.SFLR_Time_File_Mixin_w_accel, \
                    SS_model):
    def plot_model_data(self):
        SFLR_TF_models.SFLR_Time_File_Mixin_w_accel.plot_model_data(self)
        ax = self.ax
        t = self.t
        ax.plot(t, self.v, label='$v_{model}$')


    def calc_num(self):
        p_act1 = self.p_act1
        p_act2 = self.p_act2
        K_act = self.K_act
        H = self.H
        s1 = 1.0*2.0j*pi#magnitude of s at 1 Hz - fix this point for
        m1 = abs(s1*(s1+p_act1)*(s1+p_act2))
        self.num = K_act*m1#*m2


    def calc_coeffs(self):
        #This code is auto-generated by
        #/home/ryan/siue/Research/PSoC_Research/SFLR_2010/modeling/numeric_TMM_model/ss_from_ad_hoc_TF.py
        #--------------------------------
        #--------------------------------
        self.a_hat1 = self.p_act1*self.p_act2*self.wp1**2*self.wp2**2*self.wz1**2*self.wz2**2
        self.a_hat2 = self.p_act1*self.wp1**2*self.wp2**2*self.wz1**2*self.wz2**2 + self.p_act2*self.wp1**2*self.wp2**2*self.wz1**2*self.wz2**2 + 2*self.p_act1*self.p_act2*self.wp1*self.zp1*self.wp2**2*self.wz1**2*self.wz2**2 + 2*self.p_act1*self.p_act2*self.wp2*self.zp2*self.wp1**2*self.wz1**2*self.wz2**2
        self.a_hat3 = self.wp1**2*self.wp2**2*self.wz1**2*self.wz2**2 + self.p_act1*self.p_act2*self.wp1**2*self.wz1**2*self.wz2**2 + self.p_act1*self.p_act2*self.wp2**2*self.wz1**2*self.wz2**2 + 2*self.p_act1*self.wp1*self.zp1*self.wp2**2*self.wz1**2*self.wz2**2 + 2*self.p_act1*self.wp2*self.zp2*self.wp1**2*self.wz1**2*self.wz2**2 + 2*self.p_act2*self.wp1*self.zp1*self.wp2**2*self.wz1**2*self.wz2**2 + 2*self.p_act2*self.wp2*self.zp2*self.wp1**2*self.wz1**2*self.wz2**2 + 4*self.p_act1*self.p_act2*self.wp1*self.wp2*self.zp1*self.zp2*self.wz1**2*self.wz2**2
        self.a_hat4 = self.p_act1*self.wp1**2*self.wz1**2*self.wz2**2 + self.p_act1*self.wp2**2*self.wz1**2*self.wz2**2 + self.p_act2*self.wp1**2*self.wz1**2*self.wz2**2 + self.p_act2*self.wp2**2*self.wz1**2*self.wz2**2 + 2*self.wp1*self.zp1*self.wp2**2*self.wz1**2*self.wz2**2 + 2*self.wp2*self.zp2*self.wp1**2*self.wz1**2*self.wz2**2 + 2*self.p_act1*self.p_act2*self.wp1*self.zp1*self.wz1**2*self.wz2**2 + 2*self.p_act1*self.p_act2*self.wp2*self.zp2*self.wz1**2*self.wz2**2 + 4*self.p_act1*self.wp1*self.wp2*self.zp1*self.zp2*self.wz1**2*self.wz2**2 + 4*self.p_act2*self.wp1*self.wp2*self.zp1*self.zp2*self.wz1**2*self.wz2**2
        self.a_hat5 = self.wp1**2*self.wz1**2*self.wz2**2 + self.wp2**2*self.wz1**2*self.wz2**2 + self.p_act1*self.p_act2*self.wz1**2*self.wz2**2 + 2*self.p_act1*self.wp1*self.zp1*self.wz1**2*self.wz2**2 + 2*self.p_act1*self.wp2*self.zp2*self.wz1**2*self.wz2**2 + 2*self.p_act2*self.wp1*self.zp1*self.wz1**2*self.wz2**2 + 2*self.p_act2*self.wp2*self.zp2*self.wz1**2*self.wz2**2 + 4*self.wp1*self.wp2*self.zp1*self.zp2*self.wz1**2*self.wz2**2
        self.a_hat6 = self.p_act1*self.wz1**2*self.wz2**2 + self.p_act2*self.wz1**2*self.wz2**2 + 2*self.wp1*self.zp1*self.wz1**2*self.wz2**2 + 2*self.wp2*self.zp2*self.wz1**2*self.wz2**2
        self.a_hat7 = self.wz1**2*self.wz2**2
        self.b_th_hat0 = self.H*self.num*self.wp1**2*self.wp2**2*self.wz1**2*self.wz2**2
        self.b_th_hat1 = 2*self.H*self.num*self.wz1*self.zz1*self.wp1**2*self.wp2**2*self.wz2**2 + 2*self.H*self.num*self.wz2*self.zz2*self.wp1**2*self.wp2**2*self.wz1**2
        self.b_th_hat2 = self.H*self.num*self.wp1**2*self.wp2**2*self.wz1**2 + self.H*self.num*self.wp1**2*self.wp2**2*self.wz2**2 + 4*self.H*self.num*self.wz1*self.wz2*self.zz1*self.zz2*self.wp1**2*self.wp2**2
        self.b_th_hat3 = 2*self.H*self.num*self.wz1*self.zz1*self.wp1**2*self.wp2**2 + 2*self.H*self.num*self.wz2*self.zz2*self.wp1**2*self.wp2**2
        self.b_th_hat4 = self.H*self.num*self.wp1**2*self.wp2**2
        self.b_a_hat2 = self.B1*self.H*self.a_gain*self.num*self.wp1**2*self.wp2**2*self.wz2**2 + self.B2*self.H*self.a_gain*self.num*self.wp1**2*self.wp2**2*self.wz1**2
        self.b_a_hat3 = 2*self.B1*self.H*self.a_gain*self.num*self.wz2*self.zz2*self.wp1**2*self.wp2**2 + 2*self.B2*self.H*self.a_gain*self.num*self.wz1*self.zz1*self.wp1**2*self.wp2**2
        self.b_a_hat4 = self.B1*self.H*self.a_gain*self.num*self.wp1**2*self.wp2**2 + self.B2*self.H*self.a_gain*self.num*self.wp1**2*self.wp2**2

        self.a1 = self.a_hat1/self.a_hat7
        self.a2 = self.a_hat2/self.a_hat7
        self.a3 = self.a_hat3/self.a_hat7
        self.a4 = self.a_hat4/self.a_hat7
        self.a5 = self.a_hat5/self.a_hat7
        self.a6 = self.a_hat6/self.a_hat7
        self.b_th0 = self.b_th_hat0/self.a_hat7
        self.b_th1 = self.b_th_hat1/self.a_hat7
        self.b_th2 = self.b_th_hat2/self.a_hat7
        self.b_th3 = self.b_th_hat3/self.a_hat7
        self.b_th4 = self.b_th_hat4/self.a_hat7
        self.b_a2 = self.b_a_hat2/self.a_hat7
        self.b_a3 = self.b_a_hat3/self.a_hat7
        self.b_a4 = self.b_a_hat4/self.a_hat7
        #--------------------------------


    def set_zero_coeffs(self, prop_str):
        for i in range(self.N):
            cur_prop = prop_str + str(i)
            if not hasattr(self, cur_prop):
                setattr(self, cur_prop, 0.0)


    def build_matrices(self):
        self.A = zeros((self.N, self.N), dtype='float64')
        self.B = zeros((self.N,), dtype='float64')
        self.C = zeros((2, self.N), dtype='float64')

        self.B[-1] = 1.0

        for i in range(1, self.N):
            self.A[i-1,i] = 1.0

        for i in range(self.N):
            a_i_str = 'a' + str(i)
            a_i = getattr(self, a_i_str)
            self.A[-1,i] = -a_i

            b_th_i_str = 'b_th' + str(i)
            b_th_i = getattr(self, b_th_i_str)
            self.C[0,i] = b_th_i

            b_a_i_str = 'b_a' + str(i)
            b_a_i = getattr(self, b_a_i_str)
            self.C[1,i] = b_a_i


    def create_LTI(self):
        B2 = rwkmisc.colwise(self.B)
        nr, nc = self.C.shape
        self.lti1 = signal.lti(self.A, B2, self.C[0,:], 0)
        self.lti2 = signal.lti(self.A, B2, self.C[1,:], 0)



    def __init__(self, pklname, bode_opts=None, N=7):
        mydict = rwkmisc.LoadPickle(pklname)
        for key, value in mydict.iteritems():
            setattr(self, key, value)
        self.calc_num()
        self.calc_coeffs()
        self.N = N
        self.set_zero_coeffs('a')
        self.set_zero_coeffs('b_th')
        self.set_zero_coeffs('b_a')
        self.build_matrices()
        self.I = eye(N)
        self.bode_opts = bode_opts
        self.create_LTI()


    ## def calc_freq_resp_one_s(self, s_i):
    ##     mat = s_i*self.I - self.A
    ##     mati = linalg.inv(mat)
    ##     M = dot(mati, self.B)
    ##     comp = dot(self.C,M)
    ##     return comp


    ## def calc_freq_resp(self, f):
    ##     svect = 2.0j*pi*f
    ##     comp_mat = zeros((2, len(svect)), dtype='complex128')

    ##     for i, s_i in enumerate(svect):
    ##         comp_i = self.calc_freq_resp_one_s(s_i)
    ##         comp_mat[:,i] = comp_i
    ##     self.comp_mat = comp_mat
    ##     return comp_mat



def right_of_penalty(pole, line=-10.0):
    penalty = 0.0
    if real(pole) > line:
        penalty = real(pole) - line
    return penalty


def pole_zero_sorter(pole_vect, r=6.0*2*pi, imag_tol=1e-5):
    filt_poles = [pole for pole in pole_vect if imag(pole) > -1e-7]#only positive poles
    small_poles = [pole for pole in filt_poles if abs(pole) < r]
    small_imag = [pole for pole in small_poles if \
                  abs(imag(pole)) > imag_tol]
    small_real = [pole for pole in small_poles if abs(imag(pole)) <= imag_tol]
    return small_imag, small_real


def find_zetas(pole_vect):
    z = abs(real(pole_vect))/(abs(pole_vect))
    return z


class SFLR_SS_model_ABCD(SFLR_model_w_bodes, \
                         SS_model, \
                         SFLR_model_that_saves):
    def __init__(self, A=None, B=None, C=None, D=0.0):
        SS_model.__init__(self, A, B, C, D=D)


class SFLR_SS_model_GHCD(SFLR_model_w_bodes, \
                         SS_model, \
                         SFLR_model_that_saves):
    def _init2(self):
        nr,nc = self.G.shape
        self.N = nr
        nrC, ncC = self.C.shape
        self.M = nrC
        self.I = eye(self.N)


    def __init__(self, G=None, H=None, C=None, D=0.0, T=1.0/500, \
                 **kwargs):
        """Note that kwargs is just used to ignore extra keyword
        arguments from child class __init__ methods."""
        self.G = G
        self.H = H
        self.C = C
        self.D = D
        self.T = T
        self.use_dig = True
        self.mysat = 200.0
        self.accel = True
        if self.G is not None:
            self._init2()



class SS_pole_optimizer_mixin(object):
    def find_CL_poles(self, K=None):
        if K is None:
            K = self.K
        if self.use_dig:
            self.calc_Ghat(K=K)
            self.CL_poles = linalg.eigvals(self.Ghat)
        else:
            self.calc_Ahat(K=K)
            self.CL_poles = linalg.eigvals(self.Ahat)
        return self.CL_poles


    def find_OL_poles(self, K=None):
        if self.use_dig:
            self.OL_poles = linalg.eigvals(self.G)
        else:
            self.OL_poles = linalg.eigvals(self.A)
        return self.OL_poles


    def pole_sorter(self, r=6.0*2*pi, imag_tol=1e-5):
        """Find the poles whose abs is less than r and sort them into real
        and imaginary based on whether their imaginary part is greater
        than or less than imag_tol."""
        self.find_CL_poles()
        small_imag, small_real = pole_zero_sorter(self.CL_poles, \
                                                  r=r, imag_tol=imag_tol)
        self.small_real = small_real
        self.small_imag = small_imag
        #self.small_real2 = self.filter_immovable(small_real)
        return small_imag, small_real


    def second_imag_penalty(self):
        penalty = 0.0
        imag_list = self.small_imag
        if len(imag_list) > 1:
            imag_vect = array(imag_list)
            #drop the right most element
            ind = real(imag_vect).argmax()
            imag_list2 = copy.copy(imag_list)
            imag_list2.pop(ind)
            imag_vect2 = array(imag_list2)
            zetas = find_zetas(imag_vect2)
            mybool = zetas < self.second_imag_zeta
            if mybool.any():
                penalty = 1.0
        return penalty


    def unstable_check(self):
        mybool = real(self.CL_poles) > 0.0
        return mybool.any()


    def pole_cost(self, verbosity=0):
        self.pole_sorter(r=20.0)
        small_imag = self.small_imag
        small_real = self.small_real
        my_small_poles = small_imag + small_real
        cost = 0.0
        small_imag_penalties = [right_of_penalty(pole, line=self.imagline) for pole \
                                in small_imag]
        small_real_penalties = [right_of_penalty(pole, line=self.realline) for pole \
                                in  small_real]
        imag_cost = sum(small_imag_penalties)
        real_cost = sum(small_real_penalties)
        cost += imag_cost*self.imagweight + real_cost
        if self.secondweight > 0.0:
            second_imag_pen = self.second_imag_penalty()*self.secondweight
            cost += second_imag_pen
        if self.unstable_check():
            cost += 1e3
        return cost


    def mycost(self, K, verbosity=0):
        K = atleast_2d(K)
        self.K = K
        if verbosity > 10:
            print('K = ' + str(K))
        pcost = self.pole_cost()
        cost = pcost
        neg_K_penalty = 0.0
        ## if (array(K) < 0.0).any():
        ##     neg_K_penalty = 1.0e3
        ##     cost += neg_K_penalty
        if verbosity > 0:
            print('pcost = %s' % pcost)
            #print('zcost = %s' % zcost)
            print('cost = %s' % cost)
            #print('pole locs = %s' % pole_vect)
            #print('zero locs = %s' % myzeros)
            print('neg_K_penalty = %s' % neg_K_penalty)
            print('-'*30)
        if self.logger is not None:
            self.logger.log(K, cost)
        return cost


def digital_pole_zero_sorter(pole_vect, dig_radius=0.125, imag_tol=1e-5):
    filt_poles = array([pole for pole in pole_vect if imag(pole) > -1e-7])#only positive poles
    #distance from 1.0+0.0j seems to be closely related to wn
    d_filt = abs(1.0-filt_poles)
    small_poles = []
    for pole, d in zip(filt_poles, d_filt):
        if d < dig_radius:
            small_poles.append(pole)
    small_imag = [pole for pole in small_poles if \
                  abs(imag(pole)) > imag_tol]
    small_real = [pole for pole in small_poles if abs(imag(pole)) <= imag_tol]
    return small_imag, small_real


def digital_omega_penalty(pole, radius=0.98):
    penalty = 0.0
    if abs(pole) > radius:
        penalty = abs(pole) - radius
    return penalty



class SS_digital_pole_optimizer_mixin(SS_pole_optimizer_mixin):
    """This class adapts the pole optimizer approach to a digital
    system."""
    def unstable_check(self):
        mybool = abs(self.CL_poles) > 1.0
        return mybool.any()


    def _find_dig_radius(self, wn):
        """Convert continous wn into a radius from 1.0+0.0j that
        defines small digital poles."""
        zeta = arange(0,1,0.01)
        wd = sqrt(1-zeta**2)*wn
        s = -wn*zeta + 1.0j*wd
        if hasattr(self,'T'):
            dt = self.T
        else:
            dt = 1.0/500.0
        z = exp(s*dt)
        d = abs(1.0-z)
        r = d.max()
        self.dig_radius = r
        return r


    def pole_sorter(self, dig_radius=None, imag_tol=1e-5):
        """Find the poles whose abs is less than r and sort them into real
        and imaginary based on whether their imaginary part is greater
        than or less than imag_tol."""
        self.find_CL_poles()#<-- this method calculates the CL poles
        if dig_radius is None:
            if hasattr(self, 'dig_radius'):
                dig_radius = self.dig_radius
            else:
                dig_radius = self._find_dig_radius(10.0*2*pi)
        small_imag, small_real = digital_pole_zero_sorter(self.CL_poles, \
                                                          dig_radius=dig_radius, \
                                                          imag_tol=imag_tol)#<-- this function doesn't exist yet
        self.small_real = small_real
        self.small_imag = small_imag
        #self.small_real2 = self.filter_immovable(small_real)
        return small_imag, small_real




    def pole_cost(self, verbosity=0):
        self.pole_sorter()
        small_imag = self.small_imag
        small_real = self.small_real
        my_small_poles = small_imag + small_real
        cost = 0.0
        small_imag_penalties = [digital_omega_penalty(pole, radius=self.imag_radius) for pole \
                                in small_imag]
        small_real_penalties = [digital_omega_penalty(pole, radius=self.real_radius) for pole \
                                in  small_real]
        imag_cost = sum(small_imag_penalties)
        real_cost = sum(small_real_penalties)
        cost += imag_cost*self.imagweight + real_cost
        if self.secondweight > 0.0:
            raise NotImplementedError, "self.second_imag_penalty not written yet for the digital case"
            second_imag_pen = self.second_imag_penalty()*self.secondweight
            cost += second_imag_pen
        if self.unstable_check():
            cost += 1e3
        return cost



class SFLR_SS_model_ABCD_pole_opt(SS_pole_optimizer_mixin, \
                                  SFLR_SS_model_ABCD):
    def __init__(self, *args, **kwargs):
        SFLR_SS_model_ABCD.__init__(self, *args, **kwargs)
        self.realline = -12.0
        self.imagline = -12.0
        self.imagweight = 10.0
        self.secondweight = 50.0
        self.second_imag_zeta = 0.95
        self.logger = None


class SFLR_SS_model_GHCD_digial_pole_opt(SS_digital_pole_optimizer_mixin, \
                                         SFLR_SS_model_GHCD):
    def __init__(self, *args, **kwargs):
        SFLR_SS_model_GHCD.__init__(self, *args, **kwargs)
        if not hasattr(self, 'dig_radius'):
            if kwargs.has_key('wn'):
                wn = kwargs['wn']
            else:
                wn = 10.0*2*pi
            dig_radius = self._find_dig_radius(wn)
        defaults = {'imag_radius':exp(-10.0*self.T), \
                    'real_radius':exp(-10.0*self.T), \
                    'imagweight':10.0, \
                    'secondweight':0.0, \
                    'logger':None, \
                   }

        for key, value in defaults.iteritems():
            if kwargs.has_key(key):
                curval = kwargs[key]
            else:
                curval = value
            setattr(self, key, curval)


class accel_ignorer(object):
    def _get_C(self):
        """For now, I am handling multiple output systems by using
        just the first output.  So, this method just return self.C for
        systems with one output, but grabs the first row of C for
        multiple output systems."""
        temp = atleast_2d(self.C)[0,:]
        C = atleast_2d(temp)
        return C


class SFLR_SS_model_GHCD_7_pole_JCF_opt(accel_ignorer, \
                                        SFLR_SS_model_GHCD_digial_pole_opt):
    def __init__(self, *args, **kwargs):
        SFLR_SS_model_GHCD_digial_pole_opt.__init__(self, *args, **kwargs)
        self.N_conj = 5
        self.N_real = 3
        self.imag_inds = [3,5]



    def run_observer(self, i):
        if hasattr(self, 'G_ol'):
            G = self.G_ol
        else:
            G = self.G
        if hasattr(self, 'H_ol'):
            H = self.H_ol
        else:
            H = self.H

        C = self._get_C()

        self.term1[:,i] = squeeze(colwise(dot(G, self.X_tilde[:,i-1])))
        self.term2[:,i] = squeeze(colwise(H*self.vvect[i-1]))
        Y_tilde_float_i = squeeze(dot(C, self.X_tilde[:,i-1]))
        self.term3[:,i] = squeeze(colwise(dot(self.Ke, self.Y_dig[0,i-1]-Y_tilde_float_i)))
        next_x_tilde = self.term1[:,i] + self.term2[:,i] + self.term3[:,i]
        self.X_tilde[:,i] = squeeze(next_x_tilde)


    def calc_v(self, i):
        #self.Yvect[0,i-1] = self.Y_dig[i-1]
        #self.Yvect[1,i] = self.avect[i]
        self.run_observer(i)
        vtemp = self.E*self.r[i] - dot(self.K, self.X_tilde[:,i])
        if vtemp.shape == (1,):
            vtemp = vtemp[0]
        #vtemp = self.uvect[i] - self.yvect[i-1]#P control with kp=1.0
            #for debugging
        if abs(self.Y_dig[0,i-1]) > 1000:
            vtemp = 0.0#instability software limit switch

        self.vvect[i] = self.sat(real(vtemp))
        return self.vvect[i]



    def _init_sim_vectors(self):
        dtype = complex128
        N2 = len(self.r)
        self.X_tilde = zeros((self.N, N2), dtype=dtype)
        self.term1 = zeros((self.N, N2), dtype=dtype)
        self.term2 = zeros((self.N, N2), dtype=dtype)
        self.term3 = zeros((self.N, N2), dtype=dtype)
        self.vvect = zeros((N2), dtype=float64)
        Ny, Nx = self.C.shape
        self.Y_dig = zeros((Ny, N2))



    def digital_lsim_with_obs(self, r, X0=None, X0_tilde=None):
        dtype = complex128
        self.dtype = dtype
        self.r = r
        self.N2 = len(r)

        if X0 is None:
            X0 = zeros((self.N,1), dtype=dtype)
        if X0_tilde is None:
            X0_tilde = zeros((self.N,1), dtype=dtype)
        self._init_sim_vectors()
        #V = zeros_like(r)
        X = zeros((self.N, self.N2), dtype=dtype)
        X[:,0] = squeeze(X0)
        self.X_tilde[:,0] = squeeze(X0_tilde)


        if hasattr(self, 'G_ol'):
            G = self.G_ol
        else:
            G = self.G

        if hasattr(self, 'H_ol'):
            H = self.H_ol
        else:
            H = self.H
        #C = self.C
        C = self._get_C()
        Ke = self.Ke
        K = self.K

        self.Y_dig[:,0] = squeeze(dot(C, X0))

        FO = G - dot(Ke,C)
        prev_x = X0
        prev_x_tilde = X0_tilde
        debug_ind = 0
        for k in range(1,self.N2):
            #V[k] = self.E*r[k] - squeeze(dot(K, prev_x_tilde))
            self.vvect[k] = self.calc_v(k)
            ## self.term1[:,k] = squeeze(dot(FO, prev_x_tilde))
            ## self.term2[:,k] = squeeze(H*V[k-1])
            ## self.term3[:,k] = squeeze(colwise(dot(Ke, self.Y_dig[0,k-3].astype(int))))
            ## next_x_tilde =  self.term1[:,k] + self.term2[:,k] + self.term3[:,k]
            ## if term1.any() or term2.any() or term3.any():
            ##     if debug_ind < 2:
            ##         print('k = '+str(k))
            ##         print('term1 = ' + str(term1))
            ##         print('term2 = ' + str(term2))
            ##         print('term3 = ' + str(term3))
            ##         print('next_x_tilde = ' + str(next_x_tilde))
            ##         print('-'*20)
            ##         debug_ind += 1


            ## next_x_tilde = dot(FO, prev_x_tilde) + H*V[k] + \
            ##                colwise(dot(Ke, self.Y_dig[0,k-1]))
            next_x = dot(G, prev_x) + H*self.vvect[k]
            self.Y_dig[:,k] = squeeze(dot(self.C, next_x))
            X[:,k] = squeeze(next_x)
            ##self.X_tilde[:,k] = squeeze(next_x_tilde)
            prev_x = next_x
            ##prev_x_tilde = next_x_tilde

        self.X_dig = X
        #self.Y_dig = Y
        #self.v = squeeze(self.E*u + dot(self.K, self.X_dig))
        self.v = self.vvect
        self.theta = squeeze(self.Y_dig[0,:])
        return self.Y_dig


    def build_cvect_from_K(self, K):
        K = colwise(K)
        cvect = zeros((1,self.N), dtype=float64)
        cvect[0,0:self.N_real] = K[0:self.N_real,0]
        for ind in self.imag_inds:
            cvect[0,ind] = real(K[ind])
            cvect[0,ind+1] = imag(K[ind])
        return cvect


    def extract_K_i_from_C0(self):
        C0 = self.C[0,:]
        return self.build_cvect_from_K(C0)


    def build_K(self, cvect):
        cvect = squeeze(cvect)
        K = zeros((1,self.N), dtype=complex128)
        K[0,0:self.N_real] = cvect[0:self.N_real]
        for ind in self.imag_inds:
            curcomp = cvect[ind] + 1.0j*cvect[ind+1]
            K[0,ind] = curcomp
            K[0,ind+1] = conj(curcomp)
        return K


    def mycost(self, cvect, verbosity=0):
        K = self.build_K(cvect)
        self.K = K
        if verbosity > 10:
            print('K = ' + str(K))
        pcost = self.pole_cost()
        cost = pcost
        neg_K_penalty = 0.0
        ## if (array(K) < 0.0).any():
        ##     neg_K_penalty = 1.0e3
        ##     cost += neg_K_penalty
        if verbosity > 0:
            print('pcost = %s' % pcost)
            #print('zcost = %s' % zcost)
            print('cost = %s' % cost)
            #print('pole locs = %s' % pole_vect)
            #print('zero locs = %s' % myzeros)
            print('neg_K_penalty = %s' % neg_K_penalty)
            print('-'*30)
        if self.logger is not None:
            self.logger.log(K, cost)
        return cost


class SFLR_SS_model_GHCD_7_pole_JCF_opt_w_integrator(SFLR_SS_model_GHCD_7_pole_JCF_opt):
    def _init_sim_vectors(self):
        SFLR_SS_model_GHCD_7_pole_JCF_opt._init_sim_vectors(self)
        self.evect = zeros(self.N2, dtype=self.dtype)
        self.esum = zeros(self.N2, dtype=self.dtype)


    def calc_v(self, i):
        #self.Yvect[0,i-1] = self.Y_dig[i-1]
        #self.Yvect[1,i] = self.avect[i]
        self.run_observer(i)
        self.evect[i] = self.r[i] - self.Y_dig[0,i-1]
        if i > 0:
            self.esum[i] = self.esum[i-1] + self.evect[i]
        vtemp = self.E*self.r[i] - dot(self.K, self.X_tilde[:,i]) + self.Ki*self.esum[i]
        if vtemp.shape == (1,):
            vtemp = vtemp[0]
        #vtemp = self.uvect[i] - self.yvect[i-1]#P control with kp=1.0
            #for debugging
        if abs(self.Y_dig[0,i-1]) > 1000:
            vtemp = 0.0#instability software limit switch

        self.vvect[i] = self.sat(real(vtemp))
        return self.vvect[i]


    def find_CL_poles(self, K=None, Ki=None):
        if K is None:
            K = self.K
        if Ki is None:
            Ki = self.Ki
        C = self._get_C()
        G = self.G
        H = self.H

        UL = G- dot(H,K)
        UR = H*Ki
        Upper = column_stack([UL,UR])

        LL = -dot(C,G)+dot(C, dot(H,K))
        LR = 1-Ki*dot(C,H)
        Lower = column_stack([LL,LR])

        big_mat = row_stack([Upper, Lower])

        eigs = linalg.eigvals(big_mat)
        self.CL_poles = eigs
        return eigs


    def build_cvect_from_K(self, K=None, Ki=None):
        if K is None:
            K = self.K
        if Ki is None:
            Ki = self.Ki
        cvect = SFLR_SS_model_GHCD_7_pole_JCF_opt.build_cvect_from_K(self, K)
        cvect = numpy.append(cvect, [Ki])
        return cvect



    def extract_K_i_from_C0(self):
        C0 = self.C[0,:]
        cvect = zeros((1,self.N+1), dtype=float64)
        cvect[0,0:self.N_real] = C0[0:self.N_real]
        for ind in self.imag_inds:
            cvect[0,ind] = real(C0[ind])
            cvect[0,ind+1] = imag(C0[ind])
        cvect[0,self.N] = self.Ki
        return cvect


    def build_K(self, cvect):
        K = zeros((1,self.N), dtype=complex128)
        K[0,0:self.N_real] = cvect[0:self.N_real]
        for ind in self.imag_inds:
            curcomp = cvect[ind] + 1.0j*cvect[ind+1]
            K[0,ind] = curcomp
            K[0,ind+1] = conj(curcomp)
        Ki = cvect[-1]
        return K, Ki


    def mycost(self, cvect, verbosity=0):
        K, Ki = self.build_K(cvect)
        self.K = K
        self.Ki = Ki
        if verbosity > 10:
            print('K = ' + str(K))
        pcost = self.pole_cost()
        cost = pcost
        neg_K_penalty = 0.0
        ## if (array(K) < 0.0).any():
        ##     neg_K_penalty = 1.0e3
        ##     cost += neg_K_penalty
        if verbosity > 0:
            print('pcost = %s' % pcost)
            #print('zcost = %s' % zcost)
            print('cost = %s' % cost)
            #print('pole locs = %s' % pole_vect)
            #print('zero locs = %s' % myzeros)
            print('neg_K_penalty = %s' % neg_K_penalty)
            print('-'*30)
        if self.logger is not None:
            self.logger.log(K, cost)
        return cost


    def __init__(self, *args, **kwargs):
        SFLR_SS_model_GHCD_7_pole_JCF_opt.__init__(self, *args, **kwargs)
        if kwargs.has_key('Ki'):
            self.Ki = kwargs['Ki']
        else:
            self.Ki = 0.1



class SFLR_SS_model_GHCD_7_pole_JCF_opt_w_integrator2(SFLR_SS_model_GHCD_7_pole_JCF_opt_w_integrator):
    def find_pole_nearest_to_1_0j(self):
        d = abs(1.0 - self.CL_poles)
        ind = d.argmin()
        return self.CL_poles[ind]



    def mycost(self, cvect, verbosity=0):
        K, Ki = self.build_K(cvect)
        self.K = K
        self.Ki = Ki
        if verbosity > 10:
            print('K = ' + str(K))
        pcost = self.pole_cost()
        cost = pcost
        neg_K_penalty = 0.0
        slowest_pole = self.find_pole_nearest_to_1_0j()
        slowest_cost = digital_omega_penalty(slowest_pole, radius=self.imag_radius)
        cost += slowest_cost*20#triple count the slowest one
        ## if (array(K) < 0.0).any():
        ##     neg_K_penalty = 1.0e3
        ##     cost += neg_K_penalty
        if verbosity > 0:
            print('pcost = %s' % pcost)
            #print('zcost = %s' % zcost)
            print('cost = %s' % cost)
            #print('pole locs = %s' % pole_vect)
            #print('zero locs = %s' % myzeros)
            print('neg_K_penalty = %s' % neg_K_penalty)
            print('-'*30)
        if self.logger is not None:
            self.logger.log(K, cost)
        return cost


class SFLR_SS_model_GHCD_7_pole_JCF_opt_K0_max(SFLR_SS_model_GHCD_7_pole_JCF_opt):
    def mycost(self, cvect, mygoal=1.5e12, \
               verbosity=0):
        cost = SFLR_SS_model_GHCD_7_pole_JCF_opt.mycost(self, \
                                                        cvect, \
                                                        verbosity=verbosity)
        K0_penalty =0.0
        K0 = cvect[0]
        if K0 < mygoal:
            K0_penalty = (mygoal-K0)/mygoal*0.1
        cost += K0_penalty
        return cost


class closed_loop_SS_model(SFLR_SS_model):
    def lsim2(self, U, T, X0=None, returnall=False, hmax=None):
        """Simulate output of a continuous-time linear system, using ODE solver.

        Inputs:

            system -- an instance of the LTI class or a tuple describing the
            system.  The following gives the number of elements in
            the tuple and the interpretation.
            2 (num, den)
            3 (zeros, poles, gain)
            4 (A, B, C, D)
        U -- an input array describing the input at each time T
            (linear interpolation is assumed between given times).
            If there are multiple inputs, then each column of the
            rank-2 array represents an input.
        T -- the time steps at which the input is defined and at which
            the output is desired.
        X0 -- (optional, default=0) the initial conditions on the state vector.

        Outputs: (T, yout, xout)

        T -- the time values for the output.
        yout -- the response of the system.
        xout -- the time-evolution of the state-vector.
        """
        # system is an lti system or a sequence
        #  with 2 (num, den)
        #       3 (zeros, poles, gain)
        #       4 (A, B, C, D)
        #  describing the system
        #  U is an input vector at times T
        #   if system describes multiple outputs
        #   then U can be a rank-2 array with the number of columns
        #   being the number of inputs

        # rather than use lsim, use direct integration and matrix-exponential.
        if hmax is None:
            hmax = T[1]-T[0]
        U = atleast_1d(U)
        T = atleast_1d(T)
        if len(U.shape) == 1:
            U = U.reshape((U.shape[0],1))
        sU = U.shape
        if len(T.shape) != 1:
            raise ValueError, "T must be a rank-1 array."
        if sU[0] != len(T):
            raise ValueError, "U must have the same number of rows as elements in T."
        if sU[1] != self.inputs:
            raise ValueError, "System does not define that many inputs."

        if X0 is None:
            X0 = zeros(self.B.shape[0],self.A.dtype)

        # for each output point directly integrate assume zero-order hold
        #   or linear interpolation.

        ufunc = interpolate.interp1d(T, U, kind='linear', axis=0, \
                                     bounds_error=False)

        def fprime(x, t, self, ufunc):
            return dot(self.A,x) + squeeze(dot(self.B,nan_to_num(ufunc([t]))))

        xout = integrate.odeint(fprime, X0, T, args=(self, ufunc), hmax=hmax)
        yout = dot(self.C,transpose(xout)) + dot(self.D,transpose(U))
        if returnall:
            return T, squeeze(transpose(yout)), xout
        else:
            return squeeze(transpose(yout))



    def lsim(self, u, t, interp=0, \
             returnall=False, X0=None, hmax=None):
        """Find the response of the TransferFunction to the input u
        with time vector t.  Uses signal.lsim.

        return y the response of the system."""
        try:
            out1 = signal.lsim(self.lti1, u, t, interp=interp, X0=X0)
            out2 = signal.lsim(self.lti2, u, t, interp=interp, X0=X0)
        except LinAlgError:
            #if the system has a pure integrator, lsim won't work.
            #Call lsim2.
            yout = self.lsim2(u, t, X0=X0, returnall=True, hmax=hmax)
            out1 = yout[:,0]
            out2 = yout[:,1]
                #override returnall because it is handled below
        self.theta = out1[1]
        self.accel = out2[1]
        self.x_theta = out1[2]
        self.x_accel = out2[2]

        self.v = squeeze(self.E*u + dot(self.K, self.x_theta.T))


        if returnall:#most users will just want the system output y,
            #but some will need the (t, y, x) tuple that
            #signal.lsim returns
            return out1, out2
        else:
            return out1[1], out2[1]


    def calc_bodes(self, f):
        comp_mat = self.calc_freq_resp(f)
        th_u_comp = comp_mat[0,:]
        a_u_comp = comp_mat[1,:]
        th_u_opts = self.find_opt('theta','u')
        self.th_u_bode = rwkbode.rwkbode(output='theta', \
                                         input='u', \
                                         compin=th_u_comp, \
                                         seedfreq=th_u_opts.seedfreq, \
                                         seedphase=th_u_opts.seedphase)
        self.th_u_bode.PhaseMassage(f)

        a_u_opts = self.find_opt('a','u')
        self.a_u_bode = rwkbode.rwkbode(output='a', \
                                        input='u', \
                                        compin=a_u_comp, \
                                        seedfreq=a_u_opts.seedfreq, \
                                        seedphase=a_u_opts.seedphase)
        self.a_u_bode.PhaseMassage(f)
        self.bodes = [self.th_u_bode, self.a_u_bode]



class discretized_closed_loop_SS_model(closed_loop_SS_model):
    def __init__(self, pklname, bode_opts=None, N=7, \
                 T=1.0/500):
        closed_loop_SS_model.__init__(self, pklname, \
                                      bode_opts=bode_opts, \
                                      N=N)
        self.use_dig = True
        self.discretize(T)


    def lsim_from_exp_file(self, filepath, fi=1, plot=True, \
                           clear=True):
        self.load_exp_time_file(filepath)
        u = self.data_file.u
        t = self.data_file.t
        self.digital_lsim(u)


    def digital_lsim(self, u, X0=None):
        closed_loop_SS_model.digital_lsim(self, u, X0=X0)
        self.theta = squeeze(self.Y_dig[0,:])
        self.accel = squeeze(self.Y_dig[1,:])



class discretized_closed_loop_SS_model_theta_only(discretized_closed_loop_SS_model):
    def __init__(self, pklname, bode_opts=None, N=7, \
                 T=1.0/500):
        discretized_closed_loop_SS_model.__init__(self, pklname, \
                                                  bode_opts=bode_opts, \
                                                  N=N, \
                                                  T=T)
        self.C = rowwise(self.C[0,:])


    def digital_lsim(self, u, X0=None):
        closed_loop_SS_model.digital_lsim(self, u, X0=X0)
        self.theta = squeeze(self.Y_dig[0,:])


    def plot_model_data(self):
        SFLR_TF_models.SFLR_Time_File_Mixin.plot_model_data(self, accel=False)
