from pylab import *
from scipy import *
from scipy import linalg, signal, integrate, optimize
import numpy

import controls

import rwkmisc, rwkbode
from rwkmisc import colwise, rowwise
from IPython.Debugger import Pdb

import copy

import SFLR_TF_models
reload(SFLR_TF_models)

import plotting_mixin

class SS_model(plotting_mixin.item_that_plots):
    def __init__(self, A, B, C, D=0.0):
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        nr,nc = A.shape
        self.N = nr
        nrC, ncC = C.shape
        self.M = nrC
        self.I = eye(self.N)
        self.use_dig = False


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



    def calc_freq_resp_one_s(self, s_i):
        mat = s_i*self.I - self.A
        mati = linalg.inv(mat)
        M = dot(mati, self.B)
        comp = dot(self.C,M)
        return comp


    def calc_freq_resp(self, f):
        svect = 2.0j*pi*f
        comp_mat = zeros((self.M, len(svect)), dtype='complex128')

        for i, s_i in enumerate(svect):
            comp_i = self.calc_freq_resp_one_s(s_i)
            comp_mat[:,i] = squeeze(comp_i)
        self.comp_mat = comp_mat
        return comp_mat


    def save_params(self, pklpath, attrlist):
        mydict = {}
        for key in attrlist:
            val = getattr(self, key)
            mydict[key] = val
        rwkmisc.SavePickle(mydict, pklpath)
        

    def discretize(self, T):
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


    def plot_poles(self, lt='bo', fi=1, clear=True, label=None):
        if not hasattr(self, 'ax'):
            self.create_ax(fi=fi, clear=clear)
        if not hasattr(self, 'poles'):
            self.find_poles()
        self.ax.plot(real(self.poles), imag(self.poles), \
                     lt, label=label)


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
        

    def observability(self, out_ind=0):
        A, B = self._get_A_and_B()
        C = self._get_C()
        #force using only first output for now
        nrC, ncC = C.shape
        if nrC > 1:
            C = atleast_2d(C[0,:])
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
              obs=False, out_ind=0):
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
            Q = self.observability()
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


    def digital_lsim(self, u, X0=None):
        if X0 is None:
            X0 = zeros((self.N,1))
        N2 = len(u)
        X = zeros((self.N, N2))
        X[:,0] = squeeze(X0)

        G = self.G
        H = self.H
        C = self.C

        Ny, Nx = self.C.shape
        Y = zeros((Ny, N2))
        Y[:,0] = squeeze(dot(C, X0))
        
        prev_x = X0
        for k in range(1,N2):
            next_x = dot(G, prev_x) + H*u[k-1]
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


    def digital_lsim_with_obs(self, r, X0=None, X0_tilde=None):
        if X0 is None:
            X0 = zeros((self.N,1))
        if X0_tilde is None:
            X0_tilde = zeros((self.N,1))
        N2 = len(r)
        V = zeros_like(r)
        X = zeros((self.N, N2))
        X_tilde = zeros((self.N, N2))
        X[:,0] = squeeze(X0)
        X_tilde[:,0] = squeeze(X0_tilde)

        G = self.G_ol
        H = self.H_ol
        #C = self.C
        C = self._get_C()
        Ke = self.Ke
        K = self.K

        Ny, Nx = self.C.shape
        Y = zeros((Ny, N2))
        Y[:,0] = squeeze(dot(C, X0))

        FO = G - dot(Ke,C)
        prev_x = X0
        prev_x_tilde = X0_tilde
        for k in range(1,N2):
            V[k] = self.E*r[k] - squeeze(dot(K, prev_x_tilde))
            next_x_tilde = dot(FO, prev_x_tilde) + H*V[k] + \
                           colwise(dot(Ke, Y[:,k-1]))
            next_x = dot(G, prev_x) + H*V[k]
            Y[:,k] = squeeze(dot(C, next_x))
            X[:,k] = squeeze(next_x)
            X_tilde[:,k] = squeeze(next_x_tilde)
            prev_x = next_x
            prev_x_tilde = next_x_tilde

        self.X_dig = X
        self.X_tilde = X_tilde
        self.Y_dig = Y
        #self.v = squeeze(self.E*u + dot(self.K, self.X_dig))
        self.v = V
        return self.Y_dig


    def lsim_from_exp_file_w_obs(self, filepath, fi=1, plot=True, \
                                 clear=True):
        self.load_exp_time_file(filepath)
        u = self.data_file.u
        t = self.data_file.t
        self.digital_lsim_with_obs(u)


def model_from_pickle(pklpath, model_class=SS_model):
    mydict = rwkmisc.LoadPickle(pklpath)
    ss = model_class(**mydict)
    return ss


def SS_model_from_pickle(pklpath):
    return model_from_pickle(pklpath, model_class=SS_model)
    

class CCF_SS_Model_from_poles_and_zeros(SS_model):
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
        comp_mat = self.calc_freq_resp(f)
        th_v_comp = comp_mat[0,:]
        a_v_comp = comp_mat[1,:]
        th_v_opts = self.find_opt('theta','v')
        self.th_v_bode = rwkbode.rwkbode(output='theta', \
                                         input='v', \
                                         compin=th_v_comp, \
                                         seedfreq=th_v_opts.seedfreq, \
                                         seedphase=th_v_opts.seedphase)
        self.th_v_bode.PhaseMassage(f)

        a_v_opts = self.find_opt('a','v')        
        self.a_v_bode = rwkbode.rwkbode(output='a', \
                                        input='v', \
                                        compin=a_v_comp, \
                                        seedfreq=a_v_opts.seedfreq, \
                                        seedphase=a_v_opts.seedphase)
        self.a_v_bode.PhaseMassage(f)
        self.bodes = [self.th_v_bode, self.a_v_bode]

        return self.bodes
    
    

class SFLR_CCF_model(CCF_SS_Model_from_poles_and_zeros, \
                     SFLR_model_w_bodes):
    def __init__(self, poles, zeros=None, bode_opts=None):
        CCF_SS_Model_from_poles_and_zeros.__init__(self, poles, \
                                                   zeros=zeros)
        self.bode_opts = bode_opts


    def find_fit_bode(self, bode, fitbodes):
        for curfit in fitbodes:
            if curfit.input == bode.input and \
                   curfit.output == bode.output:
                return curfit

        
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


    def digital_lsim_with_obs(self, r, X0=None, X0_tilde=None):
        if X0 is None:
            X0 = zeros((self.N,1))
        if X0_tilde is None:
            X0_tilde = zeros((self.N,1))
        N2 = len(r)
        V = zeros_like(r)
        X = zeros((self.N, N2))
        X_tilde = zeros((self.N, N2))
        X[:,0] = squeeze(X0)
        X_tilde[:,0] = squeeze(X0_tilde)

        G = self.G_ol
        H = self.H_ol
        #C = self.C
        C = self._get_C()
        Ke = self.Ke
        K = self.K

        Ny, Nx = self.C.shape
        Y = zeros((Ny, N2))
        Y[:,0] = squeeze(dot(C, X0))

        FO = G - dot(Ke,C)
        prev_x = X0
        prev_x_tilde = X0_tilde
        for k in range(1,N2):
            V[k] = self.E*r[k] - squeeze(dot(K, prev_x_tilde))
            next_x_tilde = dot(FO, prev_x_tilde) + H*V[k] + \
                           colwise(dot(Ke, Y[0,k-1]))
            next_x = dot(G, prev_x) + H*V[k]
            Y[:,k] = squeeze(dot(C, next_x))
            X[:,k] = squeeze(next_x)
            X_tilde[:,k] = squeeze(next_x_tilde)
            prev_x = next_x
            prev_x_tilde = next_x_tilde

        self.X_dig = X
        self.X_tilde = X_tilde
        self.Y_dig = Y
        #self.v = squeeze(self.E*u + dot(self.K, self.X_dig))
        self.v = V
        return self.Y_dig

    
def digital_SFLR_model_from_pickle(pklpath):
    return model_from_pickle(pklpath, model_class=digital_SFLR_SS_model)
    

def digital_SFLR_model_from_pickle_ignore_accel(pklpath):
    return model_from_pickle(pklpath, digital_SFLR_SS_model_ignoring_accel)


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
