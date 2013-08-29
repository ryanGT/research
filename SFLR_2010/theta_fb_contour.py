import numpy
from numpy import arange, zeros, array, pi, log10, angle

import rwkmisc
import SFLR_TMM

from IPython.core.debugger import Pdb
import accel_fb_system

class theta_fb_sys_contour(rwkmisc.object_that_saves, \
                           accel_fb_system.generic_contour_system):
    def load_params(self):
        mydict = rwkmisc.LoadPickle(self.params_pkl_path)
        self.params = SFLR_TMM.SFLR_params(**mydict)
        

    def build_s(self):
        #f = logspace(-1, 1.5, 500)
        #maxf = 30.0
        #maxf = 5.0
        #df = 0.1
        #f = arange(-20, maxf, df)
        #f = arange(-2.0, maxf, df)
        #sr = arange(1,5,0.01)
        #f = sr/(2*pi)
        mesh_change = 3.0
        #f0 = arange(-8, -2, 0.5)
        f0 = arange(-15, -2, 0.5)
        f1 = arange(-2, mesh_change, 0.01)
        f2 = arange(mesh_change, 50, 0.5)
        f_hat = numpy.append(f1,f2)
        f = numpy.append(f0,f_hat)
        #f = f1

        #maxi = 20.0
        #maxi = 1.0
        #di = 0.1
        #si = arange(-1,1,0.01)
        #im = arange(-1, maxi, di)
        #im = si/(2*pi)
        im1 = arange(-0.5, mesh_change, 0.01)
        im2 = arange(mesh_change, 15, 0.5)
        im3 = arange(15,50,0.01)
        #im = numpy.append(im1,im2)
        im = numpy.concatenate([im1,im2,im3])
        #im = im1
        nr = len(im)
        nc = len(f)
        s = zeros((nr,nc), dtype='D')

        for i in range(nr):
            for j in range(nc):
                s[i,j] = -2.0*pi*f[j] + 2.0j*pi*im[i]

        self.f = f
        self.im = im
        self.s = s


    def calc_theta_fb_contours(self):
        self.theta_comp_mat = self.theta_u_bode_func(self.s, self.params, \
                                                     *self.bode_args, \
                                                     **self.bode_kwargs)


    def calc_a_theta_contours(self):
        self.a_comp_mat = self.accel_u_bode_func(self.s, self.params, \
                                                 *self.bode_args, \
                                                 **self.bode_kwargs)
        #self.theta_a_comp = 1.0/self.a_theta_comp


    def plot_dB_contours(self, startfi=1, titles=None, \
                         myxlim=[-20,2], myylim=[-2,20], \
                         zoomin=False):

        if titles is None:
            titles = ['Closed-Loop Contour Plot for $\\theta/u$', \
                      'Closed-Loop Contour Plot for $\\ddot{x}/u$']

        self.plot_db_mat_contour(self.theta_dB_mat, startfi, titles[0], \
                                 myxlim=myxlim, myylim=myylim, \
                                 zoomin=zoomin)
        self.plot_db_mat_contour(self.a_dB_mat, startfi+1, titles[1], \
                                 myxlim=myxlim, myylim=myylim, \
                                 zoomin=zoomin)


    def __init__(self, theta_u_bode_func, accel_u_bode_func, \
                 params_pkl_path=None, saved_pklname=None, \
                 levels=None, \
                 bode_args=(), bode_kwargs={}):
        """params_pkl_path refers to the path for loading the SFLR
        parameters.  saved_pklname refers to the path to a pickle for
        restoring Theta FB contour data without re-running the
        calculations."""
        self.theta_u_bode_func = theta_u_bode_func
        self.accel_u_bode_func = accel_u_bode_func
        self.bode_args = bode_args
        self.bode_kwargs = bode_kwargs
        self.params_pkl_path = params_pkl_path
        self.levels = levels
        self.load_params()
        
        if saved_pklname is not None:
            self.load(saved_pklname)
        else:
            self.build_s()
            self.calc_theta_fb_contours()
            self.calc_a_theta_contours()
            self.theta_abs_mat = abs(self.theta_comp_mat)
            self.theta_dB_mat = 20*log10(self.theta_abs_mat)
            self.a_abs_mat = abs(self.a_comp_mat)
            self.a_dB_mat = 20*log10(self.a_abs_mat)
        self.saveattrs = ['f','im','s','theta_comp_mat','a_comp_mat',\
                          'theta_abs_mat','theta_dB_mat',\
                          'a_abs_mat','a_dB_mat',\
                          'levels']
        self.f_contour = self.f
        self.im_contour = self.im


    def find_poles(self, dB_mat=None):
        if dB_mat is None:
            dB_mat = self.a_dB_mat
        out = accel_fb_system.generic_contour_system.find_poles(self, dB_mat)
        return out


    def find_accel_poles(self):
        dB_mat = self.a_dB_mat
        self.a_poles = accel_fb_system.generic_contour_system.find_poles(self, dB_mat)
        return self.a_poles


    def find_theta_poles(self):
        dB_mat = self.theta_dB_mat
        self.theta_poles = accel_fb_system.generic_contour_system.find_poles(self, dB_mat)
        return self.theta_poles


    def find_accel_zeros(self):
        zeros_mat = 1.0/self.a_comp_mat
        self.a_zeros_db = accel_fb_system.comp_to_db(zeros_mat)
        self.a_zeros = self.find_poles(self.a_zeros_db)
        return self.a_zeros


    def find_theta_zeros(self):
        zeros_mat = 1.0/self.theta_comp_mat
        self.theta_zeros_db = accel_fb_system.comp_to_db(zeros_mat)
        self.theta_zeros = self.find_poles(self.theta_zeros_db)
        return self.theta_zeros


    def find_all_poles(self):
        self.all_poles = self.a_poles
        return self.all_poles


    def clean_poles_and_zeros(self):
        attrlist = ['theta_poles', 'a_poles', \
                    'a_zeros', 'theta_zeros']
        for attr in attrlist:
            self._clean_one_pole_or_zero_array(attr)

    def find_phase_near_zero(self, s=1.0e-5j):
        theta_comp = self.theta_u_bode_func(s, self.params)
        self.theta_phase_origin = angle(theta_comp,1)
        accel_comp = self.accel_u_bode_func(s, self.params)
        self.accel_phase_origin = angle(accel_comp,1)


    def find_origin_power(self, s=1e-5j):
        if not hasattr(self, 'theta_phase_origin'):
            self.find_phase_near_zero(s=s)
        self.theta_origin_power = numpy.round(self.theta_phase_origin/90.0)
        self.accel_origin_power = numpy.round(self.accel_phase_origin/90.0)
        self.min_origin_power = min([self.theta_origin_power, \
                                     self.accel_origin_power])

    def append_origin_zeros(self):
        self._append_origin_zeros('theta_origin_power', 'theta_zeros')
        self._append_origin_zeros('accel_origin_power', 'a_zeros')


    def build_SS_poles_and_zeros(self):
        self.SS_poles = self._build_full_list_of_poles_or_zeros_w_conj(self.all_poles)
        self.theta_SS_zeros = self._build_full_list_of_poles_or_zeros_w_conj(self.theta_zeros)
        self.accel_SS_zeros = self._build_full_list_of_poles_or_zeros_w_conj(self.a_zeros)
        self.SS_zeros = [self.theta_SS_zeros, self.accel_SS_zeros]
        #self.SS_zeros = [self.accel_SS_zeros]


class theta_fb_sys_contour_sympy(theta_fb_sys_contour):
    def calc_theta_fb_contours(self):
        self.theta_comp, self.comp_mat = self.accel_u_bode_func(self.s, self.params, \
                                                                *self.bode_args, \
                                                                **self.bode_kwargs)


    def calc_a_theta_contours(self):
        self.a_theta_comp = self.comp_mat/self.theta_comp
        self.theta_a_comp = 1.0/self.a_theta_comp


class theta_fb_contour_CND_paper(theta_fb_sys_contour):
    """I am investigating whether or not the LQG design for the
       ROM model of the SFLR will drive higher modes unstable.  To do this,
       I need to track the third and fourth order poles, so my contour
       matrix is different than for other purposes."""
    def build_s(self):
        #f = logspace(-1, 1.5, 500)
        #maxf = 30.0
        #maxf = 5.0
        #df = 0.1
        #f = arange(-20, maxf, df)
        #f = arange(-2.0, maxf, df)
        #sr = arange(1,5,0.01)
        #f = sr/(2*pi)
        mesh_change = 3.0
        #f0 = arange(-8, -2, 0.5)
        f0 = arange(-15, -2, 0.5)
        f1 = arange(-2, mesh_change, 0.01)
        f2 = arange(mesh_change, 30, 0.5)
        f_hat = numpy.append(f1,f2)
        f = numpy.append(f0,f_hat)
        #f = f1

        #maxi = 20.0
        #maxi = 1.0
        #di = 0.1
        #si = arange(-1,1,0.01)
        #im = arange(-1, maxi, di)
        #im = si/(2*pi)
        im1 = arange(-0.5, mesh_change, 0.01)
        im2 = arange(mesh_change, 15, 0.5)
        im3 = arange(15,150,0.1)
        #im = numpy.append(im1,im2)
        im = numpy.concatenate([im1,im2,im3])
        #im = im1
        nr = len(im)
        nc = len(f)
        s = zeros((nr,nc), dtype='D')

        for i in range(nr):
            for j in range(nc):
                s[i,j] = -2.0*pi*f[j] + 2.0j*pi*im[i]

        self.f = f
        self.im = im
        self.s = s


