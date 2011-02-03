import numpy
from numpy import arange, zeros, array, pi, log10

import rwkmisc
import SFLR_TMM

from IPython.Debugger import Pdb


class theta_fb_sys_contour(rwkmisc.object_that_saves):
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
        f2 = arange(mesh_change, 20, 0.5)
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
        im2 = arange(mesh_change, 20, 0.5)
        im = numpy.append(im1,im2)
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
        self.comp_mat = self.accel_u_bode_func(self.s, self.params, \
                                               *self.bode_args, \
                                               **self.bode_kwargs)


    def calc_a_theta_contours(self):
        self.a_theta_comp = self.accel_theta_bode_func(self.s, self.params, \
                                                       *self.bode_args, \
                                                       **self.bode_kwargs)
        self.theta_a_comp = 1.0/self.a_theta_comp


    def __init__(self, params_pkl_path=None, saved_pklname=None, \
                 accel_u_bode_func=None, accel_theta_bode_func=None,
                 bode_args=(), bode_kwargs={}):
        """params_pkl_path refers to the path for loading the SFLR
        parameters.  saved_pklname refers to the path to a pickle for
        restoring Theta FB contour data without re-running the
        calculations."""
        if saved_pklname is not None:
            self.load(saved_pklname)
        else:
            self.params_pkl_path = params_pkl_path
            self.load_params()            
            self.build_s()
            self.accel_u_bode_func = accel_u_bode_func
            self.accel_theta_bode_func = accel_theta_bode_func
            self.bode_args = bode_args
            self.bode_kwargs = bode_kwargs
            self.calc_theta_fb_contours()
            self.calc_a_theta_contours()
            self.abs_mat = abs(self.comp_mat)
            self.dB_mat = 20*log10(abs(self.comp_mat))
            self.bode_args = bode_args
            self.bode_kwargs = bode_kwargs
        self.saveattrs = ['f','im','s','comp_mat','a_theta_comp',\
                          'theta_a_comp','abs_mat','dB_mat']



class theta_fb_sys_contour_sympy(theta_fb_sys_contour):
    def calc_theta_fb_contours(self):
        self.theta_comp, self.comp_mat = self.accel_u_bode_func(self.s, self.params, \
                                                                *self.bode_args, \
                                                                **self.bode_kwargs)


    def calc_a_theta_contours(self):
        self.a_theta_comp = self.comp_mat/self.theta_comp
        self.theta_a_comp = 1.0/self.a_theta_comp
