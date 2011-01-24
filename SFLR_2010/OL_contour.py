import numpy
from numpy import arange, zeros, array, pi, log10

import rwkmisc
import SFLR_TMM

from IPython.Debugger import Pdb

import theta_fb_contour
import accel_fb_system

class OL_sys_contour(theta_fb_contour.theta_fb_sys_contour, \
                     accel_fb_system.generic_contour_system):
    def calc_theta_fb_contours(self):
        #Canceling out an inherited method that doesn't make sense
        raise NotImplementedError
        ## self.comp_mat = self.accel_u_bode_func(self.s, self.params, \
        ##                                        *self.bode_args, \
        ##                                        **self.bode_kwargs)


    def calc_a_theta_contours(self):
        #Canceling out an inherited method that doesn't make sense
        raise NotImplementedError
        ## self.a_theta_comp = self.accel_theta_bode_func(self.s, self.params, \
        ##                                                *self.bode_args, \
        ##                                                **self.bode_kwargs)
        ## self.theta_a_comp = 1.0/self.a_theta_comp


    def calc_contours(self):
        self.theta_comp_mat, self.accel_comp_mat = self.bode_func(self.s, self.params)
        


    def plot_dB_contours(self, startfi=1, titles=None, \
                         myxlim=[-20,2], myylim=[-2,20], \
                         zoomin=False):

        if titles is None:
            titles = ['Open-Loop Contour Plot for $\\theta/v$', \
                      'Open-Loop Contour Plot for $\\ddot{x}/v$']
        
        self.plot_db_mat_contour(self.theta_dB_mat, startfi, titles[0], \
                                 myxlim=myxlim, myylim=myylim, \
                                 zoomin=zoomin)
        self.plot_db_mat_contour(self.accel_dB_mat, startfi+1, titles[1], \
                                 myxlim=myxlim, myylim=myylim, \
                                 zoomin=zoomin)
        
                                 

    def __init__(self, mod, params_pkl_path=None, saved_pklname=None, \
                 levels=None):
        """params_pkl_path refers to the path for loading the SFLR
        parameters.  saved_pklname refers to the path to a pickle for
        restoring the contour data without re-running the
        calculations."""
        if saved_pklname is not None:
            self.load(saved_pklname)
        else:
            self.params_pkl_path = params_pkl_path
            self.load_params()            
            self.build_s()
            self.mod = mod
            self.bode_func = mod.Bodes
            self.calc_contours()
            self.theta_abs_mat = abs(self.theta_comp_mat)
            self.theta_dB_mat = 20*log10(abs(self.theta_comp_mat))
            self.accel_abs_mat = abs(self.accel_comp_mat)
            self.accel_dB_mat = 20*log10(abs(self.accel_comp_mat))
            self.levels = levels
        self.f_contour = self.f
        self.im_contour = self.im
        self.saveattrs = ['f','im','s', \
                          'theta_comp_mat','accel_comp_mat', \
                          'theta_abs_mat','accel_abs_mat', \
                          'theta_dB_mat','accel_dB_mat', \
                          'levels']
