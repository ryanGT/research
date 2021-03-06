import numpy
from numpy import arange, zeros, array, pi, log10, real, \
     imag, angle

import rwkmisc
import SFLR_TMM
import copy
from IPython.core.debugger import Pdb

import theta_fb_contour
import accel_fb_system

mysaveattrs = ['f','im','s', \
               'theta_comp_mat','accel_comp_mat', \
               'theta_abs_mat','accel_abs_mat', \
               'theta_dB_mat','accel_dB_mat', \
               'levels']

def find_ind_of_zearest(item, vect):
    mydiffs = vect - item
    absdiffs = abs(mydiffs)
    ind = absdiffs.argmin()
    return ind


def find_zearest(item, vect):
    ind = find_ind_of_zearest(item, vect)
    return vect[ind]


def pop_nearest_from_vector(item, vect):
    ind = find_ind_of_zearest(item, vect)
    if hasattr(vect, 'tolist'):
        mylist = vect.tolist()
    else:
        mylist = vect
    mylist.pop(ind)
    array_out = array(mylist)
    return array_out


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



    def _initial_calcs(self):
        self.build_s()
        self.calc_contours()
        self.theta_abs_mat = abs(self.theta_comp_mat)
        self.theta_dB_mat = 20*log10(abs(self.theta_comp_mat))
        self.accel_abs_mat = abs(self.accel_comp_mat)
        self.accel_dB_mat = 20*log10(abs(self.accel_comp_mat))


    def __init__(self, mod, params_pkl_path=None, saved_pklname=None, \
                 levels=None, tol=1e-4):
        """params_pkl_path refers to the path for loading the SFLR
        parameters.  saved_pklname refers to the path to a pickle for
        restoring the contour data without re-running the
        calculations."""
        if saved_pklname is not None:
            self.load(saved_pklname)
        else:
            self.params_pkl_path = params_pkl_path
            self.load_params()
            self.mod = mod
            self.bode_func = mod.Bodes
            self.levels = levels
            self._initial_calcs()
        self.tol = tol
        self.f_contour = self.f
        self.im_contour = self.im
        self.saveattrs = mysaveattrs


    def find_theta_OL_poles(self):
        self.theta_poles = self.find_poles(self.theta_dB_mat)
        return self.theta_poles


    def find_accel_OL_poles(self):
        self.accel_poles = self.find_poles(self.accel_dB_mat)
        return self.accel_poles


    def find_theta_OL_zeros(self):
        if not hasattr(self, 'theta_dB_mat'):
            self._initial_calcs()
        zeros_mat = 1.0/self.theta_comp_mat
        self.theta_zeros_db = accel_fb_system.comp_to_db(zeros_mat)
        self.theta_OL_zeros = self.find_poles(self.theta_zeros_db)
        return self.theta_OL_zeros


    def find_accel_OL_zeros(self):
        if not hasattr(self, 'accel_dB_mat'):
            self._initial_calcs()
        zeros_mat = 1.0/self.accel_comp_mat
        self.accel_zeros_db = accel_fb_system.comp_to_db(zeros_mat)
        self.accel_OL_zeros = self.find_poles(self.accel_zeros_db)
        return self.accel_OL_zeros


    def clean_poles_and_zeros(self):
        attrlist = ['theta_poles', 'accel_poles', \
                    'accel_OL_zeros', 'theta_OL_zeros']
        for attr in attrlist:
            self._clean_one_pole_or_zero_array(attr)


    def pole_filter_accel(self, pole):
        return self.pole_filter(pole, attr='accel_poles')


    def find_all_poles(self):
        if not hasattr(self, 'theta_poles'):
            self.find_theta_OL_poles()
        if not hasattr(self, 'accel_poles'):
            self.find_accel_OL_poles()
        self.all_poles = self.theta_poles
        self.unique_accel_poles = filter(self.pole_filter, self.accel_poles)
        self.unique_theta_poles = filter(self.pole_filter_accel, self.theta_poles)
        self.all_poles = numpy.append(self.all_poles, self.unique_accel_poles)
        return self.all_poles
    
    def build_SS_poles_and_zeros(self):
        self.SS_poles = self._build_full_list_of_poles_or_zeros_w_conj(self.all_poles)
        self.theta_SS_zeros = self._build_full_list_of_poles_or_zeros_w_conj(self.theta_OL_zeros)
        self.accel_SS_zeros = self._build_full_list_of_poles_or_zeros_w_conj(self.accel_OL_zeros)
        self.SS_zeros = [self.theta_SS_zeros, self.accel_SS_zeros]


    def append_origin_zeros(self):
        self._append_origin_zeros('theta_origin_power', 'theta_OL_zeros')
        self._append_origin_zeros('accel_origin_power', 'accel_OL_zeros')


    def find_nearest_pole(self, loc):
        mydiffs = self.all_poles - loc
        myind = abs(mydiffs).argmin()
        return self.all_poles[myind]


    def find_matching_zeros(self, tol=2.0):
        """Find zeros that show up in both accel and theta.  These are
        candidates for removal."""
        #ind_list = []
        zero_list = []
        for i, zero in enumerate(self.theta_OL_zeros):
            #ind = find_ind_of_zearest(zero, self.accel_OL_zeros)
            closest = find_zearest(zero, self.accel_OL_zeros)
            pdiff = abs(zero-closest)/abs(zero)*100.0
            if pdiff < tol:
                #ind_list.append(ind)
                zero_list.append(closest)
        return zero_list


two_func_attrs = copy.copy(mysaveattrs)
#two_func_attrs.append('theta_bode_func')
#two_func_attrs.append('accel_bode_func')

class OL_sys_contour_two_funcs(OL_sys_contour):
    def calc_contours(self):
        self.theta_comp_mat = self.theta_bode_func(self.s, self.params)
        self.accel_comp_mat = self.accel_bode_func(self.s, self.params)


    def __init__(self, theta_bode_func, accel_bode_func, \
                 params_pkl_path=None, saved_pklname=None, \
                 levels=None, tol=1e-4):
        """params_pkl_path refers to the path for loading the SFLR
        parameters.  saved_pklname refers to the path to a pickle for
        restoring the contour data without re-running the
        calculations."""
        self.params_pkl_path = params_pkl_path
        self.load_params()
        self.theta_bode_func = theta_bode_func
        self.accel_bode_func = accel_bode_func
        if saved_pklname is not None:
            self.load(saved_pklname)
        else:
            self.levels = levels
            self._initial_calcs()
        self.tol = tol
        self.f_contour = self.f
        self.im_contour = self.im
        self.saveattrs = mysaveattrs


    def find_phase_near_zero(self, s=1.0e-5j):
        theta_comp = self.theta_bode_func(s, self.params)
        self.theta_phase_origin = angle(theta_comp,1)
        accel_comp = self.accel_bode_func(s, self.params)
        self.accel_phase_origin = angle(accel_comp,1)


    def find_origin_power(self, s=1e-5j):
        if not hasattr(self, 'theta_phase_origin'):
            self.find_phase_near_zero(s=s)
        self.theta_origin_power = numpy.round(self.theta_phase_origin/90.0)
        self.accel_origin_power = numpy.round(self.accel_phase_origin/90.0)
        self.min_origin_power = min([self.theta_origin_power, \
                                     self.accel_origin_power])



