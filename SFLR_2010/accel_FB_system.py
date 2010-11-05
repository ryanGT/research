from pylab import *
from scipy import *

import pylab as PL

import os, sys

import controls

import bode_plot_overlayer as BPO
import SFLR_bode_options

from rwkmisc import rowwise, my_import, LoadPickle

msg1 = """Cannot calculate bode response of a
Ga_ThetaFB_System unless self.thfb_a_comp
is defined."""

msg2 = """Cannot calculate the step response of
a Ga_ThetaFB_System unless self.ROM
is defined."""

def comp_to_db(comp_mat):
    dB_mat = 20*log10(abs(comp_mat))
    return dB_mat


def comp_to_zeros_db(comp_mat):
    zeros_mat = 1.0/comp_mat
    return comp_to_dB(zeros_mat)


class ga_theta_fb_system(object):
    def __init__(self, Ga, thfb_a_comp=None, \
                 ROM_model=None, substr=None, \
                 accel_comp_bode_opts=None, \
                 afb_bode_opts=None, \
                 tfb_contour_sys=None, \
                 levels=None):
        self.Ga = Ga
        self.ROM = ROM_model
        self.thfb_a_comp = thfb_a_comp
        self.substr = substr
        if accel_comp_bode_opts is None:
            accel_comp_bode_opts = SFLR_bode_options.Accel_Comp_Bode_opts[1]
        self.accel_comp_bode_opts = accel_comp_bode_opts
        if afb_bode_opts is None:
            afb_bode_opts = SFLR_bode_options.AccelFB_Bode_opts[1]
        self.afb_bode_opts = afb_bode_opts
        self.tfb_contour_sys = tfb_contour_sys
        if levels is None:
            levels = arange(-10, 50, 3)
        self.levels = levels
        if tfb_contour_sys is not None:
            self._map_from_tfb_contour_sys()
        if self.ROM is not None:
            self.update_rom(Ga)


    def _map_from_tfb_contour_sys(self):
        self.comp_mat = self.tfb_contour_sys.comp_mat
        self.s_contour = self.tfb_contour_sys.s
        self.f_contour = self.tfb_contour_sys.f
        self.im_contour = self.tfb_contour_sys.im


    def plot_dB_mat_contour(self, dB_afb, fi, title=None, \
                            myxlim=[-20,2],myylim=[-2,20],\
                            zoomin=True):
        figure(fi)
        clf()
        contour(-self.f_contour*2*pi, self.im_contour*2*pi , dB_afb, self.levels)
        if title is not None:
            PL.title(title)
        if zoomin:
            PL.xlim(myxlim)
            PL.ylim(myylim)


    def calc_ga_comp(self, s):
        Ga_comp = self.Ga(s)
        self.Ga_comp = Ga_comp


    def calc_ga_comp_from_f(self, f):
        s = 2.0j*pi*f
        self.calc_ga_comp(s)

        
    def update_rom(self, Ga=None):
        assert self.ROM is not None, msg2
        if Ga is not None:
            self.ROM.Ga = Ga
        self.ROM.build_TFs()
        self.theta_TF_simp = self.ROM.theta_TF.simplify()
        self.accel_TF_simp = self.ROM.accel_TF.simplify()


    def lsim(self, u, t, calca=True):
        assert self.ROM is not None, msg2
        self.t = t
        self.u = u
        self.theta_t = self.theta_TF_simp.lsim(u,t)
        if calca:
            self.accel_t = self.accel_TF_simp.lsim(u,t)
        else:
            self.accel_t = zeros_like(theta_t)
        return self.theta_t, self.accel_t


    def plot_time_domain(self, fig=None, fi=1, \
                         plotu=False, clear=False):
        if fig is None:
            fig = figure(1)
        if clear:
            fig.clf()
            
        if fig.get_axes():
            ax = fig.get_axes()[0]
        else:
            ax = fig.add_subplot(1,1,1)
        if plotu:
            ax.plot(self.t, self.u, label='$u$')
        plot(self.t, self.theta_t, label='$\\theta_{%s}$' % self.substr)
        plot(self.t, self.accel_t, label='$\\ddot{x}_{%s}$' % self.substr)        
        
        
    def calc_accel_comp_bode(self, f):
        assert self.thfb_a_comp is not None, msg1
        if not hasattr(self, 'Ga_comp'):
            self.calc_ga_comp_from_f(f)
        a_comp = self.Ga_comp*self.thfb_a_comp
        a_bode = BPO.comp_to_Bode(a_comp, f, \
                                  bode_opt=self.accel_comp_bode_opts)
        self.Accel_comp_Bode = a_bode


    def _set_kwargs(self, **kwargs):
        if not kwargs.has_key('label'):
            kwargs['label'] = self.substr
        return kwargs


    def plot_accel_comp_bode(self, f, fi=1, **kwargs):
        if not hasattr(self, 'Accel_comp_Bode'):
            self.calc_accel_comp_bode(f)
        kwargs = self._set_kwargs(**kwargs)
        BPO._plot_bode(self.Accel_comp_Bode, \
                       self.accel_comp_bode_opts, \
                       f, fignum=fi, **kwargs)
        

    def calc_afb_bode(self, f):
        assert self.thfb_a_comp is not None, msg1
        if not hasattr(self, 'Ga_comp'):
            self.calc_ga_comp_from_f(f)
        afb_a_comp = self.thfb_a_comp/(1.0+self.Ga_comp*self.thfb_a_comp)
        afb_a_bode = BPO.comp_to_Bode(afb_a_comp, f, \
                                      bode_opt=self.afb_bode_opts)
        self.Accel_AFB_Bode = afb_a_bode
        self.afb_a_comp = afb_a_comp


    def plot_afb_bode(self, f, fi=1, **kwargs):
        if not hasattr(self, 'Accel_AFB_Bode'):
            self.calc_afb_bode(f)
        kwargs = self._set_kwargs(**kwargs)    
        BPO._plot_bode(self.Accel_AFB_Bode, \
                       self.afb_bode_opts, \
                       f, fignum=fi, **kwargs)


    def nyquist_plot(self, f, fi=1, clear=False, **kwargs):
        if not hasattr(self, 'Ga_comp'):
            self.calc_ga_comp_from_f(f)
        Nyq = self.Ga_comp*self.thfb_a_comp
        figure(fi)
        if clear:
            clf()
        mirror = Nyq.conj()
        plot(Nyq.real, Nyq.imag)
        plot(mirror.real, mirror.imag)
        

    def calc_afb_compmat(self):
        s = self.s_contour
        Ga_comp = self.Ga(s)
        self.comp_afb = self.comp_mat/(1+Ga_comp*self.comp_mat)
        self.comp_afb_db = comp_to_db(self.comp_afb)
        return self.comp_afb


    def plot_afb_contour(self, fi, **kwargs):
        if not has_attr(self, 'comp_afb_db'):
            self.calc_afb_compmat()
        self.plot_dB_mat_contour(self.comp_afb_db, fi, **kwargs)



        

class Kp_Accel_ThetaFB_System(Ga_ThetaFB_System):
    def __init__(self, Ga, **kwargs):
        Ga = float(Ga)
        Ga_ThetaFB_System.__init__(self, Ga, **kwargs)


    def calc_Ga_comp(self, s):
        Ga_comp = self.Ga
        self.Ga_comp = Ga_comp
