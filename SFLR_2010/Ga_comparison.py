from pylab import *
from scipy import *

import os, sys

import controls

import SFLR_2010
reload(SFLR_2010)

import bode_plots
reload(bode_plots)

import add_design_dir

Accel_TMM_model = bode_plots.build_Accel_TMM_model()

import Gth
reload(Gth)

import SFLR_bode_options
reload(SFLR_bode_options)
AccelFB_bode_opts = SFLR_bode_options.AccelFB_Bode_opts
Accel_Comp_Bode_opts = SFLR_bode_options.Accel_Comp_Bode_opts

import Gth
reload(Gth)


import model_w_bm_bodes
reload(model_w_bm_bodes)

import model_w_bm_bodes_ThetaFB_editted
reload(model_w_bm_bodes_ThetaFB_editted)

import bode_plot_overlayer as BPO
reload(BPO)

from rwkmisc import rowwise, my_import, LoadPickle

import SFLR_TMM
reload(SFLR_TMM)

f = logspace(-1, 1.5, 300)

s = 2.0j*pi*f

pkl_name = 'model_w_bm_opt.pkl'
pkl_path = pkl_name

mydict = LoadPickle(pkl_path)
params = SFLR_TMM.SFLR_params(**mydict)

thfb_th_comp, thfb_a_comp = model_w_bm_bodes_ThetaFB_editted.Bodes(s, \
                                                                   params, \
                                                                   Gth.Gth)


def calc_AFB(Ga):
    Ga_comp = Ga(s)
    afb_a_comp = thfb_a_comp/(1.0+Ga_comp*thfb_a_comp)
    myopts = SFLR_bode_options.AccelFB_Bode_opts
    afb_a_bode = BPO.comp_to_Bode(afb_a_comp, f, \
                                  bode_opt=myopts[1])
    return afb_a_bode



def _plot_accel_Bode(bode, fi, linetype='-', \
                     bode_opts=None, **kwargs):
    if bode_opts is None:
        AccelFB_Bode_opts[1]
    BPO._plot_bode(bode, bode_opts, \
                   f, fignum=fi, linetype=linetype, \
                   **kwargs)


def plot_AFB(bode, fi, linetype='-', **kwargs):
    myopts = SFLR_bode_options.AccelFB_Bode_opts    
    _plot_accel_Bode(bode, fi, linetype, \
                     bode_opts=myopts[1], \
                     **kwargs)


def plot_Accel_comp(bode, fi, linetype='-', **kwargs):
    BPO._plot_bode(bode, Accel_Comp_Bode_opts[1], \
                   f, fignum=fi, linetype=linetype, **kwargs)


def calc_Accel_comp(Ga):
    Ga_comp = Ga(s)
    a_comp = Ga_comp*thfb_a_comp
    a_bode = BPO.comp_to_Bode(a_comp, f, \
                              bode_opt=Accel_Comp_Bode_opts[1])
    return a_bode



class AccelFB_System(object):
    def __init__(self, Ga, substr=''):
        self.Ga = Ga
        self.ROM = SFLR_2010.ROM_model.create_ROM_model(AccelFB_bode_opts, \
                                                        Accel_TMM_model, \
                                                        Gth.Gth, \
                                                        Ga)
        self.update_ROM(Ga)
        self.substr = substr

        
    def update_ROM(self, Ga=None):
        if Ga is not None:
            self.ROM.Ga = Ga
        self.ROM.build_TFs()
        self.theta_TF_simp = self.ROM.theta_TF.simplify()
        self.accel_TF_simp = self.ROM.accel_TF.simplify()


    def lsim(self, u, t, calca=True):
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
        
        
    def calc_Accel_comp(self):
        self.Accel_comp_Bode = calc_Accel_comp(self.Ga)


    def plot_Accel_comp(self, fi=1):
        if not hasattr(self, 'Accel_comp_Bode'):
            self.calc_Accel_comp()
        plot_Accel_comp(self.Accel_comp_Bode, fi, \
                        label=self.substr, linetype=None)


    def calc_AFB(self):
        self.Accel_AFB_Bode = calc_AFB(self.Ga)


    def plot_AFB(self, fi=1):
        if not hasattr(self, 'Accel_AFB_Bode'):
            self.calc_AFB()
        plot_AFB(self.Accel_AFB_Bode, fi, \
                 label=self.substr, linetype=None)
