"""This is a module for making Bode plots from transfer functions and
easily overlaying those Bodes with experimental data.  It should also
make curve fitting to experimental data very easy"""
from pylab import *
from scipy import *
from scipy import optimize

import pylab_util as PU
reload(PU)
import os, sys

import rwkos, rwkbode

import txt_data_processing as TDP
reload(TDP)

import rwkbode, pylab_util, rwkos
reload(rwkbode)

import controls
TF = controls.TransferFunction

import copy

import SFLR_TMM
reload(SFLR_TMM)

import bode_plot_overlayer as BPO
reload(BPO)


class DC_Motor_TF_Model(object):
    def find_bode(self, output, input):
        """This method is called by plot_bodes"""
        found = 0
        for bode in self.bodes:
            if (bode.output == output) and \
               (bode.input == input):
                found = 1
                return bode
        #if we got this far, we didn't find a match
        assert found, "Did not find a bode with output %s and input %s." % \
               (output, input)

    def plot_bodes(self, f, startfi=1, clear=False, **kwargs):
        """This method is for compatability with
        bode_plot_overlayer.py"""
        if not hasattr(self, 'bodes'):
            self.calc_bodes(f)
        if not kwargs.has_key('label'):
            kwargs['label'] = self.label
        for i, opt in enumerate(self.bode_opts):
            fignum = startfi+i
            bode = self.find_bode(opt.output_label, opt.input_label)
            BPO._plot_bode(bode, opt, f, fignum=fignum,
                           clear=clear, **kwargs)


    def build_act_iso(self):
        p_act = self.p_act1
        K_act = self.K_act
        s1 = 1.0*2.0j*pi#magnitude of s at 1 Hz - fix this point for
        m1 = abs(s1*(s1+p_act))
        num = K_act*m1#*m2
        self.G_act_iso = TF(num,[1,p_act,0])
        return self.G_act_iso


    def build_TFs(self):
        self.build_act_iso()
        self.G_act = self.G_act_iso


    def calc_bodes(self, f):
        bode = BPO.tf_to_Bode(self.G_act, f, \
                              self.bode_opts[0], PhaseMassage=True)
        self.bodes = [bode]
        self.G_act_bode = bode


    def negative_params_check(self):
        if not hasattr(self, 'pos_params'):
            self.pos_params = self.unknown_params
        for key in self.pos_params:
            val = getattr(self, key)
            if val < 0.0:
                return 10000.0
        return 0.0
    

    def act_cost(self, C):
        self.set_params(C)
        self.build_TFs()
        self.calc_bodes(self.ffit)
        cost = self.G_act_bode.dB_Phase_Error_sum(self.fit_bode)
        cost += self.negative_params_check()
        return cost


    def __init__(self, bode_opts, unknown_params=None, \
                 fit_bode=None, \
                 ffit=None, \
                 label='TF', \
                 params={'K_act':10.0, \
                         'p_act1':5*2*pi}):
        if ffit is None:
            ffit = logspace(-1,1.5,200)
        self.fit_bode = fit_bode
        self.ffit = ffit
        self.bode_opts = bode_opts
        if unknown_params is None:
            unknown_params = params.keys()
        self.unknown_params = unknown_params
        self.params = params
        for key, value in params.iteritems():
            setattr(self, key, value)

        self.label = label
        self.build_TFs()
        self.my_cost = self.act_cost
        

    def get_ig(self):
        """Build a vector of initial guesses for curve fitting based
        on the unknown parameters list self.unknown_params"""
        N = len(self.unknown_params)
        ig = zeros(N)
        for i, attr in enumerate(self.unknown_params):
            val = getattr(self, attr)
            ig[i] = val
        self.ig = ig
        return self.ig


    def set_params(self, C):
        """Set the parameters from a vector of optimization results,
        using self.unknown_params as the keys to where the values need
        to be assigned."""
        for attr, val in zip(self.unknown_params, C):
            setattr(self, attr, val)


    def fit(self):
        ig = self.get_ig()
        self.calc_TMM_bode()
        X_opt = optimize.fmin(self.my_cost, ig)
        self.X_opt = X_opt
        self.set_params(X_opt)
        self.build_TFs()
        self.calc_bodes(self.ffit)
        return X_opt
