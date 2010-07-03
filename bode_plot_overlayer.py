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

import controls

import copy

import SFLR_TMM
reload(SFLR_TMM)


def _plot_bode(bode, bode_opt, f, fignum=1, clear=False, \
               linetype='-', **kwargs):
    rwkbode.GenBodePlot(fignum, f, bode, clear=clear, \
                        linetype=linetype, **kwargs)
    PU.set_Bode_opts(fignum, bode_opt, coh=False)


def _PhaseMassage(bode, bode_opt, f):
    bode.seedphase = bode_opt.seedphase
    bode.seedfreq = bode_opt.seedfreq
    bode.PhaseMassage(f)


def tf_to_Bode(G, f, bode_opt, PhaseMassage=True):
    comp = G.FreqResp(f, fignum=None)
    bode = rwkbode.rwkbode(compin=comp, \
                           input=bode_opt.input_label, \
                           output=bode_opt.output_label)
    if PhaseMassage:
        _PhaseMassage(bode, bode_opt, f)
    return bode


def plot_bode_TMM(TMM_model, bode_opt, f, fignum=1, clear=False, \
                   PhaseMassage=False, **kwargs):
    if not hasattr(TMM_model, 'bodes'):
        print('calculating bodes')
        TMM_model.calc_bodes(f)
    bode = TMM_model.find_bode(bode_opt)
    if PhaseMassage:
        _PhaseMassage(bode, bode_opt, f)
    _plot_bode(bode, bode_opt, f, fignum=fignum, clear=clear, \
               **kwargs)


def plot_exp_bode(exp_mod, bode_opt, f=None, fignum=1, clear=False, \
                  trunc=True, PhaseMassage=False, **kwargs):
    if trunc:
        search_attr = 'trunc_avebodes'
        f = exp_mod.trunc_f
    else:
        search_attr = 'avebodes'
        f = exp_mod.f

    bode = exp_mod.find_bode(bode_opt.output_label, \
                             bode_opt.input_label, \
                             attr=search_attr)
    if PhaseMassage:
        _PhaseMassage(bode, bode_opt, f)
    _plot_bode(bode, bode_opt, f, fignum=fignum, clear=clear, \
               **kwargs)
    

def plot_different_exp_bodes_different_mods(exp_mods, bode_opts, \
                                            clear=False, **kwargs):
    first = 1
    for mod, bode_opt in zip(exp_mods, bode_opts):
        if first:
            clear = clear
            first = 0
        else:
            clear = False
        plot_exp_bode(mod, bode_opt, clear=clear, **kwargs)
    
    
def plot_multiple_TMM_models_one_Bode(model_list, bode_opt, f, \
                                      fignum=1, clear=True, \
                                      PhaseMassage=False, **kwargs):
    first = 1
    for model in model_list:
        if first:
            clear = clear
            first = 0
        else:
            clear = False
        plot_bode_TMM(model, bode_opt, f, fignum=fignum, \
                       clear=clear, PhaseMassage=PhaseMassage, **kwargs)


def plot_TMM_vs_exp_one_Bode(exp_mode, TMM_model, f, bode_opt, \
                             fignum=1, clear=True, \
                             PhaseMassage=False, \
                             labels=['exp.','TMM'], \
                             **kwargs):
    plot_exp_bode(exp_mode, bode_opt, clear=clear, \
                  fignum=fignum, PhaseMassage=PhaseMassage, \
                  label=labels[0], **kwargs)
    plot_bode_TMM(TMM_model, bode_opt, f, fignum=fignum, \
                  clear=False, PhaseMassage=PhaseMassage, \
                  label=labels[1], **kwargs)


def plot_TF_bode(G, f, bode_opt, PhaseMassage=True, \
                 **kwargs):
    bode = tf_to_Bode(G, f, bode_opt, PhaseMassage=PhaseMassage)
    _plot_bode(bode, bode_opt, f, **kwargs)
    

                             

def plot_different_TMM_models_different_Bodes(model_list, \
                                              bode_opt_list, \
                                              f, \
                                              fignum=1, clear=True, \
                                              PhaseMassage=False, \
                                              **kwargs):
    first = 1
    for model, bode_opt in zip(model_list, bode_opt_list):
        if first:
            clear = clear
            first = 0
        else:
            clear = False
        plot_bode_TMM(model, bode_opt, f, fignum=fignum, \
                       clear=clear, PhaseMassage=PhaseMassage, **kwargs)
        



class exp_bode_object(object):
    def __init__(self, modname, bode_opts, label='exp.'):
        self.mod = TDP.load_avebode_data_set(modname)
        self.bode_opts = bode_opts
        self.bode_attr = self.mod
        self.func = plot_exp_bode
        self.label = label


    def plot_bodes(self, f, startfi=1, clear=False, **kwargs):
        for i, opt in enumerate(self.bode_opts):
            fignum = startfi+i
            self.func(self.bode_attr, opt, f, \
                      fignum=fignum, clear=clear, \
                      **kwargs)


class TMM_bode_object(exp_bode_object):
    def __init__(self, TMM_model, bode_opts, label='TMM'):
        self.model = TMM_model
        self.bode_opts = bode_opts
        self.bode_attr = self.model
        self.func = plot_bode_TMM
        self.label = label


class OL_TF_bode_object(exp_bode_object):
    def __init__(self, G_th, G_a_th, bode_opts, label='TF'):
        self.G_th = G_th
        self.G_a_th = G_a_th
        self.G_a_v = G_th*G_a_th
        self.bode_opts = bode_opts
        self.func = plot_bode_TMM
        self.label = label


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


    def calc_bodes(self, f):
        th_v_comp = self.G_th.FreqResp(f, fignum=None)
        a_v_comp = self.G_a_v.FreqResp(f, fignum=None)
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


    def find_bode(self, output, input):
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
        if not hasattr(self, 'bodes'):
            self.calc_bodes(f)
        for i, opt in enumerate(self.bode_opts):
            fignum = startfi+i
            bode = self.find_bode(opt.output_label, opt.input_label)
            _plot_bode(bode, opt, f, fignum=fignum,
                       clear=clear, **kwargs)


class single_TF_bode_object(OL_TF_bode_object):
    def __init__(self, G, bode_opts, label='TF'):
        self.G = G
        self.bode_opts = bode_opts
        self.bode_opt = self.bode_opts[0]
        self.func = plot_bode_TMM
        self.label = label


    def calc_bodes(self, f):
        opt = self.bode_opt
        self.bode = tf_to_Bode(self.G, f, opt)
        self.bodes = [self.bode]

    
class Bode_Overlayer(object):
    def __init__(self, bode_obj_list, bode_opt_list):
        self.bode_obj_list = bode_obj_list
        self.bode_opt_list = bode_opt_list


    def plot_bodes(self, f, startfi=1, clear=True, **kwargs):
        first = 1
        for bode_obj in self.bode_obj_list:
            if first:
                clear = clear
                first = 0
            else:
                clear = False
            bode_obj.plot_bodes(f, startfi=startfi, clear=clear,
                                label=bode_obj.label, **kwargs)
            