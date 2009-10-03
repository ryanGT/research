from pylab import *
from scipy import *
from scipy import optimize

from rwkmisc import my_import

import time, copy

import SFLR_TMM
reload(SFLR_TMM)

import rwkbode

#import sympy_bodes
#reload(sympy_bodes)

mod_name = 'sympy_bodes_tau'
mod = my_import(mod_name)
reload(mod)

import exp_data

f2 = exp_data.f2
fexp = exp_data.f


def plot_exp(startfi=1):
    exp_data.plot_exp(startfi=startfi)
    

def calc_and_massage_a_th_Bode(f, th_bode, a_bode, \
                               **kwargs):
    a_th_bode = a_bode/th_bode
    a_th_bode.seedfreq = 10.0
    a_th_bode.seedphase = 0.0
    a_th_bode.PhaseMassage(f)
    a_th_bode = SFLR_TMM.a_massage(a_th_bode, f, **kwargs)
    return a_th_bode

def calc_and_massage_Bodes(f, params, func=mod.Bodes, \
                           **kwargs):
    s = 2.0j*pi*f
    th_comp, a_comp = func(s, params)
    a_bode = rwkbode.rwkbode(compin=a_comp)
    th_bode = rwkbode.rwkbode(compin=th_comp)
    
    a_bode.seedfreq = 2.0
    a_bode.seedphase = 30.0
    a_bode.PhaseMassage(f)
    a_bode = SFLR_TMM.a_massage(a_bode, f, **kwargs)
    th_bode.seedfreq = 1.2
    th_bode.seedphase = -90.0
    th_bode.PhaseMassage(f)

    return th_bode, a_bode


def plot_Bodes(f, th_bode, a_v_bode, a_th_bode, startfi=1, \
               clear=False, lt='-'):
    kwargs = {'clear':clear, 'linetype':lt}
    rwkbode.GenBodePlot(startfi, f, th_bode, **kwargs)
    rwkbode.GenBodePlot(startfi+1, f, a_v_bode, **kwargs)
    rwkbode.GenBodePlot(startfi+2, f, a_th_bode, **kwargs)

    
def calc_and_plot_Bodes(f, params, startfi=1, lt='-', \
                        func=mod.Bodes, clear=False, \
                        **kwargs):
    th_bode, a_bode = calc_and_massage_Bodes(f, params, func=func, \
                                             **kwargs)
    a_th_bode = calc_and_massage_a_th_Bode(f, th_bode, a_bode, \
                                           **kwargs)
    plot_Bodes(f, th_bode, a_bode, a_th_bode, startfi=startfi, \
               clear=clear, lt=lt)
    return th_bode, a_bode, a_th_bode


def run_simulation(pkl_name, mod_name, param_dict={},
                   f=f2, plot=True, startfi=1, **kwargs):
    params = SFLR_TMM.load_params('temp_params.pkl')
    params.update(param_dict)
    print('params.c_clamp = %s' % params.c_clamp)
    mymod = my_import(mod_name)
    myfunc = mymod.Bodes

    th_bode, a_bode = calc_and_massage_Bodes(f, params, func=myfunc)
    a_th_bode = calc_and_massage_a_th_Bode(f, th_bode, a_bode)
    if plot:
        #plot_exp(startfi=startfi)
        plot_Bodes(f, th_bode, a_bode, a_th_bode, startfi=startfi, \
                   **kwargs)
    return th_bode, a_bode, a_th_bode
    
