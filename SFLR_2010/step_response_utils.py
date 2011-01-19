from pylab import *
from scipy import *
import os, sys
import controls
import pylab_util as PU
import rwkos
from scipy import optimize
import bode_plot_overlayer as BPO
#reload(BPO)

import txt_data_processing as TDP
#reload(TDP)

import SFLR_TF_models

from mybodeopts import *

import bode_plots
#reload(bode_plots)

import add_design_dir

Accel_TMM_model = bode_plots.build_Accel_TMM_model()

import Gth
#reload(Gth)

import Ga5
#reload(Ga5)

z = 0.62831853071795862

def C4to5(C4):
    a2,a1,a0,gain = C4
    b1 = gain
    b0 = z*gain
    C5 = [a2,a1,a0,b1,b0]
    return C5


def Ga_to_C5(Ga):
    C5 = Ga.den.coeffs[1:].tolist() + Ga.num.coeffs.tolist()
    return C5


def Ga_to_C4(Ga):
    blist = Ga.num.coeffs.tolist()
    gain = blist[-2]
    z = blist[-1]/gain
    C4 = Ga.den.coeffs[1:].tolist() + [gain]
    return C4


def build_Ga(C):
    #a2,a1,a0,b1,b0,gain = C
    if len(C) == 4:
        a2,a1,a0,gain = C
        b1 = gain
        b0 = z*gain
        gain2 = 1.0
    elif len(C) == 5:
        a2,a1,a0,b1,b0 = C
        b2 = 0.0
        gain2 = 1.0
        ## b1 = gain
        ## b0 = z*gain
    ## elif len(C) == 6:
    ##     a2,a1,a0,b2,b1,b0 = C
    elif len(C) == 6:
        a2,a1,a0,b1,b0,gain2 = C
    Ga_opt = controls.TF([b1,b0],[1,a2,a1,a0])*gain2
    return Ga_opt


ROM_ATFB_model = SFLR_TF_models.G_th_G_a_TF(AccelFB_Bode_opts, \
                                            [], \
                                            Accel_TMM_model, \
                                            Gth=Gth.Gth, \
                                            Ga=Ga5.Ga5, \
                                            ffit=None, \
                                            label='ROM')

pkl_name = 'decent_ROM_params_09_09_10.pkl'
pkl_dir = '/home/ryan/git/research/SFLR_2010'
pkl_path = os.path.join(pkl_dir, pkl_name)

ROM_ATFB_model.load_params(pkl_path)
#ROM_ATFB_model.plot_one_bode(f, 1, fignum=1)

t = arange(0,2,0.01)
u = zeros_like(t)
u[50:] = 200.0


## th_Ga5_TF = ROM_ATFB_model.theta_TF.simplify()
## a_Ga5_TF = ROM_ATFB_model.accel_TF.simplify()

## th_Ga5 = th_Ga5_TF.lsim(u,t)
## a_Ga5 = a_Ga5_TF.lsim(u,t)

import measurement_utils
    

def step_response(C, calca=False):
    Ga_opt = build_Ga(C)
    ROM_ATFB_model.Ga = Ga_opt
    ROM_ATFB_model.build_TFs()
    t_s = ROM_ATFB_model.theta_TF.simplify()
    a_s = ROM_ATFB_model.accel_TF.simplify()

    th_a_s = t_s.lsim(u,t)
    if calca:
        y_a_s = a_s.lsim(u,t)
    else:
        y_a_s = zeros_like(th_a_s)
    return th_a_s, y_a_s


def mycost(C):
    th, accel = step_response(C)
    ts = measurement_utils.find_settling_time(th, u, t, p=0.01)
    mp = th.max()/u[-1]
    return ts + (mp-1)


if __name__ == '__main__':
    blist = Ga5.Ga5.num.coeffs.tolist()
    gain = blist[0]
    z = blist[1]/gain
    #ig = Ga5.Ga5.den.coeffs[1:].tolist() + [gain]
    ig = array([44.46407632, 258.87341724, 276.53309716, 453.43453413])
    #ig = Ga5.Ga5.den.coeffs[1:].tolist() + Ga5.Ga5.num.coeffs.tolist()
    #ig = Ga5.Ga5.den.coeffs[1:].tolist() + [0] + Ga5.Ga5.num.coeffs.tolist()
    #ig = Ga5.Ga5.den.coeffs[1:].tolist() + Ga5.Ga5.num.coeffs.tolist() + [1.01]
    ## ig = array([  1.63492735e+01,   3.10898809e+02,   9.44903871e+02,
    ##               2.48671857e+02,   3.13400963e+02,   6.62214709e-01])
    ## ig = array([  1.63492735e+01,   3.10898809e+02,   9.44903871e+02,
    ##               1.0,   240.0])


    th_Ga5, a_Ga5 = step_response(ig, calca=True)
    ts_Ga5 = measurement_utils.find_settling_time(th_Ga5, u, t, p=0.01)

    figure(1)
    clf()
    plot(t,u, label=None)
    plot(t,th_Ga5, label='$\\theta_{Ga5}$')
    plot(t,a_Ga5, label='$\\ddot{x}_{Ga5}$')



    C2 = array([  43.06911092,  317.7518689 ,  642.30362953,  250,
                     245])

    ## th2, a2 = step_response(C2, calca=True)
    ## plot(t,th2, label='$\\theta_{2}$')
    ## plot(t,a2, label='$\\ddot{x}_{2}$')
    ## ts2 = measurement_utils.find_settling_time(th2, u, t, p=0.01)

    C_opt = optimize.fmin(mycost, ig)
    th_opt, accel_opt = step_response(C_opt, calca=True)
    ts_opt = measurement_utils.find_settling_time(th_opt, u, t, p=0.01)
    plot(t,th_opt, label='$\\theta_{opt}$')
    plot(t,accel_opt, label='$\\ddot{x}_{opt}$')

    print('Ga5 t_s = %0.4f' % ts_Ga5)
    print('opt t_s = %0.4f' % ts_opt)

    show()
