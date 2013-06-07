from pylab import *
from scipy import *

import pylab_util as PU

import controls

import rwkos, rwkdataproc, rwkbode

from rwkmisc import rowwise, my_import, LoadPickle

from IPython.core.debugger import Pdb

import sys, os, copy

import SFLR_SS_model as SS
reload(SS)

import bode_plot_overlayer as BPO
reload(BPO)

import model_w_bm_bodes_ThetaFB_editted
reload(model_w_bm_bodes_ThetaFB_editted)

import SFLR_TMM
reload(SFLR_TMM)

import SFLR_bode_options
OL_Bode_opts = [SFLR_bode_options.bode_opt_th_v, \
                SFLR_bode_options.bode_opt_a_v, \
               ]

ThetaFB_Bode_opts = [SFLR_bode_options.bode_opt_th_u, \
                     SFLR_bode_options.bode_opt_a_u, \
                    ]


AccelFB_Bode_opts = copy.copy(ThetaFB_Bode_opts)
AccelFB_Bode_opts[0].phaselim = [-200,100]
AccelFB_Bode_opts[0].seedfreq = 5.0
AccelFB_Bode_opts[0].seedphase = 0.0

AccelFB_Bode_opts[1].seedfreq = 5.0
AccelFB_Bode_opts[1].seedphase = 0.0


import SFLR_data_processing
#reload(SFLR_data_processing)


#Load Exp Data




#sym_dir = rwkos.FindFullPath('siue/Research/PSoC_Research/SFLR_2010/modeling/JVC_sympy_TMM')

#pkl_dir = rwkos.FindFullPath('SFLR_2010/vibration_suppression_design/TMM')
#pkl_name = 'model_w_bm_opt_9_bit.pkl'
#pkl_path = os.path.join(pkl_dir, pkl_name)
pkl_name = 'model_w_bm_opt.pkl'
pkl_dir = '/Users/rkrauss/git/research/SFLR_2010'
pkl_path = os.path.join(pkl_dir, pkl_name)

## #the controllers Gth2 and GaX are in pkl_dir
## if pkl_dir not in sys.path:
##     sys.path.insert(1, pkl_dir)

#f = logspace(-1, 1.5, 500)
f = logspace(-1, 2, 500)

s = 2.0j*pi*f


#?SFLR_TMM.SFLR_params?
mydict = LoadPickle(pkl_path)
params = SFLR_TMM.SFLR_params(**mydict)

fig_dir = os.path.join('..','figs')

def mysave(filename, fi):
    path = os.path.join(fig_dir, filename)
    PU.mysave(path, fi)



# Create TMM models:
def build_OL_TMM_model(pklpathin=None):
    if pklpathin is None:
        pklpathin = pkl_path
    OL_TMM_model = SFLR_TMM.model_w_bm(pklpathin)
    return OL_TMM_model


def build_OL_Bode_plotter():
    relpath = 'siue/Research/PSoC_Research/SFLR_2010/data/swept_sine/August_2010/after_beam_reattachment'
    data_folder = rwkos.FindFullPath(relpath)
    if data_folder not in sys.path:
        sys.path.insert(1, data_folder)

    exp_mod_name = 'swept_sine_kp_1_after_beam_reattachment_amp_75_maxf_20_duration_40000_avebodes'

    OL_Exp_Bode_Plotter = BPO.exp_bode_object(exp_mod_name, \
                                              OL_Bode_opts)
    return OL_Exp_Bode_Plotter


def OL_bode_plots(fi=1, resave=1, TMM=True):
    OL_Exp_Bode_Plotter = build_OL_Bode_plotter()
    relpath2 = 'siue/Research/PSoC_Research/SFLR_2010/data/swept_sine/July_2010/JVC/July_02_2010'

    data_folder2 = rwkos.FindFullPath(relpath2)
    if data_folder2 not in sys.path:
        sys.path.insert(1, data_folder2)


    Kp_mod_name = 'swept_sine_kp_1_good_amp_75_maxf_20_duration_40000_avebodes'
    log_ds_mod_name = 'swept_sine_kp_1_after_beam_reattachment_amp_75_maxf_20_duration_40000_logdownsampled'

    #Kp_mod_name
    OL_Exp_Bode_Plotter2 = BPO.exp_bode_object(log_ds_mod_name, \
                                               OL_Bode_opts)

    OL_Plotters = [OL_Exp_Bode_Plotter]

    if TMM:
        OL_TMM_model = build_OL_TMM_model()
        OL_TMM_Bode_Plotter = BPO.TMM_bode_object(OL_TMM_model, \
                                                  OL_Bode_opts)
        OL_Plotters.append(OL_TMM_Bode_Plotter)

    OL_overlayer = BPO.Bode_Overlayer(OL_Plotters, \
                                      OL_Bode_opts)
    OL_overlayer.plot_bodes(f, startfi=fi)
    #OL_Exp_Bode_Plotter2.plot_bodes(f, startfi=fi, trunc=False, clear=False, linetype='o')

    PU.SetLegend(fi,axis=0,loc=3)
    PU.SetLegend(fi+1,axis=0,loc=2)


    if resave:
        mysave('ol_theta_v_bode.eps', fi)
        mysave('ol_a_v_bode.eps', fi+1)


def build_Gth_TMM_model(pklpathin=None):
    if pklpathin is None:
        pklpathin = pkl_path

    import Gth
    reload(Gth)

    Gth_TMM = SFLR_TMM.model_w_bm_with_theta_feedback_comp(pklpathin, \
                                                           p=Gth.p, \
                                                           z=Gth.z, \
                                                           gain=Gth.gain)
    return Gth_TMM



def build_ThetaFB_exp_Bode_plotter():
    thetafb_mod_name = 'notch_filtered_swept_sine_Gth_comp_no_accel_FB_amp_15_duration_40000_maxf_20_avebodes'
    relpath = 'siue/Research/PSoC_Research/SFLR_2010/data/swept_sine/August_2010/after_beam_reattachment/TMM'
    data_folder = rwkos.FindFullPath(relpath)
    if data_folder not in sys.path:
        sys.path.insert(1, data_folder)

    ThetaFB_Exp_Bode_Plotter = BPO.exp_bode_object(thetafb_mod_name, \
                                                   ThetaFB_Bode_opts)
    return ThetaFB_Exp_Bode_Plotter


def ThetaFB_bode_plots(fi=3, resave=1):
    ThetaFB_Exp_Bode_Plotter = build_ThetaFB_exp_Bode_plotter()
    design_dir = rwkos.FindFullPath('SFLR_2010/vibration_suppression_design/TMM')

    if design_dir not in sys.path:
        sys.path.append(design_dir)


    Gth_TMM = build_Gth_TMM_model()

    import Gth
    reload(Gth)

    ThetaFB_TMM_Bode_Plotter = BPO.TMM_bode_object(Gth_TMM, \
                                                   ThetaFB_Bode_opts)

    ThetaFB_Plotters =[ThetaFB_Exp_Bode_Plotter, \
                       ThetaFB_TMM_Bode_Plotter, \
                      ]

    ThetaFB_overlayer = BPO.Bode_Overlayer(ThetaFB_Plotters, \
                                           ThetaFB_Bode_opts)
    ThetaFB_overlayer.plot_bodes(f, startfi=fi)

    th_comp, a_comp = model_w_bm_bodes_ThetaFB_editted.Bodes(s, \
                                                             params, \
                                                             Gth.Gth)
    PU.SetLegend(fi,axis=0,loc=3)
    PU.SetLegend(fi+1,axis=0,loc=2)

    if resave:
        mysave('thetafb_theta_thd_bode.eps', fi)
        mysave('thetafb_a_thd_bode.eps', fi+1)


    th_bode_opt = ThetaFB_Bode_opts[0]
    a_bode_opt = ThetaFB_Bode_opts[1]
    th_bode = BPO.comp_to_Bode(th_comp, f, bode_opt=th_bode_opt)
    a_bode = BPO.comp_to_Bode(a_comp, f, bode_opt=a_bode_opt)
    BPO._plot_bode(th_bode, th_bode_opt, f, fignum=fi, linetype='k:')
    BPO._plot_bode(a_bode, a_bode_opt, f, fignum=fi+1, linetype='k:')





def build_Accel_TMM_model(pklpathin=None, Ga=None):
    if pklpathin is None:
        pklpathin = pkl_path

    import Gth
    reload(Gth)

    import Ga5
    reload(Ga5)

    import Unc_row
    reload(Unc_row)

    if Ga is None:
        Ga = Ga5.Ga5

    Gth_Ga_TMM = SFLR_TMM.model_w_bm_accel_and_theta_FB(pkl_path, \
                                                        Gth=Gth.Gth, \
                                                        Ga=Ga, \
                                                        Unc_func=Unc_row.Unc_row0)
    return Gth_Ga_TMM


def build_AccelFB_Bode_plotter():
    accelfb_mod_name = 'notch_filtered_swept_sine_Gth_comp_Ga5_amp_15_duration_40000_maxf_20_avebodes'
    relpath = 'siue/Research/PSoC_Research/SFLR_2010/data/swept_sine/August_2010/after_beam_reattachment/TMM'
    data_folder = rwkos.FindFullPath(relpath)
    if data_folder not in sys.path:
        sys.path.insert(1, data_folder)

    #I think the bode_opts for accel and theta FB are the same
    AccelFB_Exp_Bode_Plotter = BPO.exp_bode_object(accelfb_mod_name, \
                                                   AccelFB_Bode_opts)

    return AccelFB_Exp_Bode_Plotter


def build_AccelFB_TMM_plotter():
    Gth_Ga_TMM = build_Accel_TMM_model()

    AccelFB_TMM_Bode_Plotter = BPO.TMM_bode_object(Gth_Ga_TMM, \
                                                   AccelFB_Bode_opts)
    return AccelFB_TMM_Bode_Plotter


def AccelFB_bode_plots(fi=5, resave=1):
    AccelFB_Exp_Bode_Plotter = build_AccelFB_Bode_plotter()
    design_dir = rwkos.FindFullPath('SFLR_2010/vibration_suppression_design/TMM')

    if design_dir not in sys.path:
        sys.path.append(design_dir)


    import Gth
    reload(Gth)

    import Ga5
    reload(Ga5)

    import Unc_row
    reload(Unc_row)

    def myfunc(s):
        return 0.0

    Gth_Ga_TMM = build_Accel_TMM_model()

    AccelFB_TMM_Bode_Plotter = BPO.TMM_bode_object(Gth_Ga_TMM, \
                                                   AccelFB_Bode_opts)

    AccelFB_Plotters =[AccelFB_Exp_Bode_Plotter, \
                       AccelFB_TMM_Bode_Plotter, \
                      ]

    AccelFB_overlayer = BPO.Bode_Overlayer(AccelFB_Plotters, \
                                           AccelFB_Bode_opts)
    AccelFB_overlayer.plot_bodes(f, startfi=fi, PhaseMassage=True)

    PU.SetLegend(fi,axis=0,loc=3)
    PU.SetLegend(fi+1,axis=0,loc=2)
    if resave:
        mysave('accelfb_theta_thd_bode.eps', fi)
        mysave('accelfb_a_thd_bode.eps', fi+1)



if __name__ == '__main__':
    OL_bode_plots()
    ThetaFB_bode_plots()
    AccelFB_bode_plots()

    show()
