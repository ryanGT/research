from pylab import *
from scipy import *
import copy, os
import txt_data_processing, rwkbode, pylab_util, rwkos

import pylab_util as PU

#from systemid import Model, PolyHasher

data_dir = rwkos.FindFullPath('siue/Research/modeling/SLFR/data/July_07_2009')
if data_dir not in sys.path:
    sys.path.append(data_dir)

def load_data(mod_name):
    #data_mod_name = 'swept_sine_amp_75_July_07_2009_avebodes'
    #data_mod_name = 'swept_sine_amp_75_July_07_2009_log_downsampled'
    bode_data_set = txt_data_processing.load_avebode_data_set(mod_name)

    #bode_data_set.Bode_Plot2(func=rwkbode.GenBodePlot, linetype='o')
    exp_bodes = bode_data_set.avebodes[0:2]
    th_v_exp = exp_bodes[0]
    a_v_exp = exp_bodes[1]
    f = bode_data_set.f
    a_theta_exp = bode_data_set.avebodes[2]
    return f, th_v_exp, a_v_exp, a_theta_exp

## l1 = log10(0.5)
## l2 = log10(50)
## flist = arange(0.5, 1, 0.5).tolist()+f.tolist()+arange(26,30, 0.5).tolist()
## f2 = array(flist)
