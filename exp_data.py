from pylab import *
from scipy import *
import copy, os
import txt_data_processing, rwkbode, pylab_util, rwkos

import pylab_util as PU

#from systemid import Model, PolyHasher

case = 2

if case == 1:
    data_dir = rwkos.FindFullPath('siue/Research/modeling/SLFR/data/July_07_2009')
    data_mod_name = 'swept_sine_amp_75_July_07_2009_avebodes'

elif case == 2:
    relpath = 'siue/Research/SFLR_2010/data/swept_sine/August_2010/after_beam_reattachment'
    data_dir = rwkos.FindFullPath(relpath)
    data_mod_name = 'swept_sine_kp_1_after_beam_reattachment_amp_75_maxf_20_duration_40000_avebodes'

data_mod_path = os.path.join(data_dir, data_mod_name)
    
if data_dir not in sys.path:
    sys.path.append(data_dir)


#data_mod_name = 'swept_sine_amp_75_July_07_2009_log_downsampled'
#bode_data_set = txt_data_processing.load_avebode_data_set(data_mod_name)
bode_data_set = txt_data_processing.load_avebode_data_set(data_mod_path)

#bode_data_set.Bode_Plot2(func=rwkbode.GenBodePlot, linetype='o')
exp_bodes = bode_data_set.avebodes#[0:2]
th_v_exp = bode_data_set.find_bode('theta','v')
a_v_exp = bode_data_set.find_bode('a','v')
f = bode_data_set.f
a_theta_exp = bode_data_set.find_bode('a','theta')

l1 = log10(0.5)
l2 = log10(50)
flist = arange(0.5, 1, 0.5).tolist()+f.tolist()+arange(26,30, 0.5).tolist()
f2 = array(flist)

def plot_exp(startfi=1, lt='o'):
    kwargs = {'clear':True, 'linetype':lt}
    rwkbode.GenBodePlot(startfi, f, th_v_exp, **kwargs)
    rwkbode.GenBodePlot(startfi+1, f, a_v_exp, **kwargs)
    rwkbode.GenBodePlot(startfi+2, f, a_theta_exp, **kwargs)
    return startfi+3


freqlim = [0.5, 30.0]
maglims = [(-40,15),(-20,10),(-20,45)]
phaselims = [(-220,0),(-400,120),(-200,220)]
NF = 3
figdir = 'figs'

def _set_lims(startfi=1, maglims=maglims, phaselims=phaselims):
    fis = range(startfi, startfi+NF)
    for fi, maglim, phaselim in zip(fis, maglims, phaselims):
        PU.SetFreqLim(fi, freqlim)
        PU.SetMagLim(fi, maglim)
        PU.SetPhaseLim(fi, phaselim)
    

def _save_figs(filenames, startfi=1, del_eps=True):
    fis = range(startfi, startfi+NF)
    for fi, fname in zip(fis, filenames):
        epsname = fname+'.eps'
        outpath = os.path.join(figdir, epsname)
        PU.mysave(outpath, fi=fi)
        if del_eps:
            os.remove(outpath)#delete eps but leave pdf
        

def set_legends(startfi, leg_list, locs):
    fis = range(startfi, startfi+NF)
    for fi, loc in zip(fis, locs):
        PU.SetLegend(fi, leg_list, loc=loc)

    
if __name__ == '__main__':
    plot_exp(lt='-')
    olfilenames = ['th_v_ol_bode','a_v_ol_bode','a_th_ol_bode']
    _set_lims(1)
    _save_figs(olfilenames, 1)

    show()
    

