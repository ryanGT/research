from pylab import *
from scipy import *
from scipy import optimize

import pylab_util as PU
reload(PU)
import os, sys

import rwkos, rwkbode, rwkmisc

import txt_data_processing as TDP
#reload(TDP)

import rwkbode, pylab_util, rwkos
reload(rwkbode)

import controls
TF = controls.TransferFunction

import copy

import SFLR_TMM
reload(SFLR_TMM)

import bode_plot_overlayer as BPO
reload(BPO)

from IPython.Debugger import Pdb

#import bode_options

def notch_tf(C):
    wz = C[0]
    zz = C[1]
    wp = C[2]
    zp = C[3]
    notch = TF([1,2*wz*zz,wz**2],\
               [1,2*wp*zp,wp**2])*wp**2/wz**2
    return notch


class SFLR_Time_File_Mixin(object):
    def load_exp_time_file(self, filepath, \
                           col_map={0:'t', 1:'n', 2:'u', \
                                    3:'v', 4:'theta', 5:'a'}):
        self.data_file = TDP.Data_File(path=filepath, \
                                       col_map=col_map)
        self.t = self.data_file.t


    def create_ax(self, fi=1, clear=True):
        fig = figure(fi)
        if clear:
            fig.clf()
        self.ax = fig.add_subplot(1,1,1)
        return self.ax


    def plot_exp_time_data(self, accel=False):
        ax = self.ax
        t = self.t
        u = self.data_file.u
        ax.plot(t, u, label='$u$')
        ax.plot(t, self.data_file.v, label='$v_{exp}$')
        ax.plot(t, self.data_file.theta, label='$\\theta_{exp}$')
        if accel:
            ax.plot(t, self.data_file.a, label='$\\ddot{x}_{exp}$')


    def plot_model_data(self, accel=False):
        ax = self.ax
        t = self.t
        ax.plot(t, self.theta, label='$\\theta_{model}$')
        if accel:
            ax.plot(t, self.accel, label='$\\ddot{x}_{model}$')


    def label_time_domain_plot(self):
        ax = self.ax
        ax.set_xlabel('Time (sec)')
        ax.set_ylabel('Signal Amplitude (counts)')        


    def plot_time_domain_exp_vs_model(self, fi=1, legloc=5):
        self.create_ax(fi=fi)
        self.plot_exp_time_data()
        self.plot_model_data()
        self.label_time_domain_plot()
        self.ax.legend(loc=legloc)


    def lsim_from_exp_file(self, filepath, fi=1, plot=True, \
                           clear=True):
        self.load_exp_time_file(filepath)
        u = self.data_file.u
        t = self.data_file.t
        #Pdb().set_trace()
        self.lsim(u, t)
        if plot:
            self.create_ax(fi=fi, clear=clear)
            self.plot_exp_time_data()
            self.plot_model_data()
            self.label_time_domain_plot()


class SFLR_Time_File_Mixin_w_accel(SFLR_Time_File_Mixin):
    def plot_exp_time_data(self):
        SFLR_Time_File_Mixin.plot_exp_time_data(self, accel=True)


    def plot_model_data(self):
        SFLR_Time_File_Mixin.plot_model_data(self, accel=True)

    
class Rigid_Acuator_TF_Model(SFLR_Time_File_Mixin):
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

    def plot_one_bode(self, f, index, fignum=1, \
                      clear=False, **kwargs):
        """index refers to which bode option to pull from
        self.bode_opts.  This is how the method determines which bode
        to plot."""
        if not hasattr(self, 'bodes'):
            self.calc_bodes(f)
        if not kwargs.has_key('label'):
            kwargs['label'] = self.label
        opt = self.bode_opts[index]
        bode = self.find_bode(opt.output_label, opt.input_label)
        BPO._plot_bode(bode, opt, f, fignum=fignum,
                       clear=clear, **kwargs)


    def build_act_iso(self):
        p_act = self.p_act1
        K_act = self.K_act
        H = self.H
        s1 = 1.0*2.0j*pi#magnitude of s at 1 Hz - fix this point for
        m1 = abs(s1*(s1+p_act))
        num = K_act*m1#*m2
        self.G_act_iso = TF(num*H,[1,p_act,0])
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
        cost = self.G_act_bode.dB_Phase_Error_sum(self.TMM_act_bode)
        cost += self.negative_params_check()
        return cost


    def calc_TMM_bode(self):
        self.TMM_model.calc_bodes(self.ffit)
        self.TMM_act_bode = self.TMM_model.find_bode(self.bode_opts[0])


    def __init__(self, bode_opts, unknown_params, \
                 TMM_model, \
                 ffit=None, \
                 label='rigid act.', \
                 params={'K_act':(0.06727320631008038/(2*pi)), \
                         'H':180.0/pi*1024/360.0, \
                         'p_act1':120.40659975047289}):
        if ffit is None:
            ffit = logspace(-1,1.5,200)
        self.TMM_model = TMM_model
        self.ffit = ffit
        self.bode_opts = bode_opts
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


    def set_params_from_dict(self, dictin):
        """Set the parameters from a dictionary of key, value pairs."""
        for attr, val in dictin.iteritems():
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


    def lsim(self, u, t):
        self.theta = self.G_act.lsim(u, t)
        #self.accel = self.G_a_th.lsim(self.theta, t)
        return self.theta#, self.accel



    def fit_time_domain(self, filepath):
        self.load_exp_time_file(filepath)#sets self.data_file
        ig = self.get_ig()
        X_opt = optimize.fmin(self.time_domain_cost, ig)
        self.X_opt = X_opt
        self.set_params(X_opt)
        self.build_TFs()
        self.lsim(self.data_file.u, self.data_file.t)
        self.plot_time_domain_exp_vs_model()
        return X_opt


    def time_domain_cost(self, C):
        self.set_params(C)
        self.build_TFs()
        self.lsim(self.data_file.u, self.data_file.t)
        e_theta = self.theta - self.data_file.theta
        #e_accel = self.accel - self.data_file.a
        cost = sum(e_theta**2)# + sum(e_accel**2)
        #cost += self.negative_params_check()
        return cost


    def save_fit_results(self, filepath):
        self.fit_res_dict = dict(zip(self.unknown_params, self.X_opt))
        rwkmisc.SavePickle(self.fit_res_dict, filepath)


    def load_fit_results(self, filepath):
        self.fit_res_dict = rwkmisc.LoadPickle(filepath)
        self.set_params_from_dict(self.fit_res_dict)
        self.build_TFs()


    def update_params(self):
        """Update the values in self.params from getattr(self, key).
        This is necessary when you have loaded fit results and then
        want to save all the params."""
        for key in self.params.iterkeys():
            val = getattr(self, key)
            self.params[key] = val


    def fix_pole_zero_cancelation(self):
        """If someone accidentally allowed one of the poles or zeros
        that should cancel between Gth and Ga vary and not the other,
        this should fix it."""
        self.wpa1 = self.wz1
        self.wpa2 = self.wz2
        self.zpa1 = self.zz1
        self.zpa2 = self.zz2



    def save_params(self, filepath):
        rwkmisc.SavePickle(self.params, filepath)
        return self.params

    
    def load_params(self, filepath):
        mydict = rwkmisc.LoadPickle(filepath)
        self.set_params_from_dict(mydict)
        self.build_TFs()
        return mydict
        

class Second_Order_Rigid_Act_Model(Rigid_Acuator_TF_Model):
    def build_act_iso(self):
        p_act1 = self.p_act1
        p_act2 = self.p_act2
        K_act = self.K_act
        H = self.H
        s1 = 1.0*2.0j*pi#magnitude of s at 1 Hz - fix this point for
        m1 = abs(s1*(s1+p_act1)*(s1+p_act2))
        num = K_act*m1#*m2
        pole2 = TF(1,[1,p_act2])
        self.G_act_iso = TF(num*H,[1,p_act1,0])*pole2
        return self.G_act_iso


    def __init__(self, bode_opts, unknown_params, \
                 TMM_model, \
                 ffit=None, \
                 label='SO rigid act.', \
                 params = {'K_act':9.78920392e-03, \
                           'H':180.0/pi*1024/360.0, \
                           'p_act1':3.65273247e+01, \
                           'p_act2':100.0}):
        Rigid_Acuator_TF_Model.__init__(self, bode_opts, \
                                        unknown_params, \
                                        TMM_model, ffit=ffit, \
                                        label=label, \
                                        params=params)


class model_w_notch(object):
    def _build_notch_TF(self, wp, zp, wz, zz):
        """Create a second order notch filter.  This is a generic
        notch filter method called by build_notch1_TF and
        build_notch2_TF (if the class has a build_notch2_TF method)."""
        notch = TF([1,2*wz*zz,wz**2],\
                   [1,2*wp*zp,wp**2])*wp**2/wz**2
        return notch


    def build_notch1_TF(self):
        wp1 = self.wp1
        zp1 = self.zp1
        wz1 = self.wz1
        zz1 = self.zz1
        self.notch1 = self._build_notch_TF(wp1, zp1, wz1, zz1)

    
class FO_Act_w_one_notch(model_w_notch,Rigid_Acuator_TF_Model):
    def build_TFs(self):
        self.build_act_iso()
        self.build_notch1_TF()
        self.G_act = self.G_act_iso*self.notch1


    def __init__(self, bode_opts, unknown_params, \
                 TMM_model, \
                 ffit=None, \
                 label='FO 1 notch', \
                 params = {'K_act':1.12347211e-02, \
                           'H':180.0/pi*1024/360.0, \
                           'p_act1':1.38053790e+02, \
                           'wp1':1.77094515e+01, \
                           'zp1':4.14416466e-01, \
                           'wz1':1.55761111e+01, \
                           'zz1':5.10606176e-03}, \
                 ):
        Rigid_Acuator_TF_Model.__init__(self, bode_opts, \
                                        unknown_params, \
                                        TMM_model, ffit=ffit, \
                                        label=label, \
                                        params=params)

    

class SO_Act_w_one_notch(Second_Order_Rigid_Act_Model,model_w_notch):
    def build_TFs(self):
        self.build_act_iso()
        self.build_notch1_TF()
        self.G_act = self.G_act_iso*self.notch1


    def __init__(self, bode_opts, unknown_params, \
                 TMM_model, \
                 ffit=None, \
                 label='SO 1 notch', \
                 params = {'K_act':1.12347211e-02, \
                           'H':180.0/pi*1024/360.0, \
                           'p_act1':1.38053790e+02, \
                           'p_act2':3.84793971e+01, \
                           'wp1':1.77094515e+01, \
                           'zp1':4.14416466e-01, \
                           'wz1':1.55761111e+01, \
                           'zz1':5.10606176e-03}, \
                 ):
        Rigid_Acuator_TF_Model.__init__(self, bode_opts, \
                                        unknown_params, \
                                        TMM_model, ffit=ffit, \
                                        label=label, \
                                        params=params)



SO2N_keys = ['K_act','p_act1','p_act2','wp1','zp1','wz1',\
             'zz1','wp2','zp2','wz2','zz2']
SO2N_values = array([1.12144141e-02, 1.14932348e+02, 4.48157368e+01, \
                     1.74381691e+01, 4.16390661e-01, 1.55749618e+01, \
                     5.08152097e-03, 1.12829058e+02, 2.59982690e-02, \
                     1.11433819e+02, 1.63493872e-02])
SO2N_params = dict(zip(SO2N_keys, SO2N_values))
SO2N_params['H'] = 180.0/pi*1024/360.0

class SO_Act_w_two_notches(SO_Act_w_one_notch):
    def build_notch2_TF(self):
        wp2 = self.wp2
        zp2 = self.zp2
        wz2 = self.wz2
        zz2 = self.zz2
        self.notch2 = self._build_notch_TF(wp2, zp2, wz2, zz2)
        
        
    def build_TFs(self):
        self.build_act_iso()
        self.build_notch1_TF()
        self.build_notch2_TF()
        self.G_act = self.G_act_iso*self.notch1*self.notch2


    def __init__(self, bode_opts, unknown_params, \
                 TMM_model, \
                 ffit=None, \
                 label='SO 2 notch', \
                 params=SO2N_params, \
                 ):
        Rigid_Acuator_TF_Model.__init__(self, bode_opts, \
                                        unknown_params, \
                                        TMM_model, ffit=ffit, \
                                        label=label, \
                                        params=params)

         
accel_keys = ['wpa1', 'zpa1', 'B1', 'wpa2', 'zpa2', 'B2']
accel_values = [2.8*2*pi, 0.3, 1.02595313, 18.0*2*pi, 0.05, -4.93599049]
accel_dict = dict(zip(accel_keys, accel_values))
#By default, the poles of the a/theta TF will be exactly the zeros of
#the theta/v TF
accel_dict['wpa1'] = SO2N_params['wz1']
accel_dict['wpa2'] = SO2N_params['wz2']
accel_dict['zpa1'] = SO2N_params['zz1']
accel_dict['zpa2'] = SO2N_params['zz2']
accel_dict['a_gain'] = 1.55990124

accel_params = copy.copy(SO2N_params)
accel_params.update(accel_dict)

FO_accel_params = copy.copy(accel_params)
pop_params = ['p_act2','wp2','zp2','wz2','zz2','B2']
for item in pop_params:
    FO_accel_params.pop(item)
    

class TF_w_accel(object):
    def lsim(self, u, t):
        self.theta = self.G_act.lsim(u, t)
        self.accel = self.G_a_th.lsim(self.theta, t)
        return self.theta, self.accel


    def plot_exp_time_data(self):
        Rigid_Acuator_TF_Model.plot_exp_time_data(self, accel=True)


    def plot_model_data(self):
        Rigid_Acuator_TF_Model.plot_model_data(self, accel=True)


    def time_domain_cost(self, C):
        self.set_params(C)
        self.build_TFs()
        self.lsim(self.data_file.u, self.data_file.t)
        e_theta = self.theta - self.data_file.theta
        e_accel = self.accel - self.data_file.a
        cost = sum(e_theta**2) + sum(e_accel**2)
        #cost += self.negative_params_check()
        return cost


    def calc_bodes(self, f):
        bode1 = BPO.tf_to_Bode(self.G_act, f, \
                               self.bode_opts[0], PhaseMassage=True)
        bode2 = BPO.tf_to_Bode(self.G_a_v, f, \
                               self.bode_opts[1], PhaseMassage=True)

        self.bodes = [bode1, bode2]
        self.G_act_bode = bode1
        self.G_a_v_bode = bode2


    def calc_TMM_bode(self):
        self.TMM_model.calc_bodes(self.ffit)
        self.TMM_act_bode = self.TMM_model.find_bode(self.bode_opts[0])
        self.TMM_accel_v_bode = self.TMM_model.find_bode(self.bode_opts[1])


    def accel_v_cost(self, C):
        self.set_params(C)
        self.build_TFs()
        self.calc_bodes(self.ffit)
        cost = self.G_a_v_bode.dB_Phase_Error_sum(self.TMM_accel_v_bode)
        #cost += self.negative_params_check()
        return cost


class FO_1_notch_w_accel(TF_w_accel,FO_Act_w_one_notch):
    def build_a_theta_TF(self):
        wpa1 = self.wpa1
        zpa1 = self.zpa1
        B1 = self.B1
        mode1 = TF([B1,0,0],\
                   [1,2*zpa1*wpa1,wpa1**2])
        self.G_a_th = (mode1)*self.a_gain


    def build_TFs(self):
        self.build_act_iso()
        self.build_notch1_TF()
        self.build_a_theta_TF()
        self.G_act = self.G_act_iso*self.notch1
        self.G_a_v = self.G_act*self.G_a_th


    def __init__(self, bode_opts, unknown_params, \
                 TMM_model, \
                 ffit=None, \
                 label='FO 1 notch', \
                 params=FO_accel_params, \
                 ):
        Rigid_Acuator_TF_Model.__init__(self, bode_opts, \
                                        unknown_params, \
                                        TMM_model, ffit=ffit, \
                                        label=label, \
                                        params=params)
        self.my_cost = self.accel_v_cost
    
    
class Accel_w_two_notches(TF_w_accel,SO_Act_w_two_notches):
    """This model builds on SO_Act_w_two_notches by adding the
    accelerometer output."""
    def build_a_theta_TF(self):
        wpa1 = self.wpa1
        zpa1 = self.zpa1
        B1 = self.B1
        wpa2 = self.wpa2
        zpa2 = self.zpa2
        B2 = self.B2
        mode1 = TF([B1,0,0],\
                   [1,2*zpa1*wpa1,wpa1**2])
        mode2 = TF([B2,0,0],\
                   [1,2*zpa2*wpa2,wpa2**2])
        self.G_a_th = (mode1 + mode2)*self.a_gain

        
    def build_TFs(self):
        self.build_act_iso()
        self.build_notch1_TF()
        self.build_notch2_TF()
        self.build_a_theta_TF()
        self.G_act = self.G_act_iso*self.notch1*self.notch2
        self.G_a_v = self.G_act*self.G_a_th


    def __init__(self, bode_opts, unknown_params, \
                 TMM_model, \
                 ffit=None, \
                 label='ROM', \
                 params=accel_params, \
                 ):
        Rigid_Acuator_TF_Model.__init__(self, bode_opts, \
                                        unknown_params, \
                                        TMM_model, ffit=ffit, \
                                        label=label, \
                                        params=params)
        self.my_cost = self.accel_v_cost


class P_control_Theta_FB(Accel_w_two_notches):
    def build_TFs(self):
        self.build_act_iso()
        self.build_notch1_TF()
        self.build_notch2_TF()
        self.build_a_theta_TF()
        self.G_act_ol = self.G_act_iso*self.notch1*self.notch2
        self.G_act = controls.feedback(self.G_act_ol*self.kp)
        self.G_a_v = self.G_act*self.G_a_th


    def __init__(self, bode_opts, unknown_params, \
                 TMM_model, \
                 kp=1.0, \
                 ffit=None, \
                 label='ROM', \
                 params=accel_params, \
                 ):
        self.kp = kp
        Accel_w_two_notches.__init__(self, bode_opts, \
                                     unknown_params, \
                                     TMM_model, ffit=ffit, \
                                     label=label, \
                                     params=params)


p=20*2*pi
z=3.0*2*pi
gain=2.0
Gth_comp = TF([1,z], [1,p])*gain*(p/z)

class G_th_comp_Theta_FB(Accel_w_two_notches):
    def _build_TFs(self):
        self.build_act_iso()
        self.build_notch1_TF()
        self.build_notch2_TF()
        self.build_a_theta_TF()
        self.G_act_ol = self.G_act_iso*self.notch1*self.notch2
        
    def build_TFs(self):
        self._build_TFs()
        self.G_act = controls.feedback(self.G_act_ol*self.Gth)
        self.G_a_v = self.G_act*self.G_a_th


    def __init__(self, bode_opts, unknown_params, \
                 TMM_model, \
                 Gth=Gth_comp, \
                 ffit=None, \
                 label='ROM', \
                 params=accel_params, \
                 ):
        self.Gth = Gth
        Accel_w_two_notches.__init__(self, bode_opts, \
                                     unknown_params, \
                                     TMM_model, ffit=ffit, \
                                     label=label, \
                                     params=params)


class G_th_G_a_TF(G_th_comp_Theta_FB):
    def simplify(self):
        self.atf2 = self.accel_TF.simplify()
        self.thtf2 = self.theta_TF.simplify()
        
    def lsim(self, u, t, simplify=True):
        if simplify:
            self.simplify()
            self.theta = self.thtf2.lsim(u, t)
            self.accel = self.atf2.lsim(u, t)
        else:
            self.theta = self.theta_TF.lsim(u, t)
            self.accel = self.accel_TF.lsim(u, t)
        return self.theta, self.accel

    
    def calc_bodes(self, f):
        #Will need to fix this to find the bodes for a/theta_d and
        #theta/theta_d
        bode1 = BPO.tf_to_Bode(self.theta_TF, f, \
                               self.bode_opts[0], PhaseMassage=True)
        bode2 = BPO.tf_to_Bode(self.accel_TF, f, \
                               self.bode_opts[1], PhaseMassage=True)

        self.bodes = [bode1, bode2]
        self.G_act_bode = bode1#<-- this is a bad variable name
        self.G_a_v_bode = bode2#<-- this is a bad variable name


    def build_TFs(self):
        self._build_TFs()
        self.G_act_ol = self.G_act_iso*self.notch1*self.notch2
        self.G_act = controls.feedback(self.G_act_ol*self.Gth)
        self.G_a_theta_d_hat = self.G_act*self.G_a_th
        self.accel_TF = controls.feedback(self.G_a_theta_d_hat, self.Ga)
        self.theta_TF = self.accel_TF/self.G_a_th
        #Figure out how to find the a/theta_d bode here and back out
        #theta/theta_d with Ga


    def __init__(self, bode_opts, unknown_params, \
                 TMM_model, \
                 Gth=Gth_comp, \
                 Ga=None, \
                 ffit=None, \
                 label='ROM', \
                 params=accel_params, \
                 ):
        self.Ga = Ga
        G_th_comp_Theta_FB.__init__(self, bode_opts, \
                                    unknown_params, \
                                    TMM_model, \
                                    Gth=Gth, \
                                    ffit=ffit, \
                                    label=label, \
                                    params=params)
        


    
class G_th_comp_Theta_FB_no_accel(G_th_comp_Theta_FB):
    def lsim(self, u, t):
        self.theta = self.G_act.lsim(u, t)
        #self.accel = self.G_a_th.lsim(self.theta, t)
        return self.theta#, self.accel
    
    def time_domain_cost(self, C):
        self.set_params(C)
        self.build_TFs()
        self.lsim(self.data_file.u, self.data_file.t)
        e_theta = self.theta - self.data_file.theta
        #e_accel = self.accel - self.data_file.a
        cost = sum(e_theta**2)# + sum(e_accel**2)
        #cost += self.negative_params_check()
        return cost

    def plot_exp_time_data(self):
        SO_Act_w_two_notches.plot_exp_time_data(self, accel=False)


    def plot_model_data(self):
        SO_Act_w_two_notches.plot_model_data(self, accel=False)

    
def sat(vin, mymax=200):
    """The PSoC can only ouput voltages in the range +/- 2.5V.  With a
    9-bit digital to analog converter, this corresponds to +/- 255
    counts.  It is slightly nonlinear above +/- 200."""
    vmax = mymax
    vmin = -1*mymax
    if vin > vmax:
        return vmax
    elif vin < vmin:
        return vmin
    else:
        return vin

    
class G_th_comp_Theta_FB_dig_w_sat(G_th_comp_Theta_FB):
    def __init__(self, bode_opts, unknown_params, \
                 TMM_model, \
                 Gth=Gth_comp, \
                 ffit=None, \
                 label='ROM', \
                 params=accel_params, \
                 dt=1.0/500, \
                 ):
        G_th_comp_Theta_FB.__init__(self, bode_opts, \
                                    unknown_params, \
                                    TMM_model, \
                                    Gth=Gth, \
                                    ffit=ffit, \
                                    label=label, \
                                    params=params)
        self.dt = dt
        self.Gth_num, self.Gth_den = self.Gth.c2d_tustin(dt=self.dt)
        self.Gth_z = controls.Digital_Compensator(self.Gth_num, self.Gth_den)


    def plot_model_vvect(self):
        ax = self.ax
        ax.plot(self.t, self.vvect, label='$v_{model}')


    def Run_Sim(self, u, t):
        """Simulate the response of open-loop transfer function G to input
        u.  vfunc specifes how the input voltage is calculated.

        You could just find the closed-loop transfer function and simulate
        its response to u, as long as the voltage doesn't saturate."""
        N = len(u)
        self.evect = zeros((N,), dtype='int32')
        self.vvect = zeros((N,), dtype='float64')
        self.v_int = zeros((N,), dtype='int32')
        self.theta = zeros((N,), dtype='float64')
        self.accel = zeros((N,), dtype='float64')

        self.Gth_z.input = self.evect
        self.Gth_z.output = self.vvect

        dt = self.dt

        prevx = None

        for i in range(1,N):
            self.evect[i] = u[i] - self.theta[i-1]
            vtemp = self.Gth_z.calc_out(i)
            self.vvect[i] = sat(vtemp)
            self.v_int[i] = int(self.vvect[i])

            t1 = (i-1)*dt
            t2 = i*dt
            t_i, y_i, x_i = self.G_act_ol.lsim([self.v_int[i],self.v_int[i]], \
                                               [t1,t2], returnall=True, X0=prevx)
            prevx = x_i[-1]
            self.theta[i] = y_i[-1]

        self.accel = self.G_a_th.lsim(self.theta, t)


        return self.theta, self.accel, self.vvect
    

    def Run_Sim_from_exp_file(self, filepath, fi=1, plot=True, \
                              clear=True):
        self.load_exp_time_file(filepath)
        u = self.data_file.u
        t = self.data_file.t
        self.Run_Sim(u, t)
        if plot:
            self.create_ax(fi=fi, clear=clear)
            self.plot_exp_time_data()
            self.plot_model_data()
            self.plot_model_vvect()
            self.label_time_domain_plot()
