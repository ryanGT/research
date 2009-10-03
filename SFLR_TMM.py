from __future__ import division
import TMM
reload(TMM)
from scipy import pi, sqrt, arange, vectorize, exp, c_, array, transpose, real, imag, rand, cos, sin, sinh, cosh, argmax, arange, eye, zeros, shape, poly1d
import pylab
import scipy
import pdb
import time
#import SAMII
import rwkascii
from rwkmisc import prettymat, colwise, SavePickle, LoadPickle
from rwkos import FindFullPath
import shutil, os, sys, glob
import TMM.beam
reload(TMM.beam)
import TMM.rigid
reload(TMM.rigid)
import TMM.spring
reload(TMM.spring)
import TMM.velocitysource
reload(TMM.velocitysource)
import TMM.feedback
reload(TMM.feedback)
from TMM.beam import BeamElement, BeamElement_v2
from TMM.rigid import RigidMass
from TMM.spring import TorsionalSpringDamper
from TMM.velocitysource import AngularVelocitySource
from TMM.velocitysource import AVSwThetaFB
from TMM.feedback import SAMIIAccelFB
from rwkdataproc import datastruct
import re
import rwkparse
import rwkbode
reload(rwkbode)
from rwkbode import bodeout

import copy

ms=4

#beam dimensions in inches
Li = 17.0
wi = 13.0/16
ti = 1.0/32#did not use calipers
#ti = 0.029

i_to_m = 25.4/1000.0
L = Li*i_to_m
w = wi*i_to_m
t = ti*i_to_m
A = w*t

#calculating rho from total mass
Lt = (19.0+5.0/8)*25.4/1000#total length (including the part in the clamp)
m = 65.75/1000.0#total mass (65.75 grams)
mu = m/Lt#mass per unit length

E = 200*10**9#N/m**2 or 200 GPa
rho = 7850.0#kg/m**3
I = 1.0/12*w*t**3
#mu = A*rho

EI = E*I

#from Frank
#EI = 0.4167
#mu = 0.2
#L = 0.508

class SFLR_params(object):
    def __init__(self, \
                 K_act=0.05, \
                 p_act1=1.0*2*pi, \
                 p_act2=50.0*2*pi, \
                 z_act=3.0*2*pi, \
                 tau=10.0*2*pi, \
                 EI=0.172, \
                 L=0.4318, \
                 mu=0.1319, \
                 L2=None, \
                 c_beam=0.0, \
                 a_gain=1.5, \
                 k_spring=10e6, \
                 c_spring=0.0, \
                 k_clamp=10e6, \
                 c_clamp=0.0, \
                 a_m=0.0077, \
                 a_L=0.011125, \
                 a_r=0.005556, \
                 a_I=7.0627e-7, \
                 b_m=0.0077, \
                 b_L=0.011125, \
                 b_r=0.005556, \
                 b_I=7.0627e-7):
        self.K_act = K_act
        self.p_act1 = p_act1
        self.p_act2 = p_act2
        self.z_act = z_act
        self.tau = tau
        self.EI = EI
        self.L = L
        self.L2 = L2
        self.mu = mu
        self.c_beam = c_beam
        self.a_gain = a_gain
        self.k_spring = k_spring
        self.c_spring = c_spring
        self.k_clamp = k_clamp
        self.c_clamp = c_clamp
        self.a_m = a_m
        self.a_L = a_L
        self.a_r = a_r
        self.a_I = a_I
        self.b_m = b_m
        self.b_L = b_L
        self.b_r = b_r
        self.b_I = b_I


    def __repr__(self):
        outstr = ''
        for key, val in self.__dict__.iteritems():
            outstr += '%s = %s\n' % (key, val)
        return outstr


    def update(self, newdict):
        for key, val in newdict.iteritems():
            setattr(self, key, val)


    def __sub__(self, other):
        out = SFLR_params()
        for key, val in self.__dict__.iteritems():
            op = getattr(other, key)
            if (val is None) or (op is None):
                diff = None
            else:
                diff = val-op
            setattr(out, key, diff)
        return out

            
itm = 25.4/1000.0#convert inches to meters
d_a = (7.0/16.0)*itm#diameter of accel
r_a = d_a/2.0#radius of accel

Ln = 16.5*itm-r_a
L2 = 0.5*itm-r_a

class new_def_params(SFLR_params):
    def __init__(self, \
                 L=Ln, \
                 L2=L2, \
                 a_L=d_a, \
                 a_r=r_a, \
                 **kwargs):
        SFLR_params.__init__(self, L=L, L2=L2, a_L=a_L, a_r=a_r, \
                             **kwargs)


def save_params(params, filename):
    mydict = params.__dict__
    SavePickle(mydict, filename, protocol=1)

def load_params(filename):
    mydict = LoadPickle(filename)
    myparams = SFLR_params(**mydict)
    return myparams

def SFLRBeam(beamparams={'EI':EI, \
                         'mu':mu, \
                         'L':L}, \
             c=None,
             maxsize=ms, symname='Ubeam', \
             symlabel='beam', symsub=True, usez=True):
    if c is not None:
        beamparams['c'] = c
        myclass = BeamElement_v2
    else:
        myclass = BeamElement
    return myclass(beamparams, \
                   maxsize=maxsize, \
                   symlabel=symlabel, \
                   symname=symname, \
                   symsub=symsub, \
                   usez=usez)


#accel dimiensions in inches
hi = 5.0/8
di = 7.0/16

#convert to meters
h = hi*i_to_m
d = di*i_to_m
Ra = d/2.0#radius of the cylinder

#mass in kilograms
ma = 7.7/1000#i.e. 8 grams

#I about centroid
Ic = 1.0/12*ma*(3*Ra**2+h**2)
#parallel axis theorem
Ia = Ic + ma*(h/2.0)**2

#note that the accel is attached at 16.5 inches from the base
class accel_mass(RigidMass):
    def __init__(self, accel_params={'m':ma, 'L':d, \
                                     'r':Ra,'I':Ia}, \
                 maxsize=ms, symname='Uaccel', symlabel='accel', \
                 symsub=True, usez=True):
        return RigidMass.__init__(self, accel_params,\
                                  maxsize=maxsize, symlabel=symlabel, \
                                  symname=symname, symsub=symsub, \
                                  usez=usez)


class AVS_only_model(TMM.TMMSystem.ClampedFreeTMMSystem):
    def __init__(self, K, tau):
        self.avs=AngularVelocitySource({'K':K,'tau':tau}, \
                                       maxsize=ms, \
                                       symname='Uact', \
                                       symlabel='act', \
                                       unknownparams=['K','tau'])
        bodeout1 = {'input':'v', 'output':'th', 'type':'abs', \
                    'ind':self.avs, 'post':'', 'dof':1, \
                    'gain':180.0/pi*1024.0/360.0}
        self.bodeout1 = bodeout(**bodeout1)
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self,
                              [self.avs],
                              bodeouts=[self.bodeout1])
        

class SFLR_TMM_OL_model(TMM.TMMSystem.ClampedFreeTMMSystem):
    def __init__(self, xf, actvect=[], include_spring=True):
        t = xf[-1]
        A = w*t
        I = 1.0/12*w*t**3
        mu = A*rho
        EI = E*I

        self.spring = TorsionalSpringDamper({'k':xf[0],'c':xf[1]}, \
                                       maxsize=ms,symname='Usp', \
                                       symlabel='sp', \
                                       unknownparams=['k','c'])
        self.beam = SFLRBeam(beamparams={'EI':EI, 'L':L, 'mu':mu})
        if actvect:#an empty actvect means the actuator is unknown
            self.avs = AngularVelocitySource({'K':actvect[0], \
                                              'tau':actvect[1]}, \
                                             maxsize=ms, \
                                             symname='Uact', \
                                             symlabel='act')
            ind=2
        else:
            self.avs = AngularVelocitySource({'K':xf[2],'tau':xf[3]}, \
                                             maxsize=ms, \
                                             symname='Uact', \
                                             symlabel='act', \
                                             unknownparams=['K','tau'])
            ind=4
        self.accel_tip = accel_mass()#a rigid mass at the tip of the beam
        bodeout1={'input':'v', 'output':'atip', 'type':'abs', \
                  'ind':self.accel_tip, 'post':'accel', 'dof':0, \
                  'gain':xf[-2],'gainknown':False}
        if include_spring:
            b2_ind = self.spring
            my_list = [self.avs, self.spring, self.beam, self.accel_tip]
        else:
            b2_ind = self.avs
            my_list = [self.avs, self.beam, self.accel_tip]
        bodeout2={'input':'v', 'output':'th', 'type':'abs', \
                  'ind':b2_ind, 'post':'', 'dof':1, \
                  'gain':180.0/pi*1024.0/360.0}
        self.bodeout1=bodeout(**bodeout1)
        self.bodeout2=bodeout(**bodeout2)
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self,
                                                           my_list, \
                              bodeouts=[self.bodeout1, self.bodeout2])


def calc_beam_props(t=1.0/32*i_to_m, w=13.0/16*i_to_m, L=17.0*i_to_m, \
                    rho=7850.0, E=200*10**9):
    A = w*t
    I = 1.0/12*w*t**3
    mu = A*rho
    EI = E*I
    return {'EI':EI, 'mu':mu, 'L':L}

    
class Beam_Only_Model(TMM.TMMSystem.ClampedFreeTMMSystem):
    def __init__(self, beamparams=None):
        if beamparams is None:
            beamparams = calc_beam_props()
        self.beam = SFLRBeam(beamparams=beamparams)
        my_list = [self.beam]
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self,
                                                           my_list)


class system_w_beam(TMM.TMMSystem.ClampedFreeTMMSystem):
    def __init__(self, beamparams=None, c=None):
        if beamparams is None:
            beamparams = calc_beam_props()
        self.beam = SFLRBeam(beamparams=beamparams, c=c)
    
class Beam_w_Accel_Mass(TMM.TMMSystem.ClampedFreeTMMSystem):
    def __init__(self, beamparams=None):
        if beamparams is None:
            beamparams = calc_beam_props()
        self.beam = SFLRBeam(beamparams=beamparams)
        self.accel = accel_mass()
        my_list = [self.beam, self.accel]
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self, \
                                                           my_list)
    
class Beam_Base_Spring_Accel_Mass(system_w_beam):
    def __init__(self, beamparams=None, k=1000.0):
        system_w_beam.__init__(self, beamparams)
        self.clamp_spring = TorsionalSpringDamper({'k':k,'c':0}, \
                                                  maxsize=ms)
        self.accel = accel_mass()
        my_list = [self.clamp_spring, self.beam, self.accel]
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self, \
                                                           my_list)
        
class Two_Piece_Beam_w_Accel_Mass(TMM.TMMSystem.ClampedFreeTMMSystem):
    def __init__(self, beamparams=None):
        if beamparams is None:
            beamparams = calc_beam_props()
        L1 = (16.0+11.0/32)*25.4/1000
        L2 = (5.0/16)*25.4/1000
        d1 = copy.copy(beamparams)
        d1['L'] = L1
        d2 = copy.copy(beamparams)
        d2['L'] = L2
        self.beam = SFLRBeam(beamparams=d1)
        self.accel = accel_mass()
        self.tip_beam = SFLRBeam(beamparams=d2)
        my_list = [self.beam, self.accel, self.tip_beam]
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self, \
                                                           my_list)


class SFLR_TMM_OL_model_v2(system_w_beam):
    def __init__(self, xf, actvect=[], include_spring=True, \
                 beamparams=None):
        k = xf[-1]#default should be 45.0
        if beamparams is None:
            beamparams = calc_beam_props()
        beamparams['mu'] = beamparams['mu']*xf[-3]
        system_w_beam.__init__(self, beamparams)
        self.clamp_spring = TorsionalSpringDamper({'k':k,'c':0}, \
                                                  maxsize=ms)
        self.spring = TorsionalSpringDamper({'k':xf[0],'c':xf[1]}, \
                                       maxsize=ms,symname='Usp', \
                                       symlabel='sp', \
                                       unknownparams=['k','c'])
        if actvect:#an empty actvect means the actuator is unknown
            self.avs = AngularVelocitySource({'K':actvect[0], \
                                              'tau':actvect[1]}, \
                                             maxsize=ms, \
                                             symname='Uact', \
                                             symlabel='act')
            ind=2
        else:
            self.avs = AngularVelocitySource({'K':xf[2],'tau':xf[3]}, \
                                             maxsize=ms, \
                                             symname='Uact', \
                                             symlabel='act', \
                                             unknownparams=['K','tau'])
            ind=4
        self.accel_tip = accel_mass()#a rigid mass at the tip of the beam
        bodeout1={'input':'v', 'output':'atip', 'type':'abs', \
                  'ind':self.accel_tip, 'post':'accel', 'dof':0, \
                  'gain':xf[-2],'gainknown':False}
        if include_spring:
            b2_ind = self.spring
            my_list = [self.avs, self.spring, self.clamp_spring, \
                       self.beam, self.accel_tip]
        else:
            b2_ind = self.avs
            my_list = [self.avs, self.beam, self.accel_tip]
        bodeout2={'input':'v', 'output':'th', 'type':'abs', \
                  'ind':b2_ind, 'post':'', 'dof':1, \
                  'gain':180.0/pi*1024.0/360.0}
        self.bodeout1=bodeout(**bodeout1)
        self.bodeout2=bodeout(**bodeout2)
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self,
                                                           my_list, \
                              bodeouts=[self.bodeout1, self.bodeout2])


class AVS_and_Beam(system_w_beam):
    def __init__(self, K, tau, a_gain, beamparams=None, c=0.0):
        system_w_beam.__init__(self, beamparams, c=c)
        self.avs = AngularVelocitySource({'K':K,'tau':tau}, \
                                         maxsize=ms, \
                                         symname='Uact', \
                                         symlabel='act', \
                                         unknownparams=['K','tau'])
        bodeout1={'input':'v', 'output':'atip', 'type':'abs', \
                  'ind':self.beam, 'post':'accel', 'dof':0, \
                  'gain':a_gain,'gainknown':False}        
        bodeout2={'input':'v', 'output':'th', 'type':'abs', \
                  'ind':self.avs, 'post':'', 'dof':1, \
                  'gain':180.0/pi*1024.0/360.0}
        self.bodeout1=bodeout(**bodeout1)
        self.bodeout2=bodeout(**bodeout2)
        my_list = [self.avs, self.beam]
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self,
                                                           my_list, \
                              bodeouts=[self.bodeout1, self.bodeout2])
        
        
    
class SFLR_TMM_OL_model_v3(system_w_beam):
    def __init__(self, xf, actvect=[], include_spring=True, \
                 beamparams=None):
        k = xf[-1]#default should be 45.0
        c_beam = xf[-4]
        if beamparams is None:
            beamparams = calc_beam_props()
        beamparams['mu'] = beamparams['mu']*xf[-3]
        system_w_beam.__init__(self, beamparams, c=c_beam)
        self.clamp_spring = TorsionalSpringDamper({'k':k,'c':0}, \
                                                  maxsize=ms)
        self.spring = TorsionalSpringDamper({'k':xf[0],'c':xf[1]}, \
                                       maxsize=ms,symname='Usp', \
                                       symlabel='sp', \
                                       unknownparams=['k','c'])
        if actvect:#an empty actvect means the actuator is unknown
            self.avs = AngularVelocitySource({'K':actvect[0], \
                                              'tau':actvect[1]}, \
                                             maxsize=ms, \
                                             symname='Uact', \
                                             symlabel='act')
            ind=2
        else:
            self.avs = AngularVelocitySource({'K':xf[2],'tau':xf[3]}, \
                                             maxsize=ms, \
                                             symname='Uact', \
                                             symlabel='act', \
                                             unknownparams=['K','tau'])
            ind=4
        #self.accel_tip = accel_mass()#a rigid mass at the tip of the beam
        #end_ind = self.accel_tip
        end_ind = self.beam
        bodeout1={'input':'v', 'output':'atip', 'type':'abs', \
                  'ind':end_ind, 'post':'accel', 'dof':0, \
                  'gain':xf[-2],'gainknown':False}
        if include_spring:
            b2_ind = self.spring
            my_list = [self.avs, self.spring, self.clamp_spring, \
                       self.beam]#, self.accel_tip]
        else:
            b2_ind = self.avs
            my_list = [self.avs, self.beam, self.accel_tip]
        bodeout2={'input':'v', 'output':'th', 'type':'abs', \
                  'ind':b2_ind, 'post':'', 'dof':1, \
                  'gain':180.0/pi*1024.0/360.0}
        self.bodeout1=bodeout(**bodeout1)
        self.bodeout2=bodeout(**bodeout2)
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self,
                                                           my_list, \
                              bodeouts=[self.bodeout1, self.bodeout2])

def a_massage(a_bode, fvect, lim2=100.0, flow=4.0, fhigh=9.0):
    phase = a_bode.phase
    b1 = a_bode.phase < -100
    b2 = fvect < flow
    b = b1 & b2
    phase[b] += 360.0
    
    b3 = a_bode.phase > lim2
    b4 = fvect > fhigh
    b5 = b3 & b4
    phase[b5] -= 360.0
    a_bode.phase = phase
    return a_bode

def phase_shift(bode, max_phase=120):
    phase = bode.phase
    b1 = phase > max_phase
    phase[b1] -= 360.0
    bode.phase = phase
    return bode

class SFLR_TMM_OL_model_v4(TMM.TMMSystem.ClampedFreeTMMSystem):
    def _build_beam_params_dict(self, params=None, Lparam='L'):
        if params is None:
            params = self.params
        L = getattr(params, Lparam)
        beamdict = {'EI':params.EI, 'mu':params.mu, \
                      'L':L}
        return beamdict
            
    def _create_beam(self, params=None, attr='beam', Lparam='L'):
        if params is None:
            params = self.params
        beamparams = self._build_beam_params_dict(params, \
                                                  Lparam=Lparam)
        temp_beam = SFLRBeam(beamparams, c=params.c_beam)
        setattr(self, attr, copy.copy(temp_beam))

    def _create_beam2(self, params=None, attr='beam2', Lparam='L2'):
        self._create_beam(params=params, attr=attr, Lparam=Lparam)

    def _create_beams(self):
        self._create_beam()
        self._create_beam2()

    def _create_accel_mass(self):
        accel_params = {'m':self.params.a_m, \
                        'L':self.params.a_L, \
                        'r':self.params.a_r, \
                        'I':self.params.a_I}
        self.accel = accel_mass(accel_params=accel_params)

    def _system_init(self):
        """Be sure that self.list, self.bodeout1, and self.bodeout2
        are defined before calling this method - the last step in the
        init."""
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self,
                                                           self.list, \
                             bodeouts=[self.bodeout1, self.bodeout2])
        

    def _create_AVS(self):
        self.avs = AngularVelocitySource({'K':self.params.K_act, \
                                          'tau':self.params.tau}, \
                                         maxsize=ms, \
                                         symname='Uact', \
                                         symlabel='act', \
                                         unknownparams=['K','tau'])

    def _create_spring(self):
        k = self.params.k_spring
        c = self.params.c_spring
        self.spring = TorsionalSpringDamper({'k':k,'c':c}, \
                                            maxsize=ms,symname='Usp', \
                                            symlabel='sp', \
                                            unknownparams=['k','c'])

    def _create_clamp_spring(self):
        k = self.params.k_clamp
        c = self.params.c_clamp
        self.clamp_spring = TorsionalSpringDamper({'k':k,'c':c}, \
                                                  maxsize=ms)

    def _set_accel_ind(self):
        #self.accel_ind = self.list[-1]
        self.accel_ind = self.accel
        
    def _create_bode_outs(self):
        """Note that self.list and self.b2_ind must be defined before
        self._create_bode_outs is called.  self.b2_ind is the index of
        the element just before the encoder measurement.  If the
        actuator is rigid, self.b2_ind should probably be set to
        self.avs.  If there is flexibility in the actuator,
        self.b2_ind should probaby be set to the torsional
        spring/damper following the avs."""
        self._set_accel_ind()
        a_gain = self.params.a_gain
        bodeout1={'input':'v', 'output':'atip', 'type':'abs', \
                  'ind':self.accel_ind, 'post':'accel', 'dof':0, \
                  'gain':a_gain,'gainknown':False}
        bodeout2={'input':'v', 'output':'th', 'type':'abs', \
                  'ind':self.b2_ind, 'post':'', 'dof':1, \
                  'gain':180.0/pi*1024.0/360.0}
        self.bodeout1=bodeout(**bodeout1)
        self.bodeout2=bodeout(**bodeout2)
        
    def __init__(self, myparams):
        self.params = myparams

    def calc_bodes(self, fvect):
        self.fvect = fvect
        bodes = self.BodeResponse(fvect)
        theta_v_bode = bodes[1]
        theta_v_bode.seedphase = -90.0
        theta_v_bode.seedfreq = 1.2
        theta_v_bode.PhaseMassage(fvect)

        accel_v_bode = bodes[0]
        accel_v_bode.seedfreq = 2.0
        accel_v_bode.seedphase = 30.0
        accel_v_bode.PhaseMassage(fvect)
        accel_v_bode = a_massage(accel_v_bode, fvect, lim2=-90)
        accel_v_bode = phase_shift(accel_v_bode)

        a_th_bode = accel_v_bode.__div__(theta_v_bode)
        a_th_bode.seedfreq = 10.0
        a_th_bode.seedphase = 0.0
        a_th_bode.PhaseMassage(fvect)
        a_th_bode = a_massage(a_th_bode, fvect)
        self.bodes = bodes
        self.theta_v_bode = theta_v_bode
        self.accel_v_bode = accel_v_bode
        self.accel_theta_bode = a_th_bode
        return theta_v_bode, accel_v_bode, a_th_bode

    def plot_bodes(self, startfi=1, clear=False, **kwargs):
        f = self.fvect
        rwkbode.GenBodePlot(startfi, f, self.theta_v_bode, clear=clear, \
                            **kwargs)
        rwkbode.GenBodePlot(startfi+1, f, self.accel_v_bode, \
                            clear=clear, **kwargs)
        rwkbode.GenBodePlot(startfi+2, f, self.accel_theta_bode, \
                            clear=clear, **kwargs)

    def calc_and_plot_bodes(self, fvect, startfi=1, clear=False, **kwargs):
        t1 = time.time()
        self.calc_bodes(fvect)
        t2 = time.time()
        print('self.BodeResponse time='+str(t2-t1))
        self.plot_bodes(startfi, clear=clear, **kwargs)
        
class Rigid_Actuator_Model(SFLR_TMM_OL_model_v4):
    """Model for the SFLR with only the AVS and Beam (i.e. the
    actuator cannot be back-driven)."""
    def __init__(self, myparams):
        SFLR_TMM_OL_model_v4.__init__(self, myparams)
        self._create_beam()
        self._create_AVS()
        self.b2_ind = self.avs
        self.list = [self.avs, self.beam]
        self._create_bode_outs()

        return self._system_init()

class AVS_Beam_Accel_Mass(Rigid_Actuator_Model):
    def __init__(self, myparams):
        SFLR_TMM_OL_model_v4.__init__(self, myparams)
        self._create_beam()
        self._create_AVS()
        self._create_accel_mass()
        self.b2_ind = self.avs
        self.list = [self.avs, self.beam, self.accel]
        self._create_bode_outs()
        return self._system_init()


class AVS_Spring_Beam_Accel(AVS_Beam_Accel_Mass):
    def __init__(self, myparams):
        SFLR_TMM_OL_model_v4.__init__(self, myparams)
        self._create_beam()
        self._create_AVS()
        self._create_accel_mass()
        self._create_spring()
        self.b2_ind = self.spring
        self.list = [self.avs, self.spring, self.beam, self.accel]
        self._create_bode_outs()
        return self._system_init()
    

class AVS_Spring_Clamp_Spring_Beam_Accel(AVS_Spring_Beam_Accel):
    def __init__(self, myparams):
        SFLR_TMM_OL_model_v4.__init__(self, myparams)
        self._create_beam()
        self._create_AVS()
        self._create_accel_mass()
        self._create_spring()
        self._create_clamp_spring()
        self.b2_ind = self.spring
        self.list = [self.avs, self.spring, self.clamp_spring, \
                     self.beam, self.accel]
        self._create_bode_outs()
        return self._system_init()



class unactuated_beam_only(SFLR_TMM_OL_model_v4):
    def _system_init(self):
        """Be sure that self.list, self.bodeout1, and self.bodeout2
        are defined before calling this method - the last step in the
        init."""
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self,
                                                           self.list)
    
    def __init__(self, myparams):
        SFLR_TMM_OL_model_v4.__init__(self, myparams)
        self._create_beam()
        self.list = [self.beam]
        self._system_init()
        
    def find_nat_freqs(self, flist=[2.5, 17.9]):
        slist = [2.0j*pi*item for item in flist]
        nat_freqs = []
        for s in slist:
            eig_array = self.FindEig(s, disp=0)
            cur_f = eig_array[-1]/(2*pi)
            nat_freqs.append(cur_f)
        self.nat_freqs = nat_freqs
        return nat_freqs

    def print_nat_freqs(self, label=''):
        print('%s natural frequencies:' % label)
        print('   f1=%0.4f' % self.nat_freqs[0])
        print('   f2=%0.4f' % self.nat_freqs[1])
        
    def calc_det_vector(self, f):
        s = 2.0j*pi*f
        t1 = time.time()
        comp_list = [self.EigError(item) for item in s]
        t2 = time.time()
        print('calc_det_vector time=%s' % (t2-t1))
        comp_array = array(comp_list)
        mag_array = abs(comp_array)
        self.det_array = mag_array
        return mag_array
        
    def plot_det_vector(self, f, fi=1, clear=False, ylim=[0,20]):
        self.calc_det_vector(f)
        pylab.figure(fi)
        if clear:
            pylab.clf()
        pylab.plot(f, self.det_array)
        if ylim:
            pylab.ylim(ylim)
        pylab.xlabel('Freq. (Hz)')
        pylab.ylabel('Determinant')
        

class unactuated_beam_w_accel(unactuated_beam_only):
    def __init__(self, myparams):
        SFLR_TMM_OL_model_v4.__init__(self, myparams)
        self._create_beam()
        self._create_accel_mass()
        self.list = [self.beam, self.accel]
        self._system_init()
        
class unactuated_two_piece_beam_w_accel(unactuated_beam_w_accel):
    def __init__(self, myparams):
        SFLR_TMM_OL_model_v4.__init__(self, myparams)
        self._create_beams()
        self._create_accel_mass()
        self.list = [self.beam, self.accel, self.beam2]
        self._system_init()
        
class unactuated_beam_w_accel_and_clamp_spring(unactuated_beam_w_accel):
    def __init__(self, myparams):
        SFLR_TMM_OL_model_v4.__init__(self, myparams)
        self._create_beam()
        self._create_accel_mass()
        self._create_clamp_spring()
        self.list = [self.clamp_spring, self.beam, self.accel]
        self._system_init()
        

class unactuated_two_piece_beam_w_accel_and_clamp(unactuated_two_piece_beam_w_accel):
    def __init__(self, myparams):
        SFLR_TMM_OL_model_v4.__init__(self, myparams)
        self._create_beams()
        self._create_accel_mass()
        self._create_clamp_spring()
        self.list = [self.clamp_spring, self.beam, \
                     self.accel, self.beam2]
        self._system_init()
    


class model_v_2_0(AVS_Spring_Clamp_Spring_Beam_Accel):
    def __init__(self, myparams):
        SFLR_TMM_OL_model_v4.__init__(self, myparams)
        self._create_beams()
        self._create_AVS()
        self._create_accel_mass()
        self._create_spring()
        self._create_clamp_spring()
        self.b2_ind = self.spring
        self.list = [self.avs, self.spring, self.clamp_spring, \
                     self.beam, self.accel, self.beam2]
        self._create_bode_outs()
        return self._system_init()

class cheating_model(AVS_Spring_Clamp_Spring_Beam_Accel):
    def __init__(self, myparams):
        SFLR_TMM_OL_model_v4.__init__(self, myparams)
        self._create_beams()
        self._create_AVS()
        self._create_accel_mass()
        #self._create_spring()
        self._create_clamp_spring()
        self.b2_ind = self.clamp_spring
        self.list = [self.avs, self.clamp_spring, \
                     self.beam, self.accel, self.beam2]
        self._create_bode_outs()
        return self._system_init()