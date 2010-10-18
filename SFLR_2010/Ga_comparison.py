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
AccelFB_bode_opts = SFLR_bode_options.AccelFB_Bode_opts

class AccelFB_System(object):
    def __init__(self, Ga):
        self.Ga = Ga
        self.ROM = SFLR_2010.ROM_model.create_ROM_model(AccelFB_Bode_opts, \
                                                        Accel_TMM_model, \
                                                        Gth.Gth, \
                                                        Ga)
        self.update_ROM(Ga)

        
    def update_ROM(self, Ga=None):
        if Ga is not None:
            self.ROM.Ga = Ga
        self.ROM.build_TFs()
        self.t_s = self.ROM.theta_TF.simplify()
        self.a_s = self.ROM.accel_TF.simplify()


    def lsim(self, u, t, calca=True):
        self.th_a_s = self.t_s.lsim(u,t)
        if calca:
            self.y_a_s = self.a_s.lsim(u,t)
        else:
            y_a_s = zeros_like(th_a_s)
        return th_a_s, y_a_s
