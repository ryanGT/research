from __future__ import division
from scipy import sinh, cosh, sin, cos, real, imag, shape, arange, pi, zeros, array, eye, transpose, conj, c_, poly, vectorize, dot, randn, sum, squeeze, sqrt, poly1d
#from scipy import matrixmultiply
#from numpy.lib.scimath import sqrt
#from LinearAlgebra import inverse
#import pylab
#import rwkbode
#reload(rwkbode)
#from scipy.optimize import newton, fmin, fminbound
#pylab.ioff()

import sys, os, copy, time

import SFLR_TMM
reload(SFLR_TMM)

import add_design_dir

import Gth
reload(Gth)

zsin=sin
zcos=cos
zsinh=sinh
zcosh=cosh

def Bodes(s, params):
    mynum = squeeze(Gth.Gth.num)
    myden = squeeze(Gth.Gth.den)
    try:
        GthNum = mynum(s)
        GthDen = myden(s)
    except:
        mynum = poly1d(mynum)
        myden = poly1d(myden)
        GthNum = mynum(s)
        GthDen = myden(s)
        
    if type(params) == dict:
        params = SFLR_TMM.SFLR_params(**params)
    EI = params.EI
    L1 = params.L
    L2 = params.L2
    mu = params.mu
    c_beam = params.c_beam
    if c_beam > 0.0:
        EI = EI*(1.0+c_beam*s)
##     beta = pow((-1*s*s*L**4*mu/EI),0.25)
##     d1 = 0.5*(cos(beta)+cosh(beta))
##     d2 = 0.5*(sinh(beta)-sin(beta))
##     d3 = 0.5*(cosh(beta)-cos(beta))
##     d4 = 0.5*(sin(beta)+sinh(beta))
    beta1 = pow((-1*s*s*L1**4*mu/EI),0.25)
    d1_1 = 0.5*(cos(beta1)+cosh(beta1))
    d2_1 = 0.5*(sinh(beta1)-sin(beta1))
    d3_1 = 0.5*(cosh(beta1)-cos(beta1))
    d4_1 = 0.5*(sin(beta1)+sinh(beta1))
    beta2 = pow((-1*s*s*L2**4*mu/EI),0.25)
    d1_2 = 0.5*(cos(beta2)+cosh(beta2))
    d2_2 = 0.5*(sinh(beta2)-sin(beta2))
    d3_2 = 0.5*(cosh(beta2)-cos(beta2))
    d4_2 = 0.5*(sin(beta2)+sinh(beta2))
    a_m = params.a_m
    a_L = params.a_L
    a_I = params.a_I
    a_r = params.a_r
    b_m = params.b_m
    b_L = params.b_L
    b_I = params.b_I
    b_r = params.b_r
    k_spring = params.k_spring
    c_spring = params.c_spring
    k_clamp = params.k_clamp
    c_clamp = params.c_clamp
    K_act = params.K_act
    tau = params.tau
    a_gain = params.a_gain
    p_act1 = params.p_act1
    p_act2 = params.p_act2
    z_act = params.z_act
    H = params.H
    I = 1.0j
    enc_gain = 180.0/pi*1024.0/360.0
    #--------------------------------
    a_1 = s**2
    a_2 = 1/GthDen
    a_3 = sqrt(p_act1**2+4*pi**2)
    a_4 = 1/s
    a_5 = 1/(s+p_act1)
    a_6 = 1/(a_2*GthNum*K_act*a_3*a_4*a_5*H+1.0E+0)
    a_7 = 1/(c_clamp*s+k_clamp)
    a_8 = a_1*b_I-b_m*b_r*a_1*(b_L-b_r)
    a_9 = a_2*GthNum*K_act*a_3*a_4*a_5*a_7*a_8*a_6+1.0E+0*a_2*GthNum*K_act \
        *a_3*a_4*a_5*a_6
    a_10 = 1/L1
    a_11 = 1/beta1
    a_12 = 1/EI
    a_13 = beta1**2
    a_14 = 1/a_13
    a_15 = L1**2
    a_16 = -1.0E+0*a_14*b_m*b_r*d3_1*a_2*GthNum*K_act*a_3*s*a_5*a_12*a_6 \
        *a_15+1.0E+0*a_11*d4_1*a_2*GthNum*K_act*a_3*a_4*a_5*a_8*a_12* \
        a_6*L1+1.0E+0*beta1*d2_1*a_2*GthNum*K_act*a_3*a_4*a_5*b_L*a_6*a_10 \
        +d1_1*a_9
    a_17 = a_L*a_16
    a_18 = beta1**3
    a_19 = 1/a_18
    a_20 = L1**3
    a_21 = -1.0E+0*a_19*b_m*b_r*d2_1*a_2*GthNum*K_act*a_3*s*a_5*a_12*a_6 \
        *a_20+1.0E+0*a_14*d3_1*a_2*GthNum*K_act*a_3*a_4*a_5*a_8*a_12* \
        a_6*a_15+a_11*d4_1*a_9*L1+1.0E+0*d1_1*a_2*GthNum*K_act*a_3*a_4* \
        a_5*b_L*a_6
    a_22 = 1.0E+0*a_21
    a_23 = 1/(c_spring*s+k_spring)
    a_24 = a_23*a_8*a_6+1.0E+0
    a_25 = a_7*a_24+1.0E+0*a_23*a_6
    a_26 = -1.0E+0*a_14*b_m*b_r*d3_1*a_1*a_23*a_12*a_6*a_15+1.0E+0*a_11 \
        *d4_1*a_12*a_24*L1+1.0E+0*beta1*d2_1*a_23*b_L*a_6*a_10+d1_1*a_25
    a_27 = -1.0E+0*a_19*b_m*b_r*d2_1*a_1*a_23*a_12*a_6*a_20+1.0E+0*a_14 \
        *d3_1*a_12*a_24*a_15+a_11*d4_1*a_25*L1+1.0E+0*d1_1*a_23*b_L*a_6
    a_28 = 1.0E+0*a_27+a_L*a_26
    a_29 = 1/a_20
    a_30 = 1/a_15
    a_31 = -1.0E+0*beta1*d2_1*a_2*GthNum*K_act*a_3*a_4*a_5*a_8*a_6*a_10 \
        -a_13*d3_1*EI*a_9*a_30-1.0E+0*a_18*d4_1*a_2*GthNum*K_act*a_3*a_4 \
        *a_5*b_L*EI*a_6*a_29+1.0E+0*b_m*b_r*d1_1*a_2*GthNum*K_act*a_3 \
        *s*a_5*a_6
    a_32 = a_m*a_1*a_21+a_m*a_r*a_1*a_16+1.0E+0*a_31
    a_33 = beta2**3
    a_34 = a_22+a_17
    a_35 = 1/L2**3
    a_36 = beta2**2
    a_37 = 1/L2**2
    a_38 = a_L-a_r
    a_39 = a_1*a_I-a_m*a_r*a_1*a_38
    a_40 = -a_m*a_1*a_38*a_21+a_39*a_16+1.0E+0*(-1.0E+0*a_11*b_m*b_r*d4_1 \
        *a_2*GthNum*K_act*a_3*s*a_5*a_6*L1+beta1*d2_1*EI*a_9*a_10+1.0E+0 \
        *a_13*d3_1*a_2*GthNum*K_act*a_3*a_4*a_5*b_L*EI*a_6*a_30+1.0E+0 \
        *d1_1*a_2*GthNum*K_act*a_3*a_4*a_5*a_8*a_6)-a_L*a_31
    a_41 = 1/L2
    a_42 = beta2*d2_2*a_40*a_41+1.0E+0*a_36*d3_2*EI*a_16*a_37+a_33*d4_2 \
        *EI*a_34*a_35-d1_2*a_32
    a_43 = 1.0E+0*beta1*d2_1*b_L*a_10+a_13*d3_1*a_7*b_L*EI*a_30+1.0E+0 \
        *d1_1
    a_44 = -1.0E+0*a_14*d3_1*a_12*a_15-1.0E+0*a_11*d4_1*b_L*a_12*L1-d1_1 \
        *a_7*b_L
    a_45 = -1.0E+0*a_19*d2_1*a_12*a_20-1.0E+0*a_14*d3_1*b_L*a_12*a_15- \
        a_11*d4_1*a_7*b_L*L1
    a_46 = -a_m*a_1*a_38*a_45+a_39*a_44+1.0E+0*(-1.0E+0*a_11*d4_1*L1-beta1 \
        *d2_1*a_7*b_L*EI*a_10-1.0E+0*d1_1*b_L)-a_L*a_43
    a_47 = 1.0E+0*a_45+a_L*a_44
    a_48 = 1/beta2
    a_49 = a_m*a_1*a_45+a_m*a_r*a_1*a_44+1.0E+0*a_43
    a_50 = -1.0E+0*beta1*d2_1*a_24*a_10-a_13*d3_1*EI*a_25*a_30-1.0E+0* \
        a_18*d4_1*a_23*b_L*EI*a_6*a_29+1.0E+0*b_m*b_r*d1_1*a_1*a_23*a_6
    a_51 = a_m*a_1*a_27+a_m*a_r*a_1*a_26+1.0E+0*a_50
    a_52 = -a_m*a_1*a_38*a_27+a_39*a_26+1.0E+0*(-1.0E+0*a_11*b_m*b_r*d4_1 \
        *a_1*a_23*a_6*L1+beta1*d2_1*EI*a_25*a_10+1.0E+0*a_13*d3_1*a_23 \
        *b_L*EI*a_6*a_30+1.0E+0*d1_1*a_24)-a_L*a_50
    a_53 = beta2*d2_2*a_52*a_41+1.0E+0*a_36*d3_2*EI*a_26*a_37+a_33*d4_2 \
        *EI*a_28*a_35-d1_2*a_51
    a_54 = -beta2*d2_2*a_46*a_41-1.0E+0*a_36*d3_2*EI*a_44*a_37-a_33*d4_2 \
        *EI*a_47*a_35+d1_2*a_49
    a_55 = -a_48*d4_2*a_51*L2+1.0E+0*beta2*d2_2*EI*a_26*a_41+a_36*d3_2 \
        *EI*a_28*a_37+d1_2*a_52
    a_56 = 1/(a_54*a_55+a_53*(-a_48*d4_2*a_49*L2+1.0E+0*beta2*d2_2*EI* \
        a_44*a_41+a_36*d3_2*EI*a_47*a_37+d1_2*a_46))
    a_57 = a_48*d4_2*a_32*L2-1.0E+0*beta2*d2_2*EI*a_16*a_41-a_36*d3_2* \
        EI*a_34*a_37-d1_2*a_40
    RESULT = a_gain*a_1*(a_47*(a_42*a_55*a_56+a_53*a_57*a_56)+a_28*(a_54 \
        *a_57*a_56+a_42*(a_48*d4_2*a_49*L2-1.0E+0*beta2*d2_2*EI*a_44* \
        a_41-a_36*d3_2*EI*a_47*a_37-d1_2*a_46)*a_56)+a_22+a_17)
    return RESULT
