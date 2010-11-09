from __future__ import division
from scipy import sinh, cosh, sin, cos, real, imag, shape, arange, pi, zeros, array, eye, transpose, conj, c_, poly, vectorize, dot, randn, sum, squeeze, sqrt
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
    GthNum = Gth.Gth.num(s)
    GthDen = Gth.Gth.den(s)
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
    a_7 = 1.0E+0*a_2*GthNum*K_act*a_3*a_4*a_5*a_6
    a_8 = 1/(c_spring*s+k_spring)
    a_9 = beta1**3
    a_10 = L1**3
    a_11 = 1/a_10
    a_12 = beta1**2
    a_13 = 1/(c_clamp*s+k_clamp)
    a_14 = a_1*b_I-b_m*b_r*a_1*(b_L-b_r)
    a_15 = a_2*GthNum*K_act*a_3*a_4*a_5*a_13*a_14*a_6+a_7
    a_16 = L1**2
    a_17 = 1/a_16
    a_18 = 1/L1
    a_19 = -1.0E+0*beta1*d2_1*a_2*GthNum*K_act*a_3*a_4*a_5*a_14*a_6*a_18 \
        -a_12*d3_1*EI*a_15*a_17-1.0E+0*a_9*d4_1*a_2*GthNum*K_act*a_3* \
        a_4*a_5*b_L*EI*a_6*a_11+1.0E+0*b_m*b_r*d1_1*a_2*GthNum*K_act*a_3 \
        *s*a_5*a_6
    a_20 = 1/beta1
    a_21 = 1/EI
    a_22 = 1/a_12
    a_23 = -1.0E+0*a_22*b_m*b_r*d3_1*a_2*GthNum*K_act*a_3*s*a_5*a_21*a_6 \
        *a_16+1.0E+0*a_20*d4_1*a_2*GthNum*K_act*a_3*a_4*a_5*a_14*a_21 \
        *a_6*L1+1.0E+0*beta1*d2_1*a_2*GthNum*K_act*a_3*a_4*a_5*b_L*a_6* \
        a_18+d1_1*a_15
    a_24 = 1/a_9
    a_25 = -1.0E+0*a_24*b_m*b_r*d2_1*a_2*GthNum*K_act*a_3*s*a_5*a_21*a_6 \
        *a_10+1.0E+0*a_22*d3_1*a_2*GthNum*K_act*a_3*a_4*a_5*a_14*a_21 \
        *a_6*a_16+a_20*d4_1*a_15*L1+1.0E+0*d1_1*a_2*GthNum*K_act*a_3*a_4 \
        *a_5*b_L*a_6
    a_26 = a_m*a_1*a_25+a_m*a_r*a_1*a_23+1.0E+0*a_19
    a_27 = beta2**3
    a_28 = a_L*a_23
    a_29 = 1.0E+0*a_25
    a_30 = a_29+a_28
    a_31 = 1/L2**3
    a_32 = beta2**2
    a_33 = 1/L2**2
    a_34 = a_L-a_r
    a_35 = a_1*a_I-a_m*a_r*a_1*a_34
    a_36 = -a_m*a_1*a_34*a_25+a_35*a_23+1.0E+0*(-1.0E+0*a_20*b_m*b_r*d4_1 \
        *a_2*GthNum*K_act*a_3*s*a_5*a_6*L1+beta1*d2_1*EI*a_15*a_18+1.0E+0 \
        *a_12*d3_1*a_2*GthNum*K_act*a_3*a_4*a_5*b_L*EI*a_6*a_17+1.0E+0 \
        *d1_1*a_2*GthNum*K_act*a_3*a_4*a_5*a_14*a_6)-a_L*a_19
    a_37 = 1/L2
    a_38 = beta2*d2_2*a_36*a_37+1.0E+0*a_32*d3_2*EI*a_23*a_33+a_27*d4_2 \
        *EI*a_30*a_31-d1_2*a_26
    a_39 = 1.0E+0*beta1*d2_1*b_L*a_18+a_12*d3_1*a_13*b_L*EI*a_17+1.0E+0 \
        *d1_1
    a_40 = -1.0E+0*a_22*d3_1*a_21*a_16-1.0E+0*a_20*d4_1*b_L*a_21*L1-d1_1 \
        *a_13*b_L
    a_41 = -1.0E+0*a_24*d2_1*a_21*a_10-1.0E+0*a_22*d3_1*b_L*a_21*a_16- \
        a_20*d4_1*a_13*b_L*L1
    a_42 = -a_m*a_1*a_34*a_41+a_35*a_40+1.0E+0*(-1.0E+0*a_20*d4_1*L1-beta1 \
        *d2_1*a_13*b_L*EI*a_18-1.0E+0*d1_1*b_L)-a_L*a_39
    a_43 = 1.0E+0*a_41+a_L*a_40
    a_44 = 1/beta2
    a_45 = a_m*a_1*a_41+a_m*a_r*a_1*a_40+1.0E+0*a_39
    a_46 = a_8*a_14*a_6+1.0E+0
    a_47 = a_13*a_46+1.0E+0*a_8*a_6
    a_48 = -1.0E+0*beta1*d2_1*a_46*a_18-a_12*d3_1*EI*a_47*a_17-1.0E+0* \
        a_9*d4_1*a_8*b_L*EI*a_6*a_11+1.0E+0*b_m*b_r*d1_1*a_1*a_8*a_6
    a_49 = -1.0E+0*a_22*b_m*b_r*d3_1*a_1*a_8*a_21*a_6*a_16+1.0E+0*a_20 \
        *d4_1*a_21*a_46*L1+1.0E+0*beta1*d2_1*a_8*b_L*a_6*a_18+d1_1*a_47
    a_50 = -1.0E+0*a_24*b_m*b_r*d2_1*a_1*a_8*a_21*a_6*a_10+1.0E+0*a_22 \
        *d3_1*a_21*a_46*a_16+a_20*d4_1*a_47*L1+1.0E+0*d1_1*a_8*b_L*a_6
    a_51 = a_m*a_1*a_50+a_m*a_r*a_1*a_49+1.0E+0*a_48
    a_52 = 1.0E+0*a_50+a_L*a_49
    a_53 = -a_m*a_1*a_34*a_50+a_35*a_49+1.0E+0*(-1.0E+0*a_20*b_m*b_r*d4_1 \
        *a_1*a_8*a_6*L1+beta1*d2_1*EI*a_47*a_18+1.0E+0*a_12*d3_1*a_8 \
        *b_L*EI*a_6*a_17+1.0E+0*d1_1*a_46)-a_L*a_48
    a_54 = beta2*d2_2*a_53*a_37+1.0E+0*a_32*d3_2*EI*a_49*a_33+a_27*d4_2 \
        *EI*a_52*a_31-d1_2*a_51
    a_55 = -beta2*d2_2*a_42*a_37-1.0E+0*a_32*d3_2*EI*a_40*a_33-a_27*d4_2 \
        *EI*a_43*a_31+d1_2*a_45
    a_56 = -a_44*d4_2*a_51*L2+1.0E+0*beta2*d2_2*EI*a_49*a_37+a_32*d3_2 \
        *EI*a_52*a_33+d1_2*a_53
    a_57 = 1/(a_55*a_56+a_54*(-a_44*d4_2*a_45*L2+1.0E+0*beta2*d2_2*EI* \
        a_40*a_37+a_32*d3_2*EI*a_43*a_33+d1_2*a_42))
    a_58 = a_44*d4_2*a_26*L2-1.0E+0*beta2*d2_2*EI*a_23*a_37-a_32*d3_2* \
        EI*a_30*a_33-d1_2*a_36
    a_59 = a_55*a_58*a_57+a_38*(a_44*d4_2*a_45*L2-1.0E+0*beta2*d2_2*EI \
        *a_40*a_37-a_32*d3_2*EI*a_43*a_33-d1_2*a_42)*a_57
    RESULT = a_gain*a_1*(a_43*(a_38*a_56*a_57+a_54*a_58*a_57)+a_52*a_59 \
        +a_29+a_28)/(enc_gain*(1.0E+0*a_8*a_6*a_59+a_7))
    return RESULT
