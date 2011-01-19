from __future__ import division
from scipy import *

def Unc_row0(s, params, Gth):
    if type(params) == dict:
        import SFLR_TMM
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
    #--------------------------------
    x0 = s**2
    x1 = 1/EI
    x2 = 1/beta1
    x3 = c_clamp*s
    x4 = k_clamp + x3
    x5 = 1/x4
    x6 = L1**2
    x7 = x2**2
    x8 = x6**(-1/2)
    x9 = b_r**2
    x10 = x8**(-3)
    x11 = x7**(3/2)
    x12 = x0**(-1/2)
    x13 = p_act1 + s
    x14 = 1/x13
    x15 = 2*pi*I
    x16 = p_act1 + x15
    x17 = abs(x16)
    x18 = Gth(s)*H*K_act*x12*x14*x17
    x19 = 1.0 + x18
    x20 = 1/x19
    x21 = c_spring*s
    x22 = k_spring + x21
    x23 = 1/x22
    Unc00 = d1_1 + a_L*beta1*d2_1*x8 + a_L*b_m*b_r*d1_1*x0*x5 - a_L*b_L*b_m*d1_1*x0*x5 - b_m*d2_1*x0*x1*x10*x11 + L1*b_m*b_r*d4_1*x0*x2*x5 + b_m*b_r*d3_1*x0*x1*x6*x7 - L1*b_L*b_m*d4_1*x0*x2*x5 - a_L*b_m*d3_1*x0*x1*x6*x7 - b_L*b_m*d3_1*x0*x1*x6*x7 + L1*a_L*b_m*b_r*d4_1*x0*x1*x2 - L1*a_L*b_L*b_m*d4_1*x0*x1*x2
    Unc01 = a_L*d1_1 + b_L*d1_1 + L1*d4_1*x2 + a_L*b_I*d1_1*x0*x5 + a_L*b_L*beta1*d2_1*x8 + L1*b_I*d4_1*x0*x2*x5 + a_L*b_m*d1_1*x0*x5*x9 + b_I*d3_1*x0*x1*x6*x7 + L1*a_L*b_I*d4_1*x0*x1*x2 + L1*b_m*d4_1*x0*x2*x5*x9 + b_m*d3_1*x0*x1*x6*x7*x9 - a_L*b_L*b_m*b_r*d1_1*x0*x5 - b_m*b_r*d2_1*x0*x1*x10*x11 + L1*a_L*b_m*d4_1*x0*x1*x2*x9 - L1*b_L*b_m*b_r*d4_1*x0*x2*x5 - a_L*b_m*b_r*d3_1*x0*x1*x6*x7 - b_L*b_m*b_r*d3_1*x0*x1*x6*x7 - L1*a_L*b_L*b_m*b_r*d4_1*x0*x1*x2
    Unc02 = a_L*d1_1*x5 + L1*d4_1*x2*x5 + a_L*d1_1*x20*x23 + b_L*d1_1*x20*x23 + d3_1*x1*x6*x7 + L1*a_L*d4_1*x1*x2 + L1*d4_1*x2*x20*x23 + a_L*b_I*d1_1*x0*x20*x23*x5 + a_L*b_L*beta1*d2_1*x20*x23*x8 + L1*b_I*d4_1*x0*x2*x20*x23*x5 + a_L*b_m*d1_1*x0*x20*x23*x5*x9 + b_I*d3_1*x0*x1*x20*x23*x6*x7 + L1*a_L*b_I*d4_1*x0*x1*x2*x20*x23 + L1*b_m*d4_1*x0*x2*x20*x23*x5*x9 + b_m*d3_1*x0*x1*x20*x23*x6*x7*x9 - a_L*b_L*b_m*b_r*d1_1*x0*x20*x23*x5 - b_m*b_r*d2_1*x0*x1*x10*x11*x20*x23 + L1*a_L*b_m*d4_1*x0*x1*x2*x20*x23*x9 - L1*b_L*b_m*b_r*d4_1*x0*x2*x20*x23*x5 - a_L*b_m*b_r*d3_1*x0*x1*x20*x23*x6*x7 - b_L*b_m*b_r*d3_1*x0*x1*x20*x23*x6*x7 - L1*a_L*b_L*b_m*b_r*d4_1*x0*x1*x2*x20*x23
    Unc03 = -a_L*b_L*d1_1*x5 - d2_1*x1*x10*x11 - L1*b_L*d4_1*x2*x5 - a_L*d3_1*x1*x6*x7 - b_L*d3_1*x1*x6*x7 - L1*a_L*b_L*d4_1*x1*x2
    Unc04 = Gth(s)*K_act*a_L*d1_1*x12*x14*x17*x20 + Gth(s)*K_act*b_L*d1_1*x12*x14*x17*x20 + Gth(s)*K_act*L1*d4_1*x12*x14*x17*x2*x20 + Gth(s)*K_act*a_L*b_I*d1_1*s*x14*x17*x20*x5 + Gth(s)*K_act*L1*b_I*d4_1*s*x14*x17*x2*x20*x5 + Gth(s)*K_act*a_L*b_L*beta1*d2_1*x12*x14*x17*x20*x8 + Gth(s)*K_act*a_L*b_m*d1_1*s*x14*x17*x20*x5*x9 + Gth(s)*K_act*b_I*d3_1*s*x1*x14*x17*x20*x6*x7 + Gth(s)*K_act*L1*a_L*b_I*d4_1*s*x1*x14*x17*x2*x20 + Gth(s)*K_act*L1*b_m*d4_1*s*x14*x17*x2*x20*x5*x9 + Gth(s)*K_act*b_m*d3_1*s*x1*x14*x17*x20*x6*x7*x9 - Gth(s)*K_act*a_L*b_L*b_m*b_r*d1_1*s*x14*x17*x20*x5 - Gth(s)*K_act*b_m*b_r*d2_1*s*x1*x10*x11*x14*x17*x20 + Gth(s)*K_act*L1*a_L*b_m*d4_1*s*x1*x14*x17*x2*x20*x9 - Gth(s)*K_act*L1*b_L*b_m*b_r*d4_1*s*x14*x17*x2*x20*x5 - Gth(s)*K_act*a_L*b_m*b_r*d3_1*s*x1*x14*x17*x20*x6*x7 - Gth(s)*K_act*b_L*b_m*b_r*d3_1*s*x1*x14*x17*x20*x6*x7 - Gth(s)*K_act*L1*a_L*b_L*b_m*b_r*d4_1*s*x1*x14*x17*x2*x20
    return Unc00, Unc01, Unc02, Unc03, Unc04
