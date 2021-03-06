from __future__ import division
from scipy import *

def Bodes(s, params):
    if type(params) == dict:
        params = SLFR_TMM.SLFR_params(**params)
    EI = params.EI
    L1 = params.L
    L2 = params.L2
    mu = params.mu
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
    k_clamp = params.k_clamp
    c_clamp = params.c_clamp
    K_act = params.K_act
    tau = params.tau
    a_gain = params.a_gain
    p_act1 = params.p_act1
    p_act2 = params.p_act2
    z_act = params.z_act
    I = 1.0j
    #--------------------------------
    x0 = 1/s
    x1 = s + tau
    x2 = 1/x1
    x3 = x0**(-2)
    x4 = a_L - a_r
    x5 = 1/beta1
    x6 = L1**(-2)
    x7 = x5**(-2)
    x8 = L2**(-2)
    x9 = beta2**2
    x10 = 1/EI
    x11 = 1/x6
    x12 = 1/x7
    x13 = x11**(3/2)
    x14 = x12**(3/2)
    x15 = x8**(1/2)
    x16 = a_I*x3
    x17 = -a_m*a_r*x3*x4
    x18 = x16 + x17
    x19 = c_clamp*s
    x20 = k_clamp + x19
    x21 = 1/x20
    x22 = x13**(-1/3)
    x23 = d1_1*x21
    x24 = L1*d4_1*x10*x5
    x25 = x23 + x24
    x26 = d3_1*x10*x11*x12
    x27 = L1*d4_1*x21*x5
    x28 = x9**(-1/2)
    x29 = EI*d3_1*x21*x6*x7
    x30 = beta1*d2_1*x22
    x31 = x29 + x30
    x32 = -x31
    x33 = x26 + x27
    x34 = x15**3
    x35 = x28**(-3)
    x36 = d2_1*x10*x13*x14
    x37 = a_L*x26
    x38 = x36 + x37
    x39 = -x38
    x40 = -1.0*x18*x26
    x41 = -a_L*d1_1
    x42 = a_m*x3*x36*x4
    x43 = -1.0*L1*d4_1*x5
    x44 = x40 + x41 + x42 + x43
    x45 = a_m*a_r*x26*x3
    x46 = a_m*x3*x36
    x47 = x45 + x46
    x48 = d1_1 - x47
    x49 = EI*x21*x30
    x50 = x18*x25
    x51 = d1_1 + x49 + x50
    x52 = a_m*x3*x33*x4
    x53 = a_L*x32
    x54 = x52 + x53
    x55 = x51 - x54
    x56 = a_m*a_r*x25*x3
    x57 = a_m*x3*x33
    x58 = x32 + x56 + x57
    x59 = a_L*x25
    x60 = x33 + x59
    x61 = EI*d3_2*x60*x8*x9
    x62 = d1_2*x55
    x63 = EI*beta2*d2_2*x15*x25
    x64 = x61 + x62 + x63
    x65 = L2*d4_2*x28*x58
    x66 = x64 - x65
    x67 = 1/x66
    x68 = K_act*a_m*a_r*d1_1*s*tau*x2
    x69 = K_act*L1*a_m*d4_1*s*tau*x2*x5
    x70 = x68 + x69
    x71 = EI*K_act*d3_1*tau*x0*x2*x6*x7
    x72 = x70 - x71
    x73 = K_act*L1*d4_1*tau*x0*x2*x5
    x74 = K_act*a_L*d1_1*tau*x0*x2
    x75 = x73 + x74
    x76 = EI*K_act*tau*x0*x2*x30
    x77 = K_act*d1_1*tau*x0*x18*x2
    x78 = K_act*a_m*s*tau*x2*x4*x43
    x79 = a_L*x71
    x80 = x76 + x77 + x78 + x79
    x81 = EI*d3_2*x39*x8*x9
    x82 = -L2*d4_2*x28*x48
    x83 = -1.0*beta2*d2_2*d3_1*x11*x12*x15
    x84 = d1_2*x44
    x85 = x81 + x82 + x83 + x84
    x86 = -beta2*d2_2*x15*x55
    x87 = d1_2*x58
    x88 = -EI*beta2*d4_2*x34*x60*x9
    x89 = -EI*d3_2*x25*x8*x9
    x90 = x86 + x87 + x88 + x89
    x91 = d3_1*d3_2*x11*x12*x8*x9
    x92 = -EI*d4_2*x34*x35*x39
    x93 = -beta2*d2_2*x15*x44
    x94 = L2*d4_2*x28*x48
    x95 = beta2*d2_2*d3_1*x11*x12*x15
    x96 = -x81
    x97 = -x84
    x98 = x94 + x95 + x96 + x97
    x99 = x67*x90*x98
    x100 = d1_2*x48
    x101 = x100 + x91 + x92 + x93 + x99
    x102 = 1/x101
    x103 = -beta2*d2_2*x15*x80
    x104 = d1_2*x72
    x105 = -EI*beta2*d4_2*x34*x75*x9
    x106 = -EI*K_act*d1_1*d3_2*tau*x0*x2*x8*x9
    x107 = x103 + x104 + x105 + x106
    x108 = x102*x107*x67*x85
    x109 = x67**2
    x110 = x102*x109*x85*x90
    x111 = x110 + x67
    x112 = K_act*tau*x0*x2
    x113 = -1.0*d1_2*x80
    x114 = -1.0*EI*d3_2*x75*x8*x9
    x115 = -1.0*EI*beta2*d1_1*d2_2*x112*x15
    x116 = L2*d4_2*x28*x72
    x117 = x113 + x114 + x115 + x116
    x118 = x111*x117
    x119 = x108 + x118
    x120 = -x116
    x121 = EI*d3_2*x75*x8*x9
    x122 = EI*beta2*d1_1*d2_2*x112*x15
    x123 = d1_2*x80
    x124 = x120 + x121 + x122 + x123
    x125 = x119*x21
    x126 = x112 + x125
    x127 = x102*x124*x67*x90
    x128 = x102*x107
    x129 = x127 - x128
    th_out = 512.0*(x112 + x21*(x108 - x111*x124))/pi
    a_out = a_gain*x3*(a_L*(d1_1*x126 + L1*d4_1*x10*x119*x5 - d3_1*x10*x11*x12*x129) + L1*d4_1*x126*x5 + d3_1*x10*x11*x119*x12 - d2_1*x10*x129*x13*x14)
    return th_out, a_out
