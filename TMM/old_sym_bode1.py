from __future__ import division
from scipy import sinh, cosh, sin, cos, real, imag, shape, arange, pi, zeros, array, eye, transpose, dot, conj, c_, poly, vectorize, dot, randn, sum, squeeze, sqrt
#from numpy.lib.scimath import sqrt
from LinearAlgebra import inverse
import pylab
import rwkbode
reload(rwkbode)
from scipy.optimize import newton, fmin, fminbound
pylab.ioff()

import sys, os, copy, time

def old_sym_bode1(s,ucv):
	kbase=ucv[0]
	cbase=ucv[1]
	kj1=ucv[2]
	cj1=ucv[3]
	Kact=ucv[4]
	tauact=ucv[5]
	kj2=ucv[6]
	cj2=ucv[7]
	mubeam=5.7281
	EIbeam=339134.5276
	Lbeam=4.6482
	rl0=0.0902
	Ll0=0.3302
	ml0=16.032
	Il0=0.19
	rl1=0.06145
	Ll1=0.1969
	ml1=5.0264
	Il1=0.027
	rl2=0.1077
	Ll2=0.4001
	ml2=5.5799
	Il2=0.0728
	gainbode1=57.2957795131
	abeam=Lbeam**2/EIbeam
	betabeam=(-1*s**2*Lbeam**4*mubeam/EIbeam)**(0.25)
	c1beam=0.5*cos(betabeam)+0.5*cosh(betabeam)
	c2beam=-sin(betabeam)+sinh(betabeam)
	c3beam=cos(betabeam)-cosh(betabeam)
	c4beam=sin(betabeam)+sinh(betabeam)
	a_1 = 1/s
	a_2 = 1/(tauact+s)
	a_3 = 1/betabeam
	a_4 = 1/Lbeam
	a_5 = abeam*a_3*c4beam*a_4/2.0
	a_6 = 1/(cbase*s+kbase)
	a_7 = c1beam*a_6
	a_8 = 1/(cj1*s+kj1)
	a_9 = 1/abeam
	a_10 = a_9*betabeam*c2beam*Lbeam*a_6/2.0
	a_11 = s**2
	a_12 = Ll0-rl0
	a_13 = Il0*a_11-ml0*a_12*rl0*a_11
	a_14 = a_7+a_5
	a_15 = a_13*a_14
	a_16 = -betabeam*c2beam*a_4/2.0
	a_17 = betabeam**2
	a_18 = a_9*a_17*c3beam*a_6/2.0
	a_19 = -Ll0*(a_18+a_16)
	a_20 = 1/a_17
	a_21 = -abeam*a_20*c3beam/2.0
	a_22 = a_3*c4beam*Lbeam*a_6/2.0
	a_23 = a_22+a_21
	a_24 = -ml0*a_12*a_11*a_23
	a_25 = a_8*(a_24+a_19+a_15+a_10+c1beam)
	a_26 = a_25+a_7+a_5
	a_27 = Ll2-rl2
	a_28 = Il2*a_11-ml2*a_27*rl2*a_11
	a_29 = 1/betabeam**3
	a_30 = -abeam*a_29*c2beam*Lbeam*ml0*a_11/2.0
	a_31 = -abeam*a_29*c2beam*Lbeam/2.0
	a_32 = abeam*a_20*c3beam*Ll0/2.0
	a_33 = a_32+a_31
	a_34 = a_33*ml1*a_11
	a_35 = abeam*a_20*c3beam*ml0*rl0*a_11/2.0
	a_36 = abeam*a_20*c3beam/2.0
	a_37 = -a_3*c4beam*Lbeam/2.0
	a_38 = -c1beam*Ll0
	a_39 = abeam*a_29*c2beam*Lbeam*ml0*a_12*a_11/2.0
	a_40 = abeam*a_20*c3beam*a_13/2.0
	a_41 = a_8*(a_40+a_39+a_38+a_37)
	a_42 = a_41+a_36
	a_43 = ml1*rl1*a_11*a_42
	a_44 = Ll1*a_42+a_32+a_31
	a_45 = 1/(cj2*s+kj2)
	a_46 = Ll1-rl1
	a_47 = -a_33*ml1*a_46*a_11
	a_48 = a_35+a_30+c1beam
	a_49 = -Ll1*a_48
	a_50 = Il1*a_11-ml1*a_46*rl1*a_11
	a_51 = a_50*a_42
	a_52 = a_45*(a_51+a_40+a_49+a_47+a_39+a_38+a_37)+a_41+a_36
	a_53 = ml2*rl2*a_11*a_52+ml2*a_11*a_44+a_43+a_35+a_34+a_30+c1beam
	a_54 = a_43+a_35+a_34+a_30+c1beam
	a_55 = Ll0*a_14
	a_56 = a_55+a_22+a_21
	a_57 = Ll1*a_26+a_55+a_22+a_21
	a_58 = -ml1*a_46*a_11*a_56
	a_59 = ml0*rl0*a_11*a_14
	a_60 = ml0*a_11*a_23
	a_61 = -Ll1*(a_60+a_59+a_18+a_16)
	a_62 = a_50*a_26
	a_63 = a_45*(a_62+a_61+a_58+a_24+a_19+a_15+a_10+c1beam)+a_25+a_7+a_5
	a_64 = -ml2*rl2*a_11*a_63-ml2*a_11*a_57-ml1*rl1*a_11*a_26-ml1*a_11 \
		*a_56-ml0*a_11*a_23-ml0*rl0*a_11*a_14-a_9*a_17*c3beam*a_6/2.0+betabeam \
		*c2beam*a_4/2.0
	a_65 = a_28*a_63-Ll2*(ml1*rl1*a_11*a_26+ml1*a_11*a_56+a_60+a_59+a_18 \
		+a_16)-ml2*a_27*a_11*a_57+a_62+a_61+a_58+a_24+a_19+a_15+a_10+ \
		c1beam
	a_66 = 1/(a_53*a_65+(a_28*a_52-Ll2*a_54-ml2*a_27*a_11*a_44+a_51+a_40 \
		+a_49+a_47+a_39+a_38+a_37)*a_64)
	a_67 = -Kact*ml2*rl2*s*(-a_28*a_52+Ll2*a_54+ml2*a_27*a_11*a_44-a_50 \
		*a_42-abeam*a_20*c3beam*a_13/2.0+Ll1*a_48+a_33*ml1*a_46*a_11-abeam \
		*a_29*c2beam*Lbeam*ml0*a_12*a_11/2.0+c1beam*Ll0+a_3*c4beam* \
		Lbeam/2.0)*a_66*tauact*a_2-Kact*a_1*a_28*a_53*a_66*tauact*a_2
	a_68 = -Kact*ml2*rl2*s*a_65*a_66*tauact*a_2-Kact*a_1*a_28*a_64*a_66 \
		*tauact*a_2
	bode = gainbode1*(a_52*a_68-a_42*a_68+a_63*a_67-a_26*a_67+Kact*a_1 \
		*tauact*a_2)
	return bode
