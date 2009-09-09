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

def old_sym_bode0(s,ucv):
	kbase=ucv[0]
	cbase=ucv[1]
	kj1=ucv[2]
	cj1=ucv[3]
	Kact=ucv[4]
	tauact=ucv[5]
	kj2=ucv[6]
	cj2=ucv[7]
	gainbode0=ucv[8]
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
	abeam=Lbeam**2/EIbeam
	betabeam=(-1*s**2*Lbeam**4*mubeam/EIbeam)**(0.25)
	c1beam=0.5*cos(betabeam)+0.5*cosh(betabeam)
	c2beam=-sin(betabeam)+sinh(betabeam)
	c3beam=cos(betabeam)-cosh(betabeam)
	c4beam=sin(betabeam)+sinh(betabeam)
	a_1 = s**2
	a_2 = betabeam**2
	a_3 = 1/a_2
	a_4 = -abeam*a_3*c3beam/2.0
	a_5 = 1/betabeam
	a_6 = 1/(cbase*s+kbase)
	a_7 = a_5*c4beam*Lbeam*a_6/2.0
	a_8 = a_7+a_4
	a_9 = 1/s
	a_10 = Ll2-rl2
	a_11 = Il2*a_1-ml2*a_10*rl2*a_1
	a_12 = 1/betabeam**3
	a_13 = -abeam*a_12*c2beam*Lbeam*ml0*a_1/2.0
	a_14 = -abeam*a_12*c2beam*Lbeam/2.0
	a_15 = abeam*a_3*c3beam*Ll0/2.0
	a_16 = a_15+a_14
	a_17 = a_16*ml1*a_1
	a_18 = abeam*a_3*c3beam*ml0*rl0*a_1/2.0
	a_19 = abeam*a_3*c3beam/2.0
	a_20 = 1/(cj1*s+kj1)
	a_21 = -a_5*c4beam*Lbeam/2.0
	a_22 = -c1beam*Ll0
	a_23 = Ll0-rl0
	a_24 = abeam*a_12*c2beam*Lbeam*ml0*a_23*a_1/2.0
	a_25 = Il0*a_1-ml0*a_23*rl0*a_1
	a_26 = abeam*a_3*c3beam*a_25/2.0
	a_27 = a_20*(a_26+a_24+a_22+a_21)
	a_28 = a_27+a_19
	a_29 = ml1*rl1*a_1*a_28
	a_30 = Ll1*a_28+a_15+a_14
	a_31 = 1/(cj2*s+kj2)
	a_32 = Ll1-rl1
	a_33 = -a_16*ml1*a_32*a_1
	a_34 = a_18+a_13+c1beam
	a_35 = -Ll1*a_34
	a_36 = Il1*a_1-ml1*a_32*rl1*a_1
	a_37 = a_36*a_28
	a_38 = a_31*(a_37+a_26+a_35+a_33+a_24+a_22+a_21)+a_27+a_19
	a_39 = ml2*rl2*a_1*a_38+ml2*a_1*a_30+a_29+a_18+a_17+a_13+c1beam
	a_40 = a_29+a_18+a_17+a_13+c1beam
	a_41 = 1/Lbeam
	a_42 = 1/abeam
	a_43 = abeam*a_5*c4beam*a_41/2.0
	a_44 = c1beam*a_6
	a_45 = a_44+a_43
	a_46 = Ll0*a_45
	a_47 = a_46+a_7+a_4
	a_48 = a_42*betabeam*c2beam*Lbeam*a_6/2.0
	a_49 = a_25*a_45
	a_50 = -betabeam*c2beam*a_41/2.0
	a_51 = a_42*a_2*c3beam*a_6/2.0
	a_52 = -Ll0*(a_51+a_50)
	a_53 = -ml0*a_23*a_1*a_8
	a_54 = a_20*(a_53+a_52+a_49+a_48+c1beam)
	a_55 = a_54+a_44+a_43
	a_56 = Ll1*a_55+a_46+a_7+a_4
	a_57 = -ml1*a_32*a_1*a_47
	a_58 = ml0*rl0*a_1*a_45
	a_59 = ml0*a_1*a_8
	a_60 = -Ll1*(a_59+a_58+a_51+a_50)
	a_61 = a_36*a_55
	a_62 = a_31*(a_61+a_60+a_57+a_53+a_52+a_49+a_48+c1beam)+a_54+a_44+ \
		a_43
	a_63 = -ml2*rl2*a_1*a_62-ml2*a_1*a_56-ml1*rl1*a_1*a_55-ml1*a_1*a_47 \
		-ml0*a_1*a_8-ml0*rl0*a_1*a_45-a_42*a_2*c3beam*a_6/2.0+betabeam \
		*c2beam*a_41/2.0
	a_64 = a_11*a_62-Ll2*(ml1*rl1*a_1*a_55+ml1*a_1*a_47+a_59+a_58+a_51 \
		+a_50)-ml2*a_10*a_1*a_56+a_61+a_60+a_57+a_53+a_52+a_49+a_48+c1beam
	a_65 = 1/(a_39*a_64+(a_11*a_38-Ll2*a_40-ml2*a_10*a_1*a_30+a_37+a_26 \
		+a_35+a_33+a_24+a_22+a_21)*a_63)
	a_66 = 1/(tauact+s)
	bode = gainbode0*a_1*(a_8*(-Kact*ml2*rl2*s*(-a_11*a_38+Ll2*a_40+ \
		ml2*a_10*a_1*a_30-a_36*a_28-abeam*a_3*c3beam*a_25/2.0+Ll1*a_34+ \
		a_16*ml1*a_32*a_1-abeam*a_12*c2beam*Lbeam*ml0*a_23*a_1/2.0+c1beam \
		*Ll0+a_5*c4beam*Lbeam/2.0)*a_65*tauact*a_66-Kact*a_9*a_11*a_39 \
		*a_65*tauact*a_66)-abeam*a_12*c2beam*Lbeam*(-Kact*ml2*rl2*s*a_64 \
		*a_65*tauact*a_66-Kact*a_9*a_11*a_63*a_65*tauact*a_66)/2.0)
	return bode
