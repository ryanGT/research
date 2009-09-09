      subroutine bodevect(svect, ucv, outvect, n0, n1)
Cf2py integer intent(hide),depend(svect) :: n0 = len(svect)
Cf2py integer intent(hide),depend(ucv) :: n1 = len(ucv)
Cf2py intent(in) svect, ucv
Cf2py intent(out) outvect
      integer i, n0, n1
      double precision ucv(n1)
      double complex svect(n0), outvect(n0)
      double complex bode
      DO i=1,n0
         outvect(i) = bode(svect(i), ucv, n1)
      ENDDO
      END

      subroutine invbodevect(svect, ucv, outvect, n0, n1)
Cf2py integer intent(hide),depend(svect) :: n0 = len(svect)
Cf2py integer intent(hide),depend(ucv) :: n1 = len(ucv)
Cf2py intent(in) svect, ucv
Cf2py intent(out) outvect
      integer i, n0, n1
      double precision ucv(n1)
      double complex svect(n0), outvect(n0)
      double complex invbode
      DO i=1,n0
         outvect(i) = invbode(svect(i), ucv, n1)
      ENDDO
      END

      double complex function invbode(s, ucv, n1)
Cf2py integer intent(hide),depend(ucv) :: n1 = len(ucv)
Cf2py intent(in) s, ucv
      integer n1
      double precision ucv(n1)
      double complex s
Cf2py intent(out) invbode
      double complex bode
      invbode = 1.0/bode(s, ucv, n1)
      RETURN
      END

      double complex function zcosh(z)
      double complex z
      zcosh = 0.5*(exp(z)+exp(-z))
      RETURN
      END
      
      double complex function zsinh(z)
      double complex z
      zsinh = 0.5*(exp(z)-exp(-z))
      RETURN
      END

      double complex function poly(s,coeffs,n)
Cf2py intent(in) s, coeffs
Cf2py intent(out) poly
Cf2py integer intent(hide),depend(coeffs) :: n = len(coeffs)
      double complex s
      integer n
      double precision coeffs(n)
      poly=0
      DO i=1,n
!         PRINT *,'i=', i
         poly = poly+coeffs(i)*s**(n-i)
      ENDDO
      RETURN
      END

      double complex function bode(s, ucv, n1)
Cf2py integer intent(hide),depend(ucv) :: n1 = len(ucv)
Cf2py intent(in) s, ucv
      integer n1
      double precision ucv(n1)
      double complex s
Cf2py intent(out) bode
      double complex zsinh, zcosh
      double precision kbase, cbase, mubeam, EIbeam, Lbeam, rl0, Ll0,
     @    ml0, Il0, kj1, cj1, rl1, Ll1, ml1, Il1, Kact, tauact, kj2,
     @    cj2, rl2, Ll2, ml2, Il2, gainbode0, gainbode1, abeam
      double complex c1beam, c2beam, c3beam, c4beam, betabeam, a_1, a_2,
     @    a_3, a_4, a_5, a_6, a_7, a_8, a_9, a_10, a_11, a_12, a_13,
     @    a_14, a_15, a_16, a_17, a_18, a_19, a_20, a_21, a_22, a_23,
     @    a_24, a_25, a_26, a_27, a_28, a_29, a_30, a_31, a_32, a_33,
     @    a_34, a_35, a_36, a_37, a_38, a_39, a_40, a_41, a_42, a_43,
     @    a_44, a_45, a_46, a_47, a_48, a_49, a_50, a_51, a_52, a_53,
     @    a_54, a_55, a_56, a_57, a_58, a_59, a_60, a_61, a_62, a_63,
     @    a_64, a_65, a_66
      kbase=ucv(1)
      cbase=ucv(2)
      kj1=ucv(3)
      cj1=ucv(4)
      Kact=ucv(5)
      tauact=ucv(6)
      kj2=ucv(7)
      cj2=ucv(8)
      gainbode0=ucv(9)
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
      c1beam=0.5*zcos(betabeam)+0.5*zcosh(betabeam)
      c2beam=-zsin(betabeam)+zsinh(betabeam)
      c3beam=zcos(betabeam)-zcosh(betabeam)
      c4beam=zsin(betabeam)+zsinh(betabeam)
      a_1 = s**2
      a_2 = betabeam**2
      a_3 = 1/a_2
      a_4 = -0.5*abeam*a_3*c3beam
      a_5 = 1/betabeam
      a_6 = 1/(cbase*s+kbase)
      a_7 = 0.5*a_5*c4beam*Lbeam*a_6
      a_8 = a_7+a_4
      a_9 = 1/s
      a_10 = Ll2-rl2
      a_11 = Il2*a_1-ml2*a_10*rl2*a_1
      a_12 = 1/betabeam**3
      a_13 = -0.5*abeam*a_12*c2beam*Lbeam*ml0*a_1
      a_14 = -0.5*abeam*a_12*c2beam*Lbeam
      a_15 = 0.5*abeam*a_3*c3beam*Ll0
      a_16 = a_15+a_14
      a_17 = a_16*ml1*a_1
      a_18 = 0.5*abeam*a_3*c3beam*ml0*rl0*a_1
      a_19 = 0.5*abeam*a_3*c3beam
      a_20 = 1/(cj1*s+kj1)
      a_21 = -0.5*a_5*c4beam*Lbeam
      a_22 = -c1beam*Ll0
      a_23 = Ll0-rl0
      a_24 = 0.5*abeam*a_12*c2beam*Lbeam*ml0*a_23*a_1
      a_25 = Il0*a_1-ml0*a_23*rl0*a_1
      a_26 = 0.5*abeam*a_3*c3beam*a_25
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
      a_43 = 0.5*abeam*a_5*c4beam*a_41
      a_44 = c1beam*a_6
      a_45 = a_44+a_43
      a_46 = Ll0*a_45
      a_47 = a_46+a_7+a_4
      a_48 = 0.5*a_42*betabeam*c2beam*Lbeam*a_6
      a_49 = a_25*a_45
      a_50 = -0.5*betabeam*c2beam*a_41
      a_51 = 0.5*a_42*a_2*c3beam*a_6
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
      a_62 = a_31*(a_61+a_60+a_57+a_53+a_52+a_49+a_48+c1beam)+a_54+a_44+
     1   a_43
      a_63 = -ml2*rl2*a_1*a_62-ml2*a_1*a_56-ml1*rl1*a_1*a_55-ml1*a_1*a_4
     1   7-ml0*a_1*a_8-ml0*rl0*a_1*a_45-0.5*a_42*a_2*c3beam*a_6+0.5*beta
     2   beam*c2beam*a_41
      a_64 = a_11*a_62-Ll2*(ml1*rl1*a_1*a_55+ml1*a_1*a_47+a_59+a_58+a_51
     1   +a_50)-ml2*a_10*a_1*a_56+a_61+a_60+a_57+a_53+a_52+a_49+a_48+c1b
     2   eam
      a_65 = 1/(a_39*a_64+(a_11*a_38-Ll2*a_40-ml2*a_10*a_1*a_30+a_37+a_2
     1   6+a_35+a_33+a_24+a_22+a_21)*a_63)
      a_66 = 1/(tauact+s)
      bode = gainbode0*a_1*(a_8*(-Kact*ml2*rl2*s*(-a_11*a_38+Ll2*a_40+
     1   ml2*a_10*a_1*a_30-a_36*a_28-0.5*abeam*a_3*c3beam*a_25+Ll1*a_34+
     2   a_16*ml1*a_32*a_1-0.5*abeam*a_12*c2beam*Lbeam*ml0*a_23*a_1+c1be
     3   am*Ll0+0.5*a_5*c4beam*Lbeam)*a_65*tauact*a_66-Kact*a_9*a_11*a_3
     4   9*a_65*tauact*a_66)-0.5*abeam*a_12*c2beam*Lbeam*(-Kact*ml2*rl2*
     5   s*a_64*a_65*tauact*a_66-Kact*a_9*a_11*a_63*a_65*tauact*a_66))
      RETURN
      END
