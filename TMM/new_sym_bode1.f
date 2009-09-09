      a_1 = 1/s
      a_2 = 1/(tauact+s)
      a_3 = 1/betabeam
      a_4 = 1/Lbeam
      a_5 = 0.5*abeam*a_3*c4beam*a_4
      a_6 = 1/(cbase*s+kbase)
      a_7 = c1beam*a_6
      a_8 = 1/(cj1*s+kj1)
      a_9 = 1/abeam
      a_10 = 0.5*a_9*betabeam*c2beam*Lbeam*a_6
      a_11 = s**2
      a_12 = Ll0-rl0
      a_13 = Il0*a_11-ml0*a_12*rl0*a_11
      a_14 = a_7+a_5
      a_15 = a_13*a_14
      a_16 = -0.5*betabeam*c2beam*a_4
      a_17 = betabeam**2
      a_18 = 0.5*a_9*a_17*c3beam*a_6
      a_19 = -Ll0*(a_18+a_16)
      a_20 = 1/a_17
      a_21 = -0.5*abeam*a_20*c3beam
      a_22 = 0.5*a_3*c4beam*Lbeam*a_6
      a_23 = a_22+a_21
      a_24 = -ml0*a_12*a_11*a_23
      a_25 = a_8*(a_24+a_19+a_15+a_10+c1beam)
      a_26 = a_25+a_7+a_5
      a_27 = Ll2-rl2
      a_28 = Il2*a_11-ml2*a_27*rl2*a_11
      a_29 = 1/betabeam**3
      a_30 = -0.5*abeam*a_29*c2beam*Lbeam*ml0*a_11
      a_31 = -0.5*abeam*a_29*c2beam*Lbeam
      a_32 = 0.5*abeam*a_20*c3beam*Ll0
      a_33 = a_32+a_31
      a_34 = a_33*ml1*a_11
      a_35 = 0.5*abeam*a_20*c3beam*ml0*rl0*a_11
      a_36 = 0.5*abeam*a_20*c3beam
      a_37 = -0.5*a_3*c4beam*Lbeam
      a_38 = -c1beam*Ll0
      a_39 = 0.5*abeam*a_29*c2beam*Lbeam*ml0*a_12*a_11
      a_40 = 0.5*abeam*a_20*c3beam*a_13
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
      a_63 = a_45*(a_62+a_61+a_58+a_24+a_19+a_15+a_10+c1beam)+a_25+a_7+a
     1   _5
      a_64 = -ml2*rl2*a_11*a_63-ml2*a_11*a_57-ml1*rl1*a_11*a_26-ml1*a_11
     1   *a_56-ml0*a_11*a_23-ml0*rl0*a_11*a_14-0.5*a_9*a_17*c3beam*a_6+0
     2   .5*betabeam*c2beam*a_4
      a_65 = a_28*a_63-Ll2*(ml1*rl1*a_11*a_26+ml1*a_11*a_56+a_60+a_59+a_
     1   18+a_16)-ml2*a_27*a_11*a_57+a_62+a_61+a_58+a_24+a_19+a_15+a_10+
     2   c1beam
      a_66 = 1/(a_53*a_65+(a_28*a_52-Ll2*a_54-ml2*a_27*a_11*a_44+a_51+a_
     1   40+a_49+a_47+a_39+a_38+a_37)*a_64)
      a_67 = -Kact*ml2*rl2*s*(-a_28*a_52+Ll2*a_54+ml2*a_27*a_11*a_44-a_5
     1   0*a_42-0.5*abeam*a_20*c3beam*a_13+Ll1*a_48+a_33*ml1*a_46*a_11-0
     2   .5*abeam*a_29*c2beam*Lbeam*ml0*a_12*a_11+c1beam*Ll0+0.5*a_3*c4b
     3   eam*Lbeam)*a_66*tauact*a_2-Kact*a_1*a_28*a_53*a_66*tauact*a_2
      a_68 = -Kact*ml2*rl2*s*a_65*a_66*tauact*a_2-Kact*a_1*a_28*a_64*a_6
     1   6*tauact*a_2
      RESULT = gainbode1*(a_52*a_68-a_42*a_68+a_63*a_67-a_26*a_67+Kact*a
     1   _1*tauact*a_2)
