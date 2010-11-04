import txt_data_processing as TDP
#reload(TDP)

import copy

#Specify Bode Plotting Options
myfreqlim = [0.5,30]
#myfreqlim = [0.1,100]

bode_opt_th_v = TDP.Bode_Options(input_label='v', \
                                 output_label='theta', \
                                 freqlim=myfreqlim, \
                                 maglim=[-50,20], \
                                 phaselim=[-250,25], \
                                 seedfreq=1.0, \
                                 seedphase=-120.0)

bode_opt_th_u = TDP.Bode_Options(input_label='u', \
                                 output_label='theta', \
                                 freqlim=myfreqlim, \
                                 maglim=[-50,20], \
                                 phaselim=[-250,25], \
                                 seedfreq=1.0, \
                                 seedphase=-30.0)

bode_opt_a_v = TDP.Bode_Options(input_label='v', \
                                output_label='a', \
                                freqlim=myfreqlim, \
                                maglim=[-20, 20], \
                                phaselim=[-400,200], \
                                seedfreq=1.0, \
                                seedphase=90.0)

bode_opt_a_u = TDP.Bode_Options(input_label='u', \
                                output_label='a', \
                                freqlim=myfreqlim, \
                                maglim=[-20, 35], \
                                phaselim=[-400,200], \
                                seedfreq=1.0, \
                                seedphase=150.0)


bode_opt_a_theta = TDP.Bode_Options(input_label='theta', \
                                    output_label='a', \
                                    freqlim=myfreqlim, \
                                    maglim=[-20, 45], \
                                    phaselim=[-250,250], \
                                    seedfreq=1.0, \
                                    seedphase=170.0)


OL_Bode_opts = [bode_opt_th_v, bode_opt_a_v]
OL_bode_opts = OL_Bode_opts

ThetaFB_Bode_opts = [bode_opt_th_u, \
                     bode_opt_a_u, \
                    ]


AccelFB_Bode_opts = copy.copy(ThetaFB_Bode_opts)
AccelFB_Bode_opts[0].phaselim = [-200,100]
AccelFB_Bode_opts[0].seedfreq = 5.0
AccelFB_Bode_opts[0].seedphase = 0.0

AccelFB_Bode_opts[1].seedfreq = 5.0
AccelFB_Bode_opts[1].seedphase = 0.0

Accel_Comp_Bode_opts = copy.copy(AccelFB_Bode_opts)
Accel_Comp_Bode_opts[1].seedfreq = 10.0
Accel_Comp_Bode_opts[1].seedphase = -200.0
