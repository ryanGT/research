import txt_data_processing as TDP
reload(TDP)

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


OL_bode_opts = [bode_opt_th_v, bode_opt_a_v]
