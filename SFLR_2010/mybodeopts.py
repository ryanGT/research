import SFLR_bode_options

OL_Bode_opts = [SFLR_bode_options.bode_opt_th_v, \
                SFLR_bode_options.bode_opt_a_v, \
                ]

ThetaFB_Bode_opts = [SFLR_bode_options.bode_opt_th_u, \
                     SFLR_bode_options.bode_opt_a_u, \
                    ]

import copy

AccelFB_Bode_opts = copy.copy(ThetaFB_Bode_opts)
AccelFB_Bode_opts[0].phaselim = [-200,100]
AccelFB_Bode_opts[0].seedfreq = 5.0
AccelFB_Bode_opts[0].seedphase = 0.0

AccelFB_Bode_opts[1].seedfreq = 5.0
AccelFB_Bode_opts[1].seedphase = 0.0
