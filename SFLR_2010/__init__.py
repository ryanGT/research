from pylab import *
from scipy import *

import controls

#import step_response_utils, ROM_model, accel_FB_system

#'theta_fb_sympy_bodes', #<-- this module gives bad results for s with positive real part, so I removed it

modlist = ['step_response_utils', \
           'ROM_model', \
           'accel_fb_system',
           'theta_fb_contour', \
           'theta_fb_maxima_sympy_accel_u_bode', \
           'theta_fb_maxima_sympy_accel_theta_bode', \
           'bode_plots', \
           'model_w_bm_OL_bodes', \
           'OL_contour'
           ]


for modname in modlist:
    __import__('SFLR_2010.' + modname)


def DeepReload():
    for modname in modlist:
        reload(sys.modules['SFLR_2010.' + modname])
    ## reload(step_response_utils)
    ## reload(ROM_model)
    ## reload(accel_FB_system)

