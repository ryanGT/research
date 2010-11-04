from pylab import *
from scipy import *

import controls

#import step_response_utils, ROM_model, accel_FB_system

modlist = ['step_response_utils', \
           'ROM_model', \
           'accel_FB_system',
           'theta_fb_contour', \
           'theta_fb_sympy_bodes', \
           ]


for modname in modlist:
    __import__('SFLR_2010.' + modname)


def DeepReload():
    for modname in modlist:
        reload(sys.modules['SFLR_2010.' + modname])
    ## reload(step_response_utils)
    ## reload(ROM_model)
    ## reload(accel_FB_system)

