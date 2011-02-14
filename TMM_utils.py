from scipy import *


def build_s_mat(f_min=-1.0, f_max=30.0, f_step=0.1, \
                im_min=-1, im_max=30.0, im_step=0.1, Hz=True):
    """Return a matrix of s values over the grid from f_min to f_max
    and im_min to im_max, with steps of f_step and im_step.

    At each point on the grid

    s[i,j] = -f[j] + 1.0j*im[i]

    If Hz is true, then f and im are assumed to be in Hz and are each
    multiplied by 2*pi.

    The negative sign on f means that the default values focus on the
    second quadrant of the real/imaginary plane."""

    f = arange(f_min, f_max, f_step)
    im = arange(im_min, im_max, im_step)
    
    nr = len(im)
    nc = len(f)
    s = zeros((nr,nc), dtype='D')

    for i in range(nr):
        for j in range(nc):
            s[i,j] = -f[j] + 1.0j*im[i]

    if Hz:
        s = s*2.0*pi

    return s

    
