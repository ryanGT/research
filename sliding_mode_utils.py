from scipy import *

from IPython.Debugger import Pdb

def my_smooth_s(delta_x, max_v, max_a, cushion=2.0, dt=0.002):
    max_v = float(max_v)
    max_a = float(max_a)
    delta_x = float(delta_x)
    t1 = max_v/max_a
    n1 = int(t1/dt)
    t1 = n1*dt
    x1 = 0.5*max_a*t1**2

    if delta_x <= 2*x1:
        #the system won't get up to max_v
        case = 1
        x1 = delta_x/2.0
        t1 = sqrt(2*x1/max_a)
        n1 = int(t1/dt)
        t1 = dt*n1
        t2 = 2*t1
        v1 = max_a*t1
        x1 = 0.5*max_a*t1**2#recalc with actual t1
        dxa = 2*x1#actual displacment
        round_e = delta_x-dxa#displacement error do to rounding t1 and
                             #t2 to multiples of dt
    else:
        case = 2
        x12 = delta_x - 2*x1#this is the distance covered at max_v
        t12 = t1 + x12/max_v
        n12 = int(t12/dt)
        t12 = n12*dt
        x12 = t12*max_v
        t2 = t1+t12
        v1 = max_v
        dxa = 2*x1+x12
        round_e = delta_x-dxa
        
    T = t2+cushion

    t = arange(0,T,dt)
    #n1 = t1/dt
    n2 = t2/dt

    n1 += 1#trying to handle the python index conventions
    xd = zeros_like(t)
    xd_dot = zeros_like(t)
    xd_ddot = zeros_like(t)
    #part 1 - ramping up velocity
    xd_ddot[0:n1] = max_a
    xd_dot[0:n1] = max_a*t[0:n1]
    xd[0:n1] = 0.5*max_a*(t[0:n1])**2
    if case == 2:
        #part 2 - constant velocity
        n12 += 1
        n2 += 1
        xd_ddot[n1:n12] = 0
        xd_dot[n1:n12] = max_v
        xd[n1:n12] = xd[n1-1]+max_v*(t[n1:n12]-t[n1-1])
    else:
        n12 = n1
        t12 = t1
        n2 += 2
    #part 3 - slowing down
    x12 = xd[n12-1]
    dx = delta_x - x12
    dt2 = t2 - t12
    dxa = v1*dt-0.5*max_a*dt2**2
    tol = 1e-6
    if dx - dxa > tol:
        a2 = (v1*dt2-dx)*2.0/(dt2**2)
    else:
        a2 = max_a
    xd_ddot[n12:n2] = -a2
    xd_dot[n12:n2] = v1 - a2*(t[n12:n2]-t[n12-1])
    xd[n12:n2] = xd[n12-1] + v1*(t[n12:n2]-t[n12-1]) - \
                 0.5*a2*(t[n12:n2]-t[n12-1])**2
    #part 4 - done
    xd[n2:] = delta_x
    xd_dot[n2:] = 0.0
    xd_ddot[n2:] = 0.0

    return t, xd, xd_dot, xd_ddot


def sat(x):
    if abs(x) <= 1:
        return x
    elif x > 1:
        return 1.0
    elif x < -1:
        return -1.0
