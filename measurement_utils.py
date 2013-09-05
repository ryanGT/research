from pylab import *
from scipy import *

#from IPython.core.debugger import Pdb
#from RTP_utils import *
#from RTP_utils import dt

def find_overshoot(y, u):
    final_value = u[-1]#the last element in u is assumed to be the
                       #desired stopping point
    max_y = y.max()
    if max_y <= final_value:
        return 0.0
    overshoot = float(max_y - final_value)
    percent_overshoot = overshoot/final_value*100.0
    return percent_overshoot


def find_settling_index(y, u, p=0.01):
    final_value = u[-1]
    e = final_value-y
    ep = abs(p*final_value)
    if ep < 1.0:
        ep = 1.0
    bool1 = abs(e)<= ep
    ## bool0 = u == final_value#use this to ignore the portion before the
    ##                         #step input occurs where u==y==0.
    ## mybool = bool0 & bool1
    mybool = bool1
    #if the final error is not within p*final_value, then the step
    #response never settles
    if not mybool[-1]:
        #print('never settles')
        return None
    #otherwise find the highest index where the error is greater than
    #p*final_value.  The step response settles at the next index.
    not_settled = invert(mybool)
    not_settled_inds = where(not_settled)[0]
    last_not_settled = not_settled_inds.max()
    first_settled_ind = last_not_settled + 1
    return first_settled_ind


def find_settling_time(y, u, t, p=0.01):
    ind = find_settling_index(y,u,p)
    if ind is None:
        ts = None
    else:
        ts = t[ind]
    return ts

def get_t_limits(t):
    t0 = t.min()
    t1 = t.max()
    return array([t0,t1])


def plot_overshoot(u, t, Mp, fignum=1):
    final_value = u[-1]
    max_allowed = final_value*(1.0+Mp/100.0)
    t_plot = get_t_limits(t)
    figure(fignum)
    plot(t_plot, [max_allowed, max_allowed], 'r--')


def plot_settling_lines(u, t, p=0.01, fignum=1, ls='k--', linewidth=1.0):
    final_value = u[-1]
    dy = abs(final_value*p)
    if  dy < 1.0:
        dy = 1.0
    lower_limit = final_value - dy
    upper_limit = final_value + dy
    t_plot = get_t_limits(t)
    figure(fignum)
    plot(t_plot, [lower_limit, lower_limit], ls, \
         label=None, linewidth=linewidth)
    plot(t_plot, [upper_limit, upper_limit], ls, \
         label=None, linewidth=linewidth)


def plot_settling_point(y, u, t, p=0.01, fignum=1, ps='k^', **kwargs):
    ind = find_settling_index(y, u, p)
    if ind is not None:
        figure(fignum)
        plot([t[ind]], [y[ind]], ps, label=None, **kwargs)


def find_steady_state_error(y, u):
    fv = u[-1]
    y_f = y[-1]
    e_ss = 100*float(fv-y_f)/float(fv)
    return e_ss

def main(y, u, t, Mp=10.0, p=0.01, fignum=1):
    ts = find_settling_time(y, u, t, p)
    overshoot = find_overshoot(y, u)
    plot_overshoot(u, t, Mp, fignum=fignum)
    plot_settling_lines(u, t, p, fignum=fignum)
    plot_settling_point(y, u, t, p, fignum=fignum)
    return ts, overshoot


def find_and_plot_settling_time(y, u, t, p=0.01, fignum=1, \
                                ls='--', ps='k^', linewidth=1.0, \
                                markersize=10):
    ts = find_settling_time(y, u, t, p=p)
    print('settling time = %s' % ts)
    plot_settling_lines(u, t, p, fignum=fignum, ls=ls, \
                        linewidth=linewidth)
    plot_settling_point(y, u, t, p, fignum=fignum, ps=ps, \
                        markersize=markersize)
    return ts


def find_and_plot_settling_point(y, u, t, p=0.01, fignum=1, ps='k^'):
    ts = find_settling_time(y, u, t, p=p)
    print('settling time = %s' % ts)
    plot_settling_point(y, u, t, p, fignum=fignum, ps=ps)


def plot_overshoot_and_settling(y, u, t, Mp=10.0, p=0.01, fignum=1):
    plot_overshoot(u, t, Mp, fignum=fignum)
    plot_settling_lines(u, t, p, fignum=fignum)
    plot_settling_point(y, u, t, p, fignum=fignum)


def find_e_ss_penalty(y, u, t, p=0.01):
    ts = find_settling_time(y, u, t, p)
    e_ss = find_steady_state_error(y, u)
    if ts is None:
        #print('did not settle')
        penalty = abs(float(e_ss))
    else:
        #print('ts = %0.4f (before penalty)' % ts)
        penalty = 0.0
    return penalty


def find_overshoot_penalty(y, u, t, Mp=10.0):
    overshoot = find_overshoot(y, u)
    if overshoot > Mp:
        penalty = overshoot-Mp
    else:
        penalty = 0.0
    return penalty


def ts_with_penalties(y, u, t, Mp=10.0, p=0.01, fignum=1):
    #print('steady state error penalty = %0.4f' % penalty1)
    ts = find_settling_time(y, u, t, p)
    e_ss_penalty = find_e_ss_penalty(y, u, t, p=p)
    overshoot = find_overshoot(y, u)
    Mp_penalty = find_overshoot_penalty(y, u, t, Mp=Mp)

    #print('overshoot = %0.2f percent' % overshoot)

    #print('overshoot penalty = %0.4f' % penalty2)
    penalty = e_ss_penalty + Mp_penalty
    if ts is None:
        ts = penalty
    else:
        ts += penalty

    #print('total penalty = %0.4f' % penalty)
    #print('final ts = %0.4f' % ts)
    return ts, overshoot, penalty


def print_results(y, u, t, Mp=10.0, p=0.01):
    ts_raw = find_settling_time(y, u, t, p)
    if ts_raw is None:
        print('did not settle')
    else:
        print('ts = %0.4f (before penalty)' % ts_raw)
    e_ss_penalty = find_e_ss_penalty(y, u, t, p=p)
    print('steady-state eror penalty = %0.4f' % e_ss_penalty)
    ts, overshoot, penalty = ts_with_penalties(y, u, t, Mp=Mp, p=p)
    print('overshoot = %0.2f' % overshoot)
    Mp_penalty = find_overshoot_penalty(y, u, t, Mp=Mp)
    print('overshoot penalty = %0.4f' % Mp_penalty)
    print('total penalty = %0.4f' % penalty)
    print('final ts = %0.4f' % ts)


def plot_results(u, y, v, n, fignum=1, clear=True, plotu=True, \
                 legend_subscript='exp', dt=1.0/500, plot_bars=True):
    t = dt*n#dt gets imported from RTP_utils.py
    figure(fignum)
    if clear:
        clf()
    if plotu:
        plot(t, u, label='$u$')

    plot(t, y, label='$y_{%s}$' % legend_subscript)
    plot(t, v, label='$v_{%s}$' % legend_subscript)

    if plot_bars:
        plot_overshoot_and_settling(y, u, t, fignum=fignum)

    xlabel('Time (sec)')
    ylabel('Signal Amp. (counts)')
    legend(loc=5)



## class step_response_measurer(object):
##     def __init__(self, u, y, t, fignum=1, label='Experimental'):
##         self.u = u
##         self.y = y
##         self.t = t
##         self.fignum = fignum
##         self.label = label


##     def find_settling_time(self):
##         self.ts = find_settling_time(self.y, self.u, self.t, self.p)
##         return self.ts


##     def find_overshoot(self):
##         self.percent_overshoot = find_overshoot(self.y, self.u)


##     def find_settling_index(self):
##         self.index_settled = find_settling_index(self.y, \
##                                                  self.u, \
##                                                  self.p)
##         return self.index_settled


