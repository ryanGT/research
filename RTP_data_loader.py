#from pylab import *
from scipy import *

#import controls
import os
import pylab_util
import measurement_utils

class data_loader(object):
    def load_data(self):
        data = loadtxt(self.pathin, skiprows=self.skiprows, \
                       delimiter=self.delimiter)
        self.data = data
        self.t = data[:,0]
        self.n = data[:,1]
        self.u = data[:,2]
        self.v = data[:,3]
        self.y = data[:,4]

        
    def __init__(self, pathin, use_t=True, skiprows=2, delimiter='\t'):
        self.pathin = pathin
        self.use_t = use_t
        self.skiprows = skiprows
        self.delimiter = delimiter

        if os.path.exists(self.pathin):
            self.load_data()

    def _find_filename(self, pathout=None, ext='.eps'):
        if pathout is None:
            pno, old_ext = os.path.splitext(self.pathin)
            if ext[0] != '.':
                ext = '.' + ext
            pathout = pno + ext
            self.pathout = pathout
        return pathout

    def Plot(self, fignum=1, use_t=None, plotu=True, clear=True, \
             legend=None, legloc=None, xlabel=None, ylabel=None, \
             ylim=None, **kwargs):
        self.fignum = fignum
        if use_t is None:
            use_t = self.use_t
        if use_t:
            x = self.t
            if xlabel is None:
                xlabel = 'Time (sec)'
        else:
            x = self.n
            if xlabel is None:
                xlabel = 'ISR Counts'
        if ylabel is None:
            ylabel = 'Signal Amplitude (counts)'
        if plotu:
            plot_data = column_stack([self.u, self.v, self.y])
        else:
            plot_data = column_stack([self.v, self.y])

        pylab_util.plot_cols(x, plot_data, fi=fignum, \
                             leg=legend, clear=clear, ylim=ylim, \
                             legloc=legloc, \
                             xlabel=xlabel, **kwargs)
        pylab_util.set_ylabel(ylabel, fi=fignum)


    def Save(self, pathout=None, fignum=None, ext='.eps'):
        if fignum is None:
            fignum = self.fignum
        if pathout is None:
            pathout = self._find_filename(ext=ext)
        else:
            self.pathout = pathout
        pylab_util.mysave(pathout, fignum)
        

    def plot_overshoot_and_settling(self, Mp=10.0, p=0.01, fignum=None):
        if fignum is None:
            fignum = self.fignum
        measurement_utils.plot_overshoot_and_settling(self.y, self.u, self.t, \
                                                      Mp=Mp, p=p, \
                                                      fignum=fignum)
