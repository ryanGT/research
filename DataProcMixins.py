from scipy import *
from scipy import signal
from pylab import figure, cla, clf, plot, subplot, show, ylabel, xlabel, xlim, ylim, semilogx, legend, title, savefig, yticks, grid, rcParams

#from IPython.core.debugger import Pdb

import copy, os, sys, time

mt=time.time

from rwkdataproc import thresh, CalcSpectra, makefreqvect, edges
from rwkmisc import reverse

import rwkbode


def _staythresh(vectin, thresh, staypoints=20, above=True):
    if above:
        inds = where(vectin>=thresh)[0]
    else:
        inds = where(vectin<=thresh)[0]
    for n in range(staypoints):
        keepus = inds[0:-1]+1==inds[1:]#test if the current point and the next point both satisfy thresh
        inds = inds[keepus]
    return inds.min()


def staythreshfwd(vectin, thresh, staypoints=20, above=True, startind=0):
    myvect = vectin[startind:]
    ind1 = _staythresh(vectin, thresh, staypoints=staypoints, above=above)
    return ind1+startind


def staythreshbackwd(vectin, thresh, staypoints=20, above=True, startind=None):
    if startind is None:
        startind = vectin.shape[0]
    myvect = vectin[0:startind]
    myvect = reverse(myvect)
    ind1 = _staythresh(myvect, thresh, staypoints=staypoints, above=above)
    return startind-ind1



class FilterChannel:
    def __init__(self, signal, cutoff_freq, sample_freq, order=2, tvect=None):
        self.signal = signal
        self.fc = cutoff_freq
        self.fs = sample_freq
        self.order = order
        self.wn = self.fc/(0.5*self.fs)
        self.Filter()
        if tvect is not None:
            self.t = tvect
        else:
            dt = 1.0/self.fs
            self.t = dt*array(range(len(self.signal)))

    def Filter(self, cutoff_freq=None):
        if cutoff_freq is not None:
            self.fc = cutoff_freq
            self.wn = self.fc/(0.5*self.fs)
        [b,a] = signal.butter(self.order, self.wn)
        self.b = b
        self.a = a
        self.filtered = signal.lfilter(self.b, self.a, self.signal)

    def Plot(self, fig, clear=True):
        if clear:
            fig.clf()
        ax = fig.add_subplot(1,1,1)
        ax.plot(self.t, self.signal, self.t, self.filtered)
        return ax



class Trigger:
    """This class is used to create various triggers off of digital or
    analog signals."""
    def __init__(self, signal, timevector, threshlevel=None):
        self.signal = signal
        self.t = timevector
        self.threshlevel = threshlevel


    def FindTrig(self, threshlevel=None, startind=0, above=True):
        #thresh(iterin, value, startind=0, above=1)
        if threshlevel is None:
            threshlevel = self.threshlevel
        ind = thresh(self.signal, threshlevel, startind=startind, above=above)
        return ind


class DropDownLightGateTrigger(Trigger):
    def __init__(self, signal, timevector, threshlevel=1.0):
        if threshlevel is None:
            threshlevel = 1.0
        Trigger.__init__(self, signal, timevector, threshlevel=threshlevel)


    def FindTrig(self, threshlevel=None, startind=0):
        #firstedge = Trigger.FindTrig(self, startind=startind, above=False)
        #secondedge = Trigger.FindTrig(self, startind=firstedge, above=True)
        risingedges, fallingedges = edges(self.signal, edges='both', startind=startind)
        firstedge = fallingedges[0]
        self.firstedge = firstedge
        secondedge = risingedges[0]
        self.secondedge = secondedge
        self.t0 = self.t[firstedge]
        if secondedge==-1:
            self.t1 = self.t.max()
        else:
            self.t1 = self.t[secondedge]
        self.dt = self.t1-self.t0
        return secondedge


class ForwardBackwardTrigger(Trigger):
    """This class implements a trigger I found very useful in dealing
    with FMH impact data.  It first searches forward until a vector
    stays above a certain value and then searches backward until a
    vector stays below a different value."""
    def __init__(self, signal, forwardthresh, backwardthresh, staypoints=20):
        """Search through signal until it stays above forwardthresh
        for staypoints nummber of points.  Then search backward until
        vector stays below bacwardthresh for staypoints."""
        self.signal = signal
        self.forwardthresh = forwardthresh
        self.backwardthresh = backwardthresh
        self.staypoints = staypoints

    def FindTrig(self):
        ind1 = staythreshfwd(self.signal, self.forwardthresh, staypoints=self.staypoints)
        ind2 = staythreshbackwd(self.signal, self.backwardthresh, staypoints=self.staypoints, above=False, startind=ind1)
        return ind2


class StayTrigger(Trigger):
    def __init__(self, signal, thresh, staypoints=20):
        """Search through signal until it stays above or below thresh
        for staypoints nummber of points."""
        self.signal = signal
        self.thresh = thresh
        self.staypoints = staypoints

    def FindTrig(self, startind=0, above=True):
        ind1 = staythreshfwd(self.signal, self.thresh, above=above, staypoints=self.staypoints, startind=startind)
        return ind1


class AnalogTruncateChannel:
    def __init__(self, signal, timevector):
        self.signal = signal
        self.t = timevector

    def Trunc(self, trigger, threshlevel=None, duration=0.05, backup=0.01):
        """Truncate channel based on a trigger time that comes from
        trigger.FindTrig(threshlevel).  trigger is expected to be an
        instance of Trigger (or something that derives from it)."""
        ind1 = trigger.FindTrig(threshlevel)
        t0 = self.t[ind1]
        self.t0 = t0
        ind0 = thresh(self.t, t0-backup)
        try:
            ind2 = thresh(self.t, self.t[ind0]+duration)
        except:
            ind2 = len(self.t)

        self.sigtrunc = copy.copy(self.signal[ind0:ind2])
        self.ttrunc = copy.copy(self.t[ind0:ind2])-t0
        return ind0, ind2





class FFTChannel:
    def __init__(self, signal, time_vector, seedfreq=None, seedphase=None):
        self.signal = signal
        self.time = time_vector
        self.f = makefreqvect(self.time)
        N = max(shape(signal))
        self.comp = squeeze(fft(signal, None, 0)*2/N)
        self.mag = abs(self.comp)
        self.phase = arctan2(imag(self.comp),real(self.comp))*180.0/pi
        self.dBmag = 20.0*log10(self.mag)


    def Plot(self, fig, style='Bode', clear=True, label=None):
        if clear:
            fig.clf()
        ax1 = fig.add_subplot(2,1,1)
        if clear:
            ax1.cla()
        kwargs = {}
        if label:
            kwargs['label']=label
        ax1.semilogx(self.f, self.dBmag, **kwargs)
        ax2 = fig.add_subplot(2,1,2)
        if clear:
            ax2.cla()
        ax2.semilogx(self.f, self.phase, **kwargs)


    def Plot_Mag(self, ax, semilog=False, clear=True, label=None, freqlim=None):
        if clear:
            ax.cla()
        self._Plot_Mag(ax, semilog=semilog, clear=clear, label=label, freqlim=freqlim)
        ax.set_xlabel('Freq (Hz)')
        ax.set_ylabel('Linear Magnitude')
        return ax


    def Plot_dBMag(self, ax, linetype='-', semilog=True, clear=True, label=None, freqlim=None):
        if clear:
            ax.cla()
        self._Plot_prop('dBmag', ax, semilog=semilog, clear=clear, label=label, freqlim=freqlim, linetype=linetype)
        ax.set_xlabel('Freq (Hz)')
        ax.set_ylabel('Magnitude (dB)')
        return ax


    def _Plot_Mag(self, ax, semilog=False, clear=True, label=None, freqlim=None):
        self._Plot_prop('mag', ax, semilog=semilog, clear=clear, label=label, freqlim=freqlim)


    def _Plot_prop(self, prop, ax, xprop = 'f', semilog=False, clear=True, label=None, freqlim=None, linetype='-'):
        myvect = getattr(self, prop)
        myx = getattr(self, xprop)
        if clear:
            ax.cla()
        kwargs = {}
        if label:
            kwargs['label']=label
        if semilog:
            ax.semilogx(myx, myvect, linetype, **kwargs)
        else:
            ax.plot(myx, myvect, linetype, **kwargs)
        if freqlim:
            ax.set_xlim(freqlim)




class AccelMixin:
    """A mixin class intended to add accelerometer integrating
    capabilities to spreadsheet files.  The class is assumed to have
    the following properties defined:

    self.a = the accelerometer signal to integrate

    self.t = the time vector used to find the corresponding dt's

    self.ascale = the calibration factor to get the acceleration
    signal into engineering units

    self.initvel = the intial velocity - optional if v0 is specified
    in the function call"""
    def calc_velocity(self, v0=None, decel=True, tscale=1.0):
        """Calculate velocity based on trapezoidal integration of an
        accelerometer signal.

        decel specifies if the accelerometer signal is assumed to
        be in the opposite direction as self.initvel or v0

        tscale is optional and is analogous to self.ascale - it is
        used to get self.t into appropriate engineering units.  It
        would most likely be used if self.t was in milliseconds
        (tscale would then most likely = 1.0/1000.0)"""
        if v0 is None:
            v0=self.impvel
        deltav = integrate.cumtrapz(self.a*self.ascale, self.t*tscale)
        if decel:
            v = v0-deltav
        else:
            v = v0+deltav
        v=r_[v,v[-1]]
        self.v = v
        return v

    def calc_disp(self, tscale=1.0):
        """Calculate displacement based on trapezoidal integration of
        a velocity vector.  self.v and self.t are assumed to already
        be defined."""
        x = integrate.cumtrapz(self.v, self.t*tscale)
        x=r_[x,x[-1]]
        self.x = x
        return x



class BodeResponse:
    def __init__(self, inputsig, outputsig, timevector):
        self.t = timevector
        self.input = inputsig
        self.output = outputsig
        self.f = makefreqvect(timevector)
        self.spectra = CalcSpectra(inputsig, outputsig)
        self.bode = rwkbode.BodeFromSpectra(self.spectra)
        if len(self.bode.mag) < len(self.f):
            self.f = self.f[0:len(self.bode.mag)]

    def Plot(self, fig, clear=True):
        rwkbode.GenBodePlot(1, self.f, self.bode, fig=fig, clear=clear)


class AveBodeResponse(BodeResponse):
    def __init__(self, inputsig, outputsig, timevector):
        self.t = timevector
        self.input = inputsig
        self.output = outputsig
        self.f = makefreqvect(timevector)
        self.spectra = CalcSpectra(inputsig, outputsig)
        self.bode = rwkbode.AveBodeFromSpectra(self.spectra)
        if len(self.bode.mag) < len(self.f):
            self.f = self.f[0:len(self.bode.mag)]

    def PlotCoh(self, fig, clear=True):
        rwkbode.GenCohPlot(1, self.f, self.bode, fig=fig, clear=clear)

    def PlotBodeCoh(self, fig, clear=True):
        rwkbode.BodeCohPlot(1, self.f, self.bode, fig=fig, clear=clear)

