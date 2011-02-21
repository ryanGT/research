from scipy import log10, shape, zeros, c_, r_, atleast_2d, compress, imag, real, pi, cos, sin, squeeze, where, dot, arange, arctan2
#from math import atan2

#import pylab
#from pylab import figure, cla, subplot, clf, ylabel, xlabel, clf, gca, legend

#import inline_tools
#import c_spec
#from converters import blitz as cblitz
#from rwkdataproc import thresh_py as thresh
from rwkdataproc import thresh, norm2
#from rwkdataproc import mat_atan2, mat_atan2_c
from rwkdataproc import PhaseMassage, CreateDownSampleMask, colwise
from rwkmisc import my_import
import pdb
import os, copy
import mplutil
#reload(mplutil)
import rwkmisc

from  IPython.Debugger import Pdb

#print("In __init__.py of rwkbode")

from matplotlib.ticker import LogFormatterMathtext

class MyFormatter(LogFormatterMathtext):
   def __call__(self, x, pos=None):
       if pos==0: return ''  # pos=0 is the first tick
       else: return LogFormatterMathtext.__call__(self, x, pos)


class bodeout:
    """This class is used by transfer matrix models
    to define a desired bode output."""
    def __init__(self, output='',input='',ind=None,pind=None,type='abs', \
                 post='',dof=None,gain=1.0,gainknown=True, \
                 seedfreq=-1,seedphase=0.0):
        self.output=output
        self.input=input
        self.ind=ind
        self.pind=pind
        self.dof=dof
        self.post=post
        self.type=type
        self.gain=gain
        self.gainknown=gainknown
        self.seedfreq=seedfreq
        self.seedphase=seedphase

    def __repr__(self):
#        print('in bodeout.__str__')
#        attrlist=['output','input','ind','pind','dof']
        attrlist=self.__dict__.keys()
        strout='rwkbode.bodeout instance'
        for item in attrlist:
            strout+='\n'
            strout+=item+'='+str(getattr(self,item))
        return strout


bode_keys = ['input', 'output', 'mag', 'phase', 'coh', \
             'freqlim', 'maglim', 'phaselim', \
             'averaged', 'seedfreq', 'seedphase', \
             'labels', 'legloc']
   
def bode_from_dict(dictin, keys=None):
   if keys is None:
      keys = bode_keys
   my_bode = rwkbode()

   for attr in bode_keys:
      if dictin.has_key(attr):
         setattr(my_bode, attr, dictin[attr])
   return my_bode
      

class rwkbode:
    def __init__(self, output='out', input='in', \
                 mag=None, phase=None, coh=None, \
                 freqlim=[], maglim=[], phaselim=[], \
                 averaged='not specified', \
                 seedfreq=-1, seedphase=0,
                 labels=[], legloc=-1, compin=[]):
        self.output = output
        self.input = input
        if len(compin) > 0:
           if mag is None:
              self.mag = squeeze(colwise(abs(compin)))
           if phase is None:
              self.phase = squeeze(colwise(arctan2(imag(compin),real(compin))*180.0/pi))
        else:
            self.mag = squeeze(mag)
            self.phase = squeeze(phase)
        self.coh = coh
        self.averaged = averaged
        self.seedfreq = seedfreq
        self.seedphase = seedphase
        self.freqlim = freqlim
        self.maglim = maglim
        self.phaselim = phaselim
        self.labels = labels
        self.legloc = legloc

    def to_dict(self, keys=None):
       """Convert rwkbode instance into a dictionary for saving."""
       if keys is None:
          keys = bode_keys
       mydict = {}
       for key in keys:
          val = getattr(self, key)
          mydict[key] = val
       return mydict
    
       
    def __repr__(self):
#        print('in bodeout.__str__')
#        attrlist=['output','input','ind','pind','dof']
        attrlist=self.__dict__.keys()
        strout='rwkbode.rwkbode instance'
        vectlist=['mag','phase','coh']
        for item in attrlist:
            strout+='\n'
            if item in vectlist:
                strout+=item+', shape='+str(shape(getattr(self,item)))
            else:
                strout+=item+'='+str(getattr(self,item))
        return strout

    def copybodeprops(self, source, proplist=['freqlim','maglim','phaselim','averaged','seedfreq','seedphase']):
        """Copy the values of proplist from source rwkbode.  This function is used to preserve information when preforming transfer function math operations on rwkbode instances (i.e. if I add a constant to an rwkbode, I create a new rwkbode instance but may want to preserve some of the old information)."""
        for curprop in proplist:
            if hasattr(source,curprop):
                setattr(self,curprop,getattr(source,curprop))

    def __radd__(self,other):
        return self.__add__(other)

    def __add__(self,other):
        """Add two bodes or a bode and a constant.  Basically, if
        other is an rwkbode, call its ToComp method, and add the two
        complex vectors and convert back to an rwkbode, otherwise
        assume that other can simply be added to the complex vector
        (i.e. it is either a constant or a complex vector)."""
        if hasattr(other,'ToComp'):
            ocomp=other.ToComp()
        else:
            ocomp=other
        compout=self.ToComp()+ocomp
        bodeout=rwkbode(self.output,self.input,compin=compout)
        bodeout.copybodeprops(source=self)
        return bodeout

    def __div__(self,other):
        """Divide two bodes or a bode and a constant.  Basically, if
        other is an rwkbode, call its ToComp method, and divide the
        two complex vectors and convert back to an rwkbode, otherwise
        assume that self.ToComp() can simply be divided by other
        (i.e. it is either a constant or a complex vector)."""
        if hasattr(other,'ToComp'):
            ocomp=other.ToComp()
        else:
            ocomp=other
        compout=self.ToComp()/ocomp
        if hasattr(other, 'output'):
           new_input = other.output
        else:
           new_input = self.input#don't know how else to handle this
        bodeout=rwkbode(self.output,new_input,compin=compout)
        bodeout.copybodeprops(source=self)
        return bodeout

    def dB_Phase_Error_vect(self,other):
        """Subtract two bodes dBmag and phase and return a Bode with
        dBmag and phase."""
        outbode=rwkbode()
        outbode.phase=squeeze(self.phase)-squeeze(other.phase)
        outbode.dBmag=self.dBmag()-other.dBmag()
        return outbode

    def dB_Phase_Error_sum(self,other, phaseweight=0.1):
        """Subtract two bodes dBmag and phase and return the sum of
        the squared error."""
        phaseerror = squeeze(self.phase)-squeeze(other.phase)
        e_phase = (phaseerror**2).sum()
        dBmagerror = self.dBmag()-other.dBmag()
        e_dB = (dBmagerror**2).sum()
        return e_dB + e_phase*phaseweight
       
    def __mul__(self,other):
        """Multiply two bodes."""
        if type(other)==float or type(other)==int:
            myoutput=copy.deepcopy(self)
            myoutput.mag=myoutput.mag*other
            return myoutput
        myin='in'
        myout='out'
        match=1
        if self.output==other.input:
            first=self
            second=other
        elif self.input==other.output:
            first=other
            second=self
        else:
            warnme=1
            if (self.input=='in' and self.output=='out') or (other.input=='in' and other.output=='out'):
                warnme=0
            if warnme:
                print('Warning: multiplying Bodes without a matching input/output pair:\n'+self.output+'/'+self.input+ ' * ' +other.output+'/'+other.input)
            match=0
            first=self
            second=other
        if match:
            myin=first.input
            myout=first.output
        myoutput=copy.deepcopy(self)
        myoutput.input=myin
        myoutput.output=myout
        myoutput.mag=squeeze(colwise(first.mag)*colwise(second.mag))
        myoutput.phase=squeeze(colwise(first.phase)+colwise(second.phase))
        return myoutput
        
    def ToComp(self,ind=None):
        if ind!=None:
            x=self.mag[ind]*cos(pi/180.0*self.phase[ind])
            y=self.mag[ind]*sin(pi/180.0*self.phase[ind])
            return x+y*1.0j
        else:
            x=squeeze(self.mag)*squeeze(cos(pi/180.0*self.phase))
            y=squeeze(self.mag)*squeeze(sin(pi/180.0*self.phase))
        return squeeze(x+y*1.0j)

    def real(self):
        return (self.ToComp()).real

    def imag(self):
        return (self.ToComp()).imag

    def truncate(self, freq, flow, fhigh=None):
        """Truncate the mag, phase, and coherence of the rwkbode
        instance based on the indices returned by thresh(flow) and
        thresh(fhigh).  If fhigh is not given, it is assumed that flow
        is a list of [flow, fhigh]."""
        if fhigh is None:
            fhigh=flow[1]
            flow=flow[0]
        i1=thresh(freq,flow)
        i2=thresh(freq,fhigh)
        if i2-i1<max(shape(self.mag)):#test if already truncated
            self.mag=colwise(self.mag)[i1:i2,:]
            self.phase=colwise(self.phase)[i1:i2,:]
            self.coh=colwise(self.coh)[i1:i2,:]
        tfreq=freq[i1:i2]
        if self.freqlim:
            if self.freqlim[0]<min(tfreq):
                self.freqlim[0]=min(tfreq)
            if self.freqlim[1]>max(tfreq):
                self.freqlim[1]=max(tfreq)
        else:
            self.freqlim=[min(tfreq),max(tfreq)]
        return tfreq

    def downsample(self,freq,factor,ranges=[]):
        if shape(factor) and ranges:
#            pdb.set_trace()
            mask=CreateVariableMask(freq, ranges, factor)
        else:
            mask=CreateDownSampleMask(freq, factor)
#            mask=CreateLogMask(freq)
        dsfreq=compress(mask,freq)
        mylist=['mag','phase','coh']
        for item in mylist:
            tempvect=getattr(self,item)
            if tempvect is not None:
                tempvect=colwise(tempvect)
                setattr(self,item,compress(mask,tempvect,0))
        return dsfreq

    def compress(self, mask):
        mylist=['mag','phase','coh']
        for item in mylist:
            tempvect=getattr(self,item)
            if tempvect is not None:
                tempvect=colwise(tempvect)
                setattr(self,item,tempvect.compress(mask,0))

    def CrossoverFreq(self, freqin):
       t1=squeeze(self.dBmag() > 0.0)
       t2=r_[t1[1:],t1[0]]
       t3=(t1 & -t2)
       myinds=where(t3)[0]
       if not myinds.any():
          return None, []
       maxind=max(myinds)
       return freqin[maxind], maxind

    def PhaseMargin(self,freqin):
       fc,ind=self.CrossoverFreq(freqin)
       if not fc:
          return 180.0
       return 180.0+squeeze(self.phase[ind])

    def PhaseMassage(self,freqin,jump=250):
        if not self.seedfreq >0:
            print('seedfreq and seedphase must be specified before calling PhaseMassage')
            return
        self.phase=PhaseMassage(self.phase,freqin,self.seedfreq,self.seedphase,jump)

    def dBmag(self):
        return squeeze(20*log10(self.mag))

    def automag(self,freqvect,margin=0.1):
#        pdb.set_trace()
        return self.autolim("mag",freqvect,margin,db=1)

    def autophase(self,freqvect,margin=0.1):
        return self.autolim("phase",freqvect,margin)

    def autolim(self, myattr, freqvect, margin=0.1,db=0):
        if self.freqlim:
            ind1=thresh(freqvect,self.freqlim[0])
            ind2=thresh(freqvect,self.freqlim[1])
        else:
            ind1=0
            ind2=-1
        mymatrix=getattr(self,myattr)
        if len(shape(mymatrix))==1:
            submat=mymatrix[ind1:ind2]
        else:
            mymatrix=colwise(mymatrix)
            submat=mymatrix[ind1:ind2,:]
        if db:
            submat=20*log10(submat)
        # max and min need to be done columnwise
        # (maybe a Krauss matrix max)
        if len(shape(submat))==2:
            mymax=[]
            mymin=[]
            for q in range(shape(submat)[1]):
                mymax.append(max(submat[:,q]))
                mymin.append(min(submat[:,q]))
        else:
            mymax=max(submat)
            mymin=min(submat)
        if len(shape(mymax))>0:
            mymax=max(mymax)
            mymin=min(mymin)
        myspan=mymax-mymin
        mymargin=margin*myspan
        limout=[mymin-mymargin, mymax+mymargin]
        setattr(self,myattr+"lim",limout)
        return limout

def BodeFromSpectra(specin,calccoh=False):
    H=specin.Gxy/specin.Gxx
    if calccoh:
       cohnum=squeeze(norm2(specin.Gxy))
       cohden=squeeze(specin.Gxx*specin.Gyy)
       coh=cohnum/cohden
       ave=True
    else:
       coh=[]
       ave=False
    mybode=rwkbode(input=specin.input, output=specin.output, \
                   compin=squeeze(H), coh=coh, averaged=ave)
    mybode.Gxx=specin.Gxx
    mybode.Gyy=specin.Gyy
    mybode.Gxy=specin.Gxy
    return mybode

def AveBodeFromSpectra(specin):
    Gxxave=squeeze(specin.Gxx.mean(1))
    Gxyave=squeeze(specin.Gxy.mean(1))
    Gyyave=squeeze(specin.Gyy.mean(1))
    Have=Gxyave/Gxxave
    cohnum=norm2(Gxyave)
    cohden=Gxxave*Gyyave
    mybode=rwkbode(input=specin.input, output=specin.output, \
                   compin=squeeze(Have), coh=squeeze(cohnum/cohden))
    return mybode
    
   
#ttttttttttttttttttttttttttttttttttttttttttt
#
#   Temporarily commented out while switching to WX embedded compatible functions
#
#ttttttttttttttttttttttttttttttttttttttttttt
## def PreviewTrunc(startfi, freq, bodelist, freqlim):
##     """This function plots each of a list of Bodes
##     with vertical bars showing where truncation
##     would take place."""
##     pylab.ioff()
##     for x,curbode in enumerate(bodelist):
##         GenBodePlot(x+startfi,freq,curbode)
##         pylab.figure(x+startfi)
##         pylab.subplot(211)
##         AddVerticalLines(freqlim)
##         pylab.subplot(212)
##         AddVerticalLines(freqlim)
#ttttttttttttttttttttttttttttttttttttttttttt

def DownsampleList(freq,bodelist,factor,ranges=[]):
    """Downsample each bode in bode returning
    a downsampled frequency vector and a list
    of copies of the bodes that are downsampled."""
    bodesout=[]
    first=1
    for curbode in bodelist:
        tempbode=copy.deepcopy(curbode)
        if first:
            dsf=tempbode.downsample(freq,factor,ranges)
            first=0
        else:
            tempbode.downsample(freq,factor,ranges)
        bodesout.append(tempbode)
    return dsf, bodesout

def TruncList(freq,bodelist,freqlim):
    """Truncate each bode in bodelist returning
    a truncated frequency vector and a list of 
    copies of the bodes that are trucnated."""
    bodesout=[]
    first=1
    for curbode in bodelist:
        tempbode=copy.deepcopy(curbode)
        if first:
            tf=tempbode.truncate(freq,freqlim)
            first=0
        else:
            tempbode.truncate(freq,freqlim)
        bodesout.append(tempbode)
    return tf, bodesout
        
        
def AddVerticalLines(freqlim,ymin=-1000,ymax=1000):
    curylim=pylab.ylim()
    curx=pylab.xlim()
    myyvect=[ymin,ymax]
    pylab.plot([freqlim[0],freqlim[0]],myyvect)
    pylab.plot([freqlim[1],freqlim[1]],myyvect)
    pylab.ylim(curylim)
    pylab.xlim(curx)


def CombinedBodes(bodelist):
    """This function seeks to make one combined rwkbode instance from
    a list of them.  This may be useful in the case of having several
    experimental Bode data sets with various targeted frequency ranges
    that need to be treated as one data set."""
    docoh=True
    for cb in bodelist:
        if cb.coh is None:
            docoh=False
            break
    first=1
#    pdb.set_trace()
    bigcoh=[]
    for cb in bodelist:
        if first:
            first=0
            bigmag=colwise(cb.mag)
            bigphase=colwise(cb.phase)
            if docoh:
                bigcoh=colwise(cb.coh)
        else:
            bigmag=r_[bigmag,colwise(cb.mag)]
            bigphase=r_[bigphase,colwise(cb.phase)]
            if docoh:
                bigcoh=r_[bigcoh,colwise(cb.coh)]
    myoutbode=copy.deepcopy(bodelist[0])
    myoutbode.mag=squeeze(bigmag)
    myoutbode.phase=squeeze(bigphase)
    myoutbode.coh=squeeze(bigcoh)
    myoutbode.freqlim=[]
    return myoutbode

def renew(obi):
    mylist=['mag','phase','coh','freqlim','maglim','phaselim','averaged','seedfreq','seedphase','labels','legloc']
    mydict={}
    for item in mylist:
        mydict[item]=copy.deepcopy(getattr(obi,item))
    return rwkbode(obi.output, obi.input, **mydict)
#obi.mag, obi.phase, obi.coh, freqlim=obi.freqlim, maglim=obi.maglim,phaselim=obi.phaselim, averaged=obi.averaged, obi.seedfreq=obi.seedfreq, seedphase=obi.seedphase, labels=obi.labels, legloc=obi.legloc)

def concatenate(bodelist):
    bodeout=copy.deepcopy(bodelist[0])
    for curbode in bodelist[1:]:
        bodeout.mag=c_[colwise(bodeout.mag),colwise(curbode.mag)]
        bodeout.phase=c_[colwise(bodeout.phase),colwise(curbode.phase)]
        if bodeout.coh and curbode.coh:
            bodeout.coh=c_[colwise(bodeout.coh),colwise(curbode.coh)]
    return bodeout

def CreateVariableMask(vector, ranges, factors):
    mask=[]
    prev=0
    for range, factor in zip(ranges,factors):
        if type(range)==list:
            i1=thresh(vector,range[0])
            i2=thresh(vector,range[1])
        else:
            i1=thresh(vector,prev)
            i2=thresh(vector,range)
            prev=range
        tempv=vector[i1:i2]
        tmask=CreateDownSampleMask(tempv,factor)
        mask.extend(tmask)
    pad=len(vector)-len(mask)
    mypad=[0]*pad
    mask.extend(mypad)
    return mask


def CreateLogMask(vector,factor=1.1):
    N=1
    mask=[1]
    tmask=[0]*N
    while len(mask)+len(tmask)<=len(vector):
        mask.extend(tmask)
        mask.append(1)
        N=N*factor
        tmask=[0]*int(N)
    pad=len(vector)-len(mask)
    mypad=[0]*pad
    mask.extend(mypad)
    return mask

#ttttttttttttttttttttttttttttttttttttttttttt
#
#    Temporarily commented out while working
#    on the switch to WX embedded compatible functions
#
#ttttttttttttttttttttttttttttttttttttttttttt
## def GenCohPlot(fignum,freqvect,bodein,clear=True,filename='',mydpi=75,folder='figs',linestyle='',legend=[],legloc=-1,autoY=1,linewidth=None):
## #    pdb.set_trace()
##     pylab.figure(fignum)
##     if clear:
##         pylab.cla()
##     if len(shape(bodein.coh))>1:
##         bodein.coh=colwise(bodein.coh)
##         if shape(bodein.coh)[1]==1:
##             bodein.coh=bodein.coh[:,0]
##         else:
##             raise IndexError, 'Coherence must be a vector or a matrix with only one column. Got something with shape '+str(shape(bodein.coh))
##     myargs=[freqvect,bodein.coh]
##     mykwargs={}
##     if linestyle:
##         myargs.append(linestyle)
##     if linewidth:
##         mykwargs['linewidth']=linewidth
##     pylab.semilogx(*myargs,**mykwargs)
##     if bodein.freqlim:
##         pylab.axis(r_[bodein.freqlim,[0,1]])
## #    pylab.semilogx(freqvect,bodein.coh)
##     if legend:
##         legargs=[legend]
##         if legloc>=0:
##             legargs.append(legloc)
##         pylab.legend(*legargs)
##     pylab.xlabel('Freq (Hz)')
##     pylab.ylabel('Coherence')

##     if filename:
##         if not os.path.exists(folder):
##             os.mkdir(folder)
##         curpath=os.path.join(folder,filename)
##         pylab.savefig(curpath,dpi=mydpi)
#tttttttttttttttttttttttttttttttttttttttttt

def BodeFromComp(compin, expbode, freq, PhaseMassage=True):
    """Create an rwkbode instance from a complex vector and an example
    bode.  This would typically be used in conjunction with a symbolic
    bode function that outputs complex values, where an expbode is
    known to have the same input and output and desired seedfreq and
    seedphase.  The freq input is used to run PhaseMassage on the
    created rwkbode instance."""
    tbode=rwkbode(expbode.output,expbode.input,compin=compin)
    tbode.seedfreq=expbode.seedfreq
    tbode.seedphase=expbode.seedphase
    tbode.freqlim=expbode.freqlim
    if (tbode.seedfreq>0) & PhaseMassage:
#        Pdb().set_trace()
        tbode.PhaseMassage(freq)
    tbode.phase=squeeze(tbode.phase)
    tbode.mag=squeeze(tbode.mag)
    return tbode

def Bode_From_TF(TF, freq, input='', output=''):
    comp = TF.FreqResp(freq, fignum=None)
    return rwkbode(output=output, input=input, \
                   compin=comp)
   
def BodeFromModname(modname,expbode,f,funcname=None,ucv=[],optargs=(),PhaseMassage=True):
    """Generate one Bode from its module name, i.e. the name of a
    module that contains the symbolic Bode function funcname.  If
    funcname is not given, it is assumed to be the same as modname.
    expbode is used to determine the inputs and outputs as well as the
    seedfreq and seedphase."""
    s=2.0j*pi*f
    if not modname:
        bodefunc=funcname
    elif callable(modname):
        bodefunc=modname
    else:
        bodemod=my_import(modname)
        if funcname:
            bodefunc=getattr(bodemod,funcname)
        else:
            bodefunc=getattr(bodemod,modname)
#    if ucv:
    if ucv==[]:
        compv=bodefunc(s,*optargs)
    else:
        compv=bodefunc(s,ucv,*optargs)
    mybode=BodeFromComp(compv, expbode, f, PhaseMassage)
    return mybode

def BodesFromCompListbyModname(modpattern,expbodes,f,funcname=None, ucv=[]):
    """This function takes as inputs a modpattern that is a pattern
    that will be turned into modulnames by calling % x in an enumerate
    loop over the list of expbodes (assumed to be in the same order as
    _bode0, _bode1,...  when enumerating to get the modnames) and an
    input frequency vector and returns a list of model bodes.  If
    funcname is not specified, the same name as each module will be
    used.
    
    If ucv is given, it is passed to the bode function:
    bode=bodefunc(s,ucv) (ucv stands for unknown coefficent vector and
    is used in my curve fitting stuff)."""
    if type(expbodes)!=list and type(expbodes)!=tuple:
        expbodes=[expbodes]
    s=2.0j*pi*f
    bodesout=[]
    for x,curexp in enumerate(expbodes):
        curmod=modpattern % x
        mybode=BodeFromModname(curmod,curexp,f,funcname,ucv)
        bodesout.append(mybode)
    return bodesout

def PlotBodeList(startfi, freq, bodelist, **kwargs):
    """Plot a list of Bodes by calling GenBodePlot 
    for each bode in bodelist."""
    pylab.ioff()
    for x,cb in enumerate(bodelist):
        GenBodePlot(x+startfi,freq,cb,**kwargs)
        

def SetAxisLimits(startfi,freqlist,listofbodelists):
    """Take a starting figure number (startfi), a list of frequency
    vectors (freqlist), and a list of bode lists (listofbodelists)
    that have been overlayed and set sensible axis limits.  Each bode
    must have freqlim defined, or it will be skipped."""
    translist=rwkmisc.transposed(listofbodelists)
    for x,curlist in enumerate(translist):
        curfi=startfi+x
        xmin=None
        xmax=None
        magmin=None
        magmax=None
        phasemin=None
        phasemax=None
        for curf,curbode in zip(freqlist,curlist):
            if curbode.freqlim:
                if (not xmin) or min(curbode.freqlim)<xmin:
                    xmin=min(curbode.freqlim)
                if (not xmax) or max(curbode.freqlim)>xmax:
                    xmax=max(curbode.freqlim)
                curmaglim=curbode.automag(curf)
                curphaselim=curbode.autophase(curf)
                if (not magmin) or min(curmaglim)<magmin:
                    magmin=min(curmaglim)
                if (not magmax) or max(curmaglim)>magmax:
                    magmax=max(curmaglim)
                if (not phasemin) or min(curphaselim)<phasemin:
                    phasemin=min(curphaselim)
                if (not phasemax) or max(curphaselim)>phasemax:
                    phasemax=max(curphaselim)
        if xmax and xmin:
            setfreqlim=[xmin,xmax]
            rwkmplutil.SetBothXlims(curfi,setfreqlim) 
        if magmin and magmax:
            rwkmplutil.SetMagLim(curfi,[magmin,magmax])
        if phasemin and phasemax:
            rwkmplutil.SetPhaseLim(curfi,[phasemin,phasemax])

mytypes=['-','--',':','-.']
#colors='ygrbckm'
colors=['b','y','r','g','c','k']#['y','b','r','g','c','k']

def _getlinestyle(ax=None):
    if ax is None:
       import pylab
       ax = pylab.gca()
    myind=ax._get_lines.count
    return {'color':colors[myind % len(colors)],'linestyle':mytypes[myind % len(mytypes)]}

def _inccount():
    ax=gca()
    ax._get_lines.count+=1
    
def _PlotMatrixvsF(freqvect, matin, linestyle='', \
                   linewidth=None, semilogx=True, \
                   allsolid=False, axis=None, label=None):
    mykwargs={}
    usepylab = False
    if axis is None:
       import pylab
       axis = pylab.gca()
       usepylab = True
    if len(shape(matin))==1:
        myargs=[freqvect,matin]
        if linestyle:
            myargs.append(linestyle)
        else:
            mykwargs.update(_getlinestyle(axis))
        if linewidth:
            mykwargs['linewidth'] = linewidth
        if label is not None:
           mykwargs['label'] = label
        if semilogx:
            curline,=axis.semilogx(*myargs,**mykwargs)
        else:
            curline,=axis.plot(*myargs,**mykwargs)
        mylines=[curline]
#        _inccount()
    else:
        mylines=[]
        for q in range(shape(matin)[1]):
            myargs=[freqvect,matin[:,q]]
            if linestyle:
                myargs.append(linestyle)
            else:
                mykwargs.update(_getlinestyle(axis))
            if linewidth:
                mykwargs['linewidth']=linewidth
            if label is not None:
               mykwargs['label'] = label
            if semilogx:
                curline,=axis.semilogx(*myargs,**mykwargs)
            else:
                curline,=axis.plot(*myargs,**mykwargs)
            mylines.append(curline)
#            _inccount()
    return mylines


def _PlotMag(freqvect,bodein,linestyle='',linewidth=0, \
             axis=None, label=None):
    if callable(bodein.dBmag):
        myvect=bodein.dBmag()
    else:
        myvect=bodein.dBmag
    return _PlotMatrixvsF(freqvect, myvect, linestyle=linestyle, \
                          linewidth=linewidth, axis=axis, label=label)

def _PlotPhase(freqvect,bodein,linestyle='',linewidth=0, axis=None, \
               label=None):
    return _PlotMatrixvsF(freqvect,bodein.phase,linestyle=linestyle, \
                          linewidth=linewidth, axis=axis, label=label)

def _PlotCoh(freqvect,bodein,linestyle='',linewidth=0, axis=None, \
             label=None):
    return _PlotMatrixvsF(freqvect,bodein.coh, linestyle=linestyle, \
                          linewidth=linewidth, axis=axis, label=label)

#ttttttttttttttttttttttttttttttttttttttttttt
#
#    Temporarily commented out while working
#    on the switch to WX embedded compatible functions
#
#ttttttttttttttttttttttttttttttttttttttttttt
## def MagPlot(fignum,freqvect,bodein,clear=True,legendsub=111,legloc=3,**kwargs):
##     figure(fignum)
##     myargs=['linestyle','colors','linewidth']
##     subkwargs={}
##     for key in myargs:
##         if kwargs.has_key(key):
##             subkwargs[key]=kwargs[key]
##     if clear:
##         clf()
##     subplot(111)
##     if clear:
##         cla()
##     mylines=_PlotMag(freqvect,bodein,**subkwargs)
##     ylabel('Mag. Ratio (dB)')
##     xlabel('Freq. (Hz)')
##     ax=gca()
##     ax.xaxis.set_major_formatter(MyFormatter())
##     if kwargs.has_key('freqlim'):
##         rwkmplutil.SetXlims(fignum,kwargs['freqlim'],[111])
##     if kwargs.has_key('legend'):
##         myargs=(kwargs['legend'],)
##         if kwargs.has_key('legloc'):
##             myargs+=(kwargs['legloc'],)
##         legloc(*myargs)
#ttttttttttttttttttttttttttttttttttttttttttt


#ttttttttttttttttttttttttttttttttttttttttttt
#
#    Temporarily commented out while working
#    on the switch to WX embedded compatible functions
#
#ttttttttttttttttttttttttttttttttttttttttttt
## def AspectPlot(fignum, freqvect, bodein, aspectstr, semilogx=True, clear=True, legendsub=111,legloc=3,**kwargs):
##     figure(fignum)
##     myargs=['linestyle','colors','linewidth']
##     subkwargs={}
##     for key in myargs:
##         if kwargs.has_key(key):
##             subkwargs[key]=kwargs[key]
##     if clear:
##         clf()
##     subplot(111)
##     if clear:
##         cla()
##     mything=getattr(bodein, aspectstr)
##     if callable(mything):
##         mymatrix=mything()
##     else:
##         mymatrix=mything
##     mymatrix=squeeze(mymatrix)
##     mylines=_PlotMatrixvsF(freqvect, mymatrix, semilogx=semilogx, **subkwargs)
##     ylabel('Mag. Ratio (dB)')
##     xlabel('Freq. (Hz)')
## #    ax=gca()
## #    ax.xaxis.set_major_formatter(MyFormatter())
##     if kwargs.has_key('freqlim'):
##         rwkmplutil.SetXlims(fignum,kwargs['freqlim'],[111])
##     if kwargs.has_key('legend'):
##         myargs=(kwargs['legend'],)
##         if kwargs.has_key('legloc'):
##             myargs+=(kwargs['legloc'],)
##         legloc(*myargs)
#ttttttttttttttttttttttttttttttttttttttttttttttt


def BodeCohPlot(fignum, freqvect, bodein, \
                clear=True, legendsub=311, legloc=3, \
                fig=None, **kwargs):
    if fig is None:
       from pylab import figure
       fig = figure(fignum,figsize=(8,8))
    myargs=['linestyle','colors','linewidth']
    subkwargs={}
    for key in myargs:
        if kwargs.has_key(key):
            subkwargs[key]=kwargs[key]
    if clear:
        fig.clf()
    ax1 = fig.add_subplot(3, 1, 1)
    if clear:
        ax1.cla()
    _PlotCoh(freqvect, bodein, axis=ax1 , **subkwargs)
    ax1.set_ylabel('Coherence')
    ax1.xaxis.set_major_formatter(MyFormatter())
    ax2=fig.add_subplot(3, 1, 2)
    if clear:
        ax2.cla()
    mylines=_PlotMag(freqvect, bodein, axis=ax2, **subkwargs)
    ax2.set_ylabel('Mag. Ratio (dB)')
    ax2.xaxis.set_major_formatter(MyFormatter())
    ax3=fig.add_subplot(3, 1, 3)
    if clear:
        ax3.cla()
    mylines=_PlotPhase(freqvect, bodein, axis=ax3, **subkwargs)
    ax3.set_ylabel('Phase (deg.)')
    ax3.set_xlabel('Freq. (Hz)')
    ax3.xaxis.set_major_formatter(MyFormatter())
##     if kwargs.has_key('freqlim'):
##         rwkmplutil.SetXlims(fignum,kwargs['freqlim'],[311,312,313])
##     if kwargs.has_key('legend'):
##         subplot(legendsub)
##         myargs=(kwargs['legend'],)
##         if kwargs.has_key('legloc'):
##             myargs+=(kwargs['legloc'],)
##         legloc(*myargs)
    return fig


def GenBodePlot(fignum, freqvect, bodein, clear=True, \
                legend_axis=1, legloc=3, fig=None, **kwargs):
    if fig is None:
       from pylab import figure
       fig = figure(fignum)
    
    if type(bodein)==type(arange(0,1,0.01)):
        bodein=rwkbode(compin=bodein)
    #myargs=['linestyle','colors','linewidth','label']
    myargs=['linestyle','colors','linewidth','label']
    subkwargs={}
    for key in myargs:
        if kwargs.has_key(key):
            subkwargs[key]=kwargs[key]
    if clear:
        fig.clf()
    ax1 = fig.add_subplot(2,1,1)
    if clear:
        ax1.cla()
    myind=ax1._get_lines.count
    mylines=_PlotMag(freqvect,bodein, axis=ax1,**subkwargs)
    ax1.set_ylabel('Mag. Ratio (dB)')
#    ax=gca()
    ax1.xaxis.set_major_formatter(MyFormatter())
    ax2 = fig.add_subplot(2,1,2, sharex=ax1)
    if clear:
        ax2.cla()
    mylines=_PlotPhase(freqvect,bodein, axis=ax2,**subkwargs)
    ax2.set_ylabel('Phase (deg.)')
    ax2.set_xlabel('Freq. (Hz)')
    ax2.xaxis.set_major_formatter(MyFormatter())
    if kwargs.has_key('legend'):
       mplutil.Legend(kwargs['legend'], fig, axis=legend_axis, \
                      loc=legloc)
#ttttttttttttttttttttttttttttt
#    if kwargs.has_key('freqlim'):
#        rwkmplutil.SetXlims(fignum,kwargs['freqlim'],[211,212])
#        subplot(legendsub)
#        myargs=(kwargs['legend'],legloc)
#        legend(*myargs)
#ttttttttttttttttttttttttttttt
    return fig


def GenCohPlot(fignum, freqvect, bodein, \
               clear=True, legendsub=212, legloc=3, \
               fig=None, **kwargs):
    if fig is None:
       from pylab import figure
       fig = figure(fignum)
    assert isinstance(bodein, rwkbode), "GenCohPlot must be called with an rwkbode instance for bodein."

    myargs=['linestyle','colors','linewidth']
    subkwargs={}
    for key in myargs:
        if kwargs.has_key(key):
            subkwargs[key]=kwargs[key]
    if clear:
        fig.clf()
    ax = fig.add_subplot(1,1,1)
    if clear:
        ax.cla()
    myind=ax._get_lines.count
    mylines=_PlotCoh(freqvect,bodein, axis=ax,**subkwargs)
    ax.set_ylabel('Coherence')
#    ax=gca()
    ax.xaxis.set_major_formatter(MyFormatter())
    ax.set_xlabel('Freq. (Hz)')
    return fig
    
    

## def GenBodePlot(fignum,freqvect,bodein,clear=True,filename='',mydpi=75,folder='figs',linestyle='',colors=[],legend=[],legloc=-1,autoY=1,speciallegend=False,linewidth=0,legonphase=False):
##     #ticksize='large',labelsize='x-large',
##     pylab.ioff()
##     bodein.mag=colwise(bodein.mag,makecopy=False)    
##     bodein.phase=colwise(bodein.phase,makecopy=False)
##     pylab.figure(fignum)
##     pylab.subplot(211)

##     mykwargs={}
##     if clear:
##         pylab.cla()

##     if len(shape(bodein.mag))==1:
##         myargs=[freqvect,20*log10(bodein.mag)]
##         if linestyle:
##             myargs.append(linestyle)
##         elif colors:
##             myargs.append(colors[0])
##         curline,=pylab.semilogx(*myargs)
##         mylines=[curline]
##     else:
##         mylines=[]
##         for q in range(shape(bodein.mag)[1]):
##             myargs=[freqvect,20*log10(bodein.mag[:,q])]
##             if linestyle:
##                 myargs.append(linestyle)
##             elif colors:
##                 myargs.append(colors[q % len(colors)])                
##             if linewidth:
##                 mykwargs['linewidth']=linewidth
##             curline,=pylab.semilogx(*myargs,**mykwargs)
##             mylines.append(curline)
##     if bodein.freqlim:
##         if not bodein.maglim:
##             autoY=1
##         if autoY:
##             pylab.axis(r_[bodein.freqlim,bodein.automag(freqvect)])
##         else:
##             pylab.axis(r_[bodein.freqlim,bodein.maglim])
##     if legend:
##         if speciallegend:
##             legargs=[mylines]
##         else:
##             legargs=[]
##         legargs.append(legend)
##         if legloc>=0:
##             legargs.append(legloc)
##     if (legend and not legonphase):
##         pylab.legend(*legargs)
## #    locs , labels = pylab.xticks ()
## #    pylab.set(labels , size=ticksize)
## ##    xticks([1,10],[r'$10^0$',r'$10^1$'],size=ticksize,family='monospace')
## #    locs , labels = pylab.yticks()
## #    pylab.set(labels , size=ticksize)
##     pylab.ylabel('Mag. Ratio (dB)')#,fontsize=labelsize)
##     ax=pylab.gca()
##     ax.xaxis.set_major_formatter(MyFormatter())
## ##    draw()

##     pylab.subplot(212)
##     if clear:
##         pylab.cla()
##     if len(shape(bodein.phase))==1:
##         plotstr='pylab.semilogx(freqvect,bodein.phase'
##         if len(linestyle)>0:
##             plotstr+=',linestyle'
##         if len(colors)>0:
##             plotstr+=',colors[0]'
##         plotstr+=')'
##         exec(plotstr)
##     else:
##         for q in range(shape(bodein.phase)[1]):
##             myargs=[freqvect,bodein.phase[:,q]]
##             if linestyle:
##                 myargs.append(linestyle)
##             elif colors:
##                 myargs.append(colors[q % len(colors)])
##             pylab.semilogx(*myargs,**mykwargs)
##     if bodein.freqlim:
##         if autoY:
##             pylab.axis(r_[bodein.freqlim,bodein.autophase(freqvect)])
##         else:
##             pylab.axis(r_[bodein.freqlim,bodein.phaselim])
## #    locs , labels = pylab.xticks ()
## #    pylab.set(labels , size=ticksize)
## ##    xticks([1,10],[r'$10^0$',r'$10^1$'],size=ticksize,family='monospace')
## #    locs , labels = pylab.yticks()
## #    pylab.set(labels , size=ticksize)
##     pylab.xlabel('Freq (Hz)')#,fontsize=labelsize)
##     pylab.ylabel('Phase (deg.)')#,fontsize=labelsize)
##     if (legend and legonphase):
##         pylab.legend(*legargs)
##     ax=pylab.gca()
##     ax.xaxis.set_major_formatter(MyFormatter())
    
##     if len(filename)>0:#use this as the test of whether ot not to save
##         if not os.path.exists(folder):
##             os.mkdir(folder)
##         curpath=os.path.join(folder,filename)
##         pylab.savefig(curpath,dpi=mydpi)

def PrintBodeList(listin):
    for x,item in enumerate(listin):
        print(str(x)+': '+item.output + ' vs. '+item.input)

def FindMatch(bodelist, output, input):
    foundmatch=0
    for bode in bodelist:
        if bode.input==input and bode.output==output:
            bodeout=bode
            foundmatch=1
            break
    if foundmatch:
        return bodeout
    else:
        raise IndexError, 'Could not find bode with output='+output +' and input='+input

   
