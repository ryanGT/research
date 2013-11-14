#import Numeric
#from Numeric import *
from numpy import *
import numpy
from math import atan2
import sys, os, glob, time
#sys.path.insert(0,'..')
#t1=time.time()
#import inline_tools
#import c_spec
#t1b=time.time()
#from converters import blitz as cblitz
#t2=time.time()
#from scipy.io.mio import loadmat, savemat
from scipy.io import loadmat, savemat
from scipy import r_, c_, mat, shape, iscomplexobj, multiply, \
     shape, atleast_2d , fft, pi, zeros, where, column_stack, \
     arange, squeeze, array, fft

#from pylab import figure, cla, ion, ioff, show, plot, xlabel, ylabel, title, legend, semilogx, subplot
import rwkmisc, copy
import pdb
#import rwkascii#, rwkbode
import scipy
from rwkmisc import my_import
import types
#from IPython.core.debugger import Pdb
#import rwkmplutil

def is_1D(arrayin):
    shape_test = squeeze(arrayin).shape
    if len(shape_test) > 1:
        return False
    else:
        return True


def colwise(matin, makecopy=1):
    if type(matin) != ndarray:
        matin = array(matin)
    if makecopy:
        tempmat = matin.copy()
    else:
        tempmat = matin
#    t2=time.time()
#    print('copy time='+str(t2-t1))
    matout = numpy.atleast_2d(tempmat)
    myshape = matout.shape
    if myshape[0] < myshape[1]:
        matout = matout.T
    return matout

def rowwise(matin, makecopy=1):
    col_mat = colwise(matin, makecopy=makecopy)
    return col_mat.T


def python_phase_engine(vectin, seedphase, startind=1, stopind=None,
                        step=1, jump=250.0):
    vectin = colwise(vectin)
    for cur_col in vectin.T:
        prevph = seedphase
        startind -= 1#for FORTRAN compatability
        if stopind is None:
            stopind = cur_col.shape[0]
        elif stopind == 1 and step == -1:
            stopind = -1
        for i in range(startind, stopind, step):
            diff = cur_col[i] - prevph
            if (diff >= jump+360):
                cur_col[i] = cur_col[i]-720
            elif (diff >= jump):
                cur_col[i]=cur_col[i]-360
            elif (diff <= -jump-360):
                cur_col[i]=cur_col[i]+720
            elif (diff <= -jump):
                cur_col[i]=cur_col[i]+360
            prevph=cur_col[i]
    return squeeze(vectin)


try:
    from fortran_phase_engine import phaseengine
except ImportError:
    phaseengine = python_phase_engine

mt=time.time

#t3=time.time()
#t2=time.time()
#print('c import time ='+str(t2-t1))
#print('scipy + plot='+str(t3-t2))
#print('c_spec + inline_tools='+str(t1b-t1))
#print('t2-t1b='+str(t2-t1b))

def CreateDownSampleMask(vector, factor):
    tmask=[1]+[0]*factor
#    tmask.append(1)
    N=len(vector)/len(tmask)
    mask=tmask*N
    pad=len(vector)-len(mask)
    mypad=[0]*pad
    mask.extend(mypad)
    return mask

def DownSampleTimeDomain(vector, factor, axis=0):
    mask=CreateDownSampleMask(vector, factor)
    return vector.compress(mask, axis)


class Spectra:
   def __init__(self, input=None, output=None, \
                Gxx=None, Gyy=None, Gxy=None, ave=False):
      """Note that input and output are just labels."""
      self.input=input
      self.output=output
      self.Gxx=Gxx
      self.Gyy=Gyy
      self.Gxy=Gxy
      self.ave=ave

   def Average(self):
      keys=['Gxx','Gyy','Gxy']
      avedict={}
      for curkey in keys:
         curitem=getattr(self,curkey)
         avedict[curkey]=curitem.mean(1)
      return Spectra(self.input, self.output, ave=True, **avedict)

   def CutinHalf(self,ind):
      keys=['Gxx','Gyy','Gxy']
      lfdict={}
      hfdict={}
      for curkey in keys:
         curitem=getattr(self,curkey)
         lfdict[curkey]=copy.copy(curitem[0:ind])
         hfdict[curkey]=copy.copy(curitem[ind:])
      lfspec=Spectra(self.input, self.output, ave=self.ave, **lfdict)
      hfspec=Spectra(self.input, self.output, ave=self.ave, **hfdict)
      return lfspec, hfspec

   def DownSamplewithAveraging(self,dsf):
      keys=['Gxx','Gyy','Gxy']
      dsdict={}
      for curkey in keys:
         curitem=getattr(self,curkey)
         dsdict[curkey]=DownSamplewithAveraging(curitem,dsf)
      return Spectra(self.input, self.output, ave=self.ave, **dsdict)

   def Truncate(self, f, fcut):
      """Truncate self.Gxx, self.Gyy, and self.Gxy at the index corresponding to f[ind]=fcut.  return the truncated frequency vector."""
      fout=copy.copy(f)
      ind=thresh(f,fcut)
      fout=fout[0:ind]
      keys=['Gxx','Gyy','Gxy']
      for curkey in keys:
         curitem=colwise(getattr(self,curkey))
         curitem=curitem[0:ind,:]
         setattr(self,curkey,squeeze(curitem))
      return fout


def thresh_py(iterin, value, startind=0, above=1,reverse=0):
    myiter=iterin[startind:]
#    Pdb().set_trace()
    if reverse:
        myiter=myiter[::-1]
    ind=-1
    for x,item in enumerate(myiter):
        if above:
            if item>=value:
                ind=x
                break
        else:
            if item<=value:
                ind=x
                break
    ind=ind+startind
    if reverse:
        ind=len(iterin)-ind
    return ind

def _thresh(iterin, value, startind=0, above=1):
#    mt=time.time
#    t1=mt()
    if above:
#        keepinds=where(iterin>=value)[0]
        mybools=iterin>=value
    else:
#        keepinds=where(iterin<=value)[0]
        mybools=iterin<=value
#    t2=mt()
    myinds=arange(len(iterin))
#    t3=mt()
    keepinds=myinds[mybools]
#    t4=mt()
    inds2=keepinds>=startind
#    t5=mt()
    keep2=keepinds[inds2]
#    inds2=where(keepinds>=startind)[0]
#    t6=mt()
#    out=keepinds[inds2].min()
#    t7=mt()
#    print('In thresh')
#    print('total time='+str(t6-t1))
#    print('t2-t1='+str(t2-t1))
#    print('t3-t2='+str(t3-t2))
#    print('t4-t1='+str(t4-t1))
#    print('t5-t4='+str(t5-t4))
#    print('t6-t5='+str(t6-t5))
#    if inds2.any():
#        return keepinds[inds2].min()
#    else:
#        return -1
    return keep2

def edges(iterin, mythresh=2.5, startind=0, edges='rising'):
    inds = _thresh(iterin, mythresh, above=1)
    if edges=='rising' or edges=='both':
        risingdiff = inds[1:]-inds[0:-1]!=1
        risingedges = inds[risingdiff.nonzero()[0]+1]#diff.nonzero()[0] is an array of the nonzero elements of diff, similar to what would be returned by where.  The +1 is necessary because of the shifting that took place with diff = inds[1:]
        risingedges = risingedges[risingedges>=startind]
    if edges=='falling' or edges=='both':
        tc=mt()
        fallingdiff = inds[0:-1]-inds[1:]!=-1
        fallingedges = inds[fallingdiff.nonzero()[0]]
        fallingedges = fallingedges[fallingedges>=startind]
    if edges=='both':
        return risingedges, fallingedges
    elif edges=='rising':
        return risingedges
    elif edges=='falling':
        return fallingedges


def thresh(iterin, value, startind=0, above=1):
    if above:
#        keepinds=where(iterin>=value)[0]
        mybools=iterin>=value
    else:
#        keepinds=where(iterin<=value)[0]
        mybools=iterin<=value
    myinds=arange(len(iterin))
    keepinds=myinds[squeeze(mybools)]
    inds2=keepinds>=startind
    keep2=keepinds[inds2]
#    inds2=where(keepinds>=startind)[0]
    t6=mt()
#    out=keepinds[inds2].min()
#    if inds2.any():
#        return keepinds[inds2].min()
#    else:
#        return -1
    if len(keep2)==0:
        return -1
    else:
        return keep2.min()

## def thresh_old(iterin, value, startind=0, above=1):
## #        ind=-1
##     stopind=len(iterin)
## #    print('In thresh')
## #    print('type(value)='+str(type(value)))
## #    print('type(iterin[0])='+str(type(iterin[0])))
## #    print('type(startind)='+str(type(startind)))
## #    print('type(stopind)='+str(type(stopind)))
## #    print('type(above)='+str(type(above)))
##     code="""
##         #line 21 "rwkdataproc.py"
##         int ind;
##         for(int j=startind;j<stopind;j++){
##             if(above>0){
##                 if(iterin(j)>=value){
##                     ind=j;
##                     break;
##                 }
##             }
##             else{
##                 if(iterin(j)<=value){
##                     ind=j;
##                     break;
##                 }
##             }
##         }
##         return_val = ind;
##         """
##     N=inline_tools.inline(code,['iterin','value','above','startind','stopind'],
##             compiler='gcc',
##             type_converters = cblitz,
##             verbose = 1,force=1)
##     return N

def trunconesweep(rawfilename,sigdict,truncchannel,junklist=[],desL=-1,truncpath="trunc"):
    rawdata=loadmat(rawfilename)
    adict={}
    truncvar=rawdata.pop(sigdict[truncchannel])
    junklist.append(truncchannel)
    for curkey, curvalue in sigdict.iteritems():
        #exec(curkey+'_a=rawdata[\''+curvalue+'\']')
        if curkey not in junklist:
            adict[curkey]=rawdata[curvalue]
    startind=thresh_py(truncvar,0.5)
    if desL<=0:
        stopind=thresh_py(truncvar,0.5,startind+1,0)
    startind-=500
    if desL<=0:
        stopind+=500
    else:
        stopind=startind+desL
    trunc_dict={}
    for curkey, curchannel in adict.iteritems():
        trunc_dict[curkey]=adict[curkey][startind:stopind]
    trunc_dict['t']=trunc_dict['t']-trunc_dict['t'][0]
    outname,ext=os.path.splitext(rawfilename)
    outname=outname+'_trunc'+ext
#    print('pwd='+os.getcwd())
    if not os.path.exists(truncpath):
#        print('making dir '+truncpath)
        os.mkdir(truncpath)
    outpath=os.path.join(truncpath,outname)
#    print('saving to '+outpath)
    savemat(outpath,trunc_dict)
    return stopind-startind

class datastruct(dict):
    def __init__(self, **kwds):
        dict.__init__(self, **kwds)
        self.__dict__ = self

def BuildDataSet(pattern,directory='',ext='',xlabel='t'):
    extpattern=''
    if len(ext)>0 and ext[0]!='*':
        extpattern='*'+ext
    relpattern=os.path.join(directory,pattern)
    fullpattern=relpattern+extpattern
    temp=fullpattern.replace("**","*")
    while (temp!=fullpattern):
        fullpattern=temp
        temp=fullpattern.replace("**","*")
#    print('fullpattern='+fullpattern)
    filelist=glob.glob(fullpattern)

    mystruct=datastruct()
    mystruct.filenames=filelist
    mystruct.ynames=[]
    firstfile=1

    for filename in filelist:
#        print("filename="+filename)
        d=loadmat(filename)
        if firstfile==1:
            setattr(mystruct,xlabel,d[xlabel])
            mystruct.xname=xlabel
            for name, vector in d.iteritems():
                if name!=xlabel:
                    setattr(mystruct,name,c_[mat(vector).T])
                    mystruct.ynames.append(name)
            firstfile=0
        else:
            for name, vector in d.iteritems():
                if name!=xlabel:
                    tempmat=getattr(mystruct,name)
                    setattr(mystruct,name,c_[tempmat,mat(vector).T])
    return mystruct


#------------------------------------
#
#    Script Based Processing Again
#    05/09/2006
#
#------------------------------------
def PreviewOneFileForTrunc(filepath, signame, sigdict, truncsig='ce2', startcush=500, endcush=500,fi=1, fig=None):
    """Preview one untruncated signal.  signame is the truncname you
    wish to view, truncsig is the name of the channel used for
    truncation.  sigdict is a dictionary mapping long names to
    truncnames.  startcush and endcush are the cushions in data points
    to add to the data that is kept."""
    mydata=loadmat(filepath)
    trunckey=sigdict[truncsig]
    sigkey=sigdict[signame]
    truncsig=mydata[trunckey]
    signal=mydata[sigkey]
    t=mydata[sigdict['t']]
    ind1=thresh(truncsig,0.5)
    ind2=thresh(truncsig, 0.5, startind=ind1, above=0)
    if fig is None:
       from pylab import figure
       fig = figure(fi)
    ax = fig.gca()
    ax.cla()
    ax.plot(t,signal)
    i1=ind1-startcush
    i2=ind2+endcush
    ts=t[i1]
    te=t[i2]
    if fig is None:
       rwkmplutil.myvline(ts,'g:')
       rwkmplutil.myvline(te,'r:')
    return ind1,ind2

def TruneOneFile(filepath, signames, sigdict, truncsig='ce2', startcush=500, endcush=500, xlabel='t',folder='trunc'):
    from scipy.io import save_as_module
    mydata=loadmat(filepath)
    trunckey=sigdict[truncsig]
    truncsig=mydata[trunckey]
    t=mydata[sigdict[xlabel]]
    ind1=thresh(truncsig,0.5)
    ind2=thresh(truncsig, 0.5, startind=ind1, above=0)
    i1=ind1-startcush
    i2=ind2+endcush
    assert i1>0
    assert i2<max(shape(t))
    outdict={}
    outdict[xlabel]=t[i1:i2]
    for sig in signames:
       outdict[sig]=mydata[sigdict[sig]][i1:i2]
    infolder, name =os.path.split(filepath)
    fno,ext=os.path.splitext(name)
    outname=fno+'_trunc'
    #outname=fno+'_trunc'+ext
    outfolder=os.path.join(infolder,folder)
    if not os.path.exists(outfolder):
       os.mkdir(outfolder)
    outpath=os.path.join(outfolder, outname)
    #savemat(outpath,outdict)
    save_as_module(outpath,outdict)
    print('path='+outpath)
    return outdict



#=============================
#
#   Data Processing Scripts
#   intended for gui or script
#   based interface.
#   Started 10/5/05
#=============================

def mystack(listin):
    slist=[max(item.shape) for item in listin]
    myN=min(slist)
    newlist=[item[0:myN] for item in listin]
    return column_stack(newlist)


def ReadSigMap(filename):
    f = open(filename)
    mylines = f.readlines()
    f.close()
    rows = [item.split('\t') for item in mylines]
    col0 = [row[0].strip() for row in rows]
    col1 = [row[1].strip() for row in rows]
    return dict(zip(col1,col0))


def _initializedict(deschannels):
    outdict={}
    for sig in deschannels:
        outdict[sig]=[]
    return outdict

def LoadDatafromSavedModules(filelist,deschannels,xlabel='t'):
    """Builds matrices for each channel assuming that each filename
    has been saved to a module using scipy.io.save.  filelist is a
    list of the modules that will be imported using my_import, so that
    need to be somewhere on sys.path.  deschannels are the keys from
    the saved modules that will be kept."""
    outdict=_initializedict(deschannels)
    for modname in filelist:
        mod=my_import(modname)
        mydict=ConvertSaveModuletoDict(mod)
        for chn in deschannels:
           outdict[chn].append(mydict[chn])
    for key, value in outdict.iteritems():
        if key != xlabel:
#           print('key='+str(key))
#           print('type(value)='+str(type(value)))
#           print('shape(value)='+str(shape(value)))
           outdict[key]=mystack(value)
           myN=max(outdict[key].shape)
    outdict[xlabel]=mydict[xlabel][0:myN]
    return outdict

def LoadDatafromUntrunc(filelist,deschannels,namedict={},xlabel='t'):
    """Loads data from *.mat files and builds
    a datastructure with matrices of the signals.
    The files in filelist are loaded using
    loadmat(filename) without any attempt to
    search for the file.  So, filelist should either
    contain fullpaths or the files should be
    in the current directory.

    deschannels is a list of the truncnames you desire to keep.  It
    should not include 't'.

    namedict is a dictionary from a sigmap that has untruncated names
    as keys and truncnames as values."""
    outdict={}
    for sig in deschannels:
        outdict[sig]=[]
    d=loadmat(filelist[0])
    if namedict:
        xkey=namedict[xlabel]
    else:
        xkey=xlabel
    outdict[xlabel] = squeeze(d[xkey])
    for filename in filelist:
        d=loadmat(filename)
        for sig in deschannels:
            if namedict:
                curkey=namedict[sig]
            else:
                curkey=sig
            outdict[sig].append(squeeze(d[curkey]))
    for key, value in outdict.iteritems():
        if key != xlabel:
           outdict[key]=mystack(value)
           myN=max(outdict[key].shape)
    outdict[xlabel]=outdict[xlabel][0:myN]
    return outdict

def StringMatrix(matin, N):
    """reshape matrices whose columns contain individual tests to make
    tests of longer total time T.  Called by ConnectDataSets."""
    nr,nc=shape(matin)
    nr2=nr*nc/N
    nr2
    matout=matin.reshape(nr2,N)
    return matout

def MakeLongTime(t, N):
    """Create a long time vector to go along with stringing multiple
    tests together to form one longer one.  t is the time vector from
    an individual test.  N is the desired length of the new t vector."""
    dt=(t.max()-t.min())/(max(shape(t))-1)
    tout=arange(0,N,dtype='f')
    tout=tout*dt
    return tout


def ConnectDataSets(dictin, N, xlabel='t'):
    """String datasets from different tests into longer sets to
    increase total time T and reduce bias error.

    dictin is a dictionary whose keys are the truncnames for the
    signals and the values are matrices.  This is the dictionary that
    would be output by LoadDatafromUntrunc.

    N is the number of datasets to string together for each new,
    longer test."""
    outdict = {}
    for key, value in dictin.iteritems():
        if key!=xlabel:
            outdict[key]=StringMatrix(value,N)
    nro=shape(outdict[key])[0]
    outdict[xlabel]=MakeLongTime(dictin[xlabel],nro)
    return outdict

def BuildDataSetfromList(filelist, xlabel='t'):
    """Loads data from *.mat files and builds a datastructure with
    matrices of the signals.  The files in filelist are loaded using
    loadmat(filename) without any attempt to search for the file.  So,
    filelist should either contain fullpaths or the files should be in
    the current directory."""
    mydict={}
    mydict['filenames']=filelist
    mydict['ynames']=[]
    firstfile=1

    for filename in filelist:
#        print("filename="+filename)
        d=loadmat(filename)
        if firstfile==1:
            mydict[xlabel]=d[xlabel]
            mydict['xname']=xlabel
            for name, vector in d.iteritems():
                if name!=xlabel:
                    mydict[name]=c_[mat(vector).T]
                    mydict['ynames'].append(name)
            firstfile=0
        else:
            for name, vector in d.iteritems():
                if name!=xlabel:
                    tempmat=mydict[name]
                    mydict[name]=c_[tempmat,mat(vector).T]
    return mydict

def FilelistfromPattern(pattern,directory='',ext=''):
    extpattern=''
    if len(ext)>0 and ext[0]!='*':
        extpattern='*'+ext
    relpattern=os.path.join(directory,pattern)
    fullpattern=relpattern+extpattern
    temp=fullpattern.replace("**","*")
    while (temp!=fullpattern):
        fullpattern=temp
        temp=fullpattern.replace("**","*")
#    print('fullpattern='+fullpattern)
    filelist=glob.glob(fullpattern)
    return filelist

#def CalcBodeFromMatrix(time,outputmatrix,inputmatrix,truncfreq=100,avefile=True,jump=250):
#    tempfreq=makefreqvect(tddict['t'])
#    co=thresh_py(tempfreq,truncfreq)
#    freqvect=tempfreq[0:co]
#    outputmatrix=colwise(outputmatrix)
#    inputmatrix=colwise(inputmatrix)


def rwkfft(matin,cutoff=None,tvect=None):
    matin=colwise(matin)
    N=shape(matin)[0]
    fftmat=fft(matin,None,0)*2/N
    if cutoff is not None:
        tempf=makefreqvect(tvect)
        co=thresh_py(tempf,cutoff)
        fftmat=fftmat[0:co,:]
    return fftmat


def CalcSpectra(x,y, input=None, output=None):
    """Calculate Gxx, Gyy, and Gxy.  Note that input and output are
    just labels.  x and y are time domain signals."""
    N = max(shape(x))
    x_fft = squeeze(fft(x, None, 0)*2/N)
    y_fft = squeeze(fft(y, None, 0)*2/N)
    Gxx = norm2(x_fft)
    Gyy = norm2(y_fft)
    Gxy = (scipy.conj(x_fft))*y_fft
    return Spectra(input, output, Gxx, Gyy, Gxy)


def CalcBodesFromDataset(tddict, bodelist, description='', \
                         truncfreq=100, filto=2, wn=0.1,\
                         avefilt=True, jump=250):
#    pdb.set_trace()
    bodedict={}
    bodedict['description']=description

    tempfreq=makefreqvect(tddict['t'])
    co=thresh_py(tempfreq,truncfreq)
    bodedict['freq']=tempfreq[0:co]
    bodedict['xname']='freq'
    bodedict['bodes']=[]

    avedict=copy.deepcopy(bodedict)
    N=len(tddict['t'])

    fncopy=tddict['filenames']
    fnout=[]
    for fn in fncopy:
       folder,name=os.path.split(fn)
       fnout.append(name)
    bodedict['filenames']=fnout
    avedict['filenames']=fnout
    for i, curbode in enumerate(bodelist):
       outbode=copy.deepcopy(curbode)
       outbode.labels=bodedict['filenames']
       curavebode=copy.deepcopy(curbode)
       curavebode.averaged=True
       curavebode.labels=outbode.labels
       curin_t=tddict[curbode.input]
       curout_t=tddict[curbode.output]
       t1=time.time()
       curin_fft=fft(curin_t,None,0)*2/N
       curout_fft=fft(curout_t,None,0)*2/N
       t2=time.time()
       print('fft time='+str(t2-t1))
       curin_fft=curin_fft[0:co]
       curout_fft=curout_fft[0:co]
       curGxx=norm2(curin_fft)
       curGyy=norm2(curout_fft)
       curGxy=scipy.multiply(scipy.conj(curin_fft),curout_fft)
       H=scipy.divide(curGxy,curGxx)
       Hmag=abs(H)
       outbode.mag=Hmag
       Gxyave=scipy.mean(curGxy,1)
       Gxxave=scipy.mean(curGxx,1)
       Gyyave=scipy.mean(curGyy,1)
       cohnum=norm2(Gxyave)
       cohden=scipy.multiply(Gxxave,Gyyave)
       curavebode.coh=scipy.divide(cohnum,cohden)
       Have=scipy.divide(Gxyave,Gxxave)
       curavebode.mag=abs(Have)
       Hphase=mat_atan2(scipy.imag(H),scipy.real(H))
       Hphase_ave=mat_atan2(scipy.imag(Have),scipy.real(Have))
       Hphase=Hphase*180/scipy.pi
       Hphase_ave=Hphase_ave*180/scipy.pi
       if curbode.seedfreq and curbode.seedphase:
#           print('Massaging the phase')
            Hphase=PhaseMassage(Hphase,bodedict['freq'],curbode.seedfreq,curbode.seedphase)
            Hphase_ave=PhaseMassage(Hphase_ave,avedict['freq'],curbode.seedfreq,curbode.seedphase)
            if avefilt:
                Hphase=PhaseMassageFilter(Hphase,filto,wn,jump)
                Hphase_ave=PhaseMassageFilter(Hphase_ave,filto,wn,jump)
#           Hphase=AveMessage(Hphase,bodedict['freq'],curbode.seedfreq,curbode.seedphase,280.0,5)
#           Hphase_ave=AveMessage(Hphase_ave,avedict['freq'],curbode.seedfreq,curbode.seedphase,280.0,5)
       outbode.phase=Hphase
       curavebode.phase=Hphase_ave
       bodedict['bodes'].append(outbode)
       avedict['bodes'].append(curavebode)
    return bodedict,avedict

#    if len(datapath)>0:
#        os.chdir(datapath)

#    scipy.io.save('trunc/bodedict_'+endstr,bodedict)
#    scipy.io.save('trunc/avedict_'+endstr,avedict)

#    os.chdir(curdir)

def ConvertSaveModuletoDict(modulein):
#    mydict=vars(modulein)
    mydict=modulein.__dict__
    dictout={}
    for key,item in mydict.iteritems():
        if key.find('__')!=0 and type(item)!=types.ModuleType:
            dictout[key]=item
    return dictout

#=============================
#
#   End new scripts
#
#=============================

#----------------------------------
#Other newish stuff 05/09/06
#----------------------------------

def DownSamplewithAveraging(vectin,factor):
    a = shape(vectin)
    assert len(a)==1
    N = a [0]
    newN = N/factor
    f2 = vectin.reshape(newN,factor)
    f3 = f2.mean(1)
    return f3

def InterpVector(dsvect,fds,vect,f):
    dsf = len(f)/len(fds)
    myint = interpolate.interp1d(fds,squeeze(dsvect))
    n = int((dsf+0.5)/2.0)
    fb = f[n:-n]
    intv = myint(fb)
    return fb, intv

def AveSteps(dsvect, dsf):
    mymat=array([dsvect]*dsf)
    myvect=mymat.transpose().flatten()
    return myvect

def CalcBiasEror(dsvect,fds,vect,f,zeroorder=True):
    dsf = len(f)/len(fds)
    if zeroorder:
       stepvect=AveSteps(dsvect,dsf)
       evect=squeeze(vect-stepvect)
       return f, evect
    else:
       n = int((dsf+0.5)/2.0)
       fb, intv = InterpVector(dsvect,fds,vect,f)
       evect=squeeze(vect[n:-n])-squeeze(intv)
       return fb, evect/squeeze(vect[n:-n])

#-----------------------------------

def PreviewData(datastructin, fig=None):
    x=getattr(datastructin,datastructin.xname)
    usepylab = False
    if fig is None:
       from pylab import figure
       fig = figure(fi)
       usepylab = True
    ax = fig.gca()


    for key, curset in datastructin.datasets.iteritems():
        for channel in curset.ynames:
            ax.cla()
            if usepylab:
               ioff()
            curmat=getattr(curset,channel)
#            pdb.set_trace()
            for i in range(shape(curmat)[1]):
                ax.plot(x,curmat[:,i])
            ax.set_ylabel(channel)
            ax.set_xlabel(datastructin.xname)
            ax.set_title(str(key))
            ax.legend(curset.filenames)
            if usepylab:
               ion()
               show()
            raw_input('Please press return to continue...\n')

def makefreqvect(timevect):
    tspan=timevect.max()-timevect.min()
    N = timevect.shape[0]
    dt=tspan/(N-1)
    fs=1.0/dt
    T = tspan+dt
    df=1.0/T
    f=arange(0,fs,df)
    if max(shape(f))==(max(shape(timevect))+1):
        f=f[0:-1]
    return f

def FindFixedSineInputFreq(timevector,inputsignal,minfreq=0.05,maxfreq=100):
    """This function finds the input frequency
    for fixed sine data.  The timevector is used
    to calculate the frequency vector.  The return
    value is the frequency,index of the maximum
    magnitude of inputsignal between minfreq and
    maxfreq (these are required to avoid problems
    with DC values and aliasing)."""
    f=makefreqvect(timevector)
    fminind=thresh_py(f,minfreq)
    fmaxind=thresh_py(f,maxfreq)
    N=len(timevector)
    inputfft=fft(inputsignal,None,0)*2/N
    inputmag=abs(inputfft)
    inputslice=inputmag[fminind:fmaxind]
    maxtemp=argmax(inputslice)
    maxind=fminind+maxtemp
    fatmax=f[maxind]
    return fatmax,maxind

def BodefromTwoTimeDomainVectors(timevector,output,input,truncfreq=100):
    """This function calculates the Bode response between two time
    domain signals.  The timevector is used to calculate the frequency
    vector, which is then used to truncate the Bode response to reduce
    calculation time and return only useful information.  Input and
    output are time domain vectors.

    The return values are
    freq, magnitude ratio, phase, complex

    The goal of this function is to be useful for small amounts of
    data and/or as part of a routine to calculate a Bode response from
    fixed sine data."""

    N=len(timevector)
    f=makefreqvect(timevector)
    co=thresh_py(f,truncfreq)
    f=f[0:co]
    curin_fft=fft(input,None,0)*2/N
    curout_fft=fft(output,None,0)*2/N
    curin_fft=curin_fft[0:co]
    curout_fft=curout_fft[0:co]
    curGxx=norm2(curin_fft)
    curGyy=norm2(curout_fft)
    curGxy=scipy.multiply(scipy.conj(curin_fft),curout_fft)
    H=scipy.divide(curGxy,curGxx)
    Hmag=abs(H)
    Hphase=mat_atan2(scipy.imag(H),scipy.real(H))*180.0/pi
    return f,Hmag,Hphase,H

def pherrcomp(e,jump=250):
    out=0
    if e>abs(jump):
        out=-360
    elif e<-abs(jump):
        out=360
    return out

def PhaseMassageFilter(phin, N=2, Wn=0.1, jump=250,freqvect=[],chkind=0, fig=None):
    if fig is None:
       from pylab import figure
       fig = figure(fi)

    phun=colwise(phin)
    phout=zeros(shape(phun),'d')#+0.0j
    (b,a)=scipy.signal.butter(N,Wn)
    phfilt=scipy.signal.lfilter(b,a,phun,axis=0)
    for i in range(shape(phout)[1]):
        pherr=phun[:,i]-phfilt[:,i]
        phcor=map(pherrcomp,pherr)
        phout[:,i]=phun[:,i]+phcor
        if freqvect and i==chkind:
            fig.clear()
            ax1 = fig.add_subplot(2, 1, 1)
            ax1.semilogx(freqvect,phun[:,chkind])
            ax1.semilogx(freqvect,phfilt[:,chkind])
            ax1.semilogx(freqvect,phout[:,chkind])
            ax1.legend(['un','filt','out'],3)
            ax2 = fig.add_subplot(2, 1, 2)
            ax2.semilogx(freqvect,pherr)
            ax2.semilogx(freqvect,phcor)
            ax2.legend(['err','cor'],3)
    return phout

mt=time.time

def PhaseMassage(phin, freqin, seedfreq, seedphase,jump=250.0):
    if type(freqin) != ndarray:
        freqin = array(freqin)
    if type(phin) != ndarray:
        phin = array(phin)
    ind=thresh(freqin,seedfreq)
    vect1=phin[ind:]
    vect1=phaseengine(vect1,seedphase,jump=jump)
    vect2=phin[0:ind]
    vect2=phaseengine(vect2,seedphase,ind,1,-1,jump=jump)
    return r_[vect2,vect1]


## def PhaseMassage_old(phin, freqin, seedfreq, seedphase,jump=250.0):
## #    print('calling PhaseMassage')
## #    phout=copy.deepcopy(phin)
## #    ts=mt()
##     phout=colwise(phin)
##     addedcol=False
## #    ta=mt()
##     if len(shape(phout))==1:
##         phout=phout[:,NewAxis]
##         addedcol=True
##     if not hasattr(freqin, 'min'):
##         freqin=array(freqin)
##     if seedfreq < freqin.min():
##         ind=0
##     else:
##         #tb=mt()
##         #ind=thresh_py(freqin,seedfreq)
## #        tb=mt()
##         inds=where(freqin>seedfreq)
##         if type(inds)==type(('a',5)):#is where returning a tuple
##             if len(inds)==1:
##                 inds=inds[0]
##         ind=inds.min()
## #        tc=mt()
## #        print('thresh time='+str(tc-tb))
## #    td=mt()
##     ph1=phout[0:ind,:]
##     ph2=phout[ind:,:]
## #    t2=mt()
##     Pdb().set_trace()
##     ph2o=phaseengine(ph2,seedphase,1,jump)
## #    t3=mt()
##     if ind==0:
##         ph1o=[]
##         outmat=ph2o
##     else:
## #        t4=mt()
##         ph1o=phaseengine(ph1,seedphase,-1,jump)
## #        t5=mt()
##         outmat=r_[ph1o,ph2o]
## #        t6=mt()
## #    mylist=rwkascii.rwktextlist()
## #    mylist.list=[]
## #    mylist.filename='phaseenginelog_PM.txt'
## #    for r in range(shape(outmat)[0]):
## #        curline=str(r)
## #        for c in range(shape(outmat)[1]):
## #            curline+='\t'+str(outmat[r,c])
## #        mylist.append(curline)
## #    mylist.tofile(0)
##     if addedcol:
##         outmat=outmat[:,0]
## #    te=mt()
## #    tlist=[ts,t2,t3,t4,t5,t6,te]
## #    inds=range(1,len(tlist))
## #    tdiffs=[tlist[x]-tlist[x-1] for x in inds]
## #    outlist=[str(x)+':'+str(item) for x, item in enumerate(tdiffs)]
## #    outstr='\n'.join(outlist)
## #    print('-------------------')
## #    print(outstr)
## #    print('colwise time='+str(ta-ts))
## #    print('submat time='+str(t2-td))
## #    print('td-ta='+str(td-ta))
## #    print('total PhaseMassage time='+str(te-ts))
## #    print('-------------------')
##     return outmat

## def phaseengine_old(phmat,seedphase,step=1,jump=250.0):
## #    mylist=rwkascii.rwktextlist()
## #    mylist.list=[]
## #    mylist.filename='phaseenginelog.txt'
## #    print('jump='+str(jump))
##     t1=time.time()
## #    phout=copy.deepcopy(phmat)
##     phout=phmat
##     N=shape(phmat)[0]
##     nc=shape(phmat)[1]
## #    print('type (nc)='+str(type(nc)))
##     t2=time.time()
##     if step==1:
##         startind=0;
##         stopind=shape(phmat)[0]-1
##     elif step==-1:
##         startind=shape(phmat)[0]-1
##         stopind=0
## #    print('phout[startind,:]='+str(phout[startind,:]))
## #    pdb.set_trace()
##     for q in range(nc):
## #        print('phout[startind,q]='+str(phout[startind,q]))
## #        print('seedphase='+str(seedphase))
## #        print('jump='+str(jump))
## #        pdb.set_trace()
##         if (phout[startind,q]-seedphase)>jump:
##             phout[startind,q]-=360
##         elif (seedphase-phout[startind,q])>jump:
##             phout[startind,q]+=360
## #    print('phout[startind,:]='+str(phout[startind,:]))
## #----------------------------------
## #   Python Only implementation
## #----------------------------------
##     usepy=0
##     if usepy:
##         rns=range(shape(phout)[0])
##         if step==-1:
##             rns=rwkmisc.reverse(rns)
##         prevr=rns.pop(0)
## #        t3=time.time()
##     #    mylist.append('r\tc0b\tc0p\tc0a\tc1b\tc1p\tc1a\tc2b\tc2p\tc2a')
##         for r in rns:
##     #        curline=str(r)
##             for q in range(nc):
##     #            curline+='\t'+str(phout[r,q])+'\t'+str(phout[prevr,q])
##                 if (phout[r,q]-phout[prevr,q])>jump:
##                     phout[r,q]-=360
##                 elif (phout[prevr,q]-phout[r,q])>jump:
##                     phout[r,q]+=360
##     #            curline+='\t'+str(phout[r,q])
##     #        mylist.append(curline)
##             prevr=r
##     #    mylist.tofile(0)
## ###----------------------------------
## ###ccccccccccccccccccccccccccccccccc
## ###   C code implementation
## ###ccccccccccccccccccccccccccccccccc
##     else:
## #        t3=time.time()
##         startind+=step#this is for the c-code only
##         cursum=zeros((1,nc),'f')
##         curave=zeros((1,nc),'f')
##         aveout=zeros(shape(phout),'f')
##         code="""
##             #line 215 "rwkdataproc.py"
##             int previ;
##             previ=startind-step;
##             //std::cout<<"using C imp."<<std::endl;
##             //std::cout<<"I am not averaging"<<std::endl;
##             for ( int i=startind;step*i<=stopind;i+=step )
##             {
##                 //std::cout<<"i="<<i<<std::endl;
##                 for ( int m=0;m<nc;m++ )
##                 {
##                     //std::cout<<"m="<<m<<std::endl;
##                     if ((phout(i,m)-phout(previ,m))>jump){
##                         phout(i,m)-=360.0;
##                     }
##                     else if ((phout(previ,m)-phout(i,m))>jump){
##                         phout(i,m)+=360.0;
##                     }
##                 }
##                 previ=i;
##             }
##             """
##         inline_tools.inline(code,['phout','jump','step','startind','stopind','nc'],
##                 compiler='gcc',
##                 type_converters = cblitz,
##                 verbose = 1, force=1)
## ###cccccccccccccccccccccccccccccccccccccc
## #    t4=time.time()
## #    print('copy time='+str(t2-t1))
## #    print('time to start of loop='+str(t3-t1))
## #    print('weave loop time='+str(t4-t3))
## #    print('total time='+str(t4-t1))
## #    for r in range(shape(phout)[0]):
## #        curline=str(r)
## #        for c in range(shape(phout)[1]):
## #            curline+='\t'+str(phout[r,c])
## #        mylist.append(curline)
## #    mylist.tofile(0)
##     return phout

## def phaseengine_w_ave(phmat,step=1,nave=5.0,jumpave=120.0):
##     #Note that this function only does averaging, it does not use a seedphase or look for jumps from one element to the next.  It only looks at the current element compared to the average of the nave previous elements.  It is expected that this function is called only after first calling the regular phaseengine to first clean up the data.
## #    print('jumpave='+str(jumpave))
##     t1=time.time()
## #    phout=copy.deepcopy(phmat)
##     phout=phmat
##     N=shape(phmat)[0]
##     if len(shape(phmat))==1:
##         nc=1
##     else:
##         nc=shape(phmat)[1]#assume vector for now
##     t2=time.time()
##     if step==1:
##         startind=0;
##         stopind=shape(phmat)[0]-1
##     elif step==-1:
##         startind=shape(phmat)[0]-1
##         stopind=0
## ###ccccccccccccccccccccccccccccccccc
## ###   C code implementation
## ###ccccccccccccccccccccccccccccccccc
##     t3=time.time()
##     aveout=zeros(shape(phout),'f')
##     startind+=step#this is for the c-code only
##     code="""
##         #line 314 "rwkdataproc.py"
##         //std::ofstream outfile("log.txt");
##         int previ;
##         int rind;
##         int q;
##         double cursum;
##         double curave;
##         bool remove;
##         double numin;
##         //outfile<<"before loop:"<<std::endl;
##         //outfile<<"previ="<<previ;
##         //outfile<<" startind="<<startind<<std::endl;
##         //std::cout<<"Nphout[1]="<<Nphout[1]<<std::endl;
##         //std::cout<<"phout(0,0)="<<phout(0,0)<<std::endl;
##         for ( int j=0;j<nc;j++ ){
##             previ=startind-step;
##             cursum=0.0;
##             curave=0.0;
##             numin=0.0;
##             remove=false;
##             for ( int i=startind;step*i<=stopind;i+=step )
##             {
##                 //std::cout<<"i="<<i<<std::endl;
##                 if ( numin < nave )
##                 {
##                     numin+=1.0;
##                     //outfile<<"i="<<i<<" numin="<<numin<<std::endl;
##                 }
##                 else
##                 {
##                     remove=true;
##                 }
##                 cursum+=phout(previ,j);
##                 if ( remove )
##                 {
##                     //outfile<<"i="<<i<<" stopind="<<stopind<<" previ="<<previ<<" nave="<<nave<<" step="<<step;
##                     rind=static_cast<int>(previ-nave*step);
##                     //outfile<<" rind="<<rind;
##                     cursum-=phout(rind);
##                 }
##                 //outfile<<" numin="<<numin<<" cursum="<<cursum;
##                 curave=cursum/numin;
##                 aveout(previ,j)=curave;
##                 if ( (phout(i,j)-aveout(previ,j))>jumpave )
##                 {
##                     phout(i,j)-=360.0;
##                 }
##                 else if ( (aveout(previ,j)-phout(i,j))>jumpave )
##                 {
##                     phout(i,j)+=360.0;
##                 }
##                 //aaaaaaaaaaaaaaaaaaaaaaaa
##                 //outfile<<" curave="<<curave<<" phout(i,j)="<<phout(i,j)<<std::endl;
##                 previ=i;
##             }
##         //outfile<<"Exiting"<<std::endl;
##         //outfile.close();
##         }
##         """
##     inline_tools.inline(code,['phout','jumpave','step','startind','stopind','nave','aveout','nc'],
##             compiler='gcc',
##             type_converters = cblitz,
##             verbose = 1,
##             headers=["<iostream>","<fstream>"])
##     t4=time.time()
## #    print('copy time='+str(t2-t1))
## #    print('time to start of loop='+str(t3-t1))
## #    print('loop time='+str(t4-t3))
## #    print('total time='+str(t4-t1))
## #    for r in range(shape(phout)[0]):
## #        curline=str(r)
## #        for c in range(shape(phout)[1]):
## #            curline+='\t'+str(phout[r,c])
## #        mylist.append(curline)
## #    mylist.tofile(0)
##     return phout

## def AveMessage(phin, freqin, seedfreq, seedphase,jump=250.0,nave=5.0):
##     #Note: Run PhaseMassage first and then pass the data to AveMessage as a second clean up step
##     phout=copy.deepcopy(phin)
##     ind=thresh_py(freqin,seedfreq)
## #    print('ind='+str(ind))
## #    print('freq(ind)='+str(freqin[ind]))
##     ph1=phout[0:ind]
##     ph2=phout[ind:]
## #    print('shape(ph1)='+str(shape(ph1)))
## #    print('shape(ph2)='+str(shape(ph2)))
## #    pdb.set_trace()
##     ph1o=zeros(shape(ph1),'f')
##     ph2o=zeros(shape(ph2),'f')
##     ph1o=phaseengine_w_ave(ph1,-1,nave,jump)
##     ph2o=phaseengine_w_ave(ph2,1,nave,jump)
##     outmat=r_[ph1o,ph2o]
## #    print('type(ph1o)='+str(type(ph1o)))
## #    outmat=ph1o
## #    print('type(ph2o)='+str(type(ph2o)))
## #    outmat=ph2o
## #    print('last line of AveMessage')
##     return outmat

def norm2(x):
    """compute |x|^2 = x*conjugate(x)"""
    if iscomplexobj(x):
#        t1=time.time()
#        mat1=x.real**2 + x.imag**2
#        t2=time.time()
        mat2=multiply(x.real,x.real) + multiply(x.imag,x.imag)
#        t3=time.time()
#        print('---------------------------')
#        print('pow time='+str(t2-t1))
#        print('multiply time='+str(t3-t2))
#        print('pow time/multiply time='+str((t2-t1)/(t3-t2)))
#        print('shape(x)='+str(shape(x)))
#        print('x.typecode='+str(x.typecode()))
#        print('mat1.typecode='+str(mat1.typecode()))
#        print('mat2.typecode='+str(mat2.typecode()))
#        if len(shape(x))==1:
#            print('type(x[0])='+str(type(x[0])))
#            print('type(mat1[0])='+str(type(mat1[0])))
#            print('type(mat2[0])='+str(type(mat2[0])))
#        else:
#            print('type(x[0,0])='+str(type(x[0,0])))
#            print('type(mat1[0,0])='+str(type(mat1[0,0])))
#            print('type(mat2[0,0])='+str(type(mat2[0,0])))
        return mat2
    else:
        return multiply(x,x)


def mat_atan2(y,x):
    if shape(y)!=shape(x):
        raise IndexError, 'In matrix atan2, shape(x)!-shape(y)'
    x=atleast_2d(x)
    y=atleast_2d(y)
    nr=shape(y)[0]
    nc=shape(y)[1]
    outmat=zeros((nr,nc),'d')
#    pdb.set_trace()
    for i in range(nr):
        for j in range(nc):
            outmat[i,j]=atan2(y[i,j],x[i,j])
    return colwise(outmat)

## def mat_atan2_c(y,x):
## #    print('in mat_atan2_c')
##     assert(shape(y)==shape(x))
##     outmat=zeros(shape(y),'d')
##     nr=shape(y)[0]
##     if len(shape(y))==1:
##         nc=1
##     else:
##         nc=shape(y)[1]
##     code = """
##            for(int i = 0; i < nr; i++)
##                for(int j = 0; j < nc; j++)
##                    outmat(i,j) = atan2(y(i,j),x(i,j));
##            """
##     inline_tools.inline(code,['outmat','y','x','nr','nc'],
##                         type_converters = cblitz,
##                         compiler='gcc',
##                         verbose = 1)
##     return outmat

