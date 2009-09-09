from __future__ import division
from scipy import cos, cosh, sin, sinh, array, r_, c_, exp, pi, dot, real, imag, zeros, eye, shape, atleast_1d
#import MLab
import scipy
from scipy.linalg import det
import copy
import pdb

from TMM.TMMElement import TMMElement, HT4, TMMElementLHT

class ForcingElement(TMMElementLHT):
    """This class represents arbitrary forcing at any location in a
    TMM model, using a vector as a column of the augmented transfer
    matrix."""
    def __init__(self,params,**kwargs):
        """Initialize the ForcingElement.  params is a dictionary with
        only one key: 'fv', whose value is a column vector of the
        elements of the augmented column of the transfer matrix.  It
        should have only N element (the bottom 1 is assumed)."""
        TMMElementLHT.__init__(self,'forcing',params,**kwargs)
        for n, ent  in enumerate(self.params['fv']):
            if ent == 0:
                self.symparams['fv'][n]='0'
        filtfv=[item for item in self.params['fv'] if item]
        if len(filtfv)==1:
            #there is only 1 input and symparams should just have a bunch of zeros and one 1
            self.symparams['fv']=[str(item) for item in self.params['fv']]

    def GetMat(self, s, sym=False):
        """Return the element transfer matrix for the
        ForcingElement element.  If sym=True, 's' must be a
        symbolic string and a matrix of strings will be returned.
        Otherwise, 's' is a numeric value (probably complex) and the
        matrix returned will be complex.  Note that the
        ForcingElement element will return an identity matrix
        for its non-augmented transfer matrix."""
        N=self.maxsize
        if sym:
            matout=eye(N,dtype='f')
            matout=matout.astype('S1')
        else:
            matout=eye(N,dtype='D')
        return matout

    def GetAugMat(self, s, sym=False):
        """Return the augmented element transfer matrix for the
        ForcingElement element.  If sym=True, 's' must be a
        symbolic string and a matrix of strings will be returned.
        Otherwise, 's' is a numeric value (probably complex) and the
        matrix returned will be complex."""
        N=self.maxsize
        if N%2:
            N-=1
#        matout=eye(N+1,'d')+0j
        #matout=eye(N+1,dtype='D')
        if sym:
            matout = eye(N+1,dtype='f')
            matout = matout.astype('S20')
            myparams = self.symparams
        else:
            matout = eye(N+1,dtype='D')
            myparams = self.params
#        pdb.set_trace()
        fv = myparams['fv']
        matout[0:N,N]=fv
        return matout

    def GetMaximaString(self,name=None,label=None,N=None,aug=0):#omit
        if label is None:
            if self.symlabel:
                label=self.symlabel
            else:
                label='fv'
        if name is None:
            if self.symname:
                name=self.symname
            else:
                name='Uforcing'
        maximalines=[]
        if N is None:
            N=self.maxsize
        Nstr=str(N)
        maximalines=[]
        maximalines.append(name+':ident('+Nstr+')')
        return maximalines, name, [], [], label

    def GetAugMaximaString(self,name=None,label=None,N=None):#omit
        maxlines,myname,defs,params,label=self.GetMaximaString(name,label,N)
        if N is None:
            N=self.maxsize
        Nstr=str(N)
        Np1str=str(N+1)
        fv=self.params['fv']
        maxlines.append('newcol:zeromatrix('+Nstr+',1)')
        for x,val in enumerate(fv):
            if val:
                maxlines.append('newcol['+str(x+1)+',1]:'+str(real(val)))
        maxlines.append('newrow:zeromatrix(1,'+Np1str+')')
        maxlines.append('newrow[1,'+Np1str+']:1')
        maxlines.append(myname+':addcol('+myname+',newcol)')
        maxlines.append(myname+':addrow('+myname+',newrow)')
        return maxlines, myname,defs,params, label

    def GetMaximaLatexString(self,name=None,label=None,wrap=0,N=None,intro=None,aug=0):#omit
        if intro is None:
            intro='The transfer matrix for a forcing element is given by'
        return TMMElementLHT.GetMaximaLatexString(self,name=name,label=label,wrap=wrap,N=N,intro=intro,aug=aug)

#    def GetHT(self):
#        return HT4()

