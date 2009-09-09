from __future__ import division
from scipy import cos, cosh, sin, sinh, array, r_, c_, exp, pi, dot, real, imag, zeros, eye, shape, atleast_1d
import MLab
import scipy
from scipy.linalg import det
import copy
import pdb

from TMM import TMMElement, HT4

class ForcingElement(TMMElement):
    def __init__(self,params,maxsize=12,label=''):
        TMMElement.__init__(self,'forcing',params,maxsize=maxsize,label=label)

    def GetMat(self,s):
        return eye(self.maxsize,'d')+0.0j

    def GetAugMat(self,s):
        N=self.maxsize
        if N%2:
            N-=1
#        matout=eye(N+1,'d')+0j
        matout=eye(N+1,dtype='D')
        fv=self.params['fv']
        matout[0:N,N]=fv
        return matout

    def GetMaximaString(self,name=None,label=None,N=None,aug=0):
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
        return maximalines, name, [], []

    def GetAugMaximaString(self,name=None,label=None,N=None):
#        pdb.set_trace()
        maxlines,myname,defs,params=self.GetMaximaString(name,label,N)
        if N is None:
            N=self.maxsize
        Nstr=str(N)
        Np1str=str(N+1)
        fv=self.params['fv']
        maxlines.append('newcol:zeromatrix('+Nstr+',1)')
        for x,val in enumerate(fv):
            if val:
                maxlines.append('newcol['+str(x+1)+',1]:'+str(val))
        maxlines.append('newrow:zeromatrix(1,'+Np1str+')')
        maxlines.append('newrow[1,'+Np1str+']:1')
        maxlines.append(myname+':addcol('+myname+',newcol)')
        maxlines.append(myname+':addrow('+myname+',newrow)')
        return maxlines, myname,defs,params

    def GetHT(self):
        return HT4()

