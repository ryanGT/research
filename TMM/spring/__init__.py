from __future__ import division
from scipy import cos, cosh, sin, sinh, array, r_, c_, exp, pi, dot, real, imag, zeros, eye, shape, atleast_1d, shape
#import MLab
import scipy
from scipy.linalg import det
import copy
import pdb

#import TMM
#reload(TMM)
#from TMM.TMMElement import TMMElement, HT4
from TMM.TMMElement import TMMElementIHT
#import TMM
#reload(TMM)
#from TMM import HT4

from rwkmisc import symstr, SymstrMattoMaxima
import rwkmisc

from IPython.core.debugger import Pdb
#TMMElement=TMM.TMMElement

class TorsionalSpringDamper(TMMElementIHT):
    """This class represents a spring element with compliance along 1,
    2, or 3 axes (for maxsize = 4, 8, or 12)."""
    def __init__(self,params, symlabel='sd', symname='Usd',**kwargs):
        """Create a TorsionalSpringDamper instance.  params is a
        dictionary with keys of 'k' and 'c'.  'k' and 'c' are lists or
        arrays if maxsize > 4.  If maxsize==8, k = [k0,k1].  If
        maxsize==12, k=[k0, k1, k2].  'c' is optional and will be set
        to zeros(shape(k)) if it is omitted.  Otherwise c=c0, [c0],
        [c0,c1], or [c0, c1, c2]."""
        if not params.has_key('k'):
            params['k']=None
        params['k']=atleast_1d(params['k'])
        if not params.has_key('c'):
            params['c']=zeros(shape(params['k']))
        else:
            params['c']=atleast_1d(params['c'])
        TMMElementIHT.__init__(self,'spring',params,symlabel=symlabel,symname=symname,**kwargs)

    def GetMat(self,s,sym=False):
        """Return the element transfer matrix for the
        TorsionalSpringDamper element.  If sym=True, 's' must be a
        symbolic string and a matrix of strings will be returned.
        Otherwise, 's' is a numeric value (probably complex) and the
        matrix returned will be complex."""
        N=self.maxsize
        if sym:
            myparams=self.symparams
        else:
            myparams=self.params
        k=myparams['k']
        c=myparams['c']
        springterm=1/(k[0]+c[0]*s)
        if sym:
            maxlen=len(springterm)+10
            matout=eye(N,dtype='f')
            matout=matout.astype('S%d'%maxlen)
        else:
            matout=eye(N,dtype='D')
        matout[1,2]=springterm
        if max(shape(k))>1 and self.maxsize>=8:
            matout[5,6]=1/(k[1]+c[1]*s)
        if max(shape(k))>2 and self.maxsize>=12:
            matout[9,10]=1/(k[2]+c[2]*s)
        return matout

    def GetMaximaString(self,name=None,label=None,N=None,aug=0):#omit
        defs=[]
        params=[]
        if label is None:
            label=self.symlabel
##             if self.symlabel:
##                 label=self.symlabel
##             else:
##                 label='sd'
##                 self.symlabel=label
        if name is None:
            name=self.symname
##             if self.symname:
##                 name=self.symname
##             else:
##                 name='Usd'
        maximalines=[]
        if N is None:
            N=self.maxsize
        Nstr=str(N)
        maximalines.append(name+':ident('+Nstr+')')
        dof=int(N/4)
        for x in range(dof):
            if dof==1:
                rc=''
            else:
                rc=str(x+1)
            r=4*x+2
            c=r+1
            kstr='k'+label+rc
            cstr='c'+label+rc
            maximalines.append(name+'['+str(r)+','+str(c)+']:1/('+kstr+'+s*'+cstr+')')
            params.append(kstr)
            params.append(cstr)
        for line in maximalines:
            if line[-1]!=';':
                line+=';'
        return maximalines, name, defs, params, label

    def GetMaximaLatexString(self,name=None,label=None,wrap=0,N=None,intro=None,aug=0):#omit
        if intro is None:
            intro='The transfer matrix for a spring/damper element is given by'
        return TMMElementIHT.GetMaximaLatexString(self,name=name,label=label,wrap=wrap,N=N,intro=intro,aug=aug)

#    def GetHT(self):
#        return HT4()

class SpringToGround(TMMElementIHT):
    def __init__(self,params,maxsize=4,label='',):
        TMMElementIHT.__init__(self,'spring to ground',params,maxsize=maxsize,label=label)

    def GetMat(self,s,sym=False):
        """Return the element transfer matrix for the
        TorsionalSpringDamper element.  If sym=True, 's' must be a
        symbolic string and a matrix of strings will be returned.
        Otherwise, 's' is a numeric value (probably complex) and the
        matrix returned will be complex."""
        N=self.maxsize
        if sym:
            myparams=self.symparams
        else:
            myparams=self.params
        k=myparams['k']
        if sym:
            maxlen=len(k)+10
            matout=eye(N,dtype='f')
            matout=matout.astype('S%d'%maxlen)
        else:
            matout=eye(N,dtype='D')
        matout[N-1,0]=k
        return matout


class MajetteActuator(TorsionalSpringDamper):#omit me
    def __init__(self,params,maxsize=12,label='',):
        params['k']=atleast_1d(params['k'])
        if not params.has_key('c'):
            params['c']=zeros(shape(params['k']))
        else:
            params['c']=atleast_1d(params['c'])
        TMMElementIHT.__init__(self,'spring with actuator',params,maxsize=maxsize,label=label)

    def GetAugMat(self,s):
        N=self.maxsize
        tempout=eye(N+1,dtype='D')
        tempout[0:N,0:N]=self.GetMat(s)
        tempout[2,4]=1.0
        return tempout

#    def GetHT(self):
#        return HT4()


class Spring2x2(TMMElementIHT):
    """This class represents a spring element in a system where only
    linear displacement and force are used as states.  This would be
    used with mass-spring-damper simplified systems like a system of
    carts used in a lot of introductory vibrations or controls
    classes."""
    def __init__(self,params,symlabel='K',**kwargs):
        """Create an instance of the Spring2x2 class.  params is a
        dictionary for consistancy with other TMMElement's, but it has
        only one entry: 'k'."""
        TMMElementIHT.__init__(self,'spring',params,maxsize=2,symlabel=symlabel,**kwargs)

    def GetMat(self,s,sym=False):
        """Return the transfer matrix for the Spring2x2 TMMElement,
        given s as an input.  If sym=True, the transfer matrix will be
        a string."""
        N=self.maxsize
        if sym:
            myparams=self.symparams
        else:
            myparams=self.params
        k=myparams['k']
        maxlen=20
        if sym:#type(s)==type(rwkmisc.symstr('s')):
            matout=eye(N,dtype='f')
            matout=matout.astype('S%d'%maxlen)
        else:
            matout=eye(N,dtype='D')
        matout[0,1]=1/k
        return matout
