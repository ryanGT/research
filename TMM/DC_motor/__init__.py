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

def Gc_PD(s, params):
    kp = params['kp']
    kd = params['kd']
    out = kp + kd*s
    return out


def Gp_Jb(s, params):
    J = params['J_motor']
    b = params['b_motor']
    out = 1.0/(J*s**2+b*s)
    return out


class DC_motor_closed_loop_torque_dist(TMMElementIHT):
    """This class represents a DC motor under closed-loop control and
    also includes a torque disturbance input.  The transfer function
    between theta and theta_d will be

     theta         G_c*G_p
    --------- = --------------
     theta_d     1 + G_c*G_p

     
    The transfer function between theta and tau, which will be a disturbance input, will be

     theta           G_p
    --------- = --------------
      tau         1 + G_c*G_p

    """
    def __init__(self, params, Gc_func=Gc_PD, Gp_func=Gp_Jb, \
                 symlabel='dccl', symname='Udccl',**kwargs):
        """Create a DC motor instance.  params is a
        dictionary with keys required by Gc_func and Gp_func."""
        TMMElementIHT.__init__(self,74,params,symlabel=symlabel,symname=symname,maxsize=4,**kwargs)
        self.Gc_func = Gc_func
        self.Gp_func = Gp_func
        

    def GetMat(self,s,sym=False):
        """Return the element transfer matrix for the
        TorsionalSpringDamper element.  If sym=True, 's' must be a
        symbolic string and a matrix of strings will be returned.
        Otherwise, 's' is a numeric value (probably complex) and the
        matrix returned will be complex."""
        N = self.maxsize
        if sym:
            myparams = self.symparams
        else:
            myparams = self.params
            
        Gc = self.Gc_func(s, myparams)
        Gp = self.Gp_func(s, myparams)

        G_tau = Gp/(1+Gc*Gp)
        #G_theta_d = Gc*Gp/(1+Gc*Gp)
        
        if sym:
            maxlen = len(G_tau) + 100
            matout = eye(N,dtype = 'f')
            matout = matout.astype('S%d'%maxlen)
        else:
            matout = eye(N,dtype='D')
        matout[1,2] = G_tau
        #matout[1,4] = G_theta_d
        return matout


    def GetAugMat(self, s, sym=False):
        """Return the augmented element transfer matrix for the
        AngularVelocitySource element, which includes the velocity
        source portion of 1/s in the augmentend column for theta.  If
        sym=True, 's' must be a symbolic string and a matrix of
        strings will be returned.  Otherwise, 's' is a numeric value
        (probably complex) and the matrix returned will be complex."""
        N = self.maxsize
        if sym:
            myparams=self.symparams
            matout=eye(N+1,dtype='f')
            matout=matout.astype('S30')
        else:
            matout=eye(N+1,dtype='D')
            myparams=self.params
        myrow = 1# hard coding for now#(self.params['axis']-1)*4+1#axis should be 1, 2, or 3
        fourbyfour = self.GetMat(s, sym=sym)
        matout[0:4,0:4] = fourbyfour

        Gc = self.Gc_func(s, myparams)
        Gp = self.Gp_func(s, myparams)
        G_theta_d = Gc*Gp/(1+Gc*Gp)

        matout[myrow,N] = G_theta_d
        return matout


    def GetMaximaString(self,name=None,label=None,N=None,aug=0):#omit
        defs=[]
        params=[]
        if 1:
            raise NotImplementedError
        #the stuff below here is copied from another element type
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


class DC_motor_closed_loop_hold_at_zero(DC_motor_closed_loop_torque_dist):
    """This is a DC motor as well, but it is trying to hold the joint
    at 0, i.e. theta_d=0."""
    def GetAugMat(self,s,sym=False):
        outmat = DC_motor_closed_loop_torque_dist.GetAugMat(self, s, sym=sym)
        outmat[1,4] = 0
        return outmat
    
