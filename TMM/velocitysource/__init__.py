from __future__ import division
from scipy import cos, cosh, sin, sinh, array, r_, c_, exp, pi, dot, real, imag, zeros, eye, shape, atleast_1d
#import MLab
import scipy
from scipy.linalg import det
import copy
#import pdb
from  IPython.Debugger import Pdb
from TMM.TMMElement import TMMElement, HT4, TMMElementLHT
import rwkmisc
from rwkmisc import symstr, SymstrMattoMaxima

class AngularVelocitySource(TMMElementLHT):
    """This class models an angular velocity source.  This class is
    used to model rotary hydraulic actuators."""
    def __init__(self,params={},**kwargs):
        """Initialize an instance of the AngularVelocitySource class.
        params is a dictionary with keys 'K', 'tau', and 'axis'.  All
        of these keys are optional.  'K' is the gain and it defaults
        to 1.  'axis' is the axis about which the actuator rotates -
        it defaults to 1.  'tau' is the pole of the first order lag of
        the actuator (i.e. if 'tau' is given and is > 0, the trasnfer
        function of the actuator will be tau/(s*(s+tau))."""
        if not params.has_key('K'):
            params['K']=1
        if not params.has_key('axis'):
            params['axis']=1
        if shape(params['K']):
            params['K']=params['K'][0]
        if params.has_key('tau'):
            if shape(params['tau']):
                params['tau']=params['tau'][0]
            if params['tau']<=0:
                params.pop('tau')
        TMMElementLHT.__init__(self,'avs',params,**kwargs)

    def GetMat(self, s, sym=False):
        """Return the element transfer matrix for the
        AngularVelocitySource element.  If sym=True, 's' must be a
        symbolic string and a matrix of strings will be returned.
        Otherwise, 's' is a numeric value (probably complex) and the
        matrix returned will be complex.  Note that the
        AngularVelocitySource element will return an identity matrix
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
        AngularVelocitySource element, which includes the velocity
        source portion of 1/s in the augmentend column for theta.  If
        sym=True, 's' must be a symbolic string and a matrix of
        strings will be returned.  Otherwise, 's' is a numeric value
        (probably complex) and the matrix returned will be complex."""
        N=self.maxsize
        if sym:
            myparams=self.symparams
            matout=eye(N+1,dtype='f')
            matout=matout.astype('S30')
        else:
            matout=eye(N+1,dtype='D')
            myparams=self.params
        myrow=(self.params['axis']-1)*4+1#axis should be 1, 2, or 3
        K=myparams['K']
        if myparams.has_key('tau'):
            tau=myparams['tau']
            matout[myrow,N]=K*tau/(s*(s+tau))
        else:
            matout[myrow,N]=K/s
        return matout

#    def GetHT(self):
#        return HT4()

    def GetMaximaString(self,name=None,label=None,N=None,aug=0):#omit
        defs=[]
        params=[]
        if label is None:
            if self.symlabel:
                label=self.symlabel
            else:
                label='avs'
                self.symlabel=label
        if name is None:
            if self.symname:
                name=self.symname
            else:
                name='Uavs'
        maximalines=[]
        if N is None:
            N=self.maxsize
        Nstr=str(N)
        Np1str=str(N+1)
        maximalines.append(name+':ident('+Nstr+')')
        return maximalines, name,defs,params,label

    def GetAugMaximaString(self,name=None,label=None,N=None):#omit
        defs=[]
        params=[]
        if label is None:
            if self.symlabel:
                label=self.symlabel
            else:
                label='avs'
                self.symlabel=label
        if name is None:
            if self.symname:
                name=self.symname
            else:
                name='Uavs'
        maximalines=[]
        if N is None:
            N=self.maxsize
        Nstr=str(N)
        Np1str=str(N+1)
        maximalines.append(name+':ident('+Np1str+')')
        x=self.params['axis']-1
        r=4*x+2
        c=N+1
        gainstr='K'+label
        taustr='tau'+label
        if self.params.has_key('tau'):
            Gpstr=gainstr+'*'+taustr+'/(s*(s+'+taustr+'))'
        else:
            Gpstr=gainstr+'/s'
        maximalines.append(name+'['+str(r)+','+str(c)+']:'+Gpstr)
        params.append(gainstr)
        if self.params.has_key('tau'):
            params.append(taustr)
        for line in maximalines:
            if line[-1]!=';':
                line+=';'
        return maximalines, name, defs, params, label

    def GetMaximaLatexString(self,name=None,label=None,wrap=0,N=None,intro=None,aug=0):#omit
        if intro is None:
            intro='The transfer matrix for an angular velocity source element is given by'
        return TMMElementLHT.GetMaximaLatexString(self,name=name,label=label,wrap=wrap,N=N,intro=intro,aug=aug)

class AVSwThetaFB(TMMElementLHT):
    """This class models the closed-loop response of an angular
    velocity source with compliance and relative theta feedback."""
    def __init__(self,params={}, **kwargs):
        """Initialize the AVSwThetaFB instance.  params should have
        the following keys:

        'ks'   - spring constant (required)
        'c'    - damper coeffiecent (defaults to 0)
        'Ka'   - actuator gain (defaults to 1)
        'Gc'   - proportional feedback gain (defaults to 1)
                 (Gc can be a transfer function modelled as a
                 ratio of polynomials when passed to FORTRAN.)
        'axis' - axis about which the actuator rotates (defaults to 1)
        'tau'  - first order pole of actuator."""
        if not params.has_key('Ka'):
            params['Ka']=1
        if not params.has_key('axis'):
            params['axis']=1
        if not params.has_key('c'):
            params['c']=0
        if not params.has_key('Gc'):
            params['Gc']=1
        if shape(params['Ka']):
            params['Ka']=params['Ka'][0]
        if params.has_key('tau'):
            if shape(params['tau']):
                params['tau']=params['tau'][0]
        TMMElementLHT.__init__(self,'avsthfb',params,**kwargs)

    def GetMat(self,s):
        """Return the element transfer matrix for the AVSwThetaFB
        element.  If sym=True, 's' must be a symbolic string and a
        matrix of strings will be returned.  Otherwise, 's' is a
        numeric value (probably complex) and the matrix returned will
        be complex.  Note that the AngularVelocitySource element will
        return an identity matrix for its non-augmented transfer
        matrix."""
        N=self.maxsize
        tempmap=self.GetAugMat(s)
        matout=tempmap[0:N,0:N]
        return matout

    def GetAugMat(self, s, sym=False):
        """Return the augmented element transfer matrix for the
        AVSwThetaFB element.  If sym=True, 's' must be a symbolic
        string and a matrix of strings will be returned.  Otherwise,
        's' is a numeric value (probably complex) and the matrix
        returned will be complex."""
        N=self.maxsize
        if N%2:
            N=N-1
#        matout=eye(N+1,'d')+0.0j
        if sym:
            myparams=self.symparams
            matout=eye(N+1,dtype='f')
            matout=matout.astype('S30')
        else:
            matout=eye(N+1,dtype='D')
            myparams=self.params
#        matout=eye(N+1,dtype='D')
        myrow=(self.params['axis']-1)*4+1#axis should be 1, 2, or 3
        mycol=myrow+1
        Gc=myparams['Gc']
        gain=myparams['Ka']
        c=myparams['c']
        k=myparams['ks']
        if myparams.has_key('tau'):
            tau=myparams['tau']
            Gp=gain*tau/(s*(s+tau))
        else:
            Gp=gain/s
        actpart=Gc*Gp/(Gc*Gp+1.0)
        flexpart=1.0/((Gc*Gp+1)*(c*s+k))
        matout[myrow,N]=actpart
        matout[myrow,mycol]=flexpart
        return matout

    def GetMaximaString(self,name=None,label=None,N=None,aug=0):#omit
        defs=[]
        params=[]
        if label is None:
            if self.symlabel:
                label=self.symlabel
            else:
                label='avsfb'
                self.symlabel=label
        if name is None:
            if self.symname:
                name=self.symname
            else:
                name='Uavsfb'
        maximalines=[]
        if N is None:
            N=self.maxsize
        Nstr=str(N)
        Np1str=str(N+1)
        maximalines.append(name+':ident('+Nstr+')')
        useflex=0
        k=0.0
        c=0.0
        gcstr='Gc'+label
        gpstr='Gp'+label
        kastr='Ka'+label
        taustr='tau'+label
        if self.params.has_key('ks'):
            useflex=1
            kstr='ks'+label
        if self.params.has_key('c'):
            if self.params['c']!=0.0:
                useflex=1
                cstr='c'+label
        if useflex:
            x=self.params['axis']-1
            r=4*x+2
            cf=r+1
#            actpart=Gc*Gp/(Gc*Gp+1.0)
            maximalines.append('flexpart:1/(('+gcstr+'*'+gpstr+'+1)*('+cstr+'*s+'+kstr+'))')
            maximalines.append(name+'['+str(r)+','+str(cf)+']:flexpart')
            params.append(kstr)
            params.append(cstr)
        Gpdef, Kparam=self.SymGp(kastr,gpstr,taustr)
        defs.append(Gpdef)
        params.append('Ka'+label)
        params.append('Gc'+label)
        params.append(taustr)
        return maximalines, name,defs,params,label

    def GetComplexList(self,label=None): 
        """This function returns a list of all variables associated
        with this element that would need to be declared complex in a
        FORTRAN funciton that includes the element in a Bode model.
        This funcion is used in symbolic work."""
        strlist=['Gc','Gp']
        if label is None:
            if self.symlabel:
                label=self.symlabel
            else:
                label='bz'
                self.symlabel=label
        listout=[item+label for item in strlist]
        return listout

    def GetAugMaximaString(self,name=None,label=None,N=None,subGp=False):#omit
        print('subGp='+str(subGp))
        if label is None:
            if self.symlabel:
                label=self.symlabel
            else:
                label='avsfb'
                self.symlabel=label
        maxlines, myname,defs,params,dummylabel=TMMElementLHT.GetAugMaximaString(self,name,label,N)
        if N is None:
            N=self.maxsize
        Np1str=str(N+1)
        x=self.params['axis']-1
        r=4*x+2
        gcstr='Gc'+label
        gpstr='Gp'+label
        kastr='Ka'+label
        taustr='tau'+label
        maxlines.append('actpart:'+gcstr+'*'+gpstr+'/('+gcstr+'*'+gpstr+'+1.0)')
        maxlines.append(myname+'['+str(r)+','+Np1str+']:actpart')
        if subGp:
            Gpstr, gainlabel=self.SymGp(kastr,gpstr,taustr)
            maxlines.append(myname+':at('+myname+','+Gpstr+')')
#            pdb.set_trace()
            ind=-1
            for x,ent in enumerate(defs):
                if ent.find(gpstr+'=')==0:
                    ind=x
                    break
            defs.pop(x)
        return maxlines, myname,defs,params, label

    def SymGp(self,gainlabel='Kaavs',Gplabel='Gp',taulabel='tau'):#omit
#        linesout.append('\\begin{maxima}')
#        linesout.append("\t\\parseopts{lhs='G_p'}")
        if self.params.has_key('tau'):
            eqstr=Gplabel+'='+gainlabel+'*'+taulabel+'/(s*(s+'+taulabel+'))'
        else:
            eqstr=Gplabel+'='+gainlabel+'/s'
#        linesout.append(eqstr)
#        linesout.append('\\end{maxima}')
#        return linesout
        return eqstr, gainlabel

#    def GetHT(self):
#        return HT4()

class AVSwArbFB(AVSwThetaFB):#omit
    def __init__(self,params={},maxsize=12,label='',):
        if not params.has_key('K'):
            params['K']=1
        if not params.has_key('axis'):
            params['axis']=1
        if not params.has_key('c'):
            params['c']=0
        if not params.has_key('Gc'):
            params['Gc']=1
        if not params.has_key('Gx'):
            params['Gx']=0
        if not params.has_key('Gm'):
            params['Gm']=0
        if not params.has_key('Gv'):
            params['Gv']=0
        if shape(params['K']):
            params['K']=params['K'][0]
        if params.has_key('tau'):
            if shape(params['tau']):
                params['tau']=params['tau'][0]
        TMMElementLHT.__init__(self,'avsthfb',params,maxsize=maxsize,label=label)

    def GetAugMat(self,s):
        N=self.maxsize
        if N%2:
            N=N-1
#        matout=eye(N+1,'d')+0.0j
        matout=eye(N+1,dtype='D')
        myrow=(self.params['axis']-1)*4+1#axis should be 1, 2, or 3
        xcol=myrow-1
        Mcol=myrow+1
        Vcol=myrow+2
        Gc=self.params['Gc']
        Gm=self.params['Gm']
        Gx=self.params['Gx']
        Gv=self.params['Gv']
        gain=self.params['K']
        c=self.params['c']
        k=self.params['ks']
        if self.params.has_key('tau'):
            tau=self.params['tau']
            Gp=gain*tau/(s*(s+tau))
        else:
            Gp=gain/s
        dp=Gp/(1+Gc*Gp)
        actpart=Gc*Gp/(Gc*Gp+1.0)
        flexpart=1.0/((Gc*Gp+1)*(c*s+k))
        matout[myrow,N]=actpart
        matout[myrow,Mcol]=flexpart+Gm*dp
        matout[myrow,xcol]=Gx*dp
        matout[myrow,Vcol]=Gv*dp
        return matout

    def GetMaximaString(self,name=None,label=None,N=None,aug=0):
        defs=[]
        params=[]
        if label is None:
            if self.symlabel:
                label=self.symlabel
            else:
                label='avsfb'
        if name is None:
            if self.symname:
                name=self.symname
            else:
                name='Uavsfb'
        maximalines=[]
        if N is None:
            N=self.maxsize
        Nstr=str(N)
        Np1str=str(N+1)
        maximalines.append(name+':ident('+Nstr+')')
        useflex=0
        k=0.0
        c=0.0
        fblist=['Gx','Gm','Gv']
        op='Gp/(1+Gc*Gp)'
        fbparts=[item+'*'+op for item in fblist]
        if self.params.has_key('ks'):
            useflex=1
            kstr='ks'+label
        if self.params.has_key('c'):
            if self.params['c']!=0.0:
                useflex=1
                cstr='c'+label
        if useflex:
            x=self.params['axis']-1
            r=4*x+2
            cf=r+1
#            actpart=Gc*Gp/(Gc*Gp+1.0)
            maximalines.append('flexpart:1/((Gc*Gp+1)*('+cstr+'*s+'+kstr+'))')
            maximalines.append(name+'['+str(r)+','+str(cf)+']:flexpart+'+fbparts[1])
            params.append(kstr)
            params.append(cstr)
        else:
            maximalines.append(name+'['+str(r)+','+str(cf)+']:'+fbparts[1])
        maximalines.append(name+'['+str(r)+','+str(r-1)+']:'+fbparts[0])
        maximalines.append(name+'['+str(r)+','+str(cf+1)+']:'+fbparts[2])
        Gpdef, Kparam=self.SymGp()
        defs.append(Gpdef)
        params.append(Kparam)
        params.append('Gc')
        params.extend(fblist)
        return maximalines, name,defs,params

    def GetAugMaximaString(self,name=None,label=None,N=None,subGp=False):
        print('subGp='+str(subGp))
        maxlines, myname,defs,params,dummylabel=TMMElementLHT.GetAugMaximaString(self,name,label,N)
        if N is None:
            N=self.maxsize
        Np1str=str(N+1)
        x=self.params['axis']-1
        r=4*x+2
        gcstr='Gc'+label
        gpstr='Gp'+label
        kastr='Ka'+label
        taustr='tau'+label
        maxlines.append('actpart:'+gcstr+'*'+gpstr+'/('+gcstr+'*'+gpstr+'+1.0)')
#        maxlines.append('actpart:Gc*Gp/(Gc*Gp+1.0)')
        maxlines.append(myname+'['+str(r)+','+Np1str+']:actpart')
        if subGp:
            Gpstr, gainlabel=self.SymGp(kastr,gpstr,taustr)
#            Gpstr, gainlabel=self.SymGp()
            maxlines.append(myname+':at('+myname+','+Gpstr+')')
#            pdb.set_trace()
            ind=-1
            for x,ent in enumerate(defs):
                if ent.find('Gp=')==0:
                    ind=x
                    break
            defs.pop(x)
        return maxlines, myname,defs,params

