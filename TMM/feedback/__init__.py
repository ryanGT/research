from __future__ import division
from scipy import cos, cosh, sin, sinh, array, r_, c_, exp, pi, dot, real, imag, zeros, eye, shape, atleast_1d
#import MLab
import scipy
from scipy.linalg import det
import copy
import pdb

from TMM.TMMElement import TMMElement, HT4, TMMElementLHT

from IPython.core.debugger import Pdb

import rwkmisc
reload(rwkmisc)
#import ModelSpec

class SAMIIAccelFB(TMMElementLHT):
    """The class represents the accelerometer feedback element for
    SAMII.  It is a very specific class that will not be much use in
    other models."""
    def __init__(self,link0,joint1,link1,avs,Ga=1,axis=1,maxsize=4,label='',symlabel='accelfb',symname='Uaccelfb',compensators='Ga',**extraargs):
        """Initializes the accelerometer feedback element.  link0,
        joint1, link1, and avs are all TMMElement's that must be
        passed in.  Their parameters are used in creating the transfer
        matrix for acceleration feedback."""
        l0params=link0.params
#        a = -(L0*m1*r1*s^2-c1*s-k1)/(c1*s+k1);
#        b = -((L0*m1*r1^2-L0*L1*m1*r1+Iz1*L0)*s^2+(c1*L1+c1*L0)*s+k1*L1+k1*L0)/(c1*s+k1);
#        c = L0/(c1*s+k1);
#        d = L0*L1/(c1*s+k1);
        L0=l0params['L']
        l1params=link1.params
        L1=l1params['L']
        m1=l1params['m']
        r1=l1params['r']
        Iz1=l1params['I']
        j1params=joint1.params
        c1=j1params['c']
        k1=j1params['k']
        gain=avs.params['Ka']
        Gc=avs.params['Gc']
        params={'L0':L0,'L1':L1,'m1':m1,'r1':r1,'Iz1':Iz1,'c1':c1,'k1':k1,'axis':axis,'gain':gain,'Gc':Gc,'Ga':Ga}
        if avs.params.has_key('tau'):
            params['tau']=avs.params['tau']
        for curkey,curvalue in params.iteritems():
            if curkey!='Ga':
                if shape(curvalue):
                    curvalue=copy.deepcopy(curvalue)[0]
                    params[curkey]=curvalue
        TMMElementLHT.__init__(self,'accelfb',params,maxsize=maxsize,label=label,symlabel=symlabel,symname=symname,compensators=compensators,**extraargs)
        #set up symbolic params
        sL0='L'+link0.symlabel
        sL1='L'+link1.symlabel
        sm1='m'+link1.symlabel
        sr1='r'+link1.symlabel
        sIz1='I'+link1.symlabel
        sc1='c'+joint1.symlabel
        sk1='k'+joint1.symlabel
        stau='tau'+avs.symlabel
        sgain='Ka'+avs.symlabel
        sGc='Gc'+avs.symlabel
        sGa='Ga'
        values=[sL0,sL1,sm1,sr1,sIz1,sc1,sk1,stau,sgain,sGc,sGa]
        keys=['L0','L1','m1','r1','Iz1','c1','k1','tau','gain','Gc','Ga']
        symparams={}
        for key, value in zip(keys,values):
            symparams[key]=rwkmisc.symstr(value)
#            params[key]=ModelSpec.Par(value)
        symparams['axis']=axis
        self.symparams=symparams

    def GetSymAugMat(self,s):
        """Return the augmented element transfer matrix for the
        SAMIIAccelFB element.  If sym=True, 's' must be a
        symbolic string and a matrix of strings will be returned.
        Otherwise, 's' is a numeric value (probably complex) and the
        matrix returned will be complex."""
        s=rwkmisc.symstr(s)
        return self.GetAugMat(s,self.symparams)

    def GetComplexList(self):
        """This function returns a list of all variables associated
        with this element that would need to be declared complex in a
        FORTRAN funciton that includes the element in a Bode model.
        This funcion is used in symbolic work."""
        return ['Ga']

    def GetAugMaximaString(self,**kwargs):#omit
        symstrmat=self.GetSymAugMat('s')
        defs=[]
        params=[]
        maximalines=[]
        matstr=''
        first=1
        for currow in symstrmat.tolist():
            tempstr=str(currow)
            rowstr=tempstr.replace("'","")
            if first:
                first=0
            else:
                matstr+=','
            matstr+=rowstr

#        maximalines.append('The transfer matrix for the accelerometer feedback element is given by')
        maximalines.append(self.symname+':matrix('+matstr+')')
        params.append('Ga')
        return maximalines, self.symname, defs, params, self.symlabel

    def GetMaximaString(self,**kwargs):#omit
        symstrmat=self.GetSymAugMat('s')
        symstrmat=symstrmat[0:-1,0:-1]#just chop off last row and column
        defs=[]
        params=[]
        maximalines=[]
        matstr=''
        first=1
        for currow in symstrmat.tolist():
            tempstr=str(currow)
            rowstr=tempstr.replace("'","")
            if first:
                first=0
            else:
                matstr+=','
            matstr+=rowstr
#        maximalines.append('The transfer matrix for the accelerometer feedback element is given by')
        maximalines.append(self.symname+':matrix('+matstr+')')
        params.append('Ga')
        return maximalines, self.symname, defs, params, self.symlabel


    def GetMat(self,s):#omit
        raise NotImplementedError, 'AngularVelocitySource must use augmented transfer matrix'

    def GetAugMat(self,s,myparams=None):
        """Return the augmented element transfer matrix for the
        SAMIIAccelFB element.  If sym=True, 's' must be a symbolic
        string and a matrix of strings will be returned.  Otherwise,
        's' is a numeric value (probably complex) and the matrix
        returned will be complex."""
        if myparams is None:
            myparams=self.params
        L0=myparams['L0']
        L1=myparams['L1']
        m1=myparams['m1']
        r1=myparams['r1']
        Iz1=myparams['Iz1']
        c1=myparams['c1']
        k1=myparams['k1']
#        Pdb().set_trace()
#        anum=L0*m1*r1*s**2-c1*s-k1
#        aden=c1*s+k1
#        a=anum/aden*(-1)
        a=(-1)*(L0*m1*r1*s**2-c1*s-k1)/(c1*s+k1)
#        test=(L0*m1*r1*s**2-c1*s-k1)/(c1*s+k1)*(-1)
        b=(-1)*((L0*m1*r1**2-L0*L1*m1*r1+Iz1*L0)*s**2+(c1*L1+c1*L0)*s+k1*L1+k1*L0)/(c1*s+k1)
        c=L0/(c1*s+k1)
        d=L0*L1/(c1*s+k1)
        Gc=myparams['Gc']
        N=self.maxsize
        if N%2:
            N-=1
        myrow=(myparams['axis']-1)*4+1#axis should be 1, 2, or 3 ->myrow = 1, 5, or 9
        startcol=myrow-1
        if myparams.has_key('tau'):
            tau=myparams['tau']
            Gp=myparams['gain']*tau/(s*(s+tau))
        else:
            Gp=myparams['gain']/s
#        Pdb().set_trace()
        Gat=myparams['Ga']
        #Ga is assumed to either be a list or tuple that has the transfer function numerator and denominator as scipy.poly1d or it is a constant
        if type(Gat)==list or type(Gat)==tuple:
            num=Gat[0]
            den=Gat[1]
            Ga=num(s)/den(s)
        else:
            Ga=Gat
        mylist=[a,b,c,d]
#        print('just before mainpart, Ga='+str(Ga))
        mainpart=Ga*Gc*Gp*s**2/(Gc*Gp+1.0)
#        Pdb().set_trace()
        if isinstance(s,rwkmisc.symstr):
            print('Ga='+str(Ga))
            templist=[item*mainpart for item in mylist]
            lengths=[len(item) for item in templist]
            maxlen=max(lengths)
            maxlen+=100#give ourselves a little extra string room
            matout=eye(N+1,dtype='f')#I don't think any of the symbolic variables need to be complex and this will save parsing problems with Maxima and Nj-->N*%I
            matout=matout.astype('S%d'%maxlen)
        else:
    #        matout=eye(N+1,'d')+0.0j
            matout=eye(N+1,dtype='D')
        for x,co in enumerate(mylist):
            if x==1:
                matout[myrow,startcol+x]=co*mainpart+1.0
            else:
                matout[myrow,startcol+x]=co*mainpart
#        MATRIX([1,0,0,0,0],[a*Ga*Gc*Gp*s^2/(Gc*Gp+1),b*Ga*Gc*Gp*s^2/(Gc*Gp+1)+1,c*Ga*Gc*Gp*s^2/(Gc*Gp+1),d*Ga*Gc*Gp*s^2/(Gc*Gp+1),0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]);
        return matout

#    def GetHT(self):
#        return HT4()

class SAMIIAccelFBsym(SAMIIAccelFB):#omit
    def __init__(self,link0symlabel,joint1symlabel,link1symlabel,avssymlabel,Ga=1,axis=1,maxsize=4,label='',**kwargs):
        L0='L'+link0symlabel
        L1='L'+link1symlabel
        m1='m'+link1symlabel
        r1='r'+link1symlabel
        Iz1='I'+link1symlabel
        c1='c'+joint1symlabel
        k1='k'+joint1symlabel
        tau='tau'+avssymlabel
        gain='Ka'+avssymlabel
        Gc='Gc'+avssymlabel
        values=[L0,L1,m1,r1,Iz1,c1,k1,tau,gain,Gc,Ga]
        keys=['L0','L1','m1','r1','Iz1','c1','k1','tau','gain','Gc','Ga']
        params={}
        for key, value in zip(keys,values):
            params[key]=rwkmisc.symstr(value)
#            params[key]=ModelSpec.Par(value)
        params['axis']=axis
        TMMElementLHT.__init__(self,'accelfb',params,maxsize=maxsize,label=label,**kwargs)

    def GetAugMat(self,s):
        s=rwkmisc.symstr(s)
#        s=ModelSpec.Par(s)
        return SAMIIAccelFB.GetAugMat(self,s)
