from __future__ import division
from scipy import cos, cosh, sin, sinh, array, r_, c_, exp, pi, dot, real, imag, zeros, eye, shape
#import MLab
import scipy
from scipy.linalg import det
import copy, re
import pdb

from TMM.TMMElement import TMMElementLHT, HT4

from rwkmisc import symstr, SymstrMattoMaxima
import rwkmisc

class RigidMass(TMMElementLHT):
    """This class models a rigid mass element in a TMM model.  It can be 1, 2, or 3 dimensional."""
    def __init__(self,params,**kwargs):
        """Create a rigid mass element.  params should be a dictionary
         with keys 'L', 'm', 'r', and 'I'.  'L' is the length of the
         element.  'm' is the mass of the element.  'r' is the
         distance to the center of gravity of the link.  'I' is the
         second moment of inertia and it should be a vector if this is
         a 2 or 3D link.

        kwargs should contain {'usez':True} if a 4x4 rigid mass
        element should model rotation about the z-axis."""
        if kwargs.has_key('usez'):
            self.usez=kwargs['usez']
            kwargs.pop('usez')
        else:
            self.usez=False
        TMMElementLHT.__init__(self,'rigid',params,**kwargs)

    def GetMat(self, s, sym=False):
        """Return the element transfer matrix for the RigidMass
        element.  If sym=True, 's' must be a symbolic string and a
        matrix of strings will be returned.  Otherwise, 's' is a
        numeric value (probably complex) and the matrix returned will
        be complex."""
        if sym:
            myparams=self.symparams
        else:
            myparams=self.params
        if self.maxsize==4 and self.usez:
            rigidmat1=rigidmatz(s,myparams)
        else:
            rigidmat1=rigidmaty(s,myparams)
        if self.maxsize==4:
            return rigidmat1
        elif self.maxsize>4:
            rigidmat2=rigidmatz(s,myparams)
            zmat=scipy.zeros(scipy.shape(rigidmat2))
        if self.maxsize==8:
            bigmat1=c_[rigidmat1,zmat]
            bigmat2=c_[zmat,rigidmat2]
            temp=r_[bigmat1,bigmat2]
            return Transform8by8(temp)
        elif self.maxsize==12:
            rigidmat0=rigidmatx(s,myparams)
            row1=c_[rigidmat0,zmat,zmat]
            t1=c_[rigidmat1,zmat]
            t2=c_[zmat,rigidmat2]
            temp=r_[t1,t2]
            temp=Transform8by8(temp)
            part2=c_[scipy.zeros((8,4)),temp]
            return r_[row1, part2]


    def GetMaximaString(self,name=None,label=None,N=None,aug=0):#omit
        if label is None:
            if self.symlabel:
                label=self.symlabel
            else:
                label='rz'
                self.symlabel=label
        if name is None:
            if self.symname:
                name=self.symname
            else:
                name='Urigidz'
        maximalines=[]
        if N is None:
            N=self.maxsize
        Nstr=str(N)
        rzlines,defs,params=SymRigidMatz(name,label)
        return rzlines, name, defs, params, label

    def GetMaximaLatexString(self,name=None,label=None,wrap=0,N=None,intro=None,aug=0):#omit
        if intro is None:
            intro='The transfer matrix for a rigid mass element is given by'
        return TMMElementLHT.GetMaximaLatexString(self,name=name,label=label,wrap=wrap,N=N,intro=intro,aug=aug)

#    def GetHT(self):
#        return HT4(x=self.params['L'])

class RigidMass_vertg(RigidMass):#omit me
    def GetMaximaString(self,name=None,label=None,N=None,aug=0):
        if label is None:
            if self.symlabel:
                label=self.symlabel
            else:
                label='rz'
                self.symlabel=label
        if name is None:
            if self.symname:
                name=self.symname
            else:
                name='Urigidz'
        maximalines=[]
        if N is None:
            N=self.maxsize
        Nstr=str(N)
        rzlines,defs,params=SymRigidMatz_vertg(name,label)
        defs.append('g=9.81')
        return rzlines, name, defs, params, label

def rigidmatx(s,params):
    """Return the 4x4 transfer matrix for motion of a rigid mass about
    its x-axis."""
    m=params['m']
    Ix=params['I'][0]
    #this matrix assumes +w_x is the state [w_x; theta_x; M_x; V_x] and is derived in 
    #Rmatrix_3D.m in E:\GT\Research\SAmii\modeling\transfer_matrix\threeD_TMM_derivations
#    [      1,      0,      0,      0]
#    [      0,      1,      0,      0]
#    [      0, s^2*Ix,      1,      0]
#    [  m*s^2,      0,      0,      1]
    matout=array([[1.,0,0,0],[0,1.,0,0],[0,s**2*Ix,1.,0],[m*s**2,0,0,1.]])
    return matout

def SymRigidMatz(name='Urigidz',label='rz'):#omit
    linesout=[]
    rstr='r'+label
    Lstr='L'+label
    mstr='m'+label
    Izstr='I'+label
    matstr=name+':matrix([1,L,0,0],[0,1,0,0],[-m*s^2*(L-r),s^2*Iz-m*s^2*r*(L-r),1,-L],[m*s^2,m*s^2*r,0,1])'
    params=[]
    replist=[('r',rstr),('L',Lstr),('m',mstr),('Iz',Izstr)]
    for fi, ri in replist:
        params.append(ri)
        matstr=re.sub('\\b'+fi+'\\b',ri,matstr)
    linesout.append(matstr)
    return linesout, [], params

def SymRigidMatz_vertg(name='Urigidz',label='rz'):#omit
    linesout=[]
    rstr='r'+label
    Lstr='L'+label
    mstr='m'+label
    Izstr='I'+label
    matstr=name+':matrix([1,L,0,0],[0,1,0,0],[-m*s^2*(L-r),s^2*Iz-m*s^2*r*(L-r)+m*g*r,1,-L],[m*s^2,m*s^2*r,0,1])'
    params=[]
    replist=[('r',rstr),('L',Lstr),('m',mstr),('Iz',Izstr)]
    for fi, ri in replist:
        params.append(ri)
        matstr=re.sub('\\b'+fi+'\\b',ri,matstr)
    linesout.append(matstr)
    return linesout, [], params

def rigidmatz_vertg(s,params,g=9.81):#omit
    """Return the 4x4 transfer matrix for motion of a rigid mass about
    its z-axis, including the effect of gravity."""
    L=params['L']
    r=params['r']
    m=params['m']
    It=params['I']
    if shape(It):
        Iz=It[2]
    else:
        Iz=It
    #this matrix assumes +w_y is the state [w_y; theta_z; M_z; V_y] and is derived in 
    #Rmatrix_3D.m in E:\GT\Research\SAmii\modeling\transfer_matrix\threeD_TMM_derivations
#    [                    1,                    L,                    0,                    0]
#    [                    0,                    1,                    0,                    0]
#    [         -m*s^2*(L-r), s^2*Iz-m*s^2*r*(L-r)-m*g*r,                    1,                   -L]
#    [                m*s^2,              m*s^2*r,                    0,                    1]
    matout=array([[1.,L,0,0],\
                  [0,1.,0,0],\
                  [-m*s**2*(L-r),s**2*Iz-m*s**2*r*(L-r)+m*g*r,1.,-L],\
                  [m*s**2,m*s**2*r,0,1.]])
    return matout

def rigidmatz(s,params):
    """Return the 4x4 transfer matrix for motion of a rigid mass about
    its z-axis."""
    L=params['L']
    r=params['r']
    m=params['m']
    It=params['I']
    if shape(It):
        Iz=It[2]
    else:
        Iz=It
    #this matrix assumes +w_y is the state [w_y; theta_z; M_z; V_y] and is derived in 
    #Rmatrix_3D.m in E:\GT\Research\SAmii\modeling\transfer_matrix\threeD_TMM_derivations
#    [                    1,                    L,                    0,                    0]
#    [                    0,                    1,                    0,                    0]
#    [         -m*s^2*(L-r), s^2*Iz-m*s^2*r*(L-r),                    1,                   -L]
#    [                m*s^2,              m*s^2*r,                    0,                    1]
    matout=array([[1.,L,0,0],\
                  [0,1.,0,0],\
                  [-m*s**2*(L-r),s**2*Iz-m*s**2*r*(L-r),1.,-L],\
                  [m*s**2,m*s**2*r,0,1.]])
    return matout

def rigidmaty(s,params):
    """Return the 4x4 transfer matrix for motion of a rigid mass about
    its y-axis."""
    L=params['L']
    r=params['r']
    m=params['m']
    Imat=params['I']
    if shape(Imat):
        if max(shape(Imat))>=1:
            Iy=Imat[1]
        else:
            Iy=Imat[0]
    else:
        Iy=Imat
    #this matrix assumes +w_z is the state [w_z; theta_y; M_y; V_z] and is derived in 
    #Rmatrix_3D.m in E:\GT\Research\SAmii\modeling\transfer_matrix\threeD_TMM_derivations
#    [                    1,                   -L,                    0,                    0]
#    [                    0,                    1,                    0,                    0]
#    [          m*s^2*(L-r), s^2*Iz-m*s^2*r*(L-r),                    1,                    L]
#    [                m*s^2,             -m*s^2*r,                    0,                    1]
    matout=array([[1.,-L,0,0],[0,1.,0,0],[m*s**2*(L-r),s**2*Iy-m*s**2*r*(L-r),1.,L],[m*s**2,-m*s**2*r,0,1.]])
    return matout


class RigidMass2x2(TMMElementLHT):
    """The class models a 2 by 2 mass element where the states are
    displacement and force.  Used to model carts in simple
    mass-spring-damper systems."""
    def __init__(self,params,symlabel='M',**kwargs):
        """Create a 2 by 2 mass element.  params is a dictionary with
        only one key: params={'m':##.##} (the mass of the cart)."""
        TMMElementLHT.__init__(self,'rigid',params,maxsize=2,symlabel=symlabel,**kwargs)

    def GetMat(self,s,sym=False):
        """Return the 2 by 2 transfer matrix for a cart mass element."""
        N=self.maxsize
        if sym:
            myparams=self.symparams
        else:
            myparams=self.params
        m=myparams['m']
        maxlen=20
        if sym:#type(s)==type(rwkmisc.symstr('s')):
            matout=eye(N,dtype='f')
            matout=matout.astype('S%d'%maxlen)
        else:
            matout=eye(N,dtype='D')
        matout[1,0]=m*s**2
        return matout
