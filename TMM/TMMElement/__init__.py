from __future__ import division
from scipy import cos, cosh, sin, sinh, array, r_, c_, exp, pi, real, imag, zeros, eye, alltrue, shape, atleast_1d, dot, vstack, isscalar
from scipy.linalg import inv as inverse
from scipy.io import save
#import MLab
import scipy
from scipy.linalg import det
import copy, os, sys
import pdb
from rref import rrefnull
from rwkmisc import RowSwap, ColSwap, rowwise, colwise, my_import
import rwkmisc
import rwkascii
from rwkdataproc import datastruct
import rwkbode
import time, types
from textfiles.latexlist import latexlist
from textfiles.latexlist import ReplaceList
from textfiles import rwkreadfile, textlist
import rwkfortran
reload(rwkfortran)
from rwkfortran import FortranToTextList, GetLHS
#from rwkmisc import null, norm2
from  IPython.Debugger import Pdb
from rwkmisc import symstr, SymstrMattoMaxima, rwkstr
import rwkmisc

import inspect
import re
from rwkparse import GetPythonFunctionArgs


class TMMElement:
    """This is the base class for creating TMMElement models.  It is
    not intended to be used directly.  Users must derive their own
    classes from it, and derived classes must override GetMat and
    GetHT."""
    def __init__(self,elemtype,params,maxsize=12,label='',symname='',symlabel='',symsub=False,unknownparams=None,functionparams=None,compensators=None, subparams=False):
        """Create a TMMElement instance.  Since the class is not
        intended to be used without deriving from it, this code should
        be called from each derived methods __init__ method.  It
        contains code that should be common to all __init__methods for
        all TMMElement classes.

        subparams was added on 5/17/07 to fascilitate substituting in
        symstr expressions for parameters.  It is currently only
        utilized in the bendmatz part of the beam element."""
        self.elemtype=self.parsetype(elemtype)
        self.params=params
        self.maxsize=maxsize
        self.label=label
        self.symlabel=symlabel
        if not symname and symlabel:
            symname='U'+symlabel
        self.symname=symname
        self.symsub=symsub
        self.subparams=subparams
    
        if functionparams is not None:
            if not shape(functionparams):
                functionparams=[functionparams]
        self.functionparams=functionparams
        self.compensators=compensators
        if isinstance(elemtype,str):
            self.elemstr=elemtype
        else:
            self.elemstr=''
        #setup symparams
        keylist=self.params.keys()
        symvals=[]
        for key, value in self.params.iteritems():
            if not isscalar(value):
                if len(value)>1:
                    inds=range(len(value))
                    cursym=[symstr(key+str(item)+self.symlabel) for item in inds]
                else:
                    cursym=[symstr(key+self.symlabel)]
            else:
                cursym=symstr(key+self.symlabel)
            symvals.append(cursym)
        self.symparams=dict(zip(keylist,symvals))
        #if unknownparams is not None:
        #    self.unknownparams=rwkmisc.flatten([self.symparams[key] for key in unknownparams])
        #else:
        self.unknownparams = unknownparams

    def GetMat(self, s, sym=False):
        """This method must be overridden in all derived TMMElement
        classes.  It returns the element transfer matrix given 's' as
        an input.  It should be defined in a way that allows it to
        output either a numeric or symbolic transfer matrix.  For
        symbolic output, 's' should be a symstr."""
        raise NotImplementedError

    def GetAugMat(self, s, sym=False, maxlen=500):
        """This method returns the augmented element transfer matrix.
        By default, it returns an N+1 by N+1 matrix with with
        self.GetMat substituted into the upper left N by N matrix.
        This is exactly what is needed for unactuated elements.
        Unactuated elements do not need to override this method.
        Actuated elements must override this method.  Provided that
        GetMat is written in a way that is compliant with either
        numeric or symbolic output, this method will be as well."""
        N=self.maxsize
        if sym:
            tempout=eye(N+1,dtype='f')
            tempout=tempout.astype('S%d'%maxlen)
        else:
            tempout=eye(N+1,dtype='D')
        tempout[0:N,0:N]=self.GetMat(s,sym=sym)
        return tempout

    def GetHT(self):
        """This method returns the homogeneous transformation matrix
        for the TMMElement.  This matrix is used in three dimensional
        visualization of mode shapes.  This method must be overridden
        in derived classes."""
        raise NotImplementedError
    
    def SubstituteUnknowns(self, unknowndict):
        """This method is used to substitute the results from system
        identification into the proper place in the TMM model.  The
        input is a dictionary containing unknown parameters and their
        values from system i.d.  This functino looks for the unknown
        parameters of this element in the unknown dictionary and
        substitutes the appropriate values.  If all of the
        unknownparams of the element are found and substituted, then
        self.unknownparams is set to None."""
        if not self.unknownparams:
            return False
        myucv=self.unknownparams
        if myucv=='all':
            myucv=self.params.keys()
        ukeys=[item+self.symlabel for item in myucv]
        foundall=True
        for key, param in zip(ukeys, myucv):
            if unknowndict.has_key(key):
                self.params[param]=unknowndict[key]
            else:
                foundall=False
        if foundall:
            self.unknownparams=None
            return True
        else:
            return False

    def GetSymMat(self):
        """This function is called by self.GetMaximaLines to get the
        symoblic element transfer function.  By default, it simply
        calls self.GetMat(s, sym=True) where s is a symstr.  This
        function could be overridden by a derived element if for some
        reason self.GetMat cannot be written in such a way that in can
        output either a symoblic or numeric transfer matrix.  As long
        as self.GetMat can return a symbolic transfer matrix for the
        element, this function does not need to be overridden by
        derived elements."""
        s=symstr('s')
        symmat=self.GetMat(s,sym=True)
        return SymstrMattoMaxima(symmat,self.symname)

    def GetAugSymMat(self):
        """Similar to GetSymMat, this method returns a symbolic
        transfer matrix for the element.  By default, it calls
        self.GetAugSymMat(s, sym=True) with s as a symstr."""
        s=symstr('s')
        symmat=self.GetAugMat(s,sym=True)
        return SymstrMattoMaxima(symmat,self.symname)

    def GetMaximaLines(self,aug=False):
        """This is the new (as of 04/04/06) method for symbolic
        analysis based on the GetSymMat and GetAugSymMat methods.
        This function outputs a list that can be appended to a
        latexlist as part of the input to the Python-Maxima-Latex
        symbolic engine."""
        mylist=[]
        out=mylist.append
        ws='\t'
        out('\\begin{maxima}')
        if self.symname[0]=='U':
            mylabel='\\M{U}_{'+self.symname[1:]+'}'
        else:
            mylabel='\\M{'+self.symname+'}'
        out(ws+"\\parseopts{lhs='"+mylabel+"'}")
        if aug:
            out(self.GetAugSymMat())
        else:
            out(self.GetSymMat())
        out('\\end{maxima}')
        return mylist

    def parsetype(self, elemtype):
        """This method assigns a numeric value for the element type
        and is called by the __init__ method."""
        elemout=-1
        if isinstance(elemtype,str):
            teststr=elemtype.lower()
            if teststr=='beam':
                elemout=1
            elif teststr.find('rigid')==0:
                elemout=2
            elif teststr.find('rot')==0:
                elemout=3
            elif teststr.find('spring')==0:
                elemout=4
            elif teststr.find('avs')==0:
                elemout=5
            elif teststr.find('forcing')==0:
                elemout=6

        else:
            elemout=elemtype
        return elemout

    def GetComplexList(self,label=None):
        """This function returns a list of all variables associated
        with an element that would need to be declared complex in a
        FORTRAN funciton that includes the element in a Bode model.
        This funcion is used in symbolic work.
        
        This is the default function defined in the top level
        TMMElement class.  It assumes that most elements have only
        real parameters and returns an empty list."""
        return []
    
    def GetFortranDefs(self):
        """This method is used to specify definitions that need to go
        at the beginning of automatically generated FORTRAN functions,
        so that FORTRAN analysis of the element will go smoothly."""
        return []

########
#    These should be illiminated, but that would mean re-writing all existing elements to use symstr based stuff
########
    def GetMaximaString(self,name=None,label=None,N=None):#omit 
        raise NotImplementedError

    def GetAugMaximaString(self,name=None,label=None,N=None):#omit
        maxlines,myname,defs,params,label=self.GetMaximaString(name,label,N)
        if N is None:
            N=self.maxsize
        Nstr=str(N)
        Np1str=str(N+1)
        maxlines.append('newcol:zeromatrix('+Nstr+',1)')
        maxlines.append('newrow:zeromatrix(1,'+Np1str+')')
        maxlines.append('newrow[1,'+Np1str+']:1')
        maxlines.append(myname+':addcol('+myname+',newcol)')
        maxlines.append(myname+':addrow('+myname+',newrow)')
        return maxlines, myname,defs,params,label

    def GetSubLines(self):#omit
        """This function doesn't actually work yet.  Don't use it."""
        listout=[]
        out=listout.append
        out('The following substitutions will be made:')
        substr='['
        first=1
        for cp,value in self.params.iteritems():
            if first:
                first=0
            else:
                substr+=','
            cursub=cp+label+'='+str(value)
            sublines.append(cursub)
            substr+=cursub
            maximalines.append('\\begin{equation}')
            maximalines.append('\t'+cursub)
            maximalines.append('\\end{equation}')
        paramsubs=substr+']'
        maximalines.append('so that')
        if defs:
            for curdef in defs:
                curlhs,currhs=curdef.split('=',1)
                subname=curlhs+'sub'
                maximalines.append('\\begin{maxima}')
                maximalines.append('\t\\parseopts{lhs="'+curlhs+'"}')
                maximalines.append(subname+':at('+currhs+','+substr+'])')
                maximalines.append('\\end{maxima}')
                cursub=curlhs+'='+subname
                substr+=','+cursub
                sublines.append(cursub)
        substr+=']'


    def GetDefs(self, **kwargs):
        """This function is a stop gap for derived elements that don't
        have their own GetDefs method.  Ideally, each derived element
        should override this method.  This is a first step in breaking
        up the functionality of GetMaximaLatexString into more useful
        and smaller functions.  For now, this function just calls
        GetMaximaLatexString and returns one portion of its return
        arguments."""
        maximalines, defs, params, sublines = self.GetMaximaLatexString(**kwargs)
        return defs


    def GetParams(self, **kwargs):
        """Another stop gap function to break up the functionality of
        GetMaximaLatexString.  This one should work for all elements
        with good __init__ methods that correctly handle symlabel."""
        myparams = self.symparams.values()
        mylist = rwkmisc.flatten(myparams)
        badlist=['0','1']
        outlist = [item for item in mylist if item not in badlist]
        return outlist


    def GetSubs(self, **kwargs):
        """Another stop gap function to break up GetMaximaLatexString.
        This one is sort of a hack.  See GetDefs."""
        maximalines, defs, params, sublines = self.GetMaximaLatexString(**kwargs)
        return sublines
                

    def GetMaximaLatexString(self,name=None,label=None,wrap=0,N=None,intro=None,aug=0,optargs={}):#omit
        maximalines=[]
        if intro:
            maximalines.append(intro)
        maximalines.append('\\begin{maxima-noout}')
        if aug:
            templines,myname,defs,params,label=self.GetAugMaximaString(name=name,label=label,N=N,**optargs)
        else:
            templines,myname,defs,params,label=self.GetMaximaString(name=name,label=label,N=N,**optargs)
        outlines=['\t'+line for line in templines]
        maximalines.extend(outlines)
        maximalines.append('\\end{maxima-noout')
        maximalines.append('\\begin{maxima}')
        maximalines.append("\t\\parseopts{lhs='"+myname+"',wrap="+str(wrap)+"}")
        maximalines.append('\t'+myname)
        maximalines.append('\\end{maxima}')
        if defs:
            maximalines.append('where')
            for curdef in defs:
                maximalines.append('\\begin{maxima}')
                maximalines.append('\t'+curdef)
                maximalines.append('\\end{maxima}')
        if self.symsub:
            maximalines.append('The following substitutions will be made:')
            sublines=[]
            substr='['
            first=1
            for cp,value in self.params.iteritems():
                if first:
                    first=0
                else:
                    substr+=','
                cursub=cp+label+'='+str(value)
                sublines.append(cursub)
                substr+=cursub
                maximalines.append('\\begin{equation}')
                maximalines.append('\t'+cursub)
                maximalines.append('\\end{equation}')
            paramsubs=substr+']'
            maximalines.append('so that')
            if defs:
                for curdef in defs:
                    curlhs,currhs=curdef.split('=',1)
                    subname=curlhs+'sub'
                    maximalines.append('\\begin{maxima}')
                    maximalines.append('\t\\parseopts{lhs="'+curlhs+'"}')
                    maximalines.append(subname+':at('+currhs+','+substr+'])')
                    maximalines.append('\\end{maxima}')
                    cursub=curlhs+'='+subname
                    substr+=','+cursub
                    sublines.append(cursub)
            substr+=']'
            return maximalines, defs, params, sublines
        else:
            return maximalines, defs, params,[]


class TMMElementLHT(TMMElement):
    """This class is meant to be a base class for all elements whose
    homogenous transformation matrix performs only a translation along
    the x-axis by the system parameter 'L' - i.e. that elements do not
    rotate the coordinate systemused for visualization.  This class
    should be used as the base class for beam and rigid link elements.
    It is a convenience class - it saves the author of new TMMElement
    classes from having to define the GetHT method for theis new
    element."""
    def GetHT(self):
        """Return the homogenous transformation matrix for the
        TMMElementLHT class.  This matrix will represent a translation
        along the 'x' axis with no rotation."""
        return HT4(x=self.params['L'])

class TMMElementIHT(TMMElement):
    """This is another convenience class meant to be a base class for
    torsional springs and other elments whose homogenous
    transformation matrix is simpliy and identity matrix."""
    def GetHT(self):
        """Return a 4x4 identifiy matrix for the homogenous
        transformation matrix of the TMMElementIHT class."""
        return HT4()

def HT4(axis='',angle=0,x=0.,y=0.,z=0.):
    """This function returns a 4x4 homogenous transformation matrix.
    Axis is either 'x', 'y', or 'z' (strings).  If axis is empty (if
    bool(axis) if false), there is no rotation portion to the matrix
    and the upper left 3x3 will be an identity matrix.  angle is in
    degrees.  x, y, and z are the translations associated with the
    homogenous transformation matrix."""
    if axis:
        r1=rot3by3(axis,angle)
    else:
        r1=scipy.eye(3,dtype='d')
    s1=c_[r1,array([[x],[y],[z]])]
    return r_[s1,array([[0.,0.,0.,1.]])]

def DH(d,theta,a,alpha):
    """Return a Denavit-Hartenburg homogenous transformation matrix
    according to the conventions of Sciavicco and Siciliano, second
    edition, p. 45 This means translating by d about the z_{i-1} and
    rotated by theta about that same axis.  Then translate by a along
    the x' axis and rotate by alpha about the x' axis.  Alpha and
    theta are both in degrees,"""
    r1=rot3by3('z',theta)
    s1=c_[r1,array([[0.],[0.],[d]])]
    s1=r_[s1,array([[0.,0.,0.,1.]])]
    r2=rot3by3('x',alpha)
    s2=c_[r2,array([[a],[0.],[0.]])]
    s2=r_[s2,array([[0.,0.,0.,1.]])]
    return scipy.dot(s1,s2)

def rot3by3(axis,angle):
    """Generate a 3x3 rotation matrix rotating angle degrees about
    axis."""
    rad=angle*pi/180.
    if isinstance(axis,str):
        if axis.lower()=='x':
            axis=1
        elif axis.lower()=='y':
            axis=2
        elif axis.lower()=='z':
            axis=3
    if axis==1:
        R=array([[1.0,0.0,0.0],[0.0,cos(rad),-sin(rad)],[0.0,sin(rad),cos(rad)]])
    elif axis==2:
        R=array([[cos(rad),0.0,sin(rad)],[0.0,1.0,0.0],[-sin(rad),0.0,cos(rad)]])
    elif axis==3:
        R=array([[cos(rad),-sin(rad),0.0],[sin(rad),cos(rad),0.0],[0.0,0.0,1.0]])
    return R

def rotationmat(axis,angle):
    """Create a 12x12 rotation matrix for use with 3D, 12x12 TMM
    analysis."""
    R=rot3by3(axis,angle)
    firstrow=1
#    im=scipy.eye(4,'f')
    im=scipy.eye(4,dtype='d')
    for currow in R.tolist():
        firstcol=1
        for ent in currow:
            curmat=im*ent
            if firstcol:
                myrow=curmat
                firstcol=0
            else:
                myrow=c_[myrow,curmat]
        if firstrow:
            mymat=myrow
            firstrow=0
        else:
            mymat=r_[mymat,myrow]
    return mymat

def Transform8by8(matin):
    """This function takes in an 8x8 matrix and switches rows and
    columns to unshuffle the mixed states associated with 2 and 3D
    transfer matrices."""
    matin=ColSwap(matin,0,4)
    matin=ColSwap(matin,3,7)
    matin=RowSwap(matin,0,4)
    matin=RowSwap(matin,3,7)
    return matin
