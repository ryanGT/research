from __future__ import division
from scipy import cos, cosh, sin, sinh, array, r_, c_, exp, pi, dot, real, imag, zeros, eye, alltrue, shape, atleast_1d
from scipy.linalg import inv as inverse
import MLab
import scipy
from scipy.linalg import det
import copy, os
import pdb
from rref import rrefnull
from rwkmisc import RowSwap, ColSwap, rowwise, colwise
import rwkascii
from rwkdataproc import datastruct
import rwkbode
import time, types
from textfiles.latexlist import latexlist
from textfiles.latexlist import ReplaceList
from textfiles import rwkreadfile, textlist
#from rwkmisc import null, norm2

class TMMElement:
    def __init__(self,elemtype,params,maxsize=12,label='',symname='',symlabel='',symsub=False):
        self.elemtype=self.parsetype(elemtype)
        self.params=params
        self.maxsize=maxsize
        self.label=label
        self.symlabel=symlabel
        self.symname=symname
        self.symsub=symsub
        if isinstance(elemtype,str):
            self.elemstr=elemtype
        else:
            self.elemstr=''

    def parsetype(self, elemtype):
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


    def GetMat(self, s):
        raise NotImplementedError

    def GetMaximaString(self,name=None,label=None,N=None):
        raise NotImplementedError

    def GetAugMaximaString(self,name=None,label=None,N=None):
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


    def GetMaximaLatexString(self,name=None,label=None,wrap=0,N=None,intro=None,aug=0,optargs={}):
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
            print('substr='+substr)
#            maximalines.append('\\begin{maxima}')
#            maximalines.append('\t\\parseopts{lhs="'+myname+'"}')
#            maximalines.append('\t'+myname+':at('+myname+','+substr+')')
#            maximalines.append('\\end{maxima}')
            return maximalines, defs, params, sublines
        else:
            return maximalines, defs, params,[]

    def GetAugMat(self,s):
        N=self.maxsize
#        print('s='+str(s))
#        print('elemstr='+self.elemstr)
#        tempout=eye(N+1,'d')+0.0j
        tempout=eye(N+1,dtype='D')
        tempout[0:N,0:N]=self.GetMat(s)
        return tempout

    def GetHT(self):
        raise NotImplementedError

#def rigidmat(s,params,dof):
#    if isinstance(dof,str):
#        dof=dof.lower()
#        mydict={'x':0,'y':1,'z':2}
#        dof=mydict[dof]
#    L=params['EndVector'][dof]
#    h=L-params['cg'][dof]
#    m=params['m']
#    I=params['Ivect'][dof]
#    matout=array([[1,L,0,0],[0,1,0,0],[h*m*(-s**2),(h*L-h**2)*(-s**2)*m+I*s**2,1,L],[m*(-s**2),m*(-s**2)*(L-h),0,1]])
#    return matout

def rigidmatx(s,params):
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

def rigidmatz(s,params):
    L=params['L']
    r=params['r']
    m=params['m']
    Iz=params['I'][2]
    #this matrix assumes +w_y is the state [w_y; theta_z; M_z; V_y] and is derived in 
    #Rmatrix_3D.m in E:\GT\Research\SAmii\modeling\transfer_matrix\threeD_TMM_derivations
#    [                    1,                    L,                    0,                    0]
#    [                    0,                    1,                    0,                    0]
#    [         -m*s^2*(L-r), s^2*Iz-m*s^2*r*(L-r),                    1,                   -L]
#    [                m*s^2,              m*s^2*r,                    0,                    1]
    matout=array([[1.,L,0,0],[0,1.,0,0],[-m*s**2*(L-r),s**2*Iz-m*s**2*r*(L-r),1.,-L],[m*s**2,m*s**2*r,0,1.]])
    return matout

def rigidmaty(s,params):
    L=params['L']
    r=params['r']
    m=params['m']
    Iy=params['I'][1]
    #this matrix assumes +w_z is the state [w_z; theta_y; M_y; V_z] and is derived in 
    #Rmatrix_3D.m in E:\GT\Research\SAmii\modeling\transfer_matrix\threeD_TMM_derivations
#    [                    1,                   -L,                    0,                    0]
#    [                    0,                    1,                    0,                    0]
#    [          m*s^2*(L-r), s^2*Iz-m*s^2*r*(L-r),                    1,                    L]
#    [                m*s^2,             -m*s^2*r,                    0,                    1]
    matout=array([[1.,-L,0,0],[0,1.,0,0],[m*s**2*(L-r),s**2*Iy-m*s**2*r*(L-r),1.,L],[m*s**2,-m*s**2*r,0,1.]])
    return matout

def HT4(axis='',angle=0,x=0.,y=0.,z=0.):
    #4x4 homogenous transformation matrix
    if axis:
        r1=rot3by3(axis,angle)
    else:
#        r1=scipy.eye(3,'f')
        r1=scipy.eye(3,dtype='f')
    s1=c_[r1,array([[x],[y],[z]])]
    return r_[s1,array([[0.,0.,0.,1.]])]

def DH(d,theta,a,alpha):
    #Denavit-Hartenburg homogenous transformation matrix according to the conventions of Sciavicco and Siciliano, second edition, p. 45
    r1=rot3by3('z',theta)
    s1=c_[r1,array([[0.],[0.],[d]])]
    s1=r_[s1,array([[0.,0.,0.,1.]])]
    r2=rot3by3('x',alpha)
    s2=c_[r2,array([[a],[0.],[0.]])]
    s2=r_[s2,array([[0.,0.,0.,1.]])]
    return scipy.dot(s1,s2)

def rot3by3(axis,angle):
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
    R=rot3by3(axis,angle)
    firstrow=1
#    im=scipy.eye(4,'f')
    im=scipy.eye(4,dtype='f')
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


#def bendmat(s,params,EIstr='EI1'):
#    EI=params[EIstr]
#    mu=params['mu']
#    L=params['L']
#    beta=pow((-1*pow(s,2)*pow(L,4)*mu/EI),0.25)
#    a=pow(L,2)/EI
#    c0=(cosh(beta)+cos(beta))/2.0
#    c1=(sinh(beta)+sin(beta))/(2.0*beta)
#    c2=(cosh(beta)-cos(beta))/(2.0*pow(beta,2));
#    c3=(sinh(beta)-sin(beta))/(2*pow(beta,3));
#    matout=array([[c0, L*c1, a*c2, a*L*c3], [pow(beta,4)*c3/L, c0, a*c1/L, a*c2],[pow(beta,4)*c2/a, pow(beta,4)*L*c3/a, c0, L*c1],[pow(beta,4)*c1/(a*L), pow(beta,4)*c2/a, pow(beta,4)*c3/L, c0]])
#    return matout


def Transform8by8(matin):
    matin=ColSwap(matin,0,4)
    matin=ColSwap(matin,3,7)
    matin=RowSwap(matin,0,4)
    matin=RowSwap(matin,3,7)
    return matin

class TMMSystem:
    def __init__(self,elemlist,bcend=[2,3,6,7],bcbase=[0,1,4,5],bodeouts=[]):
        self.sysHT=[]
        self.elemlist=elemlist#a list of TMMElement instances
        self.maxsize=elemlist[0].maxsize
        self.bcend=self._parsebcstr(bcend)#a list of the boundary condition indices that are 0 at the end of the chain (the free end in the case of a cantilever configuration) - the defaults are for a 1D cantilever where the states are [-w, phi, M, V] and M and V are 0 at the free end
        self.bcbase=self._parsebcstr(bcbase)#a list of the boundary configuration indices that are 0 as the base of the chain (the clamped beginning of a cantilever configuration) - the defaults are for a 1D cantilever where -w and phi are 0 at the base
        templist=range(self.maxsize)
        self.bccols=[]
        for ent in templist:
            if ent not in self.bcbase:
                self.bccols.append(ent)#in evaluating the eigvalues, the columns of the system U matrix that multiply the non-zero elements of the base state vector are the ones whose sub-determinant must go to 0

#        pdb.set_trace()
        if bodeouts:
            bodeouts=atleast_1d(bodeouts)
            tempbo=[]
            for cb in bodeouts:
                if type(cb)==dict:
                    curtemp=rwkbode.bodeout(**cb)
                    tempbo.append(curtemp)
                else:
                    tempbo.append(cb)
        self.bodeouts=tempbo
#        self.bodeoutlist=[]
#        self.dofdict={}
        self.rawbodeouts=[]
        for curout in self.bodeouts:
            if not isinstance(curout,rwkbode.bodeout):
                curout=curout.toscalar()
            scalar=0
            if not shape(curout.ind):
                scalar=1
            myouts=[]
            curents=atleast_1d(curout.ind)
            for curent in curents:
                if not isinstance(curent,rwkbode.bodeout):
                    curent=curent.toscalar()
                if isinstance(curent,TMMElement):
                    ind=elemlist.index(curent)
                elif type(curent)==types.IntType:
                    ind=curent
                else:
                    raise ValueError, "Unknown type for Bode output:"+str(curent)     
                myouts.append(ind)
#                self.dofdict[ind]=curout['dof']
#            pdb.set_trace()
            for tempout in myouts:
                tempraw=copy.deepcopy(curout)
                tempraw.pind=tempout
                self.rawbodeouts.append(tempraw)
            if scalar:
                curout.pind=myouts[0]
            else:
                curout.pind=myouts
#            self.bodeoutlist.extend(myouts)
#        self.bodeoutlist.sort()

    def SymBodesAllinOne(self,texname='SymBodeDev.tex',frlist=[],basebodename='bode',basecdname='chardet',functionparams={},otherparams={},debug=False):
        print('In SymBodesAllinOne')
        #,defaultvalues={}):
        aug=1
        mylist=latexlist([],texname)
        alldefs=[]
        allparams=[]
        defaultvalues={}
        allsubs=[]
#        pdb.set_trace()
        for elem in self.elemlist:
            curlines, curdefs, curparams,cursubs=elem.GetMaximaLatexString(aug=aug)
            mylist.extend(curlines)
            alldefs+=curdefs
            allparams+=curparams
            allsubs+=cursubs
            for key,value in elem.params.iteritems():
                defaultvalues[key+elem.symlabel]=value
#        pdb.set_trace()
        print('allsubs:'+str(allsubs))
        ratparams=copy.deepcopy(allparams)
        for curdef in alldefs:
            lhs,junk=curdef.split('=',1)
            ratparams.append(lhs)
        ratstr='ratvars('
        first=1
        for ent in ratparams:
            if first:
                first=0
            else:
                ratstr+=','
            ratstr+=ent
        ratstr+=',s)'
        preamble=[]
        preamble.append('\\begin{maxima-noout}')
        preamble.append('\tshowtime:all')
        preamble.append('\tratmx:true')
        preamble.append('\tnolabels:true')
        preamble.append('\t'+ratstr)
        preamble.append('\\end{maxima-noout}')
#        pdb.set_trace()
        mylist.list=preamble+mylist.list
#        if debug:
#            mylist.tofile()
#            return mylist
        syslines=self.SymULatex()
        mylist.extend(syslines)
        submatlines=self.SymSubmat(findsubcol=True)
        mylist.extend(submatlines)
        bvlines=self.SymBaseVect()
        mylist.extend(bvlines)
        rblines=self.SymRawBodes(frlist)
        mylist.extend(rblines)
        bodelines,filenames,chardetnames=self.SymBodes(basebodename=basebodename, basecdname=basecdname,debug=debug)
        mylist.extend(bodelines)
#        return mylist,filenames,chardetnames
        mylist.tofile()
        if debug:
            mylist.tofile()
            return mylist
        print('about to ProcessMaxima')
        outname=mylist.ProcessMaxima()
        pyfilenames=[]
        for key,value in otherparams.iteritems():
            defaultvalues[key]=value
        for curfile, curcdfile in zip(filenames,chardetnames): 
            fno,junk=os.path.splitext(curfile)
            cdfno,junk=os.path.splitext(curcdfile)
            pyfilename=fno+'.py'
            pyfilenames.append(pyfilename)
            cdin=rwkreadfile(curfile)
            cdin2=rwkreadfile(curcdfile)
            cdout=[]
            cdout2=[]
            cdlines=[]
            for line in cdin:
                if line:
                    cdout.append(line)
            for line in cdin2:
                if line:
                    cdout2.append(line)
            rv,dummy=cdout[0].split('=',1)
            cdout.append('return '+rv)
            rv2,dummy=cdout2[0].split('=',1)
            cdout2.append('return '+rv2)
            pylines=textlist([],pyfilename)
            pylines.readfile('header.py')
            funcstr='s'
            for key,value in functionparams.iteritems():
                funcstr+=','+str(key)+'='+str(value)
            pylines.append('def '+fno+'('+funcstr+'):')
            cdlines.append('def '+cdfno+'('+funcstr+'):')
            for cp in allparams:
                if cp not in functionparams.keys():
                    p1='\t'+cp+'='
                    if defaultvalues.has_key(cp):
                        tempval=defaultvalues[cp]
                        if shape(tempval):
                            tempval=tempval[0]
                        p1+=str(tempval)
                    pylines.append(p1)
                    cdlines.append(p1)
            for curd in alldefs:
                curd=curd.replace('^','**')
                pylines.append('\t'+curd)
                cdlines.append('\t'+curd)
            for line in cdout:
                line=line.replace('^','**')
                pylines.append('\t'+line)
            for line in cdout2:
                line=line.replace('^','**')
                cdlines.append('\t'+line)
            pylines.append('')
            pylines.extend(cdlines)
            pylines.tofile()
        return mylist,filenames,chardetnames, pyfilenames

    def SymBodes(self,basebodename='bodeoutputs',basecdname='chardet',myext='.txt',debug=False):
        outfiles=[]
        linesout=[]
        filenames=[]
        chardetnames=[]
        for x,cb in enumerate(self.bodeouts):
            scalar=1
            bodestr='bode'+str(x)
            curname=basebodename+str(x)+myext
            curcdname=basecdname+str(x)+myext
            filenames.append(curname)
            chardetnames.append(curcdname)
            linesout.append('Bode output '+str(x) +' will be given by')
            if shape(cb.pind):
                if max(shape(cb.pind))>1:
                    scalar=0
            if scalar:
                curbodesym=self.FindSymRawBode(cb.pind,cb.dof)
                bodeeq=bodestr+'='+curbodesym
                if cb.post.lower()=='vel':
                    bodeeq+='*s'
                elif cb.post.lower()=='accel':
                    bodeeq+='*s^2'
            else:
                b1=self.FindSymRawBode(cb.pind[0],cb.dof)
                b2=self.FindSymRawBode(cb.pind[1],cb.dof)
                rhs=b1+'-'+b2
                if cb.post.lower()=='vel':
                    rhs='('+rhs+')*s'
                elif cb.post.lower()=='accel':
                    rhs='('+rhs+')*s^2'
                bodeeq=bodestr+'='+rhs
            linesout.append('\\begin{equation}')
            linesout.append('\t'+bodeeq.replace('*', ' '))
            linesout.append('\\end{equation}')
            linesout.append('\\begin{maxima-noout}')
            linesout.append(bodeeq.replace('=',':',1))
            linesout.append('\\end{maxima-noout}')
            if not debug:
                linesout.append('\\begin{maxima-noout}')
    #            linesout.append('\tstringout("'+curname+'",bode=radcan('+bodestr+'))')
                linesout.append('\tstringout("'+curname+'",bode='+bodestr+')')
                linesout.append('\tstringout("'+curcdname+'",chardet=ratdenom('+bodestr+'))')
                linesout.append('\\end{maxima-noout}')
            else:
                print('skipping stringouts')
        return linesout, filenames, chardetnames

    def FindSymRawBode(self,pind,dof):
        """This function returns the symbolic
        name for a raw bode output matching a
        specific location and dof.  It is assumed
        that you have already run SymRawBodes to
        calculate the rawbode outputs and attach
        their symbolic names."""
        for rb in self.rawbodeouts:
            if rb.pind==pind and rb.dof==dof:
                return rb.symname
        return None

    def SymRawBodes(self,frlist=[]):
        """This function calculates the raw bode outputs
        for a TMMSystem.  It assumes that the base vector bv
        is already defined, so you must have already called 
        SymBaseVect.  This function returns a LaTeX list 
        that calculates the matrices U_i that transfer 
        from the base vector to the raw bode locations 
        as well as the raw bode outputs themselves
        where

        rb_i=U_i*bv"""

        linesout=[]
#        modbodes=[]
        for x,rb in enumerate(self.rawbodeouts):
            curp=rb.pind
            pstr=str(curp)
            curintro='The transfer matrix from the base vector to model position '+pstr+' is'
            curUname='U'+pstr
            curUlines=self.SymULatex(intro=curintro,stopind=curp,sysmatname=curUname,frlist=frlist)
            linesout.extend(curUlines)  
            linesout.append('The state vector at this location is given by')
            linesout.append('\\begin{equation}')
            linesout.append('\t\\M{z}_{'+pstr+'}='+curUname+' \\M{bv}')
            linesout.append('\\end{equation}')
            cursymname='rb'+str(x)
            linesout.append('\\begin{maxima-noout}')
#            linesout.append("\t\\parseopts{lhs='"+cursymname+"'}")
            linesout.append('\ttemprow:row('+curUname+','+str(rb.dof+1)+')')
            linesout.append('\t'+cursymname+':temprow'+'.bv')
            linesout.append('\\end{maxima-noout}')
            rb.symname=cursymname
#            linesout.append('\\begin{maxima-noout}')
#            linesout.append("\t\\parseopts{lhs='"+cursymname+"',wrap=o}")
#            linesout.append('\t'+cursymname+':radcan('+cursymname+'['+str(rb.dof+1)+',1])')
#            linesout.append('\t'+cursymname+':'+cursymname+'['+str(rb.dof+1)+',1]')
#            linesout.append('\\end{maxima-noout}')
#            modbodes.append(rb)
#        self.rawbodeouts=modbodes
        return linesout

    def SymBaseVect(self,N=None,showall=False):
        """This function finds the base vector
        used to calculate the system bode response.
        It assumes that SymSubmat had already been 
        called and will append lines using submat and
        subcol, assuming they have already been defined
        by SymSubmat called with findsubcol=True.  

        This funciton returns lines to be append to a 
        LaTeX list.  These lines create a symbolic solution
        for the base vector."""
        if N is None:
            N=self.maxsize
        linesout=[]
        linesout.append("The non-zero portion of the base vector will be found using")
        linesout.append('\\begin{equation}')
        linesout.append('\\M{nzbv}=\\M{submat}^{-1}\\left(-\\M{subcol}\\right)')
        linesout.append('\\end{equation}')
        linesout.append('or')
        if showall:
            linesout.append('\\begin{maxima}')
            linesout.append("\t\\parseopts{lhs='\\M{nzbv}'}")
    #        linesout.append('\tnzbv:radcan(invert(submat).(-1*subcol))')
            linesout.append('\tnzbv:invert(submat).(-1*subcol)')
            linesout.append('\\end{maxima}')
        else:
            linesout.append('\\begin{maxima-noout}')
    #        linesout.append('\tnzbv:radcan(invert(submat).(-1*subcol))')
            linesout.append('\tnzbv:invert(submat).(-1*subcol)')
            linesout.append('\\end{maxima-noout}')
        linesout.append('\\begin{maxima-noout}')
        linesout.append('bv:zeromatrix('+str(N+1)+',1)')
        linesout.append('bv['+str(N+1)+',1]:1')
        for x,r in enumerate(self.bcend):
            linesout.append('bv['+str(r+1)+',1]:nzbv['+str(x+1)+',1]')
        linesout.append('bv:radcan(bv)')
        linesout.append('\\end{maxima-noout}')
        if showall:
            linesout.append('\\begin{maxima}')
            linesout.append("\t\\parseopts{lhs='\\M{bv}'}")
            linesout.append("\tbv")
            linesout.append('\\end{maxima}')
        return linesout

    def SymSubmat(self, sysmatname='Usys',findsubcol=False,N=None,showsub=False):
        """Note that you must call SymULatex first. 
        This function will return lines that are to
        be appended to a latexlist that already includes
        the SymULatex lines.  This function uses the 
        system boundary conditions to find a sub-matrix
        whose determinant is the characteristic equation
        of the system."""
        if N is None:
            N=self.maxsize
        keeprows=self.bcend
        keepcols=self.bccols
        delrows=[]
        delcols=[]
        if findsubcol:
            N2=N+1
        else:
            N2=N
        for x in range(N2):
            if x not in keeprows:
                delrows.append(x+1)
            if x not in keepcols:
                delcols.append(x+1)
        submatstr='submat:submatrix('
        subcolstr='subcol:submatrix('
        first=1
        for cr in delrows:
            if first:
                first=0
            else:
                submatstr+=','
                subcolstr+=','
            submatstr+=str(cr)
            subcolstr+=str(cr)
        submatstr+=','+sysmatname
        subcolstr+=','+sysmatname
        for cc in delcols:
            submatstr+=','+str(cc)
        for x in range(N):
            subcolstr+=','
            subcolstr+=str(x+1)
        submatstr+=')'
        subcolstr+=')'
        linesout=[]
        linesout.append('Based on the system boundary conditions, the submatrix whose determinant is the characteristic equation is')
        if showsub:
            linesout.append('\\begin{maxima}')
            linesout.append("\t\\parseopts{lhs='\M{subU}',wrap=0}")
            linesout.append('\t'+submatstr)
            linesout.append('\\end{maxima}')
        else:
            linesout.append('\\begin{maxima-noout}')
            linesout.append('\t'+submatstr)
            linesout.append('\\end{maxima-noout}')
            linesout.append('\\begin{equation}')
            linesout.append('\t'+submatstr.replace(':','=',1))
            linesout.append('\\end{equation}')
        if findsubcol:
            linesout.append('and the sub-column that will be used to calculated the Bode response is')
            linesout.append('\\begin{maxima}')
            linesout.append("\t\\parseopts{lhs='\M{subcol}',wrap=0}")
            linesout.append('\t'+subcolstr)
            linesout.append('\\end{maxima}')
        return linesout

    def SymCharDet(self, sysmatname='Usys',outpath='chardet.txt'):
        """This function calls SymSubmat to find the
        submatrix and then calculates its determinant
        to find the characteristic determinant of the
        system.

        Note that you must call SymULatex first."""
        linesout=self.SymSubmat(sysmatname=sysmatname)
        linesout.append('and the characteristic determinant is')
        linesout.append('\\begin{maxima}')
        linesout.append("\t\\parseopts{lhs='cd',wrap=0}")
        linesout.append('chardet:determinant(submat)')
        linesout.append('\\end{maxima}')
        linesout.append('\\begin{maxima-noout}')
        linesout.append('\tstringout("'+outpath+'",cd=chardet)')
        linesout.append('\\end{maxima-noout}')
        return linesout
                    
    def SymULatex(self,sysmatname='Usys',frlist=[],stopind=None,intro=None,showU=False):
        """This function returns lines for a LaTeX list 
        that calculate the system transfer matrix.
        By setting stopind, it can also be used to 
        calculate transfer matrices that calculate
        bode outputs at certain locations based on
        the base vector."""
        Unames=[]
        prodstr=''
        first=1
        for x,curelem in enumerate(self.elemlist):
            if first:
                first=0
            else:
                prodstr='.'+prodstr
            allouts=curelem.GetMaximaString()
            curname=allouts[1]
            Unames.append(curname)
            prodstr=curname+prodstr
            if stopind is not None:
                if x==stopind:
                    break
        latexout=[]
        if intro is None:
            intro='The system transfer matrix will be'
        latexout.append(intro)
        latexout.append('\\begin{equation}')
        templine=sysmatname+'='+prodstr.replace('.',' ')
        if frlist:
            templine=ReplaceList(templine,frlist)
        latexout.append('\t'+templine)
        latexout.append('\\end{equation}')
#        latexout.append('\\begin{maxima-noout}')
#        latexout.append('\\end{maxima-noout}')
        if showU:
            latexout.append('or')
            latexout.append('\\begin{maxima}')
    #        latexout.append('\t'+sysmatname+':radcan('+sysmatname+')')
            latexout.append("\parseopts{lhs='"+sysmatname+"',wrap=0}")
            latexout.append('\t'+sysmatname+':'+prodstr)
            latexout.append('\\end{maxima}')
        else:
            latexout.append('\\begin{maxima-noout}')
    #        latexout.append('\t'+sysmatname+':radcan('+sysmatname+')')
            latexout.append('\t'+sysmatname+':'+prodstr)
            latexout.append('\\end{maxima-noout}')
        return latexout

    def _parsebcstr(self, bcstr):
        VMlist=[[2,3],[6,7],[10,11]]
        XTlist=[[0,1],[4,5],[8,9]]
        if isinstance(bcstr,str):
            if bcstr=='free':
                bcelist=VMlist
            elif bcstr=='fixed':
                bcelist=XTlist
            else:
                raise ValueError, 'Boundary Condition String not recognized:\''+bcstr+'\''+'\nValid entries are "free" and "fixed".'
            bcout=bcelist[0]
            if self.maxsize>4:
                bcout.extend(bcelist[1])
            if self.maxsize>8:
                bcout.extend(bcelist[2])
            return bcout
        else:
            return bcstr

    def CreateSysHT(self, beammesh=10):
#        pdb.set_trace()
#        ovect=array([[0],[0],[0],[1]])
        self.sysHT=[]
        self.sysR=[]
#        prevHT=scipy.eye(4,'f')
        prevHT=scipy.eye(4,dtype='f')
        for elem in self.elemlist:
#            if elem.elemtype==1:
#                beamHT=[]
#                L=elem.params['L']
#                ce=copy.deepcopy(elem)
#                for x in range(beammesh):
#                    myL=((x+1.0)*L)/beammesh
#                    ce.params['L']=myL
#                    relHT=ce.GetHT()
#                    curBHT=dot(prevHT,relHT)
#                    beamHT.append(curBHT)
#                self.sysHT.append(beamHT)
#                prevHT=curBHT
#            else:
#                relHT=elem.GetHT()
#                curHT=dot(prevHT,relHT)
#                self.sysHT.append(curHT)
#                prevHT=curHT
            relHT=elem.GetHT()
            curHT=dot(prevHT,relHT)
            self.sysHT.append(curHT)
            self.sysR.append(curHT[0:3,0:3])
            prevHT=curHT
        return self.sysHT

    def CreateMesh(self, beammesh=10):
#        pdb.set_trace()
        self.sysmesh=[]
        basemesh=datastruct()
        basemesh['elemtype']=-1
#        basemesh['R']=eye(3,typecode='d')
        basemesh['R']=eye(3,dtype='d')
        basemesh['mesh']=zeros((1,3),'f')
        self.sysmesh.append(basemesh)
#        ne=len(self.elemlist)
        if not self.sysHT:
            self.CreateSysHT(beammesh)
        mesh=scipy.zeros((1,3),'f')
        for elem, curHT, curR in zip(self.elemlist,self.sysHT, self.sysR):
            curmesh=datastruct()
            curmesh['elemtype']=elem.elemtype
            curmesh['R']=curR
            curmesh['mesh']=[]
            cmfirst=1
            if elem.elemtype==1:
#                pdb.set_trace()
                L=elem.params['L']
                tempHT=copy.deepcopy(curHT)
                tempHT[0:3,3]=zeros((3,),'f')
                for x in range(beammesh):
                    curL=(x+1.0)*L/beammesh
                    myvect=array([curL,0.,0.,1.])
                    myvect=colwise(myvect)
                    mymesh=dot(tempHT,myvect)
                    mymesh=rowwise(mymesh[0:3])
                    mesh=r_[mesh,mymesh]
                    if cmfirst:
                        curmesh['mesh']=mymesh
                        cmfirst=0
                    else:
                        curmesh['mesh']=r_[curmesh['mesh'],mymesh]
                    
            else:
                mymesh=array([(curHT[0:3,3])])
                mesh=r_[mesh,mymesh]
                curmesh['mesh']=rowwise(mymesh)
            self.sysmesh.append(curmesh)
        return mesh

#    def MeshSystem(self, beammesh=10):
#        if not self.sysHT:
#            self.CreateSysHT(beammesh)
        
        
        

    def FindModeShape(self, eigval,beammesh=10,logtex=0,fmt='0.5g',scale=True):
        #,MV=False):
        modedict=datastruct()
        modedict['eigval']=eigval
        modedict['bodies']=[]
#        pdb.set_trace()
        if logtex:
            texlist=rwkascii.rwklatexlist([],'FindModeShape_log.tex')
            texlist.addheader()
            vectstr='\\left[\\begin{array}{c}-w_x \\\\ \\theta_x \\\\ M_x \\\\ V_x \\\\ -w_y \\\\ \\theta_y \\\\ M_y \\\\ V_y \\\\ -w_z \\\\ \\theta_z \\\\ M_z \\\\ V_z \\end{array} \\right]'
        if not self.sysHT:
            self.CreateSysHT(beammesh)
        if scipy.shape(eigval)==(2,):
            eigval=eigval[0]+eigval[1]*1j
        bv=self.FindBaseVect(eigval)
        if self.maxsize==4:
            disps=array([bv[0]])
            angles=array([bv[1]])
#            if MV:
#                moments=array([bv[2]])
#                forces=array([bv[3]])
        elif self.maxsize==8:
            disps=array([bv[0],bv[4]])
            angles=array([bv[1],bv[5]])
#            if MV:
#                moments=array([bv[2],bv[6]])
#                forces=array([bv[3],bv[7]])
        elif self.maxsize==12:
            disps=array([bv[0],bv[4],bv[8]])
            angles=array([bv[1],bv[5],bv[9]])
#            if MV:
#                moments=array([bv[2],bv[6],bv[10]])
#                forces=array([bv[3],bv[7],bv[11]])
        disps=rowwise(disps)
        angles=rowwise(angles)
#        if MV:
#            moments=rowwise(moments)
#            forces=rowwise(forces)
        curmode=datastruct()
        curmode['elemtype']=-1#this will indicate the base vector
        curmode['disps']=disps
        curmode['angles']=angles
#        if MV:
#            curmode['moments']=moments
#            curmode['forces']=forces
        modedict['bodies'].append(curmode)
        prevU=scipy.eye(self.maxsize,'f')+0j
        if logtex:
            texlist.append('Before main loop,')
            texlist.AppendEq('disps',disps,fmt=fmt)
            texlist.AppendEq('angles',angles,fmt=fmt)
#            if MV:
#                texlist.AppendEq('moments',moments,fmt=fmt)
#                texlist.AppendEq('forces',forces,fmt=fmt)
            texlist.AppendEq('basevect='+vectstr,bv,fmt=fmt)
            texlist.AppendEq('\M{U}',prevU,fmt=fmt)
        for elem, curR in zip(self.elemlist,self.sysR):
            curmode=datastruct()
            curmode['elemtype']=elem.elemtype
            curmode['disps']=[]
            curmode['angles']=[]
#            if MV:
#                curmode['moments']=[]
#                curmode['forces']=[]
            cmfirst=1
            if elem.elemtype==1:
                L=elem.params['L']
                ce=copy.deepcopy(elem)
                if logtex:
                    texlist.append('\\section{Beam}')
#                for x,curBHT in zip(range(beammesh),curHT):
                for x in range(beammesh):
#                    pdb.set_trace()
                    myL=(x+1.0)*L/beammesh
                    ce.params['L']=myL
                    elemU=ce.GetMat(eigval)
                    curU=dot(elemU,prevU)
                    curZ=dot(curU,bv)
                    curZ=curZ[:,0]
                    if self.maxsize==4:
                        tempdisps=colwise(array([curZ[0]]))
                        tempangles=array([curZ[1]])
#                        if MV:
#                            tempM=array([curZ[2]])
#                            tempV=array([curZ[3]])
                    elif self.maxsize==8:
                        tempdisps=colwise(array([curZ[0],curZ[4]]))
                        tempangles=array([curZ[1],curZ[5]])
#                        if MV:
#                            tempM=array([curZ[2],curZ[6]])
#                            tempV=array([curZ[3],curZ[7]])
                    elif self.maxsize==12:
                        tempdisps=colwise(array([curZ[0],curZ[4],curZ[8]]))
                        tempangles=array([curZ[1],curZ[5],curZ[9]])
#                        if MV:
#                            tempM=array([curZ[2],curZ[6],curZ[10]])
#                            tempV=array([curZ[3],curZ[7],curZ[11]])
                    if len(tempdisps)<3:
                        curdisps=tempdisps
                        curangles=tempangles
#                        if MV:
#                            curM=tempM
#                            curV=tempV
                    else:
                        curdisps=dot(curR,tempdisps)
                        curdisps=rowwise(curdisps)
                        tempangles=colwise(tempangles)
                        curangles=dot(curR,tempangles)
#                        if MV:
#                            tempM=colwise(tempM)
#                            curM=dot(curR,tempM)
#                            tempV=colwise(tempV)
#                            curV=dot(curR,tempV)
                    disps=r_[disps,curdisps]
                    curangles=rowwise(curangles)
                    angles=r_[angles,curangles]
                    if cmfirst:
                        cmfirst=0
                        curmode['disps']=curdisps
                        curmode['angles']=curangles
                    else:
                        curmode['disps']=r_[curmode['disps'],curdisps]
                        curmode['angles']=r_[curmode['angles'],curangles]
                    if logtex:
                        texlist.append('In beam loop,')
                        texlist.AppendEq('myL',myL)
                        texlist.AppendEq('curdisps_{body}',rowwise(tempdisps),fmt=fmt)
                        texlist.AppendEq('curdisps_{global}',rowwise(curdisps),fmt=fmt)
                        texlist.AppendEq('curangles_{body}',rowwise(tempangles),fmt=fmt)
                        texlist.AppendEq('curangles_{global}',rowwise(curangles),fmt=fmt)
#                        pdb.set_trace()
                        texlist.AppendEq('curZ='+vectstr,colwise(curZ),fmt=fmt)
                        texlist.AppendEq('elemU',elemU,fmt=fmt)
                modedict['bodies'].append(curmode)
                prevU=curU
            else:
#                pdb.set_trace()
                elemU=elem.GetMat(eigval)
                curU=dot(elemU,prevU)
                curZ=dot(curU,bv)
                curZ=curZ[:,0]
                if self.maxsize==4:
                    tempdisps=colwise(array([curZ[0]]))
                    tempangles=array([curZ[1]])
                elif self.maxsize==8:
                    tempdisps=colwise(array([curZ[0],curZ[4]]))
                    tempangles=array([curZ[1],curZ[5]])
                elif self.maxsize==12:
                    tempdisps=colwise(array([curZ[0],curZ[4],curZ[8]]))
                    tempangles=array([curZ[1],curZ[5],curZ[9]])
#                tempdisps=colwise(array([curZ[0],curZ[4],curZ[8]]))
                if len(tempdisps)<3:
                    curdisps=tempdisps
                    curangles=tempangles
                else:
                    curdisps=dot(curR,tempdisps)
                    curdisps=rowwise(curdisps[0:3,:])
                    curangles=dot(curR,tempangles)
                disps=r_[disps,curdisps]
#                tempangles=colwise(array([curZ[1],curZ[5],curZ[9]]))
                curangles=rowwise(curangles)
                angles=r_[angles,curangles]
                curmode['disps']=curdisps
                curmode['angles']=curangles
                if logtex:
                    texlist.append('\\section{element type='+str(elem.elemtype)+'}')
                    texlist.AppendEq('myL',myL)
                    texlist.AppendEq('curdisps_{body}',rowwise(tempdisps),fmt=fmt)
                    texlist.AppendEq('curdisps_{global}',rowwise(curdisps),fmt=fmt)
                    texlist.AppendEq('curangles_{body}',rowwise(tempangles),fmt=fmt)
                    texlist.AppendEq('curangles_{global}',rowwise(curangles),fmt=fmt)
#                        pdb.set_trace()
                    texlist.AppendEq('curZ='+vectstr,colwise(curZ),fmt=fmt)
                    texlist.AppendEq('curR',curR,fmt=fmt)
                    texlist.AppendEq('elemU',elemU,fmt=fmt)
                prevU=curU
                modedict['bodies'].append(curmode)               
#        pdb.set_trace()
        dmax=matmaxabs(real(disps))#*-1.#undo the -w TMM sign configuration for FEA comparison
        amax=matmaxabs(real(angles))
        dmaxi=matmaxabs(imag(disps))#*-1.#undo the -w TMM sign convention for FEA comparison
        amaxi=matmaxabs(imag(angles))
        if scale:
            tmax=matmaxabs(array([dmax,amax,dmaxi,amaxi]))
            disps=disps*0.2/tmax
            angles=angles*0.2/tmax
            for curmode in modedict['bodies']:
                curmode['disps']=curmode['disps']*0.2/tmax
                curmode['angles']=curmode['angles']*0.2/tmax
        if logtex:
            texlist.tofile()
        return disps,angles,modedict


    def FindBaseVect(self, eigval):
#        pdb.set_trace()
        submat=self.FindSubmat(eigval)
#        ns=null(submat)
        nstol=1e-8
        ns=rrefnull(submat,nstol)
        basevect=scipy.zeros((self.maxsize,1))#this will be the vector of boundary conditions at the base of the cantilever
        basevect=basevect*1j#force the vector to be cast into complex
#        print('ns[0,0]='+str(ns[0,0]))
#        print('type(ns[0,0])='+str(type(ns[0,0])))
#        print('shape(ns)='+str(scipy.shape(ns)))
        if len(scipy.shape(ns))>1:
            nsvect=ns[:,0]#*1e5#assume distinct eigenvalues/vectors and only 1D null space
        else:
            nsvect=ns
        for curind,curns in zip(self.bccols,nsvect):
            basevect[curind,0]=curns#insert non-zero values into the boundary condition vector at the base of the beam
#        endvect=scipy.dot(self.FindU(eigval),basevect)#this is the vector of boundary conditions at the free end of the cantilever
#        endvect=endvect/(scipy.linalg.norm(endvect))
        return basevect

    def BodeResponse(self,fvect):
        """Calculate the Bode response at system locations
        specified by self.bodeouts (specified during __init__
        of the TMMSystem).  fvect is a vector of frequencies
        at which the system response will be calculated.  Each
        entry in fvect is assumed to be a real number in Hz."""
#        Nout=len(self.bodeoutlist)
#        outmat=zeros((len(fvect),Nout),'d')
#        rawvect=zeros((len(fvect),Nout),'d')
#        outmat=outmat+0.j
        N=self.maxsize
#        rawbodedict={}
        for curitem in self.rawbodeouts:
            curitem.compvect=zeros(len(fvect),'d')+0.0j
        if not N%2:
            N+=1
#        pdb.set_trace()
        svect=2.0j*pi*array(fvect)
        for r,s in enumerate(svect):
#            U=scipy.eye(N,typecode='f')+0.j
            U=scipy.eye(N,dtype='D')
#            outcol=0
            bv=self.FindAugBaseVect(s)
            for x,curelem in enumerate(self.elemlist):
                tempU=curelem.GetAugMat(s)
                U=scipy.dot(tempU,U)
#                if x in self.bodeoutlist:
                curz=dot(U,bv)
                for curraw in self.rawbodeouts:
                    if curraw.pind==x:
#                        pdb.set_trace()
                        curraw.compvect[r]=curz[curraw.dof,0]#curz is 2d but has only 1 column - the 0 is to get a scalar from the only column
#                    outcol+=1
#        return outmat
#        return rawbodedict
        bodes=[]
        for bodedict in self.bodeouts:
            if not isinstance(bodedict,rwkbode.bodeout):
                bodedict=bodedict.toscalar()
            if bodedict.type.lower()=='diff':
#                col1=rawbodedict[bodedict.pind[0]]
#                col2=rawbodedict[bodedict.pind[1]]
                col1=self.FindRawBodeout(bodedict.pind[0],bodedict.dof)
                col2=self.FindRawBodeout(bodedict.pind[1],bodedict.dof)
                tempcomp=col1-col2
            else:
#                tempcomp=rawbodedict[bodedict.pind]
                tempcomp=self.FindRawBodeout(bodedict.pind,bodedict.dof)
            if bodedict.post.lower()=='accel':
                tempcomp=tempcomp*svect**2
            elif bodedict.post.lower()=='vel':
                tempcomp=tempcomp*svect
            bodes.append(rwkbode.rwkbode(bodedict.output,bodedict.input,compin=tempcomp))
        return bodes

    def FindRawBodeout(self,pind,dof):
        for curraw in self.rawbodeouts:
            if curraw.pind==pind and curraw.dof==dof:
                return curraw.compvect
        return -1

    def FindAugBaseVect(self,s):
        submat,augcol=self.FindAugSubmat(s)
        N=self.maxsize
        if not N%2:#assume an even maxsize means that the maxsize doesn't include the augmented column
            N+=1
        basevect=scipy.zeros((N,1))
        basevect=basevect+0.0j
        basevect[-1]=1
        tempvect=dot(inverse(submat),-augcol)
        for curind,curns in zip(self.bccols,tempvect):
            basevect[curind,0]=curns
        return basevect

    def FindEig(self, guess,mytol=1e-14,maxfun=1e4,maxiter=1e3,returncomplex=False,useabs=True):
#        t1=time.time()
#        eigout=scipy.optimize.fmin(self.EigError,guess,xtol=mytol,ftol=mytol,maxfun=maxfun,maxiter=maxiter)
        if useabs:
            if scipy.shape(guess)==():
                guess=scipy.array([scipy.real(guess),scipy.imag(guess)])
            eigout=scipy.optimize.fmin(self.EigError,guess,xtol=mytol,ftol=mytol,maxfun=maxfun,maxiter=maxiter)
        else:
            eigout=scipy.optimize.newton(self.EigError,guess,tol=mytol,maxiter=maxiter,args=(False,))
    #        eigout=scipy.optimize.fmin_bfgs(self.EigError,guess,gtol=1e-20)#,epsilon=1e-16)#,norm=2)
#            t2=time.time()
#            print('FindEig time='+str(t2-t1))
        if returncomplex:
            if shape(eigout):
                eigout=eigout[0]+eigout[1]*1.0j
        return eigout


    def EigError(self, value,useabs=True):
#        pdb.set_trace()
        if not shape(value):
            value=complex(value)
        submat=self.FindSubmat(value)
#        s = scipy.linalg.svdvals(submat)
#        return s[-1]/s[0] 
        shapevect=scipy.shape(submat)
        if shapevect:
            if len(shapevect)>1 and max(shapevect)>1:
                chardet=scipy.linalg.det(submat)
                if useabs:
                    return abs(chardet)
                else:
                    return chardet
#                s = scipy.linalg.svdvals(submat)
#                return s[-1]#/s[0] 
#                e=scipy.linalg.eigvals(submat)
#                e=abs(e)
#                return min(e)
            elif max(shapevect)==1:
                if len(shapevect)==2:
                    if useabs:
                        return abs(submat[0,0])
                    else:
                        return submat[0,0]
                elif len(shapevect)==1:
                    if useabs:
                        return abs(submat[0])
                    else:
                        return submat[0]
        else:
            return submat

#        mydet=det(submat)
###        print('mydet='+'%1.8e'+' '+'%1.8e'+'j') % (scipy.real(mydet),scipy.imag(mydet))
###        mye=scipy.sqrt(norm2(mydet))
###        print('mye='+str(mye))
#        mye=abs(mydet)
###        mye=pow(mye,0.1)
###        return norm2(mydet)
###        absval=scipy.absolute(mydet)
###        print('abs(mydet)='+str(absval))
#        return mye

    def FindU(self,value):
#        tu1=time.time()
        if scipy.shape(value):
            value=value[0]+value[1]*1j
#        templist=copy.deepcopy(self.elemlist)
#        curelem=templist.pop(0)
#        U=curelem.GetMat(value)
#        U=scipy.eye(self.maxsize,typecode='f')+0.j
        U=scipy.eye(self.maxsize,dtype='D')
        for curelem in self.elemlist:
            tempU=curelem.GetMat(value)
            U=scipy.dot(tempU,U)
#        tu2=time.time()
#        print('FindU time='+str(tu2-tu1))
        return U
    
    def FindAugU(self,value):
        if scipy.shape(value):
            value=value[0]+value[1]*1j
#        U=scipy.eye(self.maxsize+1,typecode='f')+0.j
        U=scipy.eye(self.maxsize+1,dtype='D')
        for curelem in self.elemlist:
            tempU=curelem.GetAugMat(value)
            U=scipy.dot(tempU,U)
        return U

    def FindSubmat(self,value,eps=1e-14):
        U=self.FindU(value)
        N=self.maxsize
        n=N/2
        matout=scipy.zeros((n,n),'f')+0j
        for r,curri in enumerate(self.bcend):
            for c,curci in enumerate(self.bccols):
                matout[r,c]=U[curri,curci]
#        pdb.set_trace()
        testmat=max(abs(imag(matout)))
        if shape(testmat):
            testmat=max(testmat)
        if testmat<=eps:
            matout=real(matout)
        return matout

    def FindSubmatwAugcol(self,value):
#        pdb.set_trace()
        U=self.FindAugU(value)
        N=self.maxsize
        n=N/2
        matout=scipy.zeros((n,n),'f')+0j
        colout=scipy.zeros((n,1),'f')+0j
        for r,curri in enumerate(self.bcend):
            for c,curci in enumerate(self.bccols):
                matout[r,c]=U[curri,curci]
            colout[r,0]=U[curri,N]
        return matout, colout

    def FindAugSubmat(self,value):
#        pdb.set_trace()
        U=self.FindAugU(value)
        N=self.maxsize
        n=int(N/2)
        matout=scipy.zeros((n,n),'d')+0j
        vectout=scipy.zeros((n,),'d')+0j
        for r,curri in enumerate(self.bcend):
            vectout[r]=U[curri,N]
            for c,curci in enumerate(self.bccols):
                matout[r,c]=U[curri,curci]
        return matout, vectout

    def BaseVectBode(self, wvect):
        if isinstance(wvect,list):
            wvect=array(wvect)
        if alltrue(wvect==real(wvect)):
            svect=wvect*1.0j
        else:
            svect=wvect
        nw=max(shape(wvect))
        matout=zeros((self.maxsize+1,nw),'d')+0j
        matout[self.maxsize,:]=1
        for n,s in enumerate(svect):
            subU,vect=self.FindAugSubmat(s)
            curb=dot(inverse(subU),-vect)
            for curcol,curent in zip(self.bccols,curb.tolist()):
                matout[curcol,n]=curent
        return matout


def null(A, eps=1e-10):
#    e,w=scipy.linalg.eig(A)
#    e=abs(e)
#    ind=minindex(e)
#    return scipy.conj(w[:,ind])
    u, s, vh = scipy.linalg.svd(A)
    mask = (s/s[0] <= eps) 
    null_space = scipy.compress(mask, vh, axis=0) 
    return scipy.transpose(scipy.conj(null_space)) 

def minindex(listin):
    mymin=min(listin)
    for x, ent in enumerate(listin):
        if ent==mymin:
            return x

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

def matmax(matin):
#    import MLab
    maxout=MLab.max(matin)
    if scipy.shape(maxout):
        maxout=max(maxout)
    return maxout

def matmaxi(matin):
    rowinds=scipy.argmax(matin,0)
    maxperrow=MLab.max(matin)
    colind=scipy.argmax(maxperrow)
    return [rowinds[colind],colind]

def matmaxabs(matin):
    if len(scipy.shape(matin))==1:
        ind=scipy.argmax(abs(matin))
        return matin[ind]
    else:
        ind=matmaxi(abs(matin))
        return matin[ind]
