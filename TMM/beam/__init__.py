from __future__ import division
from scipy import array, r_, c_, exp, pi, dot, real, imag, zeros, eye
#import MLab
import scipy
from scipy.linalg import det
import copy, re
import pdb

from  IPython.Debugger import Pdb

from TMM.TMMElement import TMMElement, HT4, Transform8by8, TMMElementLHT
##from TMM import HT4, Transform8by8

from rwkmisc import symstr, SymstrMattoMaxima

##import TMM
#reload(TMM)
##from TMM import TMMElement, HT4, Transform8by8
#TMMElement=TMM.TMMElement
#from TMM import HT4

def cos(ent):
    """Polymorphic cos function that adapts to either a numeric or
    symbolic string input."""
    if isinstance(ent,str):
        return symstr('cos('+ent+')')
    else:
        return scipy.cos(ent)

def cosh(ent):
    """Polymorphic cosh function that adapts to either a numeric or
    symbolic string input."""
    if isinstance(ent,str):
        return symstr('cosh('+ent+')')
    else:
        return scipy.cosh(ent)
   
def sin(ent):
    """Polymorphic sin function that adapts to either a numeric or
    symbolic string input."""
    if isinstance(ent,str):
        return symstr('sin('+ent+')')
    else:
        return scipy.sin(ent)

def sinh(ent):
    """Polymorphic sinh function that adapts to either a numeric or
    symbolic string input."""
    if isinstance(ent,str):
        return symstr('sinh('+ent+')')
    else:
        return scipy.sinh(ent)

class BeamElement(TMMElementLHT):
    """This class represents a continuous beam element.  It can be 1,
    2, or 3 dimensional.  In the 3D case, the x-axis includes axial
    and torsional vibration."""
    def __init__(self,params,symname='Ubeam',symlabel='beam', **kwargs):
        """Create a continous beam element.  params should be a
        dictionary with keys 'mu', 'L', and 'EI'.  If the beam is more
        than 1 dimensional (i.e. if maxsize!=4), than params should
        specify 'EI1' and 'EI2' as well as 'EA', 'm11', 'GJ' if axial
        and torsional analysis is to be performed.

        kwargs must include {'usez':True} if a 4x4 beam element should
        model bending about the z-axis.  The y-axis is used if 'usez'
        is False or not specified."""
        if kwargs.has_key('usez'):
            self.usez=kwargs['usez']
            kwargs.pop('usez')
        else:
            self.usez=False
        TMMElementLHT.__init__(self,'beam',params,symname=symname,symlabel=symlabel,**kwargs)
        #print('in BeamElement.__init__'+str(isinstance(self,TMM.TMMElementLHT)))

    def GetDefs(self):
        mydefs=['beta%s = (-1*s^2*L%s^4*mu%s/EI%s)^(0.25)']
        if not self.subparams:
            mydefs.append('a%s = L%s^2/EI%s')
            mydefs.append('c1%s = 0.5*cos(beta%s)+0.5*cosh(beta%s)')
            mydefs.append('c2%s = -sin(beta%s)+sinh(beta%s)')
            if self.usez:
                mydefs.append('c3%s = cos(beta%s)-cosh(beta%s)')
            else:
                mydefs.append('c3%s = -cos(beta%s)+cosh(beta%s)')
            mydefs.append('c4%s = sin(beta%s)+sinh(beta%s)')
        outdefs = [line.replace('%s',self.symlabel) for line in mydefs]
        return outdefs


    def GetMat(self,s,sym=False,debug=0):
        """Return the element transfer matrix for the BeamElement
        element.  If sym=True, 's' must be a symbolic string and a
        matrix of strings will be returned.  Otherwise, 's' is a
        numeric value (probably complex) and the matrix returned will
        be complex."""
#        Pdb().set_trace()
        if sym:
            myparams=self.symparams
        else:
            myparams=self.params
        if self.maxsize==4 and self.usez:
            bendmat1=bendmatz(s,myparams,'EI',symlabel=self.symlabel, subparams=self.subparams, debug=debug)
        elif myparams.has_key('EI1'):
            bendmat1=bendmaty(s,myparams,'EI1')
        elif myparams.has_key('EI'):
            bendmat1=bendmaty(s,myparams,'EI')
        else:
            raise KeyError, 'Neither EI nor EI1 found in params'
        if self.maxsize==4:
            return bendmat1
        elif self.maxsize>4:
            bendmat2=bendmatz(s,myparams,'EI2',subparams=self.subparams)
            zmat=scipy.zeros(scipy.shape(bendmat2))
         
        if self.maxsize==8:
            bigmat1=c_[bendmat1,zmat]
            bigmat2=c_[zmat,bendmat2]
            temp=r_[bigmat1,bigmat2]
            return Transform8by8(temp)
        elif self.maxsize==12:
            atmat=axialtorsionmat(s,myparams)
            t1=c_[bendmat1,zmat]
            t2=c_[zmat,bendmat2]
            temp=r_[t1,t2]
            temp=Transform8by8(temp)
            row1=c_[atmat,zmat,zmat]
            part2=c_[scipy.zeros((8,4)),temp]
#                row3=c_[zmat,zmat,bendmat2]
            return r_[row1, part2]

#    def GetHT(self):
#        return HT4(x=self.params['L'])

    def GetComplexList(self,label=None): 
        """This function returns a list of all variables associated
        with a beam bending about one axis that would need to be
        declared complex in a FORTRAN funciton that includes a beam in
        a Bode model.  This funcion is used in symbolic work."""
        strlist = ['beta']
        if not self.subparams:
            strlist.extend(['c1','c2','c3','c4'])
        if label is None:
            label=self.symlabel
##             if self.symlabel:
##                 label=self.symlabel
##             else:
##                 label='bz'
##                 self.symlabel=label
        listout=[item+label for item in strlist]
        return listout

    def GetMaximaString(self,name=None,label=None,N=None,aug=0):#omit
        if label is None:
            label=self.symlabel
##             if self.symlabel:
##                 label=self.symlabel
##             else:
##                 label='bz'
##                 self.symlabel=label
        if name is None:
            name=self.symname
##             if self.symname:
##                 name=self.symname
##             else:
##                 name='Ubeamz'
        maximalines=[]
        if N is None:
            N=self.maxsize
        Nstr=str(N)
        bzlines,defs,params=bendzsym(name,label)
        return bzlines, name, defs, params, label

    def GetMaximaSubstitutions(self,name=None,label=None,aug=0):
        """This function is used to substitute sinh and cosh terms
        into the symbolic matrix so that round off problems can be
        avoided by symbolic manipulation."""
        lines, defs, params, subs=self.GetMaximaLatexString(aug=aug)
#        Pdb().set_trace()
        filtdefs=[curdef for curdef in defs if curdef.find('beta')!=0]
        return filtdefs


    def GetMaximaLatexString(self,name=None,label=None,wrap=0,N=None,intro=None,aug=0):#omit
        if intro is None:
            intro='The transfer matrix for a beam element is given by'
        return TMMElementLHT.GetMaximaLatexString(self,name=name,label=label,wrap=wrap,N=N,intro=intro,aug=aug)
        
    def __GetMaximaLatexString__(self,name=None,label=None,wrap=0,N=None,intro=None,aug=0):#omit
        maximalines=[]
        if N is None:
            N=self.maxsize
        if aug:
            bzlines,myname, defs, params=self.GetAugMaximaString(name=name,label=label)
        else:
            bzlines,myname, defs, params=self.GetMaximaString(name=name,label=label)
        latexlines=[]
        if intro is None:
            latexlines.append('The transfer matrix for a beam elment is given by')
        else:
            latexlines.append(intro)
        if aug:
            newcolind=-1
            for x,line in enumerate(bzlines):
                if line.find('newcol:')==0:
                    newcolind=x
                    break
            wherepart=bzlines[1:newcolind]
            bzlines[1:newcolind]=[]
            preind=len(bzlines)
            bzlines.extend(wherepart)
        else:
            preind=1
        prelines=bzlines[0:preind]
        for line in prelines:
            latexlines.append('\\begin{maxima-noout}')
            latexlines.append(line)
            latexlines.append('\\end{maxima-noout}')
        latexlines.append('\\begin{maxima}')
        latexlines.append("\t\\parseopts{lhs='"+myname+"',wrap="+str(wrap)+"}")
        latexlines.append('\t'+myname)
        latexlines.append('\\end{maxima}')
        latexlines.append('where')
        deflines=bzlines[preind:]
        for dl in deflines:
            latexlines.append('\\begin{maxima}')
            latexlines.append('\t'+dl)
            latexlines.append('\\end{maxima}')
        return latexlines, defs, params

def bendzsym(name='Ubeamz',label='bz'):#omit
    linesout=[]
    params=[]
    defs=[]
    matstr=name+':matrix([c1,1/2*L*c4/beta,-1/2*a*c3/beta^2,-1/2*L*a*c2/beta^3],[1/2*beta*c2/L,c1,1/2*a*c4/beta/L,1/2*a*c3/beta^2],[-1/2*beta^2*c3/a,1/2*beta*L*c2/a,c1,-1/2*L*c4/beta],[-1/2*beta^3*c4/L/a,1/2*beta^2*c3/a,-1/2*beta*c2/L,c1])'
    linesout.append(matstr)
    astr='a'+label
    bstr='beta'+label
    mustr='mu'+label
    EIstr='EI'+label
    Lstr='L'+label
#    linesout.append('eq'+astr+':'+astr+'='+Lstr+'^2/'+EIstr)
#    linesout.append('eq'+bstr+':'+bstr+'=(-1*s^2*'+Lstr+'^4*'+mustr+'/'+EIstr+')^(0.25)')
    defs.append(astr+'='+Lstr+'^2/'+EIstr)
    defs.append(bstr+'=(-1*s^2*'+Lstr+'^4*'+mustr+'/'+EIstr+')^(0.25)')
#    defs.append(bstr+'=(-1*s^2)^0.25*('+Lstr+'^4*'+mustr+'/'+EIstr+')^(0.25)')
    clist=['0.5*cos(beta)+0.5*cosh(beta)','-sin(beta)+sinh(beta)','cos(beta)-cosh(beta)','sin(beta)+sinh(beta)']
    params=[mustr,EIstr,Lstr]
    for x,curc in enumerate(clist):
        cstr='c'+str(x+1)
#        cline='eq'+cstr+':'+cstr+'='+curc
        cline=cstr+'='+curc
#        junk,rp=cline.split(':',1)
#        defs.append(rp)
#        linesout.append(cline)
        defs.append(cline)
    replines=[]
    strlist=['c1','c2','c3','c4','beta','EI','a','L']
    for line in linesout:
#        pdb.set_trace()
        for ent in strlist:
#            line=line.replace(ent,ent+label)
            line=re.sub('\\b'+ent+'\\b',ent+label,line)
        replines.append(line)
    tempdef=[]
    for line in defs:
        for ent in strlist:
            line=re.sub('\\b'+ent+'\\b',ent+label,line)
        tempdef.append(line)
    return replines, tempdef, params#[0],replines[1:]


def bendmaty(s,params,EIstr='EI1',symlabel=''):
    """Return a 4x4 transfer matrix for a beam in bending about the
    y-axis.  This function is used for part of the GetMat function of
    an 8x8 or 12x12 beam.  It will be the entire 4x4 transfer matrix
    if usez=False for the beam element."""
    EI=params[EIstr]
    mu=params['mu']
    L=params['L']
#[                 c(1),     -1/2*L*c(4)/beta,   -1/2*a*c(3)/beta^2, -1/2*L*a*c(2)/beta^3]
#[     -1/2*beta*c(2)/L,                 c(1),    1/2*a*c(4)/beta/L,    1/2*a*c(3)/beta^2]
#[   -1/2*beta^2*c(3)/a,    1/2*beta*L*c(2)/a,                 c(1),      1/2*L*c(4)/beta]
#[ -1/2*beta^3*c(4)/L/a,    1/2*beta^2*c(3)/a,      1/2*beta*c(2)/L,                 c(1)]
#
#   This matrix is derived in bending_about_y_axis.m originally in 
#   E:\GT\Research\SAmii\modeling\transfer_matrix\threeD_beam_derivations
#
#   The states that multipy the matrix are [w_z;theta_y;M_y;V_z] and I have gone away from using -w as the state i.e. +w_z is the first state.
    if isinstance(s,str):
        beta=symstr('beta'+symlabel)
        a=symstr('a'+symlabel)
        c1=symstr('c1'+symlabel)
        c2=symstr('c2'+symlabel)
        c3=symstr('c3'+symlabel)
        c4=symstr('c4'+symlabel)
    else:
        beta=pow((-1*s*s*L**4*mu/EI),0.25)
        a=L*L/EI
        c1=0.5*cos(beta)+0.5*cosh(beta)
        c2=-sin(beta)+sinh(beta)
        c3=-cos(beta)+cosh(beta)
        c4=sin(beta)+sinh(beta)
    outmat=array([[c1,-0.5*L*c4/beta,-0.5*a*c3/(beta*beta),-0.5*L*a*c2/beta**3],[-0.5*beta*c2/L,c1,0.5*a*c4/beta/L,0.5*a*c3/(beta*beta)],[-0.5*(beta*beta)*c3/a,0.5*beta*L*c2/a,c1,0.5*L*c4/beta],[-0.5*beta**3*c4/L/a,0.5*(beta*beta)*c3/a,0.5*beta*c2/L,c1]])
    return outmat

def bendmatz(s, params, EIstr='EI2', symlabel='', subparams=False, debug=0):
    """Return a 4x4 transfer matrix for a beam in bending about the
    z-axis.  This function is used for part of the GetMat function of
    an 8x8 or 12x12 beam.  It will be the entire 4x4 transfer matrix
    if usez=True for the beam element."""
    if debug>0:
        print('In bendmatz')
    EI=params[EIstr]
    mu=params['mu']
    L=params['L']
#    [                 mys(1),      1/2*L*mys(4)/beta,   -1/2*a*mys(3)/beta^2, -1/2*L*a*mys(2)/beta^3]
#    [      1/2*beta*mys(2)/L,                 mys(1),    1/2*a*mys(4)/beta/L,    1/2*a*mys(3)/beta^2]
#    [   -1/2*beta^2*mys(3)/a,    1/2*beta*L*mys(2)/a,                 mys(1),     -1/2*L*mys(4)/beta]
#    [ -1/2*beta^3*mys(4)/L/a,    1/2*beta^2*mys(3)/a,     -1/2*beta*mys(2)/L,                 mys(1)]
#
#   This matrix is derived in bending_about_z_axis.m originally in 
#   E:\GT\Research\SAmii\modeling\transfer_matrix\threeD_beam_derivations
#
#   The states that multipy the matrix are [w_y;theta_z;M_z;V_y] and I have gone away from using -w as the state i.e. +w_y is the first state.
    if isinstance(s,str):
        beta=symstr('beta'+symlabel)
        if subparams:
            a=L*L/EI
            c1=0.5*cos(beta)+0.5*cosh(beta)
            c2=-sin(beta)+sinh(beta)
            c3=cos(beta)-cosh(beta)
            c4=sin(beta)+sinh(beta)
        else:       
            a=symstr('a'+symlabel)
            c1=symstr('c1'+symlabel)
            c2=symstr('c2'+symlabel)
            c3=symstr('c3'+symlabel)
            c4=symstr('c4'+symlabel)
    else:
        beta=pow((-1*s*s*L**4*mu/EI),0.25)
        # there is an error in my thesis and this is it: (eqn 285)
        #a = beta/(L**2)#to check Brian Posts work  
        a=L*L/EI
        if debug>0:
            print('thesis='+str(thesis))
            print('a='+str(a))
            print('beta='+str(beta))
        c1=0.5*cos(beta)+0.5*cosh(beta)
        c2=-sin(beta)+sinh(beta)
        c3=cos(beta)-cosh(beta)
        c4=sin(beta)+sinh(beta)
        clist = [c1,c2,c3,c4]
        if debug>0:
            for ind, c in enumerate(clist):
                print('c'+str(ind+1)+'='+str(c))
#    outmat=array([[c1,0.5*L*c4/beta,-0.5*a*c3/(beta*beta),-0.5*L*a*c2/beta**3],[0.5*beta*c2/L,c1,0.5*a*c4/beta/L,0.5*a*c3/(beta*beta)],[-0.5*(beta*beta)*c3/a,0.5*beta*L*c2/a,c1,-0.5*L*c4/beta],[-0.5*beta**3*c4/L/a,0.5*(beta*beta)*c3/a,-0.5*beta*c2/L,c1]])
#    outmat=zeros((4,4),'D')
#    outmat[0,:]=[c1,0.5*L*c4/beta,-0.5*a*c3/(beta*beta),-0.5*L*a*c2/beta**3]
#    outmat[1,:]=[0.5*beta*c2/L,c1,0.5*a*c4/beta/L,0.5*a*c3/(beta*beta)]
#    outmat[2,:]=[-0.5*(beta*beta)*c3/a,0.5*beta*L*c2/a,c1,-0.5*L*c4/beta]
#    outmat[3,:]=[-0.5*beta**3*c4/L/a,0.5*(beta*beta)*c3/a,-0.5*beta*c2/L,c1]
    row1=[c1,0.5*L*c4/beta,-0.5*a*c3/(beta*beta),-0.5*L*a*c2/beta**3]
    row2=[0.5*beta*c2/L,c1,0.5*a*c4/beta/L,0.5*a*c3/(beta*beta)]
    row3=[-0.5*(beta*beta)*c3/a,0.5*beta*L*c2/a,c1,-0.5*L*c4/beta]
    row4=[-0.5*beta**3*c4/L/a,0.5*(beta*beta)*c3/a,-0.5*beta*c2/L,c1]
    outmat=array([row1,row2,row3,row4])
    return outmat

def axialtorsionmat(s,params):#L,mu,EA,m11,GJ):
    """Return a 4x4 transfer matrix for axial and torsional vibration
    of a beam about its x-axis.  This function is used for part of the
    GetMat function of a 12x12 beam."""
    axmat=axialmat(s,params)
    tormat=torsionmat(s,params)
#    matout=scipy.zeros((2*scipy.shape(axmat)[0],2*scipy.shape(axmat)[1]))
    matout=scipy.zeros((4,4))
    matout=matout*1j
    matout[1:3,1:3]=tormat#python slice assignment goes to 1 less than the last index
    matout[0,0]=axmat[0,0]
    matout[0,3]=axmat[0,1]
    matout[3,0]=axmat[1,0]
    matout[3,3]=axmat[1,1]
    return matout
            
def torsionmat(s,params):#L,m11,GJ):
    """Return the torsional vibration 2 by 2 transfer matrix used as
    part of axialtorsionmat."""
    #the states for this matrix are [[theta],[tau]]
#    [                     1/2*exp(sigma*L)+1/2*exp(-sigma*L), 1/2*exp(sigma*L)/G/J/sigma-1/2*exp(-sigma*L)/G/J/sigma]
#    [ 1/2*G*J*sigma*exp(sigma*L)-1/2*G*J*sigma*exp(-sigma*L),                     1/2*exp(sigma*L)+1/2*exp(-sigma*L)]
    L=params['L']
    m11=params['m11']
    GJ=params['GJ']
    sigma=s*scipy.sqrt(m11/GJ)
    matout=scipy.array([[0.5*exp(sigma*L)+0.5*exp(-sigma*L), 0.5*exp(sigma*L)/GJ/sigma-0.5*exp(-sigma*L)/GJ/sigma],[0.5*GJ*sigma*exp(sigma*L)-0.5*GJ*sigma*exp(-sigma*L),0.5*exp(sigma*L)+0.5*exp(-sigma*L)]])
    return matout


def axialmat(s,params):#L,mu,EA):
    """Return the axial vibration 2 by 2 transfer matrix used as
    part of axialtorsionmat."""
    #axial
    #[1/2*exp(sigma*L)+1/2*exp(-sigma*L),1/2*exp(sigma*L)/E/A/sigma-1/2*exp(-sigma*L)/E/A/sigma]
    #[1/2*E*A*sigma*exp(sigma*L)-1/2*E*A*sigma*exp(-sigma*L),1/2*exp(sigma*L)+1/2*exp(-sigma*L)]
    L=params['L']
    mu=params['mu']
    EA=params['EA']
    sigma=s*scipy.sqrt(mu/EA)
    #
    #   This matrix is derived in axial_rod.m which was originally in
    #   E:\GT\Research\Samii\modeling\transfer_matrix\torsional_shaft_derivation
    #
    #   This matrix assumes +w_x and V_x are the states
    #
    matout=array([[0.5*exp(sigma*L)+0.5*exp(-sigma*L),0.5*exp(sigma*L)/(EA*sigma)-0.5*exp(-sigma*L)/(EA*sigma)],[0.5*EA*sigma*exp(sigma*L)-0.5*EA*sigma*exp(-sigma*L),0.5*exp(sigma*L)+0.5*exp(-sigma*L)]])
    return matout


def bendmatz_comp(s, params, EIstr='EI2', symlabel='', \
                  subparams=False, debug=0):
    """Return a 4x4 transfer matrix for a beam in bending about the
    z-axis.  This function is used for part of the GetMat function of
    an 8x8 or 12x12 beam.  It will be the entire 4x4 transfer matrix
    if usez=True for the beam element."""
    if debug>0:
        print('In bendmatz')
    EI = params[EIstr]
    mu = params['mu']
    L = params['L']
    if not params.has_key('c'):
        c = 0.0#no damping
    elif params['c'] == 0.0:
        c = 0.0
    else:
        w = imag(s)
        c = params['c']/w

##     if abs(s) > 5*2*pi:
##         c = 0.0
    beta=pow((-1*s*s*L**4*mu/(EI*(c*s+1))),0.25)
    
    d1 = 0.5*(cos(beta)+cosh(beta))
    d2 = 0.5*(sinh(beta)-sin(beta))
    d3 = 0.5*(cosh(beta)-cos(beta))
    d4 = 0.5*(sin(beta)+sinh(beta))

    # there is an error in my thesis and this is it: (eqn 285)
    #a = beta/(L**2)#to check Brian Posts work  
    a=L*L/EI
    
    outmat = array([[d1, L*d4/beta, a*d3/(beta**2*(1 + c*s)), \
                     -L*a*d2/(beta**3*(1 + c*s))], \
                    [beta*d2/L, d1, a*d4/(L*beta*(1 + c*s)), \
                     -a*d3/(beta**2*(1 + c*s))], \
                    [d3*beta**2*(1 + c*s)/a, L*beta*d2*(1 + c*s)/a, \
                     d1, -L*d4/beta], \
                    [-d4*beta**3*(1 + c*s)/(L*a), \
                     -d3*beta**2*(1 + c*s)/a, -beta*d2/L, d1]])
    return outmat



class BeamElement_v2(BeamElement):
    """This class is being used to test my newly derived beam matrix
    with damping (complex stiffness).  For now, it only overrides the
    GetMat method of BeamElement."""

    def GetMat(self,s,sym=False,debug=0):
        """Return the element transfer matrix for the BeamElement
        element.  If sym=True, 's' must be a symbolic string and a
        matrix of strings will be returned.  Otherwise, 's' is a
        numeric value (probably complex) and the matrix returned will
        be complex."""
        #Pdb().set_trace()
        if sym:
            myparams=self.symparams
        else:
            myparams=self.params
        if self.maxsize==4 and self.usez:
            bendmat1=bendmatz_comp(s, myparams, 'EI',\
                                   symlabel=self.symlabel, \
                                   subparams=self.subparams, \
                                   debug=debug)
        elif myparams.has_key('EI1'):
            bendmat1=bendmaty(s,myparams,'EI1')
        elif myparams.has_key('EI'):
            bendmat1=bendmaty(s,myparams,'EI')
        else:
            raise KeyError, 'Neither EI nor EI1 found in params'
        if self.maxsize==4:
            return bendmat1
        elif self.maxsize>4:
            bendmat2=bendmatz_comp(s, myparams, 'EI2', \
                                   subparams=self.subparams)
            zmat=scipy.zeros(scipy.shape(bendmat2))
         
        if self.maxsize==8:
            bigmat1=c_[bendmat1,zmat]
            bigmat2=c_[zmat,bendmat2]
            temp=r_[bigmat1,bigmat2]
            return Transform8by8(temp)
        elif self.maxsize==12:
            atmat=axialtorsionmat(s,myparams)
            t1=c_[bendmat1,zmat]
            t2=c_[zmat,bendmat2]
            temp=r_[t1,t2]
            temp=Transform8by8(temp)
            row1=c_[atmat,zmat,zmat]
            part2=c_[scipy.zeros((8,4)),temp]
#                row3=c_[zmat,zmat,bendmat2]
            return r_[row1, part2]
