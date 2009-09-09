from textfiles import rwkreadfile, textlist
from textfiles.latexlist import ReplaceList
import os
import pdb
import copy
from scipy import shape

from  IPython.Debugger import Pdb
#mytrace=Pdb().set_trace

rwkpydir='/home/ryan/rwkpython'
ws=' '*6

def _CleanList(listin):
    listin.replace('%PI','pi')
    return listin
    
def GetHeader(filename='header.f'):
    return rwkreadfile(os.path.join(rwkpydir,filename))

def PythonFileFromMaximaFortran(fortranname,headername, pythonname,newretname='',ws='\t'):
    mylist=textlist([],pythonname)
    mylist.list=[]
    mylist.readfile(headername)
    flist,retname=FortranToTextList(fortranname,newretname,ws)
    listout=[ws+item for item in flist]
    mylist.list.extend(listout)
    mylist=_CleanList(mylist)
    mylist.tofile()
    return mylist

def UnWrapFortranLine(filename, whitespace='\t'):
    operators=['+','-','/','*']
    flist=rwkreadfile(filename)
    outlist=[line for line in flist if line]
    outlist=textlist(outlist)
    x=0
    while x < len(outlist.list):
        curline=outlist[x]
        if x>0 and curline[5]!=' ':
            prevline=outlist[x-1]
            curline=curline[6:]
            curline=curline.strip()
            if prevline[-1]=='*' and curline[0]=='*':#then the line split on a ** and this needs to be fixed
                prevline=prevline[0:-1]
                curline='*'+curline
            if (prevline[-1] not in operators) and (curline[0] not in operators):#then we need to move somethings around
                optloc=len(curline)
                for op in operators:
                    curloc=curline.find(op)
                    if op=='*' and curloc<(len(curline)-1):
                        leftpart=''
                        rightpart=curline
                        while rightpart[curloc+1]=='*':#we really found ** which cannot be split
                            leftpart+=rightpart[0:curloc+2]
                            rightpart=rightpart[curloc+2:]
                            curloc=rightpart.find(op)
                        if curloc>-1:
                            curloc+=len(leftpart)
                    if curloc>-1:
                        if curloc<optloc:
                            optloc=curloc
                movepart=curline[0:optloc]
                leavepart=curline[optloc:]
                prevline=prevline+movepart
                curline=leavepart
            if curline:
                curline=whitespace+curline
                prevline+=' \\'
            if prevline[-1]=='*' and curline[0]=='*':
                prevline=prevline[0:-1]
                curline='*'+curline
            outlist[x-1]=prevline
        else:
            curline=curline.strip()
        if curline:
            outlist[x]=curline
            x+=1
        else:
            outlist.pop(x)
    return outlist.list

def FortranToTextList(filename,newretval='',whitespace='\t'):
    operators=['+','-','/','*']
    flist=rwkreadfile(filename)
    outlist=[line for line in flist if line]
    outlist=textlist(outlist)
#    pdb.set_trace()
    eqind=outlist.findall('=')
#    eqind=filter(eqlinefilter,eqind)
    valeqinds=[]
    for ln in eqind:
        curline=outlist[ln]
        if curline[5]==' ':
            valeqinds.append(ln)
    retind=max(valeqinds)
    retline=outlist[retind]#assume that the value to return is the lefthand side of the last equation
    retval,rest=retline.split('=',1)
    retval=retval.strip()
    if newretval:
        retline=retline.replace(retval,newretval)
        outlist[retind]=retline
        retval=newretval
    x=0
    while x < len(outlist.list):
#        if x==1072:
#            pdb.set_trace()
#        print('x='+str(x))
        curline=outlist[x]
        if x>0 and curline[5]!=' ':
            prevline=outlist[x-1]
            curline=curline[6:]
            curline=curline.strip()
            if prevline[-1]=='*' and curline[0]=='*':#then the line split on a ** and this needs to be fixed
                prevline=prevline[0:-1]
                curline='*'+curline
            if (prevline[-1] not in operators) and (curline[0] not in operators):#then we need to move somethings around
                optloc=len(curline)
                for op in operators:
                    curloc=curline.find(op)
                    if op=='*' and curloc<(len(curline)-1):
                        leftpart=''
                        rightpart=curline
                        while rightpart[curloc+1]=='*':#we really found ** which cannot be split
                            leftpart+=rightpart[0:curloc+2]
                            rightpart=rightpart[curloc+2:]
                            curloc=rightpart.find(op)
                        if curloc>-1:
                            curloc+=len(leftpart)
                    if curloc>-1:
                        if curloc<optloc:
                            optloc=curloc
                movepart=curline[0:optloc]
                leavepart=curline[optloc:]
                prevline=prevline+movepart
                curline=leavepart
            if curline:
                curline=whitespace+curline
                prevline+=' \\'
            if prevline[-1]=='*' and curline[0]=='*':
                prevline=prevline[0:-1]
                curline='*'+curline
            outlist[x-1]=prevline
        else:
            curline=curline.strip()
        if curline:
            outlist[x]=curline
            x+=1
        else:
#            pdb.set_trace()
            outlist.pop(x)
    outlist.append('return '+retval)
    return outlist.list,retval

def GetOptimizeVariableNames(filename,varprefix='a_'):
    """This function searched through a FORTRAN 
    file that was created by Maxima and the 
    Maxima function fortran_optimize, and 
    creates a list of all of the intermediate
    variables created by the fortran_optimize
    function.  These variables will all be 
    on the left hand side of equations at the
    begining of the file and will all start
    with varprefix.
    
    The list of variables is returned."""
    flist=rwkreadfile(filename)
    outlist=[line for line in flist if line]
    outlist=textlist(outlist)
    valeqlines=filter(nonextlinefilter,outlist.list)
    valeqlines=[line.strip() for line in valeqlines]
    allvars=map(GetLHS,valeqlines)
    optvars=[curvar for curvar in allvars if curvar.find(varprefix)==0]
    return optvars

def MakeDeclarationLine(varlist, vartype='double precision'):
    """This file takes a list of variables (varlist)
    and returns one FORTRAN line (starting with 6 spaces)
    that is made up of vartype var1, var2, ....

    This line may be lond and may need to be wrap across
    multiple lines using the function WrapFortranLine."""
    strout=' '*6+vartype+' '
    strout+=', '.join(varlist)
    return strout

def WrapLines(listin):
    linesout = []
    for line in listin:
        if (not line) or (line[0]=='C'):
            linesout.append(line)
        else:
            curwrap = WrapFortranLine(line)
            linesout.extend(curwrap)
    return linesout

def WrapFortranLine(linein, operators=['+','-','**',',','*','/']):
    """This function wraps lines longer than
    72 characters by looking for the right 
    most occurance of an operator in curline[0:73].
    
    Even if len(linein)<=72, this function will 
    still return a list of lines - but in that
    case the list will have only one line in it."""
    linesout=[]
    remainder=linein
    prefix=' '*5+'@'+' '*3
    first=1
    while len(remainder)>72:
        firstpart=remainder[0:72]
        curind=-1
        for opt in operators:
            myind=firstpart.rfind(opt)
            if myind>curind:
                curind=myind
        if curind==-1:
            raise StandardError, "could not find an operator to break the line at:\n"+remainder
        curline=remainder[0:curind+1]
        linesout.append(curline)
        remainder=prefix+remainder[curind+1:]
    linesout.append(remainder)
    return linesout

def GetAllRealVariables(allparams, alldefs, allcomp=[],unknownparams=[]):
    """alldef is assumed to be a list of equation strings whose left
    hand sides are variables that will need to be declared as either
    double precision or double complex.  allparams is a list of
    variable names. allcomp is a list of complex variables.  This is
    used to filter out complex values from allparams+GetLHS(alldefs).
    
    The return value is a list of real variables."""
    defvars=map(GetLHS,alldefs)
    myparams=allparams+defvars
#    Pdb().set_trace()
    notincluded=[item for item in unknownparams if item not in myparams]
    myparams+=notincluded
    realout=[item for item in myparams if item not in allcomp]
    return realout

class FortranVector:
    def __init__(self,name,real,length=None):
        self.name=name
        self.real=real
        self.length=length

def GetLength(listin, namein):
    for item in listin:
        if item.name==namein:
            return item.length
    return None

def _F2pyVectorDefs(listin):
    listout=[]
    defline='Cf2py integer intent(hide),depend(%s) :: %s = len(%s)'
    for item in listin:
        curdef=defline%(item.name, item.length, item.name)
        listout.append(curdef)
    return listout

def _GetNames(varlist):
    mylist=[]
    for item in varlist:
        if hasattr(item, 'name'):
            mylist.append(item.name)
        else:
            mylist.append(item)
    return mylist

def _GetNameStr(listin):
    return ', '.join(_GetNames(listin))

def _CreateInputLine(varlist):
    mylist=[]
    for item in varlist:
        if hasattr(item, 'name'):
            mylist.append(item.name)
        else:
            mylist.append(item)
    return 'Cf2py intent(in) '+', '.join(mylist)


def _ReplaceTagwithList(listin, tag, replist):
    listout=copy.copy(listin)
    while listout.findall(tag):
        inds=listout.findall(tag)
        curind=inds[0]
        listout.pop(curind)
        listout[curind:curind]=replist
    return listout
        
def _GetLengths(listin):
    strout=''
    for item in listin:
        if strout:
            strout+=', '
        if hasattr(item,'length'):
            strout+=item.length
        else:
            strout+=item
    return strout

def _InputNames(listin):
    return ', '.join([item.name for item in listin])

def _DeclareVectors(listin):
    imagline=''
    realline=''
    for item in listin:
        if item.real:
            if realline:
                realline+=', '+item.name+'('+item.length+')'
            else:
                realline=ws+'double precision '+item.name+'('+item.length+')'
        else:
            if imagline:
                imagline+=', '+item.name+'('+item.length+')'
            else:
                imagline=ws+'double complex '+item.name+'('+item.length+')'
    outlist=[]
    if realline:
        outlist.append(realline)
    if imagline:
        outlist.append(imagline)
    return outlist



def CreateHeader(TMMsysmodel,poly=True):
    mylist=textlist()
    mylist.list=[]
    mylist.readfile('newheader.f')
    ucv=TMMsysmodel.GetUnkownParams()
    vectinputs=[FortranVector('svect',False)]
    realvects=[]
    complexvects=[]
    externalinputs=[]
    if ucv:
        vectinputs.append(FortranVector('ucv',True))
    comps=TMMsysmodel.GetCompensators()
    if comps:
        if not poly:
            externalinputs.extend(comps)
        else:
            for comp in comps:
                vectinputs.append(FortranVector(comp+'_num',True))
                vectinputs.append(FortranVector(comp+'_den',True))
##    defvects=copy.copy(vectinputs)+copy.copy(realvects)+copy.copy(complexvects)
    for x,item in enumerate(vectinputs):
        item.length='n'+str(x)
##    vectinputs.extend(externalinputs)
    vectdefs=_F2pyVectorDefs(vectinputs)
    vectdefs.append(_CreateInputLine(vectinputs))
    vectdefs.append('Cf2py intent(out) outvect')
##    vectinputs.append('outvect')
    vints=_GetLengths(['i']+vectinputs)
    vectdefs.append(ws+'integer '+vints)
    outlength=GetLength(vectinputs,'svect')
    vectallvectors=vectinputs+[FortranVector('outvect',False,outlength)]
    vectdefs.extend(_DeclareVectors(vectallvectors))
    vectinstr=_InputNames(vectinputs)
    if externalinputs:
        vectinstr+=', '+', '.join(externalinputs)
    vectinstr+=', outvect, '+_GetLengths(vectinputs)
    mylist.replace('%%vectinputs%%',vectinstr)
    scalarinputs=vectinputs[1:]
    scalestr='s'
    if scalarinputs:
        scalardefs=_F2pyVectorDefs(scalarinputs)
        scalestr+=', '+_InputNames(scalarinputs)
        scalarintstr=_GetLengths(scalarinputs)
        scalestr+=', '+_GetLengths(scalarinputs)
    else:
        scalardefs=[]
    scalardefs.append(_CreateInputLine(['s']+scalarinputs))
    if scalarinputs:
        scalardefs.append(ws+'integer '+scalarintstr)
        scalardefs.extend(_DeclareVectors(scalarinputs))
    scalardefs.append(ws+'double complex s')
    mylist.replace('%%scalarinputs%%',scalestr)
    mylist.replace('%%scalarcallingscalar%%',scalestr)
    vcs=copy.copy(scalarinputs)
    vcs.insert(0,'svect(i)')
    vcsstr=_GetNameStr(vcs)
    if scalarinputs:
        vcsstr+=', '+_GetLengths(scalarinputs)
    mylist.replace('%%vectcallingscalar%%',vcsstr)
    mylist=_ReplaceTagwithList(mylist, '%%vectdefs%%', vectdefs)
    mylist=_ReplaceTagwithList(mylist, '%%scalardefs%%', scalardefs)
    templist = WrapLines(mylist.list)
    mylist.list = templist
    return mylist, scalarinputs

def _CompensatorLines(complist, inputlist):
    linesout=[]
    deflist=[]
    if not complist:
        return []
    for comp in complist:
        numstr=comp+'_num'
        numlen=GetLength(inputlist,numstr)
        numline=ws+comp+'N=poly(s,'+numstr+','+numlen+')'
        denstr=comp+'_den'
        denlen=GetLength(inputlist,denstr)
        denline=ws+comp+'D=poly(s,'+denstr+','+denlen+')'
        linesout.append(numline)
        linesout.append(denline)
        linesout.append(ws+comp+'='+comp+'N/'+comp+'D')
        deflist.append(comp+'N')
        deflist.append(comp+'D')
    defline=ws+'double complex poly,'+','.join(deflist)
    deflines=WrapFortranLine(defline)
    return deflines+linesout
                    
    
##       subroutine bodevect(svect,ucv,outvect,n1,n2)
## Cf2py integer intent(hide),depend(svect) :: n1 = len(svect)
## Cf2py integer intent(hide),depend(ucv) :: n2 = len(ucv)
## Cf2py intent(in) svect, ucv
## Cf2py intent(out) outvect
##       integer n1, n2, i
##       double complex svect(n1),outvect(n1), bode
##       double precision ucv(n2)

def MakeFortranFunction(filename, TMMsysmodel, curvefit=1, \
                        headername=None, newoutname='bode', aug=1, \
                        **unsedargs):
    """This function takes a filename that refers to a raw FORTRAN
    file output by Maxima and makes it into a file ready for f2py and
    Bode analysis."""
    
    print('curvefit='+str(curvefit))
    print('headername='+str(headername))
#    Pdb().set_trace()
#    mytrace()
    alldefs, allparams, defaultvalues, unknownparams, allsubs=TMMsysmodel.GetParamsandDefs(aug=aug)
    compvars=TMMsysmodel.GetAllComplexVariables()
    mycomp=copy.copy(compvars)
    if headername is None:
        linesout, scalarinputs = CreateHeader(TMMsysmodel)
        ## if ucv:
##             linesout=GetHeader()
##         else:
##             linesout=GetHeader('header_noucv.f')
    else:
        linesout=rwkreadfile(headername)
        scalarinputs=[]
##     if curvefit and ucv:
##         linesout.append(ws+'double precision ucv(n)')
    dpvars=GetAllRealVariables(allparams, alldefs, mycomp,unknownparams)
    dpline=MakeDeclarationLine(dpvars)
    dplines=WrapFortranLine(dpline)
    linesout.extend(dplines)
    optvars=GetOptimizeVariableNames(filename)
    mycomp.extend(optvars)
    filtdict={}
    for key,value in defaultvalues.iteritems():
        if shape(value):
            value=value[0]
        filtdict[key]=value
    if mycomp:
        optline=MakeDeclarationLine(mycomp,'double complex')
        optlines=WrapFortranLine(optline)
        linesout.extend(optlines)
    deflines=TMMsysmodel.GetAllFortranDefs()
#    Pdb().set_trace()
    defout=[ws+line.strip() for line in deflines]
    linesout.extend(defout)
    compensatorlines=_CompensatorLines(TMMsysmodel.GetCompensators(), scalarinputs)
    linesout.extend(compensatorlines)
    if curvefit and unknownparams:
        uplines=MakeUnknownParamLines(unknownparams)
        linesout.extend(uplines)
        myparams=[item for item in allparams if item not in unknownparams]
    else:
        myparams=allparams
    fplist=TMMsysmodel.GetAllFunctionParams()
    complist = TMMsysmodel.GetCompensators()
    excludelist = fplist+complist
    filtparams=[item for item in myparams if (item in filtdict.keys()) and (item not in excludelist)]
    paramlines=[ws+item +' = '+ str(filtdict[item]) for item in filtparams]
    linesout.extend(paramlines)    
    mydefs=map(MyFortranFilter,alldefs)
    mydefs=[ws+item for item in mydefs]
    linesout.extend(mydefs)
    newlines=ReplaceReturnName(filename,newoutname)
    linesout.extend(newlines)
    linesout.append(ws+'RETURN')
    linesout.append(ws+'END')
    fno,ext=os.path.splitext(filename)
    outname=fno+'_out'+ext
    mytextlist=textlist(linesout,outname)
    mytextlist.tofile()
    return outname

def MakeUnknownParamLines(unknownparams,python=0):
    """Create the lines used in the unknown parameter
    section of a FORTRAN or python bode function.
    If python=1, the lines start with a tab, square
    brakets are used instead of parenthesis, and
    the array index starts from 0."""
    linesout=[]
    if python:
        startind=0
        lp='['
        rp=']'
        myws='\t'
    else:
        startind=1
        lp='('
        rp=')'
        myws=' '*6
    for x,cp in enumerate(unknownparams):
        curline=myws+cp+'=ucv'+lp+str(x+startind)+rp
        linesout.append(curline)
    return linesout

def MyFortranFilter(linein):
    frlist=[('\\^','**'),('cos','zcos'),('sin','zsin'),('cosh','zcosh'),('sinh','zsinh')]
    strout=ReplaceList(linein,frlist)
    strout=strout.replace('^','**')#there may be problems with ')^('
    return strout


def ReplaceReturnName(filename,newname, oldname='RESULT',makecopy=False):
    """Replace the return value of a FORTRAN
    program (presumably one created by Maxima),
    by searching for a line that starts with
    6 spaces and then 'RESULT =' (where the 
    space is optional) and replacing RESULT on
    that line with newname.  

    Returns a textlist.  filename can also be a
    textlist or a list of text lines.
    
    If makecopy=True, a copy of the list is made
    before it is changed.  Otherwise, both the
    the list in and the list out will reflect
    the substitution of newname for oldname."""
    if type(filename)==str:
        mylist=rwkreadfile(filename)
    elif isinstance(filename, textlist):
        mylist=textlist.list
    elif type(filename)==list:
        mylist=filename
    else:
        raise StandardError, "filename must be a string, a textlist, or a list."
    if makecopy:
        mylist=copy.copy(mylist)
    myind=FindReturnLine(mylist,oldname)
    newline=mylist[myind]
    newline=newline.replace(oldname,newname)
    mylist[myind]=newline
    return mylist

def nonextlinefilter(line):
    return line[5]==' '

def GetLHS(line):
    """This function takes as an input a string
    that represents an equation and returns the
    left hand side."""
    strout,rest=line.split('=',1)
    return strout.strip()

def FindReturnLine(listin, retname='RESULT'):
    import re
    p=re.compile('      '+retname+' ?=')
    outind=-1
    for x,line in enumerate(listin):
        q=p.match(line)
        if q:
            outind=x
            break
    return outind

