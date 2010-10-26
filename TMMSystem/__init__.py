from __future__ import division
from scipy import cos, cosh, sin, sinh, array, r_, c_, exp, pi, dot, real, imag, zeros, eye, alltrue, shape, atleast_1d, dot, vstack, isscalar
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
from rwkmisc import symstr, SymstrMattoMaxima, rwkstr, reverse

from TMM.TMMElement import TMMElement

import inspect
import re
from rwkparse import GetPythonFunctionArgs

class TMMSubSystem(list):
    def invert(self):
        print("This method should return an inverted transfer matrix for non-collocated feedback.")

        
class TMMSystem:
    def __init__(self,elemlist,bcend=[2,3,6,7],bcbase=[0,1,4,5],bodeouts=[]):
        self.sysHT=[]
        myelemlist=[]
        self.classlist=[item.__class__ for item in elemlist]
        self.inputlist=elemlist
        for elem in elemlist:
            if isinstance(elem,TMMSubSystem):
                myelemlist.extend(elem)
            else:
                myelemlist.append(elem)
        self.elemlist=myelemlist#a list of TMMElement instances
        self.maxsize=self.elemlist[0].maxsize
        self.bcend=self._parsebcstr(bcend)#a list of the boundary condition indices that are 0 at the end of the chain (the free end in the case of a cantilever configuration) - the defaults are for a 1D cantilever where the states are [-w, phi, M, V] and M and V are 0 at the free end
        self.bcbase=self._parsebcstr(bcbase)#a list of the boundary configuration indices that are 0 as the base of the chain (the clamped beginning of a cantilever configuration) - the defaults are for a 1D cantilever where -w and phi are 0 at the base
        templist=range(self.maxsize)
        self.bccols=[]
        for ent in templist:
            if ent not in self.bcbase:
                self.bccols.append(ent)#in evaluating the eigvalues, the columns of the system U matrix that multiply the non-zero elements of the base state vector are the ones whose sub-determinant must go to 0
        self.bodeouts=[]
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
        self.rawbodeouts=[]
        for curout in self.bodeouts:
            scalar=0
            if not shape(curout.ind):
                scalar=1
            myouts=[]
            curents=atleast_1d(curout.ind)
            for curent in curents:
                if isinstance(curent,TMMElement):
                    ind=self.elemlist.index(curent)
                elif type(curent)==types.IntType:
                    ind=curent
                else:
                    raise ValueError, "Unknown type for Bode output:"+str(curent)     
                myouts.append(ind)
            for tempout in myouts:
                tempraw=copy.deepcopy(curout)
                tempraw.pind=tempout
                self.rawbodeouts.append(tempraw)
            if scalar:
                curout.pind=myouts[0]
            else:
                curout.pind=myouts

    def PreparePythonFiles(self, filenames,chardetnames=[],ratmx=False,curvefit=True):
        """Similar to PrepareFortranFiles, but this function generates
        python files instead of FORTRAN files.  The input is a list of
        auto-generated FORTRAN files (raw FORTRAN files output by
        Maxima).  If these were processed with ratmx=True, then a list
        of chardetnames may also be passed.  These functions would
        contain expressions for the characteristic determinant of the
        transfer function.  The output is a list of python filenames
        and a vector of the initial guesses, which come from the
        parameter values declared in the model for the unknownparams."""
        if not ratmx:
            chardetnames=[item+'_cd' for item in filenames]
        if curvefit:
            pynames,igv=self.SymBodesParseMaxima(filenames,chardetnames,curvefit=curvefit)
            return pynames,igv
        else:
            pynames=self.SymBodesParseMaxima(filenames,chardetnames,curvefit=curvefit)
            return pynames

    def CleanPythonFiles(self, filenames,frlistin=[]):
        frlist=[('SINH','sinh'),('SIN','sin'),('COSH','cosh'),('COS','cos'),('SQRT','sqrt')]
        frlist.extend(frlistin)
        for filename in filenames:
            mylist=textlist([],filename)
            mylist.list=[]
            mylist.readfile(filename)
            mylist.replacev(frlist)
            mylist.tofile()

    def CallF2py(self, finalfnames):
        """This function takes as input a list of FORTRAN files called
        finalfnames which are ready to compile FORTRAN files
        presumably generated by PrepareFortranFiles.  the '_out.f' is
        stripped off of each filename and an 'f' is appended to create
        the module name.  f2py is then called with

        os.system('f2py -c -m '+modulename+' '+filename)

        for each file in finalfnames.  The list of modulenames is
        returned."""
        modnames=[]
        for filename in finalfnames:
            fno,ext=os.path.splitext(filename)
            ind=filename.find('_out')
            if ind>-1:
                fno=fno[0:ind]
            modname=fno+'f'
            cmd='f2py -c -m '+modname+' '+filename
            os.system(cmd)
            modnames.append(modname)
        return modnames

    def GetAllFortranDefs(self):
        outlist=[]
        for elem in self.elemlist:
            curflist=elem.GetFortranDefs()
            outlist.extend(curflist)
        return outlist
                        
    def PrepareFortranFiles(self, filenames,headername=None,**kwargs):
        """This function takes a list of FORTRAN filenames which are
        files that were auto-generated by Maxima (raw FORTRAN files).
        It gets parameters and other definitions from the TMMSystem
        model and calls rwkfortran.MakeFortranFunction to create
        FORTRAN files that are ready to be compiled with f2py -
        including adding headers.  The output is a list of the
        ready-to-compile filenames."""        
#        alldefs, allparams, defaultvalues, unknownparams, allsubs=self.GetParamsandDefs()
#        compvars=self.GetAllComplexVariables()
        print('In PrepareFortranFiles')
        for key,value in kwargs.iteritems():
            print(key+'='+str(value))
        outfiles=[]
        for curfile in filenames:
#            outname=rwkfortran.MakeFortranFunction(curfile,allparams,alldefs,compvars,defaultvalues,unknownparams,headername=headername)
            outname=rwkfortran.MakeFortranFunction(curfile,self,headername=headername,**kwargs)
            outfiles.append(outname)
        return outfiles

    def GetAllComplexVariables(self):
        """This function calls GetComplexList for
        each element in self.elemlist in order to
        return a list of all the complex variables
        that will need to be declared in a FORTRAN
        function generated by Maxima."""
        allcomp=[]
        for elem in self.elemlist:
            curcomp=elem.GetComplexList()
            allcomp.extend(curcomp)
        return allcomp

    def GenMaxima(self, prefix='', skipMaxima=True, \
                  showbode=False, newsym=True, **kwargs):
        """This file generates a LaTeX file that describes the
        TMMSystem.  This file can be used as an input to the
        Python -> Maxima -> LaTeX symoblic engine.  Optionally, if
        skipMaxima=False, Maxima will be called, Python or FORTRAN
        files will be generated and an output LaTeX file will be
        created that contains the symoblic results from Maxima.  The
        return values are outname, bodenames, chardetnames."""
        if prefix:
            basename=prefix+'_bode'
            basecdname=prefix+'_chardet'
            kwargs['texname']=prefix+'_dev.tex'
        else:
            basename='bode'
            basecdname='chardet'
        kwargs['basebodename']=basename
        kwargs['basecdname']=basecdname
        kwargs['skipMaxima']=skipMaxima
        kwargs['showbode']=showbode
        if newsym:
            myfunc=self.SymBodeAll
        else:
            myfunc=self.SymBodeMaximaAll
        outname,bodenames,chardetnames=myfunc(**kwargs)
        return outname,bodenames,chardetnames

    def RunMaxima(self, texname,frlist=[]):
        """Takes as an input the name of a LaTeX file and runs the
        symoblic Python->Maxima->LaTeX engine on that file.  The
        output of that file is saved to raw FORTRAN files which will
        need to have headers added to them and some other parsing done
        before they can be compiled (FORTRAN_ or imported (Python).
        PreparePythonFiles and PrepareFortranFiles take care of this
        additional parsing.

        The resulting LaTeX outputs from Maxima are used to
        replace the maxima environments in the input LaTeX file and an
        output LaTeX file is generated.  The return value is the name
        of that output LaTeX file."""
        mylist=latexlist([])
        mylist.readfile(texname)
        mylist.filename=texname
        outname=mylist.ProcessMaxima(frlist)
        if frlist:
            mylist.LatexSubs(frlist)
            mylist.filename=outname
            mylist.tofile()
        return outname
    
    def GetAllFunctionParams(self):
        fps=[]
        for elem in self.elemlist:
            if hasattr(elem,'functionparams'):
                cfp=elem.functionparams
                if cfp:
                    fps.extend(cfp)
        return fps

    def IntegratedCurveFit(self,expmodname='',prefix='',fitfortran=True,skipMaxima=False,**kwargs):
        """This function runs an automated system identification process for TMMSystem objects in three steps: 
        1. Create Python and/or FORTRAN functions for the Bode responses of the system.  These functions will have the names prefix+'bode'+n, where the integer n is from iterating over the defined bodeouts of the TMMSystem. 
        2. Prepare the experimental data for easy use with the cost function of the curve-fitting routine.  This is done by taking the name of a module contaiing experimental rwkbode instances, searching for matches to the input/output pairs defined in the bodeouts of the model, and creating a list of exp. bodes from these matches and saving the list to a new module.  Saved and loaded modules come from the scipy.io.save method that shelves a dictionary. 
        3. Autogenerate a script that will actually run the the minimization that does the curve fitting. 
        
        prefix will be used in a number of other places for filenames and such. 
        
        kwargs will be passed to SymBodeMaximaAll and consist of the following:
            texname='SymBodeDev.tex'
            frlist=[]
            basebodename='bode'
            basecdname='chardet'
            debug=False
            grind=True
            optimize=True
            skipMaxima=False
        """
#        pdb.set_trace()
        if prefix:
            basename=prefix+'_bode'
            basecdname=prefix+'_chardet'
        else:
            basename='bode'
            basecdname='chardet'
        kwargs['basebodename']=basename
        kwargs['basecdname']=basecdname
        oldway=0
        if oldway:
            outname,bodenames,chardetnames=self.SymBodeMaximaAll(**kwargs)
    #        pdb.set_trace()
            parseoptslist=['debug','grind','optimize','ratmx']
            symparseopts=Getkwargs(kwargs,parseoptslist)
            symparseopts['curvefit']=True
            pyfiles,ig=self.SymBodesParseMaxima(bodenames,chardetnames,**symparseopts)
            if expmodname:
                savedname=self.PrepareExpData(expmodname,prefix)
        #        pdb.set_trace()
                autogenname=self.AutoGenerateFitFunction(savedname,pyfiles,ig,prefix)
                return autogenname
            else:
                return None
        else:
            outname,bodenames,chardetnames=self.GenMaxima(prefix)
            print('outname='+outname)
            if not skipMaxima:
                latexoutname=self.RunMaxima(outname)
                print('latexoutname='+latexoutname)
            myfinalfnames=self.PrepareFortranFiles(bodenames,**kwargs)
            mypynames=self.PreparePythonFiles(bodenames)
            savedname=self.PrepareExpData(expmodname,prefix)
            ig=self.ExtractInitialGuesses()
            if fitfortran:
                fmodnames=self.CallF2py(myfinalfnames)
                autoname=self.AutoGenerateFitFunction(savedname,fmodnames,functioname='bodevect',initialguesses=ig,prefix=prefix,fortran=True)
            else:
                autoname=self.AutoGenerateFitFunction(savedname,bodenames,initialguesses=ig,prefix=prefix)
            return autoname

    def CreateFortranandPythonModules(self,prefix,runcompile=True, newsym=True):
        """This function generates the Maxima code, runs Maxima,
        prepares the Fortran files, and optionally calls F2Py.  It
        essentially does the initial steps of IntegratedCurveFit, but
        does not step up the experimental data nor create the cost
        functions."""
        outname,bodenames,chardetnames=self.GenMaxima(prefix, newsym=newsym)
        latexoutname=self.RunMaxima(outname)
        myfinalfnames=self.PrepareFortranFiles(bodenames)
        mypynames=self.PreparePythonFiles(bodenames)
        if runcompile:
            fmodnames=self.CallF2py(myfinalfnames)

        
    def FindExpData(self,modulename,freqnames=['freq','expfreq'],bodenames=['bodes','bode','expbodes','expbode'],datadir=None,prefix='exp'):
        """Import the module named modulename and find the 
        Bodes in it that match the bodeouts of the TMMSystem 
        model.  modulename is assumed to be a module saved
        using the scipy.io.save method that contains a 
        frequency vector modulename.freqname.  A list of 
        Bodes is expected to be found at modulename.bodename
        where freqname and bodename are one of the values
        in freqnames and bodenames respectively.  They will
        be tried in order and a break will be executed once
        one is found using hasattr.

        The method returns a dictionary whose keys are 
        prefix+'freq' (the frequency vector) and 
        prefix+'bodes' (a list of Bodes in the same order 
        as self.bodeouts).
        
        If datadir is specified, it will be appended to
        sys.path if it is not already in it.
        
        Alternately, modulename can be a dictionary instead
        of a module name.  If so, it should have keys named
        freqname and bodename1 or bodename2."""
        if datadir:
            if datadir not in sys.path:
                sys.path.append(datadir)
        if type(modulename)==dict:
            expdata=rwkmisc.dictobject(**modulename)
        else:
            expdata=my_import(modulename)
        for freqname in freqnames:
            if hasattr(expdata,freqname):
                f=getattr(expdata,freqname)
                break
        else:
            raise IndexError, 'Cannot find the frequency variable in expdata imported from module '+modulename
        expbodes=[]
        for bode in self.bodeouts:
            for bodename in bodenames:
                if hasattr(expdata,bodename):
                    tempbodes=getattr(expdata,bodename)
                    if not (type(tempbodes)==list or type(tempbodes)==tuple):
                        tempbodes=[tempbodes]
                    curbode=rwkbode.FindMatch(tempbodes,bode.output, bode.input)
                    break
            else:
                raise IndexError, 'Cannot find '+bodename1+' or '+bodename2+' in expdata imported from module '+modulename
            curbode=rwkbode.renew(curbode)
            expbodes.append(curbode)
        mydict={prefix+'freq':f,prefix+'bodes':expbodes}
        return mydict

    def PrepareExpData(self,modulename,prefix='',**otherargs):
        if prefix:
            savename=prefix+'_expdata'
        else:
            savename='expdata'
        mydict=self.FindExpData(modulename,**otherargs)
        save(savename,mydict)
        return savename

    def AutoGenerateFitFunction(self,savedname,bodenames,functioname=None,initialguesses=[],prefix='',fortran=False):
        ws=' '*4
        if prefix:
            outname=prefix+'_autofitfunc'
        else:
            outname='autofitfunc'
        if fortran:
            outname+='f'
        freqname='expfreq'
        bodename='expbodes'
#        bodenames=['bode0','bode1']
        header='header.py'
        #def GenFitFun(outname,savename,freqname,bodename):
        fno,ext=os.path.splitext(outname)
        outname=fno+'.py'
        pylist=textlist([],outname)
        pylist.prependfile(header)
        out=pylist.append
        out('from scipy.io import write_array, read_array')
        out('mt=time.time')
        out('import '+savedname)
        funcstr='funclist=['
        for n,cf in enumerate(bodenames):
            cfno,ext=os.path.splitext(cf)
            out('import '+cfno)
            if not functioname:
                out('func'+str(n)+'='+cfno+'.'+cfno)
            else:
                out('func'+str(n)+'='+cfno+'.'+functioname)
            funcstr+='func'+str(n)+','
        funcstr=funcstr[0:-1]
        funcstr=funcstr+']'
        out(funcstr)
        out('global fitf,s,mybodes')
        out('fitf='+savedname+'.'+freqname)
        out('fitf=array(fitf)')
        out('s=2.0j*pi*fitf')
        out('mybodes='+savedname+'.'+bodename)
        out('bodenames='+str(bodenames))
        out('assert len(bodenames)==len(mybodes)')
        #ng=[(rand()-0.5)/5.0*cg+cg for cg in ig]
        out('initialguesses='+str(initialguesses))
        #out('rg='+str(initialguesses))
        out('')
        out('def mycost(ucv,phaseweight=0.1,returnbl=False):')
        out(ws+'totale=0.0')
        out(ws+'if returnbl:')
        out(2*ws+'bl=[]')
        out(ws+'for cf,cb in zip(funclist,mybodes):')
        out(2*ws+'curc=cf(s,ucv)')
        out(2*ws+'curb=rwkbode.rwkbode(cb.output,cb.input,compin=curc)')
        out(2*ws+'curb.seedfreq=cb.seedfreq')
        out(2*ws+'curb.seedphase=cb.seedphase')
        out(2*ws+'curb.PhaseMassage(fitf)')
        out(2*ws+'magE=squeeze(cb.dBmag())-squeeze(curb.dBmag())')
        out(2*ws+'phaseE=squeeze(cb.phase)-squeeze(curb.phase)')
        out(2*ws+'totale+=sum(magE**2)+phaseweight*sum(phaseE**2)')
        out(2*ws+'if returnbl:')
        out(2*ws+ws+'bl.append(curb)')
        out(ws+'if returnbl:')
        out(2*ws+'return totale,bl')
        out(ws+'else:')
        out(2*ws+'return totale')
        out('')
        out('def RunFit(phaseweight=0.1,ig=None, maxiter=None, maxfun=None):')
        out(ws+'t1=mt()')
        out(ws+'if ig is None:')
        out(2*ws+'ig=initialguesses')
        out(ws+'fitres=fmin(mycost,ig,(phaseweight,), maxiter=maxiter, maxfun=maxfun)')
        out(ws+'t2=mt()')
        out(ws+"print('RunFit time='+str(t2-t1))")
        out(ws+'return fitres')
        out('')
#        out('def RunandPlot(startfig=1,phaseweight=0.1,ig=None):')
#        out(ws+'if ig is None:')
#        out(2*ws+'ig=initialguesses')
#        out(ws+'fitres,junk,finalbodelist, finale=RunFit(phaseweight,ig)')
#        out(ws+'be,initialbodelist=mycost(ig,phaseweight,True)')
#        out(ws+'fi=startfig')
#        out(ws+'for cb,fcb,scb in zip(mybodes,finalbodelist,initialbodelist):')
#        out(2*ws+'rwkbode.GenBodePlot(fi,fitf,cb)')
#        out(2*ws+'rwkbode.GenBodePlot(fi,fitf,scb,clear=False)')
#        out(2*ws+'rwkbode.GenBodePlot(fi,fitf,fcb,clear=False,legend=("Data","Initial","Final",))')
#        out(2*ws+'fi+=1')
#        out(ws+'return fitres, fitf, finalbodelist, finale')
        out('def PlotExp(startfig=1,clear=True):')
        out(ws+'fi=startfig')
        out(ws+'for cb in mybodes:')
        out(2*ws+'rwkbode.GenBodePlot(fi,fitf,cb,clear=clear)')
        out(2*ws+'fi+=1')
        out('')
        out('def PlotRes(ucv,startfig=1,clear=False,legend=None):')
        out(ws+'e, bodelist = mycost(ucv, returnbl=True)')
        out(ws+'fi=startfig')
        out(ws+'plotargs={"clear":clear}')
        out(ws+'if legend is not None:')
        out(2*ws+'plotargs["legend"]=legend')
        out(ws+'for cb in bodelist:')
        out(2*ws+'rwkbode.GenBodePlot(fi,fitf,cb,**plotargs)')
        out(2*ws+'fi+=1')
        out('')
        out('def PlotAll(startfig=1):')
        out(ws+'PlotExp()')
        out(ws+'PlotRes(initialguesses)')
        out(ws+'PlotRes(fitres,legend=["Exp. Data","Initial","Final"])')
        out('')
        out('def RunLongandSave(filename):')
        out(ws+'fitres=RunFit(maxiter=5000,maxfun=5000)')
        out(ws+'write_array(filename,fitres)')
        out(ws+'return fitres')
        out('')
        out('## if __name__=="__main__":')
        out('## '+ws+'fitres = RunFit()')
        out('## '+ws+'PlotAll()')
        out('## '+ws+'pylab.show()')
        pylist.tofile()

    def SymCharDetAll(self,aug=0,sub=1,outpath='chardet.f',**kwargs):
        myoptions=['texname','frlist','basebodename','basecdname','debug','grind','optimize','skipMaxima','ratmx','radexpand']
        texname,frlist,basebodename,basecdname,debug,grind,optimize,skipMaxima,ratmx,radexpand=self.GetValuesorDefaults(kwargs,myoptions)
        mylist=latexlist([],texname)
        ulist=self.SymDefineAllUs(aug,**kwargs)
        preamble=[]
        preamble.append('\\begin{maxima-noout}')
        preamble.append('\tshowtime:all')
        preamble.append('\tnolabels:true')
        preamble.append('\tdeclare (s, complex)')
        if not radexpand:
            preamble.append('\tradexpand:false')
        if grind:
            preamble.append('\tgrind:true')
        preamble.append('\\end{maxima-noout}')
        mylist.extend(preamble)
        mylist.extend(ulist)
        Usyslist=self.SymULatex()
        mylist.extend(Usyslist)
        cdlist=self.SymCharDet(sub=sub,outpath=outpath)
        mylist.extend(cdlist)
        mylist.tofile()
        if not skipMaxima:
            latexoutname=self.RunMaxima(mylist.filename,frlist=frlist)
            print('latexoutname='+latexoutname)
        if not kwargs.has_key('headername'):
            kwargs['headername']='/home/ryan/rwkpython/chardet_header.f'
            kwargs['newoutname']='chard'
        myfinalfnames=self.PrepareFortranFiles([outpath],**kwargs)
        mykwargs={}
        if kwargs.has_key('curvefit'):
            mykwargs['curvefit']=kwargs['curvefit']
            curvefit=kwargs['curvefit']
        else:
            curvefit=False
        if curvefit:
            mypynames,igv=self.PreparePythonFiles([outpath],**mykwargs)
        else:
            mypynames=self.PreparePythonFiles([outpath],**mykwargs)
        self.CleanPythonFiles(mypynames)
        return mylist, latexoutname

    def SymDefineAllUs(self,aug=0,**kwargs):
        myoptions=['texname','frlist','basebodename','basecdname','debug','grind','optimize','skipMaxima','ratmx']
        texname,frlist,basebodename,basecdname,debug,grind,optimize,skipMaxima,ratmx=self.GetValuesorDefaults(kwargs,myoptions)
        mylist=latexlist([],texname)
        for elem in self.elemlist:
            curlines, curdefs, curparams,cursubs=elem.GetMaximaLatexString(aug=aug)
            mylist.extend(curlines)
        return mylist

    def SymBodeAll(self,**kwargs):
        """This method is the new, improved approach to symbolic Bode
        analysis as of April 4, 2006.  It relies on the new symstr
        based GetAugSymMat methods which are typically inherited from
        the base TMMElement class.
        
        This method generates a *.tex file that is an input file for
        the Python+Maxima+LaTeX symbolic engine.  This input file
        defines a transfer matrix for each element in the system as
        well as a system transfer matrix and the bode outputs.  If
        (not skipMaxima), Maxima is called and files will be created
        for each bode output and characteristic determinant.  If
        optimize=True, these output files will be FORTRAN files,
        otherwise they will be *.txt files.

        The output the filename of a latexlist containing the Maxima
        lines, a list of the raw FORTRAN files for the Bodes, and a
        list of raw FORTAN files for the characteristic determinants.""" 
        myoptions=['texname','frlist','basebodename','basecdname','debug','grind','optimize','skipMaxima','ratmx']
        texname,frlist,basebodename,basecdname,debug,grind,optimize,skipMaxima,ratmx=self.GetValuesorDefaults(kwargs,myoptions)
        print('In new sym')
        mylist=latexlist([],texname)
        alldefs=[]
        allparams=[]
        defaultvalues={}
        allsubs=[]
        out=mylist.append
        for elem in self.elemlist:
            curUlines=elem.GetMaximaLines(aug=True)
            mylist.extend(curUlines)
        mylist.extend(self.SymUsys(**kwargs))
        mylist.extend(self.SymSubmat(findsubcol=True))
        mylist.extend(self.SymBaseVect())
        mylist.extend(self.SymRawBodes(frlist, funcname='SymUsys'))
        bodelines,filenames,chardetnames=self.SymBodes(basebodename=basebodename, basecdname=basecdname,debug=debug)
        mylist.extend(bodelines)
        mylist.tofile()
        return mylist.filename, filenames, chardetnames

    def SymUsys(self,sysmatname='Usys',frlist=[],stopind=None,intro=None,showU=False,**unusedargs):
        """This is the new method for outputing the Maxima lines that
        calculate the system transfer matrix, assuming the element
        matrices have already been defined in Maxima/LaTeX.  This
        method replaces SymULatex.  It returns a list of lines to be
        appended to a Maxima/Latex list."""
        Unames=[item.symname for item in self.elemlist]
        if stopind is not None:
            keeplist=Unames[0:stopind+1]
        else:
            keeplist=Unames
        prodstr='.'.join(reverse(keeplist))
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
        if showU:
            latexout.append('or')
            latexout.append('\\begin{maxima}')
            latexout.append("\parseopts{lhs='"+sysmatname+"',wrap=0}")
            latexout.append('\t'+sysmatname+':'+prodstr)
            latexout.append('\\end{maxima}')
        else:
            latexout.append('\\begin{maxima-noout}')
            latexout.append('\t'+sysmatname+':'+prodstr)
            latexout.append('\\end{maxima-noout}')
        return latexout
    
    def SymBodeMaximaAll(self,**kwargs):
#    def SymBodeMaximaAll(self,texname='SymBodeDev.tex',frlist=[],basebodename='bode',basecdname='chardet',debug=False,grind=True,optimize=True,skipMaxima=False,ratmx=False,**unusedargs):
        """This method generates a *.tex file that is an input file
        for the Python+Maxima+LaTeX symbolic engine.  This input file
        defines a transfer matrix for each element in the system as
        well as a system transfer matrix and the bode outputs.  If
        (not skipMaxima), Maxima is called and files will be created
        for each bode output and characteristic determinant.  If
        optimize=True, these output files will be FORTRAN files,
        otherwise they will be *.txt files."""
        myoptions=['texname','frlist','basebodename','basecdname','debug','grind','optimize','skipMaxima','ratmx']
        texname,frlist,basebodename,basecdname,debug,grind,optimize,skipMaxima,ratmx=self.GetValuesorDefaults(kwargs,myoptions)
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
#                print(key+elem.symlabel+'='+str(value))
#        pdb.set_trace()
#        print('allsubs:'+str(allsubs))
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
        if ratmx:
            preamble.append('\tratmx:true')
        preamble.append('\tnolabels:true')
        preamble.append('\t'+ratstr)
        if grind:
            preamble.append('\tgrind:true')
        preamble.append('\\end{maxima-noout}')
#        pdb.set_trace()
        mylist.list=preamble+mylist.list
#        if debug:
#            mylist.tofile()
#            return mylist
        syslines=self.SymULatex(**kwargs)
        mylist.extend(syslines)
        submatlines=self.SymSubmat(findsubcol=True)
        mylist.extend(submatlines)
#        pdb.set_trace()
        bvlines=self.SymBaseVect(**kwargs)
#        bvlines=self.SymBaseVect()
        mylist.extend(bvlines)
        rblines=self.SymRawBodes(frlist)
        mylist.extend(rblines)
#        bodelines,filenames,chardetnames=self.SymBodes(basebodename=basebodename, basecdname=basecdname,debug=debug,optimize=optimize,ratmx=ratmx)
        bodelines,filenames,chardetnames=self.SymBodes(**kwargs)
        mylist.extend(bodelines)
#        return mylist,filenames,chardetnames
        mylist.tofile()
        if debug:
            mylist.tofile()
            return mylist
        if not skipMaxima:
            print('about to ProcessMaxima')
            outname=mylist.ProcessMaxima()
            return outname,filenames,chardetnames
        else:
            return mylist.filename,filenames,chardetnames

    def LoadExpBodeData(self,modulename,datadir):
        if datadir not in sys.path:
            sys.path.insert(0,datadir)
        myexpdata=my_import(modulename)
        expfreq=myexpdata.freq
        allavedata=myexpdata.bodes_ave.bodes
        expbodes=[]
        for curbode in self.bodeouts:
            curexpbode=rwkbode.FindMatch(allavedata,curbode.output,curbode.input)
            expbodes.append(curexpbode)
        self.expbodes=expbodes
        myind=sys.path.index(datadir)
        sys.path.pop(myind)
        return expfreq, expbodes

    def ExtractInitialGuesses(self):
        initalguesses=[]
        alldefs, allparams, defaultvalues, unknownparams, allsubs=self.GetParamsandDefs()
        for ukp in unknownparams:
            myig=defaultvalues[ukp]
            if shape(myig):
                if max(shape(myig))==1:
                    myig=myig[0]
                if shape(myig):
                    if max(shape(myig))==1:
                        myig=myig[0] 
            initalguesses.append(myig)
        return initalguesses

    def GetUnkownParams(self):
        unknownparams=[]
        for elem in self.elemlist:
            if elem.unknownparams:
                curucv=elem.unknownparams
                curlabel=elem.symlabel
                if curucv=='all':
                    curucv=elem.params.keys()
                for uc in curucv:
                    unknownparams.append(uc+curlabel)
        for x,cbo in enumerate(self.bodeouts):
            curname='gainbode'+str(x)
            if not cbo.gainknown:
                unknownparams.append(curname)
        return unknownparams

    def GetCompensators(self):
        comps=[]
        for elem in self.elemlist:
            curcomp=elem.compensators
            curdict=elem.symparams
            if curcomp:
                if type(curcomp)==list:
                    curents = [str(curdict[key]) for key in curcomp]
                    comps.extend(curents)
                else:
                    comps.append(str(curdict[curcomp]))
        return comps
                

    def SubstituteUnknowns(self, unknowndict):
        sublist=[]
        for elem in self.elemlist:
            if elem.SubstituteUnknowns(unknowndict):
                sublist.append(elem)
        for x,cbo in enumerate(self.bodeouts):
            curname='gainbode'+str(x)
            if not cbo.gainknown and unknowndict.has_key(curname):
                cbo.gain=unknowndict[curname]
                cbo.gainknown=True
        return sublist

    def BuildUnknownDict(self, ucv):
        unknownparams=self.GetUnkownParams()
        return dict(zip(unknownparams,ucv))

    def NewPythonFilewithSubs(self, filname, unknowndict, sublist, prefix='my',headername='myheader.py'):
        outlist=[]
        subsysnames=[]
        symlabellist=[item.symlabel for item in sublist]
#        Pdb().set_trace()
        for item in self.inputlist:
            if isinstance(item, TMMSubSystem):
                mylines, junk=inspect.getsourcelines(item.__class__)
                classname=GetClassNamefromInspectList(mylines)
                newname=prefix+classname
                subsysnames.append((classname,newname))
                mylines=ReplaceClassName(mylines,classname,newname)
                newlines=[ProcessOneSubstitutionLine(line, symlabellist, unknowndict) for line in mylines]
                newlines=[CleanInitLine(line) for line in newlines]
                outlist.extend(newlines)
                outlist.append('\n')
        maincode, junk=inspect.getsourcelines(self.__class__)
        classname=GetClassNamefromInspectList(maincode)
        newname=prefix+classname
        maincode=ReplaceClassName(maincode,classname,newname)
        newlines=[ProcessOneSubstitutionLine(line, symlabellist, unknowndict) for line in maincode]
        newlines=[CleanInitLine(line) for line in newlines]
        newlines=[ReplaceSubSystemCalls(line, subsysnames) for line in newlines]
        gaindict=GetUnknownBodeouts(newlines,unknowndict)
        newlines=[ProcessBodeoutLine(line,gaindict) for line in newlines]
        outlist.extend(newlines)
        f=open(headername,'r')
        hlines=f.readlines()
        hlines.append('\n')
        outlines=hlines+outlist
        f.close()
        f=open(filname,'w')
        f.writelines(outlines)
        f.close()
        return outlines

    def GetParamsandDefs(self,curvefit=1,otherparams={},prefersubs=0,aug=0):
        alldefs=[]
        allparams=[]
        defaultvalues={}
        unknownparams=[]
#        initalguesses=[]
        allsubs=[]
#        pdb.set_trace()
        for elem in self.elemlist:
            curlines, curdefs, curparams,cursubs=elem.GetMaximaLatexString(aug=aug)
            if prefersubs and hasattr(elem,'GetMaximaSubstitutions'):
                curdefs=elem.GetMaximaSubstitutions(aug=aug)
#            mylist.extend(curlines)
            if curvefit and elem.unknownparams:
                curucv=elem.unknownparams
                curlabel=elem.symlabel
                if curucv=='all':
#                    pdb.set_trace()
                    curucv=elem.params.keys()
                for uc in curucv:
                    unknownparams.append(uc+curlabel)
            alldefs+=curdefs
            allparams+=curparams
            allsubs+=cursubs
            for key,value in elem.params.iteritems():
                defaultvalues[key+elem.symlabel]=value
        for x,cbo in enumerate(self.bodeouts):
            curname='gainbode'+str(x)
            allparams.append(curname)
            defaultvalues[curname]=cbo.gain
            if not cbo.gainknown:
                unknownparams.append(curname)
#                initalguesses.append(cbo.gain)
        for key,value in otherparams.iteritems():
            defaultvalues[key]=value
#        return alldefs, allparams, defaultvalues, unknownparams, initalguesses, allsubs
#        for curdef in alldefs:
#            print('curdef:'+str(curdef))
        return alldefs, allparams, defaultvalues, unknownparams, allsubs

    def SymBodesParseMaxima(self,filenames,chardetnames,functionparams={},otherparams={},debug=False,grind=True,optimize=True,curvefit=False,ratmx=False,**unusedargs):
        aug=1
        pyfilenames=[]
        alldefs=[]
        allparams=[]
        defaultvalues={}
        unknownparams=[]
        initalguesses=[]
        allsubs=[]
#        pdb.set_trace()
        for elem in self.elemlist:
            curlines, curdefs, curparams,cursubs=elem.GetMaximaLatexString(aug=aug)
#            mylist.extend(curlines)
            if curvefit and elem.unknownparams:
                curucv=elem.unknownparams
                curlabel=elem.symlabel
                if curucv=='all':
#                    pdb.set_trace()
                    curucv=elem.params.keys()
                for uc in curucv:
                    unknownparams.append(uc+curlabel)
                    myig=elem.params[uc]
                    if shape(myig):
                        if max(shape(myig))==1:
                            myig=myig[0]
                        if shape(myig):
                            if max(shape(myig))==1:
                                myig=myig[0] 
                    initalguesses.append(myig)
            alldefs+=curdefs
            allparams+=curparams
            allsubs+=cursubs
            for key,value in elem.params.iteritems():
                defaultvalues[key+elem.symlabel]=value
        for x,cbo in enumerate(self.bodeouts):
            curname='gainbode'+str(x)
            if cbo.gainknown:
                allparams.append(curname)
                defaultvalues[curname]=cbo.gain
            else:
                unknownparams.append(curname)
                initalguesses.append(cbo.gain)
        for key,value in otherparams.iteritems():
            defaultvalues[key]=value
#        pdb.set_trace()
        for curfile, curcdfile in zip(filenames,chardetnames): 
            fno,junk=os.path.splitext(curfile)
            cdfno,junk=os.path.splitext(curcdfile)
            bodenum=-1
            if fno.find('bode')>=0:
                start=fno.find('bode')
                bodenum=int(fno[start+4:])
            pyfilename=fno+'.py'
            pyfilenames.append(pyfilename)
            if optimize:
#                pdb.set_trace()
                cdout,retval=FortranToTextList(curfile,'bode')
                if ratmx:
                    cdout2,retval2=FortranToTextList(curcdfile,'chardet')
            else:
                cdin=rwkreadfile(curfile)
                cdout=[]
#                inputs=[cdin,cdin2]
                inputs=[cdin]
#                outputs=[cdout,cdout2]
                outputs=[cdout]
                if ratmx:
                    cdin2=rwkreadfile(curcdfile)
                    cdout2=[]
                    inputs.append(cdin2)
                    outputs.append(cdout2)
                for curin,curout in zip(inputs,outputs):
                    for x,line in enumerate(curin):
                        if line:
                            if x<(len(curin)-1):
                                curout.append(line+' \\')
                            else:
                                curout.append(line)
                rv,dummy=cdout[0].split('=',1)
                cdout.append('return '+rv)
                if ratmx:
                    rv2,dummy=cdout2[0].split('=',1)
                    cdout2.append('return '+rv2)
            pylines=textlist([],pyfilename)
            cdlines=[]
            pylines.readfile('header.py')
            funcstr='s'
            if curvefit:
                funcstr+=',ucv'
            for key,value in functionparams.iteritems():
                if key not in unknownparams:
                    funcstr+=','+str(key)+'='+str(value)
            pylines.append('def '+fno+'('+funcstr+'):')
            cdlines.append('def '+cdfno+'('+funcstr+'):')
#            pdb.set_trace()
            for x,uc in enumerate(unknownparams):
                gainbodenum=-1
                skip=False
                if uc.find('gainbode')==0:
                    gainbodenum=int(uc[8:])
                if bodenum>-1 and gainbodenum>-1:
                    if bodenum!=gainbodenum:
                        skip=True
                if not skip:
                    curline='\t'+uc+'='+'ucv['+str(x)+']'
                    pylines.append(curline)
                    cdlines.append(curline)
            for cp in allparams:
                if (cp not in functionparams.keys()) and (cp not in unknownparams):
                    gainbodenum=-1
                    skip=False
                    if cp.find('gainbode')==0:
                        gainbodenum=int(cp[8:])
                    if bodenum>-1 and gainbodenum>-1:
                        if bodenum!=gainbodenum:
                            skip=True
                    if not skip:
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
            if ratmx:
                for line in cdout2:
                    line=line.replace('^','**')
                    cdlines.append('\t'+line)
                pylines.append('')
                pylines.extend(cdlines)
            pylines.tofile()
        if curvefit:
            return pyfilenames, initalguesses
        else:
            return pyfilenames

    def SymBodesAllinOne(self,texname='SymBodeDev.tex',frlist=[],basebodename='bode',basecdname='chardet',functionparams={},otherparams={},debug=False,grind=True,optimize=True,skipMaxima=False,ratmx=False):
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
#                print(key+elem.symlabel+'='+str(value))
#        pdb.set_trace()
#        print('allsubs:'+str(allsubs))
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
        if grind:
            preamble.append('\tgrind:true')
        preamble.append('\\end{maxima-noout}')
        mylist.list=preamble+mylist.list
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
        mylist.tofile()
        if debug:
            mylist.tofile()
            return mylist
        if not skipMaxima:
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
            if optimize:
                cdout,retval=FortranToTextList(curfile,'bode')
                cdout2,retval2=FortranToTextList(curcdfile,'chardet')
            else:
                cdin=rwkreadfile(curfile)
                cdin2=rwkreadfile(curcdfile)
                cdout=[]
                cdout2=[]
                inputs=[cdin,cdin2]
                outputs=[cdout,cdout2]
                for curin,curout in zip(inputs,outputs):
                    for x,line in enumerate(curin):
                        if line:
                            if x<(len(curin)-1):
                                curout.append(line+' \\')
                            else:
                                curout.append(line)
                rv,dummy=cdout[0].split('=',1)
                cdout.append('return '+rv)
                rv2,dummy=cdout2[0].split('=',1)
                cdout2.append('return '+rv2)
            pylines=textlist([],pyfilename)
            cdlines=[]
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

    def SymBodes(self,basebodename='bodeoutputs',basecdname='chardet',myext='.txt',debug=False,optimize=True,ratmx=False,showbode=False,sub=True,**unusedargs):
        """This method converts rawbodeouts (which are just the states
        at specified model locations) into actual model simulations of
        measured signals.  This is done by calculating the differnce
        between to rawbodeouts for differential sensors, by
        multiplying by gains when needed, and by multiplying by s or
        s**s for velocity or acceleration sensors."""
        if optimize:
            myext='.f'
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
            linesout.append('A gain is introduced to scale the bode output to engineering units')
            linesout.append('\\begin{equation}')
            linesout.append('\t'+bodestr+'='+bodestr+'*gain'+bodestr)
            linesout.append('\\end{equation}')
            linesout.append('\\begin{maxima-noout}')
            linesout.append('\t'+bodestr+':'+bodestr+'*gain'+bodestr)
            linesout.append('\\end{maxima-noout}')
            if showbode:
                substr=bodestr+'_sub'
                if sub:
                    sublines=self.MaximaSubstitute(bodestr,substr)
                    linesout.extend(sublines)
                else:
                    linesout.append('\\begin{maxima-out}')
                    linesout.append('\t'+substr+':'+bodestr)
                    linesout.append('\\end{maxima-out}')
                linesout.append('\\begin{maxima-out}')
                linesout.append('\t'+substr+':radcan('+substr+')')
                linesout.append('\\end{maxima-out}')
                linesout.append('\\begin{maxima-out}')
                linesout.append('\t'+substr+':factor('+substr+')')
                linesout.append('\\end{maxima-out}')
                linesout.append('\\begin{maxima}')
                linesout.append("\t\\parseopts{lhs='"+bodestr+"'}")
                linesout.append("\t"+substr)
                linesout.append('\\end{maxima-noout}')
            if not debug:
                linesout.append('\\begin{maxima-noout}')
    #            linesout.append('\tstringout("'+curname+'",bode=radcan('+bodestr+'))')
                if optimize:
                    linesout.append('\twith_stdout ("'+curname+'", fortran_optimize ('+bodestr+'))')
                    if ratmx:
                        linesout.append('\twith_stdout ("'+curcdname+'", fortran_optimize (ratdenom('+bodestr+')))')
                else:
                    linesout.append('\tstringout("'+curname+'",bode='+bodestr+')')
                    if ratmx:
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

    def SymRawBodes(self, frlist=[], funcname='SymULatex'):
        """This function calculates the raw bode outputs for a
        TMMSystem.  It assumes that the base vector bv is already
        defined, so you must have already called SymBaseVect.  This
        function returns a LaTeX list that calculates the matrices U_i
        that transfer from the base vector to the raw bode locations
        as well as the raw bode outputs themselves where

        rb_i=U_i*bv

        By specifying funcname='SymUsys', this function is compatible
        with the new symbolic functions."""

        linesout=[]
        myfunc=getattr(self,funcname)
        for x,rb in enumerate(self.rawbodeouts):
            curp=rb.pind
            pstr=str(curp)
            curintro='The transfer matrix from the base vector to model position '+pstr+' is'
            curUname='U'+pstr
            curUlines=myfunc(intro=curintro,stopind=curp,sysmatname=curUname,frlist=frlist)
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

    def SymBaseVect(self,N=None,showall=False,ratmx=False,**unusedargs):
        """This function finds the base vector used to calculate the
        system bode response.  It assumes that SymSubmat has already
        been called and will append lines using submat and subcol,
        assuming they have already been defined by SymSubmat called
        with findsubcol=True.

        This funciton returns lines to be append to a LaTeX list.
        These lines create a symbolic solution for the base vector."""
        if N is None:
            N=self.maxsize
        linesout=[]
        linesout.append("The non-zero portion of the base vector will be found using")
        linesout.append('\\begin{equation}')
        linesout.append('\\M{nzbv}=\\M{submat}^{-1}\\left(-\\M{subcol}\\right)')
        linesout.append('\\end{equation}')
        if showall:
            linesout.append('or')
            linesout.append('\\begin{maxima}')
            linesout.append("\t\\parseopts{lhs='\\M{nzbv}'}")
            linesout.append('\tnzbv:invert(submat).(-1*subcol)')
            linesout.append('\\end{maxima}')
        else:
            linesout.append('\\begin{maxima-noout}')
            linesout.append('\tnzbv:invert(submat).(-1*subcol)')
            linesout.append('\\end{maxima-noout}')
        linesout.append('\\begin{maxima-noout}')
        linesout.append('bv:zeromatrix('+str(N+1)+',1)')
        linesout.append('bv['+str(N+1)+',1]:1')
        for x,r in enumerate(self.bcend):
            linesout.append('bv['+str(r+1)+',1]:nzbv['+str(x+1)+',1]')
        if ratmx:
            linesout.append('bv:radcan(bv)')
        linesout.append('\\end{maxima-noout}')
        if showall:
            linesout.append('\\begin{maxima}')
            linesout.append("\t\\parseopts{lhs='\\M{bv}'}")
            linesout.append("\tbv")
            linesout.append('\\end{maxima}')
        return linesout

    def SymSubmat(self, sysmatname='Usys',findsubcol=False,N=None,showsub=False,sub=False):
        """Note that you must call SymUsys or SymULatex first.  This
        function will return lines that are to be appended to a
        latexlist that already includes the SymUsys or SymULatex
        lines.  This function uses the system boundary conditions to
        find a sub-matrix whose determinant is the characteristic
        equation of the system."""
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
        linesout.append('\\begin{maxima-noout}')
        linesout.append('\t'+submatstr)
        linesout.append('\\end{maxima-noout}')
        if showsub:
            if sub:
                sublines=self.MaximaSubstitute('submat')
                linesout.extend(sublines)
            linesout.append('\\begin{maxima}')
            linesout.append("\t\\parseopts{lhs='\M{subU}',wrap=0}")
            linesout.append('\tsubmat')
            linesout.append('\\end{maxima}')
        else:
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

    def SymCharDet(self, sysmatname='Usys',outpath='chardet.f',sub=True,showsub=True,radcan=False):
        """This function calls SymSubmat to find the submatrix and
        then calculates its determinant to find the characteristic
        determinant of the system.

        Note that you must call SymULatex first."""
        linesout=self.SymSubmat(sysmatname=sysmatname,showsub=showsub,sub=sub)
        linesout.append('and the characteristic determinant is')
        linesout.append('\\begin{maxima-noout}')
        linesout.append('\tchardet:determinant(submat)')
        linesout.append('\\end{maxima-noout}')   
        if sub:
            sublines=self.MaximaSubstitute('chardet')
            linesout.extend(sublines)
        linesout.append('\\begin{maxima}')
        linesout.append("\t\\parseopts{lhs='cd',wrap=0}")
        linesout.append('\tchardet')
        linesout.append('\\end{maxima}')
        if radcan:
            linesout.append('Simplfying gives')
            linesout.append('\\begin{maxima}')
            linesout.append("\t\\parseopts{lhs='cd',wrap=0}")
            linesout.append('\tchardet:radcan(chardet)')
            linesout.append('\\end{maxima}')
        linesout.append('\\begin{maxima-noout}')
#        linesout.append('\tstringout("'+outpath+'",cd=chardet)')
        linesout.append('\twith_stdout ("'+outpath+'", fortran_optimize (chardet))')
        linesout.append('\\end{maxima-noout}')   
        return linesout
                    
    def MaximaSubstitute(self,namebefore,nameafter=''):
        linesout=[]
        if not nameafter:
            nameafter=namebefore
        defs, allparams, defaultvalues, unknownparams, allsubs =self.GetParamsandDefs(curvefit=0,prefersubs=1)
        substr=', '.join(defs)
        substr='['+substr+']'
        linesout.append('\\begin{maxima-noout}')
        linesout.append(nameafter+':at('+namebefore+','+substr+')')
        linesout.append('\\end{maxima-noout}')   
        return linesout

    def SymULatex(self,sysmatname='Usys',frlist=[],stopind=None,intro=None,showU=False,**unusedargs):#omit me
        """This function returns lines for a LaTeX list that calculate
        the system transfer matrix.  By setting stopind, it can also
        be used to calculate transfer matrices that calculate bode
        outputs at certain locations based on the base vector."""
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
        XMlist=[[0,2],[4,6],[8,10]]
        if isinstance(bcstr,str):
            if bcstr=='free':
                bcelist=VMlist
            elif bcstr=='fixed':
                bcelist=XTlist
            elif bcstr=='pinned':
                bcelist=XMlist
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
        prevHT=scipy.eye(4,dtype='d')
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
        basemesh['mesh']=zeros((1,3),'d')
        self.sysmesh.append(basemesh)
#        ne=len(self.elemlist)
        if not self.sysHT:
            self.CreateSysHT(beammesh)
        mesh=scipy.zeros((1,3),'d')
#        Pdb().set_trace()
        prevHT=eye(4,dtype='d')
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
                #tempHT[0:3,3]=zeros((3,),'d')
                tempHT[0:3,3]=prevHT[0:3,3]
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
            prevHT=curHT
        return mesh

    def FindModeShape(self, eigval,beammesh=10,logtex=0,fmt='0.5g',scale=True,modenum=0):
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
        nc=shape(bv)[1]
        if nc==1:
            bv=bv[:,0]
        else:
            bv=bv[:,modenum]
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
#        prevU=scipy.eye(self.maxsize,'f')+0j
        prevU=scipy.eye(self.maxsize,dtype='D')
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
#                    curZ=curZ[:,0]
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
                    curdisps=rowwise(curdisps)
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
#                curZ=curZ[:,0]
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
#        dmax=matmaxabs(real(disps))#*-1.#undo the -w TMM sign configuration for FEA comparison
#        amax=matmaxabs(real(angles))
#        dmaxi=matmaxabs(imag(disps))#*-1.#undo the -w TMM sign convention for FEA comparison
#        amaxi=matmaxabs(imag(angles))
        scalemat=vstack([real(disps),imag(disps),real(angles),imag(angles)])
        myindex=(abs(scalemat)).argmax()
        maxelem=scalemat.flat[myindex]
#        dmax=real(disps).max()
#        amax=real(angles).max()
#        dmaxi=imag(disps).max()
#        amaxi=imag(angles).max()
        if scale:
#            tmax=matmaxabs(array([dmax,amax,dmaxi,amaxi]))
            tmax=maxelem
            disps=disps*0.2/tmax
            angles=angles*0.2/tmax
            for curmode in modedict['bodies']:
                curmode['disps']=curmode['disps']*0.2/tmax
                curmode['angles']=curmode['angles']*0.2/tmax
        if logtex:
            texlist.tofile()
        return disps,angles,modedict

    def FindBaseVect(self, eigval,usesvd=True):
#        pdb.set_trace()
        submat=self.FindSubmat(eigval)
#        ns=null(submat)
        nstol=1e-12
        if usesvd:
            ns=null(submat)
        else:
            ns=rrefnull(submat,nstol)
#        if len(scipy.shape(ns))>1:
#            nsvect=ns[:,0]#*1e5#assume distinct eigenvalues/vectors and only 1D null space
#        else:
#            nsvect=ns
        basevect=zeros((self.maxsize,shape(ns)[1]),'D')
        ns=colwise(ns)
        for c,nsvect in enumerate(ns.transpose()):#this is broke because it returns individual elements if ns is a column vector
            for curind,curns in zip(self.bccols,nsvect):
                basevect[curind,c]=curns#insert non-zero values into the boundary condition vector at the base of the beam
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
            curitem.compvect=zeros(len(fvect),dtype='D')
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
#                curz=dot(U,bv)
                curz=dot(U,bv)
                for curraw in self.rawbodeouts:
                    if curraw.pind==x:
#                        pdb.set_trace()
#                        curraw.compvect[r]=curz[curraw.dof,0]#curz is 2d but has only 1 column - the 0 is to get a scalar from the only column
                        curz=curz.squeeze()
                        curraw.compvect[r]=curz[curraw.dof]
#                    outcol+=1
#        return outmat
#        return rawbodedict
        bodes=[]
        for bodedict in self.bodeouts:
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
            tempcomp=tempcomp*bodedict.gain
#            print('gain='+str(bodedict.gain))
            bodes.append(rwkbode.rwkbode(bodedict.output,bodedict.input,compin=tempcomp,seedfreq=bodedict.seedfreq,seedphase=bodedict.seedphase))
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
        basevect=scipy.zeros((N,1),dtype='D')
#        basevect=basevect+0.0j
        basevect[-1]=1
        tempvect=dot(inverse(submat),-augcol)
        for curind,curns in zip(self.bccols,tempvect):
            basevect[curind,0]=curns
        return basevect

    def FindEig(self, guess,mytol=1e-14,maxfun=1e4,maxiter=1e3,returncomplex=False,useabs=True):
        if useabs:
            if scipy.shape(guess)==():
                guess=scipy.array([scipy.real(guess),scipy.imag(guess)])
            eigout=scipy.optimize.fmin(self.EigError,guess,xtol=mytol,ftol=mytol,maxfun=maxfun,maxiter=maxiter)
        else:
            eigout=scipy.optimize.newton(self.EigError,guess,tol=mytol,maxiter=maxiter,args=(False,))
        if returncomplex:
            if shape(eigout):
                eigout=eigout[0]+eigout[1]*1.0j
        return eigout

    def EigError(self, value,useabs=True):
        if not shape(value):
            value=complex(value)
        submat=self.FindSubmat(value)
        shapevect=scipy.shape(submat)
        if shapevect:
            if len(shapevect)>1 and max(shapevect)>1:
                chardet=scipy.linalg.det(submat)
                if useabs:
                    return abs(chardet)
                else:
                    return chardet
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

    def FindU(self,value):
        if scipy.shape(value):
            value=value[0]+value[1]*1.0j
        U=scipy.eye(self.maxsize,dtype='D')
        for curelem in self.elemlist:
            tempU=curelem.GetMat(value)
            U=scipy.dot(tempU,U)
        return U
    
    def FindAugU(self,value):
        if scipy.shape(value):
            value=value[0]+value[1]*1j
        U=scipy.eye(self.maxsize+1,dtype='D')
        for curelem in self.elemlist:
            tempU=curelem.GetAugMat(value)
            U=scipy.dot(tempU,U)
        return U

    def FindSubmat(self,value,eps=1e-14):
        U=self.FindU(value)
        N=self.maxsize
        n=N/2
        matout=scipy.zeros((n,n),'D')
        for r,curri in enumerate(self.bcend):
            for c,curci in enumerate(self.bccols):
                matout[r,c]=U[curri,curci]
        testmat=imag(matout)<=eps
        if testmat.all():
            return real(matout)
        else:
            return matout

    def FindSubmatwAugcol(self,value):
#        pdb.set_trace()
        U=self.FindAugU(value)
        N=self.maxsize
        n=N/2
        matout=scipy.zeros((n,n),'D')#+0j
        colout=scipy.zeros((n,1),'D')#+0j
        for r,curri in enumerate(self.bcend):
            for c,curci in enumerate(self.bccols):
                matout[r,c]=U[curri,curci]
            colout[r,0]=U[curri,N]
        return matout, colout

    def FindAugSubmat(self,value):
        U=self.FindAugU(value)
        N=self.maxsize
        n=int(N/2)
        matout=scipy.zeros((n,n),'D')#+0j
        vectout=scipy.zeros((n,),'D')#+0j
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
        matout=zeros((self.maxsize+1,nw),'D')#+0j
        matout[self.maxsize,:]=1
        for n,s in enumerate(svect):
            subU,vect=self.FindAugSubmat(s)
            curb=dot(inverse(subU),-vect)
            for curcol,curent in zip(self.bccols,curb.tolist()):
                matout[curcol,n]=curent
        return matout

    def GetValuesorDefaults(self,dictin,keylist):
        """This is a convience function for passing values
        to lots of different subfunctions.  dictin is used
        to update the defaults dict specified in this function.
        keylist is then used to retreive the desired values
        from the update dict.  The return value is the list
        of values."""
        defaults={'N':self.maxsize,'texname':'SymBodeDev.tex','frlist':[],'basebodename':'bode','basecdname':'chardet','debug':False,'grind':True,'optimize':True,'skipMaxima':False,'ratmx':False,'radexpand':False}
        defaults.update(dictin)
        argsout=[defaults[key] for key in keylist]
        return argsout

class ClampedFreeTMMSystem(TMMSystem):
    def __init__(self, elemlist, bcend='free', bcbase='fixed', \
                 bodeouts=[]):
        TMMSystem.__init__(self, elemlist, bcend=bcend, \
                           bcbase=bcbase, bodeouts=bodeouts)


def Getkwargs(dictin,keylist):
    dictout={}
    for curkey in keylist:
        if dictin.has_key(curkey):
            dictout[curkey]=dictin[curkey]
    return dictout

def null(A, eps=1e-10):
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
        mat2=multiply(x.real,x.real) + multiply(x.imag,x.imag)
        return mat2
    else:
        return multiply(x,x)

def matmax(matin):
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

def CleanInitLine(linein):
    if linein.find('def __init__')==-1:
        return linein
    p1,rest=linein.split('(',1)
    middle,p2=rest.split(')')
    return p1+'(self)'+p2

def ReplaceSubSystemCalls(linein, subsysnames):
    lineout=copy.copy(linein)
    for oldname,newname in subsysnames:
        if lineout.find(oldname)>-1:
            lineout=re.sub(oldname+'\(.*?\)',newname+'()',lineout)
    return lineout
                        

def GetClassNamefromInspectList(listin):
    line0=listin[0]
    ind=line0.find('class')
    i2=ind+len('class')+1
    temp=line0[i2:]
    nameout,junk=temp.split('(',1)
    return nameout

def ReplaceClassName(listin, oldname, newname):
    listout=copy.copy(listin)
    line0=listout[0]
    newline=line0.replace(oldname,newname)
    listout[0]=newline
    return listout

def ProcessOneSubstitutionLine(linein, symsublist, unknowndict):
    p=re.compile("symlabel='(.*?)'")
    q=p.search(linein)
    if not q:
        return linein
    mysymlabel=q.group(1)
    if mysymlabel not in symsublist:
        return linein
    print('linein='+str(linein))
    myopts=GetPythonFunctionArgs(linein)
    filtopts=[item for item in myopts if item.find('unknownparams')==-1]
    myparams=filtopts.pop(0)
    assert myparams[0]=='{' and myparams[-1]=='}'
    myparams=myparams[1:-1]
    templist=myparams.split(',')
    keylist=[]
    for item in templist:
        ind=item.find(':')
        key=item[0:ind]
        key=key.strip()
        assert key[0]=="'" and key[-1]=="'"
        key=key[1:-1]
        keylist.append(key)
    paramstr=''
    for key in keylist:
        if paramstr:
            paramstr+=', '
        curkey=key+mysymlabel
        curval=unknowndict[curkey]
        paramstr+="'"+key+"':"+str(curval)
    filtopts.insert(0,'{'+paramstr+'}')
    optstr=', '.join(filtopts)
    p1,rest=linein.split('(',1)
    junk,p2=rest.split(')',1)
    lineout=p1+'('+optstr+')'+p2
    print('lineout='+str(lineout))
    return lineout


def GetBodeouts(listin):
    p=re.compile('__init__.*bodeouts=\[(.*?)\]')
    for line in listin:
        q=p.search(line)
        if q:
            mystr=q.group(1)
            mylist=mystr.split(',')
            outlist=[item.strip() for item in mylist]
            return outlist
    return []

def GetUnknownBodeouts(listin, unknowndict):
    bodeouts=GetBodeouts(listin)
    ubodeouts={}
    for x, bode in enumerate(bodeouts):
        curkey='gainbode'+str(x)
        if unknowndict.has_key(curkey):
            ubodeouts[bode]=unknowndict[curkey]
    return ubodeouts
    
def ProcessBodeoutLine(linein, gaindict):
    temp=rwkstr(linein)
    if not (temp.contains('bodeout') and (temp.contains('gainknown=False') or temp.contains('gainknown=0'))):
        return linein
    curbodeout,rest=linein.split('=',1)
    curbodeout=curbodeout.strip()
    if not gaindict.has_key(curbodeout):
        return linein
    p1,rest=linein.split('(',1)
    middle,p2=rest.split(')',1)
    myargs=GetPythonFunctionArgs(linein)
    p=re.compile('gain=.*')
    myargs2=[p.sub('gain='+str(gaindict[curbodeout]),line) for line in myargs]
    
    myargs3=[item for item in myargs2 if item.find('gainknown=')==-1]
    newargstr=', '.join(myargs3)
    print('linein='+str(linein))
    lineout=p1+'('+newargstr+')'+p2
    print('lineout='+str(lineout))
    return lineout
