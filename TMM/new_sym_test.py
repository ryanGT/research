from __future__ import division
#import TMM
#reload(TMM)
from scipy import sinh, cosh, sin, cos, real, imag, shape, arange, pi, zeros, array, eye, transpose, dot, conj, c_, poly, sqrt, vectorize, dot, randn, rand, squeeze, poly1d, logspace, log10
from LinearAlgebra import inverse
#import pylab
from pylab import show, ioff, figure, xlabel, ylabel
import rootlocus
reload(rootlocus)
import rwkbode
#reload(rwkbode)
from scipy.optimize import newton, fmin, fminbound
from scipy.io import read_array, save
ioff()
from textfiles.latexlist import RunLatex, latexlist
from textfiles import textlist
from rwkmisc import my_import
import sys, os, copy, time

from rwkmplutil import SetPlotLims
from rwkdataproc import thresh

from IPython.core.debugger import Pdb

import rwkreportgen
reload(rwkreportgen)

import glob, re

ioff()

modeldir='/home/ryan/thesis/python_files'
if modeldir not in sys.path:
    sys.path.insert(1,modeldir)

import curvefit_samii_model as cfsm
#reload(cfsm)
olmodel=cfsm.olsamiimodel_withig()

olddir='/home/ryan/thesis/sym_control_design'
if olddir not in sys.path:
    sys.path.insert(1,olddir)

longtdir='/home/ryan/thesis/statistical_analysis/longT'
if longtdir not in sys.path:
    sys.path.insert(2,longtdir)

fitdir='/home/ryan/thesis/swdesign/integrated_id/April_26_redo'
if fitdir not in sys.path:
    sys.path.insert(3,fitdir)

myreport=rwkreportgen.Report('dumb_fit.tex')
myreport.SetModelFreqVect(0.1,50,0.01)
myreport.SetFreqLim([0.1,30])

mypairs=[('a1','j2v'),('j2a','j2v')]

myreport.LoadBodeDataSet('no_afb_1_5deg_ave','fullsweep','Exp.',pairstokeep=mypairs)

myreport.LoadBodeDataSet('stat_downsample','statexp','Statistical Analysis Data',pairstokeep=mypairs)

ioff()

#mylist, bodenames, chardetnames = olmodel.SymBodeAll(texname='new_sym.tex',basebodename='new_sym_bode')
#mylist.tofile()
#outname=olmodel.RunMaxima('new_sym.tex')
#outname=olmodel.RunMaxima('new_sym_test.tex')
#RunLatex(outname)
runf=1
if runf:
    olmodel.CreateFortranandPythonModules('new_sym')
    olmodel.CreateFortranandPythonModules('old_sym',newsym=False)
#    myfinalfnames=olmodel.PrepareFortranFiles(bodenames)
#    mypynames=olmodel.PreparePythonFiles(bodenames)
#    fmodnames=olmodel.CallF2py(myfinalfnames)



#olmodel.SymBodeMaximaAll(texname='old_sym.tex',basebodename='old_sym_bode')

dfcomp=read_array(os.path.join(fitdir,'dumbfit_comp_w1.txt'))
dfdb=read_array(os.path.join(fitdir,'dumbfit_db.txt'))

nsp='new_sym'
osp='old_sym'
nl='New Sym'
ol='Old Sym'
myreport.NewBodeDataSetFromFortranModels(nsp, 'new_sym1', nl ,bodenums=[0,1],ucv=dfcomp,expkey='fullsweep')
myreport.NewBodeDataSetFromFortranModels(osp, 'old_sym1', ol ,bodenums=[0,1],ucv=dfcomp,expkey='fullsweep')

myreport.NewBodeDataSetFromFortranModels(nsp, 'new_sym2', nl ,bodenums=[0,1],ucv=dfdb,expkey='fullsweep')
myreport.NewBodeDataSetFromFortranModels(osp, 'old_sym2', ol ,bodenums=[0,1],ucv=dfdb,expkey='fullsweep')

key1='test'
myreport.GenListofBodeFigs(key1,['fullsweep','new_sym1','old_sym1','new_sym2','old_sym2'],[3,3],[212,212],maglims=[[-40,5],[-25,25]],phaselims=[[-400,200],[-220,-50]],magticks=[arange(-40,5,10),arange(-20,25,10)],phaseticks=[arange(-360,220,90),arange(-220,-55,40)])

myreport.PlotListofBodeFigs(1, key1)

show()
