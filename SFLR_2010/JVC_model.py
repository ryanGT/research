from pylab import *
from scipy import *

import pylab_util as PU

import exp_data
reload(exp_data)
import sympy_bode_analysis
reload(sympy_bode_analysis)
import rwkbode, txt_mixin
from rwkmisc import my_import
#from sympy_bode_analysis import calc_and_plot_Bodes
from sympy_optimize_utils import myoptimize, _cost, fexp
from rwkdataproc import thresh

import sympy_utils
reload(sympy_utils)

ind1 = thresh(fexp,5)

from sympy import Symbol
import sympy_TMM
reload(sympy_TMM)

import os, copy, time

from IPython.Debugger import Pdb

mu = Symbol('mu')
EI = Symbol('EI')
L1 = Symbol('L1')
beta1 = Symbol('beta1')
beam_params1 = {'mu':mu, 'EI':EI, 'L':L1, 'beta':beta1}
beam1 = sympy_TMM.Sympy_Beam_Element(beam_params1, label='_1')

L2 = Symbol('L2')
beta2 = Symbol('beta2')
beam_params2 = {'mu':mu, 'EI':EI, 'L':L2, 'beta':beta2}
beam2 = sympy_TMM.Sympy_Beam_Element(beam_params2, label='_2')

a_m = Symbol('a_m')
a_L = Symbol('a_L')
a_r = Symbol('a_r')
a_I = Symbol('a_I')
a_gain = Symbol('a_gain')
am_params = {'m':a_m, 'L':a_L, 'r':a_r, 'I':a_I}
Accel_Mass = sympy_TMM.Sympy_Rigid_Mass_Element(am_params)

K_act = Symbol('K_act')
p_act1 = Symbol('p_act1')
AVS_params = {'K_act':K_act, 'p_act1':p_act1}
AVS = sympy_TMM.Sympy_AVS1_Element(AVS_params)


k_spring = Symbol('k_spring')
c_spring = Symbol('c_spring')
TSD_params = {'k':k_spring, 'c':c_spring}
TSD = sympy_TMM.Sympy_TSD_Element(TSD_params)

H = Symbol('H')
AVS_TFB_params = {'K_act':K_act, 'p_act1':p_act1, \
                  'H':H, 'k':k_spring, 'c':c_spring}
AVS_ThetaFB = sympy_TMM.Sympy_AVS_ThetaFB_Element(AVS_TFB_params)
GthNum = Symbol('GthNum')
GthDen = Symbol('GthDen')
AVS_ThetaFB2 = sympy_TMM.Sympy_AVS_ThetaFB_Element(AVS_TFB_params, \
                                                   Gth=GthNum/GthDen)

k_clamp = Symbol('k_clamp')
c_clamp = Symbol('c_clamp')
TSD_clamp_params = {'k':k_clamp, 'c':c_clamp}
TSD_clamp = sympy_TMM.Sympy_TSD_Element(TSD_clamp_params)

b_m = Symbol('b_m')
b_L = Symbol('b_L')
b_r = Symbol('b_r')
b_I = Symbol('b_I')
bm_params = {'m':b_m, 'L':b_L, 'r':b_r, 'I':b_I}
Base_Mass = sympy_TMM.Sympy_Rigid_Mass_Element(bm_params)

a_gain = Symbol('a_gain')

s = Symbol('s')


import SFLR_TMM

freqlim = [0.5, 30.0]
maglims = [(-40,15),(-20,10),(-20,45)]
phaselims = [(-220,0),(-400,120),(-200,220)]
NF = 3
figdir = 'figs'

maxima_bv = ['submat:submatrix(1,2,5,Usys,1,2,5)$', \
             'subcol:submatrix(1,2,5,Usys,1,2,3,4)$', \
             'nzbv:invert(submat).(-1*subcol)$', \
             'bv:zeromatrix(5,1)$', \
             'bv[5,1]:1$', \
             'bv[3,1]:nzbv[1,1]$', \
             'bv[4,1]:nzbv[2,1]$']


class JVC_model(object):
    def build_Maxima_Umat_str(self, maxind):
        """Create a Maxima line to multiply matrices together from 0
        to maxind."""
        U_str = None
        Uname = 'Ubode%i' % maxind
        for ind in range(maxind+1):
            curU = 'U%i' % ind
            if U_str is None:
                U_str = curU
            else:
                U_str = curU + '.' + U_str
        U_str = Uname + ':' + U_str + '$'
        return U_str


    def Maxima_bodes(self, sensor_inds=[], dofs=[], sensor_post=[]):
        """Find the final bode outputs for each sensor location,
        assuming Usys and the basevector bv have already been
        calculated at this point in the Maxima script."""
        bode_inds = range(len(sensor_inds))
        outlines = None
        for bi, ind, dof, post in zip(bode_inds, sensor_inds, \
                                      dofs, sensor_post):
            Uline = self.build_Maxima_Umat_str(ind)
            zline = 'z_bode%i:Ubode%i.bv$' % (bi, ind)
            bodeline1 = 'bode%i:z_bode%i[%i,1]$' % (bi, bi, dof+1)#z_bodei
            #is a 2D column matrix.  We want to extract a scalar, so
            #we need to specify column 1
            curlist = [Uline, zline, bodeline1]
            if outlines is None:
                outlines = curlist
            else:
                outlines.extend(curlist)
            if post:
                bodeline2 = 'bode%i:bode%i*%s$' % (bi, bi, post)
                outlines.append(bodeline2)
            bodeNline = 'bode%iN:ratnumer(bode%i)$' % (bi, bi)
            outlines.append(bodeNline)
            bodeDline = 'bode%iD:ratdenom(bode%i)$' % (bi, bi)
            outlines.append(bodeDline)
        return outlines
            
    
    def Usys_Maxima(self, attrlist):
        Usys_str = None
        for attr in attrlist:
            if Usys_str is None:
                Usys_str = attr
            else:
                Usys_str = attr + '.' + Usys_str
        Usys_str = 'Usys:' + Usys_str + '$'
        return Usys_str

    def save_Maxima_bodes_to_Fortran(self,  num_bodes, base_mod_name):
        outlines = None
        for bi in range(num_bodes):
            filename = base_mod_name + str(bi) + '.f'
            Nfilename = base_mod_name + str(bi) + '_num.f'
            Dfilename = base_mod_name + str(bi) + '_den.f'
            outline = 'with_stdout ("%s", fortran_optimize (bode%i))$' % \
                      (filename, bi)
            if outlines is None:
                outlines = [outline]
            else:
                outlines.append(outline)
            Nline = 'with_stdout ("%s", fortran_optimize (bode%iN))$' % \
                    (Nfilename, bi)
            outlines.append(Nline)
            Dline = 'with_stdout ("%s", fortran_optimize (bode%iD))$' % \
                    (Dfilename, bi)
            outlines.append(Dline)
        return outlines

        
    def to_Maxima(self, pathout, attrlist=['U0'], num_bodes=2, \
                  base_mod_name='maxima_bode',  **kwargs):
        """Create a Maxima batch file for the system.  attrlist is a
        list of strings referring to the element matrices.  The
        elements of attrlist should be in order, starting with U0 and
        stopping at U_n."""
        #consider adding ratdenom and ratnumer
        mylist = txt_mixin.txt_list()
        out = mylist.append
        out('showtime:all$')
        out('nolabels:true$')
        #ratvars(mubz,EIbz,Lbz,abz,betabz,c1bz,c2bz,c3bz,c4bz,s)$
        out('grind:true$')
        for attr in attrlist:
            U = getattr(self, attr)
            Uline = sympy_utils.matrix_to_Maxima_string(U,attr)
            out(Uline)
        Usys_line = self.Usys_Maxima(attrlist)
        out(Usys_line)
        mylist.extend(maxima_bv)
        bode_lines = self.Maxima_bodes(**kwargs)
        mylist.extend(bode_lines)
        save_lines = self.save_Maxima_bodes_to_Fortran(num_bodes, \
                                                       base_mod_name)
        mylist.extend(save_lines)
        self.maxima_list = mylist
        txt_mixin.dump(pathout, mylist)

        

        
    def _def_params(self):
        raise NotImplementedError

    def _load_params(self):
        myparams = self._def_params()
        self.params = myparams
    
    def _load_mod(self):
        try:
            self.mod = my_import(self.mod_name)
            reload(self.mod)
            self.func = self.mod.Bodes
        except ImportError:
            self.mod = None
            self.func = None

    def _set_a_lims(self):
        """These limits are for phase massaging of the accel Bodes"""
        self.flow = 3.44
        self.fhigh = 21.3
        self.lim2 = -10.0
        
    def __init__(self, mod_name, con_dict={}, f=exp_data.f, **kwargs):
        self.mod_name = mod_name
        self.con_dict = con_dict
        self.f = f
        for key, val in kwargs.iteritems():
            setattr(self, key, val)
        self._load_mod()
        self.plist = self.con_dict.keys()
        self._load_params()
        self._set_a_lims()

    def plot_exp(self, startfi=1, lt='y-'):
        exp_data.plot_exp(startfi=startfi, lt=lt)

    def _set_lims(self, startfi=1):
        if hasattr(self, 'maglims'):
            mymaglims = self.maglims
        else:
            mymaglims = maglims
        fis = range(startfi, startfi+NF)
        for fi, maglim, phaselim in zip(fis, mymaglims, phaselims):
            PU.SetFreqLim(fi, freqlim)
            PU.SetMagLim(fi, maglim)
            PU.SetPhaseLim(fi, phaselim)

    def save_figs(self, basename, startfi=1, del_eps=True):
        endings = ['th_v_bode','a_v_bode','a_th_bode']
        filenames = [basename +'_'+item for item in endings]
        fis = range(startfi, startfi+NF)
        for fi, fname in zip(fis, filenames):
            epsname = fname+'.eps'
            outpath = os.path.join(figdir, epsname)
            PU.mysave(outpath, fi=fi)
            if del_eps:
                os.remove(outpath)#delete eps but leave pdf

    def set_legends(self, leg_list, startfi=1, locs=None, axes=None):
        if locs is None:
            locs = [1,3,3]
        if axes is None:
            axes = [0,1,1]
        fis = range(startfi, startfi+NF)
        for fi, loc, ax in zip(fis, locs, axes):
            PU.SetLegend(fi, leg_list, loc=loc, axis=ax)

    def calc_and_massage_a_th_Bode(self,  **kwargs):
        a_th_bode = self.a_v_bode/self.th_bode
        a_th_bode.seedfreq = 10.0
        a_th_bode.seedphase = 0.0
        a_th_bode.PhaseMassage(self.f)
        #a_th_bode = SFLR_TMM.a_massage(a_th_bode, f, **kwargs)
        self.a_th_bode = a_th_bode
        return a_th_bode

    def calc_and_massage_Bodes(self, params, **kwargs):
        s = 2.0j*pi*self.f
        th_comp, a_comp = self.func(s, params)
        a_v_bode = rwkbode.rwkbode(compin=a_comp)
        th_bode = rwkbode.rwkbode(compin=th_comp)

        a_v_bode.seedfreq = 2.0
        a_v_bode.seedphase = 30.0
        a_v_bode.PhaseMassage(self.f)
        #a_v_bode = SFLR_TMM.a_massage(a_v_bode, f, **kwargs)
        th_bode.seedfreq = 1.2
        th_bode.seedphase = -90.0
        th_bode.PhaseMassage(self.f)

        self.th_bode = th_bode
        self.a_v_bode = a_v_bode
        return th_bode, a_v_bode


    def massage_a_v_bode(self):
        pass

    def massage_a_th_bode(self):
        pass
    
    def plot_Bodes(self, startfi=1, clear=False, lt='k-.'):
        kwargs = {'clear':clear, 'linetype':lt}
        rwkbode.GenBodePlot(startfi, self.f, self.th_bode, **kwargs)
        rwkbode.GenBodePlot(startfi+1, self.f, self.a_v_bode, **kwargs)
        rwkbode.GenBodePlot(startfi+2, self.f, self.a_th_bode, **kwargs)


    def calc_and_plot_Bodes(self, params=None, startfi=1, lt='k-.', \
                            clear=False, **kwargs):
        if params is None:
            params = self.params
        th_bode, a_v_bode = self.calc_and_massage_Bodes(params, **kwargs)
        a_th_bode = self.calc_and_massage_a_th_Bode(**kwargs)
        self.massage_a_v_bode()
        self.massage_a_th_bode()
        self.plot_Bodes(startfi=startfi, clear=clear, lt=lt)
        self._set_lims(startfi)
        return th_bode, a_v_bode, a_th_bode

##     def _calc_and_plot_bodes(self, params=None, startfi=1):
##         self.calc_and_plot_Bodes(params, func=self.func, \
##                                  startfi=startfi, flow=self.flow, \
##                                  fhigh=self.fhigh, lim2=self.lim2)


    def plot_start(self, startfi=1, plot_exp=True, **kwargs):
        """Overlay the experimental Bode plot with the results from
        self.params.  These are not the optimized params.
        self._load_params must be overridden for this to work."""
        if plot_exp:
            self.plot_exp(startfi)
        self.calc_and_plot_Bodes(startfi=startfi, **kwargs)

    def _save_params(self, params=None, pkl_name=None):
        if pkl_name is None:
            pkl_name = self.pkl_name
        if params is None:
            params = self.opt_params
        SFLR_TMM.save_params(params, pkl_name)

    def optimize(self, save=1, plot=1, startfi=1):
        opt_params = myoptimize(self.params, param_names=self.plist, \
                                func=self.func, con_dict=self.con_dict)
        self.opt_params = opt_params
        if save:
            self._save_params(opt_params)
        if plot:
            self.calc_and_plot_Bodes(params=opt_params, \
                                      startfi=startfi)

    def cse_to_file(self):
        sympy_TMM.cse_to_file(self.sym_bodes, self.mod_name+'.py',\
                              ['th_out','a_out'],'Bodes',\
                              inputs=['s','params'], \
                              headerfile='header.py')


    def build_param_mat(self, keys):
        """Build the matrix used for brute force parameter variation.
        There is one row per parameter and the row contains all the
        variations of that parameter that will be tried.

        keys is a list of the parameters to vary."""
        self.bfkeys = keys
        percents = array([0.9, 1.0, 1.1])
        #percents = array([0.95, 1.0, 1.05])

        N1 = len(keys)
        N2 = len(percents)

        param_mat = zeros((N2, N1))
        myparams = self.params
        
        for n, key in enumerate(keys):
            nom = getattr(myparams, key)
            vals = percents*nom
            if self.con_dict.has_key(key):
                limits = self.con_dict[key]
                if vals.min() < min(limits):
                    vals[vals.argmin()] = min(limits)
                if vals.max() > max(limits):
                    vals[vals.argmax()] = max(limits)
            param_mat[:,n] = vals

        self.param_mat = param_mat
        return self.param_mat

    def build_ind_mat(self):
        """Build the matrix of indices for the brute force parameter
        variation."""
        N2, N1 = self.param_mat.shape
        nr = N2**N1
        nc = N1

        inds = zeros((nr, nc))

        for c in range(nc):
            n = N2**(N1-c-1)
            for q in range(nr/n):
                ind = q % N2
                inds[q*n:q*n+n,c] = ind
        self.bfinds = inds

    def set_params(self, row, params=None):
        if params is None:
            params = self.bfparams
        for n, key in enumerate(self.bfkeys):
            val = self.param_mat[row[n], n]
            setattr(params, key, val)
        return params


######################################
#
# Brute Force Stuff
#
######################################
    def clear_log(self):
        f = open(self.logfile,'wb')
        f.close()
    
    def str_out(self, string):
        f = open(self.logfile,'ab')
        f.write(string)
        f.close()

    def row_of_strings_out(self, row):
    ##     assert(len(row)==N), "Got a row of bad length:len(row)=" + \
    ##                          str(len(row))
        row_str = '\t'.join(row) +'\n'
        self.str_out(row_str)
        return row_str

    def row_of_floats_out(self, row, fmt='%0.8g'):
        string_list = [fmt % item for item in row]
        row_str = self.row_of_strings_out(string_list)
        return row_str

    def build_param_list(self, params):
        mylist = [getattr(params, item) for item in self.bfkeys]
        return mylist

    def save_row(self, params, e, z1, z2):
        param_list = self.build_param_list(params)
        row = param_list + [e, z1, z2]
        self.row_of_floats_out(row)
        return row

    def find_peaks(self, a_v_bode):
        dB = a_v_bode.dBmag()
        i1 = argmax(dB[0:ind1])
        i2 = argmax(dB[ind1:])+ind1
        p1 = fexp[i1]
        p2 = fexp[i2]
        return p1, p2

    def cost(self, params=None, extra=0):
        if params is None:
            params = self.params
        return _cost(params, extra=extra, func=self.func, \
                     con_dict=self.con_dict)
    
    def run_brute_force(self, keys):
        """You must set the parameter self.logfile before calling this
        method.  This should probably happen in the __init__ method of
        a derived class."""
        self.bfparams = copy.copy(self.params)
        self.build_param_mat()
        self.build_ind_mat()

        self.bflabels = self.bfkeys+['e','p1','p2']

        self.clear_log()

        start_line = 'start time: ' + time.ctime()
        print(start_line)
        self.str_out(start_line+'\n')

        for row in self.bfinds:
            self.bfparams = self.set_params(row, self.bfparams)
            e, th_bode, a_v_bode = self.cost(self.bfparams, extra=1)
##             e, th_bode, a_v_bode = _cost(self.bfparams,
##                                        extra=1, \
##                                        func=self.func, \
##                                        con_dict=self.con_dict)
            p1, p2 = self.find_peaks(a_v_bode)
            self.save_row(self.bfparams, e, p1, p2)


        print('end time: ' + time.ctime())


    def load_bf_results(self, filename=None):
        if filename is None:
            filename = self.logfile
        self.bfdata = loadtxt(filename, skiprows=2)
        self.error = self.bfdata[:,-3]
        self.p1 = self.bfdata[:,-2]
        self.p2 = self.bfdata[:,-1]
        self.minind = self.error.argmin()
        self.minrow = self.bfdata[self.minind,:]
        return self.minrow

    def row_to_params(self, rowin=None, set_params=False):
        if rowin is None:
            rowin = self.minrow
        myparams = self._def_params()
        for key, val in zip(self.bfkeys, rowin[0:-3]):
            setattr(myparams, key, val)
        if set_params:
            self.params = myparams
        return myparams

    def save_minrow(self, pklname):
        bfparams = self.row_to_params()
        SFLR_TMM.save_params(bfparams, pklname)

    def save_params(self, pklname, params=None):
        if params is None:
            params = self.params
        SFLR_TMM.save_params(params, pklname)

        
class model1(JVC_model):
    def _def_params(self):
        myparams = SFLR_TMM.SFLR_params()
        myparams.p_act1 = 15*2*pi
        return myparams
        
    def __init__(self, mod_name='model1_bodes', \
                 pkl_name='model1_opt.pkl'):
        JVC_model.__init__(self, mod_name=mod_name, pkl_name=pkl_name, \
                           con_dict={'K_act':(0,10), \
                                     'p_act1':(0,1000), \
                                     'a_gain':(0,100)}, \
                           maglims=[(-40,10),(-25,50),(-20,70)])


class model1b(JVC_model):
    def _def_params(self):
        if self.start_pkl_name is not None:
            myparams = SFLR_TMM.load_params(self.start_pkl_name)
        elif os.path.exists(self.pkl_name):
            myparams = SFLR_TMM.load_params(self.pkl_name)
        else:
            myparams = SFLR_TMM.new_def_params()
            myparams.p_act1 = 15*2*pi
        return myparams
        
    def __init__(self, mod_name='model1b_bodes', \
                 pkl_name='model1b_opt.pkl', start_pkl_name=None):
        JVC_model.__init__(self, mod_name=mod_name, pkl_name=pkl_name, \
                           con_dict={'K_act':(0,10), \
                                     'p_act1':(0,1000), \
                                     'a_gain':(0,100), \
                                     'mu':(0.1, 0.16), \
                                     'EI':(0.12, 0.21)}, \
                           maglims=[(-40,10),(-25,50),(-20,70)], \
                           start_pkl_name=start_pkl_name)
        self.bfkeys = ['mu','EI','a_gain','K_act','p_act1']
        self.logfile = 'model1b_brute_force_params.txt'

    def myphase_massage(self, attr, lim1, lim2):
        bode = getattr(self, attr)
        phase = bode.phase
        inds1 = self.f > 10.0
        inds2 = phase > lim1
        inds = inds1 & inds2
        phase[inds] -= 360.0

        indsa = self.f < 5.0
        indsb = phase < lim2
        indsc = indsa & indsb
        phase[indsc] += 360.0
        bode.phase = phase
        setattr(self, attr, bode)
        

    def massage_a_v_bode(self):
        self.myphase_massage('a_v_bode', -40.0, -130.0)

    def massage_a_th_bode(self):
        self.myphase_massage('a_th_bode', 30.0, -30.0)

    def find_symbolic_bodes(self, save=1):
        U0 = AVS.Get_Aug_Mat(s)
        U1 = beam1.Get_Aug_Mat(s)
        U2 = Accel_Mass.Get_Aug_Mat(s)
        U3 = beam2.Get_Aug_Mat(s)
        Usys = U3*(U2*(U1*U0))
        z_b = sympy_TMM.find_base_vector(Usys)
        z_enc = U0*z_b
        enc_gain = 180.0/pi*1024.0/360.0
        theta = z_enc[1]*enc_gain

        z_a = U2*(U1*z_enc)
        atip = s**2*z_a[0]*a_gain
        self.sym_bodes = [theta, atip]
        self.Usys = Usys
        if save:
            self.cse_to_file()
        return self.sym_bodes

    def build_param_mat(self):
        JVC_model.build_param_mat(self, self.bfkeys)

    def run_brute_force(self):
        return JVC_model.run_brute_force(self, self.bfkeys)

class model2(model1):
    """This model includes the actuator, a torisonal spring/damper and
    the beam.  Accel mass is not included.  The beam is one piece"""
    def _def_params(self):
        if self.start_pkl_name is not None:
            myparams = SFLR_TMM.load_params(self.start_pkl_name)
        elif os.path.exists(self.pkl_name):
            myparams = SFLR_TMM.load_params(self.pkl_name)
        else:
            myparams = SFLR_TMM.new_def_params()
            myparams.k_spring = 10.0
            myparams.c_spring = 0.1
        return myparams
        
    def __init__(self, mod_name='model2_bodes', \
                 pkl_name='model2_opt.pkl', start_pkl_name=None):
        JVC_model.__init__(self, mod_name=mod_name, pkl_name=pkl_name, \
                           con_dict={'K_act':(0,10), \
                                     'p_act1':(0,1000), \
                                     'a_gain':(0,100), \
                                     'k_spring':(0.1, 1000), \
                                     'c_spring':(0,100), \
                                     'mu':(0.1, 0.16), \
                                     'EI':(0.12, 0.21)}, \
                           maglims=[(-40,10),(-25,20),(-20,70)], \
                           start_pkl_name=start_pkl_name)
                        
        self.bfkeys = ['k_spring','c_spring','mu','EI', \
                       'K_act','p_act1']
        self.logfile = 'model2_brute_force_params.txt'


    def find_symbolic_bodes(self, save=1):
        U0 = AVS.Get_Aug_Mat(s)
        U1 = TSD.Get_Aug_Mat(s)
        U2 = beam.Get_Aug_Mat(s)
        Usys = U2*(U1*U0)
        z_b = sympy_TMM.find_base_vector(Usys)
        z_enc = U1*(U0*z_b)
        enc_gain = 180.0/pi*1024.0/360.0
        theta = z_enc[1]*enc_gain

        z_a = U2*z_enc
        atip = s**2*z_a[0]*a_gain
        self.sym_bodes = [theta, atip]
        self.Usys = Usys
        if save:
            self.cse_to_file()
        return self.sym_bodes


    def build_param_mat(self):
        JVC_model.build_param_mat(self, self.bfkeys)
        ###override the k_spring row
        self.param_mat[:,0] = array([0.1, 0.2, 0.5])        
        ###override the c_spring row
        self.param_mat[:,1] = array([0, 0.1, 1])
        return self.param_mat


    def run_brute_force(self):
        return JVC_model.run_brute_force(self, self.bfkeys)
        
    
class model2b(model1b):
    """This model includes the accelerometer mass and a two piece
    beam."""
    def _def_params(self):
        if self.start_pkl_name is not None:
            return SFLR_TMM.load_params(self.start_pkl_name)
        case = 1
        if case == 1:
            myparams = SFLR_TMM.load_params('model2b_bf_min.pkl')
        elif case == 2 and os.path.exists(self.pkl_name):
            myparams = SFLR_TMM.load_params(self.pkl_name)
        else:
            myparams = SFLR_TMM.new_def_params()
            myparams.k_spring = 10.0
            myparams.c_spring = 0.1
        return myparams
        
    def __init__(self, mod_name='model2b_bodes', \
                 pkl_name='model2b_opt.pkl', start_pkl_name=None):
        JVC_model.__init__(self, mod_name=mod_name, pkl_name=pkl_name, \
                           con_dict={'K_act':(0,10), \
                                     'p_act1':(0,1000), \
                                     'a_gain':(0,100), \
                                     'k_spring':(0.1, 1000), \
                                     'c_spring':(0,100), \
                                     'mu':(0.1, 0.16), \
                                     'a_m':(0.006, 0.01), \
                                     'EI':(0.12, 0.21)}, \
                           maglims=[(-40,15),(-25,40),(-20,70)], \
                           start_pkl_name=start_pkl_name)
                        
        self.bfkeys = ['k_spring','c_spring','mu','EI','a_m',\
                       'K_act','p_act1']
        self.logfile = 'model2b_brute_force_params.txt'


    def find_symbolic_bodes(self, save=1):
        U0 = AVS.Get_Aug_Mat(s)
        U1 = TSD.Get_Aug_Mat(s)
        U2 = beam1.Get_Aug_Mat(s)
        U3 = Accel_Mass.Get_Aug_Mat(s)
        U4 = beam2.Get_Aug_Mat(s)
        Usys = U4*(U3*(U2*(U1*U0)))
        z_b = sympy_TMM.find_base_vector(Usys)
        z_enc = U1*(U0*z_b)
        enc_gain = 180.0/pi*1024.0/360.0
        theta = z_enc[1]*enc_gain

        z_a = U3*(U2*z_enc)
        atip = s**2*z_a[0]*a_gain
        self.sym_bodes = [theta, atip]
        self.Usys = Usys
        if save:
            self.cse_to_file()
        return self.sym_bodes



    def build_param_mat(self):
        JVC_model.build_param_mat(self, self.bfkeys)
        ###override the k_spring row
        #self.param_mat[:,0] = array([0.1, 1, 10.0])
        ###override the c_spring row
        #self.param_mat[:,1] = array([0, 0.1, 1])
        return self.param_mat


    def run_brute_force(self):
        return JVC_model.run_brute_force(self, self.bfkeys)



class model2bd(model2b):
    """This model includes the accelerometer mass and a two piece
    beam and complex stiffness for the beam."""
    def _def_params(self):
        if self.start_pkl_name is not None:
            print('loading '+self.start_pkl_name)
            myparams = SFLR_TMM.load_params(self.start_pkl_name)
        elif  os.path.exists(self.pkl_name):
            print('loading '+self.pkl_name)
            myparams = SFLR_TMM.load_params(self.pkl_name)
        else:
            myparams = SFLR_TMM.new_def_params()
            myparams.k_spring = 10.0
            myparams.c_spring = 0.1
        return myparams
        
    def __init__(self, mod_name='model2b_bodes', \
                 pkl_name='model2bd_opt.pkl', start_pkl_name=None):
        JVC_model.__init__(self, mod_name=mod_name, pkl_name=pkl_name, \
                           con_dict={'K_act':(0,10), \
                                     'p_act1':(0,1000), \
                                     'a_gain':(0,100), \
                                     'k_spring':(0.1, 1000), \
                                     'c_spring':(0,100), \
                                     'mu':(0.1, 0.16), \
                                     'a_m':(0.006, 0.01), \
                                     'EI':(0.12, 0.21), \
                                     'c_beam':(0,1.0)}, \
                           maglims=[(-40,15),(-25,20),(-20,70)], \
                           start_pkl_name=start_pkl_name)
                        
        self.bfkeys = ['k_spring','c_spring','c_beam','mu','EI','a_m',\
                       'K_act','p_act1']
        self.logfile = 'model2bd_brute_force_params.txt'


    def build_param_mat(self):
        JVC_model.build_param_mat(self, self.bfkeys)
        ###override the k_spring row
        self.param_mat[:,0] = array([0.1, 0.2, 0.5])
        ###override the c_spring row
        self.param_mat[:,1] = array([0.05, 0.12, 0.2])
        ###override the c_beam row
        self.param_mat[:,2] = array([1e-5, 1e-4, 1e-3])
        return self.param_mat

##     def build_param_mat(self):
##         """Build the matrix used for brute force parameter variation.
##         There is one row per parameter and the row contains all the
##         variations of that parameter that will be tried.

##         keys is a list of the parameters to vary."""
##         keys = self.bfkeys
##         N1 = len(keys)
##         N2 = 5

##         param_mat = zeros((N2, N1))
##         myparams = self.params
##         param_mat[:,0] = array([0.1, 0.15, 0.2, 0.5, 1.0])
##         param_mat[:,1] = array([0.05, 0.1, 0.128, 0.15, 0.2])
##         param_mat[:,2] = array([1.0e-6, 0.00001, 0.0001, 0.001, 0.01])

##         self.param_mat = param_mat




class model2cs(model2bd):
    def __init__(self, mod_name='model2cs_bodes', \
                 pkl_name='model2cs_opt.pkl', start_pkl_name=None, \
                 **kwargs):
        model2bd.__init__(self, mod_name=mod_name, pkl_name=pkl_name, \
                          start_pkl_name=start_pkl_name, \
                          **kwargs)
        self.maglims = [(-40,15),(-25,15),(-20,45)]
        self.bfkeys = ['k_spring','c_spring','c_beam']#,'mu','EI','a_m',\
                       #'K_act','p_act1']
        self.logfile = 'model2cs_brute_force_params.txt'

    
class model_w_bm(model2b):
    """This model includes the base mass as well as the accelerometer
    mass and a two piece beam."""
    def _def_params(self):
        if self.start_pkl_name is not None:
            myparams = SFLR_TMM.load_params(self.start_pkl_name)
        elif os.path.exists(self.pkl_name):
            myparams = SFLR_TMM.load_params(self.pkl_name)
        else:
            myparams = SFLR_TMM.new_def_params()
            myparams.k_spring = 10.0
            myparams.c_spring = 0.1
        return myparams

    def AVS_mat(self, s):
        U0 = AVS.Get_Aug_Mat(s)
        self.U0 = U0
        return U0


    def find_sym_matrices(self):
        self.U0 = AVS.Get_Aug_Mat(s)
        self.U1 = TSD.Get_Aug_Mat(s)
        self.U2 = Base_Mass.Get_Aug_Mat(s)
        self.U3 = TSD_clamp.Get_Aug_Mat(s)
        self.U4 = beam1.Get_Aug_Mat(s)
        self.U5 = Accel_Mass.Get_Aug_Mat(s)
        self.U6 = beam2.Get_Aug_Mat(s)


    def to_Maxima(self, pathout, num_bodes=2, \
                  base_mod_name='maxima_bode'):
        attrlist = ['U0','U1','U2','U3','U4','U5','U6']
        JVC_model.to_Maxima(self, pathout, attrlist, \
                            num_bodes=num_bodes, \
                            base_mod_name=base_mod_name)

        
    def find_symbolic_bodes(self, save=1):
        self.Usys = self.U6*(self.U5*(self.U4*(self.U3*\
                                               (self.U2*(self.U1*self.U0)))))
        z_b = sympy_TMM.find_base_vector(self.Usys)
        z_enc = self.U2*(self.U1*(self.U0*z_b))
        enc_gain = 180.0/pi*1024.0/360.0
        theta = z_enc[1]*enc_gain

        z_a = self.U5*(self.U4*(self.U3*z_enc))
        atip = s**2*z_a[0]*a_gain
        self.sym_bodes = [theta, atip]
        self.Usys = self.Usys
        if save:
            self.cse_to_file()
        return self.sym_bodes


    def Maxima_bodes(self):
        sensor_inds = [2,5]
        dofs = [1,0]
        #enc_gain = 180.0/pi*1024.0/360.0
        #enc_gain_str = str(enc_gain)
        sensor_post = ['enc_gain', 's^2*a_gain']
        return JVC_model.Maxima_bodes(self, sensor_inds, dofs, sensor_post)

        
    def __init__(self, mod_name='model_w_bm_bodes', \
                 pkl_name='model_w_bm_opt.pkl', start_pkl_name=None):
        JVC_model.__init__(self, mod_name=mod_name, pkl_name=pkl_name, \
                           con_dict={'K_act':(0,10), \
                                     'p_act1':(0,1000), \
                                     'a_gain':(0,100), \
                                     'k_spring':(0.1, 1000), \
                                     'c_spring':(0,100), \
                                     'mu':(0.1, 0.16), \
                                     'a_m':(0.006, 0.01), \
                                     'EI':(0.12, 0.21)},
                                     #'c_beam':(0,1.0)}, \
                           maglims=[(-40,15),(-25,15),(-20,45)], \
                           start_pkl_name=start_pkl_name)
                        
        self.bfkeys = ['k_spring','c_spring','mu','EI','a_m',\
                       'K_act','p_act1']
        self.logfile = 'model_w_bm_brute_force_params.txt'


    def build_param_mat(self):
        JVC_model.build_param_mat(self, self.bfkeys)
        ###override the k_spring row
        self.param_mat[:,0] = array([0.1, 0.2, 0.5])
        ###override the c_spring row
        self.param_mat[:,1] = array([0.05, 0.12, 0.2])
        ###override the c_beam row
        ###self.param_mat[:,2] = array([0.01, 0.05, 0.1])
        return self.param_mat


    def Unc(self):
        U0 = AVS.Get_Aug_Mat(s)
        U1 = TSD.Get_Aug_Mat(s)
        Unc = U1*U0
        self.Unc = Unc
        return Unc


class model_w_bm_theta_FB(model_w_bm):
    def __init__(self, mod_name='model_w_bm_bodes_ThetaFB', \
                 pkl_name='model_w_bm_opt.pkl', start_pkl_name=None):
        JVC_model.__init__(self, mod_name=mod_name, pkl_name=pkl_name, \
                           con_dict={'K_act':(0,10), \
                                     'p_act1':(0,1000), \
                                     'a_gain':(0,100), \
                                     'k_spring':(0.1, 1000), \
                                     'c_spring':(0,100), \
                                     'mu':(0.1, 0.16), \
                                     'a_m':(0.006, 0.01), \
                                     'EI':(0.12, 0.21)},
                                     #'c_beam':(0,1.0)}, \
                           maglims=[(-40,15),(-25,15),(-20,45)], \
                           start_pkl_name=start_pkl_name)

        self.bfkeys = ['k_spring','c_spring','mu','EI','a_m',\
                       'K_act','p_act1']
        self.logfile = 'model_w_bm_brute_force_params.txt'


    def find_symbolic_bodes(self, save=1):
        ## self.list = [self.avs, self.base_mass, \
        ##              self.clamp_spring, \
        ##              self.beam, self.accel, self.beam2]
        U0 = AVS_ThetaFB.Get_Aug_Mat(s)
        U1 = Base_Mass.Get_Aug_Mat(s)
        U2 = TSD_clamp.Get_Aug_Mat(s)
        U3 = beam1.Get_Aug_Mat(s)
        U4 = Accel_Mass.Get_Aug_Mat(s)
        U5 = beam2.Get_Aug_Mat(s)
        Usys = U5*(U4*(U3*(U2*(U1*U0))))
        z_b = sympy_TMM.find_base_vector(Usys)
        z_enc = U1*(U0*z_b)
        enc_gain = 180.0/pi*1024.0/360.0
        theta = z_enc[1]*enc_gain

        z_a = U4*(U3*(U2*z_enc))
        atip = s**2*z_a[0]*a_gain
        self.sym_bodes = [theta, atip]
        self.Usys = Usys
        if save:
            self.cse_to_file()
        return self.sym_bodes
    

    def Unc(self):
        U0 = AVS_ThetaFB.Get_Aug_Mat(s)
        U1 = Base_Mass.Get_Aug_Mat(s)
        U2 = TSD_clamp.Get_Aug_Mat(s)
        U3 = beam1.Get_Aug_Mat(s)
        U4 = Accel_Mass.Get_Aug_Mat(s)
        Unc = U4*(U3*(U2*(U1*U0)))
        self.Unc = Unc
        return Unc

    
    def Unc_row_to_file(self, pathout):
        Unc = self.Unc
        Unc00 = Unc[0,0]
        Unc01 = Unc[0,1]
        Unc02 = Unc[0,2]
        Unc03 = Unc[0,3]
        Unc04 = Unc[0,4]
        outlabels = ['Unc00', 'Unc01', 'Unc02', 'Unc03', 'Unc04']
        funcname = 'Unc_row0'
        expr_list = [Unc00, Unc01, Unc02, Unc03, Unc04]
        inputs = ['s', 'params', 'Gth']
        replace_dict = {'Gth':'Gth(s)'}
        sympy_TMM.cse_to_file(expr_list, \
                              pathout, \
                              outlabels, \
                              funcname, \
                              inputs, \
                              headerfile='header.py', \
                              replace_dict=replace_dict)
        


class model_w_bm_theta_FB_ND(model_w_bm_theta_FB):
    """This class is setup to find the numerator and denominator of
    the theta feedback transfer functions.  As such, it is important
    that Gth is made up of GthNum/GthDen."""
    ## def find_symbolic_bodes(self, save=1):
    ##     ## self.list = [self.avs, self.base_mass, \
    ##     ##              self.clamp_spring, \
    ##     ##              self.beam, self.accel, self.beam2]
    ##     U0 = AVS_ThetaFB2.Get_Aug_Mat(s)
    ##     U1 = Base_Mass.Get_Aug_Mat(s)
    ##     U2 = TSD_clamp.Get_Aug_Mat(s)
    ##     U3 = beam1.Get_Aug_Mat(s)
    ##     U4 = Accel_Mass.Get_Aug_Mat(s)
    ##     U5 = beam2.Get_Aug_Mat(s)
    ##     Usys = U5*(U4*(U3*(U2*(U1*U0))))
    ##     z_b = sympy_TMM.find_base_vector(Usys)
    ##     z_enc = U1*(U0*z_b)
    ##     enc_gain = 180.0/pi*1024.0/360.0
    ##     theta = z_enc[1]*enc_gain

    ##     z_a = U4*(U3*(U2*z_enc))
    ##     atip = s**2*z_a[0]*a_gain
    ##     self.sym_bodes = [theta, atip]
    ##     self.Usys = Usys
    ##     if save:
    ##         self.cse_to_file()
    ##     return self.sym_bodes
    


    def find_sym_matrices(self):
        self.U0 = AVS_ThetaFB2.Get_Aug_Mat(s)
        self.U1 = Base_Mass.Get_Aug_Mat(s)
        self.U2 = TSD_clamp.Get_Aug_Mat(s)
        self.U3 = beam1.Get_Aug_Mat(s)
        self.U4 = Accel_Mass.Get_Aug_Mat(s)
        self.U5 = beam2.Get_Aug_Mat(s)


    def to_Maxima(self, pathout, num_bodes=2, \
                  base_mod_name='maxima_bode'):
        attrlist = ['U0','U1','U2','U3','U4','U5']
        JVC_model.to_Maxima(self, pathout, attrlist, \
                            num_bodes=num_bodes, \
                            base_mod_name=base_mod_name)


    def find_symbolic_bodes(self, save=1):
        self.Usys = self.U5*(self.U4*(self.U3*\
                                (self.U2*(self.U1*self.U0))))
        z_b = sympy_TMM.find_base_vector(self.Usys)
        z_enc = self.U1*(self.U0*z_b)
        enc_gain = 180.0/pi*1024.0/360.0
        theta = z_enc[1]*enc_gain

        z_a = self.U4*(self.U3*(self.U2*z_enc))
        atip = s**2*z_a[0]*a_gain
        self.sym_bodes = [theta, atip]
        self.Usys = self.Usys
        if save:
            self.cse_to_file()
        return self.sym_bodes


    def Maxima_bodes(self):
        sensor_inds = [1,4]
        dofs = [1,0]
        #enc_gain = 180.0/pi*1024.0/360.0
        #enc_gain_str = str(enc_gain)
        sensor_post = ['enc_gain', 's^2*a_gain']
        return JVC_model.Maxima_bodes(self, sensor_inds, dofs, sensor_post)
