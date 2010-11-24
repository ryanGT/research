from pylab import *
from scipy import *
from scipy import optimize

import copy
import pylab as PL
import pylab_util as PU

from IPython.Debugger import Pdb

import os, sys

import controls, rwkos, txt_mixin

import bode_plot_overlayer as BPO
import SFLR_bode_options

from rwkmisc import rowwise, my_import, LoadPickle

msg1 = """Cannot calculate bode response of a
Ga_ThetaFB_System unless self.thfb_a_comp
is defined."""

msg2 = """Cannot calculate the step response of
a Ga_ThetaFB_System unless self.ROM
is defined."""

def comp_to_db(comp_mat):
    dB_mat = 20*log10(abs(comp_mat))
    return dB_mat


def comp_to_zeros_db(comp_mat):
    zeros_mat = 1.0/comp_mat
    return comp_to_db(zeros_mat)


def pole_locs(dB_mat):
    bool1 = dB_mat[1:-1,1:-1] > dB_mat[2:,1:-1]
    bool2 = dB_mat[1:-1,1:-1] > dB_mat[0:-2,1:-1]
    bool3 = dB_mat[1:-1,1:-1] > dB_mat[1:-1,2:]
    bool4 = dB_mat[1:-1,1:-1] > dB_mat[1:-1,0:-2]
    pole_bool = bool1 & bool2 & bool3 & bool4
    indr, indc = where(pole_bool)
    return indr, indc



class ga_theta_fb_system(object):
    def __init__(self, Ga, thfb_a_comp=None, \
                 ROM_model=None, substr=None, \
                 accel_comp_bode_opts=None, \
                 afb_bode_opts=None, \
                 tfb_contour_sys=None, \
                 levels=None):
        self.Ga = Ga
        self.ROM = ROM_model
        self.thfb_a_comp = thfb_a_comp
        self.substr = substr
        if accel_comp_bode_opts is None:
            accel_comp_bode_opts = SFLR_bode_options.Accel_Comp_Bode_opts[1]
        self.accel_comp_bode_opts = accel_comp_bode_opts
        if afb_bode_opts is None:
            afb_bode_opts = SFLR_bode_options.AccelFB_Bode_opts[1]
        self.afb_bode_opts = afb_bode_opts
        self.tfb_contour_sys = tfb_contour_sys
        if levels is None:
            levels = arange(-10, 50, 3)
        self.levels = levels
        if tfb_contour_sys is not None:
            self._map_from_tfb_contour_sys()
        if self.ROM is not None:
            self.update_rom(Ga)


    def substr_to_filename(self):
        filename = self.substr
        replist = [' ','.','=']
        for item in replist:
            filename = filename.replace(item, '_')
        filename = filename.replace('__','_')
        filename = filename.replace('__','_')
        self.subfilename = filename
        return filename
        

    def build_contour_filename(self, pre='', post='', ext='.eps'):
        self.substr_to_filename()
        self.contour_name = 'contour_plot_' + pre + self.subfilename + post + ext
        return self.contour_name


    def build_contour_path(self, folder='', pre='', post='', ext='.eps'):
        self.build_contour_filename(pre=pre, post=post)
        self.contour_path = os.path.join(folder, self.contour_name)
        return self.contour_path
    

    def _map_from_tfb_contour_sys(self):
        self.comp_mat = self.tfb_contour_sys.comp_mat
        self.ol_db_mat = comp_to_db(self.comp_mat)
        self.s_contour = self.tfb_contour_sys.s
        self.f_contour = self.tfb_contour_sys.f
        self.im_contour = self.tfb_contour_sys.im
        self.a_theta_comp = self.tfb_contour_sys.a_theta_comp 


    def plot_db_mat_contour(self, dB_afb, fi, title=None, \
                            myxlim=[-20,2],myylim=[-2,20],\
                            zoomin=False):
        figure(fi)
        clf()
        contour(-self.f_contour*2*pi, self.im_contour*2*pi , dB_afb, self.levels)
        PL.xlabel('real($s$)')
        PL.ylabel('imag($s$)')
        if hasattr(self, 'contour_fignums'):
            self.contour_fignums.append(fi)
        else:
            self.contour_fignums = [fi]
        if title is not None:
            PL.title(title)
        if zoomin:
            PL.xlim(myxlim)
            PL.ylim(myylim)


    def find_poles(self, dB_mat):
        indr, indc = pole_locs(dB_mat)

        pole_vect = zeros((len(indr)), dtype='D')
        nc = len(self.f_contour)
        nr = len(self.im_contour)

        for i in range(len(indr)):
            ir = indr[i]
            ic = indc[i]
            if ir < nr - 1:
                ir += 1
            if ic < nc - 1:
                ic += 1
            pole_i = -self.f_contour[ic]*2*pi + self.im_contour[ir]*2.0j*pi
            pole_vect[i] = pole_i

        return pole_vect


    def set_contour_lims(self, xlim=None, ylim=None):
        """Set the limits of all figures in self.contour_fignums.  If
        xlim or ylim is None, use the max and min values of
        self.s_contour.real and self.s_contour.imag."""
        if xlim is None:
            xlim = [self.s_contour.real.min(), self.s_contour.real.max()]
        if ylim is None:
            ylim = [self.s_contour.imag.min(), self.s_contour.imag.max()]
        for fi in self.contour_fignums:
            PU.SetXlim(fi, xlim)
            PU.SetYlim(fi, ylim)
        


    def calc_ga_comp(self, s):
        Ga_comp = self.Ga(s)
        self.Ga_comp = Ga_comp


    def calc_ga_comp_from_f(self, f):
        s = 2.0j*pi*f
        self.calc_ga_comp(s)

        
    def update_rom(self, Ga=None):
        assert self.ROM is not None, msg2
        if Ga is not None:
            self.ROM.Ga = Ga
        self.ROM.build_TFs()
        self.theta_TF_simp = self.ROM.theta_TF.simplify()
        self.accel_TF_simp = self.ROM.accel_TF.simplify()


    def lsim(self, u, t, calca=True):
        assert self.ROM is not None, msg2
        self.t = t
        self.u = u
        self.theta_t = self.theta_TF_simp.lsim(u,t)
        if calca:
            self.accel_t = self.accel_TF_simp.lsim(u,t)
        else:
            self.accel_t = zeros_like(theta_t)
        return self.theta_t, self.accel_t


    def plot_time_domain(self, fig=None, fi=1, \
                         plotu=False, clear=False):
        if fig is None:
            fig = figure(1)
        if clear:
            fig.clf()
            
        if fig.get_axes():
            ax = fig.get_axes()[0]
        else:
            ax = fig.add_subplot(1,1,1)
        if plotu:
            ax.plot(self.t, self.u, label='$u$')
        plot(self.t, self.theta_t, label='$\\theta_{%s}$' % self.substr)
        plot(self.t, self.accel_t, label='$\\ddot{x}_{%s}$' % self.substr)        
        
        
    def calc_accel_comp_bode(self, f):
        assert self.thfb_a_comp is not None, msg1
        if not hasattr(self, 'Ga_comp'):
            self.calc_ga_comp_from_f(f)
        a_comp = self.Ga_comp*self.thfb_a_comp
        a_bode = BPO.comp_to_Bode(a_comp, f, \
                                  bode_opt=self.accel_comp_bode_opts)
        self.Accel_comp_Bode = a_bode


    def _set_kwargs(self, **kwargs):
        if not kwargs.has_key('label'):
            kwargs['label'] = self.substr
        return kwargs


    def plot_accel_comp_bode(self, f, fi=1, **kwargs):
        if not hasattr(self, 'Accel_comp_Bode'):
            self.calc_accel_comp_bode(f)
        kwargs = self._set_kwargs(**kwargs)
        BPO._plot_bode(self.Accel_comp_Bode, \
                       self.accel_comp_bode_opts, \
                       f, fignum=fi, **kwargs)
        

    def calc_afb_bode(self, f):
        assert self.thfb_a_comp is not None, msg1
        if not hasattr(self, 'Ga_comp'):
            self.calc_ga_comp_from_f(f)
        afb_a_comp = self.thfb_a_comp/(1.0+self.Ga_comp*self.thfb_a_comp)
        afb_a_bode = BPO.comp_to_Bode(afb_a_comp, f, \
                                      bode_opt=self.afb_bode_opts)
        self.Accel_AFB_Bode = afb_a_bode
        self.afb_a_comp = afb_a_comp


    def plot_afb_bode(self, f, fi=1, **kwargs):
        if not hasattr(self, 'Accel_AFB_Bode'):
            self.calc_afb_bode(f)
        kwargs = self._set_kwargs(**kwargs)    
        BPO._plot_bode(self.Accel_AFB_Bode, \
                       self.afb_bode_opts, \
                       f, fignum=fi, **kwargs)


    def nyquist_plot(self, f, fi=1, clear=False, **kwargs):
        if not hasattr(self, 'Ga_comp'):
            self.calc_ga_comp_from_f(f)
        Nyq = self.Ga_comp*self.thfb_a_comp
        figure(fi)
        if clear:
            clf()
        mirror = Nyq.conj()
        plot(Nyq.real, Nyq.imag)
        plot(mirror.real, mirror.imag)


    def calc_ga_comp_contour(self):
        s = self.s_contour
        Ga_comp = self.Ga(s)
        self.Ga_comp_contour = Ga_comp
    

    def calc_afb_compmat(self):
        self.calc_ga_comp_contour()
        self.comp_afb = self.comp_mat/(1 + self.Ga_comp_contour*self.comp_mat)
        self.comp_afb_db = comp_to_db(self.comp_afb)
        return self.comp_afb_db


    def find_afb_zeros(self):
        if not hasattr(self, 'comp_afb'):
            self.calc_afb_compmat()
        zeros_mat = 1.0/self.comp_afb
        self.afb_zeros_db = comp_to_db(zeros_mat)
        return self.afb_zeros_db


    def calc_theta_comp_contour(self):
        if not hasattr(self, 'comp_afb'):
            self.calc_afb_compmat()
        self.theta_comp = self.comp_afb/self.a_theta_comp
        return self.theta_comp


    def find_afb_theta_zeros(self):
        if not hasattr(self, 'theta_comp'):
            self.calc_theta_comp_contour()
        zeros_mat = 1.0/self.theta_comp
        self.afb_theta_zeros_db = comp_to_db(zeros_mat)
        self.afb_theta_zeros = self.find_poles(self.afb_theta_zeros_db)
        return self.afb_theta_zeros


    def plot_afb_contour(self, fi, **kwargs):
        if not hasattr(self, 'comp_afb_db'):
            self.calc_afb_compmat()
        self.plot_db_mat_contour(self.comp_afb_db, fi, **kwargs)


    def plot_ol_contour(self, fi, **kwargs):
        self.plot_db_mat_contour(self.ol_db_mat, fi, **kwargs)


    def find_afb_poles(self):
        if not hasattr(self, 'comp_afb_db'):
            self.calc_afb_compmat()
        self.afb_poles = self.find_poles(self.comp_afb_db)
        return self.afb_poles


    def find_ol_poles(self):
        self.ol_poles = self.find_poles(self.ol_db_mat)
        return self.ol_poles

        

class kp_accel_theta_fb_system(ga_theta_fb_system):
    def __init__(self, Ga, **kwargs):
        Ga = float(Ga)
        ga_theta_fb_system.__init__(self, Ga, **kwargs)


    def calc_ga_comp(self, s):
        self.Ga_comp = self.Ga


    def calc_ga_comp_contour(self):
        self.Ga_comp_contour = self.Ga


def find_zetas(pole_vect):
    z = abs(real(pole_vect))/(abs(pole_vect))
    return z


def right_of_penalty(pole, line=-10.0):
    penalty = 0.0
    if real(pole) > line:
        penalty = real(pole) - line
    return penalty


def print_pole_locs(pole_vect, fmt='%0.3f%+0.3fj'):
    for pole in pole_vect:
        print(fmt % (real(pole), imag(pole)))


def pole_locs_to_list(pole_vect, fmt='%0.3f%+0.3fj'):
    outlist = []
    for pole in pole_vect:
        curstr = fmt % (real(pole), imag(pole))
        outlist.append(curstr)
    return outlist


def Ga_to_C5(Ga):
    C5 = Ga.den.coeffs[1:].tolist() + Ga.num.coeffs.tolist()
    return C5


def Ga_to_C4(Ga):
    blist = Ga.num.coeffs.tolist()
    gain = blist[-2]
    z = blist[-1]/gain
    C4 = Ga.den.coeffs[1:].tolist() + [gain]
    return C4


def pole_zero_sorter(pole_vect, r=6.0*2*pi, imag_tol=1e-5):
    filt_poles = [pole for pole in pole_vect if imag(pole) > -1e-7]#only positive poles
    small_poles = [pole for pole in filt_poles if abs(pole) < r]
    small_imag = [pole for pole in small_poles if \
                  abs(imag(pole)) > imag_tol]
    small_real = [pole for pole in small_poles if abs(imag(pole)) <= imag_tol]
    return small_imag, small_real
    

class ga_pole_optimizer(ga_theta_fb_system):
    """This is a class for optimizing the poles of a
    ga_theta_fb_system based on contour plots and closed-loop pole
    locations."""
    def __init__(self, *args, **kwargs):
        ga_theta_fb_system.__init__(self, *args, **kwargs)
        self.logger = None
        self.optcase = 30
        self.second_imag_zeta = 0.95


    def _find_poles_if_necessary(self):
        if not hasattr(self, 'afb_poles'):
            self.find_afb_poles()
        

    def unstable_check(self, tol=1e-7):
        self._find_poles_if_necessary()
        bool_vect = real(self.afb_poles) > tol
        return bool_vect.any()


    def filter_immovable(self, pole_vect=None, \
                         zeros=False, \
                         loc=None, eps=1.0e-7):
        """There seems to be a small pole near -1.3194689145077254 that is
        present for all choices of Ga.  I conclude it is immovable and
        will ignore it in my cost functions."""
        if pole_vect is None:
            self._find_poles_if_necessary()
            pole_vect = self.afb_poles

        if loc is None:
            if zeros:
                loc = -1.38230077
            else:
                loc = -1.3194689145077254
        
        def myfunc(pole):
            rp = real(pole)
            ignore = loc
            if (ignore - eps) < rp < (ignore + eps):
                return False
            else:
                return True
        filtered = filter(myfunc, pole_vect)

        if zeros:
            self.filt_afb_theta_zeros = filtered
        else:
            self.filt_afb_poles = filtered
        return filtered


    def find_immovable_pole(self):
        self.immovable_small = [pole for pole in self.small_real
                                if pole not in self.small_real2]
        return self.immovable_small
    

    def pole_sorter(self, r=6.0*2*pi, imag_tol=1e-5):
        """Find the poles whose abs is less than r and sort them into real
        and imaginary based on whether their imaginary part is greater
        than or less than imag_tol."""
        self._find_poles_if_necessary()
        small_imag, small_real = pole_zero_sorter(self.afb_poles, \
                                                  r=r, imag_tol=imag_tol)
        self.small_real = small_real
        self.small_imag = small_imag
        self.small_real2 = self.filter_immovable(small_real)
        return small_imag, small_real


    def second_imag_penalty(self):
        penalty = 0.0
        imag_list = self.small_imag
        if len(imag_list) > 1:
            imag_vect = array(imag_list)
            #drop the right most element
            ind = real(imag_vect).argmax()
            imag_list2 = copy.copy(imag_list)
            imag_list2.pop(ind)
            imag_vect2 = array(imag_list2)
            zetas = find_zetas(imag_vect2)
            mybool = zetas < self.second_imag_zeta
            if mybool.any():
                penalty = 1.0
        return penalty


    def print_small_poles(self):
        self.pole_sorter(r=20.0)
        print('small_imag poles = ')
        print_pole_locs(self.small_imag)
        print('')
        print('small_real2 poles = ')
        print_pole_locs(self.small_real2)
        print('')
        self.find_immovable_pole()
        print('immovable_small pole(s) = ')
        print_pole_locs(self.immovable_small)
        print('')


    def find_other_poles(self):
        already_included = self.small_imag + self.small_real2 + self.immovable_small
        other_poles = [pole for pole in self.afb_poles if pole not in already_included]
        return other_poles
    

    def find_other_zeros(self):
        mylist = self.afb_theta_zeros.tolist()
        already_included = self.small_imag_zeros + self.small_real_zeros
        other_zeros = [zero for zero in  mylist if zero not in already_included]
        return other_zeros
    
        
    def pole_zero_report(self):
        self.pole_sorter(r=20.0)
        mylist = []
        out = mylist.append
        out('small imag poles = ')
        small_imag_list = pole_locs_to_list(self.small_imag)
        mylist.extend(small_imag_list)
        out('small real2 poles = ')
        small_real2_list = pole_locs_to_list(self.small_real2)
        mylist.extend(small_real2_list)
        self.find_immovable_pole()
        out('immovable_small pole(s) = ')
        immovable_list = pole_locs_to_list(self.immovable_small)
        mylist.extend(immovable_list)
        other_poles = self.find_other_poles()
        out('other poles = ')
        mylist.extend(pole_locs_to_list(other_poles))
        self.find_afb_theta_zeros()
        self.theta_zero_sorter()
        out('small imag zeros = ')
        mylist.extend(pole_locs_to_list(self.small_imag_zeros))
        out('small real zeros = ')
        mylist.extend(pole_locs_to_list(self.small_real_zeros))
        other_zeros = self.find_other_zeros()
        out('other zeros = ')
        mylist.extend(pole_locs_to_list(other_zeros))
        self.pz_report = mylist
        return mylist


    def pole_cost(self, verbosity=0):
        self.pole_sorter(r=20.0)
        small_imag = self.small_imag
        small_real2 = self.small_real2
        ## if small_real2:
        ##     cost += max(real(small_real2))
        my_small_poles = small_imag + small_real2
        cost = 0.0
        ## ind_small_imag = real(array(small_imag)).argmax()
        ## smallest_imag = small_imag[ind_small_imag]
        ## if real(smallest_imag) > -10.0:
        ##     si_penalty = 10.0*(real(smallest_imag)+10)
        ##     cost += si_penalty
        ##     if verbosity > 0:
        ##         print('smallest_imag penalty = %0.4g' % si_penalty)
        myline = -8.0#push polls to the left of this line if possible
        small_imag_penalties = [right_of_penalty(pole, line=myline) for pole \
                                in small_imag]
        small_real_penalties = [right_of_penalty(pole) for pole in  small_real2]
        imag_cost = sum(small_imag_penalties)
        real_cost = sum(small_real_penalties)
        cost += 10*imag_cost + real_cost
        second_imag_pen = self.second_imag_penalty()*50.0
        cost += second_imag_pen
        ## for pole in my_small_poles:
        ##     if real(pole) > -10.0:
        ##         cost += real(pole) + 10.0
        if self.unstable_check():
            cost += 1e3
        if verbosity > 0:
            print('small_imag poles = ' + str(small_imag))
            print('small_real2 poles = ' + str(small_real2))
            print('small_imag_penalties = %s' % small_imag_penalties)
            print('small_real_penalties = %s' % small_real_penalties)
            print('imag_cost = %0.4g' % imag_cost)
            print('real_cost = %0.4g' % real_cost)
            print('second_imag_pen = %0.4g' % second_imag_pen)

        return cost


    def theta_zero_sorter(self, r=6.0*2*pi, imag_tol=1e-5):
        """Find the poles whose abs is less than r and sort them into real
        and imaginary based on whether their imaginary part is greater
        than or less than imag_tol."""
        self.find_afb_theta_zeros()
        small_imag, small_real = pole_zero_sorter(self.afb_theta_zeros, \
                                                  r=r, imag_tol=imag_tol)
        self.small_real_zeros = small_real
        self.small_imag_zeros = small_imag
        self.small_real2_zeros = self.filter_immovable(small_real)
        return small_imag, small_real


    def zero_cost(self, verbosity=0):
        small_imag, small_real = self.theta_zero_sorter(r=10)
        small_real2 = self.filter_immovable(small_real, zeros=True)
        cost = 0.0
        if small_imag:
            cost += sum(imag(small_imag))
        if verbosity > 0:
            print('small_imag zeros = ' + str(small_imag))
            print('small_real2 zeros = ' + str(small_real2))
        return cost


    def build_Ga(self, C):
        #a2,a1,a0,b1,b0,gain = C
        if len(C) == 4:
            a2,a1,a0,gain = C
            b1 = gain
            b0 = z*gain
            gain2 = 1.0
        elif len(C) == 5:
            a2,a1,a0,b1,b0 = C
            b2 = 0.0
            gain2 = 1.0
            ## b1 = gain
            ## b0 = z*gain
        ## elif len(C) == 6:
        ##     a2,a1,a0,b2,b1,b0 = C
        elif len(C) == 6:
            a2,a1,a0,b1,b0,gain2 = C
        Ga = controls.TF([b1,b0],[1,a2,a1,a0])*gain2
        self.Ga = Ga
        return Ga
        

    def mycost(self, C, verbosity=0):
        if verbosity > 10:
            print('C = ' + str(C))
        Ga = self.build_Ga(C)
        self.calc_afb_compmat()
        self.calc_theta_comp_contour()
        self.find_afb_poles()
        pcost = self.pole_cost()
        zcost = self.zero_cost(verbosity=verbosity)
        cost = pcost + zcost*10.0
        neg_C_penalty = 0.0
        if (array(C) < 0.0).any():
            neg_C_penalty = 1.0e3
            cost += neg_C_penalty
        if verbosity > 0:
            print('pcost = %s' % pcost)
            print('zcost = %s' % zcost)
            print('cost = %s' % cost)
            #print('pole locs = %s' % pole_vect)
            #print('zero locs = %s' % myzeros)
            print('neg_C_penalty = %s' % neg_C_penalty)
            print('-'*30)
        if self.logger is not None:
            self.logger.log(C, cost)
        return cost
    

    def create_ig(self):
        return Ga_to_C5(self.Ga)


    def set_C_from_Ga(self):
        self.C_opt = Ga_to_C5(self.Ga)


    def optimize_Ga(self, C_ig=None, logger=None):
        self.logger = logger
        if C_ig is None:
            C_ig = self.create_ig()
        C_opt = optimize.fmin(self.mycost, C_ig)
        self.C_opt = C_opt


    def save_opt_case(self, optcase=None, folder=None):
        if optcase is None:
            optcase = self.optcase
        if folder is None:
            folder = rwkos.FindFullPath('~/siue/Research/PSoC_Research/SFLR_2010/vibration_suppression_design/TMM/')
        ws = ' '*4
        outname = 'Ga_opt_optcase_%i.txt' % optcase
        outpath = os.path.join(folder, outname)
        outlist = ['elif case == %i:' % optcase]
        outlist.append(ws + '# contour optimization story')
        outlist.append(ws + '# ' + str(self.__class__))
        cline1 = ws + 'C_opt = array([ %0.8f, %0.8f, %0.8f, \\' % tuple(self.C_opt[0:3])
        cline2 = ws + ' '*16 + '%0.8f, %0.8f])' % tuple(self.C_opt[3:])
        outlist.append(cline1)
        outlist.append(cline2)
        savelines = [ws + item for item in outlist]
        txt_mixin.dump(outpath, savelines)
        return savelines


    def header_list(self):
        header = ['contour optimization story', \
                  'optcase = %i' % self.optcase, \
                  str(self.__class__), \
                  ]
        cline1 = 'C_opt = array([ %0.8f, %0.8f, %0.8f, \\' % tuple(self.C_opt[0:3])
        cline2 = ' '*16 + '%0.8f, %0.8f])' % tuple(self.C_opt[3:])
        header.append(cline1)
        header.append(cline2)
        self.header = header
        return header
        

    def report_list(self):
        report = []
        report.extend(self.header_list())
        report.extend(self.pole_zero_report())
        self.report = report
        return self.report


    def print_report(self):
        report = self.report_list()
        for line in report:
            print(line)
        
            

class first_pole_optimizer(ga_pole_optimizer):
    """I am attempting to tell a story in the SFLR_2010 paper about
    how I developed a good cost function for this problem.  I am
    trying to do this by going back retracing the steps I took to get
    to my final cost function.

    What was the fisrt thing I tried?  If I was starting fresh, what
    is the first thing I would try?

    It seems like the main issue is driving the mode 1 pole to the
    left.

    This class will also serve as the base class for optimizers
    without a zero cost."""
    def __init__(self, *args, **kwargs):
        ga_pole_optimizer.__init__(self, *args, **kwargs)
        self.optcase = 31
        self.use_zcost = False
        
        
    def pole_cost(self, verbosity=0):
        self.pole_sorter(r=20.0)
        small_imag = self.small_imag
        small_real2 = self.small_real2
        cost = real(small_imag).max()
        if self.unstable_check():
            cost += 1e3
        return cost


    def mycost(self, C, verbosity=0):
        if verbosity > 10:
            print('C = ' + str(C))
        Ga = self.build_Ga(C)
        self.calc_afb_compmat()
        self.find_afb_poles()
        pcost = self.pole_cost()
        if self.use_zcost:
            self.calc_theta_comp_contour()
            zcost = self.zero_cost(verbosity=verbosity)
        else:
            zcost = 0.0
        cost = pcost + zcost*10.0
        neg_C_penalty = 0.0
        if (array(C) < 0.0).any():
            neg_C_penalty = 1.0e3
            cost += neg_C_penalty
        if verbosity > 0:
            print('pcost = %s' % pcost)
            #print('zcost = %s' % zcost)
            print('cost = %s' % cost)
            #print('pole locs = %s' % pole_vect)
            #print('zero locs = %s' % myzeros)
            print('neg_C_penalty = %s' % neg_C_penalty)
            print('-'*30)
        if self.logger is not None:
            self.logger.log(C, cost)
        return cost


class second_pole_optimizer(first_pole_optimizer):
    """drive all the small poles as far to the left as possible"""
    def __init__(self, *args, **kwargs):
        ga_pole_optimizer.__init__(self, *args, **kwargs)
        self.optcase = 32
        self.use_zcost = False


    def pole_cost(self, verbosity=0):
        self.pole_sorter(r=20.0)
        small_imag = self.small_imag
        small_real2 = self.small_real2
        cost = sum(real(small_imag)) + sum(real(small_real2))
        if self.unstable_check():
            cost += 1e3
        return cost


class third_pole_optimizer(second_pole_optimizer):
    """drive all the small poles as far to the left as possible, but
    with a smaller radius to define the small poles"""
    def __init__(self, *args, **kwargs):
        ga_pole_optimizer.__init__(self, *args, **kwargs)
        self.optcase = 33
        self.use_zcost = False


    def pole_cost(self, verbosity=0):
        self.pole_sorter(r=15.0)
        small_imag = self.small_imag
        small_real2 = self.small_real2
        cost = sum(real(small_imag)) + sum(real(small_real2))
        if self.unstable_check():
            cost += 1e3
        return cost


class fourth_pole_optimizer(third_pole_optimizer):
    """drive all small poles to the left of -10"""
    def __init__(self, *args, **kwargs):
        ga_pole_optimizer.__init__(self, *args, **kwargs)
        self.use_zcost = False
        self.optcase = 34
        self.imagline = -10
        self.realline = -10
        self.imagweight = 1.0
        self.secondweight = 0.0


    def pole_cost(self, verbosity=0):
        self.pole_sorter(r=20.0)
        small_imag = self.small_imag
        small_real2 = self.small_real2
        my_small_poles = small_imag + small_real2
        cost = 0.0
        small_imag_penalties = [right_of_penalty(pole, line=self.imagline) for pole \
                                in small_imag]
        small_real_penalties = [right_of_penalty(pole, line=self.realline) for pole \
                                in  small_real2]
        imag_cost = sum(small_imag_penalties)
        real_cost = sum(small_real_penalties)
        cost += imag_cost*self.imagweight + real_cost
        if self.secondweight > 0.0:
            second_imag_pen = self.second_imag_penalty()*self.secondweight
            cost += second_imag_pen
        if self.unstable_check():
            cost += 1e3
        return cost


class fifth_pole_optimizer(fourth_pole_optimizer):
    """drive all small poles to the left of -10.  Penalize a second
    imag pole"""
    def __init__(self, *args, **kwargs):
        ga_pole_optimizer.__init__(self, *args, **kwargs)
        self.use_zcost = False
        self.optcase = 35
        self.imagline = -10
        self.realline = -10
        self.imagweight = 1.0
        self.secondweight = 10.0


class sixth_pole_optimizer(fourth_pole_optimizer):
    """drive all small poles to the left of -10.  Penalize a second
    imag pole with weight of 10; weight image poles more"""
    def __init__(self, *args, **kwargs):
        ga_pole_optimizer.__init__(self, *args, **kwargs)
        self.use_zcost = False
        self.optcase = 36
        self.imagline = -10
        self.realline = -10
        self.imagweight = 10.0
        self.secondweight = 10.0


class seventh_pole_optimizer(fourth_pole_optimizer):
    """drive all small poles to the left of -10.  Penalize a second
    imag pole with weight of 50; weight image poles more"""
    def __init__(self, *args, **kwargs):
        ga_pole_optimizer.__init__(self, *args, **kwargs)
        self.use_zcost = False
        self.optcase = 37
        self.imagline = -10
        self.realline = -10
        self.imagweight = 10.0
        self.secondweight = 50.0


class eighth_pole_optimizer(fourth_pole_optimizer):
    """drive all small real poles to the left of -10. drive small
    image poles to the left of -8.  Penalize a second imag pole with
    weight of 50; weight image poles more"""
    def __init__(self, *args, **kwargs):
        ga_pole_optimizer.__init__(self, *args, **kwargs)
        self.use_zcost = False
        self.optcase = 38
        self.imagline = -8
        self.realline = -10
        self.imagweight = 10.0
        self.secondweight = 50.0


class pole_optimizer_with_zcost(fourth_pole_optimizer):
    def __init__(self, Ga, optcase=None, imagline=-8, \
                 realline=-10, imagweight=10, \
                 secondweight=50.0, **kwargs):
        ga_pole_optimizer.__init__(self, Ga, **kwargs)
        self.use_zcost = True
        self.optcase = optcase
        self.imagline = imagline
        self.realline = realline
        self.imagweight = imagweight
        self.secondweight = secondweight


class pole_optimizer9(pole_optimizer_with_zcost):
    def __init__(self, Ga, **kwargs):
        pole_optimizer_with_zcost.__init__(self, Ga, \
                                           optcase=39, \
                                           imagline=-8, \
                                           realline=-10, \
                                           secondweight=50.0, \
                                           **kwargs)


class pole_optimizer10(pole_optimizer_with_zcost):
    def __init__(self, Ga, **kwargs):
        pole_optimizer_with_zcost.__init__(self, Ga, \
                                           optcase=40, \
                                           imagline=-10, \
                                           realline=-10, \
                                           secondweight=50.0, \
                                           **kwargs)

    
class pole_optimizer11(pole_optimizer_with_zcost):
    def __init__(self, Ga, **kwargs):
        pole_optimizer_with_zcost.__init__(self, Ga, \
                                           optcase=41, \
                                           imagline=-12, \
                                           realline=-12, \
                                           secondweight=50.0, \
                                           **kwargs)


class pole_optimizer12(pole_optimizer_with_zcost):
    def __init__(self, Ga, **kwargs):
        pole_optimizer_with_zcost.__init__(self, Ga, \
                                           optcase=42, \
                                           imagline=-10, \
                                           realline=-12, \
                                           secondweight=50.0, \
                                           **kwargs)


class pole_optimizer13(pole_optimizer_with_zcost):
    def __init__(self, Ga, **kwargs):
        pole_optimizer_with_zcost.__init__(self, Ga, \
                                           optcase=43, \
                                           imagline=-12, \
                                           realline=-12, \
                                           secondweight=50.0, \
                                           **kwargs)
        self.second_imag_zeta = 0.975
