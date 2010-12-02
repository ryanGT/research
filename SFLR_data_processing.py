from scipy import *
import numpy

import re, os

#import controls

import txt_data_processing
#reload(txt_data_processing)

import txt_mixin

import plotting_mixin

from rwkmisc import rowwise, colwise

from IPython.Debugger import Pdb


class SFLR_Data_File(txt_data_processing.Data_File, \
                     plotting_mixin.item_that_plots):
    def clean_mat_rows(self):
        """This is one of several functions used to load the system
        matrices from the comments at the top of a saved data file.
        See self._find_matrices_in_header for the full process.

        This function eliminates the comments sign and any leading
        spaces from each line as well as spaces around square
        brackets.  Any remaining spaces are assumed to be between
        elements of a row and are replaced with commas."""
        pat1 = '^# +'
        self.mat_rows.replaceallre(pat1,'')
        pat2 = ' *\[ +'
        self.mat_rows.replaceallre(pat2,'[')
        pat3 = ' +\] *'
        self.mat_rows.replaceallre(pat3,']')
        pat4 = '[ \t]+'
        self.mat_rows.replaceallre(pat4,', ')


    def find_ind_mat_lines(self):
        """One of several functions to extract the system matrices
        from the data file comments.  See
        self._find_matrices_in_header for the full process.

        Each line with an equal sign is assumed to be the first line
        of a new system matrix."""
        inds = self.mat_rows.findall('=')
        prev_ind = inds.pop(0)
        inds.append(None)
        mat_code = []
        for ind in inds:
            cur_lines = self.mat_rows[prev_ind:ind]
            mat_code.append(cur_lines)
            prev_ind = ind
        self.mat_chunks = mat_code
        

    def chunk_to_key_and_code(self, chunk):
        """One of several functions to extract the system matrices
        from the data file comments.  See
        self._find_matrices_in_header for the full process.

        Each line with an equal sign is assumed to be of the form
        matname = [[1, 2, 3, 4].

        This code grabs matname, drops the equal sign and puts the
        rest of the code in array(...).  The code lines are joined
        with a comma.  Each line is assumed to contain one row.
        """
        first_line = chunk.pop(0)
        key, rest = first_line.split('=',1)
        key = key.strip()
        rest = rest.strip()
        if key == 'E':
            return key, rest
        assert rest[0] == '[', "expected mat chunk to start with ["
        rest = 'array(' + rest
        if chunk:#there was more than two lines of chunk before popping the first line
            tail = ', ' + ', '.join(chunk)
        else:
            tail = ''
        code = rest + tail + ')'
        return key, code


    def build_mat_code_dict(self):
        """One of several functions to extract the system matrices
        from the data file comments.  See
        self._find_matrices_in_header for the full process."""
        mydict = {}
        for chunk in self.mat_chunks:
            key, code = self.chunk_to_key_and_code(chunk)
            mydict[key] = code
        self.mat_code_dict = mydict
        

    def assign_mat_code_to_attrs(self):
        """One of several functions to extract the system matrices
        from the data file comments.  See
        self._find_matrices_in_header for the full process."""
        for key, code in self.mat_code_dict.iteritems():
            temp = eval(code)
            setattr(self, key, temp)
            
        
    def _find_matrices_in_header(self):
        """This is the main method for extracting the system matrices
        from the comments at the beginning of the data file."""
        if not hasattr(self, 'header'):
            self._get_header()
        self.mat_rows = txt_mixin.txt_list(self.header[1:-1])
        self.clean_mat_rows()
        self.find_ind_mat_lines()
        self.build_mat_code_dict()
        self.assign_mat_code_to_attrs()
        

    def __init__(self, path=None, **kwargs):
        txt_data_processing.Data_File.__init__(self, path, **kwargs)
        if hasattr(self, 'tildey'):
            self.y_tilde = self.tildey
        self.plot_vars = ['u','theta','a','v', \
                          'y_tilde','v_obs']
        self._find_matrices_in_header()


    def _get_label(self, key):
        found = 0
        if hasattr(self, 'plot_labels'):
            if self.plot_labels.has_key(key):
                found = 1
                return self.plot_labels[key]
        if not found:
            return key
        
        
    def plot(self, fi=1, fig=None, clear=True, \
             plot_vars=None, linestyles=None, \
             **kwargs):
        if plot_vars is None:
            plot_vars = self.plot_vars

        if linestyles is None:
            linestyles = ['-']*len(plot_vars)
            
        ax = self._prep_ax(fi=fi, fig=fig, clear=clear)
        t = self.t
        for key, ls in zip(plot_vars, linestyles):
            if hasattr(self, key):
                vect = getattr(self, key)
                label = self._get_label(key)
                ax.plot(t, vect, label=label, \
                        linestyle=ls, **kwargs)

        self.label_axis()


    def set_legend(self, loc=1):
        self.ax.legend(loc=loc)

        
    def plot_states_all_on_one(self, fi=2, fig=None, clear=True):
        ax = self._prep_ax(fi=fi, fig=fig, clear=clear)
        keep_going = True
        n = 0
        while keep_going:
            attr = 'X_tilde%i' % n
            if hasattr(self, attr):
                vect = getattr(self, attr)
                ax.plot(self.t, vect, label=attr)
                n += 1
            else:
                keep_going = False
                break
        self.label_axis()


    def plot_states_many_figs(self, fi=3, fig=None, clear=True, \
                              label_fmt=None):
        if label_fmt is None:
            label_fmt = 'X_tilde%i'
        keep_going = True
        n = 0
        while keep_going:
            attr = 'X_tilde%i' % n
            label = label_fmt % n
            if hasattr(self, attr):
                ax = self._prep_ax(fi=fi, fig=fig, clear=clear)                
                vect = getattr(self, attr)
                ax.plot(self.t, vect, label=label)
                ax.set_title(attr)
                n += 1
                fi += 1
            else:
                keep_going = False
                break
        self.label_axis()



    def plot_states_from_verification(self, fi=3, clear=False, \
                                      attr='X_tilde', \
                                      label_attr=None):
        if label_attr is None:
            label_attr = attr
        mat = getattr(self, attr)
        mat = colwise(mat)
        nr, nc = mat.shape
        for n in range(nc):
            cur_col = mat[:,n]
            label = label_attr + '%i' % n
            ax = self._prep_ax(fi=fi, clear=clear)                
            ax.plot(self.t, cur_col, label=label)
            #ax.set_title(label)
            fi += 1


    def plot_obs_terms_from_verification(self, fi=300, clear=True, \
                                         ind=0, \
                                         label_suffix='', \
                                         many_plots=True):
        keep_going = True
        n = 1
        while keep_going:
            attr = 'term%i' % n
            label_attr = attr + label_suffix
            if hasattr(self, attr):
                ax = self._prep_ax(fi=fi, clear=clear)                
                vect = getattr(self, attr)
                ax.plot(self.t, vect[ind,:], label=label_attr)
                ax.set_title(attr)
                self.label_axis()
                n += 1
                if many_plots:
                    fi += 1
                else:
                    clear = False
            else:
                keep_going = False
                break
        if many_plots:
            ax = self._prep_ax(fi=fi, clear=clear)
            #vect = self.term1[ind, :] + self.term2[ind,:] + self.term3[ind,:]        
            vect = self.term1[ind, :] + self.term3[ind,:]
            ax.plot(self.t, vect, label='all terms '+label_suffix)
            ax.plot(self.t, self.X_tilde[ind, :], label='X_tile_' + str(ind) + ' ' + label_suffix)

            self.label_axis()


    def run_observer(self, i):
        #self.term1[:,i] = squeeze(colwise(dot(self.FO, self.X_tilde[:,i-1])))
        self.term1[:,i] = squeeze(colwise(dot(self.G, self.X_tilde[:,i-1])))
        self.term2[:,i] = squeeze(colwise(self.H*self.vvect[i-1]))
        Y_tilde_float_i = squeeze(dot(self.C, self.X_tilde[:,i-1]))
        #Y_tilde_int_i = Y_tilde_float_i.astype(int)
        #self.term3[:,i] = squeeze(colwise(dot(self.Ke, self.Yvect[:,i]-Y_tilde_int_i)))
        self.term3[:,i] = squeeze(colwise(dot(self.Ke, self.Yvect[:,i-1]-Y_tilde_float_i)))

        ## if term1.any() or term2.any() or term3.any():
        ##     Pdb().set_trace()
        next_x_tilde = self.term1[:,i] + self.term2[:,i] + self.term3[:,i]
            #colwise(dot(self.Ke, self.Yvect[:,i-1]))
        ## if term1.any() or term2.any() or term3.any():
        ##     if self.debug_ind < 2:
        ##         print('i = '+str(i))
        ##         print('term1 = ' + str(term1))
        ##         print('term2 = ' + str(term2))
        ##         print('term3 = ' + str(term3))
        ##         print('next_x_tilde = ' + str(next_x_tilde))
        ##         print('-'*20)
        ##         self.debug_ind += 1

        self.X_tilde[:,i] = squeeze(next_x_tilde)
    

    def check_observer(self):
        nr, nc = self.G.shape
        n = len(self.t)
        self.n = n
        self.X_tilde = zeros((nr, n))
        self.term1 = zeros((nr, n))
        self.term2 = zeros((nr, n))
        self.term3 = zeros((nr, n))

        self.Yvect = rowwise(self.theta)
        self.vvect = self.v
        #self.FO = self.G - dot(self.Ke, self.C)

        self.debug_ind = 0
        
        for i in range(1,n):
            self.run_observer(i)
        
        
#collables:
#t	n	u	v	$\theta$	a	$\hat{\theta}_d$
SFLR_col_map =  {0:'t', 1:'n', 2:'u', 3:'v', \
                 4:'theta', 5:'a', 6:'theta_d_hat'}


class SFLR_Exp_Data_File(SFLR_Data_File):
    def __init__(self, path=None, substr=None, \
                 col_map=SFLR_col_map, plot_vars=None, \
                 **kwargs):        
        txt_data_processing.Data_File.__init__(self, path, \
                                               col_map=col_map, \
                                               **kwargs)
        if plot_vars is None:
            plot_vars = ['u','theta','a','v', \
                         'theta_d_hat']
        self.plot_vars = plot_vars
        labels = ['\\theta_d','\\theta','\\ddot{x}','v','\\hat{\\theta}_d']
        if substr is not None:
            labels[1:] = [item + '_{%s}' % substr for item in labels[1:]]
        labels = ['$%s$' % item for item in labels]
        self.plot_labels = dict(zip(self.plot_vars, labels))
        

    

class SFLR_Exp_Step_Response_Set(txt_data_processing.Data_Set):
    """A class for loading a group of experimental SFLR step response
    txt data files."""
    def __init__(self, filepaths, substrs=None):
        self.filepaths = filepaths
        self.col_map = SFLR_col_map
        if substrs is None:
            substrs = [None]*len(filepaths)
        self.substrs = substrs
        #self.title_dict = title_dict
        self.folder, filename = os.path.split(self.filepaths[0])
        self.Load_Data()

    
    def Load_One_File(self, path, substr=None):
        curfile = SFLR_Exp_Data_File(path, col_map=self.col_map,\
                                     substr=substr)
        return curfile


    def Load_Data(self):
        data_files = []
        for filename, substr in zip(self.filepaths, self.substrs):
            curfile = self.Load_One_File(filename, substr=substr)
            data_files.append(curfile)
        self.data_files = data_files
    

    def Overlay_Step_Responses(self, fi=1, \
                               plot_vars=['theta','a'], \
                               linestyles=None):
        first = 1
        if linestyles is None:
            linestyles = ['-']*len(self.data_files)
        for df,lt in zip(self.data_files, linestyles):
            if first:
                first = 0
                cur_plot_vars = ['u'] + plot_vars
                df.plot(fi=fi, fig=None, clear=True, \
                        plot_vars=cur_plot_vars, linestyles=[lt])
            else:
                df.plot(fi=fi, fig=None, clear=False, \
                        plot_vars=plot_vars, linestyles=[lt])

        
