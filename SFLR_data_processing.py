from scipy import *
import numpy

import re

#import controls

import txt_data_processing
reload(txt_data_processing)

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
        self.y_tilde = self.tildey
        self.plot_vars = ['u','theta','a','v', \
                          'y_tilde','v_obs']
        self._find_matrices_in_header()
        
        
    def plot(self, fi=1, fig=None, clear=True):
        ax = self._prep_ax(fi=fi, fig=fig, clear=clear)
        t = self.t
        for key in self.plot_vars:
            if hasattr(self, key):
                vect = getattr(self, key)
                ax.plot(t, vect, label=key)

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


    def plot_states_many_figs(self, fi=3, fig=None, clear=True):
        keep_going = True
        n = 0
        while keep_going:
            attr = 'X_tilde%i' % n
            if hasattr(self, attr):
                ax = self._prep_ax(fi=fi, fig=fig, clear=clear)                
                vect = getattr(self, attr)
                ax.plot(self.t, vect, label=attr)
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
        Y_tilde_int_i = Y_tilde_float_i.astype(int)
        self.term3[:,i] = squeeze(colwise(dot(self.Ke, self.Yvect[:,i]-Y_tilde_int_i)))

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
        
        
