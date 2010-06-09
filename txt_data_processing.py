from scipy import *
from scipy import io
import glob, os, copy, re, sys
import mplutil
reload(mplutil)
import time

#import rst_creator
#reload(rst_creator)

import txt_mixin
import rwkdataproc
import rwkbode

from IPython.Debugger import Pdb

pat = '\$\\\\(.*)\$'
p = re.compile(pat)


td_map_attrs = ['filepaths', 'col_map', 'time_label', 'title_dict']

bode_keys = ['f', 'bode_list', 'col_map', \
             'filepaths', 'title_dict']


def thresh(iterin, value, startind=0, above=1):
    if above:
        mybools = iterin >= value
    else:
        mybools = iterin <= value
    myinds = arange(len(iterin))
    keepinds = myinds[squeeze(mybools)]
    inds2 = keepinds >= startind
    keep2 = keepinds[inds2]

    if len(keep2)==0:
        return -1
    else:
        return keep2.min()


###########################
#
#  Variable Down Sampling of data for rwkbode.compress
#
###########################

def FindLogInds(f, minf, maxf, N=50):
    logf = logspace(log10(minf), log10(maxf), N)
    logfinds = [thresh(f, item) for item in logf]
    return logfinds


def FindVarLogInds(f, minf, maxf, N):
    """N is points per decade."""
    nd = log10(maxf)-log10(minf)
    NN = int(nd*N+0.5)
    return FindLogInds(f, minf, maxf, NN)


def BuildMask(inds, vect):
    mymask = [item in inds for item in range(len(vect))]
    return mymask

###########################

def clean_key(keyin):
    if keyin[0] == '#':
        keyin = keyin[1:]
    temp = p.sub('\\1', keyin)
    return temp


from rwkmisc import my_import

def import_mod(mod_path):
    folder, mod_name = os.path.split(mod_path)
    if folder:
        abspath = os.path.abspath(folder)
        if abspath not in sys.path:
            sys.path.append(abspath)
    my_mod = my_import(mod_name)
    return my_mod


def filterx(label):
    xlabels = ['t','n']
    if label in xlabels:
        return False
    else:
        return True
    
    
class Data_File(object):
    """A class for representing one txt data file assumed to be
    composed of data in columns with a certain number of initial rows
    at the top of the file that should be skipped.

    col_map is a dictionary with column numbers as the keys at the
    names of attributes to store those columns as as the values, i.e.

    colmap = {0:'t', 1:'u', 2:'y'}

    would lead to the equivalent of

    self.t = data[:,0]
    self.u = data[:,1]
    self.y = data[:,2]

    (the code actually uses setattr(self, key, data[:,ind])).
    """
    def __init__(self, path, col_map={}, delim='\t', skiprows=None):
        self.path = path
        self.col_map = col_map
        self.delim = delim
        self.skiprows = skiprows
        if os.path.exists(self.path):
            self.Load_Data()


    def _load_raw(self):
        temp = txt_mixin.txt_file_with_list(self.path)
        self.raw_list = temp.list
        self.N = len(self.raw_list)
        

    def sniff(self):
        """Attempt to determine the correct value for self.skiprows by
        looking for the first row that doesn't start with a #.

        If the file doesn't start with lines that start with #, try
        loading data and increment skiprows each time scipy.loadtxt
        throws a ValueError."""
        skip = 0
        self._load_raw()
        while skip < self.N:
            if self.raw_list[skip][0] == '#':
                skip += 1
            else:
                if skip > 0:
                    self.skiprows = skip
                    return self.skiprows
                else:
                    break
        #if we got to this point, we didn't find any lines starting
        #with a #
        skip = 0
        while skip < self.N:
            try:
                data = loadtxt(self.path, delimiter=self.delim, \
                               skiprows=skip)
                self.skiprows = skip
                return self.skiprows
            except ValueError:
                skip += 1


    def get_labels(self):
        """Read the labels from the data file, assuming they are in
        the row just above the data itself."""
        if self.skiprows is None:
            self.sniff()
        if not hasattr(self, 'raw_list'):
            self._load_raw()
        labels = self.raw_list[self.skiprows-1]
        raw_labels = labels.split(self.delim)
        self.labels = map(clean_key, raw_labels)
        


    def _load_data(self):
        if self.skiprows is None:
            self.sniff()
        skip = self.skiprows
        while skip < self.N:
            try:
                self.data = loadtxt(self.path, delimiter=self.delim, \
                                    skiprows=skip)
                self.skiprows = skip
                return self.data
            except ValueError:
                skip += 1
        

    def Load_Data(self):
        """Load the data and store the columns in appropriate
        attributes."""
        if not self.col_map:
            self.get_labels()
            N = len(self.labels)
            keys = range(N)
            self.col_map = dict(zip(keys, self.labels))
        data = self._load_data()
        for ind, attr in self.col_map.iteritems():
            setattr(self, attr, data[:,ind])


    def get_time_axis(self, fignum):
        from pylab import figure
        fig = figure(fignum)
        ax = fig.add_subplot(111)
        return ax
    
        
    def Time_Plot(self, labels=None, ax=None, fignum=1,
                  clear=True, legloc=4, legend_dict={}, \
                  ylabel='Voltage (counts)', \
                  basename=None, save=False, \
                  ext='.png', fig_dir='', **plot_opts):                  
        if ax is None:
            ax = self.get_time_axis(fignum)
        if clear:
            ax.clear()
        if labels is None:
            labels = filter(filterx, self.labels)
        for label in labels:
            curvect = getattr(self, label)
            mplutil.plot_vect(ax, self.t, curvect, clear=False, \
                              ylabel=ylabel, **plot_opts)
        #ax.legend(labels, loc=legloc)
##         if basename and savefigs:
##             fig_name = basename+'_time_plot'
##             fig_name += '_%s' % '_'.join(labels)
##             fig_path = os.path.join(fig_dir, fig_name)
##             mplutil.mysave(fig_path, fig, ext=ext)
        return ax


class Data_Set(object):
    """A class for loading a group of related txt data files and
    making all the data columns into one matrix per
    variable, so that each column of each matrix represents a
    different experimental test.

    If pattern is None, no attempt will be made to load the data.
    This is useful only if you are loading a saved data set from a
    module, otherwise, you want pattern to be a glob pattern for your
    txt files."""
    def __init__(self, pattern=None, time_col=0, time_label='t', \
                 col_map={}, title_dict={}, \
                 delim='\t', skiprows=None):
        self.pattern = pattern
        self.time_col = time_col
        self.time_label = time_label
        self.col_map = col_map
        self.title_dict = title_dict
        self.delim = delim
        self.skiprows = skiprows
        if pattern is not None:
            self.filepaths = glob.glob(pattern)
            self.filepaths.sort()
            assert len(self.filepaths) > 0, \
                   "Could not find any files matching pattern: %s" % \
                   pattern
            self.folder, filename = os.path.split(self.filepaths[0])
            self.Load_Data()


    def _get_list_of_figs(self, N, fignum=1, size=None):
        from pylab import figure
        figs = []
        for i in range(N):
            fig = figure(fignum+i, size)
            figs.append(fig)
        return figs
            


    def Append_One_File(self, path):
        curfile = Data_File(path, col_map=self.col_map, \
                            delim=self.delim, \
                            skiprows=self.skiprows)
        if not self.col_map:
            self.col_map = curfile.col_map
        #copy the time vector if one does exist
        if not hasattr(self, self.time_label):
            if hasattr(curfile, self.time_label):
                t = getattr(curfile, self.time_label)
                setattr(self, self.time_label, t)
        #copy each column of data in col_map, either adding it to an
        #existing matrix (if it exists), or creating a new matrix from
        #the column
        attrs = self.col_map.values()
        attrs.remove(self.time_label)
        for attr in attrs:
            new_col = getattr(curfile, attr)
            if hasattr(self, attr):
                existing_data = getattr(self, attr)
                cur_mat = column_stack([existing_data, new_col])
            else:
                cur_mat = column_stack([new_col])
            setattr(self, attr, cur_mat)



    def Load_Data(self):
        """Load the data columns of the data files and store them as
        matrices in the corresponding attributes of self."""
        for filename in self.filepaths:
            self.Append_One_File(filename)
                
            
    def Time_Plot(self, labels=['y'], ax=None, fignum=1,
                  legend_dict={}, ylabel=None, \
                  basename=None, savefigs=True, \
                  ext='.png', fig_dir=''):                  
        if ax is None:
            from pylab import figure
            fig = figure(fignum)
            fig.clf()
            ax = fig.add_subplot(111)
        for label in labels:
            cur_mat = getattr(self, label)
            curlabel = ylabel
            if curlabel is None:
                if self.title_dict.has_key(label):
                    curlabel = self.title_dict[label]
            mplutil.plot_cols(ax, self.t, cur_mat, ylabel=curlabel, \
                              clear=False, labels=curlabel)
        if basename and savefigs:
            fig_name = basename+'_time_plot'
            fig_name += '_%s' % '_'.join(labels)
            fig_path = os.path.join(fig_dir, fig_name)
            mplutil.mysave(fig_path, fig, ext=ext)


    def _get_labels(self, exclude=['n']):
        labels = self.col_map.values()
        labels.remove(self.time_label)
        labels = [item for item in labels if item not in exclude]
        return labels


    def _build_data_matrix(self, labels):
        col_list = []
        for label in labels:
            cur_mat = getattr(self, label)
            col_list.append(cur_mat)
        data = column_stack(col_list)
        return data
    
        
    def Time_Plots(self, labels=None, figs=None, fignum=1, \
                   exclude=['n'], overlay=False, **kwargs):
        if labels is None:
            labels = self._get_labels(exclude=exclude)
        N = len(labels)
        if figs is None:
            figs = self._get_list_of_figs(N=N, fignum=fignum, size=None)
        for item, fig in zip(labels, figs):
            if not overlay:
                ax = fig.add_subplot(111)
            if self.title_dict.has_key(item):
                my_ylabel = self.title_dict[item]
            else:
                my_ylabel = item
            self.Time_Plot(labels=[item], ax=ax, ylabel=my_ylabel,
                           **kwargs)


    def Overlay_Time_Plots(self, labels=None, fig=None, fignum=1, \
                           exclude=['n'], clear=True, **kwargs):
        if labels is None:
            labels = self._get_labels(exclude=exclude)
        if fig is None:
            import pylab
            fig = pylab.figure(fignum)
        ax = fig.add_subplot(111)
        if clear:
            ax.clear()
        data = self._build_data_matrix(labels)
        mplutil.plot_cols(ax, self.t, data)
        



    def save(self, module_name, signal_names=None, overwrite=True):
        """Create a dictionary containing the time vector and the
        matrices associated with signal_names (a list of attrs of
        self) and save that dictionary to the module module_name using
        scipy.io.save_as_module"""
        if signal_names is None:
            signal_names = self.col_map.values()
        mydict = {}
        t = getattr(self, self.time_label)
        mydict[self.time_label] = t
        mydict['signal_names'] = signal_names
        for attr in td_map_attrs:
            mydict[attr] = getattr(self, attr)
        for key in signal_names:
            mat = getattr(self, key)
            mydict[key] = mat

        if overwrite:
            old_files = glob.glob(module_name+'.*')
            for curfile in old_files:
                os.remove(curfile)
        io.save_as_module(module_name, mydict)



def load_time_domain_data_set(module_name):
    """Load a time domain data set that was saved to a module using
    the save method."""
    my_mod = import_mod(module_name)
    my_data_set = Data_Set(pattern=None)
    t = getattr(my_mod, my_mod.time_label)
    setattr(my_data_set, my_data_set.time_label, t)
    my_data_set.filepaths = my_mod.filepaths
    all_attrs = my_mod.signal_names + td_map_attrs 
    for attr in all_attrs:
        mat = getattr(my_mod, attr)
        setattr(my_data_set, attr, mat)
    return my_data_set
        
    
    
        

class Bode_Options(object):
    """This class exists to encapsulate the options associated with
    creating a Bode plot from to columns of a text file.
    Bode_Data_Set will be passed a list of these options to create its
    Bodes.

    input_label and output_label refer to valid values of the col_map
    dictionary of the Bode_Data_Set, i.e. the data set will have an
    attribute of that name that is a matrix of time domain data so
    that something like

    input_matrix = getattr(self, input_label)
    output_matrix = getattr(self, output_label)

    will retreive the appropriate data.

    seedphase and seedfreq are used to massage the phase of the
    experimental Bode data.

    freqlim, phaselim, maglim, and cohlim are used to zoom the final
    plots correctly."""
    def __init__(self, input_label, output_label, \
                 seedfreq=None, seedphase=None, \
                 freqlim=None, phaselim=None, maglim=None, \
                 cohlim=None):
        self.input_label = input_label
        self.output_label = output_label
        self.seedfreq = seedfreq
        self.seedphase = seedphase
        self.freqlim = freqlim
        self.phaselim = phaselim
        self.maglim = maglim
        self.cohlim = cohlim
                 
                 

class Bode_Data_Set(Data_Set):
    """A class for loading a group of related txt data files and
    making Bode plots, including averaging and coherence.  A data
    channel from each test becomes a column in a matrix of the data
    set, so that each column of a data set matrix represents a
    different experimental test.

    bode_list is a list of Bode_Options instances specifying which
    channels should be used for each Bode plot, along with other
    optional fine tuning parameters."""
    def __init__(self, pattern, bode_list, \
                 time_col=0, time_label='t', \
                 col_map={}, coh=False, coh_size=(6,9), \
                 #coh_size=(8,12), \
                 title_dict={}, \
                 delim='\t', skiprows=None):
        Data_Set.__init__(self, pattern, time_col=time_col, \
                          time_label=time_label, col_map=col_map, \
                          title_dict=title_dict, \
                          delim=delim, skiprows=skiprows)
        self.bode_list = bode_list
        self.N = len(self.bode_list)
        self.coh_size = coh_size
        self.coh = coh


    def _ave_dict(self, attr, dict_in={}, save_attr=None):
        """Make a list of bode dicts from the bode list returned by
        getattr(self, attr) and place that dict list in dict_in.  This
        is part of saving Bode_Data_Set instances using scip.io.save_as_module.

        save_attr refers to the keys used to insert the dict list into
        dict_in.  If it is not given, attr is used."""
        if hasattr(self, attr):
            my_bode_list = getattr(self, attr)
            if save_attr is None:
                save_attr = attr
            dict_list = []
            for bode in my_bode_list:
                cur_dict = bode.to_dict()
                dict_list.append(cur_dict)
            dict_in[save_attr] = dict_list
        return dict_in
        

    def _f_into_dict(self, f_attr, dict_in={}, save_f_attr=None):
        """Similar to _ave_dict above, get the attribute f_attr and
        put it in dict_in using the key save_f_attr or f_attr if
        save_f_attr is not given."""
        if hasattr(self, f_attr):
            if save_f_attr is None:
                save_f_attr = f_attr
            fvect = getattr(self, f_attr)
            dict_in[save_f_attr] = fvect
        return dict_in
        

    def _copy_bode_keys(self, dict_in={}, key_list=None):
        if key_list is None:
            key_list = bode_keys
        for key in key_list:
            dict_in[key] = getattr(self, key)
            
        
    def build_ave_dict(self, attrs=['avebodes', 'trunc_avebodes', \
                                    'compressed_avebodes'],
                       f_attrs=['f','trunc_f','compressed_f']):
        """Build a dictionary that can be used with scipy.io.save_as_module.
        This also includes making a list of dictionaries for each bode
        in self.avebodes or whatever attributes are listed in attrs."""
        if not hasattr(self, 'avebodes'):
            self.Calc_Ave_Bodes()
        if not hasattr(self, 'f'):
            self.Make_Freq_Vect()
        keys = bode_keys
        mydict = {}
        for key in keys:
            mydict[key] = getattr(self, key)
        for attr in attrs:
            mydict = self._ave_dict(attr, mydict)
        for f_attr in f_attrs:
            mydict = self._f_into_dict(f_attr, mydict)
        return mydict


    def _delete_old(self, mod_name):
        old_files = glob.glob(mod_name+'.*')
        for curfile in old_files:
            os.remove(curfile)
        

    def save_ave(self, mod_name, **kwargs):
        mydict = self.build_ave_dict(**kwargs)
        self._delete_old(mod_name)
        io.save_as_module(mod_name, mydict)


    def save_as(self, mod_name, attr='compressed_avebodes', \
                f_attr='compressed_f', \
                save_attr='avebodes', save_f='f'):
        """This function is used to save one attr, such as compressed
        or trunc avebodes as another, such as avebodes, for the sake
        of using the saved data with curve fitting."""
        mydict = {}
        self._ave_dict(attr, dict_in=mydict, save_attr=save_attr)
        self._f_into_dict(f_attr, dict_in=mydict, \
                          save_f_attr=save_f)
        my_keys = copy.copy(bode_keys)
        my_keys.remove('f')
        self._copy_bode_keys(dict_in=mydict, key_list=my_keys)
        self._delete_old(mod_name)        
        io.save_as_module(mod_name, mydict)
    

    def Calc_Spectra(self):
        self.spectra = []
        for cur_bode in self.bode_list:
            in_label = cur_bode.input_label
            out_label = cur_bode.output_label
            in_mat = getattr(self, in_label)
            out_mat = getattr(self, out_label)
            spec = rwkdataproc.CalcSpectra(in_mat, out_mat, \
                                           input=in_label, \
                                           output=out_label)
            self.spectra.append(spec)


    def PhaseMassage(self, attr):
        """Massage the phase of the bodes in list attr if
        self.seed_freqs and self.seed_phases are set."""
        if not hasattr(self, 'f'):
            self.Make_Freq_Vect()
        bode_list = getattr(self, attr)
        for bode in bode_list:
            bode.PhaseMassage(self.f)


    def _calc_bodes(self, attr, func):
        if not hasattr(self, 'spectra') or (not self.spectra):
            self.Calc_Spectra()
        outlist = []
        copylist = ['seedfreq', 'seedphase', 'freqlim', 'phaselim', \
                    'maglim', 'cohlim']
        for spec, bode_opts in zip(self.spectra, self.bode_list):
            cur_bode = func(spec)
            for copy_attr in copylist:
                cur_val = getattr(bode_opts, copy_attr)
                setattr(cur_bode, copy_attr, cur_val)
            outlist.append(cur_bode)

        setattr(self, attr, outlist)
        self.PhaseMassage(attr)#this won't do anything if
                                  #self.seed_freqs and
                                  #self.seed_phases aren't set
        

    def Calc_Bodes(self):
        myfunc = rwkbode.BodeFromSpectra
        myattr = 'bodes'
        self._calc_bodes(myattr, myfunc)


    def Calc_Ave_Bodes(self):
        if len(self.filepaths) == 1:
            #don't average one file, instead set self.avebodes with a
            #placeholder (a list of Nones) to make Bode_Plot happy,
            #but to also throw errors if the script tries to do
            #anything with avebodes.
            self.avebodes = [None]*len(self.bode_list)
            return
        myfunc = rwkbode.AveBodeFromSpectra
        myattr = 'avebodes'
        self._calc_bodes(myattr, myfunc)


    def Make_Freq_Vect(self):
        t = getattr(self, self.time_label)
        self.f = rwkdataproc.makefreqvect(t)


    def _get_list_of_figs(self, N=None, fignum=1, size=None):
        if N is None:
            N = len(self.bode_list)
        figs = Data_Set._get_list_of_figs(self, N, fignum=fignum, \
                                          size=size)
        return figs


    def _build_auto_title(self, bode, title_dict):
        if title_dict.has_key(bode.output):
            outstr = title_dict[bode.output]
        else:
            outstr = bode.output 
        if title_dict.has_key(bode.input):
            instr = title_dict[bode.input]
        else:
            instr = bode.input
        titlestr = '%s/%s' % (outstr, instr)
        return titlestr
    

    def _set_title(self, fig, bode, title=None, autotitle=True,
                   title_dict=None):
        if title_dict is None:
            title_dict = self.title_dict
        if title is not None:
            mytitle = title
        elif autotitle:
            mytitle = self._build_auto_title(bode, \
                                             title_dict=title_dict)
        else:
            mytitle = ''
        if mytitle:
            mplutil.SetTitle(fig, mytitle)


    def _check_bode_existance(self):
        if not hasattr(self, 'avebodes') or (not self.avebodes):
            self.Calc_Ave_Bodes()
        if not hasattr(self, 'bodes') or (not self.bodes):
            self.Calc_Bodes()
        if not hasattr(self, 'f'):
            self.Make_Freq_Vect()


    def _bode_plot(self, bode_attr, f_attr='f', figs=None, \
                   fignum=1, size=None, clear=True, \
                   func=rwkbode.BodeCohPlot, **plotargs):
        if func == rwkbode.BodeCohPlot and size is None:
            size = self.coh_size
            self.coh = True
        if figs is None:
            kwargs = {}
            if size is not None:
                kwargs['size'] = size
            figs = self._get_list_of_figs(fignum=fignum, **kwargs)
        bode_list = getattr(self, bode_attr)
        f = getattr(self, f_attr)
        for bode, fig in zip(bode_list, figs):
            func(100, f, bode, fig=fig, clear=clear, **plotargs)
            #note that the 100 is a dummy, it won't be used since fig
            #is passed in
        return figs
        


    def _set_plot_opts(self, bode_attr, figs):
        bode_list = getattr(self, bode_attr)
        #Pdb().set_trace()
        for bode, fig in zip(bode_list, figs):
            mplutil.set_Bode_opts(fig, bode, coh=self.coh)  
        

    def _set_titles(self, bode_attr, figs, \
                    title=None, autotitle=True, title_dict=None):
        bode_list = getattr(self, bode_attr)
        for bode, fig in zip(bode_list, figs):
            self._set_title(fig, bode, title, autotitle, \
                            title_dict=title_dict)


    def find_bode(self, output_label, input_label, attr='avebodes'):
        bode_list = getattr(self, attr)
        for item in bode_list:
            if (item.output == output_label) and (item.input == input_label):
                return item


    def Bode_Plot2(self, attr='avebodes', f_attr='f', \
                   figs=None, fignum=1, size=None, \
                   func=rwkbode.BodeCohPlot, **plotargs):
        figs = self._bode_plot(attr, f_attr=f_attr, \
                               figs=figs, fignum=fignum, \
                               size=size, func=func, **plotargs)
        self._set_plot_opts(attr, figs)
        self._set_titles(attr, figs)
        return figs
        

    def Bode_Plot(self, overlay_ave=True, figs=None, fignum=1, \
                  title=None, autotitle=True, coh=False, \
                  title_dict={}, size=(9,15), \
                  basename=None, savefigs=True, \
                  ext='.png', fig_dir='', rst_file=None):
        """Generate Bode plots for the Bode_Data_Set.  If coh is True,
        call rwkbode.BodeCohPlot, otherwise call rwkbode.GenBodePlot.

        If title is None and autotitle is True, then the title of the
        Bode plot will be generated automatically as output/input,
        where output and input are the strings associated with the
        input and output of the Bodes.  title_dict is a dictionary
        that can be used to map those attribute strings to something
        else, including LaTeX strings.

        Note that if you want a figure to be saved, you must specify
        basename.  savefigs is used to turn off saving when basename
        is specified."""
        if coh:
            myfunc = rwkbode.BodeCohPlot
        else:
            myfunc = rwkbode.GenBodePlot

        self._check_bode_existance()

        if figs is None:
            kwargs = {}
            if coh:
                kwargs['size'] = size
            figs = self._get_list_of_figs(fignum=fignum, **kwargs)
        for bode, avebode, fig in zip(self.bodes, self.avebodes, figs):
            if coh:
                #do only the average Bode on a BodeCohPlot:
                myfunc(100, self.f, avebode, fig=fig)
            else:
                #otherwise, plot the regular Bodes using GenBodePlot
                #and overlay the average if overlay_ave is True
                myfunc(100, self.f, bode, fig=fig)
                if overlay_ave:
                    myfunc(100, self.f, avebode, fig=fig, \
                           clear=False, linetype='k-')
            mplutil.set_Bode_opts(fig, bode)
            self._set_title(fig, bode, title, autotitle, \
                            title_dict=title_dict)
            if basename:
                fig_name = basename+'_bode'
                if coh:
                    fig_name +='_coh'
                fig_name += '_%s_vs_%s' % (bode.output, bode.input)
                fig_path = os.path.join(fig_dir, fig_name)
                if savefigs:
                    mplutil.mysave(fig_path, fig, ext=ext)
                if rst_file:
                    caption = 'Bode plot for %s/%s.' % (bode.output, \
                                                        bode.input)
                    relpath = os.path.join('figs', fig_name)+ext
                    rst_file.add_figure(relpath, caption=caption)
        return figs

                
            

##     def Report(self, startfi=1, timedomain=True, bode=True, \
##                overlay_ave=True, coh=True, bode_coh_only=False, \
##                time_signals=None, title_dict={},
##                report_dir=None, basename=None, rst_file=None, \
##                savefigs=True, ext='.png'):
##         """Generate a report for the Bode_Data_Set.  The report may
##         include time domain plots for each signal if timedomain is
##         True as well as Bodes with the optional overlays of the
##         averaged Bodes from each individual data file/test if bode is
##         True (and overlay_ave is True).  bode_coh_only = True
##         overrides both of these options and plots just the Bode plots
##         with coherence.

##         Note that here, bode refers to a Bode plot with only magintude
##         and phase; coh refers to Bode plots with coherence (coh, mag,
##         and phase).

##         This function will not really generate a report unless
##         report_dir and basename are specified.  If they are not
##         specified, figures will be created, but not saved and no rst
##         will be generated.  If report_dir and basename are specified,
##         then figures will be saved to report_dir/figs (provided
##         savefigs is True), and an rst file called basename.rst will be
##         created.  basename will also be the first portion of the
##         figure file names.  rst_file can contain an introduction to
##         the report and will append to if given."""
##         gen_report = False
##         if report_dir and basename:
##             gen_report = True
##             if not os.path.exists(report_dir):
##                 os.mkdir(report_dir)
##             fig_dir = os.path.join(report_dir, 'figs')
##             if not os.path.exists(fig_dir):
##                 os.mkdir(fig_dir)
##             if rst_file is None:
##                 rst_name = basename+'.rst'
##                 rst_file = rst_creator.rst_file(rst_name)
##         if not gen_report:
##             savefigs = False
##             basename = None
##             fig_dir = None
##         if bode_coh_only:
##             timedomain = False
##             bode = False
##             coh = True
##         if time_signals is None:
##             time_signals = ['theta','u','v']
##         time_plots = []
##         for key in time_signals:
##             if title_dict.has_key(key):
##                 label = title_dict[key]
##             else:
##                 label = key
##             time_plots.append((key, label))
##         fi = startfi
##         if coh:
##             if rst_file:
##                 rst_file.add_section('Bode Plots with Coherence')
##             self.Bode_Plot(fignum=fi, coh=True, title_dict=title_dict, \
##                            basename=basename, savefigs=savefigs, \
##                            fig_dir=fig_dir, ext=ext, \
##                            rst_file=rst_file)
##             fi += self.N
##         if bode:
##             self.Bode_Plot(fignum=fi, coh=False, \
##                            overlay_ave=overlay_ave, \
##                            title_dict=title_dict, \
##                            basename=basename, savefigs=savefigs, \
##                            fig_dir=fig_dir, ext=ext)
##             fi += self.N
##         if timedomain:
##             for attr, label in time_plots:
##                 self.Time_Plot(labels=[attr], ylabel=label, fignum=fi, \
##                            basename=basename, savefigs=savefigs, \
##                            fig_dir=fig_dir, ext=ext)
##                 fi += 1
##         if gen_report:
##             rst_file.save()
##             rst_file.to_html()
            
        

        
    def Truncate(self, flow, fhigh):
        self._check_bode_existance()
        self.trunc_bodes = copy.deepcopy(self.bodes)
        self.trunc_avebodes = copy.deepcopy(self.avebodes)
        for bode in self.trunc_bodes:
            self.trunc_f = bode.truncate(self.f, flow=flow, fhigh=fhigh)
        for bode in self.trunc_avebodes:
            bode.truncate(self.f, flow=flow, fhigh=fhigh)
        


    def Overlay_Trunc_Bodes(self, ave=True, figs=None, fignum=1, \
                            title_dict=None):
        if ave:
            attr = 'avebodes'
            trunc_attr = 'trunc_avebodes'
            size = self.coh_size
            func = rwkbode.BodeCohPlot
        else:
            attr = 'bodes'
            trunc_attr = 'trunc_bodes'
            size = None
            func = rwkbode.GenBodePlot
        assert hasattr(self, trunc_attr), \
               "You must call self.Truncate before calling\n " + \
               "self.Overlay_Trunc_Bodes."
            
        figs = self._bode_plot(attr, f_attr='f', figs=figs, \
                   fignum=fignum, size=size, clear=True, \
                   func=func)
        self._bode_plot(trunc_attr, f_attr='trunc_f', figs=figs, \
                        clear=False, func=func)
        self._set_plot_opts(attr, figs)
        self._set_titles(attr, figs, title=None, autotitle=True, \
                         title_dict=title_dict)
        return figs
    


    def Truncate_and_Review(self, flow, fhigh, figs=None, fignum=1):
        self.Truncate(flow, fhigh)
        self.Overlay_Trunc_Bodes(ave=True, figs=figs, fignum=fignum)


    def save_bodes_and_time_domain(self, basename):
        td_name = basename + '_time_domain'
        avename = basename + '_avebodes'
        self.save(td_name)
        self.save_ave(avename)


    def Log_Compress_Data(self, freqs, Ns, attr='avebodes', fattr='f'):
        """Logarithmically downsample the data in attr by using the
        compress method of an rwkbode instance.  This is done to
        improve the effectiveness of curve fitting swept sine data
        where higher modes are overemphasized on Bode plots due to the
        semi-logx nature of Bodes.

        freqs is a list of frequencies that mark the boundaries of
        different values for downsampling.  Ns is a list of the N
        values to use for each frequency range where N is points per
        decade.  Note that freqs should be a list with one more entry
        than Ns (each frequency range starts at the end of the
        previous one, so there are len(freqs)-1 ranges)."""
        inds_log = []
        prev_f = freqs.pop(0)
        fvect = getattr(self, fattr)
        for f, N in zip(freqs, Ns):
            cur_inds = FindVarLogInds(fvect, prev_f, f, N)
            inds_log += cur_inds
            prev_f = f
        varinds = unique(inds_log)#unique is from scipy
        varmask = BuildMask(varinds, fvect)
        fvarlog = fvect.compress(varmask)
        bodes = copy.deepcopy(getattr(self, attr))
        for bode in bodes:
            bode.compress(varmask)
        attr_out = 'compressed_'+attr
        self.compressed_f = fvarlog
        setattr(self, attr_out, bodes)
        return varinds, varmask, fvarlog

        
    
allkeys = ['delim', \
           'trunc_avebodes', \
           'trunc_bodes', \
           'pattern', \
           'bode_list', \
           'filepaths', \
           'test', \
           'theta', \
           'folder', \
           'col_map', \
           'coh_size', \
           'N', \
           'time_col', \
           'time_label', \
           'a', \
           'avebodes', \
           'f', \
           'spectra', \
           'trunc_f', \
           'skiprows', \
           'n', \
           'data_store', \
           'u', \
           't', \
           'v', \
           'bodes']


def list_of_dicts_to_bodes(mod_in, attr='avebodes'):
    dict_list = getattr(mod_in, attr)
    bodelist = []
    for cur_dict in dict_list:
        curbode = rwkbode.bode_from_dict(cur_dict)
        bodelist.append(curbode)
    return bodelist

        
def load_avebode_data_set(module_name):
    """Load an avebodes data set that was saved to a module using
    the save_ave method of Bode_Data_Set."""
    my_mod = import_mod(module_name)
    my_data_set = Bode_Data_Set(pattern=None, \
                                bode_list=my_mod.bode_list)
    for key in bode_keys:
        val = getattr(my_mod, key)
        setattr(my_data_set, key, val)
    my_data_set.avebodes = list_of_dicts_to_bodes(my_mod, 'avebodes')
    if hasattr(my_mod, 'trunc_avebodes'):
        my_data_set.trunc_avebodes = list_of_dicts_to_bodes(my_mod, \
                                                            attr='trunc_avebodes')
        my_data_set.trunc_f = my_mod.trunc_f
    return my_data_set


def find_matching_bode(bode1, bode_list):
    for bode in bode_list:
        if bode.input == bode1.input:
            if bode.output == bode1.output:
                return bode
            

def merge_trunc_ave_data_sets(data_set1, data_set2):
    """Merge data_set1 and data_set2, which are presumed to be
    Bode_Data_Set instances presumably from doing focused swept sine
    tests with different frequency ranges."""
    my_keys = ['bode_list', 'col_map', 'title_dict']
    bode_merge_keys = ['mag', 'phase', 'coh']
    bode_copy_keys = ['input', 'output', 'phaselim', \
                      'maglim', 'seedfreq', 'seedphase']
    f1 = column_stack([data_set1.trunc_f])
    f2 = column_stack([data_set2.trunc_f])
    f = row_stack([f1,f2])
    my_data_set = Bode_Data_Set(pattern=None, \
                                bode_list=data_set1.bode_list)
    for key in my_keys:
        val = getattr(data_set1, key)
        setattr(my_data_set, key, val)

    bodes_out = []
    search_list = data_set2.trunc_avebodes
    for bode in data_set1.trunc_avebodes:
        curbode = rwkbode.rwkbode()
        bode2 = find_matching_bode(bode, search_list)
        for key in bode_copy_keys:
            val = getattr(bode, key)
            setattr(curbode, key, val)
        for key in bode_merge_keys:
            vect1 = getattr(bode, key)
            mat1 = column_stack([vect1])
            vect2 = getattr(bode2, key)
            mat2 = column_stack([vect2])
            mat_out = row_stack([mat1, mat2])
            setattr(curbode, key, squeeze(mat_out))
        bodes_out.append(curbode)
    my_data_set.avebodes = bodes_out
    my_data_set.f = f
    my_data_set.filepaths = data_set1.filepaths+data_set2.filepaths
    return my_data_set
        

