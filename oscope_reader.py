from numpy import *
import glob, os
import txt_mixin

class oscope_folder(object):
    """A class for easily loading data from my oscope (Tektronix).
    The data is assumed to be in .csv files in a certain folder."""
    def __init__(self, folder_path, col_map={}, delimiter=',', skiprows=18):
        self.folder_path = folder_path
        self.delimiter = delimiter
        self.skiprows = skiprows


    def find_csvs(self):
        pat1 = os.path.join(self.folder_path, '*.CSV')
        pat2 = os.path.join(self.folder_path, '*.csv')
        files1 = glob.glob(pat1)
        files2 = glob.glob(pat2)
        files3 = [item for item in files2 if item not in files1]
        self.csv_files = files1+files3


    def _load_raw(self):
        raw_data = []
        for path in self.csv_files:
            temp = txt_mixin.txt_file_with_list(path)
            raw_data.append(temp.list)
        self.raw_data = raw_data
        return raw_data


    def convert_row(self, rowin):
        curlist = rowin.split(self.delimiter)
        while not curlist[-1]:
            curlist.pop(-1)
        t_i = curlist[-2]
        v_i = curlist[-1]
        return t_i, v_i


    def zero_out_t(self):
        for t in self.tlist:
            t0 = t[0]
            t -= t0

    def shift_t(self, shift):
        for t in self.tlist:
            t += shift

    def raw_data_to_floats(self):
        if not hasattr(self, 'raw_data'):
            self._load_raw()
        self.str_list = []
        self.tlist = []
        self.vlist = []
        for col in self.raw_data:
            cur_list = col[self.skiprows:]
            self.str_list.append(cur_list)
            N = len(cur_list)
            t = zeros(N)
            v = zeros(N)
            for i, row in enumerate(cur_list):
                t_i, v_i = self.convert_row(row)
                t[i] = float(t_i)
                v[i] = float(v_i)
            self.tlist.append(t)
            self.vlist.append(v)


    def Plot(self, fig=None, fignum=1, clear=True, \
             ylabel='Voltage', xlabel=None, use_ms=0):
        if xlabel is None:
            if use_ms:
                xlabel = 'Time (ms)'
            else:
                xlabel = 'Time (sec)'

        if fig is None:
            import pylab as PL
            fig = PL.figure(fignum)
        if clear:
            fig.clf()
        ax = fig.add_subplot(111)
        for t, v in zip(self.tlist, self.vlist):
            if use_ms:
                t*=1000.0
            ax.plot(t,v)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        
