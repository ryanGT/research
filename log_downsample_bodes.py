#!/usr/bin/env python
from pylab import *
from scipy import *

import os, pdb#, re

import txt_data_processing
#reload(txt_data_processing)

import txt_mixin

from calc_and_save_bodes import _clean, get_data_set_options

def load_and_parse_opts(pathin):
    myfile = txt_mixin.txt_file_with_list(pathin)
    mylist = myfile.list
    mylist = _clean(mylist)
    my_dict = get_data_set_options(mylist, key_list=['load_name', \
                                                     'save_name', \
                                                     'freqs', \
                                                     'points', \
                                                     'plot'])
    exec_list = ['freqs','points']
    for key in exec_list:
        val_str = my_dict[key]
        exec('val='+val_str)
        my_dict[key] = val
    return my_dict


if __name__ == '__main__':
    from optparse import OptionParser

    usage = 'usage: %prog [options] log_opts_file.txt'
    parser = OptionParser(usage)

    parser.add_option("-g", "--plot", dest="plot", \
                      help="should the Bode plots be shown on a graph.",\
                      default=1, type="int")

##     parser.add_option("-o", "--output", dest="output", \
##                       help="Output name for saving the data.",\
##                       default="", type="str")

    (options, args) = parser.parse_args()

    print('options='+str(options))
    print('args='+str(args))

    opts_file = args[0]
    folder, opts_name = os.path.split(opts_file)
    log_opt_dict = load_and_parse_opts(opts_file)
    load_path = os.path.join(folder, log_opt_dict['load_name'])
    data_set = txt_data_processing.load_avebode_data_set(load_path)
    freqs = log_opt_dict['freqs']
    points = log_opt_dict['points']
    inds, mask, flog = data_set.Log_Compress_Data(freqs, points)
    outpath = os.path.join(folder, log_opt_dict['save_name'])
    data_set.save_as(outpath)
    
    if options.plot:
        figs = data_set.Bode_Plot2()

        data_set.Bode_Plot2(attr='compressed_avebodes', \
                    f_attr='compressed_f', figs=figs, \
                    clear=False, linetype='o')

        show()
        
