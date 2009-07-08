#!/usr/bin/env python
from pylab import *
from scipy import *

import os, pdb#, re

import txt_data_processing
reload(txt_data_processing)

import txt_mixin

from txt_data_processing import Bode_Options, Bode_Data_Set
#import rst_creator
#reload(rst_creator)
#label_re = re.compile('(\\w+) *?: *?(.*)')

def _clean_comments(listin):
    listout = txt_mixin.txt_list([])
    for line in listin:
        ind = line.find('#')
        if ind == 0:
            line = ''
        elif ind > 0:
            line = line[0:ind]
        listout.append(line)
    return listout

def _clean_blanks(listin):
    listout = filter(None, listin)
    return listout

def _clean(listin):
    temp = _clean_comments(listin)
    listout = _clean_blanks(temp)
    return listout

def get_data_set_options(listin, key_list=None):
    if key_list is None:
        key_list = ['pat','truncstart','truncend', 'output']
    N = len(listin)
    i = 0
    my_dict = {}
    while i < N:
        curline = listin[i]
        curline = curline.strip()
        curline = curline.replace(' =','=')
        curline = curline.replace('= ','=')
        found_any = 0
        for key in key_list:
            if curline.find(key+'=')==0:
                listin.pop(i)
                key_str, value_str = curline.split('=',1)
                try:
                    value = float(value_str)
                except ValueError:
                    value = value_str
                my_dict[key] = value
                found_any = 1
                break
        if found_any:
            N -= 1
        else:
            i += 1
    return my_dict

def _fix_continued_lines(cleanlist):
    """By the time you pass a list to this function, it should be
    clean, i.e. it should have no comments and no blank lines.  So,
    each line is either the start of a new bode definition or it is a
    continuatin of a previous line.  The start of a new definition
    begins with 'label:' where label is composed of alphanumeric
    characters and there may be spaces around the colon."""
    N = len(cleanlist)
    i = 0
    while i < N:
        line = cleanlist[i]
        if line.find(':') == -1:
            #this is a continuation line
            line = line.strip()
            prevline = cleanlist[i-1].rstrip()
            if prevline[-1] == '\\':
                prevline = prevline[0:-1]
                prevline = prevline.rstrip()
            if prevline[-1] != ',':
                prevline += ','
            cleanlist.pop(i)
            N -= 1
            cleanlist[i-1] = prevline + line
        else:
            i += 1
    return cleanlist
    
def load_bode_options(filepath):
    myfile = txt_mixin.txt_file_with_list(filepath)
    mylist = myfile.list
    mylist = _clean(mylist)
    ds_opts = get_data_set_options(mylist)
    mylist = _fix_continued_lines(mylist)
    return mylist, ds_opts

def _fix_bad_splits(listin):
    """Spliting the key=value pairs string at the commas causes
    problems with freqlim, phaselim, and so on where there is a list
    for the value, i.e. phaselim=[-200,200] has a comma in the middle
    of the value.  For now, the fix is to search each item after doing
    a str.split(',') and checking that each item has an '='.  Any item
    that doesn't should be merged with the previous item in listin:"""
    N = len(listin)
    i = 0
    while i < N:
        item = listin[i]
        if item.find('=') == -1:
            #this is a bad split victim
            listin.pop(i)
            previtem = listin[i-1]
            listin[i-1] = previtem + ','+ item
            N -= 1
        else:
            i += 1
    return listin
    
def _parse(str_in):
    """Take a string of key=value pairs separated by commas and turn
    it into a dictionary."""
    my_list = str_in.split(',')
    my_list = _fix_bad_splits(my_list)
    my_dict = {}
    for item in my_list:
        key, val_str = item.split('=')
        key = key.strip()
        val_str = val_str.strip()
        try:
            val = float(val_str)
        except ValueError:
            val = val_str
        exec_list = ['phaselim','freqlim','maglim','cohlim']
        if key in exec_list:
            exec('val='+val_str)
        my_dict[key] = val
    return my_dict
    
def parse_bode_options(listin):
    """Each line of listin should start with 'label:' followed by a
    list of key=value pairs separated by commas."""
    labels = []
    bode_opt_list = []
    for line in listin:
        label_str, rest = line.split(':',1)
        label_str = label_str.strip()
        try:
            label = int(label_str)
        except ValueError:
            label = label_str
        labels.append(label)
        opt_dict = _parse(rest)
        in_str = opt_dict.pop('input')
        out_str = opt_dict.pop('output')
        cur_opts = Bode_Options(in_str, out_str, **opt_dict)
        bode_opt_list.append(cur_opts)
    return dict(zip(labels, bode_opt_list))
    
def create_data_set(pat, bode_opt_dict):
    #since a Bode_Data_Set has a find_bode option and is built around
    #Bode lists rather than dictionaries, I am going to turn the dict
    #into a list.
    my_keys = opt_dict.keys()
    my_keys.sort()
    bode_list = []
    for key in my_keys:
        bode_list.append(opt_dict[key])

    data_set = Bode_Data_Set(pat, bode_list)
    data_set.Calc_Bodes()
    data_set.Calc_Ave_Bodes()

    return data_set


if __name__ == '__main__':
    from optparse import OptionParser

    usage = 'usage: %prog [options] opts_file.txt'
    parser = OptionParser(usage)

    parser.add_option("-s", "--truncstart", dest="truncstart", \
                      help="The frequency to use as the start of the data truncation.\n" + \
                      "Data will be kept from truncstart to truncend.", \
                      default=None, type="float")

    parser.add_option("-e", "--truncend", dest="truncend", \
                      help="The frequency to use as the end of the data truncation.\n" + \
                      "Data will be kept from truncstart to truncstop.", \
                      default=None, type="float")

    parser.add_option("-p", "--pat", dest="pat", \
                      help="glob pattern for txt data files to include in the data set", \
                      default="", type="str")

    parser.add_option("-g", "--plot", dest="plot", \
                      help="should the Bode plots be shown on a graph.",\
                      default=1, type="int")

    parser.add_option("-o", "--output", dest="output", \
                      help="Output name for saving the data.",\
                      default="", type="str")

    (options, args) = parser.parse_args()

    print('options='+str(options))
    print('args='+str(args))

    opts_file = args[0]
    folder, opts_name = os.path.split(opts_file)
    list_in, ds_opts = load_bode_options(opts_file)
    opt_dict = parse_bode_options(list_in)
    if options.pat:
        pat = options.pat
    else:
        pat = ds_opts['pat']
    fullpat = os.path.join(folder, pat)
    data_set = create_data_set(fullpat, opt_dict)
    if options.truncstart:
        ds_opts['truncstart'] = options.truncstart
    if options.truncend:
        ds_opts['truncend'] = options.truncend
    save_name = ''
    if options.output:
        save_name = options.output
    elif ds_opts.has_key('output'):
        save_name = ds_opts['output']
    if save_name:
        save_name, ext = os.path.splitext(save_name)
        save_path = os.path.join(folder, save_name)
        data_set.save_bodes_and_time_domain(save_path)
    
    call_show = 0
    if ds_opts.has_key('truncstart') and ds_opts.has_key('truncend'):
        data_set.Truncate_and_Review(ds_opts['truncstart'], ds_opts['truncend'])
        call_show = 1
    elif options.plot:
        data_set.Bode_Plot2()
        call_show = 1

    if call_show:
        show()
        
