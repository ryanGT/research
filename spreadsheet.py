"""spreadsheet module spreadsheet.py

Basic usage example:

colmap = {'Student':'student','Team #':'team'}
mysheet = spreadsheet.SpreadsheetFromPath('class_list_w_group_nums.csv', \
                                          colmap=colmap)

mysheet.FindLabelRow('Student')#this is the upper left column label
#note you could  also pass in a list
mysheet.FindDataColumns()
mysheet.MapCols()


after running the above code, mysheet will have attributes student and team.

"""

import csv, pdb, time, os
from scipy import array, column_stack, integrate, r_, fft, \
     arctan2, pi, shape, log10, squeeze, imag, real, signal, \
     all, fromstring, zeros, io, dtype, row_stack, arange, \
     atleast_2d

import re

import numpy

from numpy import float64, int32

from IPython.core.debugger import Pdb

import xlrd, copy

import cPickle#, dumb_shelve
from scipy.io import dumb_shelve

import rwkmisc, mplutil

import txt_mixin
reload(txt_mixin)

from rwkdataproc import thresh, CalcSpectra, makefreqvect

import DataProcMixins
reload(DataProcMixins)

def col_from_nested(nested, col):
    return [row[col] for row in nested]

def myfloat(stringin):
    """Try converting stringin to a float.  Return 0 for empty
    strings.  Convert the '#DIV/0!' to 0.  On all other float
    conversions failures, just return stringin."""
    if stringin:
        try:
            out = float(stringin)
        except ValueError:
            #apparently, this is supposed to fail gracefully for string data
            if stringin == '#DIV/0!':
                out = 0.0
            else:
                out = stringin
        return out
    else:
        return 0.0


def split_names(namelist, delim=' ', reverse=False):
    lastnames = []
    firstnames = []
    for curname in namelist:
        first, last = curname.split(delim, 1)
        if reverse:
            temp = first
            first = last
            last = temp
        last = last.strip()
        if last[-1] == ',':
            last = last[0:-1]
        first = first.strip()
        if first.find(' ') > -1:
            first, middle = first.split(' ',1)
            first = first.strip()
        lastnames.append(last)
        firstnames.append(first)
    return firstnames, lastnames


class tabdelim(csv.Dialect):
    # placeholders
    delimiter = '\t'
    quotechar = '"'
    doublequote = False
    skipinitialspace = False
    lineterminator = '\r\n'
    quoting = csv.QUOTE_MINIMAL

csv.register_dialect('tabdelim', tabdelim)

class mycsv(csv.Dialect):
    # placeholders
    delimiter = ','
    quotechar = '"'
    doublequote = True
    skipinitialspace = False
    lineterminator = '\r\n'
    quoting = csv.QUOTE_MINIMAL

csv.register_dialect('mycsv', mycsv)

class lsdcsv(csv.Dialect):
    # placeholders
    delimiter = ','
    quotechar = '"'
    doublequote = True
    skipinitialspace = False
    lineterminator = ',\r\n'
    quoting = csv.QUOTE_MINIMAL

csv.register_dialect('lsdcsv', lsdcsv)

def LoadFromShelf(pathin):
    dir,filename = os.path.split(pathin)
    filebase = filename.split('.')[0]
    fn = os.path.join(dir, filebase)
    f = dumb_shelve.open(fn, "r")
    return dict(f)

def ObjectFromShelf(pathin):
    mydict = LoadFromShelf(pathin)
    return rwkmisc.dictobject(**mydict)

def all(listin):
    fitems = [item for item in listin if not item]
    return not fitems


def get_lineterminator(linein):
    lt=''
    while linein[-1]=='\r' or linein[-1]=='\n':
        lt=linein[-1]+lt
        linein = linein[0:-1]
        if not linein:
            break
    return lt


def clean_line(linein):
#    lt=''
    while linein[-1]=='\r' or linein[-1]=='\n':
#        lt=linein[-1]+lt
        linein = linein[0:-1]
        if not linein:
            break
    return linein


def determine_lt(linesin, dialectin):
        """Verify the line terminator of a dialect from my sniff function."""
        lts=[get_lineterminator(line) for line in linesin]
        tests=[item == lts[0] for item in lts]
        if all(tests):
            dialectin.lineterminator = lts[0]
        return dialectin

def open_workbook(filepath):
    wb=xlrd.open_workbook(filepath)
    return wb

def get_workbook_sheet(filepath, sheet=0):
    wb=open_workbook(filepath)
    sh1=wb.sheet_by_index(sheet)
    return sh1

def mycomp(strin, pattern, exact=False):
    """check pattern against strin.  If exact, then retrun
    pattern==strin, else return strin.find(pattern)>-1."""
    if exact:
        return bool(strin==pattern)
    else:
        return bool(strin.find(pattern)>-1)

def transpose(m):
    """transpose(m): transposes a 2D matrix, made of tuples or lists of tuples or lists,
    keeping their type.

    >>> transpose([])
    Traceback (most recent call last):
      ...
    IndexError: list index out of range
    >>> transpose([[]])
    []
    >>> transpose([1,2,3])
    Traceback (most recent call last):
      ...
    TypeError: zip argument #1 must support iteration
    >>> transpose([[1,2,3]])
    [[1], [2], [3]]
    >>> transpose( [[2, 2, 2], [2, 2, 2]] )
    [[2, 2], [2, 2], [2, 2]]
    >>> transpose( [(2, 2, 2), (2, 2, 2)] )
    [(2, 2), (2, 2), (2, 2)]
    >>> transpose( ([2, 2, 2], [2, 2, 2]) )
    ([2, 2], [2, 2], [2, 2])
    >>> transpose( ((2, 2, 2), (2, 2, 2)) )
    ((2, 2), (2, 2), (2, 2))
    >>> t = [[[1], [2]], [[3], [4]], [[5], [6]]]
    >>> transpose(t)
    [[[1], [3], [5]], [[2], [4], [6]]]
    """
    if isinstance(m, list):
        if isinstance(m[0], list):
            return map(list, zip(*m))
        else:
            return zip(*m) # faster
    else:
        if isinstance(m[0], list):
            return tuple(map(list, zip(*m)))
        else:
            return tuple( zip(*m) )

def get_col(nestedlist, index):
    return transpose(nestedlist)[index]

def interp(x, x1, x2, y1, y2, eps=1e-16):
    """Return a linear interpolation of the y that corresponds to x
    where x is between x1 and x2.  The end points of the line segment
    which is being interpolated is (x1,y1) and (x2,y2)."""
    assert x1-eps <= x <= x2+eps, "The point you are interpolating must be between x1 and x2."
    if x1-eps <= x <= x1+eps:
        return y1
    elif x2-eps <= x <= x2+eps:
        return y2
    else:
        return y1 + ((x - x1) * (y2 - y1))/(x2 - x1)


class SpreadSheet(object):
    def __init__(self, pathin=None, skiprows=0, \
                 collabels=None, colmap=None, \
                 datafunc=float, picklekeys=[], \
                 labelsin=None, datain=None):
#        print('pathin='+pathin)
        self.path = pathin
        self.labelrow = -1
        self.datacolumns = None
        self.header = []
        self.alldata = []
        self.labels = []
        self.data = []
        self.tunits = 'sec'
        self.skiprows = skiprows
        self.picklekeys = picklekeys
        self.colmap = colmap
        self.collabels = collabels
        self.datafunc = datafunc
        if (labelsin is not None) and (datain is not None):
            self._init_from_data_and_labels(labelsin, datain)


    def _init_from_data_and_labels(self, labelsin, datain):
        self.labels = labelsin
        self.alldata = datain
        self.data = copy.copy(datain)
        self.labelrow = 0
        if self.colmap is not None:
            self.MapCols()


    def FindColLabel(self, label):
        """Search self.labels for label and return the index if found.
        Return -1 if not found."""
        ind = -1
        for n, curlabel in enumerate(self.labels):
            if curlabel == label:
                return n
        return ind


    def PopList(self, indlist):
        if type(indlist)== numpy.ndarray:
            indlist = indlist.tolist()
        indlist.sort(reverse=True)
        return [self.Pop(ind) for ind in indlist]


    def Pop(self, ind):
        """Remove the row with index ind from self.data"""
        if type(self.data)==list:
            return self.data.pop(ind)
        else:
            poprow = self.data[ind]
            newmat = numpy.delete(self.data, ind, axis=0)
            self.data = newmat
            return poprow


    def MapOut(self, outpath, colmap, labels=None, dtype='|S100'):
        """Output selected attributes to a spreadsheet file using
        colmap to map attributes to specific columns of the output
        file.  The keys in colmap are the column labels and the values
        are attributes that will be retrieved via getattr.  These
        should return arrays or lists of the same length."""
        if labels is None:
            labels = colmap.keys()
        NC = len(labels)
        firstcol = getattr(self, colmap[labels[0]])
        NR = len(firstcol)
        data = zeros((NR,NC), dtype=dtype)
        for i, label in enumerate(labels):
            curcol = array(getattr(self, colmap[label]))
            data[:,i] = curcol.astype(dtype)
        self.WriteCSV(outpath, labels, data)


    def MapCols(self, colmap = None, minrows=-1):
        if colmap is None:
            colmap = self.colmap
        else:
            self.colmap = colmap
        assert colmap is not None, "Cannot map columns without a colmap - either pass one in, or make sure self.colmap is set."
        if self.labelrow == -1:
            self.FindLabelRow()
        if not self.datacolumns:
            self.FindDataColumns(self.collabels)
        if self.alldata == []:
            self.ReadData(minrows=minrows)
        if self.data == []:
            self.ReadDataColumns()
        for key, val in colmap.iteritems():
            curcol = self.ReadDataColumn(key, parsefunc=self.datafunc)#, exact=exact, removeempty=removeempty)
            if key.lower().find('time')==0 and key.lower().find('(ms)')>-1:
                print('converting from ms')
                curcol=curcol/1000.0#convert from ms
            setattr(self, val, curcol)
        self.keys = colmap.values()


    def GetNull(self):
        try:
            if isinstance((self.data[0,0]), str):
                null=''
            elif isinstance((self.data[0,0]), float):
                null=0.0
            else:
                null=None
        except:
            null=''
        return null


    def GetDataFromRow(self, rowin, labels):
        """Passing in a row from self.data and desired labels, extract
        the corresponding elements from the row and return those
        elements as a list."""
        inds = [self.collabels.index(item) for item in labels]
        return [rowin[ind] for ind in inds]


    def Get_Dict_from_Row(self, rowind, labels):
        """Passing in a row index and desired labels, extract
        the corresponding elements from the row and return those
        elements as a dict with labels as the keys."""
        row = self.alldata[rowind]
        mylist = self.GetDataFromRow(row, labels)
        return dict(zip(labels, mylist))


    def GetDataFromRowIndex(self, index, labels):
        myrow = self.data[index]
        return self.GetDataFromRow(myrow, labels)

    def _search_vect(self, value, myvect):
        if type(myvect)==numpy.ndarray:
            return myvect == value
        else:
            boolist = [item==value for item in myvect]
            return array(boolist)


    def SearchCol(self, value, collabel):
        """Return a boolean vector corresponding to whether or not
        self[collabel] == value."""
        myvect = self[collabel]
        return self._search_vect(value, myvect)


    def SearchAttr(self, value, attr):
        myvect = getattr(self, attr)
        return self._search_vect(value, myvect)


    def GetMatches(self, vallist, collist):
        """Return a boolean vector corresponding to rows of self.data
        that match vallist.  The items of collist will be used to
        determine which columns of self.data should be searched."""
        first = 1
        assert vallist!=[] and collist!=[], "You cannot call GetMatches with vallist or collist empty:\n vallist="+str(vallist)+"\n collist="+str(collist)
        for curval, curcol in zip(vallist, collist):
            if first:
                boolvect = self.SearchCol(curval, curcol)
            else:
                curbool = self.SearchCol(curval, curcol)
                boolvect *= curbool
            first = 0
            if not boolvect.any():#exit as soon as the first test returns all falses - everything will stay false after that
                return boolvect
        return boolvect


    def GetMatchInds(self, vallist, collist):
        boolvect = self.GetMatches(vallist, collist)
        indvect = arange(self.data.shape[0])
        matchinds = indvect[boolvect]
        return matchinds


    def FilterMatches(self, vallist, collist):
        """Call self.GetMatches to find the boolean vector for vallist
        and collist and then use that boolvect to filter self.data.

        Returns self.data[boolvect]"""
        boolvect = self.GetMatches(vallist, collist)
        return self.data[boolvect]


    def SpreadSheetFromData(self, datain):
        """This function would be used in conjunction with
        self.FilterMatches to return a Spreadsheet object that
        contains only the fitlered data."""
        sheetout = copy.deepcopy(self)
        sheetout.data = datain
        return sheetout


    def SearchForMatch(self, vallist, collist):
        """Determine whether or not the items in vallist correspond to
        a row already in self.data.  The items of collist will be used
        to determine which columns of self.data should be searched."""
        boolvect = self.GetMatches(vallist, collist)
        return boolvect.any()


    def BuildIndList(self, collabels):
        inds = []
        for item in collabels:
            try:
                curind = self.collabels.index(item)
                inds.append(self.datacolumns[curind])
            except ValueError:
                inds.append(-1)
        return inds


    def Search_Down_Col(self, colindex, value):
        """Search down one column of self.alldata to find value.
        Return the index of the first == match, or -1 if no match
        occurs."""
        collist = self.get_col(colindex)
        for n, item in enumerate(collist):
            if item == value:
                return n
        return -1


    def AppendRow(self, row):
        assert len(row)==len(self.datacolumns)
        newrow = [row[ind] for ind in self.datacolumns]
        self.data = row_stack([self.data, newrow])


    def AppendSpreadSheet(self, othersheet, checkcols, defaults={}):
        """Append the rows of othersheet to self.  A row of othersheet
        whose values corresponding to checkcols match a row already in
        self will be considered a duplicated and not added.

        defaults is a dictionary whose keys correspond to collabels
        and whose values are the corresponding row entries.

        othersheet must already be opened and have collabels and data
        which are not empty.

        The columns of othersheet corresponding to self.collabels are
        the ones that will be copied and they will be copied into
        self.data."""
        msg = 'othersheet must already be opened and have collabels and data which are not empty'
        assert othersheet.collabels is not None, msg
        assert othersheet.data is not None, msg
        assert len(othersheet.collabels) > 0, msg
        assert len(othersheet.data) > 0, msg
        mynull = self.GetNull()
        mylabels = self.collabels
        inds = othersheet.BuildIndList(mylabels)
        checkinds = othersheet.BuildIndList(checkcols)
        for row in othersheet.data:
            vallist = BuildListforMatching(row, checkinds, checkcols, defaults=defaults, null=mynull)
            curmatch = self.SearchForMatch(vallist, checkcols)
            if not curmatch:#then the current row of othersheet is NOT already in self.data
                vallist = BuildListforMatching(row, inds, mylabels, defaults=defaults, null=mynull)
                self.AppendRow(vallist)



    def Filter(self, boolvect):
        """Return self.data[boolvect] where boolvect is the output of
        boolfunc.

        So, the return value is a numpy array."""
        return self.data[boolvect]


    def SpreadSheetFromData(self, datain):
        """This function would be used in conjunction with self.Filter
        to return a Spreadsheet object that contains only the fitlered
        data."""
        sheetout = copy.deepcopy(self)
        sheetout.data = datain
        return sheetout


    def __getitem__(self, index):
        if type(index)==str:
            return self.ReadDataColumn(index)
        elif type(index)==int:
            return self.get_col(index)
        else:
            print('SpreadSheet.__getitem__ cannot understand index of '+str(index)+', expecting a string or integer.')


##     def __setitem__(self, index):
##         if type(index)==str:
##             return self.ReadDataColumn(index)
##         elif type(index)==int:
##             return self.get_col(index)
##         else:
##             print('SpreadSheet.__getitem__ cannot understand index of '+str(index)+', expecting a string or integer.')


    def _BuildDict(self, picklekeys=[]):
        if not picklekeys:
            picklekeys = self.picklekeys
        assert len(picklekeys) > 0, "_BuildDict called with no picklekeys on a SpreadSheet instance with no self.picklekeys."
        mydict={}
        for key in picklekeys:
            curitem = getattr(self, key)
            mydict[key] = curitem
        return mydict


    def Pickle(self, pklpath=None, picklekeys=[]):
        if pklpath is None:
            mypath, myext = os.path.splitext(self.path)
            pklpath = mypath+'.pkl'
        mypkl = open(pklpath,'wb')
        mydict = self._BuildDict(picklekeys)
        cPickle.dump(mydict, mypkl, protocol=2)
        mypkl.close()
        self.pklpath = pklpath
        return pklpath


    def Shelve(self, shelfpath=None, picklekeys=[]):
        if shelfpath is None:
            shelfpath, myext = os.path.splitext(self.path)
        mydict = self._BuildDict(picklekeys)
        io.save(shelfpath, mydict)
        self.shelfpath = shelfpath


    def DropEmptyRows(self, cols=[0,1]):
        """Delete rows from self.alldata whose cols have a False
        boolean value."""
        n=0
        while n <= (len(self.alldata)-1):
            empty = True
            currow = self.alldata[n]
            if currow:
                for ind in cols:
                    if currow[ind]:
                        empty = False
                        break
            if empty:
                self.alldata.pop(n)
            else:
                n+=1

    def ReadRows(self, maxrows=None, startrow=0, parsefunc=None, \
                 stop_on_blank=True, minrows=-1):
        i=0
        dataout = []
        first = 1
#        floattime = 0.0
#        appendtime = 0.0
        for row in self.iterrows():
            if i>=startrow:
                if stop_on_blank and ((i-startrow) > minrows) and (row[0] == ''):
                    print('in SpreadSheet.ReadRows, i=%s, startrow=%i, minrows=%s' % (i, startrow, minrows))
                    print('stopping based on blank cell in first column, i=%i' % i)
                    break

                if first and (parsefunc is None):#try float anyways if parsefunc is None
                    try:
                        temprow = map(float, row)
                        parsefunc = float
                    except ValueError:
                        parsefunc = None
                first = 0
                if parsefunc is not None:
#                    t0 = time.time()
                    temprow = map(parsefunc, row)
#                    t1 = time.time()
                    dataout.append(temprow)
#                    t2 = time.time()
#                    floattime+=t1-t0
#                    appendtime+=t2-t1
                else:
                    dataout.append(row)
            if maxrows is not None:
                if i >= maxrows+startrow-1:
                    break
            i+=1
#        print('ReadRows: floattime='+str(floattime))
#        print('ReadRows: appendtime='+str(appendtime))
        if type(dataout[0][0])==float:
            return array(dataout)
        else:
            return dataout


    def ReadHeader(self, headerrows=100):
        header = self.ReadRows(maxrows=headerrows)
        self.header = header


    def GetLabelRow(self):
        if self.labelrow < 0:
            raise ValueError, "self.labelrow must be defined.  You must call FindLabelRow before you call GetLabelRow."
        if not self.header:
            self.ReadHeader(headerrows=self.labelrow+1)
        self.labels = self.header[self.labelrow]
        self.strip_labels()
        return self.labels


    def GetHeaderInfo(self, searchlist, searchcol = 0, usere = False, startrow = 0, exact=False):
        """Search for info that is stored in the first lines of the
           spreadsheet.  The labels are assumed to be stored in
           searchcol and the values are assumed to be stored in the
           next column of the same row.  For example, if the file were
           CSV, the format would be:

           Date, 12/25/06
           Test Description, "Blah, Blah, Blah"
           next label, next value
           ...
           """
        valuesout = []
        labellist = []
        for item in searchlist:
            foundit = False
            for row in self.iterrows():
                if len(row)>searchcol:
                    curlabel = row[searchcol]
                    if mycomp(curlabel,item,exact):
                        labellist.append(curlabel)
                        valuesout.append(row[searchcol+1])
                        foundit = True
                        break
            if not foundit:
                labellist.append('Did not find:'+item)
                valuesout.append(None)
        return valuesout, labellist


    def FindLabelRow(self, labels = None, startcol=0, exact=False, removeemptylabels=True):
        """Search for a row of column labels, where labels is a list
        of the labels to search for (automatically turned into a list
        of one item if passed a string).

        The algorithm searches down column startcol for the first item
        in labels and then continues searching across the row for the
        other items in labels until one fails to be found.

        If all items in labels are found in order in one row, that is
        the labelrow and the index of that row is returned."""
        if labels is None:
            labels = self.collabels
        if not (type(labels)==list or type(labels)==tuple):
            labels=[labels]
        labelrow = -1
        for ind, row in enumerate(self.iterrows()):
            if len(row)>startcol:
                foundall = True
                for n, item in enumerate(labels):
                    if not mycomp(row[startcol+n],item,exact):
                        foundall = False
                        break
                if foundall:
                    labelrow = ind
                    self.labelrow = labelrow
                    if removeemptylabels:
                        self.labels = [item for item in row if item]
                    else:
                        self.labels = row
                    return labelrow
        return labelrow


    def FindDataColumns(self, collabels=None, exact=False):
        """Search self.labels for collabels and return the indices as
        a list.  self.labels contains all the labels in self.labelrow.
        collabels is a list of the column labels from which data is to
        be extracted."""
        if not self.labels:
            raise ValueError, "self.labels must be defined.  You must call FindLabelRow before you call FindDataColumns."
        colnums = []
        if not collabels:
            assert self.colmap, "You cannot call FindDataColumns without either passing collabels or having self.colmap already defined."
            collabels = [item for item in self.labels if item in self.colmap.keys()]
        for item in collabels:
            foundit = False
            for n, curlabel in enumerate(self.labels):
                if mycomp(curlabel, item, exact):
                    colnums.append(n)
                    foundit = True
                    break
            if not foundit:
                print('Cound not find column with label '+item+'.')
        self.datacolumns = colnums
        self.collabels = collabels
        return colnums


    def ReadData(self, skiprows=None, maxrows=None, parsefunc=None, **kwargs):
        """Read in all the data below labelrow starting at row
        labelrow+skiprows+1.  Stop after maxrows or read all if
        maxrows is None.  If parsefunc is not None, map parsefunc onto
        each row."""
        if skiprows is None:
            if hasattr(self, 'skiprows'):
                skiprows = self.skiprows
            else:
                skiprows = 0
        data = self.ReadRows(startrow=(self.labelrow+1+skiprows), \
                             maxrows=maxrows, parsefunc=parsefunc, \
                             **kwargs)
        self.alldata = data


    def get_col_old(self, ind):
        if self.alldata == []:
            self.ReadData()
        return transpose(self.alldata)[ind]

    def get_col(self, ind):
        col_out = None
        for row in self.alldata:
            ent = row[ind]
            if col_out is None:
                col_out = [ent]
            else:
                col_out.append(ent)
        return col_out

    def get_data_col(self, ind):
        if hasattr(self.data, 'transpose'):
            if callable(self.data.transpose):
                return self.data.transpose()[ind]
            else:
                return self.data.transpose[ind]
        else:
            return transpose(self.data)[ind]


    def GetColFromArrayorNestedList(self, label, arrayorlist, labellist, exact=False):
        ind = rwkmisc.searchlist(labellist, label, exact=exact)
        assert ind > -1, label +" not found in self.labels."
        if type(arrayorlist)==type(array([1.5,2.5])):
            curcol = arrayorlist[:,ind]
        else:
            curcol = col_from_nested(arrayorlist, ind)
        return curcol


    def GetColFromAllData(self, label, exact=False):
        return self.GetColFromArrayorNestedList(label, self.alldata, self.labels, exact=exact)


    def GetColFromData(self, label, exact=False):
        return self.GetColFromArrayorNestedList(label, self.data, self.collabels, exact=exact)


    def ReadDataColumn(self, label, skiprows=None, maxrows=None, \
                       parsefunc=None, exact=False, \
                       removeempty=False, **kwargs):
        """Read in one column of data with the label label and return
        the data, passing it through parsefunc if it is not None."""
        if self.alldata==[]:
            print('reading data')
            self.ReadData(skiprows=skiprows, **kwargs)
        assert self.labels, "You must call self.FindLabelRow before trying to read a data column."
        ind = -1
        #search first in collabels and return the appropriate column of
        #self.data, if that fails, search in self.labels and return
        #the appropriate column of self.alldata
        if self.collabels:
            if len(self.collabels)>0 and len(self.data)>0:
                return self.GetColFromData(label, exact=exact)
        return self.GetColFromAllData(label, exact=exact)


    def SetDataColumn(self, label, listin):
        myarray = array(listin)
        myindex = self.collabels.index(label)
        self.data[:,myindex] = myarray


    def ReadRecArray(self, parsefuncs, colmap=None):
        """After calling FindDataColumns and ReadDataColumns with
        parsefunc=None, call this function to create a numpy recarray
        from self.data.

        colmap is a dictionary that maps self.collabels to the labels
        desired for the recarray columns, similar to the colmap used
        else where in this module."""
        assert hasattr(self, 'data'), "You must call self.ReadDataColumns(parsefunc=None) before calling self.ReadRecArray"
        assert len(parsefuncs)==len(self.data[0]), "len(parsefuncs) does not match len(self.data[0]) - i.e. the number of columns in self.data"
        deftype = self.data.dtype#the default type will be the type used by self.data, which is most likely a string of the size needed to contain the largest element of self.data
        typelist = []
        for curfunc in parsefuncs:
            if curfunc == int:
                typelist.append(dtype(int32))
            elif curfunc == float:
                typelist.append(dtype(float64))
            else:
                typelist.append(deftype)
        if colmap is not None:
            mylabels = [colmap[curlabel] for curlabel in self.collabels]
        else:
            mylabels = self.collabels
        mytups = [(curlabel, curtype) for curlabel, curtype in zip(mylabels, typelist)]
        dt = dtype(mytups)
        biglist = []
        for row in self.data:
            rowout = []
            for elem, curfunc in zip(row, parsefuncs):
                if curfunc:
                    curelem = curfunc(elem)
                    rowout.append(curelem)
                else:
                    rowout.append(elem)
            biglist.append(tuple(rowout))
        myarray = array(biglist, dtype=dt)
        self.recarray = myarray
        self.parsefuncs = parsefuncs
        self.colmap = colmap
        return myarray

    def ReadDataColumns(self, skiprows=None, maxrows=None, parsefunc=myfloat):
        """Similar to ReadData, but reads only the columns in
        self.datacolumns and attempts to parse the columns using
        parsefunc which defaults to float.  Assigns an array of the
        columns to self.data."""
        N = len(self.datacolumns)
        if type(parsefunc)!=list and type(parsefunc)!=tuple:
            parsefuncs = N*[parsefunc]
        else:
            parsefuncs = parsefunc

        if type(self.alldata)==type(array([1.5,2.5])):
            self.data = self.alldata[:,self.datacolumns]
        else:
            data = []
            if not self.alldata:
#                t0 = time.time()
                self.ReadData(skiprows=skiprows)
#                t1 = time.time()
#                print('ReadDataColumns ReadData time = '+str(t1-t0))
            for ind, curfunc in zip(self.datacolumns, parsefuncs):
                curcol = self.get_col(ind)
                if curfunc is not None:
                    try:
                        curcol = map(curfunc, curcol)
                    except ValueError:
                        print('ValueError in ReadDataColumns')
                        curcol = curcol
                data.append(curcol)
            #print('type(data[0][0])='+str(type(data[0][0])))
            myarray = column_stack(data)
            #print('type(myarray[0,0])='+str(type(myarray[0,0])))
            self.data = myarray


    def _Plot_prop_vs_t(self, prop, fig, clear=True, tprop='t'):
        """use getattr on self to get prop and tprop and then plot
        prop vs. tprop."""
        if clear:
            fig.clf()
        ax=fig.add_subplot(111)
        vect = getattr(self,prop)
        myt = getattr(self, tprop)
        ax.plot(myt, vect)
        ax.set_xlabel('Time (' + self.tunits+')')
        return ax


    def WriteAllDataCSV(self, outpath, dialect=None, append=False):
        if dialect is None:
            if hasattr(self,'dialect'):
                dialect=self.dialect
            else:
                dialect=mycsv
        if append:
            f=open(outpath,'ab')
        else:
            f=open(outpath, 'wb')
        writer = csv.writer(f,dialect)
        if self.labels:
            writer.writerow(self.labels)
        writer.writerows(self.alldata)
        f.close()


    def WriteCSV(self, outpath, labels, data, \
                 dialect=None, append=False):
        if dialect is None:
            if hasattr(self,'dialect'):
                dialect=self.dialect
            else:
                dialect=mycsv
        if append:
            f=open(outpath,'ab')
        else:
            f=open(outpath, 'wb')
        writer = csv.writer(f,dialect)
        writer.writerow(labels)
        writer.writerows(data)
        f.close()


    def WriteDataCSV(self, outpath, dialect=None, append=False):
        data = self.data
        labels = self.collabels
        self.WriteCSV(outpath, labels, data, dialect=dialect, \
                      append=append)


##     def WriteDataCSV(self, outpath, dialect=None, append=False):
##         if dialect is None:
##             if hasattr(self,'dialect'):
##                 dialect=self.dialect
##             else:
##                 dialect=mycsv
##         if append:
##             f=open(outpath,'ab')
##         else:
##             f=open(outpath, 'wb')
##         writer = csv.writer(f,dialect)
##         writer.writerow(self.collabels)
##         writer.writerows(self.data)
##         f.close()

    def InsertColFromList(self, matchlabel, matchlist, destlabel, valuelist):
        """This method is for copying a column from one spreadsheet to
        another in an intelligent way.  I need this a lot when I have
        a spreadsheet for calculating a project or final exam grade
        and I need to copy the final score into an overall spreadsheet
        for the class, but I cannot be sure that the names are exactly
        in the same order for various reasons (students dropping the
        course for example).

        The method searches down the column with label matchlabel for
        the items in matchlist and assuming it finds a match, it
        copies the corresponding element from valuelist into the
        column with label destlabel.

        For example, assuming matchlist was a list of student names,
        matchlabel would be the column label for names in the main
        spreadsheet - probably 'Name'.  If valuelist contained a list
        of scores for project 1, destlabel might be 'Project 1'.  Note
        that destlabel must exist in the labelrow of self."""
        destcol = self.labels.index(destlabel)
        searchcol = self.labels.index(matchlabel)

        for row in self.alldata:
            cursearch = row[searchcol]
            if cursearch:
                curind = search_list(matchlist, cursearch)
                if curind==-1:
                    print('could not find '+cursearch)
                else:
                    if valuelist[curind]:
                        if destcol>len(row)-1:
                            row.append(valuelist[curind])
                        else:
                            row[destcol] = valuelist[curind]


    def AppendColFromList(self, matchlabel, matchlist, destlabel, valuelist):
        self.labels.append(destlabel)
        self.InsertColFromList(matchlabel, matchlist, destlabel, valuelist)


    def AppendCol(self, collabel, listin=None, parsefunc=None):
        if collabel not in self.labels:
            self.labels.append(collabel)
        self.collabels.append(collabel)
        if listin is None:
            mydtype = self.data.dtype
            N = self.data.shape[0]
            newcol = zeros((N,), dtype=mydtype)
            listin = newcol
            #self.data = column_stack([self.data, newcol])

        if (self.data is None) or (len(self.data)==0):
            self.data = column_stack([listin])
        else:
            self.data = column_stack([self.data, listin])
        if hasattr(self, 'parsefuncs'):
            self.parsefuncs.append(parsefunc)
        self.FindDataColumns(self.collabels)


    def InsertCol(self, index, collabel, listin=None, parsefunc=None):
        self.collabels.insert(index,collabel)
        colsbefore = self.data[:,0:index]
        colsafter = self.data[:,index:]
        if listin is None:
            mydtype = self.data.dtype
            N = self.data.shape[0]
            listin = zeros((N,), dtype=mydtype)
        if index == 0:
            self.data = column_stack([listin, colsafter])
        else:
            self.data = column_stack([colsbefore, listin, colsafter])
        if hasattr(self, 'parsefuncs'):
            self.parsefuncs.insert(index, parsefunc)


    def SetDataCol(self, label, listin):
        index = rwkmisc.searchlist(self.collabels, label)
        assert index > -1, label +" not found in self.collabels."
        colsbefore = self.data[:,0:index]
        colsafter = self.data[:,index+1:]
        if index == 0:
            self.data = column_stack([listin, colsafter])
        else:
            self.data = column_stack([colsbefore, listin, colsafter])

    def strip_labels(self):
        self.labels = [item.strip() for item in self.labels]


## class DataSpreadSheet:
##     """This mix in class overrides the ReadData, ReadDataColumn, and ReadDataColumns methods for Spreadsheets that contain only numeric data."""

class spreadsheet_from_data(SpreadSheet):
    """This class seeks to make it possible to have spreadsheet
    functionality even if the data isn't read from a spreadsheet file.
    This might be usefull if data is pulled from a larger database
    file."""
    def __init__(self, labels, alldata, data=[], colmap={}, datafunc=float):
        self.labels = labels
        self.alldata = alldata
        self.colmap = colmap
        self.data = data
        self.labelrow = 0
        self.datacolumns = None
        self.datafunc = datafunc
        if colmap:
            self.collabels = [item for item in self.labels if item in self.colmap.keys()]



class CSVSpreadSheet(SpreadSheet):
    """Note that following the example of OpenOffice and the Python
    csv module, CSV does not have to mean comma seperated value, but
    means only delimitted text.  If the csv dialect is not specified,
    this class will attempt to sniff it out.  TabDelimSpreadSheet does
    specify the delimitter explicitly.  So, if you are not sure about
    your delimitter, use CSVSpreadSheet."""
    def __init__(self, pathin=None, dialect=None, skiprows=0, \
                 collabels=None, colmap=None, **kwargs):
        self.dialect = dialect
        SpreadSheet.__init__(self, pathin=pathin, skiprows=skiprows, \
                             colmap=colmap, collabels=collabels, **kwargs)


    def sniff(self, sniffbytes=1000, resniff=False):
        """Try and determine the dialect of a text file by reading in
        the first sniffbytes bytes of data from filepath.  My slightly
        more intelligent version than the csv module.

        If self.dialect is not None, this function will not attempt to
        re-sniff unless resniff is True"""
        if (self.dialect is not None) and not resniff:
            return
        else:
            f=open(self.path,'r')
            mylist=f.read(sniffbytes)
            mysniff=csv.Sniffer()
            cleanlines = [clean_line(item) for item in mylist]
            sniffstr=''.join(cleanlines)
            mydialect=mysniff.sniff(sniffstr)
            mydialect=determine_lt(mylist,mydialect)
            self.dialect = mydialect


    def iterrows(self):
        if self.dialect is None:
            self.sniff()
        reader = csv.reader(open(self.path,'rb'), self.dialect)
        for row in reader:
            yield row


    def ReadRows(self, maxrows=None, startrow=0, parsefunc=None, \
                 tryfast=True, **kwargs):
        """Attempt to read floating point text data files very
        quickly, if tryfast is True by calling ReadFloatRows."""
        if tryfast:
            return self.ReadFloatRows(startrow=startrow, \
                                      maxrows=maxrows, **kwargs)
        else:
            return SpreadSheet.ReadRows(self, startrow=startrow, \
                                        **kwargs)


    def ReadFloatRows(self, startrow=0, maxrows=None, **kwargs):
        """After some tinkering, this is the fastest approach I have
        come up with for reading large data files containing ascii
        representations of floating point data."""
        f=open(self.path)
        contents = f.readlines()
        f.close()

        if self.dialect is None:
            self.sniff()

        delim=self.dialect.delimiter

        contents[0:startrow]=[]
        if maxrows is not None:
            contents=contents[0:maxrows]

        testrow = contents[0].split(delim)
        try:
            float(testrow[0])
        except ValueError:
#            print('In ReadFloatRows, the first element of testrow could not be made into a float:'+testrow[0])
            return SpreadSheet.ReadRows(self, startrow=startrow, \
                                        maxrows=maxrows, **kwargs)
        try:
            map(float,testrow)
        except ValueError:
            return SpreadSheet.ReadRows(self, startrow=startrow, \
                                        maxrows=maxrows, **kwargs)
        row1 = fromstring(contents[0],sep=delim)
        mymat = zeros((len(contents), len(row1)))

        for n, line in enumerate(contents):
            currow = fromstring(line, sep=delim)
            mymat[n,:] = currow
        return mymat


class TrueCSVSpreadSheet(CSVSpreadSheet):
    def __init__(self, pathin=None, dialect=mycsv, **kwargs):
        CSVSpreadSheet.__init__(self, pathin=pathin, \
                                dialect=dialect, \
                                **kwargs)


class ExcelSpreadSheet(SpreadSheet):
    def __init__(self, pathin=None, datasheet=0, skiprows=0):
        """Initialize an ExcelSpreadSheet derived from SpreadSheet.
        datasheet is the index of the sheet containing the data that
        will be searched, loaded, or otherwise analyzed."""
        self.datasheet = datasheet
        self.activesheet = None
        SpreadSheet.__init__(self, pathin=pathin, skiprows=skiprows)
        self.tunits = 'ms'

    def SetActiveSheet(self):
        self.activesheet = get_workbook_sheet(self.path, self.datasheet)

    def iterrows(self):
        if self.activesheet is None:
            self.SetActiveSheet()
        for n in range(self.activesheet.nrows):
            curslice=self.activesheet.row_slice(n)
            currow=[str(item.value) for item in curslice]
            yield currow

class DynatupSpreadSheet(CSVSpreadSheet):
    def __init__(self, pathin=None, dialect=mycsv, collabels=['Time','Ch1Load/Curve','Ch1Raw']):
        self.dialect = dialect
        SpreadSheet.__init__(self, pathin = pathin)
        self.collabels = collabels
        self.impvel = -1
        if pathin is not None:
            #print pathin
            self.GetImpactVelocity()
            self.FindLabelRow()
            self.ReadDataColumns()

    def FFT(self):
        if not hasattr(self, 't'):
            self.t = self.tms/1000.0
        self.Ffft = DataProcMixins.FFTChannel(self.FkN, self.t)

    #mycsv = spreadsheet.CSVSpreadSheet('12210601.exp')
    def GetImpactVelocity(self):
        vals, labels = self.GetHeaderInfo(['Impact velocity (1)'])
        self.impvel = float(vals[0])


    def FindLabelRow(self):
        mylabels=['Index','Time']
        CSVSpreadSheet.FindLabelRow(self, labels=mylabels, exact=True)

    def FindDataColumns(self):
        CSVSpreadSheet.FindDataColumns(self, self.collabels, exact=True)

    def ReadDataColumns(self, **kwargs):
        if not self.datacolumns:
            self.FindDataColumns()
        CSVSpreadSheet.ReadDataColumns(self, skiprows=1, **kwargs)
        self.tms = self.data[:,0]
        self.FkN = self.data[:,1]
        self.Fraw = self.data[:,2]


    def Plot_F_vs_t(self, fig, clear=True):
        if clear:
            fig.clf()
        ax=fig.add_subplot(111)
        ax.plot(self.tms,self.FkN)
        ax.set_xlabel('Time (ms)')
        ax.set_ylabel('Force (kN)')
        #Do something with units here

    def Plot_FkN_vs_raw(self, fig, clear=True):
        if clear:
            fig.clf()
        ax=fig.add_subplot(111)
        ax.plot(self.Fraw,self.FkN)
        ax.set_xlabel('Counts')
        ax.set_ylabel('Force (kN)')

    def CheckCal(self, zeroc, maxv, calf, nbits):
        numcounts = 2**nbits
        self.cal1 = (self.Fraw+0.0)/numcounts*maxv*calf
        self.cal2 = (self.Fraw+0.0-zeroc)/numcounts*maxv*calf

    def PlotCal(self, fig, clear=True):
        if clear:
            fig.clf()
        ax=fig.add_subplot(111)
        ax.plot(self.FkN, self.cal1, self.FkN, self.cal2, self.FkN, self.FkN)
        ax.set_ylabel('Force (kN) - ??')
        ax.set_xlabel('Force (kN)')
        ax.legend(['Cal 1','Cal 2','Check'])

    def Plot_x_vs_t(self, fig, clear=True):
        myax = self._Plot_prop_vs_t('x', fig, clear=clear)
        myax.set_ylabel('Displacement (m)')

    def Plot_v_vs_t(self, fig, clear=True):
        myax = self._Plot_prop_vs_t('v', fig, clear=clear)
        myax.set_ylabel('Velocity (m/s)')


    def _Plot_prop_vs_t(self, prop, fig, clear=True):
        if clear:
            fig.clf()
        ax=fig.add_subplot(111)
        vect = getattr(self,prop)
        ax.plot(self.tms, vect)
        ax.set_xlabel('Time (ms)')
        return ax


    def cumtrapztest(self, mass, scale=1.0):
        myt = self.tms/1000.0
        v = self.impvel - integrate.cumtrapz(self.FkN*1000.0/mass*scale, myt)
        v=r_[v,v[-1]]
        x = integrate.cumtrapz(v, myt)
        x=r_[x,x[-1]]
        self.x = x
        self.v = v
        return x, v

class TabDelimSpreadSheet(CSVSpreadSheet):
    def __init__(self, pathin=None, dialect=tabdelim, skiprows=0, collabels=None, colmap=None, datafunc=float, picklekeys=[]):
        self.dialect = dialect
        SpreadSheet.__init__(self, pathin=pathin, skiprows=skiprows, collabels=collabels, colmap=colmap, datafunc=datafunc, picklekeys=picklekeys)


class LabviewSpreadSheet(TabDelimSpreadSheet):
    def __init__(self, pathin=None, tlabel=None, collabels=None, dialect=tabdelim, colmap=None, datafunc=float, exact=False, removeempty=False):
        """Create an instance of the LabviewSpreadSheet class, where
        the data is assumed to be in a tab delimited text file.

        If colmap is not None, it is expected to be a dictionary with
        column labels as keys and the corresponding attributes as
        values.  For example, if a data file contained 3 columns with
        labels 't', 'Light Gate', and 'Accel (V)', colmap could be
        {'Light Gate':'lg', 'Accel':a}.  This would result in self.lg
        and self.a having the appropriate data columns in them.
        datafunc is the parsefunc passed to self.ReadDataColumn.
        exact and removeempty are keyword arguments passed to
        self.ReadDataColumn - which is only called if colmap is not
        None."""
        self.datafunc = datafunc
        TabDelimSpreadSheet.__init__(self, pathin=pathin, dialect=dialect)
        if (pathin is not None) and tlabel and collabels:
            self.FindLabelRow([tlabel])
            self.FindDataColumns(collabels)
            self.ReadData()
            self.ReadDataColumns()
            self.t = self.ReadDataColumn(tlabel, parsefunc=datafunc)
            if colmap is not None:
                for key, val in colmap.iteritems():
                    curcol = self.ReadDataColumn(key, parsefunc=self.datafunc, exact=exact, removeempty=removeempty)
                    setattr(self, val, curcol)
                self.colmap = colmap
            else:
                self.lg = self.get_data_col(1)
                self.a = self.get_data_col(2)


    def __call__(self, t):
        ind1 = thresh(self.t, t)
        if ind1 > 0:
            ind2 = ind1-1
        elif t==0.0 and self.t[0]==0.0:
            return self.a[0]
        else:
            raise ValueError, "Called LabviewSpreadSheet with t < 0."
        a1 = self.a[ind1]
        a0 = self.a[ind2]
        t1 = self.t[ind1]
        t0 = self.t[ind2]
        return interp(t, t0, t1, a0, a1)


class TXTDataFile(TabDelimSpreadSheet):
    def __init__(self, pathin=None, tlabel=None, collabels=None, dialect=tabdelim, colmap=None, datafunc=float):#, exact=False, removeempty=False):
        TabDelimSpreadSheet.__init__(self, pathin=pathin, dialect=dialect, collabels=collabels, colmap=colmap, datafunc=datafunc)
        if pathin:
            if tlabel is None:
                tlabel = 'time'
            else:
                tlabel = tlabel.lower()
            self.ReadHeader(headerrows=20)
            A1 = self.header[0][0].lower()
            if A1.find(tlabel)>-1:
                self.labelrow = 0
                self.labels = self.header[0]
                if not colmap:
                    colmap = {}
                    for item in self.labels:
                        prop = copy.copy(item)
                        k = prop.find('(')
                        if k > 0:
                            prop = prop[0:k]
                        prop = prop.strip()
                        prop = prop.replace(' ','_')
                        if prop == 'Time':
                            prop='t'
                        colmap[item] = prop
                self.colmap = colmap
                self.MapCols()


class LabeledDataFile(TXTDataFile):
    """A class for files whose columns are well labeled (may have
    issues with spaces).  Automatically creates a column map from the
    labels."""
    def _make_clean_colmap(self):
        """Assuming self.labels is already defined, use it to create a
        colmap."""
        colmap = None
        for item in self.labels:
            prop = copy.copy(item)
            k = prop.find('(')
            if k > 0:
                prop = prop[0:k]
            prop = prop.strip()
            prop = prop.replace(' ','_')
            if prop == 'Time':
                prop='t'
            if colmap is None:
                colmap = {item:prop}
            else:
                colmap[item] = prop
        self.colmap = colmap

    def __init__(self, pathin, skiprows=0, poprows=0, dialect=None, \
                 colmap=None):
        TabDelimSpreadSheet.__init__(self, pathin=pathin, \
                                     dialect=dialect)
        self.ReadData()
        for i in range(skiprows):
            self.alldata.pop(0)
        for i in range(poprows):
            self.alldata.pop()#from the end
        self.labelrow = 0
        self.labels = self.alldata.pop(0)
        if colmap is None:
            self._make_clean_colmap()
        else:
            self.colmap = colmap
        self.MapCols()


class email_list(CSVSpreadSheet):
    """A spreadsheet class for looking up emails in a class list.
    Typically, there is one row of labels and the first two columns
    contain student names and emails."""
    def __init__(self, pathin, name_col=0, email_col=1, \
                 labelrow=0, dialect=None):
        CSVSpreadSheet.__init__(self, pathin=pathin, dialect=dialect)
        self.labelrow=labelrow
        self.GetLabelRow()
        self.ReadData()
        if self.labels[1].lower().find('first') > -1:
            #the last names are in column 0 and the first names are in
            #column 1
            self.first_names = txt_mixin.txt_list(self.get_col(1))
            self.last_names = txt_mixin.txt_list(self.get_col(0))
            names = []
            for first, last in zip(self.first_names, self.last_names):
                curname = last + ', ' + first
                names.append(curname)
            self.names = txt_mixin.txt_list(copy.copy(names))

            for i, item in enumerate(self.labels):
                if item.lower().find('email') > -1:
                    email_col = i
        else:
            self.names = txt_mixin.txt_list(self.get_col(name_col))
        self.emails = txt_mixin.txt_list(self.get_col(email_col))


    def get_email(self, lastname, firstname=None, aslist=False, \
                  fail_quietly=False):
        """Search for email in list using lastname or lastname,
        firstname.  aslist=True forces email to be returned in list
        even if it is only one address.  Students can have more than
        one list if the entry in that column is delimited by [,; ]."""
        p = re.compile('[,; ]+')
        name = lastname
        if firstname is not None:
            name += ', '+firstname
        inds = self.names.findall(name)
        if not fail_quietly:
            assert len(inds)==1, "Did not find exactly one match for %s, %s.  len(inds)=%i" % \
                  (lastname, firstname, len(inds))
        if len(inds) == 0:
            print('did not find %s in self.names' % name)
        if len(inds) > 1:
            print('found more than 1 %s in self.names' % name)
        email = self.emails[inds[0]]
        email = email.strip()
        q = p.search(email)
        if q or aslist:
            return p.split(email)
        else:
            return email

    def get_emails(self, lastnames, firstnames=None):
        if firstnames is None:
            firstnames = [None]*len(lastnames)
        email_list = None
        for lastname, firstname in zip(lastnames, firstnames):
            curemails = self.get_email(lastname, firstname, aslist=1)
            if email_list is None:
                email_list = curemails
            else:
                email_list.extend(curemails)
        return email_list


class group_list(LabeledDataFile):
    """A class list to find the members of a group given the group or
    team name.  The spreadsheet has one label row, group names in the
    first column, and team members in the second column."""
    def __init__(self, pathin, team_name_col=0, members_col=1, \
                 labelrow=0, dialect=None):
        CSVSpreadSheet.__init__(self, pathin=pathin, dialect=dialect)
        self.labelrow=labelrow
        self.GetLabelRow()
        if self.labels[0] == "Team #":
            team_name_col = 1
            members_col = 2
        self.ReadData()
        self.Project_Name = txt_mixin.txt_list(self.get_col(team_name_col))
        self.Group_Members = txt_mixin.txt_list(self.get_col(members_col))


    def find_group(self, group_name):
        ind = self.Project_Name.findall(group_name)
        assert len(ind)==1, 'Did not find exactly 1 group called ' + \
               group_name + ' len(ind) = ' + str(len(ind))
        return ind[0]

    def get_team_members(self, group_name):
        ind = self.find_group(group_name)
        members = self.Group_Members[ind]
        return members

    def _get_names(self, member_string):
        name_list = member_string.split(',')
        name_list = [item.strip() for item in name_list]
        lastnames = None
        firstnames = None
        for name in name_list:
            if name.find('and ') == 0:
                name = name.replace('and ','')
                name = name.strip()
            first, last = name.split(' ',1)
            last = last.strip()
            first = first.strip()
            if lastnames is None:
                lastnames = [last]
                firstnames = [first]
            else:
                lastnames.append(last)
                firstnames.append(first)
        return lastnames, firstnames

    def get_last_names(self, group_name):
        members = self.get_team_members(group_name)
        last, first = self._get_names(members)
        return last

    def get_names(self, group_name):
        members = self.get_team_members(group_name)
        return self._get_names(members)


class group_list_2010(group_list):
    """A class list to find the members of a group given the group or
    team name.  The spreadsheet has one label row, group names in the
    first column, and team members in the second column."""
    def __init__(self, pathin, team_name_col=1, members_col=2, \
                 labelrow=0, dialect=None, team_number_col=0):
        group_list.__init__(self, pathin, team_name_col=team_name_col, \
                            members_col=members_col, labelrow=labelrow, \
                            dialect=dialect)
        self.Team_Number = txt_mixin.txt_list(self.get_col(team_number_col))


class mini_project_group_list(group_list):
    """A class list to find the members of a group given the group or
    team name.  The spreadsheet has one label row, group numbers in
    the first column, and team members in the second column. The file
    is probably tab delimited"""
    def __init__(self, pathin, team_num_col=0, members_col=1, \
                 labelrow=0, dialect=None):
        f = open(pathin, 'rb')
        first_row = f.readline()

        if first_row.find('\t') > -1:
            TabDelimSpreadSheet.__init__(self, pathin=pathin, \
                                         dialect=tabdelim)
        else:
            CSVSpreadSheet.__init__(self, pathin=pathin, dialect=dialect)
        self.labelrow=labelrow
        self.GetLabelRow()
        self.ReadData()
        self.Project_Name = txt_mixin.txt_list(self.get_col(team_num_col))
        self.team_number = map(int, self.Project_Name)
        self.Group_Members = txt_mixin.txt_list(self.get_col(members_col))


    def get_names(self, group_name):
        if group_name.find('Team') == 0:
            group_name = group_name[5:]
        members = self.get_team_members(group_name)
        return self._get_names(members)


class JCILabviewSpreadSheet(LabviewSpreadSheet, DataProcMixins.AccelMixin):
    def __init__(self, pathin=None, tlabel='Time', \
                 collabels=['Time','Light Gate','Accel'], \
                 dialect=tabdelim, ascale=9.81*1000.0/5.0):
        """Create an instance of the JCILabviewSpreadSheet class,
        where the data is assumed to be in a tab delimited text file.
        This data file was used in Dynatup testing with an
        accelerometer on the back of the flat impactor."""
        TabDelimSpreadSheet.__init__(self, pathin=pathin, dialect=dialect)
        if (pathin is not None) and tlabel and collabels:
            self.FindLabelRow([tlabel])
            self.FindDataColumns(collabels)
            self.ReadDataColumns()
            self.t = self.get_data_col(0)
            self.dt = (self.t.max()-self.t.min())/(len(self.t)-1)
            self.fs = 1.0/self.dt
            self.lg = self.get_data_col(1)
            self.a = self.get_data_col(2)
            self.ascale = ascale
            self.afft = FFTChannel(self.a, self.t)
            self.afilter = FilterChannel(self.a, 5000.0, self.fs, tvect=self.t)

    def Plot_a_vs_t(self, fig, clear=True):
        myax = self._Plot_prop_vs_t('a', fig, clear=clear, tprop='t')
        myax.set_ylabel('Acceleration')
        return myax

    def SWTrunc(self, athresh=0.1, duration = 0.05, backup=0.01):
        mythresh = athresh+self.a.mean()
        ind1 = thresh(self.a, athresh)
        t0 = self.t[ind1]
        self.t0 = t0
        try:
            ind2 = thresh(self.t, t0+duration)
        except:
            ind2 = len(self.t)
        ind0 = thresh(self.t, t0-backup)
        self.atrunc = copy.copy(self.a[ind0:ind2])
        self.ttrunc = copy.copy(self.t[ind0:ind2])-t0

    def Plot_a_vs_t_trunc(self, fig, clear=True):
        myax = self._Plot_prop_vs_t('atrunc', fig, clear=clear, tprop='ttrunc')
        myax.set_ylabel('Acceleration')
        return myax


    def FilterTrunc(self, cutoff = 1000.0):
        self.afilterch = FilterChannel(self.atrunc, cutoff_freq=cutoff, sample_freq=self.fs, order=2, tvect=self.ttrunc)
        self.afilterch.Filter()


class JCILabviewHammerTest(LabviewSpreadSheet):
    def __init__(self, pathin=None, tlabel='Time', collabels=['Time','Hammer','Accel'], dialect=tabdelim, colmap={'Hammer':'hammer','Accel':'a'}):
        LabviewSpreadSheet.__init__(self, pathin=pathin, tlabel=tlabel, collabels=collabels, dialect=dialect, colmap=colmap)
        self.picklekeys = colmap.values()
        self.picklekeys.append('t')


class LSDynaCSVFile(CSVSpreadSheet):
    def ReadData(self, skiprows=0, maxrows=None, parsefunc=float):
        CSVSpreadSheet.ReadData(self, skiprows=skiprows, maxrows=maxrows, parsefunc=parsefunc)
        self.data=array(self.alldata)
        self.t = self.data[:,0]
        if self.params and self.datacols:
            for curparam, curcol in zip(self.params, self.datacols):
                setattr(self, curparam, self.data[:,curcol])

    def __init__(self, pathin=None, datacols=[], params=[],  dialect=lsdcsv):
        CSVSpreadSheet.__init__(self, pathin=pathin, dialect=dialect)
        self.datacols = datacols
        self.params = params
        self.labelrow = 0
        self.GetLabelRow()
        self.ReadData()


    def iterrows(self):
        """The LS-Dyna csv files have an extra comma at the end of a
        row that the Python CSV module interprets as an empty cell.
        So, this iterrows method chops off the last element in each
        row."""
        if self.dialect is None:
            self.sniff()
        reader = csv.reader(open(self.path,'rb'), self.dialect)
        for row in reader:
            rowout = row[0:-1]
            yield rowout


class BlackBoardGBFile(CSVSpreadSheet):
    def __init__(self, pathin='gb_export.csv', dialect=mycsv, **kwargs):
        CSVSpreadSheet.__init__(self, pathin=pathin, dialect=dialect, **kwargs)
        self.labelrow=0
        self.GetLabelRow()
        self.ReadData()
        if self.alldata[0][0].find('AAStudent') > -1:
            self.alldata.pop(0)
        self.clean_empty_rows()
        self.ParseNames()

    def assign_col_to_attr(self, collabel, attr, array=1, dtype=float):
        ind = self.labels.index(collabel)
        col = self.get_col(ind)
        if array:
            col = numpy.array(col, dtype=dtype)

        setattr(self, attr, col)

    def assign_letter_grades(self, attr='course_grades', \
                             cutoffs=None, append=True):
        N = len(self.lastnames)
        #letter_grades = numpy.zeros(N, dtype='S1')
        letter_grades = array(['F']*N)
        num_grades = numpy.zeros(N)
        if cutoffs is None:
            cutoffs = [59.5, 69.5, 79.5, 89.5]
        grades = getattr(self, attr)
        gpas = [1.0, 2.0, 3.0, 4.0]
        letters = ['D','C','B','A']
        for cutoff, gpa, letter in zip(cutoffs, gpas, letters):
            num_grades = numpy.where(grades > cutoff, \
                                     gpa, num_grades)
            letter_grades = numpy.where(grades > cutoff, \
                                        letter, letter_grades)
        self.num_grades = num_grades
        self.letter_grades = letter_grades
        self.course_gpa = self.num_grades.mean()
        if append:
            self.AppendColFromList(self.lastnames, \
                                   'Letter Grade', \
                                   self.letter_grades, \
                                   splitnames=False)
            self.AppendColFromList(self.lastnames, \
                                   'GPA', \
                                   self.num_grades, \
                                   splitnames=False)


    def save(self, csvpathout):
        self.WriteAllDataCSV(csvpathout)


    def clean_empty_rows(self):
        """Clean rows at the bottom of self.alldata that may be used
        for averaging or other analysis.  Any row that doesn't start
        with a lastname, firstname entry is popped (actually, we just
        search for the comma for now)."""
        N = len(self.alldata)
        i = 0
        while i < N:
            i += 1
            lastrow = self.alldata.pop()
            if lastrow[0].find(',') > -1:
                self.alldata.append(lastrow)
                break


    def Lastname_from_First(self, firstname, alternates={}):
        indlist = []
        return NotImplementedError


    def Firstname_from_Last(self, lastname):
        lastname = lastname.strip()
        lastname = lastname.capitalize()
        lastname = lastname.replace(' ','')
        try:
            ind = self.lastnames.index(lastname)
        except ValueError:
            ind = self.upperlast.index(lastname.upper())
        return self.firstnames[ind].capitalize()


    def ParseNames(self):
        lastnames = []
        firstnames = []
        uids = []
        for row in self.alldata:
            namemess = row[0]
            if namemess.find(',') > -1:
                lastname, rest = namemess.split(',',1)
                lastname = lastname.strip()
                lastnames.append(lastname)
                if rest.find('(') > -1:
                    firstname, rest2 = rest.split('(',1)
                    uid, rest3 = rest2.split(')',1)
                    uid=uid.strip()
                else:
                    firstname = rest
                    uid = None
                firstname = firstname.strip()
                if firstname.find(' ') > -1:
                    firstname, mi = firstname.split(' ',1)
                    firstname = firstname.strip()
                uids.append(uid)
                firstnames.append(firstname)
        self.lastnames = lastnames
        self.firstnames = firstnames
        temp = [last.upper() for last in self.lastnames]
        self.upperlast =[last.replace(' ','') for last in temp]
        self.uids = uids


    def InsertColFromList(self, namelist, destlabel, valuelist, \
                          splitnames=True, verbosity=2):
        """This is an attempt to mimic the InsertCol method of the
        Spreadsheet class, but taking into account that the
        BlackBoardGBFile has one column containg the last name, first
        name (uid).

        This method assumes we are matching student names.  If
        splitnames=True, then it is further assumed that namelist
        contains a list of 'FirstName LastName' with a space
        seperating them.

        Last Name and First Name are now in two separate columns, but
        I have an issue for 482, Fall 2010: two Miller's."""
        destcol = self.labels.index(destlabel)
        if splitnames:
            lastnames = []
            for curname in namelist:
                #handling two cases for now:
                # 1. last, first
                # 2. first last
                if curname.find(',') > -1:
                    last, first = curname.split(',')
                else:
                    first, last = curname.split(' ')
                last = last.strip()
                lastnames.append(last)
            namelist = lastnames
        found = 0
        for curname, row in zip(self.lastnames, self.alldata):
            curind = search_list(namelist, curname, match=1)
            #print('curind = ' + str(curind))
            if curind == -1:
                if verbosity > 1:
                    print('could not find ' + curname)
                    print('namelist = ' + str(namelist))
                    print('destlabel = ' + destlabel)
                    raise ValueError
            else:
                #if valuelist[curind]:#I don't know why this was here,
                #I guess to allow empty cells not mess up averages
                found += 1
                N = len(row)-1
                if verbosity > 10:
                    print('found: ' + curname)
                    print('curind = '+ str(curind))
                    print('lastnames[curind] = ' +self.lastnames[curind])
                if destcol > N:
                    row += [None]*(destcol-N)
                    #row.append(valuelist[curind])
                row[destcol] = valuelist[curind]
        if (found != len(namelist)) and (verbosity >0):
            print('found = ' + str(found))
            print('len(namelist) = '+str(len(namelist)))


    def find_student(self, lastname, firstname, verbosity=0):
        last_inds = self.lastnames.findall(lastname)
        assert len(last_inds) > 0, "Did not find a student with the lastname " + str(lastname)
        if len(last_inds) == 1:
            return last_inds[0]
        else:
            possible_firstnames = []
            for ind in last_inds:
                possible_firstnames.append(self.firstnames[ind])
            possible_firstnames = txt_mixin.txt_list(possible_firstnames)
            first_inds = possible_firstnames.findall(firstname)
            if verbosity > 0:
                print('firstname = ' + firstname)
                print('found: ' + str(self.alldata[last_inds[first_inds[0]]]))
            assert len(first_inds) == 1, "Did not find exacly one firstname for lastname " + str(lastname)
            return last_inds[first_inds[0]]


    def InsertColFromList_v2(self, lastnames, firstnames, \
                             destlabel, valuelist, \
                             verbosity=2):
        """Fall 2010 is the first time I had two students with the
        same last name in one class.  I am writing this new method to
        address this issue."""
        destcol = self.labels.index(destlabel)

        for last, first, value in zip(lastnames, firstnames, valuelist):
            rowind = self.find_student(last, first)
            row = self.alldata[rowind]
            N = len(row)-1
            if destcol > N:
                row += [None]*(destcol-N)
                    #row.append(valuelist[curind])
            row[destcol] = value


    def AppendColFromList(self, namelist, destlabel, valuelist, \
                          splitnames=True):
        self.labels.append(destlabel)
        self.InsertColFromList(namelist, destlabel, valuelist, \
                               splitnames=splitnames)


    def Append_Col_of_Zeros(self, val='0.0'):
        for row in self.alldata:
            row.append(val)


    def AppendColFromList_v2(self, lastnames, firstnames, \
                             destlabel, valuelist, verbosity=0):
        self.labels.append(destlabel)
        self.InsertColFromList_v2(lastnames, firstnames, \
                                  destlabel, valuelist)


    def Append_From_GradeSpreadSheet(self, gradesheet, labels=None, \
                                     parsefunc=myfloat, names_in=None):
        if labels is None:
            labels = gradesheet.valuelabels
        if hasattr(gradesheet, 'values'):
            names = gradesheet.names
            values = gradesheet.values
        else:
            names, values = gradesheet.ReadNamesandValues()
        if names_in is None:
            lastnames = gradesheet.lastnames
        else:
            lastnames = names_in
        if (len(labels)==1) and type(values)==list:
            self.AppendColFromList(lastnames, labels[0], \
                                   values, splitnames=False)
        else:
            for label, col in zip(labels, values.T):
                self.AppendColFromList(lastnames, label, \
                                       col, splitnames=False)


class BlackBoardGBFile_v_8_0(BlackBoardGBFile):
    def __init__(self, pathin='gb_export.csv', dialect=mycsv, **kwargs):
        CSVSpreadSheet.__init__(self, pathin=pathin, dialect=dialect, **kwargs)
        self.labelrow=0
        self.GetLabelRow()
        self.ReadData()
        if self.alldata[0][0].find('AAStudent') > -1:
            self.alldata.pop(0)
        self.ParseNames()


    def ParseNames(self):
        self.MapCols({'Last Name':'lastnames', \
                      'First Name':'firstnames', \
                      'Username':'uids'})
        self.lastnames = txt_mixin.txt_list(self.lastnames)
        self.firstnames = txt_mixin.txt_list(self.firstnames)


class Fake_BlackBoard_File(BlackBoardGBFile_v_8_0):
    def __init__(self, lastnames, firstnames=None):
        if firstnames is None:
            firstnames, lastnames = split_names(lastnames, reverse=1)
        self.lastnames = lastnames
        self.firstnames = firstnames
        self.uids = [None]*len(self.lastnames)
        self.labels = ['Last Name', 'First Name', 'Username']
        alldata_tups = zip(self.lastnames, self.firstnames, \
                           self.uids)
        self.alldata = [list(item) for item in alldata_tups]


class GradeSpreadSheet(CSVSpreadSheet):
    def __init__(self, pathin=None, namelabel='Name',valuelabel='Total', \
                 dialect=None, skiprows=0):
        CSVSpreadSheet.__init__(self, pathin=pathin, dialect=dialect, \
                                skiprows=skiprows)
        self.namelabel = namelabel
        self.valuelabel = valuelabel
        self.valuelabels = [valuelabel]#for Append_From_GradeSpreadSheet
        self.FindLabelRow([namelabel])
        self.FindDataColumns([namelabel,valuelabel],exact=True)
        self.ReadData()
        namecol = self.labels.index(namelabel)
        self.DropEmptyRows(cols=[namecol])#drop all rows that do not
                                          #have a name in them
        self.names = self.ReadDataColumn(self.namelabel, exact=True)
        self.split_names()

    def split_names(self):
        """Take self.names (which is assumed to already exist) and
        split into self.lastnames and self.firstnames.  If a name in
        self.names has a comma, it is assumed to be last, first (with
        possible initials).  If the name contains no commas, it is
        assumed to be first (middle) last, where (MI) is optional."""
        lastnames = []
        firstnames = []
        for name in self.names:
            if name.find(',') > -1:
                #assume last, first (middle)
                last, rest = name.split(',', 1)
                last = last.strip()
                rest = rest.strip()
                if rest.find(' ') > -1:
                    first, middle = rest.split(' ', 1)
                    first = first.strip()
                else:
                    first = rest
            else:
                if name.find(' ') > -1:
                    #assume first (middle) last
                    first, rest = name.split(' ')
                    first = first.strip()
                    rest = rest.strip()
                    if rest.find(' ') > -1:
                        #assume middle last
                        middle, last = rest.split(' ', 1)
                        last = last.strip()
                    else:
                        last = rest
                else:
                    #with no comma and no space, assume just a last name
                    last = name
                    first = ''
            last = last.replace(' ','')#Van Pelt
            last = last.replace("'",' ')#O'Donnel
            lastnames.append(last)
            firstnames.append(first)
        self.lastnames = lastnames
        self.firstnames = firstnames
        return self.lastnames, self.firstnames

    def ReadNamesandValues(self, exact=True, parsefunc=None):
        names=self.ReadDataColumn(self.namelabel, exact=exact)
        values=self.ReadDataColumn(self.valuelabel, exact=exact)
        if parsefunc is not None:
            values = map(parsefunc, values)
        return names, values


class GradeSpreadSheetMany(GradeSpreadSheet):
    """A class for extracting many (or at least more than one) grades
    from a spreadsheet.

    valuelabels contains a list of the column labels you want to
    extract."""
    def __init__(self, pathin=None, namelabel='Name',valuelabels=[], \
                 dialect=None, skiprows=0, split_names=None):
        CSVSpreadSheet.__init__(self, pathin=pathin, dialect=dialect, \
                                skiprows=skiprows)
        self.namelabel = namelabel
        self.valuelabels = valuelabels
        self.FindLabelRow([namelabel])
        self.FindDataColumns([namelabel]+valuelabels, exact=True)
        self.ReadData()
        namecol = self.labels.index(namelabel)

        self.DropEmptyRows(cols=[namecol])#drop all rows that do not
                                          #have a name in them
        self.names = self.ReadDataColumn(self.namelabel, exact=True)
        if split_names is None:
            if namelabel.lower().find('last') > -1:
                split_names = 0
            else:
                split_names = 1
        if split_names:
            self.split_names()
        else:
            self.lastnames = self.names

    def get_col(self, label):
        index = self.valuelabels.index(label)
        col = self.values[:,index]
        return col

    def append_col(self, label, col):
        self.valuelabels.append(label)
        self.values = numpy.append(self.values, col, -1)

    def _prep_average_or_total(self, new_label):
        self.valuelabels.append(new_label)
        nr, nc = self.values.shape
        zero_col = zeros((nr,1))
        self.values = numpy.append(self.values, zero_col, -1)

    def average_sheet(self, new_label, drop_lowest=False):
        self._prep_average_or_total(new_label)
        for n, row in enumerate(self.values):
            all_scores = row[0:-1]#there is a new 0 at the end for the ave
            if drop_lowest:
                imin = all_scores.argmin()
                top_scores = numpy.delete(all_scores, imin)
                ave = top_scores.mean()
            else:
                ave = all_scores.mean()
            self.values[n, -1] = ave

    def total_sheet(self, new_label, drop_lowest=False):
        self._prep_average_or_total(new_label)
        for n, row in enumerate(self.values):
            all_scores = row[0:-1]#there is a new 0 at the end for the ave
            if drop_lowest:
                imin = all_scores.argmin()
                top_scores = numpy.delete(all_scores, imin)
                total = top_scores.sum()
            else:
                total = all_scores.sum()
            self.values[n, -1] = total

    def ReadNamesandValues(self, exact=True, parsefunc=myfloat):
        nc = len(self.valuelabels)
        nr = len(self.names)
        values = zeros((nr,nc))
        for n, label in enumerate(self.valuelabels):
            curcol = self.ReadDataColumn(label, exact=exact)
            curcol = map(parsefunc, curcol)
            values[:,n] = array(curcol)
        self.values = values
        return self.names, self.values


class QuizScoreSpreadSheet(GradeSpreadSheetMany):
    """A class for reading in quiz scores and dropping the N lowest."""
    def __init__(self, *args, **kwargs):
        GradeSpreadSheetMany.__init__(self, *args, **kwargs)
        self.names, self.quiz_scores = self.ReadNamesandValues()
        self.ReadDataColumns()


    def Drop_Scores_and_Average(self, N):
        """Drop the lowest N quiz scores."""
        self.scores_kept = []
        sorted_scores = copy.copy(self.quiz_scores)
        for row in sorted_scores:
            row.sort()
            cur_keep = row[N:]
            self.scores_kept.append(cur_keep)
        self.scores_kept = array(self.scores_kept)
        self.quiz_averages = self.scores_kept.mean(axis=1)
        self.AppendCol('Quiz Average', self.quiz_averages)





class Survey_Answer(object):
    def count_answers(self):
        N = len(self.choices)
        self.histogram = numpy.zeros(N)
        for i, choice in enumerate(self.choices):
            count = self.answers.count(choice)
            self.histogram[i] = count
        self.percentages = self.histogram/self.N*100.0


    def _get_fig(self, fig=None, fignum=1, figsize=None):
        if fig is None:
            import pylab as P
            fig = P.figure(fignum, figsize=figsize)
        return fig

    def build_caption(self, caption=None, add_key=False):
        if caption is None:
            caption = ['Responses to survey question number %i.' % \
                      self.number]
        if add_key:
            caption.append('`\\newline`')

            caption.append('Key: `\\newline`')
            for n, answer in enumerate(self.choices):
                curline ='%i = "%s"' % (n+1, answer)
                if n < self.NC-1:
                    curline += ' `\\newline`'
##                 caption +=', '
                elif n == self.NC-1:
                    curline += ' `\\label{fig:surveyQ%i}`' % self.number
                    caption.append(curline)
        else:
            caption.append('`\\label{fig:surveyQ%i}`' % self.number)

        self.caption = caption


    def get_fig_path(self, name=None, folder='figs', ext='.eps'):
        if name is None:
            name = 'question_%0.2i' % self.number

        name_only, old_ext = os.path.splitext(name)
        if not old_ext:
            name = name_only + ext
        if not os.path.exists(folder):
            os.mkdir(folder)
        path = os.path.join(folder, name)
        return path


    def save_fig(self, name=None, fig=None, fignum=1, folder='figs', \
                 ext='.eps'):
        path = self.get_fig_path(name, folder, ext)
        fig = self._get_fig(fig, fignum)
        mplutil.mysave(path, fig)

    def save_rst(self, caption=None, name=None, folder='figs', \
                 relfolder=None):
        #Pdb().set_trace()
        path = self.get_fig_path(name, folder, ext='.rst')
        self.build_caption(caption=caption)
        junk, nameonly = os.path.split(path)
        fno, old_ext = os.path.splitext(nameonly)
        pdfname = fno + '.pdf'
        if relfolder is None:
            relfolder = folder
        relpath = os.path.join(relfolder, pdfname)
        ws = ' '*4
        rstlist = ['.. figure:: ' + relpath]
        rstlist.append(ws + ':width: 4.5in')
        rstlist.append('')
        for curline in self.caption:
            rstlist.append(ws + curline)
        rstlist.append('')
        txt_mixin.dump(path, rstlist)



    def plot_histogram(self, fig=None, fignum=1, clear=True, \
                       title=None, xlim=None, use_percentages=False, \
                       figsize=None):
        fig = self._get_fig(fig, fignum, figsize=figsize)

        if title is None:
            title = self.question
        if clear:
            fig.clf()
        if xlim is None:
            low = 0.5
            high = self.NC + 0.5
            xlim = [low, high]
        ax = fig.add_subplot(1,1,1)
        #ax = fig.add_axes([0.1,0.1,0.8,0.8])
        ind = range(1,self.NC+1)
        if use_percentages:
            ax.bar(ind, self.percentages, align='center')
            ax.set_ylabel('Percentage of Students')
        else:
            ax.bar(ind, self.histogram, align='center')
            ax.set_ylabel('Number of Students')
        ax.set_title(title)
        ax.set_xticks(ind)
        ax.set_xlabel('Response')
        fig.autofmt_xdate(rotation=0.0)
        if xlim:
            ax.set_xlim(xlim)


    def add_percentages_to_labels_for_1_student(self, N):
        """If only 1 student gave a certain answer, move the apct with
        percentage out to the label."""

        for i, label in enumerate(self.nonempty_histogram):
            cur_n = float(label)
            if cur_n == 1:
                cur_label = self.clean_pie_labels[i]
                p = 100.0*cur_n/N
                p_label = '%i (%0.4g' % (cur_n, p) + '%)'
                new_label = cur_label + ':\n ' + p_label
                self.clean_pie_labels[i] = new_label


    def plot_pie_chart(self, fig=None, fignum=1, clear=True, title=None, \
                       figsize=(11,8), use_percentages=False, \
                       use_title=True, labeldistance=1.1, shadow=True, N=None, \
                       empty_x=None, empty_y=None):
        fig = self._get_fig(fig, fignum, figsize=figsize)
        if title is None:
            title = self.question
        if clear:
            fig.clf()

        #ax = fig.add_subplot(1,1,1)
        ax = fig.add_axes([0.25,0.15,0.45,0.625])
        #Pdb().set_trace()
        mygrays = []
        for value in arange(1,0,-0.15):
            mytup = (value, value, value)
            mygrays.append(mytup)

        if use_percentages:
            slices = self.nonempty_percentages
        else:
            slices = self.nonempty_histogram

        if (N is not None) and (not use_percentages):
            self.add_percentages_to_labels_for_1_student(N)


        ret = ax.pie(slices, labels=self.clean_pie_labels, \
                     autopct='%i', colors=mygrays, shadow=shadow, \
                     labeldistance=labeldistance)

        wedges = ret[0]
        apcts = ret[-1]
        if not use_percentages:
            for a, w, label in zip(apcts, wedges, self.nonempty_histogram):
                if N is not None:
                    val = float(label)
                    p = 100.0*val/N
                    mylabel = '%i\n(%0.4g' % (label, p) + '%)'
                    if val > 1:
                        a.set_text(mylabel)
                    else:
                        a.set_text('')
                else:
                    a.set_text('%i' % label)

        print('ret = ' + str(ret))
        print('len(ret) = ' + str(len(ret)))

        ## ax.pie(self.nonempty_percentages, labels=self.clean_pie_labels, \
        ##        autopct='%1.1f%%', colors=mygrays, shadow=True)
        if use_title:
            t = ax.set_title(title)
            t.set_position((0.5, 1.05))

        if empty_x is None:
            empty_x = -2.0
            #empty_x = 2.0
        if empty_y is None:
            empty_y = -1.5
        dy = 0.15
        if hasattr(self, 'empty_x'):
            empty_x = self.empty_x
        if hasattr(self, 'empty_y'):
            empty_y = self.empty_y
        if hasattr(self, 'dy'):
            dy = self.dy

        NE = len(self.empty_percentages)
        if (NE > 0) and (NE < 3):
            if use_percentages:
                values = self.empty_percentages
            else:
                values = self.empty_histogram
            for label, value in zip(self.empty_labels, values):
                if use_percentages:
                    msg = '%s: %0.1f' % (label, value) + '%'
                else:
                    msg = '%s: %i' % (label, value)
                ax.text(empty_x, empty_y, msg)
                empty_y += dy


    def print_histogram(self):
        print('Histogram:')
        print('-'*15)
        for answer, num, percent in \
                zip(self.choices, self.histogram, self.percentages):
            print('%s : %i (%0.2f percent)' % (answer, num, percent))


    def build_table_row(self, question=None):
        if question is None:
            question = self.question

        histfloatlist = self.histogram.tolist()
        histstrlist = ['%i' % item for item in histfloatlist]
        mylist = [question] + histstrlist
        myrow = ' & '.join(mylist) + ' \\\\'
        return myrow


    def __init__(self, number, question, answers, choices=None, pretty_labels={}):
        self.number = number
        self.question = question
        self.answers = txt_mixin.txt_list(answers)
        self.pretty_labels = pretty_labels
        self.N = len(self.answers)
        if choices is None:
            self.choices = self.answers.find_unique()
        else:
            self.choices = choices
        self.NC = len(self.choices)
        self.count_answers()
        self.clean_pie_labels = copy.copy(self.choices)
        for i, key in enumerate(self.clean_pie_labels):
            if self.pretty_labels.has_key(key):
                self.clean_pie_labels[i] = self.pretty_labels[key]
        self.empty_percentages = []
        self.empty_labels = []
        self.nonempty_percentages = copy.copy(self.percentages)
        self.nonempty_histogram = copy.copy(self.histogram)
        self.empty_histogram = []
        N = len(self.nonempty_percentages)
        i = 0
        while i < N:
            value = self.nonempty_percentages[i]
            if value == 0:
                self.nonempty_percentages = numpy.delete(self.nonempty_percentages, i)
                self.nonempty_histogram = numpy.delete(self.nonempty_histogram, i)
                self.empty_percentages.append(value)
                self.empty_histogram.append(value)
                label = self.clean_pie_labels[i]
                self.clean_pie_labels = numpy.delete(self.clean_pie_labels, i)
                label = label.replace('\n',' ')
                label = label.replace('  ',' ')
                label = label.replace('  ',' ')
                self.empty_labels.append(label)
            else:
                i += 1
            N = len(self.nonempty_percentages)



class BlackBoard_Survey_File(BlackBoardGBFile):
    """Note: if you are getting NULL byte errors from the csv reader,
    you haven't cleaned the file of yucky xls characters.  Saving as
    ods, closing oocalc, then reopening the ods and saving as csv
    worked for me."""
    def __init__(self, pathin, dialect=mycsv, pretty_labels={}, \
                 **kwargs):
        CSVSpreadSheet.__init__(self, pathin=pathin, dialect=dialect, **kwargs)
        self.labelrow=0
        self.GetLabelRow()
        self.ReadData()
        self.data = {}
        self.pretty_labels = pretty_labels

    def Find_Question_Col(self, number):
        ind = self.FindColLabel('Question %i' % number)
        return ind


    def Find_Answer_Col(self, number):
        ind = self.FindColLabel('Answer %i' % number)
        return ind


    def Parse_One_Question(self, number, choices=None):
        qind = self.Find_Question_Col(number)
        aind = self.Find_Answer_Col(number)
        qdata = self.get_col(qind)
        question = qdata[0]
        answers = self.get_col(aind)
        self.data[number] = Survey_Answer(number, question, answers, \
                                          choices=choices, \
                                          pretty_labels=self.pretty_labels)



class PSoCLabviewSpreadSheet(LabviewSpreadSheet, DataProcMixins.AccelMixin):
    def __init__(self, pathin=None, tlabel='Time', \
                 collabels=['Time','u (V)','y (V)','v (V)','loop timer'], \
                 dialect=tabdelim, \
                 colmap={'u (V)':'u','y (V)':'y','v (V)':'v',\
                         'loop timer':'lt'}):
        """Create an instance of the PSoCLabviewSpreadSheet class,
        where the data is assumed to be in a tab delimited text file."""
        TabDelimSpreadSheet.__init__(self, pathin=pathin, dialect=dialect)
        self.colmap = colmap
        if (pathin is not None) and tlabel and collabels:
            self.FindLabelRow([tlabel])
            self.FindDataColumns(collabels)
            self.ReadDataColumns()
            self.t = self.get_data_col(0)
            self.dt = (self.t.max()-self.t.min())/(len(self.t)-1)
            self.fs = 1.0/self.dt
            self.MapCols()


def dump(matrix, pathout, labels=[], delim='\t', append=False):
    """Dump matrix to a file whose path is pathout, labeling the
    columns with labels and using delim as the delimiter.  delim
    should either be '\t' or ','.  If delim is not '\t', then a comma
    is used."""
    if append:
        f=open(pathout,'ab')
    else:
        f=open(pathout, 'wb')
    if delim=='\t':
        dialect=tabdelim
    else:
        dialect=mycsv
    writer = csv.writer(f,dialect, escapechar='~')
    if labels:
        writer.writerow(labels)
    writer.writerows(matrix)
    f.close()


def dumpcsv(matrix, pathout, labels=[], delim=',', append=False):
    dump(matrix, pathout, labels=labels, delim=delim, append=append)


def AddColswithKeys(vals1, labels1, vals2, labels2):
    """Intelligently add to columns, vals1 and vals2, by searching
    labels1 and labels2 to find out which rows should be added."""
    #This could probably be done cleaner with dictionaries
    mylabels = copy.copy(labels1)
    myvals = copy.copy(vals1)
    for curlabel, curval in zip(labels2, vals2):
        curind = search_list(mylabels, curlabel)
        if curind==-1:
            print('could not find '+curlabel)
        else:
            myvals[curind] += curval
    return mylabels, myvals



def BuildListforMatching(row, inds, collabels, defaults={}, null=''):
    """Build a list of values to be used along with
    SpreadSheet.SearchForMatch by finding the items of row
    corresponding to inds or subsituting from defaults if an index in
    inds is -1.  The last option is to use null if defaults does not
    have a key corresponding to the current label."""
    listout = []
    for n, curlabel in zip(inds, collabels):
        if n == -1:
            if defaults.has_key(curlabel):
                listout.append(defaults[curlabel])
            else:
                listout.append(null)
        else:
            listout.append(row[n])
    return listout

def SpreadsheetFromPath(pathin, **kwargs):
    pathnoext, ext = os.path.splitext(pathin)
    if ext == '.xls':
        return ExcelSpreadSheet(pathin, **kwargs)
    else:
        return CSVSpreadSheet(pathin, **kwargs)


def Dict_List_From_Path(pathin, labellist, collabels=None, **kwargs):
    """Read in a spreadsheet file where each row contains one record.
    Return a list of dictionaries with collabels as the keys and the
    corresponding values from each row as the values.

    labellist is a list of values to help the spreadsheet code find
    the labelrow.  This can also just be a string of the top left
    column label.

    collabels is a list of the columns you want extracted from the
    file and put into the dictionaries.  If not passed in, all the
    columns will be read."""
    mysheet = SpreadsheetFromPath(pathin, **kwargs)
    mysheet.FindLabelRow(labellist)
    if collabels is None:
        mysheet.collabels = mysheet.GetLabelRow()
    mysheet.ReadData()
    dict_list = []
    N = len(mysheet.alldata)
    for i in range(N):
        curdict = mysheet.Get_Dict_from_Row(i, mysheet.collabels)
        dict_list.append(curdict)
    return dict_list


def NewCSVWriter_f(fin, dialect=mycsv):
    writer = csv.writer(fin, dialect=dialect)
    return writer


def NewCSVWriter(pathin, dialect=mycsv):
    writer = csv.writer(open(pathin, "wb"), dialect=dialect)
    return writer

def AppendableCSVWriter(pathin, dialect=mycsv):
    writer = csv.writer(open(pathin, "ab"), dialect=dialect)
    return writer

def determine_lt(linesin, dialectin):
    """Verify the line terminator of a dialect from my sniff function."""
    lts=[get_lineterminator(line) for line in linesin]
    tests=[item == lts[0] for item in lts]
    if all(tests):
        dialectin.lineterminator = lts[0]
    return dialectin


def sniff(filepath,sniffbytes=1000):
    """Try and determine the dialect of a text file by reading in the
    first sniffbytes bytes of data from filepath.  My slightly more
    intelligent version than the csv module."""
    f=open(filepath,'r')
    mylist=f.read(sniffbytes)
    mysniff=csv.Sniffer()
    cleanlines = [clean_line(item) for item in mylist]
    sniffstr=''.join(cleanlines)
    mydialect=mysniff.sniff(sniffstr)
    mydialect=determine_lt(mylist,mydialect)
    return mydialect

def getrows(filepath, maxrows=None, dialect=None, startrow=0):
    """Read data from filepath from startrow to startrow+maxrows.  If
    dialect is not specified, it is determined using sniff."""
#    pdb.set_trace()
    if dialect is None:
        dialect = sniff(filepath)
    reader = csv.reader(open(filepath,'rb'),dialect)
    i=0
    dataout = []
    for row in reader:
        if i>=startrow:
            dataout.append(row)
        if maxrows is not None:
            if i >= maxrows+startrow-1:
                break
        i+=1
    return dataout

def NestedStringListtoFloatArray(listin):
    """The csv module reads all data in as strings.  This simple
    function calls float for every entry in listin, assuming listing
    is a nested list of spreadsheet strings (i.e. a list of the rows
    of a text file)."""
    floatlist = [map(float, row) for row in listin]
    return array(floatlist)


def ArrayFromFile(filepath, startrow, dialect=None, maxrows=None):
    """Given a filepath and a startrow, this function calls getrows
    and NestedStringListtoFloatArray to read in a chunk of data and
    convert it all to an array of floating point numbers."""
    strlist = getrows(filepath, maxrows, dialect, startrow)
    return NestedStringListtoFloatArray(strlist)

def StrArrayFromFile(filepath, startrow, dialect=None, maxrows=None):
    """Given a filepath and a startrow, this function calls getrows
    to read in a chunk of data."""
    strlist = getrows(filepath, maxrows, dialect, startrow)
    return strlist

def WriteMatrixtoText(matrix, filepath, dialect, labels=None, append=False):
    if append:
        f=open(filepath,'ab')
    else:
        f=open(filepath, "wb")
    writer = csv.writer(f,dialect)
    if labels and (not append):
        writer.writerow(labels)
    writer.writerows(matrix)
    f.close()

def WriteMatrixtoTabDelim(matrix, filepath, labels=None, append=False):
    WriteMatrixtoText(matrix, filepath, tabdelim, labels)

def WriteMatrixtoCSV(matrix, filepath, labels=None, append=False):
    WriteMatrixtoText(matrix, filepath, mycsv, labels)



################################
#
#   Excel
#
################################

def search_list(mylist, ent, check_case=False, match=False):
    ind=-1
    for x, item in enumerate(mylist):
        ci = -1
        if match and check_case:
            if item == ent:
                ci=1
        elif match and not check_case:
            if item.lower() == ent.lower():
                ci = 1
        elif check_case:
            ci=item.find(ent)
        else:
            ci=item.lower().find(ent.lower())
        if ci > -1:
            ind=x
            break
    return ind

def find_first_empty_row(nested_list_in,searchcol=0):
    searchdata=zip(*nested_list_in)[searchcol]
    for n,item in enumerate(searchdata):
        if not item:
            break
    return n



def NestedListfromPath(filepath, sheet=0):
    mysheet=get_workbook_sheet(filepath, sheet)
    return ExcelSheettoNestedList(mysheet)


def float_from_string_list(listin):
    return map(float,listin)


def get_column_data(nestedlist, cols, startrow, endrow=None, makefloat=True):
    """Extract data in columns cols from nestedlist starting at
    startrow and ending at endrow.  Attempts to intelligently read all
    the data if endrow is None.

    Keep in mind that python slices stop one short of the end, so
    endrow should probably one more than you think."""
    nestout = []
    for col in cols:
        curcol = get_col(nestedlist, col)
        if endrow is not None:
            curcol = curcol[startrow:endrow]
        else:
            curcol = curcol[startrow:]
        curcol = clean_empty_string_from_bottom_of_column(curcol)
        if makefloat:
            try:
                curcol = map(float, curcol)
            except ValueError:
                curcol = curcol
        nestout.append(curcol)
    return nestout

def clean_empty_string_from_bottom_of_column(colin, copy=False):
    if copy:
        colout=copy.deepcopy(colin)
    else:
        colout=colin
    N=range(len(colout)-1,0,-1)
    for ind in N:
        if not colout[ind]:
            colout.pop(ind)
        else:
            break
    return colout



def ExcelSheettoNestedList(sheetin):
    listout=[]
    for n in range(sheetin.nrows):
        curslice=sheetin.row_slice(n)
        currow=[str(item.value) for item in curslice]
        listout.append(currow)
    return listout

def find_header_row(p, nestedlist):
    indr=-1
    p=re.compile(r'[Ff]irst.*[Nn]ame')

    for x in range(100):
	temp=str(sh1.cell_value(x,0))
	if p.match(temp):
		indr=x
		break


def search_uids(ent, uids):
    ind=-1
    ind=search_list(uids, ent)
    return ind

def load_cols(filename, collabels, headingrow=0, worksheet=0):
    """Open spreadsheet filename and search for collabels in
    headingrow.  Returna nested list of the corresponding columns."""
    wb=xlrd.open_workbook(filename)
    ws=wb.sheet_by_index(worksheet)
    headersc=ws.row_slice(headingrow)
    headers=[item.value for item in headersc]
    nestedlist=[]
    for col in collabels:
        ind=-1
        ind=search_list(headers,col)
        if ind>-1:
            curcol=ws.col_slice(ind)
        else:
            curcol=[]
        nestedlist.append(curcol)
    return nestedlist

def get_values_for_col(colin):
    colout=[item.value for item in colin]
    return colout

def get_values_nested(nestin):
    nestout=[get_values_for_col(col) for col in nestin]
    return nestout

def truncate_nested(nestin,startind=0,endind=None):
    if endind is None:
        nestout=[col[startind:] for col in nestin]
    else:
        nestout=[col[startind:endind] for col in nestin]
    return nestout

def clean_strings(colin):
    temp=colin[0]
    if type(temp)==type(u'temp') or type(temp)==type('temp'):
        colout=[str(item).strip() for item in colin]
    else:
        colout=colin
    return colout

def clean_nested(nestin):
    nestout=[clean_strings(col) for col in nestin]
    return nestout

def find_first_empty_row(col):
    end_row=None
    for x, ent in enumerate(col):
	if not ent:
            end_row=x#this is the first empty row
            break
    return end_row

def sort_col(colin, indmap):
    inds=range(len(colin))
    colout=[colin[indmap[ind]] for ind in inds]
    return colout

def sort_nested(nestin, indmap):
    nestout=[sort_col(col, indmap) for col in nestin]
    return nestout

def trunc_nested_based_on_col(nestin, trunccol):
    end_row=find_first_empty_row(trunccol)
    nestout=truncate_nested(nestin,endind=end_row)
    return nestout

## def split_names(names):
##     first_names=[]
##     last_names=[]
##     for name in names:
##         firstn,lastn=name.split(' ',1)
##         firstn=firstn.strip()
##         lastn=lastn.strip()
##         first_names.append(firstn)
##         last_names.append(lastn)
##     return first_names, last_names

