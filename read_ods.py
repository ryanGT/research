#!/usr/bin/env python

# import the needed modules
import zipfile
import xml.parsers.expat

from IPython.core.debugger import Pdb
import copy
# get content xml data from OpenDocument file
mypath = '/home/ryan/siue/classes/484/2010/final_report_grades_and_assessment.ods'
#mypath = 'test_sheet.ods'

class Element(list):
    def __init__(self, name, attrs):
        self.name = name
        self.attrs = attrs

class TreeBuilder:
    def __init__(self):
        self.root = Element("root", None)
        self.path = [self.root]

    def start_element(self, name, attrs):
        element = Element(name, attrs)
        self.path[-1].append(element)
        self.path.append(element)

    def end_element(self, name):
        assert name == self.path[-1].name
        self.path.pop()

    def char_data(self, data):
        self.path[-1].append(data)



def showtree(node, prefix=""):
    print prefix, node.name
    for e in node:
        if isinstance(e, Element):
            showtree(e, prefix + "  ")
        else:
            print prefix + "  ", e


def find_body(root):
    N1 = len(root)
    found_content = 0
    for i in range(N1):
        if root[i].name == 'office:document-content':
            found_content = 1
            content = root[i]
    assert found_content, 'Did not find an element in root with name office:document-content'


    N2 = len(content)
    found_body = 0
    for i in range(N2):
        if content[i].name == 'office:body':
            found_body = 1
            body = content[i]
    assert found_body, 'Did not find an element in content with name office:body'
    return body


def find_spreadsheet(body):
    N1 = len(body)
    found_sheet = 0
    for i in range(N1):
        if body[i].name == 'office:spreadsheet':
            found_sheet = 1
            sheet = body[i]
    assert found_sheet, 'Did not find an element in body with name office:spreadsheet'
    return sheet

def nonempty_row_or_col(row_or_col):
    if row_or_col == []:
        return False
    elif row_or_col == [[]]:
        return False
    else:
        return True


def nonempty_table(table):
    found_any = False
    for elem in table:
        if nonempty_row_or_col(elem):
            found_any = True
            break
    return found_any

def nonempty_first_cell(row):
    return bool(row[0])

#showtree(treebuilder.root)

def process_nested(nested_in):
    """Take a nested list of table data and create a nested list of
    strings and floats."""
    for i, row in enumerate(nested_in):
        for j, elem in enumerate(row):
            value = elem[0][0].encode()
            if elem.attrs['office:value-type'] == 'float':
                value = float(value)
            if j == 0:
                row_out = [value]
            else:
                row_out.append(value)
        if i == 0:
            nested_out = [row_out]
        else:
            nested_out.append(row_out)
    return nested_out


def read_ods(pathin):
    #open zip archive
    ziparchive = zipfile.ZipFile(pathin, "r")
    xmldata = ziparchive.read("content.xml")
    ziparchive.close()

    # create parser and parsehandler
    parser = xml.parsers.expat.ParserCreate()
    treebuilder = TreeBuilder()
    # assign the handler functions
    parser.StartElementHandler  = treebuilder.start_element
    parser.EndElementHandler    = treebuilder.end_element
    parser.CharacterDataHandler = treebuilder.char_data

    # parse the data
    parser.Parse(xmldata, True)

    b = find_body(treebuilder.root)
    s = find_spreadsheet(b)
    mytables = filter(nonempty_table, s)
    assert len(mytables) == 1, 'Did not find exactly one non-empty table in the sheet, len(s) = '+str(len(s))
    mytable = mytables[0]
    alldata = filter(nonempty_row_or_col, mytable)
    cleandata = filter(nonempty_first_cell, alldata)
    processed_data = process_nested(cleandata)
    return processed_data


mydata = read_ods(mypath)
#there are 3 tables in the sheet s, two of them are empty, but proving
#that isn't simple:

## In [128]: s
## Out[128]:
## [[[],
##   [[[u'Col1']], [[u'Col2']], [[u'Col3']], [[u'Col4']]],
##   [[[u'team1']], [[u'1']], [[u'2']], [[u'3']]],
##   [[[u'team2']], [[u'3']], [[u'4']], [[u'7']]],
##   [[[u'team2']], [[u'5']], [[u'6']], [[u'5.5']]],
##   [[], [[u'3']], [[u'4']], [[u'5.17']]]],
##  [[], [[]]],
##  [[], [[]]]]

#the empty tables have an empty column that is a single list and an
#empty row that is a nested list:

# In [129]: len(s)
# Out[129]: 3

## In [130]: s[1]
## Out[130]: [[], [[]]]

## In [131]: s[2]
## Out[131]: [[], [[]]]

#I need a function to tell if a table is empty to not.  I could pop
#empty lists, but the bool of an empty nested list is trickier.
