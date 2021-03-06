import sympy, copy, re

int_div_pat = re.compile('/([0-9])+([,\\]])')

def int_divide_replace(match):
    """Replace a string in generated code that might lead to a divid
    by zero.  The pattern above is based on a divide by an integer
    followed by a comma or ] as part of a matrxi.  I was seeing this
    in my autogenerated Lagrangian models for a beam broken into rigid
    masses and torsional spring dampers."""
    outstr = '/' + match.group(1) + '.0' + match.group(2)
    return outstr


def replace_int_div_list(listin):
    listout = []

    for line in listin:
        lineout = int_div_pat.sub(int_divide_replace, line)
        listout.append(lineout)

    return listout

    
def matrix_subs(matin, subs_dict):
    mat_out = copy.copy(matin)

    for i, row in enumerate(matin):
        for j, item in enumerate(row):
            item_out = item.subs(subs_dict)
            mat_out[i,j] = item_out

    return mat_out

    
def Symbol_filt(linein):
    if linein.find('Symbol(') > -1:
        return False
    else:
        return True


def _clean_Symbol(listin, make_copy=False):
    if make_copy:
        out_list = copy.copy(listin)
        out_list = filter(Symbol_filt, out_list)
    else:
        out_list = filter(Symbol_filt, listin)
    return out_list


def sage_poly_to_python_list(sage_poly, s=None):
    """Convert a symbolic sage expression in s to a python list of
    coeffs."""
    if s is None:
        s = sympy.Symbol('s')
    sympy_poly = sympy.Poly(sage_poly,s)

    if callable(sympy_poly.coeffs):
        mycoeffs = sympy_poly.coeffs()
    else:
        mycoeffs = sympy_poly.coeffs
    coeff_str = sympy.python(mycoeffs)
    coeff_list = coeff_str.split('\n')
    coeff_list = _clean_Symbol(coeff_list)
    return coeff_list


def sympy_code_list_to_useable_code(codelist, outname):
    assert len(codelist) == 1, "got codelist whose length != 1"
    codeline = codelist[0]
    assert codeline.find('e = ') == 0, "codeline doesn't start with 'e = '"
    outline = outname + ' = ' + codeline[4:]
    return outline


def sage_to_python(expr):
    sympy_expr = sympy.sympify(expr)
    python_str = sympy.python(sympy_expr)
    python_list = python_str.split('\n')
    last_line = python_list[-1]
    assert last_line.find('e = ') == 0, "last_line doesn't start with 'e = ':" + \
           last_line
    return last_line[4:]

    
    
def sage_tf_to_python(tf, name_out='tf', s=None, substr=''):
    N = tf.numerator()
    D = tf.denominator()
    Nsympy = sympy.sympify(N)
    Dsympy = sympy.sympify(D)
    num_list = sage_poly_to_python_list(Nsympy)
    den_list = sage_poly_to_python_list(Dsympy)

    num_name = 'num'
    den_name = 'den'
    if substr:
        num_name += '_' + substr
        den_name += '_' + substr
        name_out += '_' + substr

    out_list = []
    num_line = sympy_code_list_to_useable_code(num_list, num_name)
    out_list.append(num_line)
    den_line = sympy_code_list_to_useable_code(den_list, den_name)
    out_list.append(den_line)
    last_line = '%s = controls.TF(%s, %s)' % (name_out, num_name, den_name)
    out_list.append(last_line)
    
    return out_list


def sage_matrix_to_nested_list(matin):
    nr, nc = matin.dimensions()

    nested_list = []

    for i in range(nr):
        row_list = []
        for j in range(nc):
            ent = matin[i,j]
            row_list.append(ent)

        nested_list.append(row_list)

    return nested_list


def sage_matrix_to_sympy(matin):
    matlist = sage_matrix_to_nested_list(matin)
    sympy_list = []
    for row in matlist:
        row_list = []
        for item in row:
            item_s = sympy.sympify(item)
            row_list.append(item_s)
        sympy_list.append(row_list)
    matout = sympy.Matrix(sympy_list)
    return matout


def cse_sage_sympy(exprs):
    if str(type(exprs)).find("'sympy.") > -1:
        expr_sympy = exprs
    elif (exprs.__class__.__name__ == 'Matrix_symbolic_dense') and \
             (str(type(exprs)).find("'sage.") > -1):
        expr_sympy = sage_matrix_to_sympy(exprs)
    else:
        expr_sympy = sympy.sympify(exprs)
    out = sympy.cse(expr_sympy)
    outlist = []
    mydefs = out[0]
    result_list = out[1]
    def_fmt = '%s = %s'
    for line in mydefs:
        line_str = def_fmt % line
        outlist.append(line_str)

    assert len(out[1]) == 1, "Not sure what to do with cse results that do not have exactly one results: %s" % out[1]
    out_str = 'result = %s' % out[1][0]
    outlist.append(out_str)
    outlist.append('return result')
    return outlist


def _row_to_str(row_in):
    str_out = ''

    for elem in row_in:
        cur_str = sage_to_python(elem)
        if str_out:
            str_out += ', '
        str_out += cur_str

    return '[' + str_out + ']'


def sage_array_to_python(array_in, lhs_str):
    """Convert an array that is most likely a numpy object array to
    code that can be pasted into a numeric python file to be executed
    in IPython.  Basically, use the sage_to_python function above on
    each element of the array and output a string."""

    if len(array_in.shape) == 1:
        #1D array
        row_str = _row_to_str(array_in)
        out_str = '%s = array(%s)' % (lhs_str, row_str)
        return [out_str]

    #if we get to this point, we have a 2D array
    line1 = '%s = array([\\' % lhs_str

    ws = ' '*(len(line1) -1)
    mylines = [line1]

    nr, nc = array_in.shape
    
    for i, row in enumerate(array_in):
        row_str = _row_to_str(row)
        if i < (nr-1):
            #append comma to all rows except the last one
            row_str += ','
        mylines.append(ws + row_str + ' \\')

    mylines.append(ws[0:-1] + '])')

    return mylines
