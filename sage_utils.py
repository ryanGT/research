import sympy, copy

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
    
    coeff_str = sympy.python(sympy_poly.coeffs)
    coeff_list = coeff_str.split('\n')
    coeff_list = _clean_Symbol(coeff_list)
    return coeff_list


def sympy_code_list_to_useable_code(codelist, outname):
    assert len(codelist) == 1, "got codelist whose length != 1"
    codeline = codelist[0]
    assert codeline.find('e = ') == 0, "codeline doesn't start with 'e = '"
    outline = outname + ' = ' + codeline[4:]
    return outline

    
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


def cse_sage_sympy(exprs):
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
