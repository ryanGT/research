import numpy

def expr_to_dttmm_row(expr, z_list, zero_dict):
    """Find one row of a DT-TMM transfer matrix corresponding to expr
    by substituting one var in vars as one and the rest are zero.
    Do this for each element in vars to find the row so that when we
    are done dot(mat,vars) would reproduce expr.

    This is for use in deriving DT-TMM transfer matrices.
    """
    remainder = expr.subs(zero_dict)
    
    n = len(z_list)
    row_list = []
    
    for var in z_list:
        seom = expr.subs({var:1}) - remainder
        seom = seom.subs(zero_dict)
        row_list.append(seom)

    row_list.append(remainder)
    return numpy.array(row_list)

