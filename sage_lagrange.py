import numpy

def eom_one_row(T,V,W,vel,pos,diff_dict):
    """Find one row of the equations of motion by taking
    d/dt(T.diff(vel)) and dV/d pos and d W/d pos"""
    Tpartdiff = T.diff(vel)
    Tderiv = Tpartdiff.subs(diff_dict)#<-- eventually I should really
                                      #do d/dt here, for now I am just
                                      #substituteing theta_ddot for
                                      #theta_dot and such
    Vderiv = V.diff(pos)
    Q = W.diff(pos)
    eom = Tderiv+Vderiv == Q
    return eom



def eom_to_one_matrix_row(eom, vars, zero_dict):
    """Find one row of one matrix of the EOMs by substituting that one
    var in vars is one and the rest are zero.  Do this for each
    element in vars to find the row so that when we are done
    dot(mat,vars) would reproduce eom.

    zeros_dict is used to zero out all velocity, position, and
    acceleration variable after 1 is subsituted for the element who
    coefficient we are trying to find.  So, zero_dict must have keys
    for every position, velocity, and acceleration variable and all
    the values must be 0."""
    n = len(vars)
    row_list = []
    lhs = eom.lhs()
    
    for var in vars:
        seom = lhs.subs({var:1})
        seom = seom.subs(zero_dict)
        row_list.append(seom)

    return numpy.array(row_list)


def eom_list_to_one_matrix(eom_list, vars, zero_dict):
    big_list = []

    for eom in eom_list:
        cur_row = eom_to_one_matrix_row(eom, vars, zero_dict)
        big_list.append(cur_row)

    return numpy.array(big_list)


    
def find_eom_list(T,V,W,list_of_tuples,diff_dict):
    """Create a list of equations of motion where each eom is one row
    of Lagrange's equations.  list_of_tuples should contain a list of
    (vel,pos) tuples that correspond to q and q_dot for each row."""
    eom_list = []

    for vel, pos in list_of_tuples:
        eom_i = eom_one_row(T,V,W,vel,pos,diff_dict)
        eom_list.append(eom_i)

    return eom_list
