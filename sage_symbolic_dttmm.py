"""This module will be used from within sage to generate symbolic
expressions for various aspects of DT-TMM analysis.  There will be a
system class, :py:class:`symbolic_dttmm_system` and at least two
element classes, :py:class:`symbolic_rigid_mass_element` and
:py:class:`symbolic_tsd_element`."""
from sage.all import *

from numpy import array

import sage_utils, os, sympy

class symbolic_dttmm_system(object):
    def __init__(self):
        A, D = var('A','D')
        self.A = A
        self.D = D


    def assign_element_list(self, element_list):
        self.element_list = element_list


    def _assign_parent(self):
        for element in self.element_list:
            element.parent = self


    def _assign_element_indices(self):
        for i, element in enumerate(self.element_list):
            element.i = i

            
    def _assign_prev_elements(self):
        for i in range(1,len(self.element_list)):
            self.element_list[i].prev_element = self.element_list[i-1]


    def prep(self):
        self._assign_parent()
        self._assign_element_indices()
        self._assign_prev_elements()

        for element in self.element_list:
            element._create_and_assign_parameters()


    def find_Ui(self):
        Ui_list = []

        for element in self.element_list:
            Ui = element.Get_Mat()
            Ui_list.append(Ui)

        self.Ui_list = Ui_list
        return Ui_list


    def multiply_Ui(self):
        U0i_list = []

        first = 1

        for Ui in self.Ui_list:
            if first:
                U0i = Ui
                first = 0
            else:
                temp = Ui*U0i
                U0i = temp
                
            U0i_list.append(U0i)

        self.U_sys = U0i#last one
        self.U0i_list = U0i_list
            

    def solve_boundary_conditions(self):
        raise NotImplementedError


    def find_zi(self):
        """Find symbolic expressions for the state vector at each
        location in the model."""
        zi_list = []

        for U0i in self.U0i_list:
            zi = U0i*self.z0
            zi_list.append(zi)

        self.zi_list = zi_list


    def gen_list_of_all_params(self):
        param_list = ['A','D']

        for element in self.element_list:
            cur_params = element._gen_list_of_element_params()
            param_list.extend(cur_params)

        self.all_params = param_list
        return param_list


    def gen_assignment_lines(self, input_dict_name='params'):
        """Create lines of code for assigning the variables within a
        function from params or whatever the input dictionary is
        called."""

        all_params = self.gen_list_of_all_params()

        outlines = []

        fmt = '%s = ' + input_dict_name + '["%s"]'
        
        for param in all_params:
            curline = fmt % (param, param)
            outlines.append(curline)

        return outlines


    def gen_z_mat(self):
        """Create a nested list where the first row contains the x for
        each station and the secone row contains the thetas.  Then
        pass the nested list to sympy.Matrix."""
        x_list = []
        theta_list = []
        for z_i in self.zi_list:
            x_i = sympy.sympify(z_i[0])
            th_i = sympy.sympify(z_i[1])
            x_list.append(x_i)
            theta_list.append(th_i)

        mat_list = [x_list, theta_list]
        self.z_mat = sympy.Matrix(mat_list)
        return self.z_mat
    
        
    def gen_all_module_lines(self, input_dict_name='params', filename=None):
        outlines = ['from __future__ import division', \
                    'from numpy import *', \
                    '']

        ## func_fmt = 'def z%i(' + input_dict_name + '):'
        func_line = 'def z_mat(%s):' % input_dict_name
        ws = ' '*4
        assignment_lines_raw = self.gen_assignment_lines()
        assignment_lines = [ws + line for line in assignment_lines_raw]

        func_lines = [func_line]
        func_lines.extend(assignment_lines)
        self.gen_z_mat()
        cse_lines = sage_utils.cse_sage_sympy(self.z_mat)
        cse_lines_clean = [ws + line.replace('Matrix([','array([') for line in cse_lines]
        func_lines.extend(cse_lines_clean)
        outlines.extend(func_lines)


        ## for i, z in enumerate(self.zi_list):
        ##     j = i + 1
        ##     func_line = func_fmt % j
        ##     func_lines = [func_line]
        ##     func_lines.extend(assignment_lines)
        ##     self.gen_z_mat()
        ##     #cse_lines = sage_utils.cse_sage_sympy(z)
        ##     cse_lines = sage_utils.cse_sage_sympy(self.z_mat)
        ##     cse_lines_clean = [ws + line.replace('Matrix([','array([') for line in cse_lines]
        ##     func_lines.extend(cse_lines_clean)
        ##     func_lines.append('')
        ##     func_lines.append('')

        ##     outlines.extend(func_lines)

        while not outlines[-1]:
            outlines.pop(-1)


        if filename is not None:
            import txt_mixin
            txt_mixin.dump(filename, outlines)
            
        return outlines
                               
            
        

class symbolic_dttmm_system_clamped_free(symbolic_dttmm_system):
    def solve_boundary_conditions(self):
        submat = self.U_sys[2:4, 2:4]
        last_col = self.U_sys[2:4,4]
        submat_inv = submat.inverse()
        MV_base = -submat_inv*last_col
        #z0 = matrix(SR, [[0.0],[0.0],[MV_base[0,0]], [MV_base[1,0]],[1.0]])
        z0 = matrix(SR, [[0.0],[0.0],[MV_base[0,0]], [MV_base[1,0]],[1.0]])
        self.z0 = z0
        return self.z0

        

common_params = ['Bx','Bth','Ex','Eth']
rigid_param_strings = ['m','L','I','r','F'] + common_params
tsd_param_strings = ['k','b'] + common_params

    
class symbolic_rigid_mass_element(object):
    def __init__(self, i=None, parent=None, prev_element=None):
        """i is the index for the element.  It will be used to ensure
        unique names for the elements parameters.

        I am allowing a default of None so that a symbolic DT-TMM
        system could assign the indices later.  But almost nothing
        will work before the index i is assigned a value.

        parent will be a DT-TMM symbolic system that will at least
        contain the symbolic A and D """
        self.i = i
        self.parent = parent
        self.prev_element = prev_element
        self.param_strings = rigid_param_strings


    def _check_i(self):
        assert self.i is not None, "You must assign an index number i to each element before calling almost any other method."


    def _gen_list_of_element_params(self):
        self._check_i()
        int_str = '%i' % self.i
        self.my_param_strings = [item + int_str for item in self.param_strings]
        return self.my_param_strings


    def _create_and_assign_parameters(self):
        my_param_strings = self._gen_list_of_element_params()
        lhs = ', '.join(my_param_strings)
        quoted_params = ['"%s"' % item for item in my_param_strings]
        rhs_p1 = ', '.join(quoted_params)
        exec_line = '%s = var(%s)' % (lhs, rhs_p1)
        exec(exec_line)
        for attr, var in zip(self.param_strings, my_param_strings):
            setattr(self, attr, eval(var))
        return exec_line


    def Get_Mat(self):
        #system parameter
        A = self.parent.A

        #previous element parameter
        if self.prev_element is None:
            Bxp = 0
        else:
            Bxp = self.prev_element.Bx

        #num. int. parameter
        Bth = self.Bth

        #my parameters
        m = self.m
        I = self.I
        L = self.L
        r = self.r
        F = self.F
        
        ## U_RM = array([\
        ##              [1.0, L, 0, 0, 0], \
        ##              [0, 1, 0, 0, 0], \
        ##              [-A*L*m + A*m*r, A*I - A*L*m*r + A*m*r**2, 1, -L, Bth*I - Bth*L*m*r + Bth*m*r**2 - Bxp*L*m + Bxp*m*r], \
        ##              [A*m, A*m*r, 0, 1, Bth*m*r + Bxp*m - F], \
        ##              [0, 0, 0, 0, 1] \
        ##             ])
        mylist = [\
                 [1.0, L, 0, 0, 0], \
                 [0, 1, 0, 0, 0], \
                 [-A*L*m + A*m*r, A*I - A*L*m*r + A*m*r**2, 1, -L, Bth*I - Bth*L*m*r + Bth*m*r**2 - Bxp*L*m + Bxp*m*r], \
                 [A*m, A*m*r, 0, 1, Bth*m*r + Bxp*m - F], \
                 [0, 0, 0, 0, 1] \
                 ]

        U_RM = matrix(SR,mylist)
        self.U = U_RM
        return U_RM


#tsd_param_strings = ['Bx','Bth','k','b']

class symbolic_tsd_element(symbolic_rigid_mass_element):
    def __init__(self, i=None, parent=None, prev_element=None):
        symbolic_rigid_mass_element.__init__(self, i=i, \
                                             parent=parent, \
                                             prev_element=prev_element)
        self.param_strings = tsd_param_strings
        

    def Get_Mat(self):
        #system parameter
        D = self.parent.D

        #my parameters
        k = self.k
        b = self.b
        Eth = self.Eth
        
        den = (k + b*D)
        
        if self.prev_element is None:
            E_prev = 0
        else:
            E_prev = self.prev_element.Eth


        ## self.U = array([[1.0, 0.0, 0.0, 0.0, 0.0], \
        ##                 [0.0, 1.0,  1.0/den,  0.0,\
        ##                  -b*(self.Eth - E_prev)/den], \
        ##                 [0.0, 0.0, 1.0, 0.0, 0.0], \
        ##                 [0.0, 0.0, 0.0, 1.0, 0.0], \
        ##                 [0.0, 0.0, 0.0, 0.0, 1.0]])
        mylist = [[1.0, 0.0, 0.0, 0.0, 0.0], \
                  [0.0, 1.0,  1.0/den,  0.0,\
                   -b*(self.Eth - E_prev)/den], \
                  [0.0, 0.0, 1.0, 0.0, 0.0], \
                  [0.0, 0.0, 0.0, 1.0, 0.0], \
                  [0.0, 0.0, 0.0, 0.0, 1.0]]
        U_TSD = matrix(SR,mylist)
        self.U = U_TSD
        return U_TSD



if __name__ == '__main__':
    mysys = symbolic_dttmm_system_clamped_free()
    mylist = []

    N = 10
    
    for i in range(N):
        tsd_i = symbolic_tsd_element()
        mylist.append(tsd_i)
        rm_i = symbolic_rigid_mass_element()
        mylist.append(rm_i)


    mysys.assign_element_list(mylist)
    mysys.prep()
    mysys.find_Ui()
    mysys.multiply_Ui()
    z0 = mysys.solve_boundary_conditions()
    mysys.find_zi()

    import rwkos
    outfolder = rwkos.FindFullPath('siue/Research/work/2013/sabbatical_summer_fall_2013/DTTMM_vs_FEA_vs_Lagrange_TSDs')
    modname = 'sym_tsd_rigid_mass_dttmm_N_%i.py' % N
    outpath = os.path.join(outfolder, modname)
    #mysys.gen_all_module_lines(filename=outpath)
