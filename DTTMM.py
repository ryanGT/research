"""This module is a first attempt to create an object-oriented
framework for DT-TMM analysis."""
# Authors: Ryan Krauss and Mohamed Okasha
# Started: 07/20/2011
# Version: 0.1
from scipy import *
from scipy import optimize
import numpy

import time
import rwkmisc

beta = 1.0/6.0
gamma = 0.5
theta_W = 1.5

import pdb
from IPython.core.debugger import Pdb

class DT_TMM_System(object):
    """A DT-TMM System will consist of a connection of elements and a
    specification of boundary conditions.  """
    def __init__(self, element_list, N_states, \
                 N=None, dt=None, int_case=1, \
                 initial_conditions=None, \
                 unknown_params=[], \
                 param_limits={}):
        self.element_list = element_list
        self.N_states = N_states
        self.N = N
        self.dt = dt
        self.int_case = int_case
        self.initial_conditions = initial_conditions
        self.unknown_params = unknown_params
        self.param_limits = param_limits


    def _initialize_vectors(self, N):
        for element in self.element_list:
            element._initialize_vectors(N)
            
        
    def calculate_ABDE(self, i, dt, int_case=2, **kwargs):
        for element in self.element_list:
            element.calculate_ABDE(i, dt, int_case=int_case, **kwargs)

            
    def calculate_transfer_matrices(self, i):
        for element in self.element_list:
            element.calculate_transfer_matrix(i)#must set self.U


    def calculate_system_transfer_matrix(self, i=None):
        #I don't think we want to keep the transfer matrices for each
        #time step, so I am ignoring i for now
        U = eye(self.N_states + 1)
        for element in self.element_list:
            U_j = element.U
            U = dot(U_j, U)
        self.U_sys = U


    def calculate_state_vectors(self, i):
        prev_z = self.z0
        for element in self.element_list:
            cur_z = element.calculate_state_vector(prev_z, i)#z_j = dot(U_j, z_(j-1)),
                                                             #store x[i] for each element
            prev_z = cur_z#set up for next element


    def calculate_velocity_and_accel(self, i):
        for element in self.element_list:
            element.calculate_velocity_and_accel(i)
            

    def Run_Simulation(self, N, dt, int_case=2):
        self._initialize_vectors(N)
        for i in range(1,N):    # Time loop
            self.calculate_ABDE(i, dt, int_case=int_case)
            self.calculate_transfer_matrices(i)
            self.calculate_system_transfer_matrix(i)
            self.solve_boundary_conditions(i)#must set self.z0
            self.calculate_state_vectors(i)
            self.calculate_velocity_and_accel(i)


    def solve_boundary_conditions(self, i=None):
        raise NotImplementedError


    def get_param_by_name(self, name):
        """Assuming name is element_name.param_name, find the index of
        element_name and then call getattr for that element."""
        elem_name, param_name = name.split('.',1)
        elem = self._get_element_by_name(elem_name)
        param = getattr(elem, param_name)
        return param


    def set_param_by_name(self, name, param):
        """Assuming name is element_name.param_name, find the index of
        element_name and then call setattr for that element."""
        elem_name, param_name = name.split('.',1)
        elem = self._get_element_by_name(elem_name)
        setattr(elem, param_name, param)
        
    

    def get_ig(self):
        """Return the initial guess vector for use with optimize.fmin
        system identification by using getattr for the parameters in
        self.unknown_params."""
        N = len(self.unknown_params)
        C_vect = zeros(N)
        for i, attr in enumerate(self.unknown_params):
            C_vect[i] = self.get_param_by_name(attr)

        return C_vect


    def _find_element_index_by_name(self, name):
        """This is a helper method for _find_param_element_ids that
        searches the elements in self.element_list to find one whose
        name matches name."""
        index = -1
        for i, elem in enumerate(self.element_list):
            if hasattr(elem, 'name'):
                if elem.name == name:
                    index = i
                    
        assert index != -1, "did not find element with name %s" % name
        
        return index


    def _get_element_by_name(self, name):
        ind = self._find_element_index_by_name(name)
        return self.element_list[ind]
            

    def _find_param_element_ids(self):
        """The keys of self.unknown_params are
        element_name.param_name.  In order to get or set the
        unknown_params efficiently, we need a list of the element
        indices."""
        element_ids = []
        self._find_element_and_parameter_names()
        for elem_name in self.unknown_element_names:
            ind = self._find_element_index_by_name(elem_name)
            element_ids.append(ind)

        self.unknown_element_ids = element_ids
        return element_ids


    def _find_element_and_parameter_names(self):
        """Separate the names of the unknown_params into the element
        names and parameter names, i.e. split the unknown_params at
        the dot: element_name.param_name"""
        element_names = []
        param_names = []
        
        for key in self.unknown_params:
            elem_name, param_name = key.split('.',1)
            element_names.append(elem_name)
            param_names.append(param_name)

        self.unknown_element_names = element_names
        self.unknown_parameter_names = param_names
    

    def set_params(self, C):
        """Set the parameters listed in self.unknown_params with the
        elements of C, which are assumed to be in the same order as
        self.unknown_params."""
        if not hasattr(self, 'unknown_element_ids'):
            self._find_param_element_ids()
        for C_i, ind, attr in zip(C, self.unknown_element_ids, \
                                  self.unknown_parameter_names):
            elem = self.element_list[ind]
            setattr(elem, attr, C_i)


    def set_params_from_dict(self, mydict):
        for attr, value in mydict.iteritems():
            self.set_param_by_name(attr, value)
        
            
        
    def build_fit_res_dict(self):
        """Create dictionary with self.unknown_params as the keys and
        the current parameter values from getattr as the values."""
        mydict = {}
        for attr in self.unknown_params:
            val = self.get_param_by_name(attr)
            mydict[attr] = val
        return mydict
    
        
    def save_fit_res_to_dict(self, filename):
        mydict = self.build_fit_res_dict()
        rwkmisc.SavePickle(mydict, filename)


    def load_fit_res_dict(self, filename):
        mydict = rwkmisc.LoadPickle(filename)
        self.set_params_from_dict(mydict)
    
        
    def check_limits(self, C):
        """Used to apply penalties to a cost in a cost function if a
        parameter is outside the limits specified in
        self.param_limits.  The elements of C are assumed to be in the
        same order as self.unknown_params.  The elements of
        self.unknown_params are the keys used to look up the limits in
        the dictionary self.param_limits.  If self.param_limits does
        not have the corresponding key, no penalty is applied for that
        parameter."""
        penalty = 1.0e4
        out = 0.0

        for C_i, attr in zip(C, self.unknown_params):
            if self.param_limits.has_key(attr):
                limits = self.param_limits[attr]
                if C_i < limits[0]:
                    out += penalty
                elif C_i > limits[1]:
                    out += penalty
            elif C_i < 0.0:
                #penalize for negative parameters even if they aren't
                #in the param_limits dictionary
                out += penalty

        return out


    def mymodel(self, C):
        raise NotImplementedError


    def mycost(self, C):
        raise NotImplementedError
        

    def run_fit(self):
        ig = self.get_ig()
        C_opt = optimize.fmin(self.mycost, ig)
        self.set_params(C_opt)
        return C_opt


## def set_params2(C):
##     base_mass.m = C[0]
##     base_mass.I = C[1]
##     accel_mass.I = C[2]

## def mymodel(C):
##     set_params2(C)
##     sys.Run_Simulation(N, dt, initial_conditions=initial_conditions, \
##                        int_case=1)

##     #x_tip_DTTMM = last_mass_DTTMM.x
##     #a_tip_DTTMM = last_mass_DTTMM.xddot
##     #x_tip_DTTMM = accel_mass.x
##     a_tip_DTTMM = accel_mass.xddot
##     theta_mass0 = base_mass.theta*mydict['H']
##     return theta_mass0, a_tip_DTTMM

## ## ig = array([k_spring, c_spring, \
## ##             k_clamp, c_clamp])

## ig = array([base_mass.m, base_mass.I, accel_mass.I])

## plot_ig = 1

## if plot_ig:
##     #ig = array([  1.00000000e+01, 1.00000000e-03, 3.46849426e+00, 8.07734624e-12])

##     theta_ic, a_tip_DTTMM_ic = mymodel(ig)

##     plot(t, theta_ic, ':', linewidth=2)
##     plot(t, a_tip_DTTMM_ic, '-.', linewidth=2)


## myfunc = mymodel

## def mycost(C):
##     theta_DTTMM, a_tip_DTTMM = myfunc(C)
##     e_theta_vect = theta_exp - theta_DTTMM
##     e_theta = ((e_theta_vect)**2).sum()
##     e_a_vect = a_exp - a_tip_DTTMM
##     e_a = ((e_a_vect)**2).sum()

##     penalty = 0.0
##     if (C < 0.0).any():
##         penalty = 1.0e6
##     return e_theta + e_a + penalty



## #first opt result:
## #C_opt = array([  3.28962591e-01,   2.02157293e-01,   3.46849426e+00, 8.07734624e-12])

## # opt results with int_case = 1:
## #C_opt = array([  8.96725909e-01,   3.63151344e-02,   4.00905919e+00,
## #                 8.12525142e-14])

## run_fit = 1

## if run_fit:
##     C_opt = optimize.fmin(mycost, ig, maxiter=1000, maxfun=10000)

##     theta_opt, a_tip_DTTMM_opt = mymodel(C_opt)
        


class DT_TMM_System_clamped_free_two_states(DT_TMM_System):
    """This class models a DT-TMM System with clamped/free boundary
    conditions and two states (displacement and force - i.e. 1 DOF
    mass/spring/damper carts vibrating sideways)."""
    def solve_boundary_conditions(self, i=None):
        F_base = -self.U_sys[1,2]/self.U_sys[1,1]
        z0 = array([[0.0],[F_base],[1.0]])
        self.z0 = z0
        


class DT_TMM_Element(object):
    def __init__(self, **kwargs):
        for key, value in kwargs.iteritems():
            setattr(self, key, value)
            

    def _get_ax(self, fig=None, fignum=1, clear=False, subplot=111):
        if fig is None:
            import pylab
            fig = pylab.figure(fignum)
        if clear:
            fig.clf()
            ax = fig.add_subplot(subplot)
        else:
            ax = fig.axes[0]
        return ax


    def plot_attrs(self, attr_list, fig=None, fignum=1, clear=True, \
                   labels=None, **kwargs):
        first = True
        if labels is None:
            labels = attr_list
        ax = self._get_ax(fig=fig, fignum=fignum, clear=clear, **kwargs)

        for attr, label in zip(attr_list, labels):
            myvect = getattr(self, attr)
            ax.plot(self.t, myvect, label=label)
            


    def _initialize_vectors(self, N):
        self.x = zeros(N)
        self.xdot = zeros(N)
        self.xddot = zeros(N)
        
        

    def calculate_velocity_and_accel(self, i):
        #to find x at every station
        self.xdot[i] = self.D*self.x[i] + self.E
        self.xddot[i] = self.A*self.x[i] + self.B
        

    def calculate_ABDE(self, i, dt, int_case=2):
        #print('int_case = %i' % int_case)
        if int_case == 1:
            #Finite difference 
            A = 4.0/(dt**2)
            B = -A*(self.x[i-1] + dt*self.xdot[i-1] + (dt**2/4.0)*self.xddot[i-1])
            D = 2.0/dt
            E = -D*(self.x[i-1] + 0.5*dt*self.xdot[i-1])
        elif int_case == 2:
            #Newmark Beta
            A = 1.0/(beta*dt**2)
            B = -1.0/(beta*dt**2)*(self.x[i-1] + dt*self.xdot[i-1] + (0.5-beta)*dt**2*self.xddot[i-1])
            D = gamma/(beta*dt)
            E = self.xdot[i-1] + dt*((1.0-gamma)*self.xddot[i-1]+gamma*B)      
        elif int_case== 3:
            #Fox-Euler
            print('Fox-Euler is badddddd!!!!!!!!!!!!!!!')
            A=2.0/dt**2
            B=-A*(self.x[i-1]+dt*self.xdot[i-1])
            D=2.0/dt
            E=-(D*self.x[i-1]+self.xdot[i-1])
        self.A = A
        self.B = B
        self.D = D
        self.E = E


    def calculate_transfer_matrix(self, i=None):
        #This method must be overwritten for any derived element class
        #to be valid.  It must set self.U
        raise NotImplementedError


    def calculate_state_vector(self, prev_z, i, xdof=0):
        self.z = dot(self.U, prev_z)
        self.x[i] = self.z[xdof]
        return self.z



class DT_TMM_mass_2_states(DT_TMM_Element):
    def __init__(self, m, f=None, **kwargs):
        """f is the applied external force acting on the mass.  The f
        vector can be passed in here.  If it is, f must be left as
        None in calculate_transfer_matrix.  Then f[i] will be applied
        to the mass at each step in the simulation."""
        DT_TMM_Element.__init__(self, **kwargs)
        self.m = m
        self.f = f


    def calculate_transfer_matrix(self, i=None, f=None):
        """If f is passed in here, it overrides self.f.  If f is a
        vector, f[i] is applied to the mass (then i must not be None)."""
        if f is None:
            f = self.f
        if f is None:
            f_i = 0.0
        elif (not isscalar(f)):
            f_i = f[i]
        else:
            f_i = f
        self.U = array([[1.0,0.0,0.0], \
                        [self.m*self.A,1.0,self.m*self.B - f_i], \
                        [0.0,0.0,1.0]])
        return self.U




class DT_TMM_spring_2_states(DT_TMM_Element):
    def __init__(self, k, b, prev_element=None, **kwargs):
        """A spring/damper element.

        k is the stiffness coefficient (N/m).  b is the damping
        coefficient (N*s/m).

        The spring needs the displacement and velocity of the previous
        element to calculate the relative displacement and velocity of
        its spring and damper.  If prev_element is None, it is assumed
        the spring is connected to a wall and the D and E parameters
        of the previous element are zero."""
        DT_TMM_Element.__init__(self, **kwargs)
        self.k = k
        self.b = b
        self.prev_element = prev_element
    

    def calculate_transfer_matrix(self, i=None):
        den = (self.k+self.b*self.D)
        if self.prev_element is None:
            D_prev = 0.0
            E_prev = 0.0
        else:
            D_prev = self.prev_element.D
            E_prev = self.prev_element.E
            
        self.U = array([[(self.k+self.b*D_prev)/den, \
                         1.0/den, -(self.b*(self.E-E_prev))/den], \
                        [0.0, 1.0 ,0.0], \
                        [0.0, 0.0, 1.0]])
        return self.U



class DT_TMM_mass_spring_2_states(DT_TMM_spring_2_states):
    def __init__(self, m, f, k, b, prev_element=None, **kwargs):
        """An element that represents a spring/damper with a mass
        attached to it right side.

        m is the mass of the element
        f is the force applied to the element

        k is the stiffness coefficient (N/m).  b is the damping
        coefficient (N*s/m).

        The spring needs the displacement and velocity of the previous
        element to calculate the relative displacement and velocity of
        its spring and damper.  If prev_element is None, it is assumed
        the spring is connected to a wall and the D and E parameters
        of the previous element are zero."""
        DT_TMM_Element.__init__(self, **kwargs)
        self.m = m
        self.f = f
        self.k = k
        self.b = b
        self.prev_element = prev_element

        
    def calculate_transfer_matrix(self, i=None, f=None):
        if f is None:
            f = self.f
        if f is None:
            f_i = 0.0
        elif (not isscalar(f)):
            f_i = f[i]
        else:
            f_i = f
            
        sden = (self.k+self.b*self.D)
        
        if self.prev_element is None:
            D_prev = 0.0
            E_prev = 0.0
        else:
            D_prev = self.prev_element.D
            E_prev = self.prev_element.E

        snum11 = self.k+self.b*D_prev
        snum13 = -(self.b*(self.E-E_prev))
        
        A = self.A
        B = self.B
        m = self.m
        
        self.U = array([[snum11/sden, 1.0/sden, snum13/sden], \
                        [A*m*snum11/sden, A*m/sden + 1, A*m*snum13/sden + B*m - f],
                        [0.0, 0.0, 1.0]])
        return self.U


###############################################
#
#  Four States:
#
#    z = [x, theta, M, V, 1]
#
###############################################


class DT_TMM_Element_4_states(DT_TMM_Element):
    def _initialize_vectors(self, N):
        self.x = zeros(N)
        self.xdot = zeros(N)
        self.xddot = zeros(N)
        self.theta = zeros(N)
        self.thetadot = zeros(N)
        self.thetaddot = zeros(N)
        self.z_mat = zeros((N,4))
        #self.Ax_vect = zeros(N)
        self.Bx_vect = zeros(N)
        self.Dx_vect = zeros(N)
        self.Ex_vect = zeros(N)
        #self.Ath_vect = zeros(N)
        self.Bth_vect = zeros(N)
        self.Dth_vect = zeros(N)
        self.Eth_vect = zeros(N)
        self.int_det_x = zeros(N)
        self.int_det_th = zeros(N)
        self.int_cond_x = zeros(N)
        self.int_cond_th = zeros(N)
        
        
    def calculate_velocity_and_accel(self, i):
        #to find x at every station
        ## self.xdot[i] = self.Dx*self.x[i] + self.Ex
        ## self.xddot[i] = self.Ax*self.x[i] + self.Bx
        ## self.thetadot[i] = self.Dth*self.theta[i] + self.Eth
        ## self.thetaddot[i] = self.Ath*self.theta[i] + self.Bth

        self.xdot[i] = self.Dx_vect[i]*self.x[i] + self.Ex_vect[i]
        self.xddot[i] = self.A*self.x[i] + self.Bx_vect[i]
        self.thetadot[i] = self.Dth_vect[i]*self.theta[i] + self.Eth_vect[i]
        self.thetaddot[i] = self.A*self.theta[i] + self.Bth_vect[i]


    def calculate_ABDE(self, i, dt, int_case=2, debug=0, print_times=0):
        #print('int_case = %i' % int_case)
        #t0 = time.time()
        #if print_times:
        #    print('int_case = %i' % int_case)
        if int_case == 1:
            #Finite difference 
            A = 4.0/(dt**2)
            Bx = -A*(self.x[i-1] + dt*self.xdot[i-1] + \
                     (dt**2/4.0)*self.xddot[i-1])
            Bth = -A*(self.theta[i-1] + dt*self.thetadot[i-1] + \
                      (dt**2/4.0)*self.thetaddot[i-1])
            D = 2.0/dt
            Ex = -D*(self.x[i-1] + 0.5*dt*self.xdot[i-1])
            Eth = -D*(self.theta[i-1] + 0.5*dt*self.thetadot[i-1])
        elif int_case == 2:
            #Newmark Beta
            # B = -1.0/(beta*dt**2)*(self.x[i-1] + dt*self.xdot[i-1] + (0.5-beta)*dt**2*self.xddot[i-1])
            # E = self.xdot[i-1] + dt*((1.0-gamma)*self.xddot[i-1]+gamma*B)      

            ta = time.time()
            A = 1.0/(beta*dt**2)
            tb = time.time()
            Bx = (-1.0/(beta*dt**2))*(self.x[i-1] + \
                                      dt*self.xdot[i-1] + \
                                      (0.5-beta)*dt**2*self.xddot[i-1])
            tc = time.time()
            Bth = (-1.0/(beta*dt**2))*(self.theta[i-1] + \
                                       dt*self.thetadot[i-1] + \
                                       (0.5-beta)*dt**2*self.thetaddot[i-1])
            td = time.time()
            D = gamma/(beta*dt)
            te = time.time()
            Ex = self.xdot[i-1] + \
                 dt*((1.0-gamma)*self.xddot[i-1]+gamma*Bx)
            tf = time.time()
            Eth = self.thetadot[i-1] + \
                  dt*((1.0-gamma)*self.thetaddot[i-1]+gamma*Bth)
            tg = time.time()

            if print_times:
                ta_vect = array([ta, tb, tc, td, te, tf, tg])
                dta_vect = ta_vect[1:] - ta_vect[0:-1]
                print('dta_vect = ' + str(dta_vect))

        elif int_case == 3:
            #Fox-Euler
            #print('Fox-Euler is badddddd!!!!!!!!!!!!!!!')

            #ta = time.time()
            A = 2.0/dt**2
            #tb = time.time()
            Bx = -A*(self.x[i-1]+dt*self.xdot[i-1])
            #tc = time.time()
            Bth = -A*(self.theta[i-1]+dt*self.thetadot[i-1])
            #td = time.time()
            D = 2.0/dt
            #te = time.time()
            Ex = -(D*self.x[i-1]+self.xdot[i-1])
            #tf = time.time()
            Eth = -(D*self.theta[i-1]+self.thetadot[i-1])
            #tg = time.time()

            ## if print_times:
            ##     ta_vect = array([ta, tb, tc, td, te, tf, tg])
            ##     dta_vect = ta_vect[1:] - ta_vect[0:-1]
            ##     print('dta_vect = ' + str(dta_vect))

        elif int_case == 4:
            #Wilson Theta
            A = 6.0/((theta_W*dt)**2)
            Bx = -A*(self.x[i-1]+theta_W*dt*self.xdot[i-1] + \
                     self.xddot[i-1]*(theta_W*dt)**2/3.0)
            D = 3.0/(theta_W*dt)
            Ex = -D*(self.x[i-1]+self.xdot[i-1]*(2.0*theta_W*dt)/3.0 + \
                 self.xddot[i-1]*(theta_W*dt)**2/3.0)
            Bth = -A*(self.theta[i-1]+theta_W*dt*self.thetadot[i-1] + \
                      self.thetaddot[i-1]*(theta_W*dt)**2/3.0)
            Eth = -D*(self.theta[i-1] + \
                      self.thetadot[i-1]*(2.0*theta_W*dt)/3.0 + \
                      self.thetaddot[i-1]*(theta_W*dt)**2/3.0)

        elif int_case == 5:
            #Houbolt
            if i == 1:
                A = 6.0/(dt**2)
                Bx = -2.0/(dt**2)*(3.0*self.x[i-1] + 3.0*dt*self.xdot[i-1] + \
                                   dt**2*self.xddot[i-1])
                D = 3.0/dt
                Ex = -1.0/(2*dt)*(6*self.x[i-1] + 4.0*dt*self.xdot[i-1] + \
                                  dt**2*self.xddot[i-1])
                Bth = -2.0/(dt**2)*(3.0*self.theta[i-1] + 3.0*dt*self.thetadot[i-1] + \
                                    dt**2*self.thetaddot[i-1])
                Eth = -1.0/(2*dt)*(6*self.theta[i-1] + 4.0*dt*self.thetadot[i-1] + \
                                   dt**2*self.thetaddot[i-1])
            elif i == 2:
                A = 2.0/(dt**2)
                Bx = -1.0/(dt**2)*(4.0*self.x[i-1] - 2.0*self.x[i-2] + \
                                   2.0*dt**2*self.xddot[i-2])
                D = 11.0/6.0*dt
                Ex = -1.0/(6*dt)*(16*self.x[i-1] - 5.0*self.x[i-2] + \
                                  dt**2*self.xddot[i-2])
                Bth = -1.0/(dt**2)*(4.0*self.theta[i-1] - 2.0*self.theta[i-2] + \
                                   2.0*dt**2*self.thetaddot[i-2])
                Eth = -1.0/(6*dt)*(16*self.theta[i-1] - 5.0*self.theta[i-2] + \
                                   dt**2*self.thetaddot[i-2])
            else: 
                A = 2.0/(dt**2)
                Bx = -1.0/(dt**2)*(5.0*self.x[i-1] - 4.0*self.x[i-2] + \
                                   self.x[i-3])
                D = 11.0/6.0*dt
                Ex = -1.0/(6*dt)*(18*self.x[i-1] - 9.0*self.x[i-2] + \
                                  2*self.x[i-3])
                Bth = -1.0/(dt**2)*(5.0*self.theta[i-1] - 4.0*self.theta[i-2] + \
                                    self.theta[i-3])
                Eth = -1.0/(6*dt)*(18*self.theta[i-1] - 9.0*self.theta[i-2] + \
                                   2*self.theta[i-3])


        #t1 = time.time()
            
        self.A = A
        self.Bx = Bx
        self.Bth = Bth
        self.Dx = D
        self.Dth = D
        self.Ex = Ex
        self.Eth = Eth
        #save for debugging purposes

        #t2 = time.time()

        #self.Ax_vect[i] = A
        #self.Ath_vect[i] = A
        self.Bx_vect[i] = Bx
        self.Bth_vect[i] = Bth
        self.Dx_vect[i] = D
        self.Dth_vect[i] = D
        self.Ex_vect[i] = Ex
        self.Eth_vect[i] = Eth

        #t3 = time.time()

        ## if print_times:
        ##     t_vect = array([t0, t1, t2, t3])
        ##     dt_vect = t_vect[1:] - t_vect[0:-1]
        ##     print('dt_vect = ' + str(dt_vect))

        if debug:
            x_mat = array([[A, Bx],[D, Ex]])
            self.int_det_x[i] = numpy.linalg.det(x_mat)
            self.int_cond_x[i] = numpy.linalg.cond(x_mat)

            th_mat = array([[A, Bth],[D, Eth]])
            self.int_det_th[i] = numpy.linalg.det(th_mat)
            self.int_cond_th[i] = numpy.linalg.cond(th_mat)

        
        

    def calculate_state_vector(self, prev_z, i, xdof=0, thetadof=1):
        self.z = dot(self.U, prev_z)
        self.x[i] = self.z[xdof]
        self.theta[i] = self.z[thetadof]
        self.z_mat[i,:] = squeeze(self.z[0:4])
        return self.z
        


class DT_TMM_System_clamped_free_four_states(DT_TMM_System):
    """This class models a DT-TMM System with clamped/free boundary
    conditions and two states (displacement and force - i.e. 1 DOF
    mass/spring/damper carts vibrating sideways)."""
    def __init__(self, element_list, N_states=4, \
                 **kwargs):
        DT_TMM_System.__init__(self, element_list, N_states=N_states, **kwargs)

    
    def solve_boundary_conditions(self, i=None):
        ###########################################
        #
        # [   x_tip   ]          [    0    ]
        # [ theta_tip ]          [    0    ]
        # [     0     ]  = U_sys [  M_base ]
        # [     0     ]          [  V_base ]
        # [     1     ]          [    1    ]
        #
        ###########################################
        submat = self.U_sys[2:4, 2:4]
        last_col = self.U_sys[2:4,4]
        cond = numpy.linalg.cond(submat)
        #print('i=%i, cond=%s' % (i, cond))
        assert cond is not nan, 'singular matrix'
        submat_inv = numpy.linalg.inv(submat)
        MV_base = squeeze(dot(submat_inv, -last_col))
        z0 = array([[0.0],[0.0],[MV_base[0]], [MV_base[1]],[1.0]])
        self.z0 = z0


    def _set_initial_conditions(self, initial_conditions):
        assert len(initial_conditions) == len(self.element_list), "The length of the outer list of " + \
               "initial_conditions must have the same length as self.element_list:\n" + \
               "len(self.element_list) = %i" % len(self.element_list) + \
               "len(initial_conditions) = %i" % len(initial_conditions)
        for element, ic in zip(self.element_list, initial_conditions):
            assert len(ic) in [2,4], "each element in initial_conditions must be a list of length 2 or 4:\n" + \
                   "ic = %s" % ic
            element.x[0] = ic[0]
            element.theta[0] = ic[1]
            if len(ic) > 2:
                element.xdot[0] = ic[2]
                element.thetadot[0] = ic[3]


    def _assign_t_vect(self):
        for element in self.element_list:
            element.t = self.t


    def _set_dt(self, dt):
        self.dt = dt
        for element in self.element_list:
            element.dt = dt
            
        
    def Run_Simulation(self, N, dt, initial_conditions=None, \
                       int_case=2):
        """If initial_conditions is not None, then it must be a list
        of lists.  The outer list must have the same length as
        self.element_list.  Each inner list must be either 2 or 4
        elements containing either [x0, th0] or [x0, th0, x0_dot,
        th0_dot] for each element."""
        self._initialize_vectors(N)
        self.t = arange(0, N*dt, dt)
        self._assign_t_vect()
        self._set_dt(dt)
        if initial_conditions is not None:
            self._set_initial_conditions(initial_conditions)
        for i in range(1,N):    # Time loop
            self.calculate_ABDE(i, dt, int_case=int_case)
            self.calculate_transfer_matrices(i)
            self.calculate_system_transfer_matrix(i)
            self.solve_boundary_conditions(i)#must set self.z0
            self.calculate_state_vectors(i)
            self.calculate_velocity_and_accel(i)


    def plot_element_attr(self, attr, fig=None, fignum=1, clear=True):
        if fig is None:
            import pylab
            fig = pylab.figure(fignum)
        if clear:
            fig.clf()
            ax = fig.add_subplot(111)
        else:
            ax = fig.axes[0]
        for i, element in enumerate(self.element_list):
            myvect = getattr(element, attr)
            mylabel = attr + '_' + str(i)
            ax.plot(self.t, myvect, label=mylabel)
        

class DT_TMM_Sensor(object):
    def __init__(self, name, signal, sensor_type, elem1, elem2=None, gain=1.0):
        self.name = name
        self.signal = signal
        self.sensor_type = sensor_type
        self.elem1 = elem1
        self.elem2 = elem2
        self.gain = gain
        

    def initialize_vector(self, N):
        self.signal_vect = zeros(N)
        

    def read(self, i):
        vect1 = getattr(self.elem1, self.signal)
        val1 = vect1[i]
        if self.sensor_type.lower() == 'diff':
            vect2 = getattr(self.elem2, signal)
            val2 = vect2[i]
            out = val1 - val2
        else:
            out = val1
        self.signal_vect[i] = out*self.gain
        return self.signal_vect[i]
    

class DT_TMM_Actuator(object):
    def __init__(self, element, name, attr):
        self.element = element
        self.name = name
        self.attr = attr


    def initialize_vector(self, N):
        in_vect = getattr(self.element, self.attr)
        if in_vect is None:
            in_vect = zeros(N)
            setattr(self.element, self.attr, in_vect)
            

    def set_input(self, i, val):
        in_vect = getattr(self.element, self.attr)
        in_vect[i] = val
        

class DT_TMM_System_from_XML_four_states(DT_TMM_System_clamped_free_four_states):
    def __init__(self, element_list, sensors=[], actuators=[], \
                 **kwargs):
        DT_TMM_System_clamped_free_four_states.__init__(self, element_list, **kwargs)
        self.sensors = sensors
        self.actuators = actuators


    def init_sim(self, N, initial_conditions=None):
        self._initialize_vectors(N)
        self.t = arange(0, N*self.dt, self.dt)
        self._assign_t_vect()
        for sensor in self.sensors:
            sensor.initialize_vector(N)

        for act in self.actuators:
            act.initialize_vector(N)
            
        if initial_conditions is not None:
            self._set_initial_conditions(initial_conditions)
        

    def run_sim_one_step(self, i, inputs=[], \
                         int_case=2):
        """If initial_conditions is not None, then it must be a list
        of lists.  The outer list must have the same length as
        self.element_list.  Each inner list must be either 2 or 4
        elements containing either [x0, th0] or [x0, th0, x0_dot,
        th0_dot] for each element."""
        for act, val in zip(self.actuators, inputs):
            act.set_input(i, val)
            
        self.calculate_ABDE(i, self.dt, int_case=int_case)
        self.calculate_transfer_matrices(i)
        self.calculate_system_transfer_matrix(i)
        self.solve_boundary_conditions(i)#must set self.z0
        self.calculate_state_vectors(i)
        self.calculate_velocity_and_accel(i)

        for sensor in self.sensors:
            sensor.read(i)
        

    
###########################################
#
# Next steps:
#
# - create TSD and rigid mass 3 state elements
#
# - create system with one TSD and one 4 state rigid mass and find its
#   step response
#
#    > compare this step response to the known symbolic solution (or
#      odeint results)
#
# - create cantilevered beam model based on a chain of TSD's and rigid
#   masses
#
#    > compare the Bode plot of the beam from tradional TMM with that
#      based on FFT of DT-TMM impulse response
#
# - develop DT-TMM actuator model
#
#    > predict the open-loop impulse response
#
###########################################


class DT_TMM_TSD_4_states(DT_TMM_Element_4_states):
    """This class models a torsional spring/damper with 4 states."""
    def __init__(self, k=0.0, b=0.0, prev_element=None, **kwargs):
        """A torsional spring/damper element.

        k is the stiffness coefficient (N/m).  b is the damping
        coefficient (N*s/m).

        The spring needs the displacement and velocity of the previous
        element to calculate the relative displacement and velocity of
        its spring and damper.  If prev_element is None, it is assumed
        the spring is connected to a wall and the D and E parameters
        of the previous element are zero."""
        DT_TMM_Element_4_states.__init__(self, **kwargs)
        self.k = k
        self.b = b
        self.prev_element = prev_element


    def calculate_transfer_matrix(self, i=None):
        den = (self.k+self.b*self.Dth)
        if self.prev_element is None:
            D_prev = 0.0
            E_prev = 0.0
        else:
            ## D_prev = self.prev_element.Dth
            ## E_prev = self.prev_element.Eth
            D_prev = self.prev_element.Dth_vect[i]
            E_prev = self.prev_element.Eth_vect[i]

        k = self.k
        b = self.b
        self.U = array([[1.0, 0.0, 0.0, 0.0, 0.0], \
                        [0.0, 1.0,  1.0/den,  0.0,\
                         -b*(self.Eth_vect[i] - E_prev)/den], \
                        [0.0, 0.0, 1.0, 0.0, 0.0], \
                        [0.0, 0.0, 0.0, 1.0, 0.0], \
                        [0.0, 0.0, 0.0, 0.0, 1.0]])
        return self.U



class DT_TMM_rigid_mass_4_states(DT_TMM_Element_4_states):
    def __init__(self, m, L, r, I, f=None, prev_element=None, **kwargs):
        """m is the mass of the rigid element and L is its length.  r
        is the distance from the previous element to the center of
        gravity of rigid element.  I is the second moment of inertia
        about the center of gravity.

        f is the applied external force acting on the mass.  The f
        vector can be passed in here.  If it is, f must be left as
        None in calculate_transfer_matrix.  Then f[i] will be applied
        to the mass at each step in the simulation."""
        DT_TMM_Element_4_states.__init__(self, **kwargs)
        self.m = m
        self.L = L
        self.r = r
        self.I = I
        self.f = f
        self.prev_element = prev_element


    def calculate_transfer_matrix(self, i=None, f=None):
        """If f is passed in here, it overrides self.f.  If f is a
        vector, f[i] is applied to the mass (then i must not be None)."""
        if f is None:
            f = self.f
        if f is None:
            f_i = 0.0
        elif (not isscalar(f)):
            f_i = f[i]
            
        if self.prev_element is None:
            A = 0.0#<--- this seems risky, it isn't really true, but it should only multiply things that are 0.0
            Bx_prev = 0.0
            Ath = 0.0
            Bth = 0.0
        else:
            ## Ax_prev = self.prev_element.Ax
            ## Bx_prev = self.prev_element.Bx
            A = self.prev_element.A
            Bx_prev = self.prev_element.Bx_vect[i]
            Bth = self.prev_element.Bth_vect[i]

        L = self.L
        m = self.m
        I = self.I
        r = self.r
        #Ath = self.Ath
        #Bth = self.Bth
        ## Ath = self.prev_element.Ath
        ## Bth = self.prev_element.Bth

        
        self.U = array([[1.0, L, 0.0, 0.0, 0.0], \
                        [0.0, 1.0, 0.0, 0.0, 0.0], \
                        [-(L-r)*A*m, -((L*r-r**2)*m-I)*A, 1.0, -L, \
                         -(L*m-m*r)*Bx_prev-(L*m*r-m*r**2-I)*Bth], \
                        [A*m, A*m*r, 0.0, 1.0, Bth*m*r + Bx_prev*m-f_i], \
                        [0.0, 0.0, 0.0, 0.0, 1.0]])
        return self.U
    


class DT_TMM_rigid_mass_4_states_Taylor(DT_TMM_rigid_mass_4_states):
    def calculate_transfer_matrix(self, i=None, f=None, debug=1):
        if f is None:
            f = self.f
        if f is None:
            f_i = 0.0
        elif (not isscalar(f)):
            f_i = f[i]

        if self.prev_element is None:
            A = 0.0
            Bx_prev = 0.0
            Bth = 0.0
        else:
            ## Ax_prev = self.prev_element.Ax
            ## Bx_prev = self.prev_element.Bx
            A = self.prev_element.A
            Bxp = self.prev_element.Bx_vect[i]
            Bth = self.prev_element.Bth_vect[i]

        L = self.L
        m = self.m
        I = self.I
        r = self.r
        F = f_i
        #Ath = self.Ath
        #Bth = self.Bth
        ## Ath = self.prev_element.Ath
        ## Bth = self.prev_element.Bth

        if i > 0:
            th0m1 = self.theta[i-1]
            thdotm1 = self.thetadot[i-1]
            thddotm1 = self.thetaddot[i-1]
        else:
            th0m1 = 0.0#<-- should probably do an initial condition here
            thdotm1 = 0.0
            thddotm1 = 0.0

        dt = self.dt
        Sbar = sin(th0m1) - th0m1*cos(th0m1) - 0.5*sin(th0m1)*(thdotm1*dt)**2
        cbar = cos(th0m1)*(1.0-0.5*(thdotm1*dt)**2) - sin(th0m1)*(thdotm1*dt+0.5*thddotm1*dt**2)
        u35 = -Bth*L*cbar**2*m*r + Bth*cbar**2*m*r**2 \
              -Bxp*L*cbar*m + Bxp*cbar*m*r + Bth*I

        if debug > 0:
            def myprint(lhs, var):
                fmt = '%s = %0.6g'
                print(fmt % (lhs, var))
            myprint('th0m1', th0m1)
            myprint('Sbar', Sbar)
            myprint('cbar', cbar)
        
        self.U = array([[1, L*cos(th0m1), 0, 0, L*Sbar], \
                        [0, 1, 0, 0, 0,], \
                        [-A*L*cbar*m + A*cbar*m*r, -A*L*cbar**2*m*r + A*cbar**2*m*r**2 + A*I, \
                         1, -L*cbar, u35], \
                        [A*m, A*cbar*m*r, 0, 1, Bth*cbar*m*r + Bxp*m - F], \
                        [0, 0, 0, 0, 1]])


class DT_TMM_DC_motor_4_states(DT_TMM_Element_4_states):
    def c2d(self):
        p = self.p
        g = self.g
        dt = self.dt
        
        # From Dorsey, pages 472-473, adapted from Ex. 14.17.1:
        b1 = p*dt-1.0+exp(-p*dt)
        b2 = 1.0-exp(-p*dt)-p*dt*exp(-p*dt)
        b_Dorsey = g*array([b1,b2])
        a1 = 1.0
        a2 = -(1.0+exp(-p*dt))
        a3 = exp(-p*dt)
        a_Dorsey = p*array([a1,a2,a3])

        a_vect = a_Dorsey/a_Dorsey[0]
        b_vect = b_Dorsey/a_Dorsey[0]

        self.num_z = numpy.append([0.0],b_vect)
        self.den_z = a_vect
        self.Nden = len(self.den_z)
        
        
    def __init__(self, g, p, dt, v=None, **kwargs):
        """Model a DC motor as an actuator for a flexible robot.  The
        motor transfer function is

          theta        g*p
        --------- = --------
         voltage     s(s+p)

        The motor transfer function is coverted to digital using ZOH
        and the formula borrowed from Dorsey (pages 472-473).
         
        v is the motor voltage.  The v vector can be passed in here.
        If it is, v must be left as None in calculate_transfer_matrix.
        Then v[i] will be applied to the motor at each step in the
        simulation."""
        DT_TMM_Element_4_states.__init__(self, **kwargs)
        self.g = g
        self.p = p
        self.dt = dt

        self.v = v

        self.c2d()

    def _initialize_vectors(self, N):
        DT_TMM_Element_4_states._initialize_vectors(self, N)
        self.theta_act = zeros(N)
        


    def calculate_theta_act(self, i):
        out = 0.0
        for n, bn in enumerate(self.num_z):
            if i >= n:
                newpart = self.v[i-n]*bn
                ## print('i, n, bn, newpart = %i, %i, %0.4g, %0.4g' % \
                ##       (i, n, bn, newpart))
                out += newpart

        for n in range(1, self.Nden):
            if i >= n:
                an = self.den_z[n]
                newpart = self.theta_act[i-n]*an
                ## print('i, n, an, newpart = %i, %i, %0.4g, %0.4g' % \
                ##       (i, n, an, newpart))
                out -= newpart
        out = out/self.den_z[0]
        self.theta_act[i] = out
        return out
        

    def calculate_transfer_matrix(self, i, v=None):
        """If v is passed in here, it overrides self.v.  If v is a
        vector, v[i] is applied to the mass (then i must not be None)."""
        if v is None:
            v = self.v
        if v is None:
            v_i = 0.0
        elif (not isscalar(v)):
            v_i = v[i]

        self.v[i] = v_i

        theta_act = self.calculate_theta_act(i)
        self.U = eye(5)
        self.U[1,4] = theta_act
        return self.U



class forcing_element_4_states(DT_TMM_Element_4_states):
    """This is a base class for external forcing inputs to a DT-TMM
    model.  The plan is to try and just treat them like pure TMM
    forcing elements: generate the transfer matrix directly without
    worrying about ABDE, velocity or acceleration.  Hopefully this
    doesn't mess up other elements or the system.  So, several methods
    that are called for each element in a system are overridden to
    simply pass.  This is a base class and all derived classes must
    override the method """
    def calculate_velocity_and_accel(self, i):
        pass
    
    def calculate_ABDE(self, i, dt, int_case=2, debug=0, print_times=0):
        pass


    def calculate_transfer_matrix(self, i, f=None):
        raise NotImplementedError, "classes derived from forcing_element_4_states must override calculate_transfer_matrix"
    

class transverse_forcing_element_4_states(forcing_element_4_states):
    def calculate_transfer_matrix(self, i, f=None):
        """If f is passed in here, it overrides self.f.  If f is a
        vector, f[i] is applied to the mass (then i must not be None)."""
        if f is None:
            f = self.f
        if f is None:
            f_i = 0.0
        elif (not isscalar(f)):
            f_i = f[i]

        self.f[i] = f_i

        self.U = eye(5)
        self.U[3,4] = f_i
        return self.U


    def _initialize_vectors(self, N):
        DT_TMM_Element_4_states._initialize_vectors(self, N)
        if not hasattr(self, 'f'):
            self.f = zeros(N)



# Rui 7x7 analysis
# =============================

class DT_TMM_Element_6_states(DT_TMM_Element):
    def _initialize_vectors(self, N):
        self.x = zeros(N)
        self.xdot = zeros(N)
        self.xddot = zeros(N)        
        self.y = zeros(N)
        self.ydot = zeros(N)
        self.yddot = zeros(N)
        self.theta = zeros(N)
        self.thetadot = zeros(N)
        self.thetaddot = zeros(N)
        self.z_mat = zeros((N,6))
        #self.Ax_vect = zeros(N)
        self.Bx_vect = zeros(N)
        #self.Dx_vect = zeros(N)
        self.Ex_vect = zeros(N)
        self.By_vect = zeros(N)
        self.Ey_vect = zeros(N)
        #self.Ath_vect = zeros(N)
        self.Bth_vect = zeros(N)
        #self.Dth_vect = zeros(N)
        self.Eth_vect = zeros(N)
        self.int_det_x = zeros(N)
        self.int_det_y = zeros(N)
        self.int_det_th = zeros(N)
        self.int_cond_x = zeros(N)
        self.int_cond_y = zeros(N)        
        self.int_cond_th = zeros(N)


    def calculate_velocity_and_accel(self, i):
        #to find x at every station
        ## self.xdot[i] = self.Dx*self.x[i] + self.Ex
        ## self.xddot[i] = self.Ax*self.x[i] + self.Bx
        ## self.thetadot[i] = self.Dth*self.theta[i] + self.Eth
        ## self.thetaddot[i] = self.Ath*self.theta[i] + self.Bth

        self.xdot[i] = self.D*self.x[i] + self.Ex_vect[i]
        self.xddot[i] = self.A*self.x[i] + self.Bx_vect[i]
        self.ydot[i] = self.D*self.y[i] + self.Ey_vect[i]
        self.yddot[i] = self.A*self.y[i] + self.By_vect[i]
        self.thetadot[i] = self.D*self.theta[i] + self.Eth_vect[i]
        self.thetaddot[i] = self.A*self.theta[i] + self.Bth_vect[i]


    def calculate_ABDE(self, i, dt, int_case=2, debug=0, print_times=0):
        #print('int_case = %i' % int_case)
        #t0 = time.time()
        #if print_times:
        #    print('int_case = %i' % int_case)
        if int_case == 1:
            #Finite difference 
            A = 4.0/(dt**2)
            Bx = -A*(self.x[i-1] + dt*self.xdot[i-1] + \
                     (dt**2/4.0)*self.xddot[i-1])
            By = -A*(self.y[i-1] + dt*self.ydot[i-1] + \
                     (dt**2/4.0)*self.yddot[i-1])
            Bth = -A*(self.theta[i-1] + dt*self.thetadot[i-1] + \
                      (dt**2/4.0)*self.thetaddot[i-1])
            D = 2.0/dt
            Ex = -D*(self.x[i-1] + 0.5*dt*self.xdot[i-1])
            Ey = -D*(self.y[i-1] + 0.5*dt*self.ydot[i-1])            
            Eth = -D*(self.theta[i-1] + 0.5*dt*self.thetadot[i-1])
        elif int_case == 2:
            #Newmark Beta
            # B = -1.0/(beta*dt**2)*(self.x[i-1] + dt*self.xdot[i-1] + (0.5-beta)*dt**2*self.xddot[i-1])
            # E = self.xdot[i-1] + dt*((1.0-gamma)*self.xddot[i-1]+gamma*B)      

            ta = time.time()
            A = 1.0/(beta*dt**2)
            tb = time.time()
            Bx = -1.0/(beta*dt**2)*(self.x[i-1] + \
                                    dt*self.xdot[i-1] + \
                                    (0.5-beta)*dt**2*self.xddot[i-1])
            By = -1.0/(beta*dt**2)*(self.y[i-1] + \
                                    dt*self.ydot[i-1] + \
                                    (0.5-beta)*dt**2*self.yddot[i-1])
            tc = time.time()
            Bth = -1.0/(beta*dt**2)*(self.theta[i-1] + \
                                     dt*self.thetadot[i-1] + \
                                     (0.5-beta)*dt**2*self.thetaddot[i-1])
            td = time.time()
            D = gamma/(beta*dt)
            te = time.time()
            Ex = self.xdot[i-1] + \
                 dt*((1.0-gamma)*self.xddot[i-1]+gamma*Bx)
            Ey = self.ydot[i-1] + \
                 dt*((1.0-gamma)*self.yddot[i-1]+gamma*By)
            tf = time.time()
            Eth = self.thetadot[i-1] + \
                  dt*((1.0-gamma)*self.thetaddot[i-1]+gamma*Bth)
            tg = time.time()

            if print_times:
                ta_vect = array([ta, tb, tc, td, te, tf, tg])
                dta_vect = ta_vect[1:] - ta_vect[0:-1]
                print('dta_vect = ' + str(dta_vect))

        elif int_case == 3:
            #Fox-Euler
            #print('Fox-Euler is badddddd!!!!!!!!!!!!!!!')

            #ta = time.time()
            A = 2.0/dt**2
            #tb = time.time()
            Bx = -A*(self.x[i-1]+dt*self.xdot[i-1])
            By = -A*(self.y[i-1]+dt*self.ydot[i-1])            
            #tc = time.time()
            Bth = -A*(self.theta[i-1]+dt*self.thetadot[i-1])
            #td = time.time()
            D = 2.0/dt
            #te = time.time()
            Ex = -(D*self.x[i-1]+self.xdot[i-1])
            Ey = -(D*self.y[i-1]+self.ydot[i-1])            
            #tf = time.time()
            Eth = -(D*self.theta[i-1]+self.thetadot[i-1])
            #tg = time.time()

            ## if print_times:
            ##     ta_vect = array([ta, tb, tc, td, te, tf, tg])
            ##     dta_vect = ta_vect[1:] - ta_vect[0:-1]
            ##     print('dta_vect = ' + str(dta_vect))

        elif int_case == 4:
            #Wilson Theta
            A = 6.0/((theta_W*dt)**2)
            Bx = -A*(self.x[i-1]+theta_W*dt*self.xdot[i-1] + \
                     self.xddot[i-1]*(theta_W*dt)**2/3.0)
            By = -A*(self.y[i-1]+theta_W*dt*self.ydot[i-1] + \
                     self.yddot[i-1]*(theta_W*dt)**2/3.0)
            D = 3.0/(theta_W*dt)
            Ex = -D*(self.x[i-1]+self.xdot[i-1]*(2.0*theta_W*dt)/3.0 + \
                 self.xddot[i-1]*(theta_W*dt)**2/3.0)
            Ey = -D*(self.y[i-1]+self.ydot[i-1]*(2.0*theta_W*dt)/3.0 + \
                 self.yddot[i-1]*(theta_W*dt)**2/3.0)
            Bth = -A*(self.theta[i-1]+theta_W*dt*self.thetadot[i-1] + \
                      self.thetaddot[i-1]*(theta_W*dt)**2/3.0)
            Eth = -D*(self.theta[i-1] + \
                      self.thetadot[i-1]*(2.0*theta_W*dt)/3.0 + \
                      self.thetaddot[i-1]*(theta_W*dt)**2/3.0)

        elif int_case == 5:
            #Houbolt
            if i == 1:
                A = 6.0/(dt**2)
                Bx = -2.0/(dt**2)*(3.0*self.x[i-1] + 3.0*dt*self.xdot[i-1] + \
                                   dt**2*self.xddot[i-1])
                By = -2.0/(dt**2)*(3.0*self.y[i-1] + 3.0*dt*self.ydot[i-1] + \
                                   dt**2*self.yddot[i-1])
                D = 3.0/dt
                Ex = -1.0/(2*dt)*(6*self.x[i-1] + 4.0*dt*self.xdot[i-1] + \
                                  dt**2*self.xddot[i-1])
                Ey = -1.0/(2*dt)*(6*self.y[i-1] + 4.0*dt*self.ydot[i-1] + \
                                  dt**2*self.yddot[i-1])
                Bth = -2.0/(dt**2)*(3.0*self.theta[i-1] + 3.0*dt*self.thetadot[i-1] + \
                                    dt**2*self.thetaddot[i-1])
                Eth = -1.0/(2*dt)*(6*self.theta[i-1] + 4.0*dt*self.thetadot[i-1] + \
                                   dt**2*self.thetaddot[i-1])
            elif i == 2:
                A = 2.0/(dt**2)
                Bx = -1.0/(dt**2)*(4.0*self.x[i-1] - 2.0*self.x[i-2] + \
                                   2.0*dt**2*self.xddot[i-2])
                By = -1.0/(dt**2)*(4.0*self.y[i-1] - 2.0*self.y[i-2] + \
                                   2.0*dt**2*self.yddot[i-2])
                D = 11.0/6.0*dt
                Ex = -1.0/(6*dt)*(16*self.x[i-1] - 5.0*self.x[i-2] + \
                                  dt**2*self.xddot[i-2])
                Ey = -1.0/(6*dt)*(16*self.y[i-1] - 5.0*self.y[i-2] + \
                                  dt**2*self.yddot[i-2])
                Bth = -1.0/(dt**2)*(4.0*self.theta[i-1] - 2.0*self.theta[i-2] + \
                                   2.0*dt**2*self.thetaddot[i-2])
                Eth = -1.0/(6*dt)*(16*self.theta[i-1] - 5.0*self.theta[i-2] + \
                                   dt**2*self.thetaddot[i-2])
            else: 
                A = 2.0/(dt**2)
                Bx = -1.0/(dt**2)*(5.0*self.x[i-1] - 4.0*self.x[i-2] + \
                                   self.x[i-3])
                By = -1.0/(dt**2)*(5.0*self.y[i-1] - 4.0*self.y[i-2] + \
                                   self.y[i-3])
                D = 11.0/6.0*dt
                Ex = -1.0/(6*dt)*(18*self.x[i-1] - 9.0*self.x[i-2] + \
                                  2*self.x[i-3])
                Ey = -1.0/(6*dt)*(18*self.y[i-1] - 9.0*self.y[i-2] + \
                                  2*self.y[i-3])
                Bth = -1.0/(dt**2)*(5.0*self.theta[i-1] - 4.0*self.theta[i-2] + \
                                    self.theta[i-3])
                Eth = -1.0/(6*dt)*(18*self.theta[i-1] - 9.0*self.theta[i-2] + \
                                   2*self.theta[i-3])


        #t1 = time.time()

        self.A = A
        self.Bx = Bx
        self.By = By
        self.Bth = Bth
        self.D = D
        self.Ex = Ex
        self.Ey = Ey
        self.Eth = Eth
        #save for debugging purposes

        #t2 = time.time()

        self.Bx_vect[i] = Bx
        self.By_vect[i] = By        
        self.Bth_vect[i] = Bth
        self.Ex_vect[i] = Ex
        self.Ey_vect[i] = Ey        
        self.Eth_vect[i] = Eth

        #t3 = time.time()

        ## if print_times:
        ##     t_vect = array([t0, t1, t2, t3])
        ##     dt_vect = t_vect[1:] - t_vect[0:-1]
        ##     print('dt_vect = ' + str(dt_vect))

        if debug:
            x_mat = array([[A, Bx],[D, Ex]])
            self.int_det_x[i] = numpy.linalg.det(x_mat)
            self.int_cond_x[i] = numpy.linalg.cond(x_mat)

            th_mat = array([[A, Bth],[D, Eth]])
            self.int_det_th[i] = numpy.linalg.det(th_mat)
            self.int_cond_th[i] = numpy.linalg.cond(th_mat)




    def calculate_state_vector(self, prev_z, i, xdof=0, ydof=1, thetadof=2):
        self.z = dot(self.U, prev_z)
        self.x[i] = self.z[xdof]
        self.y[i] = self.z[ydof]        
        self.theta[i] = self.z[thetadof]
        self.z_mat[i,:] = squeeze(self.z[0:6])
        return self.z



class DT_TMM_TSD_6_states(DT_TMM_Element_6_states):
    """This class models a torsional spring/damper with 6 states.  The
    state vector is

    z = [x, y, theta3, mz, qx, qy, 1]^T

    where theta3 is the angular displacement about the z axis."""
    def __init__(self, k=0.0, b=0.0, prev_element=None, **kwargs):
        """A torsional spring/damper element.

        k is the stiffness coefficient (N/m).  b is the damping
        coefficient (N*s/m).

        The spring needs the displacement and velocity of the previous
        element to calculate the relative displacement and velocity of
        its spring and damper.  If prev_element is None, it is assumed
        the spring is connected to a wall and the D and E parameters
        of the previous element are zero."""
        DT_TMM_Element_6_states.__init__(self, **kwargs)
        self.k = k
        self.b = b
        self.prev_element = prev_element


    def calculate_transfer_matrix(self, i=None):
        den = (self.k+self.b*self.D)
        if self.prev_element is None:
            D_prev = 0.0
            E_prev = 0.0
        else:
            ## D_prev = self.prev_element.Dth
            ## E_prev = self.prev_element.Eth
            D_prev = self.prev_element.D
            E_prev = self.prev_element.Eth_vect[i]

        k = self.k
        b = self.b
        self.U = eye(7)

        #z = [x, y, theta3, mz, qx, qy, 1]^T
        u34 = 1.0/den#starting with index 1
        self.U[2,3] = u34#zero index
        u37 = -b*(self.Eth_vect[i] - E_prev)/den
        self.U[2,6] = u37
        return self.U


class DT_TMM_rigid_mass_6_states(DT_TMM_Element_6_states):
    def __init__(self, m, L, r, I, f=None, prev_element=None, **kwargs):
        """This is a one-dimensional rigid body for use in
        verification of Rui's 2005 paper and investigation into
        numeric stability.  So, the y dimensions of the element
        (y_{2,o}, y_{IO}, y_{IC},...) will be 0. 

        m is the mass of the rigid element and L is its length.  r
        is the distance from the previous element to the center of
        gravity of rigid element.  I is the second moment of inertia
        about the center of gravity.

        f is the applied external force acting on the mass.  The f
        vector can be passed in here.  If it is, f must be left as
        None in calculate_transfer_matrix.  Then f[i] will be applied
        to the mass at each step in the simulation."""
        DT_TMM_Element_6_states.__init__(self, **kwargs)
        self.m = m
        self.L = L
        self.r = r
        self.I = I
        self.f = f
        self.x2o = L
        self.y2o = 0.0
        self.x2c = r
        self.y2c = 0.0
        self.prev_element = prev_element


    def calculate_transfer_matrix(self, i=None, f=None):
        """If f is passed in here, it overrides self.f.  If f is a
        vector, f[i] is applied to the mass (then i must not be None).

        Note that I am not sure what to do with f on a mass element
        currently because Rui applied f to the center of mass."""
        if f is None:
            f = self.f
        if f is None:
            f_i = 0.0
        elif (not isscalar(f)):
            f_i = f[i]

        if self.prev_element is None:
            A_prev = 0.0
            Bx_prev = 0.0
            By_prev = 0.0
            Bth = 0.0
        else:
            ## Ax_prev = self.prev_element.Ax
            ## Bx_prev = self.prev_element.Bx
            A_prev = self.prev_element.A
            Bx_prev = self.prev_element.Bx_vect[i]
            By_prev = self.prev_element.By_vect[i]            
            Bth = self.prev_element.Bth_vect[i]


        # hard coding zero forces for now
        #
        # - this will only work if I use forcing elements
        fxc = 0.0
        fyc = 0.0
        mc = 0.0
        
        L = self.L
        m = self.m
        I = self.I
        r = self.r
        self.U = eye(7)

        x2o = self.x2o
        y2o = self.y2o
        x2c = self.x2c
        y2c = self.y2c

        dt = self.dt
        
        s = sin(self.theta[i-1])
        c = cos(self.theta[i-1])

        chi1 = self.A
        
        # based on my understanding of the Rui 2005 paper (eqn 32)
        # and the Taylor series expansion that appears in the 2010 JSV paper
        # (eqn 8), I am assuming that s and c in u_{1,3} really mean
        # the sin and cos of theta[i-1]
        u13 = -x2o*s - y2o*c
        self.U[0,2] = u13#zero index

        # comparing the 2005 and 2010 papers, it seems clear that G1 =
        # capital C bar and G2 = capital S bar, defined in the JSV
        # 2010 eqn (8), so the quanties in G1 and G2 in the 2005 paper
        # are all i-1
        thm1 = self.theta[i-1]
        thdotm1 = self.thetadot[i-1]
        thddotm1 = self.thetaddot[i-1]

        G1 = c + thm1*s - 0.5*c*(thdotm1*dt)**2
        G2 = s + thm1*c - 0.5*s*(thdotm1*dt)**2        
        u17 = x2o*G1 - y2o*G2
        self.U[0,6] = u17

        u23 = x2o*c - y2o*s
        self.U[1,2] = u23

        u27 = x2o*G2 + y2o*G1
        self.U[1,6] = u27

        sbar = sin(thm1)*(1.-0.5*(thdotm1*dt)**2) + \
               cos(thm1)*(thdotm1*dt + 0.5*thddotm1*dt**2)
        cbar = cos(thm1)*(1.-0.5*(thdotm1*dt)**2) + \
               -sin(thm1)*(thdotm1*dt + 0.5*thddotm1*dt**2)
        
        xic = x2c*cbar - y2c*sbar
        yic = x2c*sbar + y2c*cbar
        xio = x2o*cbar - y2o*sbar
        yio = x2o*sbar + y2o*cbar

        u51 = -m*chi1
        u53 = m*chi1*(x2c*s+y2c*c)
        u57 = fxc - m*chi1*(x2c*G1-y2c*G2)-m*self.Bx_vect[i]
        print('u57 = ' + str(u57))
        if i == 5:
            Pdb().set_trace()
        self.U[4,0] = u51
        self.U[4,2] = u53
        self.U[4,6] = u57

        u62 = -m*chi1
        u63 = -m*chi1*(x2c*c-y2c*s)
        u67 = fyc - m*chi1*(x2c*G2 + y2c*G2) - m*self.By_vect[i]#<-- the double G2 in this line is
                                                                #    suspicious to me
        self.U[5,1] = u62
        self.U[5,2] = u63
        self.U[5,6] = u67
        
        u41 = m*chi1*(yio-yic)
        u42 = m*chi1*(xic-xio)
        u43 = u63*xio-u53*yio+I*chi1
        u45 = -yio
        u46 = xio
        u47 = -mc + u67*xio - u57*yio + I*self.Bth_vect[i] + \
              xic*(m*self.By_vect[i] - fyc) + \
              (fxc - m*self.Bx_vect[i])*yic
        self.U[3,:] = array([u41,u42,u43,1.0,u45,u46,u47])
        
        return self.U



class forcing_element_6_states(DT_TMM_Element_6_states):
    """This is a base class for external forcing inputs to a DT-TMM
    model.  The plan is to try and just treat them like pure TMM
    forcing elements: generate the transfer matrix directly without
    worrying about ABDE, velocity or acceleration.  Hopefully this
    doesn't mess up other elements or the system.  So, several methods
    that are called for each element in a system are overridden to
    simply pass.  This is a base class and all derived classes must
    override the method """
    def calculate_velocity_and_accel(self, i):
        pass

    def calculate_ABDE(self, i, dt, int_case=2, debug=0, print_times=0):
        pass


    def calculate_transfer_matrix(self, i, f=None):
        raise NotImplementedError, "classes derived from forcing_element_6_states must override calculate_transfer_matrix"


class transverse_forcing_element_6_states(forcing_element_6_states):
    def calculate_transfer_matrix(self, i, f=None):
        """If f is passed in here, it overrides self.f.  If f is a
        vector, f[i] is applied to the mass (then i must not be None)."""
        if f is None:
            f = self.f
        if f is None:
            f_i = 0.0
        elif (not isscalar(f)):
            f_i = f[i]

        self.f[i] = f_i

        self.U = eye(7)
        self.U[5,6] = f_i
        return self.U


    def _initialize_vectors(self, N):
        DT_TMM_Element_6_states._initialize_vectors(self, N)
        if not hasattr(self, 'f'):
            self.f = zeros(N)


class DT_TMM_System_clamped_free_six_states(DT_TMM_System_clamped_free_four_states):
    """This class models a DT-TMM System with clamped/free boundary
    conditions and six states (planar motion)"""
    def __init__(self, element_list, N_states=6, \
                 **kwargs):
        DT_TMM_System.__init__(self, element_list, N_states=N_states, **kwargs)


    def solve_boundary_conditions(self, i=None):
        ###########################################
        #    z = [x, y, theta3, mz, qx, qy, 1]^T
        #
        # [   x_tip   ]          [    0    ]
        # [   y_tip   ]          [    0    ]
        # [ theta_tip ]          [    0    ]
        # [     0     ]  = U_sys [  M_base ]
        # [     0     ]          [ qx_base ]
        # [     0     ]          [ qy_base ]
        # [     1     ]          [    1    ]
        #
        ###########################################
        submat = self.U_sys[3:6, 3:6]
        last_col = self.U_sys[3:6,6]
        cond = numpy.linalg.cond(submat)
        #print('i=%i, cond=%s' % (i, cond))
        assert cond is not nan, 'singular matrix'
        submat_inv = numpy.linalg.inv(submat)
        MV_base = squeeze(dot(submat_inv, -last_col))
        z0 = array([[0.0], \
                    [0.0], \
                    [0.0], \
                    [MV_base[0]], \
                    [MV_base[1]], \
                    [MV_base[2]], \
                    [1.0]])
        self.z0 = z0




if __name__ == "__main__":
    from pylab import figure, show, clf, plot, xlabel, ylabel, title, legend, savefig
    import controls
    
    T = 10.0    # Run time
    dt = 2.0/1000 # Time step
    t = arange(0,T,dt)
    N = len(t)

    M = 2 # Number of stations

    m = 2.0 #mass
    k = 100.0 #spring stiffness
    b = 1.0 # coefficient of damping

    f = zeros(N)
    f[10:] = 1.0 #unit step force input that turns on at i=10

    case = 1
    
    if case == 1:
        f1 = f
        f2 = None
    elif case == 2:
        f1 = None
        f2 = f

    k1 = DT_TMM_spring_2_states(k, b, prev_element=None)
    m1 = DT_TMM_mass_2_states(m, f=f1)
    k2 = DT_TMM_spring_2_states(k, b, prev_element=m1)
    m2 = DT_TMM_mass_2_states(m, f=f2)
    elem_list = [k1, m1, k2, m2]

    sys = DT_TMM_System_clamped_free_two_states(elem_list, N_states=2)
    #N2 = 15
    N2 = N
    sys.Run_Simulation(N2, dt)

    x1 = sys.element_list[1].x
    x2 = sys.element_list[3].x
    

    figure(1)
    clf()
    plot(t[0:N2], x1, label='DT-TMM $x_1$') # Plotting x1
    plot(t[0:N2], x2, label='DT-TMM $x_2$') # Plotting x2


    run_lsim = 1
    lw = 1.5

    if run_lsim:
        m1 = m
        m2 = m
        c1 = b
        c2 = b
        k1 = k
        k2 = k

        #The denominator is the same for both cases:
        Den = [m1*m2, \
               (c1*m2 + c2*m1 + c2*m2), \
               (c1*c2 + k1*m2 + k2*m1 + k2*m2), \
               (c1*k2 + c2*k1), \
               k1*k2]

        if case == 1:
            # From Sage/Maxima:

            # Den = m1*m2*s^4 + (c1*m2 + c2*m1 + c2*m2)*s^3 + (c1*c2 +
            # k1*m2 + k2*m1 + k2*m2)*s^2 + (c1*k2 + c2*k1)*s + k1*k2
            Num1 = [m2, c2, k2]
            Num2 = [c2, k2]

        elif case == 2:
            Num1 = [c2, k2]
            Num2 = [m1, c1+c2, k1+k2]


        G1 = controls.TF(Num1, Den)
        G2 = controls.TF(Num2, Den)

        y1 = G1.lsim(f,t)
        y2 = G2.lsim(f,t)

        plot(t,y1,'--', linewidth=lw, label='lsim $x_1$')
        plot(t,y2,'k-.', linewidth=lw, label='lsim $x_2$')


    xlabel('Time (sec)')
    ylabel('Displacement')
    legend()
    #ylim([-2,2])
    mytitle = 'Case %i' % case
    title(mytitle)

    filename = 'twoDOF_verification_case%i.png' % case
    #savefig(filename, dpi=125)

    show()




