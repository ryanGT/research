from numpy import *

import copy, time

beta = 1.0/6
gamma = 0.5
theta_W = 1.5


class dttmm_sys_from_sim(object):
    def __init__(self, z_func, num_elements, params_in, vector_props=[]):
        self.z_func = z_func
        self.num_elements = num_elements
        self.params_in = params_in
        self.params = copy.copy(params_in)
        self.vector_props = vector_props
        

    def _initialize_vectors(self, N):
        ne = self.num_elements
        myshape = (N,ne)
        self.Bx_mat = zeros(myshape)
        self.Bth_mat = zeros(myshape)
        self.Ex_mat = zeros(myshape)
        self.Eth_mat = zeros(myshape)
        self.x = zeros(myshape)
        self.xdot = zeros(myshape)
        self.xddot = zeros(myshape)
        self.theta = zeros(myshape)
        self.thetadot = zeros(myshape)
        self.thetaddot = zeros(myshape)
        #self.z_array = zeros((N,ne,5))
        

    def calculate_state_vectors(self, i):
        # - pass params into each function in self.z_func_list and get a z vect back
        # - store the z vetors in a 3D array
        # - pull out x and theta for each station and store in self.x[i,j] and
        #   self.theta[i,j]
        # - once you have self.x and self.theta for i, you are done

        ## for j, func in enumerate(self.z_func_list):
        ##     z_j = squeeze(func(self.params))
        ##     #self.z_array[i,j,:] = z_j #<--- I don't really need the other states.
        ##     self.x[i,j] = z_j[0]
        ##     self.theta[i,j] = z_j[1]
        z_mat = self.z_func(self.params)
        self.x[i,:] = z_mat[0,:]
        self.theta[i,:] = z_mat[1,:]
        

    def calculate_velocity_and_accel(self, i):
        #to find x at every station
        self.xdot[i,:] = self.D*self.x[i,:] + self.Ex_mat[i,:]
        self.xddot[i,:] = self.A*self.x[i,:] + self.Bx_mat[i,:]
        self.thetadot[i,:] = self.D*self.theta[i,:] + self.Eth_mat[i,:]
        self.thetaddot[i,:] = self.A*self.theta[i,:] + self.Bth_mat[i,:]
        

        
    def update_params(self,i):
        for prop in self.vector_props:
            vect = self.params_in[prop]
            self.params[prop] = vect[i]


        #Get Bx props
        vect_list = ['Bx','Bth','Ex','Eth']
        for name in vect_list:
            mat_name = name + '_mat'
            mymatrix = getattr(self, mat_name)
            for j in range(self.num_elements):
                key = name + str(j)
                self.params[key] = mymatrix[i,j]

        
    def Run_Simulation(self, N, dt, initial_conditions=None, \
                       int_case=1):
        """If initial_conditions is not None, then it must be a list
        of lists.  The outer list must have the same length as
        self.element_list.  Each inner list must be either 2 or 4
        elements containing either [x0, th0] or [x0, th0, x0_dot,
        th0_dot] for each element."""
        self._initialize_vectors(N)
        self.t = arange(0, N*dt, dt)
        #self._assign_t_vect()
        #self._set_dt(dt)
        self.dt = dt
        self.calculate_A_and_D(dt, int_case)
        if initial_conditions is not None:
            self._set_initial_conditions(initial_conditions)
        dt1 = 0.0
        dt2 = 0.0
        dt3 = 0.0
        dt4 = 0.0
        for i in range(1,N):    # Time loop
            t0 = time.time()
            self.calculate_B_and_E(i, dt, int_case=int_case)
            t1 = time.time()
            self.update_params(i)
            t2 = time.time()
            self.calculate_state_vectors(i)
            t3 = time.time()
            self.calculate_velocity_and_accel(i)
            t4 = time.time()
            dt1 += t1-t0
            dt2 += t2-t1
            dt3 += t3-t2
            dt4 += t4-t3

        print('dt1 = %0.4g' % dt1)
        print('dt2 = %0.4g' % dt2)
        print('dt3 = %0.4g' % dt3)
        print('dt4 = %0.4g' % dt4)
        


    def calculate_A_and_D(self, dt, int_case):
        if int_case == 1:
            #Finite difference 
            A = 4.0/(dt**2)
            D = 2.0/dt        

        elif int_case == 2:
            #Newmark Beta
            A = 1.0/(beta*dt**2)
            D = gamma/(beta*dt)

        elif int_case == 3:
            #Fox-Euler
            A = 2.0/dt**2
            D = 2.0/dt

        elif int_case == 4:
            #Wilson Theta
            A = 6.0/((theta_W*dt)**2)
            D = 3.0/(theta_W*dt)

        self.A = A
        self.D = D
        self.params['A'] = A
        self.params['D'] = D
        

    def calculate_B_and_E(self, i, dt, int_case=2, debug=0, print_times=0):
        A = self.A
        D = self.D
        
        if int_case == 1:
            #Finite difference 
            Bx = -A*(self.x[i-1,:] + dt*self.xdot[i-1,:] + \
                     (dt**2/4.0)*self.xddot[i-1,:])
            Bth = -A*(self.theta[i-1,:] + dt*self.thetadot[i-1,:] + \
                      (dt**2/4.0)*self.thetaddot[i-1,:])
            Ex = -D*(self.x[i-1,:] + 0.5*dt*self.xdot[i-1,:])
            Eth = -D*(self.theta[i-1,:] + 0.5*dt*self.thetadot[i-1,:])
        elif int_case == 2:
            #Newmark Beta
            Bx = (-1.0/(beta*dt**2))*(self.x[i-1,:] + \
                                      dt*self.xdot[i-1,:] + \
                                      (0.5-beta)*dt**2*self.xddot[i-1,:])
            Bth = (-1.0/(beta*dt**2))*(self.theta[i-1,:] + \
                                       dt*self.thetadot[i-1,:] + \
                                       (0.5-beta)*dt**2*self.thetaddot[i-1,:])
            Ex = self.xdot[i-1,:] + \
                 dt*((1.0-gamma)*self.xddot[i-1,:]+gamma*Bx)
            Eth = self.thetadot[i-1,:] + \
                  dt*((1.0-gamma)*self.thetaddot[i-1,:]+gamma*Bth)

        elif int_case == 3:
            #Fox-Euler

            Bx = -A*(self.x[i-1,:]+dt*self.xdot[i-1,:])
            Bth = -A*(self.theta[i-1,:]+dt*self.thetadot[i-1,:])
            Ex = -(D*self.x[i-1,:]+self.xdot[i-1,:])
            Eth = -(D*self.theta[i-1,:]+self.thetadot[i-1,:])

        elif int_case == 4:
            #Wilson Theta
            Bx = -A*(self.x[i-1,:]+theta_W*dt*self.xdot[i-1,:] + \
                     self.xddot[i-1,:]*(theta_W*dt)**2/3.0)
            Ex = -D*(self.x[i-1,:]+self.xdot[i-1,:]*(2.0*theta_W*dt)/3.0 + \
                 self.xddot[i-1,:]*(theta_W*dt)**2/3.0)
            Bth = -A*(self.theta[i-1,:]+theta_W*dt*self.thetadot[i-1,:] + \
                      self.thetaddot[i-1,:]*(theta_W*dt)**2/3.0)
            Eth = -D*(self.theta[i-1,:] + \
                      self.thetadot[i-1,:]*(2.0*theta_W*dt)/3.0 + \
                      self.thetaddot[i-1,:]*(theta_W*dt)**2/3.0)

        ## elif int_case == 5:
        ##     #Houbolt
        ##     if i == 1:
        ##         A = 6.0/(dt**2)
        ##         Bx = -2.0/(dt**2)*(3.0*self.x[i-1,:] + 3.0*dt*self.xdot[i-1,:] + \
        ##                            dt**2*self.xddot[i-1,:])
        ##         D = 3.0/dt
        ##         Ex = -1.0/(2*dt)*(6*self.x[i-1,:] + 4.0*dt*self.xdot[i-1,:] + \
        ##                           dt**2*self.xddot[i-1,:])
        ##         Bth = -2.0/(dt**2)*(3.0*self.theta[i-1,:] + 3.0*dt*self.thetadot[i-1,:] + \
        ##                             dt**2*self.thetaddot[i-1,:])
        ##         Eth = -1.0/(2*dt)*(6*self.theta[i-1,:] + 4.0*dt*self.thetadot[i-1,:] + \
        ##                            dt**2*self.thetaddot[i-1,:])
        ##     elif i == 2:
        ##         A = 2.0/(dt**2)
        ##         Bx = -1.0/(dt**2)*(4.0*self.x[i-1,:] - 2.0*self.x[i-2] + \
        ##                            2.0*dt**2*self.xddot[i-2])
        ##         D = 11.0/6.0*dt
        ##         Ex = -1.0/(6*dt)*(16*self.x[i-1,:] - 5.0*self.x[i-2] + \
        ##                           dt**2*self.xddot[i-2])
        ##         Bth = -1.0/(dt**2)*(4.0*self.theta[i-1,:] - 2.0*self.theta[i-2] + \
        ##                            2.0*dt**2*self.thetaddot[i-2])
        ##         Eth = -1.0/(6*dt)*(16*self.theta[i-1,:] - 5.0*self.theta[i-2] + \
        ##                            dt**2*self.thetaddot[i-2])
        ##     else: 
        ##         A = 2.0/(dt**2)
        ##         Bx = -1.0/(dt**2)*(5.0*self.x[i-1,:] - 4.0*self.x[i-2] + \
        ##                            self.x[i-3])
        ##         D = 11.0/6.0*dt
        ##         Ex = -1.0/(6*dt)*(18*self.x[i-1,:] - 9.0*self.x[i-2] + \
        ##                           2*self.x[i-3])
        ##         Bth = -1.0/(dt**2)*(5.0*self.theta[i-1,:] - 4.0*self.theta[i-2] + \
        ##                             self.theta[i-3])
        ##         Eth = -1.0/(6*dt)*(18*self.theta[i-1,:] - 9.0*self.theta[i-2] + \
        ##                            2*self.theta[i-3])


        self.Bx = Bx
        self.Bth = Bth
        self.Ex = Ex
        self.Eth = Eth

        self.Bx_mat[i,:] = Bx
        self.Bth_mat[i,:] = Bth
        self.Ex_mat[i,:] = Ex
        self.Eth_mat[i,:] = Eth

        ## if debug:
        ##     x_mat = array([[A, Bx],[D, Ex]])
        ##     self.int_det_x[i] = numpy.linalg.det(x_mat)
        ##     self.int_cond_x[i] = numpy.linalg.cond(x_mat)

        ##     th_mat = array([[A, Bth],[D, Eth]])
        ##     self.int_det_th[i] = numpy.linalg.det(th_mat)
        ##     self.int_cond_th[i] = numpy.linalg.cond(th_mat)

