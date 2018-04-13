import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
import numpy as np
fft = np.fft.fft
import time, os, glob
from numpy import sin, cos, pi, arange, zeros, zeros_like
import bode_utils

import serial_utils, txt_mixin
import control
import control_utils
import digcomp

class beam_serial_test(serial_utils.serial_test):
    def print_all_wait(self):
        n = 0
        test = ''
            
        while not test:
            out = self.read_all()
            test = out.strip()
            n += 1;
            if n > 100:
                break
            if not test:
                time.sleep(0.5)

        # at this point, we have a non-blank response
        # now read until it is blank

        for i in range(20):
            out2 = self.read_all()
            out += out2
            if not out2:
                break
            else:
                time.sleep(0.1)

        print(out)


    def create_step_u(self, amp, N=500, start_ind=50):
        if start_ind is None:
            start_ind=50
        u_step = np.zeros(N,dtype=int)
        u_step[start_ind:] = amp
        return u_step

        
    def start_inv_db(self):
        self.flush()
        self.write_byte(5)
        
    
    def get_inv_db_res(self):
        self.write_byte(6)
        pos_break = self.read_two_bytes()
        neg_break = self.read_two_bytes_twos_comp()
        return pos_break, neg_break
    
    
    def get_enc_and_accel(self):
        self.write_byte(7)
        enc = self.read_two_bytes_twos_comp()
        accel = self.read_two_bytes_twos_comp()
        return enc, accel


    def zero_accel(self, N=10):
        accel_vect = np.zeros(N)

        for i in range(N):
            enc, a_i = self.get_enc_and_accel()
            accel_vect[i] = a_i

        a_ave = accel_vect.mean()
        self.accel_zero = a_ave
        self.a_zero_vect = accel_vect
        return self.accel_zero
        

    def set_gains(self, kp=3.0, kd=0.2):
        self.kp = kp
        self.kd = kd
        
        
    def sat(self, vin, mymax=255, mymin=-255):
        if vin > mymax:
            vout = mymax
        elif vin < mymin:
            vout = mymin
        else:
            vout = vin
        return vout


    def initialize(self):
        self.open()
        self.print_all_wait()
        self.start_inv_db()#inv deadband
        time.sleep(2)
        pb, nb = self.get_inv_db_res()
        return pb, nb
        
        
    def run_test(self, u, OL=False):
        N = len(u)
        
        self.flush()

        self.write_byte(2)#start new test
        check_byte = self.read_two_bytes()
        print('check_byte = %s' % check_byte)

        dt = 0.01
        e_prev = 0.0

        nvect = np.zeros(N,dtype=int)
        accel = np.zeros_like(nvect)
        enc = np.zeros_like(nvect)
        v = np.zeros_like(nvect)
        evect = np.zeros_like(nvect)

        self.nvect = nvect
        self.enc = enc
        self.v = v
        
        #TFz.input = evect
        #TFz.output = v

        msb_list = []
        lsb_list = []

        t0 = time.time()

        for i in range(N):
            if OL:
                vout = u[i]
            else:
                e = u[i] - enc[i-1]
                evect[i] = e
                e_dot = (e-e_prev)/dt
                e_prev = e#save for next time
                ### handle digcomp later
                #vout = int(TFz.calc_out(i))
                vout = int(e*self.kp+e_dot*self.kd)#PD calculation

            #time.sleep(0.0001)
            vout = self.sat(vout)
            v[i] = vout
            #print('%i, %f' % (i,vout))
            self.write_byte(1)#new n and voltage are coming
            self.write_two_byte_int(i)
            self.write_two_byte_int(int(vout))
            
            #time.sleep(0.0001)
            
            nvect[i] = self.read_two_bytes()
            #nvect[i] = serial_utils.Read_Two_Bytes(self.ser)
            #v_echo[i] = serial_utils.Read_Two_Bytes(ser)
            enc[i] = self.read_two_bytes_twos_comp()
            accel[i] = self.read_two_bytes_twos_comp()
            
            nl_check = self.read_byte()
            assert nl_check == 10, "newline problem"


        t1 = time.time()
        print('Total time = %0.4g' % (t1-t0))
        dt_ave = (t1-t0)/N
        print('average dt = %0.4g' % (dt_ave))

        time.sleep(0.1)
        self.write_byte(3)#stop test
        time.sleep(0.1)
        self.write_byte(3)#stop test
        
        anom = accel[0:10].mean()
        self.accel_shifted = accel - anom
        
        self.nvect = nvect
        self.accel = accel
        self.enc = enc
        self.u = u
        self.v = v
        self.evect = evect
        self.dt = dt
        self.t = dt*nvect
        
        
    def run_OL_test(self, amp=200, width=10, start=20, N=500):
        v_ol = np.zeros(N,dtype=int)
        v_ol[start:start+width] = amp
        self.run_test(v_ol, OL=True)
        return v_ol


    def run_double_pulse_test(self, amp=200, width=10, start=20, N=500,
                              start2=300):
        v_ol = np.zeros(N,dtype=int)
        v_ol[start:start+width] = amp
        v_ol[start2:start2+width] = -amp
        self.run_test(v_ol, OL=True)
        return v_ol
    

    def run_accel_fb_p_control(self, u, kp):

        if not hasattr(self, 'accel_zero'):
            self.zero_accel()
            
        N = len(u)
        
        self.flush()

        self.write_byte(2)#start new test
        check_byte = self.read_two_bytes()
        print('check_byte = %s' % check_byte)

        dt = 0.01
        e_prev = 0.0

        nvect = np.zeros(N,dtype=int)
        accel = np.zeros_like(nvect)
        enc = np.zeros_like(nvect)
        v = np.zeros_like(nvect)
        evect = np.zeros_like(nvect)
        theta_d = np.zeros_like(nvect)

        #accel_comp = digcomp.Digital_Compensator(b, a, accel, theta_d)
                                                 
        self.nvect = nvect
        self.enc = enc
        self.v = v
        
        #TFz.input = evect
        #TFz.output = v

        msb_list = []
        lsb_list = []

        t0 = time.time()

        for i in range(N):
            e = u[i] - enc[i-1] - theta_d[i-1]
            evect[i] = e
            e_dot = (e-e_prev)/dt
            e_prev = e#save for next time
            ### handle digcomp later
            #vout = int(TFz.calc_out(i))
            vout = int(e*self.kp+e_dot*self.kd)#PD calculation

            #time.sleep(0.0001)
            vout = self.sat(vout)
            v[i] = vout
            #print('%i, %f' % (i,vout))
            self.write_byte(1)#new n and voltage are coming
            self.write_two_byte_int(i)
            self.write_two_byte_int(int(vout))
            
            #time.sleep(0.0001)
            
            nvect[i] = self.read_two_bytes()
            #nvect[i] = serial_utils.Read_Two_Bytes(self.ser)
            #v_echo[i] = serial_utils.Read_Two_Bytes(ser)
            enc[i] = self.read_two_bytes_twos_comp()
            accel[i] = self.read_two_bytes_twos_comp() - self.accel_zero
            
            nl_check = self.read_byte()
            assert nl_check == 10, "newline problem"
            theta_d[i] = kp*accel[i]
            

        t1 = time.time()
        print('Total time = %0.4g' % (t1-t0))
        dt_ave = (t1-t0)/N
        print('average dt = %0.4g' % (dt_ave))

        time.sleep(0.1)
        self.write_byte(3)#stop test
        time.sleep(0.1)
        self.write_byte(3)#stop test
        
        anom = accel[0:10].mean()
        self.accel_shifted = accel - anom
        
        self.nvect = nvect
        self.accel = accel
        self.enc = enc
        self.u = u
        self.v = v
        self.evect = evect
        self.dt = dt
        self.t = dt*nvect
        self.theta_d = theta_d


    def run_accel_fb_test(self, u, a, b):
    
        if not hasattr(self, 'accel_zero'):
            self.zero_accel()

        N = len(u)

        self.flush()
        
        self.write_byte(2)#start new test
        check_byte = self.read_two_bytes()
        print('check_byte = %s' % check_byte)

        dt = 0.01
        e_prev = 0.0

        nvect = np.zeros(N,dtype=int)
        accel = np.zeros_like(nvect)
        enc = np.zeros_like(nvect)
        v = np.zeros_like(nvect)
        evect = np.zeros_like(nvect)
        theta_d = np.zeros_like(nvect)
        
        accel_comp = digcomp.Digital_Compensator(b, a, accel, theta_d)

        self.nvect = nvect
        self.enc = enc
        self.v = v
            
        #TFz.input = evect
        #TFz.output = v

        msb_list = []
        lsb_list = []

        t0 = time.time()

        for i in range(N):
            e = u[i] - enc[i-1] - theta_d[i-1]
            evect[i] = e
            e_dot = (e-e_prev)/dt
            e_prev = e#save for next time
            ### handle digcomp later
            #vout = int(TFz.calc_out(i))
            vout = int(e*self.kp+e_dot*self.kd)#PD calculation

            #time.sleep(0.0001)
            vout = self.sat(vout)
            v[i] = vout
            #print('%i, %f' % (i,vout))
            self.write_byte(1)#new n and voltage are coming
            self.write_two_byte_int(i)
            self.write_two_byte_int(int(vout))

            #time.sleep(0.0001)

            nvect[i] = self.read_two_bytes()
            #nvect[i] = serial_utils.Read_Two_Bytes(self.ser)
            #v_echo[i] = serial_utils.Read_Two_Bytes(ser)
            enc[i] = self.read_two_bytes_twos_comp()
            accel[i] = self.read_two_bytes_twos_comp() - self.accel_zero

            nl_check = self.read_byte()
            assert nl_check == 10, "newline problem"
            theta_d[i] = accel_comp.calc_out(i)


        t1 = time.time()
        print('Total time = %0.4g' % (t1-t0))
        dt_ave = (t1-t0)/N
        print('average dt = %0.4g' % (dt_ave))

        time.sleep(0.1)
        self.write_byte(3)#stop test
        time.sleep(0.1)
        self.write_byte(3)#stop test

        anom = accel[0:10].mean()
        self.accel_shifted = accel - anom

        self.nvect = nvect
        self.accel = accel
        self.enc = enc
        self.u = u
        self.v = v
        self.evect = evect
        self.dt = dt
        self.t = dt*nvect
        self.theta_d = theta_d


    
    def _find_max_and_min(self, pos_margin=10, neg_margin=10):
        ymax = np.max([self.v.max(),self.enc.max(), \
                       self.accel_shifted.max()]) + pos_margin
        ymin = np.min([self.v.min(),self.enc.min(), \
                       self.accel_shifted.min()]) - neg_margin
        return [ymin, ymax]

    
    def set_labels(self):
        plt.xlabel('Time (sec.)')
        plt.ylabel('Signal Amplitude (counts)')
        
    
    def plot_OL(self,legloc=1):
        plt.figure()
        t = self.t
        kwargs = {'linewidth':2}
        plt.plot(t, self.v, t, self.enc, t, self.accel_shifted, **kwargs)
        plt.legend(['v','$\\theta$','$\\ddot{x}$'],loc=legloc)
        ylims = self._find_max_and_min()
        plt.ylim(ylims)
        self.set_labels()
        
    
    def run_step_test(self, amp=50, start=20, N=500):
        u_step = np.zeros(N,dtype=int)
        u_step[start:] = amp
        self.run_test(u_step)
        return u_step
    
    
    def plot_CL(self, legloc=1):
        plt.figure()
        t = self.t
        kwargs = {'linewidth':2}
        plt.plot(t, self.u, t, self.v, t, self.enc, \
                 t, self.accel_shifted, **kwargs)
        plt.legend(['$u$','$v$','$\\theta$','$\\ddot{x}$'], \
                   loc=legloc)
        ylims = self._find_max_and_min()
        plt.ylim(ylims)
        self.set_labels()
        
        
    def _get_attrs(self, mylist):
        """Look up the possibly abbreviated attrs in mylist
        and get the full attribute name to use with getattr.  
        If an element of mylist is not in attr_dict, it is
        assumed to be a valid, full attribute name.
        attr_dict is a dictionary of abbreviated attr names."""
        attr_dict = {'a':'accel_shifted'}

        attr_list = []
        
        for key in mylist:
            if key in attr_dict:
                attr_list.append(attr_dict[key])
            else:
                attr_list.append(key)
                
        return attr_list
    

    def plot_list(self, mylist, legloc=1):
        """Convience function for plotting a list of attrs, 
        possibly using short cut names for the attributes."""
        label_dict = {'u':'$u$', 'v':'$v$', 'enc':'$\\theta$', \
                      'a':'$\\ddot{x}$'}
        
        attr_list = self._get_attrs(mylist)
        label_list = [label_dict[key] for key in mylist]
        
        plt.figure()
        t = self.t
        kwargs = {'linewidth':2}
        
        for attr, label in zip(attr_list, label_list):
            vect = getattr(self, attr)
            plt.plot(t, vect, label=label, **kwargs)
            
        plt.legend(loc=legloc)
        #ylims = self._find_max_and_min()
        #plt.ylim(ylims)
        self.set_labels()
        
        
    def gen_u_swept_sine(self,amp=10,fmax=7,fmin=0,tspan=20,\
                         dt=0.01,deadtime=0.2):
        u_ss,t = control_utils.create_swept_sine_signal(fmax=fmax, \
                                                        fmin=fmin, \
                                                        tspan=tspan, \
                                                        dt=dt, \
                                                        deadtime=deadtime)
        u_ss *= amp
        u_int = u_ss.astype(int)
        return u_int
    
        
    def gen_u_fixed_sine(self, f_fs, amp=5, num_cycles=10, T=None, dt=0.01):
        if T is None:
            T = (1.0/f_fs)*(num_cycles+0.25)
            
        t_fs = np.arange(0,T,dt)
        u_fs = amp*np.sin(2*np.pi*f_fs*t_fs)
        u_int = u_fs.astype(int)
        return u_int


    def save_data(self, filename=None):
        if filename is None:
            filename = self.find_filename()

        if os.path.exists(filename):
            print("exists: %s" % filename)
            return
        else:
            print("saving: %s" % filename)
            
        data = np.column_stack([self.nvect, self.t, self.u, \
                                self.v, self.enc, self.accel_shifted])
        mylabels = ['n','t','u','v','enc','accel (shifted)']
        txt_mixin.dump_delimited(filename, data, labels=mylabels)


    def save_find_filename(self, base_name):
        filename = self.find_filename(base_name)
        self.save_data(filename)
        

    def find_filename(self, base_name):
        pat = base_name + '_test_%0.2i.csv'
        for i in range(1,100):
            cur_name = pat % i
            if not os.path.exists(cur_name):
                break

        next_name = cur_name
        return next_name
    

class beam_digcomp_theta_fb(beam_serial_test):        
    def run_test(self, a, b, u=None, amp=None, start_ind=None, N=500):
        if u is None:
            u = self.create_step_u(amp=amp, N=N, start_ind=start_ind)
            
        N = len(u)
        
        self.flush()

        self.write_byte(2)#start new test
        check_byte = self.read_two_bytes()
        print('check_byte = %s' % check_byte)

        dt = 0.01
        e_prev = 0.0

        nvect = np.zeros(N,dtype=int)
        accel = np.zeros_like(nvect)
        enc = np.zeros_like(nvect)
        v = np.zeros_like(nvect)
        evect = np.zeros_like(nvect)

        theta_comp = digcomp.Digital_Compensator(b, a, evect, v)
            
        self.nvect = nvect
        self.enc = enc
        self.v = v
        
        #TFz.input = evect
        #TFz.output = v

        msb_list = []
        lsb_list = []

        t0 = time.time()

        for i in range(N):
            e = u[i] - enc[i-1]
            evect[i] = e
            e_dot = (e-e_prev)/dt
            e_prev = e#save for next time
            ### handle digcomp later
            vout = int(theta_comp.calc_out(i))
            #vout = int(e*self.kp+e_dot*self.kd)#PD calculation

            #time.sleep(0.0001)
            vout = self.sat(vout)
            v[i] = vout#<-- do I want my compensator to know saturation happend?
            #print('%i, %f' % (i,vout))
            self.write_byte(1)#new n and voltage are coming
            self.write_two_byte_int(i)
            self.write_two_byte_int(int(vout))
            
            #time.sleep(0.0001)
            
            nvect[i] = self.read_two_bytes()
            #nvect[i] = serial_utils.Read_Two_Bytes(self.ser)
            #v_echo[i] = serial_utils.Read_Two_Bytes(ser)
            enc[i] = self.read_two_bytes_twos_comp()
            accel[i] = self.read_two_bytes_twos_comp()
            
            nl_check = self.read_byte()
            assert nl_check == 10, "newline problem"


        t1 = time.time()
        print('Total time = %0.4g' % (t1-t0))
        dt_ave = (t1-t0)/N
        print('average dt = %0.4g' % (dt_ave))

        time.sleep(0.1)
        self.write_byte(3)#stop test
        time.sleep(0.1)
        self.write_byte(3)#stop test
        
        anom = accel[0:10].mean()
        self.accel_shifted = accel - anom
        
        self.nvect = nvect
        self.accel = accel
        self.enc = enc
        self.u = u
        self.v = v
        self.evect = evect
        self.dt = dt
        self.t = dt*nvect
