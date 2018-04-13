import platform, re, os
import serial_utils
import time
import txt_mixin
import matplotlib.pyplot as plt
import numpy as np

def data_to_float_array(data):
    mylist = data.split('\n')
    clean_list = mylist#[3:]
    clean_2 = list(filter(None,clean_list))
    
    start_ind = None
    
    for i, row in enumerate(clean_2):
        if row.find("1,") == 0:
            start_ind = i
            
    end_ind = None
    
    for i, row in enumerate(clean_2):
        if row.find("input") > -1:
            end_ind = i
            
    clean_3 = clean_2[start_ind:end_ind]
    nested_list = [row.split(',') for row in clean_3]
    str_array = np.array(nested_list)
    float_array = str_array.astype(float)
    return float_array
#return str_array


class SBR_serial_test(serial_utils.serial_test):
    def guess_os(self):
        plat_str = platform.platform().lower()

        myos = 'win'

        if 'linux' in plat_str:
            myos = 'linux'
        elif 'darwin' in plat_str:
            myos = 'mac'

        self.os = myos
        return self.os

    
    def guess_portname(self):
        if not hasattr(self, 'os'):
            self.guess_os()

        # I could try to search for existing ports using os.path.exists
        # or something on Mac and Windows
        if self.os == 'linux':
            portname = serial_utils.find_portname_RPi()
        elif self.os == 'mac':
            portname = '/dev/tty.usbmodem1421'
            #portname = '/dev/cu.usbmodem1411'
        else:
            portname = 'COM3'

        self.portname = portname
        return self.portname
        
    
    def __init__(self, portname=''):
        if not portname:
            self.guess_portname()
        else:
            self.portname = portname
        serial_utils.serial_test.__init__(self, \
                                          self.portname, \
                                          baudrate=115200, \
                                          timeout=1)
        


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


    def send_char_int_wait(self, mychar, myint):
        """Send a char for menu input followed by an integer that
        sets some parameter; then wait for a response"""
        self.write_char(mychar)
        self.write_int(myint)
        self.print_all_wait()


    def send_char_float_wait(self, mychar, myfloat, fmt="%0.4f"):
        """Send a char for menu input followed by an integer that
        sets some parameter; then wait for a response"""
        self.write_char(mychar)
        self.write_float(myfloat, fmt=fmt)
        self.print_all_wait()
        
        
    def run_test(self, char="1", N=5):
        self.write_char(char)
        time.sleep(0.2)
        data = self.get_data()

        for i in range(N):
            print('getting data, i = %i' % i)
            new_data = self.get_data()
            data += new_data
            time.sleep(0.1)
            
        
        self.data = data
        self.data_array = data_to_float_array(data)


    def get_data_labels(self, arduino_path):
        """Extract the column labels for an Arduino test based
        on the assumption that the printing is handled by a 
        function called print_line_serial"""
        myfile = txt_mixin.txt_file_with_list(arduino_path)
        ind1 = myfile.findnext('void print_line_serial')
        ind2 = myfile.findnext('}',ind1)
        self.func_list = myfile.list[ind1:ind2]

        assert (self.func_list[0].find('void') > -1), "list did not start with void"
        self.func_list.pop(0)
        p_com = re.compile('[ \t]*//')
        filt1 = [item for item in self.func_list \
                          if not p_com.search(item)]

        filt2 = [item for item in filt1 if item.strip()]

        for i in range(-1,-5,-1):
            curline = filt2[i].strip()
            if curline.find("Serial.print('\\n')") == 0:
                filt2.pop(i)

        # get rid of printing commas
        filt3 = []
        for line in filt2:
            lineout = line.strip()
            lineout = lineout.replace('Serial.print(",");','')
            filt3.append(lineout)
        self.filt_list = filt3

        # get variable names
        labels = []
        p_var_name = re.compile("Serial.print\((.*)\);")

        for line in filt3:
            q = p_var_name.match(line)
            labels.append(q.group(1))

        self.labels = labels
        return self.labels


    def assign_cols_to_attrs(self):
        if not hasattr(self, 'labels'):
            self.get_data_labels()

        msg = "This method can only be called after a test has been run"
        assert hasattr(self, 'data_array'), msg

        for i, label in enumerate(self.labels):
            curcol = self.data_array[:,i]
            setattr(self, label, curcol)

        # assign t if it does not exist and t_ms does
        if 't' not in self.labels:
            if hasattr(self, 't_ms'):
                self.t = self.t_ms/1000.0
        


    def plot(self, attr_list, ax=None, figsize=(9,6), set_xlabel=True):
        # need to handle ax = None for plot2 and plot3
        # to use this
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            
        if type(attr_list) == str:
            attr_list = [attr_list]
            plt.figure()

        t = self.t

        for attr in attr_list:
            cur_vect = getattr(self, attr)
            ax.plot(t,cur_vect, label=attr)

            if set_xlabel:
                ax.set_xlabel('Time (sec.)')

            ax.set_ylabel('Signal Amplitude')


    def plot_seperate(self, attr_list):
        N = len(attr_list)
        fig = plt.figure()

        t = self.t
        
        for i in range(N):
            attr = attr_list[i]
            cur_vect = getattr(self, attr)
            v = i+1
            ax = fig.add_subplot(N, 1, v)
            ax.plot(t, cur_vect, label=attr)

            
    def plot2(self, attr_list1, attr_list2, figsize=(9,6)):
        """Plot two groups of attributes/signals using subplot(211)
        and subplot(212)"""
        fig = plt.figure(figsize=figsize)
        ax1 = fig.add_subplot(211)
        self.plot(attr_list1, ax1, set_xlabel=False)
        ax2 = fig.add_subplot(212)
        self.plot(attr_list2, ax2)


    
    def plot3(self, attr_list1, attr_list2, attr_list3, figsize=(9,6)):
        fig = plt.figure(figsize=figsize)
        ax1 = fig.add_subplot(311)
        self.plot(attr_list1, ax1, set_xlabel=False)
        ax2 = fig.add_subplot(312)
        self.plot(attr_list2, ax2, set_xlabel=False)
        ax3 = fig.add_subplot(313)
        self.plot(attr_list3, ax3)


    def save(self, filename, folder=''):
        outpath = os.path.join(folder,filename)
        if not os.path.exists(outpath):
            txt_mixin.dump_delimited(outpath,self.data_array,labels=self.labels)
            print('saved: %s' % outpath)
        else:
            print('already exists: %s' % outpath)
