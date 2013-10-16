"""This module is being developed in conjunction with BD-MIL and the
paper I am writing for the IL/IN ASEE conference in Valparaiso.  For
now, the main idea is to create a block diagram system that can draw
its own TiKZ block diagram.  Eventually it will also generate Python
code for BD-MIL execution and simulation as well as drawing a png
block diagram probably using matplotlib with invisible axes."""
from numpy import *
from scipy import signal

import copy, basic_file_ops, time
import copy, basic_file_ops, time
import DTTMM, control_utils

import system_with_serial

header_str = r"""\documentclass[landscape,letterpaper,11pt]{article}
\usepackage[utf8x]{inputenc} % utf8 encoding
\usepackage[T1]{fontenc} % use T1 fonts
\usepackage{amsmath} % nice math symbols
\usepackage{bm} % bold math
\usepackage{color} % change text color

\usepackage{tikz}
\usetikzlibrary{calc,patterns,decorations.pathmorphing,decorations.markings}
\tikzstyle{emptynode}=[minimum width=0pt, inner sep=0pt, minimum height=0pt, draw=none]
%%%<
\usepackage{verbatim}
\usepackage[active,tightpage]{preview}
\PreviewEnvironment{tikzpicture}
\setlength\PreviewBorder{5pt}%
%%%>
\def \springlength {2.0cm}
\pgfmathparse{\springlength*3}
\let\damperlength\pgfmathresult
\def \groundX {0.0cm}
\def \groundwidth {4cm}
\def \masswidth {2.5cm}
\pgfmathparse{\masswidth/2}
\let\halfmasswidth\pgfmathresult
\def \wallwidth {0.35cm}
\pgfmathparse{\wallwidth/2}
\let\halfwallwidth\pgfmathresult
\pgfmathparse{\masswidth+0.7cm}
\let\wallheight\pgfmathresult

\usetikzlibrary{shapes,arrows}
\tikzstyle{block} = [draw, fill=blue!10, rectangle,
    minimum height=1.0cm, minimum width=1.0cm]
\tikzstyle{multilineblock} = [draw, fill=blue!10, rectangle,
    minimum height=1.25cm, minimum width=1.0cm,
    text width=2cm,text centered,midway]
\tikzstyle{sum} = [draw, fill=blue!20, circle, node distance=1.5cm]
\tikzstyle{input} = [coordinate]
\tikzstyle{output} = [emptynode]
\tikzstyle{myarrow} = [coordinate, node distance=1.5cm]
\tikzstyle{pinstyle} = [pin edge={to-,thin,black}]
\tikzstyle{serialnode} = [inner sep=0.5mm,rectangle,draw=black, fill=black]
\tikzstyle{serialline} = [draw, ->, ultra thick, densely dashed]

\def \capdist {0.75cm}
"""

## \node (input) {$u$};
## \node [sum, right of=input] (sum) {};
## \node [block, right of=sum, node distance=1.75cm] (controller)
##     {$G_c(s)$};
## \node [emptynode, below of=controller, node distance=\capdist]
##     (caption1) {Controller};
## \node [block, right of=controller] (plant)
##     {$G(s)$};
## \node [emptynode, below of=plant, node distance=\capdist]
##     (caption2) {Plant};
## \node (output) [right of=plant] {$y$};

## \draw [->] (input) -- (sum) node[pos=0.9, yshift=0.25cm] {\small{$+$}};
## \draw [->] (sum) -- (controller);
## \draw [->] (controller) -- (plant);
## \draw [->] (plant) -- (output) node [emptynode] (outarrow) [pos=0.5] {};

## \coordinate [below of=plant, node distance=2.0cm] (tmp);
## \draw [->] (outarrow) |- (tmp) -| (sum) node[pos=0.9, xshift=0.2cm] {{\small $-$}};

from IPython.core.debugger import Pdb

fmt_dict = {float:'%0.10g'}

import parse_DTTMM_xml


def val_to_str(val):
    if type(val) in fmt_dict.keys():
        fmt = fmt_dict[type(val)]
    else:
        fmt = '%s'
    str_out = fmt % val
    return str_out


class block_diagram_system(object):
    def __init__(self, blocks=None, wires=None, \
                 intermediate_points=None, \
                 annotations=None, \
                 default_node_distance='2.5cm'):
        self.blocks = blocks
        self.wires = wires
        self.intermediate_points = intermediate_points
        self.annotations = annotations
        self.default_node_distance = default_node_distance


    def get_block(self, name):
        matches = []

        for block in self.blocks:
            if block.name == name:
                matches.append(block)

        assert len(matches) > 0, "did not find a match for name %s" % name
        assert len(matches) == 1, "found more than one match for name %s" % name

        return matches[0]



    def get_blocks_by_type(self, type_name):
        matches = []

        for block in self.blocks:
            if block.blocktype == type_name:
                matches.append(block)

        return matches

        
    def sort_blocks(self):
        sorted_list = []
        unsorted_list = copy.copy(self.blocks)

        i = 0
        while i < len(unsorted_list):
            block = unsorted_list[i]
            if isinstance(block, source_block):
                unsorted_list.remove(block)
                sorted_list.append(block)
            else:
                i += 1


        i = 0
        while i < len(unsorted_list):
            block = unsorted_list[i]
            if isinstance(block, zoh_block):
                unsorted_list.remove(block)
                sorted_list.append(block)
            else:
                i += 1


        i = 0
        while i < len(unsorted_list):
            block = unsorted_list[i]
            if block.input in sorted_list:
                unsorted_list.remove(block)
                sorted_list.append(block)
            else:
                i += 1

        self.sorted_list = sorted_list
        self.unsorted_list = unsorted_list

        return self.sorted_list, self.unsorted_list
    
            

    def prep_for_sim(self, N, t, dt):
        self.N = N
        self.t = t
        self.dt = dt
        
        for block in self.blocks:
            block.prep_for_sim(N, t, dt)


    def Run_Simulation(self):
        for i in range(1,self.N):
            for block in self.sorted_list:
                block.sim_one_step(i)
        
            

    def add_simple_wires(self):
        """Search for blocks where the input and rel_block are the
        same.  If they are, assume a simple wire between the input and
        the block is needed."""
        simple_wires = []
        #Pdb().set_trace()
        for block in self.blocks:
            if block.input is not None:
                if block.input == block.rel_block:
                    curwire = simple_wire(block.input, block)
                    simple_wires.append(curwire)

        if self.wires is None:
            self.wires = simple_wires
        else:
            self.wires.extend(simple_wires)


    def append_wire(self, wire):
        if self.wires is None:
            self.wires = [wire]
        else:
            self.wires.append(wire)


    def to_tikz(self):
        opt_str =  "[every node/.style={font=\large}, node distance=%s,>=latex']" % self.default_node_distance
        header_list = header_str.split('\n')
        self.tikz_list = copy.copy(header_list)
        out = self.tikz_list.append
        extend = self.tikz_list.extend
        out('\\begin{document}')
        out('\\begin{tikzpicture}%s' % opt_str)

        for block in self.blocks:
            curline = block.to_tikz()
            out(curline)

        if self.intermediate_points:
            for point in self.intermediate_points:
                curline = point.to_tikz()
                out(curline)

        for wire in self.wires:
            curline = wire.to_tikz()
            out(curline)


        for annotation in annotations:
            curline = annotation.to_tikz()
            out(curline)

        out('\\end{tikzpicture}')
        out('\\end{document}')


    def save_tikz(self, filename):
        basic_file_ops.writefile(filename, self.tikz_list, append=False)


class exp_block_diagram_system(block_diagram_system, \
                               system_with_serial.system_with_serial):
    def __init__(self, *args, **kwargs):
        block_diagram_system.__init__(self, *args, **kwargs)
        self.ser = None


    def prep_for_exp(self, N, t, dt):
        self.N = N
        self.t = t
        self.dt = dt

        for block in self.blocks:
            block.prep_for_exp(N, t, dt)


    def Run_Exp(self):
        for i in range(1,self.N):
            for block in self.sorted_list:
                block.exp_one_step(i)


    def get_serial_block(self):
        ser_blocks = self.get_blocks_by_type('serial_plant')#<-- this is tricky if psoc and arduino blocks aren't exactly this type
        #                                                   #    I guess they are for now,just be careful
        assert len(ser_blocks) > 0, "did not find any serial_plant blocks"
        assert len(ser_blocks) == 1, "found more than one serial_plant blocks"
        ser_block = ser_blocks[0]
        return ser_block
        
    
    def _open_ser(self):
        system_with_serial.system_with_serial._open_ser(self)
        ser_block = self.get_serial_block()
        ser_block.set_ser_sys(self)
        
    def attach_ser(self, ser):
        self.ser = ser
        ser_block = self.get_serial_block()
        ser_block.set_ser_sys(self)
        
        
class psoc_exp_block_diagram_system(exp_block_diagram_system):
    def Get_IC(self):
        #self._open_ser()#I guess this is ok for PSoC, but makes me nervous with arduino
        self.WriteByte(113)
        self.IC = self.Read_Two_Bytes_Twos_Comp()
        return self.IC
        
    def Run_Exp(self):
        self.flush_ser()
        #Reset Theta? #<-- probably
        self.Get_IC()
        exp_block_diagram_system.Run_Exp(self)
        self.Stop_Test()
        
        
    def Stop_Test(self):
        #send a 0 voltage
        self.WriteByte(47)
        self.WriteInt(0)
        #stop serial on the psoc
        self.WriteByte(115)


    def Reset_Theta(self):
        ##self._open_ser()
        self.WriteByte(55)
        

class block(object):
    def __init__(self, name, label='', caption='', \
                 input=None, output=None, \
                 input_name=None, \
                 input_ind=None, \
                 position='abs', coordinates=(0,0), \
                 rel_block=None, node_distance=None, \
                 tikz_style=None, yshift=None, xshift=None,\
                 xml_attrs=None, blocktype=None, **kwargs):
        """Create a new block.  position must be one of 'abs', 'right
        of', \ 'left of', 'above of', or 'below of' (following TiKZ
        arguments).  coordinates is an (x,y) pair that is only used if
        position is 'abs'.  Similarly, rel_block and node_distance are
        only used if position is relative.  If node_distance is not
        None, it must be a string containing a valid latex distance
        such as '2.5cm'.  rel_block must be set before calling to_tikz
        if the position is relative.  rel_block should be set to a
        string containing the name of the reference block.

        input_name refers to the name of the input block.  This will
        be used along with xml system descriptions to later find the
        input block instance.

        input_ind refers to the index of the input block's outputs,
        for systems with more than one output."""
        self.name = name
        self.label = label
        self.caption = caption
        self.input = input
        self.output = output
        self.input_name = input_name
        self.input_ind = input_ind
        self.position = position.lower()
        self.coordinates = coordinates
        self.rel_block = rel_block
        self.node_distance = node_distance
        self.tikz_style = tikz_style
        self.yshift = yshift
        self.xshift = xshift
        self.xml_attrs = xml_attrs
        self.blocktype = blocktype

        for key, val in kwargs.iteritems():
            if not hasattr(self, key):#don't overwrite existing properties
                setattr(self, key, val)


    def _append_option(self, string):
        """append string to self.opt_str correctly handling whether or
        not a comma is needed before hand, i.e. if self.opt_str is not ''"""
        if self.opt_str:
            self.opt_str += ', ' + string
        else:
            self.opt_str = string
        return self.opt_str



    def get_input(self, i):
        if self.input_ind is not None:
            my_input = self.input.output[i][self.input_ind]
        else:
            my_input = self.input.output[i]
            
        return my_input


    def prep_for_sim(self, N, t=None, dt=None):
        if t is not None:
            self.t = t
        if dt is not None:
            self.dt = dt
        self.output = zeros(N)


    def prep_for_exp(self, N, t=None, dt=None):
        """Unless a child block overrides this method, assume that
        preparing for exp is the same as preparing for sim"""
        return self.prep_for_sim(N, t=t, dt=dt)
        

    def get_output(self, i):
        if hasattr(self, 'output'):
            return self.output[i]
        else:
            raise ValueError, 'do not know what to do about my output: %s' % self


    def exp_one_step(self, i):
        """For blocks that don't need to do anything different between
        simulation and experiment, call the sim method whenever this
        is called.

        Blocks whose exp method is different will need to override
        this method."""
        self.sim_one_step(i)
        

    def _build_opt_str(self):
        ## \node [block, right of=sum, node distance=1.75cm] (controller)
        ##     {$G_c(s)$};
        ## \node [emptynode, below of=controller, node distance=\capdist]
        ##     (caption1) {Controller};
        ## \node [block, right of=controller] (plant)
        ##     {$G(s)$};

        if self.tikz_style:
            self.opt_str = self.tikz_style
        else:
            self.opt_str = ''
        if self.position != 'abs':
            pos_str = self.position + '=' + self.rel_block.name
            self._append_option(pos_str)
        if self.node_distance is not None:
            nd_str = 'node distance=%s' % self.node_distance
            self._append_option(nd_str)
        not_none_list = ['yshift','xshift']
        for attr in not_none_list:
            value = getattr(self, attr)
            if value is not None:
                attr_str = '%s=%s' % (attr, value)
                self._append_option(attr_str)
        return self.opt_str


    def to_tikz(self):
        #\node (input) {$u$};
        self._build_opt_str()
        if self.opt_str:
            tikz_str = '\\node [%s] (%s)' % (self.opt_str, self.name)
        else:
            tikz_str = '\\node (%s)' % self.name
        if self.position == 'abs':
            tikz_str += ' at (%s,%s)' % self.coordinates
        label_str = ' {%s}' % self.label
        tikz_str += label_str
        tikz_str += ';'
        self.tikz_str = tikz_str
        return self.tikz_str
        #\node [sum, right of=input] (sum) {};


    def _get_input_name(self):
        if self.input is not None:
            return self.input.name


    def build_XML_dict(self):
        assert self.xml_attrs is not None, "xml_attrs is None, cannot build dict"
        mydict = {}
        for attr in self.xml_attrs:
            if attr == 'input':
                val = self._get_input_name()
            else:
                val = getattr(self, attr)
            val_str = val_to_str(val)
            mydict[attr] = val_str

        if self.input is not None:
            if 'input' not in self.xml_attrs:
                mydict['input'] = self._get_input_name()

        self.xml_dict = mydict
        return mydict


class source_block(block):
    def __init__(self, name, **kwargs):
        if kwargs.has_key('input'):
            msg = "The input for a source_block must be None."
            assert kwargs['input'] is None, msg
        if not kwargs.has_key('blocktype'):
            kwargs['blocktype'] = 'source'
        block.__init__(self, name, **kwargs)


    def prep_for_sim(self, N, t=None, dt=None):
        block.prep_for_sim(self, N, t, dt)
        self.build_u(t)


    def sim_one_step(self, i):
        self.output[i] = self.u[i]


class arbitrary_input(source_block):
    def __init__(self, name, u=None, **kwargs):
        source_block.__init__(self, name, blocktype='arbitrary_input', \
                              **kwargs)

    def build_u(self, *args, **kwargs):
        pass

    def set_u(self, u):
        self.u = u



class swept_sine(source_block):
    def __init__(self, name, fmin=0.0, fmax=20.0, \
                 dt=1.0/500, tspan=10.0, deadtime=0.1, \
                 t=None, amp=1.0, **kwargs):
        block.__init__(self, name, blocktype='swept_sine', **kwargs)
        self.fmin = fmin
        self.fmax = fmax
        self.dt = dt
        self.tspan = tspan
        self.deadtime = deadtime
        self.t = t
        self.amp = amp
        self.xml_attrs = ['fmin','fmax','dt','tspan','deadtime','t','amp']


    def build_u(self, t):
        self.u, t2 = control_utils.create_swept_sine_signal(fmax=self.fmax, \
                                                            fmin=self.fmin, \
                                                            t=t, \
                                                            deadtime=self.deadtime)
        return self.u



class finite_width_pulse(source_block):
    def __init__(self, name, t_on=0.0, t_off=1.0, amp=1.0, **kwargs):
        block.__init__(self, name, blocktype='finite_width_pulse', **kwargs)
        self.t_on = t_on
        self.t_off = t_off
        self.amp = amp
        

    def build_u(self, t):
        N = len(t)
        dt = t[1]-t[0]
        u = zeros(N)
        self.ind_on = self.t_on/dt
        self.ind_off = self.t_off/dt + 1
        u[self.ind_on:self.ind_off] = self.amp

        self.u = u
        return self.u


class step_input(finite_width_pulse):
    def __init__(self, name, t_on=0.0, amp=1.0, **kwargs):
        block.__init__(self, name, blocktype='step_input', **kwargs)
        self.t_on = t_on
        self.amp = amp


    def build_u(self, t):
        N = len(t)
        dt = t[1]-t[0]
        u = zeros(N)
        self.ind_on = self.t_on/dt
        u[self.ind_on:] = self.amp

        self.u = u
        return self.u


class summing_block(block):
    def __init__(self, name, label='', caption='', \
                 input=None, input2=None, signs=[1.0,-1.0], **kwargs):
        block.__init__(self, name, label=label, \
                       caption=caption, input=input, \
                       blocktype='summing_block', \
                       **kwargs)

        self.input2 = input2
        self.tikz_style = 'sum'
        self.signs = signs
        #\node [sum, right of=input] (sum) {};


    def sim_one_step(self, i):
        cur_out = 0.0
        self.output[i] = self.input.output[i]*self.signs[0] + \
                         self.input2.output[i]*self.signs[1]


class zoh_block(block):
    def __init__(self, name, label='', caption='', \
                 input=None, input_ind=None, **kwargs):
        block.__init__(self, name, label=label, \
                       caption=caption, input=input, \
                       blocktype='zoh_block', \
                       input_ind=input_ind, \
                       **kwargs)


    def sim_one_step(self, i):
        if self.input_ind is not None:
            self.output[i] = self.input.output[i-1][self.input_ind]
        else:
            self.output[i] = self.input.output[i-1]


class gain_block(block):
    def __init__(self, name, gain=1.0, **kwargs):
        block.__init__(self, name, blocktype='gain_block', \
                       **kwargs)
        self.gain = gain


    def sim_one_step(self, i):
        my_input = self.get_input(i)
        out = my_input*self.gain
        self.output[i] = out

        return out


class PID_block(block):
    def __init__(self, name, kp=1.0, kd=0.0, ki=0.0, **kwargs):
        block.__init__(self, name, blocktype='gain_block', \
                       **kwargs)
        self.kp = kp
        self.kd = kd
        self.ki = ki


    def prep_for_sim(self, N, t=None, dt=None):
        if t is not None:
            self.t = t
        if dt is not None:
            self.dt = dt
        self.output = zeros(N)
        self.evect = zeros(N)
        self.esum = zeros(N)
        self.kd_hat = self.kd/dt
        self.ki_hat = self.ki*dt
        

    def sim_one_step(self, i):
        e_i = self.get_input(i)
        self.evect[i] = e_i
        self.esum[i] = self.esum[i-1] + e_i
        ediff = e_i - self.evect[i-1]
        out = self.kp*e_i + self.ki_hat*self.esum[i] + self.kd_hat*ediff
        self.output[i] = out

        return out
    

class saturation_block(block):
    def __init__(self, name, max=200.0, min=None, **kwargs):
        block.__init__(self, name, blocktype='saturation_block', \
                       **kwargs)
        self.max = max
        if min is None:
            self.min = -1.0*self.max
        else:
            self.min = min

            
    def sim_one_step(self, i):
        my_input = self.get_input(i)
        if my_input > self.max:
            out = self.max
        elif my_input < self.min:
            out = self.min
        else:
            out = my_input

        self.output[i] = out

        return out


class digital_TF_block(block):
    def __init__(self, name, numz=[], denz=[], \
                 **kwargs):
        """The class represents a digital transfer function.  numz and
        denz are the coefficients of z.  It is assumed that the
        digital sampling time to generate numz and denz are the same
        as whatever is used in the system simulation."""
        block.__init__(self, name, blocktype='digital_TF_block', \
                       **kwargs)
        self.tikz_style = 'block'
        self.numz = numz
        self.denz = denz
        self.Nden = len(self.denz)


    def prep_for_sim(self, N, t=None, dt=None):
        if t is not None:
            self.t = t
        if dt is not None:
            self.dt = dt
        self.output = zeros(N)
        self.input_vector = zeros(N)


    def sim_one_step(self, i):
        out = 0.0
        self.input_vector[i] = self.get_input(i)

        for n, bn in enumerate(self.numz):
            if (i-n) > 0:
                out += self.input_vector[i-n]*bn

        for n in range(1, self.Nden):
            if (i-n) > 0:
                out -= self.output[i-n]*self.denz[n]

        out = out/self.denz[0]
        self.output[i] = out
        return out


class TF_block(digital_TF_block):
    def __init__(self, name, num=[], den=[], c2dmethod='tustin', \
                 **kwargs):
        """The class represents a continuous transfer function that
        will be converted to a digital one for simulation purposes.
        num and den are coefficients vectors or lists for the
        continous laplace polynomials in s.  These vectors will be
        converted to the coefficients of a digital transfer function
        using the c2d method c2dmethod (i.e. tustin, zoh, ...) when
        prep_for_sim is called with a specified dt."""
        block.__init__(self, name, blocktype='TF_block', \
                       **kwargs)
        self.tikz_style = 'block'
        self.num = num
        self.den = den
        self.c2dmethod = c2dmethod
        

    def c2d(self, dt):
        out = signal.cont2discrete((squeeze(self.num), squeeze(self.den)),dt,self.c2dmethod)
        self.numz = squeeze(out[0])
        self.denz = squeeze(out[1])
        self.Nden = len(self.denz)
        

    def prep_for_sim(self, N, t=None, dt=None):
        digital_TF_block.prep_for_sim(self, N, t=t, dt=dt)
        self.c2d(dt)
                
        ## \node [block, right of=sum, node distance=1.75cm] (controller)
        ##     {$G_c(s)$};


class DTTMM_block(block):
    """This block will be used to create a DT-TMM system from an XML
    file.  The XML file will (hopefully fully) specify the options for
    the block (i.e. everything needed to perform a simulation for one
    time step)."""
    def _get_sensor_list(self):
        sensors = []
        for sensor in self.dttmm_parser.parsed_sensor_list:
            sensors.append(sensor.name)
        self.sensors = sensors


    def _get_actuator_list(self):
        actuators = []
        for actuator in self.dttmm_parser.parsed_actuator_list:
            actuators.append(actuator.name)
        self.actuators = actuators


    def convert(self, numeric_parameters={}, dt=None):
        elemlist = self.dttmm_parser.convert(numeric_parameters)
        sys = DTTMM.DT_TMM_System_from_XML_four_states(self.dttmm_parser.elemlist, \
                                                       sensors=self.dttmm_parser.sensor_list, \
                                                       actuators=self.dttmm_parser.act_list, \
                                                       dt=dt)
        self.sys = sys


    def sim_one_step(self, i):
        cur_input = self.input.get_output(i)
        self.sys.run_sim_one_step(i, inputs=[cur_input], \
                                  int_case=self.int_case)

        for j, sensor in enumerate(self.sys.sensors):
            self.output[i,j] = sensor.signal_vect[i]



    def prep_for_sim(self, N, t=None, dt=None):
        block.prep_for_sim(self, N, t, dt)
        NS = len(self.sensors)
        self.output = zeros((N,NS))
        if dt is not None:
            self.sys.dt = dt
        self.sys.init_sim(N)
        

    def __init__(self, name, xmlpath, int_case=1, **kwargs):
        block.__init__(self, name, blocktype='DTTMM_block', **kwargs)
        self.tikz_style = 'block'
        self.xmlpath = xmlpath
        self.dttmm_parser = parse_DTTMM_xml.DTTMM_xml_parser(xmlpath)
        self._get_sensor_list()
        self._get_actuator_list()
        self.xml_attrs = ['xmlpath','sensors','actuators']
        self.int_case = int_case
        ## \node [block, right of=sum, node distance=1.75cm] (controller)
        ##     {$G_c(s)$};


class serial_plant_block_arduino(block):
    def read_serial(self, i):
        self.nvect[i] = self.ser_sys.Read_Two_Bytes_Twos_Comp()
        for j in self.num_sensors:
            self.output[i,j] = self.ser_sys.Read_Two_Bytes_Twos_Comp()

        newline = self.ser_sys.Read_Byte()
        assert newline == 10, "newline problem"


    def write_serial(self, i, v):
        #I guess I am assuming single input for now
        #
        #  - that seems unnecessarily limiting, but I don't have a
        #  multi-input system to test on right now
        #
        self.ser_sys.WriteByte(1)
        self.ser_sys.WriteInt(i)
        self.ser_sys.WriteInt(v)


    def exp_one_step(self, i):
        cur_input = self.input.get_output(i)
        self.write_serial(i, int(cur_input))
        self.read_serial(i)


    def prep_for_exp(self, N, t=None, dt=None):
        #do I really want to call this prep_for_exp?
        block.prep_for_sim(self, N, t, dt)
        NS = len(self.sensors)
        self.output = zeros((N,NS))
        self.num_sensors = NS
        self.nvect = zeros(N)
        

    def set_ser_sys(self, ser_sys):
        self.ser_sys = ser_sys

        
    def __init__(self, name, ser_sys=None, actuators=[], \
                 sensors=[], **kwargs):
        block.__init__(self, name, blocktype='serial_plant', **kwargs)
        self.name = name
        self.actuators = actuators
        self.sensors = sensors
        self.ser_sys = ser_sys


class serial_plant_block_psoc(serial_plant_block_arduino):
    def read_serial(self, i):
        self.nvect[i] = self.ser_sys.Read_Two_Bytes_Twos_Comp()
        for j in range(self.num_sensors):
            self.output[i,j] = self.ser_sys.Read_Two_Bytes_Twos_Comp()

        #### arduino does newline verification each time
        # newline = self.ser_sys.Read_Byte()
        # assert newline == 10, "newline problem"


    def write_serial(self, i, v):
        #I guess I am assuming single input for now
        #
        #  - that seems unnecessarily limiting, but I don't have a
        #  multi-input system to test on right now
        #
        self.ser_sys.WriteByte(47)#arduino uses 1
        #self.ser_sys.WriteInt(i)#arduino sends i here as a check
        self.ser_sys.WriteInt(v)


class output_block(source_block):
    def __init__(self, name, **kwargs):
        if kwargs.has_key('output'):
            msg = "The output for a source_block must be None."
            assert kwargs['output'] is None, msg
        block.__init__(self, name, **kwargs)


class intermediate_point(block):
    """This class is used to create coordinates for tikz drawings that
    will be part of more complicated wires."""
    def __init__(self, name, position='abs', coordinates=(0,0), \
                 rel_block=None, node_distance=None, \
                 tikz_style=None):
        block.__init__(self, name, position=position, \
                       coordinates=coordinates, \
                       rel_block=rel_block, \
                       node_distance=node_distance, \
                       tikz_style=tikz_style)


    def to_tikz(self):
        #\node (input) {$u$};
        self._build_opt_str()
        if self.opt_str:
            tikz_str = '\\coordinate [%s] (%s)' % (self.opt_str, self.name)
        else:
            tikz_str = '\\node (%s)' % self.name
        if self.position == 'abs':
            tikz_str += ' at (%s,%s)' % self.coordinates
        tikz_str += ';'
        self.tikz_str = tikz_str
        return self.tikz_str
        ## \coordinate [below of=plant, node distance=2.0cm] (tmp);


class simple_wire(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end


    def to_tikz(self):
        tikz_str = '\\draw [->] (%s) -- (%s);' % (self.start.name, self.end.name)
        self.tikz_str = tikz_str
        return self.tikz_str


class feedback_wire(simple_wire):
    def __init__(self, start, end, intermediate_point):
        self.start = start
        self.end = end
        self.intermediate_point = intermediate_point


    def to_tikz(self):
        tikz_str = '\\draw [->] (%s) |- (%s) -| (%s);' % (self.start.name, \
                                                          self.intermediate_point.name, \
                                                          self.end.name)
        self.tikz_str = tikz_str
        return self.tikz_str


        ## \draw [->] (outarrow) |- (tmp) -| (sum) node[pos=0.9, xshift=0.2cm] {{\small $-$}};


class annotation(block):
    """This class allows blocks to be annotated with text relative to them."""
    def __init__(self, myblock, text, position='below of', \
                 node_distance='0.5cm', \
                 xshift=None, yshift=None):
        self.block = myblock
        self.text = text
        block.__init__(self, text, label=text, rel_block=myblock, \
                       position=position, node_distance=node_distance, \
                       yshift=yshift, xshift=xshift)


    def to_tikz(self):
        #\node[below of=sum1, node distance=0.4cm, xshift=0.25cm] {{\small $-$}};
        #\node[above of=sum1, node distance=0.25cm, xshift=-0.4cm] {{\small $+$}};
        self._build_opt_str()
        if self.opt_str:
            tikz_str = '\\node [%s]' % self.opt_str
        else:
            tikz_str = '\\node '#this seems like it should never happen
        if self.position == 'abs':
                tikz_str += ' at (%s,%s)' % self.coordinates
        label_str = ' {%s}' % self.label
        tikz_str += label_str
        tikz_str += ';'
        self.tikz_str = tikz_str
        return self.tikz_str



# below is the code needed to draw the wires for my simple closed-loop
# block diagram

## \draw [->] (input) -- (sum) node[pos=0.9, yshift=0.25cm] {\small{$+$}};
## \draw [->] (sum) -- (controller);
## \draw [->] (controller) -- (plant);
## \draw [->] (plant) -- (output) node [emptynode] (outarrow) [pos=0.5] {};


## \draw [->] (outarrow) |- (tmp) -| (sum) node[pos=0.9, xshift=0.2cm] {{\small $-$}};


if __name__ == '__main__':
    input = source_block('input', label='$u$')
    sum1 = summing_block('sum1', position='right of', \
                           rel_block=input, input=input)
    controller = TF_block('controller', label='$G_c(s)$', \
                          position='right of', \
                          rel_block=sum1, input=sum1, node_distance='2cm')
    plant = TF_block('plant', label='$G(s)$', \
                     position='right of', \
                     rel_block=controller, input=controller)
    output = output_block('output', label='$y$', \
                          position='right of', \
                          rel_block=plant, input=plant)
    feedback_point = intermediate_point('fb1', position='below of', \
                                        rel_block=plant, \
                                        node_distance='1.5cm')
    take_off_point = intermediate_point('out_hat', position='left of', \
                                        rel_block=output, \
                                        node_distance='0.75cm')
    block_list = [input, sum1, controller, plant, output]
    intermediate_points = [feedback_point, take_off_point]
    fb_wire = feedback_wire(take_off_point, sum1, feedback_point)
    plus_sign = annotation(sum1, '{\small $+$}', \
                           position='above of', xshift='-0.4cm', \
                           node_distance='0.2cm')
    minus_sign = annotation(sum1, '{\small $-$}', \
                            position='below of', xshift='0.2cm', \
                            node_distance='0.4cm')
    Gc_caption = annotation(controller, 'Controller', \
                            position='below of', node_distance='0.7cm')
    plant_caption = annotation(plant, 'Plant', \
                               position='below of', node_distance='0.7cm')
    annotations = [plus_sign, minus_sign, Gc_caption, plant_caption]
    sys = block_diagram_system(block_list, \
                               intermediate_points=intermediate_points, \
                               annotations=annotations)
    sys.add_simple_wires()
    sys.append_wire(fb_wire)
    sys.to_tikz()
    import os
    outdir = '/home/ryan/siue/Research/papers/ASEE_IL_IN_BD_MIL/'
    outname = 'tikz_sys_test.tex'
    outpath = os.path.join(outdir, outname)
    sys.save_tikz(outpath)

