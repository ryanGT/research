import xml.etree.ElementTree as ET
from xml.dom import minidom

import pdb

from DTTMM_xml import prettify

import copy

from IPython.core.debugger import Pdb

## element_params = {'dc_motor':['g','p','dt'], \
##                   'torsional_spring_damper':['k','c'], \
##                   'rigid_mass':['M','L','R','I'], \
##                   'beam':['mu','L','EI','t','N'], \
##                   }

## sorted_elements = sorted(element_params.iterkeys())

#print('sorted_elements = %s' % sorted_elements)

class bd_XML_element(object):
    def __init__(self, name=None, blocktype=None, params={}, \
                 pullist=['input']):
        """pullist refers to a list of parameters that will be pulled
        from params before converting the xml representation to an
        actual elememt.  The input param will be pulled to process at
        a later time. DT-TMM system elements will have their sensors
        and actuators params pulled because those were only included
        for information purposes (the full actuator and sensor
        specifications are in the DT-TMM xml file).
        """
        if not bool(name) and bool(blocktype) and bool(params):
            print('all inputs must be non-empty:')
            print('name: %s' % name)
            print('blocktype: %s' % blocktype)
            print('params: %s' % params)
            raise ValueError, "problems with non-empty inputs"
        self.name = name
        self.blocktype = blocktype
        self.params = params
        self.pullist = pullist


    def clean_params(self):
        clean_params = copy.copy(self.params)
        for attr in self.pullist:
            if clean_params.has_key(attr):
                val = clean_params.pop(attr)
                setattr(self, attr, val)
                

        #attempt converting params to float
        clean2 = {}
        for key, val in clean_params.iteritems():
            try:
                val_out = float(val)
            except:
                val_out = val
            if val_out == 'None':
                val_out = None
            clean2[key] = val_out

        self.clean_params = clean2


class DTTMM_block(bd_XML_element):
    def __init__(self, name=None, blocktype=None, params={}):
        pullist = ['input','sensors','actuators']
        bd_XML_element.__init__(self, name=name, blocktype=blocktype, \
                                params=params, pullist=pullist)


                 
                 

## class sensor_XML_element(object):
##     def __init__(self, name=None, sensor_type=None, signal=None, \
##                  elem1=None, elem2=None):
##         """name is the name of the sensor; sensor_type is either abs or diff;
##         signal should be one of x, xdot, xddot, theta, thetadot,
##         thetaddot, V or M; elem1 is the only element used for sensor_type
##         abs; if sensor_type is diff, the sensor signal will be calculated
##         from elem1-elem2"""
##         if not bool(name) and bool(sensor_type) and bool(signal) \
##                and bool(elem1):
##             print('only elem2 is optional; all other inputs must be non-empty:')
##             print('name: %s' % name)
##             print('sensor_type: %s' % sensor_type)
##             print('signal: %s' % signal)
##             print('elem1: %s' % elem1)
##             raise ValueError, "problems with non-empty inputs"

##         if sensor_type == 'diff' and not bool(elem2):
##             raise ValueError, "if sensor sensor_type is diff, elem2 must not be empty"

##         self.name = name
##         self.sensor_type = sensor_type
##         self.signal = signal
##         self.elem1 = elem1
##         self.elem2 = elem2
