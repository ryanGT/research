import xml.etree.ElementTree as ET
from xml.dom import minidom

import pdb

ERR_TOL = 1e-5 # floating point slop for peak-detection


def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")


element_params = {'dc_motor':['g','p','dt'], \
                  'torsional_spring_damper':['k','c'], \
                  'rigid_mass':['M','L','R','I'], \
                  'beam':['mu','L','EI','t','N'], \
                  }

sorted_elements = sorted(element_params.iterkeys())

#print('sorted_elements = %s' % sorted_elements)

class DTTMM_XML_element(object):
    def __init__(self, name=None, elemtype=None, params={}):
        if not bool(name) and bool(elemtype) and bool(params):
            print('all inputs must be non-empty:')
            print('name: %s' % name)
            print('elemtype: %s' % elemtype)
            print('params: %s' % params)
            raise ValueError, "problems with non-empty inputs"
        self.name = name
        self.elemtype = elemtype
        self.params = params


class sensor_XML_element(object):
    def __init__(self, name=None, sensor_type=None, signal=None, \
                 elem1=None, elem2=None):
        """name is the name of the sensor; sensor_type is either abs or diff;
        signal should be one of x, xdot, xddot, theta, thetadot,
        thetaddot, V or M; elem1 is the only element used for sensor_type
        abs; if sensor_type is diff, the sensor signal will be calculated
        from elem1-elem2"""
        if not bool(name) and bool(sensor_type) and bool(signal) \
               and bool(elem1):
            print('only elem2 is optional; all other inputs must be non-empty:')
            print('name: %s' % name)
            print('sensor_type: %s' % sensor_type)
            print('signal: %s' % signal)
            print('elem1: %s' % elem1)
            raise ValueError, "problems with non-empty inputs"

        if sensor_type == 'diff' and not bool(elem2):
            raise ValueError, "if sensor sensor_type is diff, elem2 must not be empty"

        self.name = name
        self.sensor_type = sensor_type
        self.signal = signal
        self.elem1 = elem1
        self.elem2 = elem2


class actuator_XML_element(object):
    def __init__(self, name=None, signal=None):
        """name is the name of the actuator element and signal is the
        variable that will be set at each step through the simulation
        loop."""
        if not bool(name) and bool(signal):
            print('name and signal must both be non-empty:')
            print('name: %s' % name)
            print('signal: %s' % signal)
        self.name = name
        self.signal = signal
