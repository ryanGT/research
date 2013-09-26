import xml.etree.ElementTree as ET
from xml.dom import minidom

import pdb

import xml_utils
from xml_utils import prettify

import copy

from IPython.core.debugger import Pdb

## element_params = {'dc_motor':['g','p','dt'], \
##                   'torsional_spring_damper':['k','c'], \
##                   'rigid_mass':['M','L','R','I'], \
##                   'beam':['mu','L','EI','t','N'], \
##                   }

## sorted_elements = sorted(element_params.iterkeys())

#print('sorted_elements = %s' % sorted_elements)

class bd_XML_element(xml_utils.xml_writer):
    def __init__(self, name=None, blocktype=None, params={}, \
                 pullist=[]):
        """pullist refers to a list of parameters that will be pulled
        from params before converting the xml representation to an
        actual elememt.  For the DT-TMM blocks, input was pulled to
        process later. I am not sure if there are any params to pull
        for block diagram objects.
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
        assert len(pullist) == 0, "pullist is not actually implemented yet"
        self.xml_attrs = ['name','blocktype'] + pullist#?how to handle params?
        self.xml_tag_name = 'block'


    def create_xml(self, root):
        my_elem = xml_utils.xml_writer.create_xml(self, root)
        params_xml = ET.SubElement(my_elem, 'params')
        for attr, val in self.params.iteritems():
            cur_xml = ET.SubElement(params_xml, attr)
            try:
                cur_xml.text = val.encode()
            except:
                cur_xml.text = str(val)
                


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


class block_diagram_system_parser(xml_utils.xml_parser):
    def get_blocks(self):
        self.block_xml_list = xml_utils.find_child_if_it_exists(self.root, \
                                                                'blocks')

    def parse(self):
        """Parse the block diagram system"""
        assert self.root.tag == 'block_diagram_system', \
               "This does not appear to be a valid block_diagram_system_parser xml file"
        self.get_blocks()

        block_parsers = []

        for block_xml in self.block_xml_list:
            cur_parser = block_parser(filename=None)
            cur_parser.set_root(block_xml)
            cur_parser.parse()
            block_parsers.append(cur_parser)
        
        self.block_parsers = block_parsers
        return self.block_parsers
    
        
    def convert(self):
        block_list = []

        for block_parser in self.block_parsers:
            cur_inst = block_parser.convert()
            block_list.append(cur_inst)

        self.block_list = block_list
        return block_list

    


class block_parser(xml_utils.xml_parser):
    def validate_and_get_body(self):
        """I don't know if I would ever really try to pass in an xml
        file that contains a single block, but this is what I did for
        individual figures and full GUI state xml files in
        data_vis_gui"""
        if self.root.tag == 'block':
            body = self.root
            return body
        elif self.root.tag == 'blocks':
            children = self.root.getchildren()
            #a figure file should have one child that is either a
            #bode_figure or a time_domain_figure
            assert len(children) == 1, "problem with the children in my xml file"
            body = children[0]
            return body
        else:
            raise ValueError, "Not sure how to proceed for a figure with tag %s" % self.root.tag


    def parse(self):
        """convert the XML associated with self to a list something
        ready to be converted to a bd_XML_element instance"""
        body = self.validate_and_get_body()
        self.name = xml_utils.get_child_value(body, 'name')
        self.blocktype = xml_utils.get_child_value(body, 'blocktype')
        self.params = xml_utils.get_params(body)


    def convert(self):
        """convert the list of parse dictionaries in
        :py:attr:`self.parsed_dicts to :py:class:`plot_description`
        instances and then also create and return a :py:class:`figure`
        instance"""
        bd_instance = bd_XML_element(name=self.name, \
                                     blocktype=self.blocktype, \
                                     params=self.params)
        return bd_instance

                 

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
