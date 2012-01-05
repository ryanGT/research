"""This module is being developed in conjunction with BD-MIL and the
paper I am writing for the IL/IN ASEE conference in Valparaiso.  For
now, the main idea is to create a block diagram system that can draw
its own TiKZ block diagram.  Eventually it will also generate Python
code for BD-MIL execution and simulation as well as drawing a png
block diagram probably using matplotlib with invisible axes."""

import copy, basic_file_ops

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

from IPython.Debugger import Pdb

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



class block(object):
    def __init__(self, name, label='', caption='', \
                 input=None, output=None, \
                 position='abs', coordinates=(0,0), \
                 rel_block=None, node_distance=None, \
                 tikz_style=None, yshift=None, xshift=None):
        """Create a new block.  position must be one of 'abs', 'right
        of', \ 'left of', 'above of', or 'below of' (following TiKZ
        arguments).  coordinates is an (x,y) pair that is only used if
        position is 'abs'.  Similarly, rel_block and node_distance are
        only used if position is relative.  If node_distance is not
        None, it must be a string containing a valid latex distance
        such as '2.5cm'.  rel_block must be set before calling to_tikz
        if the position is relative.  rel_block should be set to a
        string containing the name of the reference block."""
        self.name = name
        self.label = label
        self.caption = caption
        self.input = input
        self.output = output
        self.position = position.lower()
        self.coordinates = coordinates
        self.rel_block = rel_block
        self.node_distance = node_distance
        self.tikz_style = tikz_style
        self.yshift = yshift
        self.xshift = xshift


    def _append_option(self, string):
        """append string to self.opt_str correctly handling whether or
        not a comma is needed before hand, i.e. if self.opt_str is not ''"""
        if self.opt_str:
            self.opt_str += ', ' + string
        else:
            self.opt_str = string
        return self.opt_str


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



class source_block(block):
    def __init__(self, name, **kwargs):
        if kwargs.has_key('input'):
            msg = "The input for a source_block must be None."
            assert kwargs['input'] is None, msg
        block.__init__(self, name, **kwargs)



class summing_block(block):
    def __init__(self, name, label='', caption='', \
                 input=None, input2=None, **kwargs):
        block.__init__(self, name, label=label, \
                       caption=caption, input=input, \
                       **kwargs)
        
        self.input2 = input2
        self.tikz_style = 'sum'
        #\node [sum, right of=input] (sum) {};


class TF_block(block):
    def __init__(self, name, **kwargs):
        block.__init__(self, name, **kwargs)
        self.tikz_style = 'block'
        ## \node [block, right of=sum, node distance=1.75cm] (controller) 
        ##     {$G_c(s)$};


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
    
