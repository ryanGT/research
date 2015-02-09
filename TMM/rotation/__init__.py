from __future__ import division
from scipy import cos, cosh, sin, sinh, array, r_, c_, exp, pi, dot, real, imag, zeros, eye, shape
#import MLab
import scipy
from scipy.linalg import det
import copy, re
import pdb

from TMM.TMMElement import TMMElementIHT, HT4, Transform8by8, rotationmat

from rwkmisc import symstr, SymstrMattoMaxima, get_first_key_that_exists
import rwkmisc



class rotation_matrix(TMMElementIHT):
    def __init__(self, axis, angle, **kwargs):
        elemtype='rot'
        params={'axis':axis, \
                'angle':angle}
        TMMElementIHT.__init__(self, elemtype, params, **kwargs)
        self.axis = axis
        self.angle = angle
        self.rotmat = rotationmat(axis, angle)


    def GetMat(self, s):
        return self.rotmat
