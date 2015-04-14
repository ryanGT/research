import control
import numpy

class Digital_Compensator(object):
    def __init__(self, num, den, input_vect=None, output_vect=None):
        self.num = num
        self.den = den
        self.input = input_vect
        self.output = output_vect
        self.Nnum = len(self.num)
        self.Nden = len(self.den)


    def calc_out(self, i):
        out = 0.0
        for n, bn in enumerate(self.num):
            out += self.input[i-n]*bn

        for n in range(1, self.Nden):
            out -= self.output[i-n]*self.den[n]
        out = out/self.den[0]
        return out


class Dig_Comp_from_ctime(Digital_Compensator):
    def __init__(self, TF_ctime, dt, method='tustin', **kwargs):
        TFz = control.matlab.c2d(TF_ctime, dt, method=method)
        num = numpy.squeeze(TFz.num)
        den = numpy.squeeze(TFz.den)
        Digital_Compensator.__init__(self, num, den, **kwargs)
        
        
