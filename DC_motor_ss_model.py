from scipy import *
from numpy.linalg import det, inv


class DC_motor_observer(object):
    def __init__(self, s0=-75, m=0.00103322, b=0.0164207, dt=1.0/500):
        """Create a digital state-space observer for the DC motor.  s0 is
        the observer pole in rad/sec, m is the mass of the motor model and
        b is the damping coefficient.  dt is the digital time step."""
        self.m = m
        self.b = b
        self.dt = dt
        self.s0 = s0
        self.p0 = exp(s0*dt)

        
    def find_F_and_G(self):
        m = self.m
        b = self.b
        dt = self.dt
        F = array([[1.0, m*(1.0 - exp(-dt*b/m))/b], [0, exp(-dt*b/m)]])
        G = array([[m*(1.0 - exp(dt*b/m) + \
                       dt*b*exp(dt*b/m)/m)*exp(-dt*b/m)/b**2],
                   [(1.0 - exp(-dt*b/m))/b]])
        self.F = F
        self.G = G
        return F, G


    def num_observer(self):
        if not hasattr(self, 'F'):
            self.find_F_and_G()
        h = atleast_2d(array([1,0]))#note that this is really h^T as far as
                                #Dorsey is concerned
        F = self.F
        G = self.G
        p0 = self.p0

        h1 = h[0,0]
        h2 = h[0,1]
        f11 = F[0,0]
        f12 = F[0,1]
        f21 = F[1,0]
        f22 = F[1,1]
        lhs = array([[h1,h2],[f22*h1-f21*h2,f11*h2-f12*h1]])
        rhs = array([F.trace()-2*p0,det(F)-p0**2])

        L = atleast_2d(dot(inv(lhs),rhs)).T
        self.L = L
        return L
