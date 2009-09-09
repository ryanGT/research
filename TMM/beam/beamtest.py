from __future__ import division
import TMM
reload(TMM)
import TMM.beam
reload(TMM.beam)
from TMM.beam import BeamElement
from scipy import sqrt, c_, rand, imag, array, cos, cosh, real,argmax, arange, sin, sinh
import pdb
import pylab

#One Dimensional (4x4) test
EI=30000.
L=4.75
mu=0.25
mybeam=BeamElement({'EI':EI,'L':L,'mu':mu},maxsize=4)

#pdb.set_trace()
mysys=TMM.TMMSystem([mybeam],bcbase='fixed',bcend='free')



#Analytic solution
const=sqrt(EI/(mu*L**4))
myvect=c_[1.875,4.694,7.855]
myvect=myvect**2*const
#w1a=1.875**2*const
#w2a=4.694**2*const
#w3a=7.855**2*const
print('Analytic Eigenvalues:')
for k,ent in enumerate(myvect):
    q=k+1
    print('w%da='%q+str(ent))

evect=rand(3)
evect+=0.5-1
evect=evect*10
ivect=evect+myvect

tmmevect=array(map(mysys.FindEig,ivect*1.0j))
for q,cureig in enumerate(tmmevect):
    print('TMM eig %d= '%(q+1)+str(cureig[1])) 
errors=myvect-tmmevect[:,1]
perrors=abs(errors/myvect*100.)

print('max. percent error (Analytic-TMM)='+str(max(perrors)))

beta=(myvect**2*mu/EI)**0.25
beta2=(tmmevect[:,1]**2*mu/EI)**0.25
test=cos(beta*L)*cosh(beta*L)+1
test2=cos(beta2*L)*cosh(beta2*L)+1
xvect=arange(0,1.01,0.01)
xvect=xvect*L
for n,cureig,curbeta in zip(range(len(tmmevect)),tmmevect,beta):
    curAY=(sin(curbeta*L)-sinh(curbeta*L))*(sin(curbeta*xvect)-sinh(curbeta*xvect))+(cos(curbeta*L)+cosh(curbeta*L))*(cos(curbeta*xvect)-cosh(curbeta*xvect))
    mind=argmax(abs(curAY))
    mymax=curAY[mind]
    curAY=curAY*0.2/mymax
    disps,angles,modedict=mysys.FindModeShape(cureig)
    disps=real(disps)[:,0]
    mind2=argmax(abs(disps))
    mymax2=disps[mind2]
    disps=disps*0.2/mymax2
    mymesh=mysys.CreateMesh()
    myx=mymesh[:,0]
    pylab.ioff()
    pylab.figure(n+1)
    pylab.cla()
    pylab.plot(myx,disps)
    pylab.plot(xvect,curAY)
pylab.show()
