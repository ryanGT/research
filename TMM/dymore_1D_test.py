import TMM, Dymore
from scipy import pi, sqrt, arange, vectorize, exp, c_, array, transpose, real, imag, rand, cos, sin, sinh, cosh, argmax, arange, eye, zeros
import pylab
import scipy
import pdb
import time
import SAMII
import rwkascii
from rwkmisc import prettymat, colwise
from rwkos import FindFullPath
import shutil, os, sys, glob
import TMM.beam
reload(TMM.beam)
import TMM.rigid
reload(TMM.rigid)
from TMM.beam import BeamElement
from TMM.rigid import RigidMass
from rwkdataproc import datastruct
import re
import rwkparse
sys.path.append('E:\\pythonscripts\\dymore_parsing\\mode_shape_parser\\')
ts=time.time()

beamonly=0
fi=1
#======================
#
#   Beam Properties
#
#======================
EI=3.39137e+005
mu=5.72815e+000
L=7.2
mybeam=BeamElement({'EI':EI,'L':L,'mu':mu},maxsize=4)

#pdb.set_trace()
mysys=TMM.TMMSystem([mybeam],bcbase='fixed',bcend='free')

if not beamonly:
    m=20.
    a=0.5
    b=0.2
    c=0.1
    r2=0.1
    d2=a-r2
    I11=1./12.*m*(a**2+c**2)
    print('I11='+str(I11))
    I22=1./12.*m*(b**2+c**2)
    print('I22='+str(I22))
    I33=1./12.*m*(b**2+a**2)
    print('I33='+str(I33))
    I33o=I33+m*d2**2
    print('I33o='+str(I33o))
    I11o=I11+m*d2**2
    print('I11o='+str(I11o))
    I22o=I22+m*d2**2

#dddddddddddddddddddddddddddddddd
#
#   Dymore
#
#dddddddddddddddddddddddddddddddd
if beamonly:
    dymorename='beam_only.dym'
    pointnames=['PointA','PointB']
    coords=[[0.,0.,0.],[L,0.,0.]]
else:
    dymorename='beam_w_rigid_link.dym'
    pointnames=['PointA','PointB','PointC']
    coords=[[0.,0.,0.],[L,0.,0.],[L+a,0.,0.]]
problemname,ext=os.path.splitext(dymorename)
dymorepath=os.path.join(os.getcwd(),problemname)
if not os.path.exists(dymorepath):
    os.mkdir(dymorepath)
dymorefullpath=os.path.join(dymorepath,dymorename)
mylist=rwkascii.rwkdymorelist([],dymorefullpath)
mylist.readfile('dymore_header.txt')
mylist.AppendPoints(pointnames,coords)
if not beamonly:
    RigidLinkparams={'m':m,'L':a,'r':d2,'I':[I11,I22,I33]}
    RigidLink=RigidMass(RigidLinkparams,maxsize=4)
    mylist.AppendMesh('RigidLinkMesh',2,3)
    mylist.AppendLine('RigidLinkCurve','PointB','PointC','RigidLinkMesh')
    mylist.AppendMassProp('RigidLinkMass',m,[d2,0.,0.],[I11,I22o,I33o])
    mylist.AppendRigidBody('RigidLink','Beam1','PointB','FreeC','PointC','RigidLinkMass','RigidLinkCurve','RigidLinkGraphParams')
    mylist.AppendGraphParams('RigidLinkGraphParams','Symbol','Yes')

if beamonly:
    mylist.readfile('beam_only_file_1D.txt')
else:
    mylist.readfile('beam_file_1D.txt')
mylist.append('@BOUNDARY_CONDITION_DEFINITION{')
clampA=Dymore.BoundaryCondition('ClampA','PointA','Beam1',[1,1,1],[1,1,1])
mylist.extend(clampA.ToList())
if beamonly:
    freeB=Dymore.BoundaryCondition('FreeB','PointB','Beam1',[0,0,0],[0,0,0])
    mylist.extend(freeB.ToList())
else:
    freeC=Dymore.BoundaryCondition('FreeC','PointC','RigidLink',[0,0,0],[0,0,0])
    mylist.extend(freeC.ToList())
mylist.append('}')
mylist.readfile('dymore_tail.txt')
mylist.tofile()
dymfigspath=os.path.join(dymorepath,'FIGURES')
if not os.path.exists(dymfigspath):
    os.mkdir(dymfigspath)
#os.startfile(os.path.join(dymorefullpath))
#pdb.set_trace()
curdir=os.getcwd()
os.chdir(dymorepath)
print('before startfile')
os.startfile(dymorename)
print('after startfile')
pathwoext,ext=os.path.splitext(dymorefullpath)
dymputlist=rwkascii.rwkdymorelist([],pathwoext+'.dyp')
dymputlist.append('')
dymputlist.append(' Problem name: '+pathwoext)
dymputlist.tofile()
os.chdir(curdir)
#dddddddddddddddddddddddddddd
#
#   Parse Dymore
#
#dddddddddddddddddddddddddddd
from dymore_mode_shape_parser import MainParsingFunction
outpath=os.path.join(dymorepath,'parsed')
if not os.path.exists(outpath):
    os.mkdir(outpath)
meshname=problemname+'~MSH~.grf'
evcname=problemname+'~EVC~.grf'
meshpath=os.path.join(dymorepath,meshname)
evcpath=os.path.join(dymorepath,evcname)
basedymname=problemname+'_dym.txt'
basepath=os.path.join(outpath,basedymname)
indepind=0#x
myoutdofs=range(6)#all
MainParsingFunction(meshpath,basepath,evcpath,indepind,myoutdofs)

rot_b0=TMM.TMMElement('rot',{'axis':'z','angle':90})
if beamonly:
    mysys=TMM.TMMSystem([mybeam],'free','fixed')
else:
    mysys=TMM.TMMSystem([mybeam,RigidLink],'free','fixed')
mesh=mysys.CreateMesh()
basename=problemname+'_TMM_'
meshdict={'beammesh':mysys.sysmesh}
meshfiles=glob.glob(basename+'mesh'+'*')
for file in meshfiles:
        os.remove(file)

scipy.io.save(basename+'mesh',meshdict)
print('mesh=')
#--------------------------
#
#   Dymore Mesh
#--------------------------
dymdir=dymorepath
destdir=dymorepath
meshname=problemname+'~MSH~.grf'
evcname=problemname+'~EVC~.grf'

dymmeshlist=[]
basemesh=datastruct()
basemesh['label']='base'
basemesh['R']=eye(3,typecode='d')
basemesh['mesh']=zeros((1,3),'f')
dymmeshlist.append(basemesh)

meshpath=os.path.join(dymdir,meshname)
meshlines=rwkascii.rwkreadfile(meshpath)
p=re.compile('[\s]*-*\d\.[\d]+')
l=re.compile('[\s]*[a-zA-z]+[\w]*')
first=1
for x,line in enumerate(meshlines):
    tempout=p.match(line)
    if tempout:
#            print("Keeping:"+line)
        curmesh['mesh'].append(line)
    labeltest=l.match(line)
#    print('x='+str(x))
#    print('n-1='+str(len(meshlines)-1))
#    print('test='+str(x==len(meshlines)-1))
    if labeltest or x==len(meshlines)-1:
        if not first:
            curmesh['mesh']=rwkparse.matlines2nestedlist(curmesh['mesh'])
            curmesh['mesh']=array(curmesh['mesh'])
            dymmeshlist.append(curmesh)
        else:
            first=0
        curmesh=datastruct()
        curmesh['label']=line.strip()
        curmesh['mesh']=[]

itmesh=meshdict['beammesh']
for x,cm in enumerate(itmesh):
    if cm.elemtype==3:
        rotmesh=datastruct()
        rotmesh['label']='rotation'
        rotmesh['mesh']=rowwise(dymmeshlist[x-1]['mesh'][-1,:])
        dymmeshlist.insert(x,rotmesh)

for dm, tm in zip(dymmeshlist,itmesh):
    dm['elemtype']=tm['elemtype']
    dm['R']=tm['R']
    dm['mesh']=dm['mesh'][:,0:3]
dymmesh=dm['mesh']*L
dymx=dymmesh[:,0]

if beamonly:
    fn=os.path.join(outpath,'beam_only_dym_z.txt')
else:
    fn=os.path.join(outpath,'beam_w_rigid_link_dym_z.txt')
dymres=scipy.io.read_array(fn)
if beamonly:
    dymx=dymres[:,0]*L
else:
    dymx=dymres[:,0]*(L+a)
Ntoplot=3
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
    pylab.figure(fi)
    pylab.cla()
    m=n+1
    curdym=dymres[:,m]
    dind=argmax(abs(curdym))
    dmax=curdym[dind]
    curdym=curdym*0.2/dmax
    pylab.plot(dymx,curdym,'r--')
    pylab.plot(myx,disps)
    pylab.plot(xvect,curAY)
    pylab.legend(('FEA','TMM','Analytic'))
    fi+=1

print(str(mesh))
pylab.show()
