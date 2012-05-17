from scipy import shape,argmax, c_, zeros, dot
from rwkmisc import reverse
from rwkmisc import RowSwap
import pdb
from IPython.core.debugger import Pdb
import copy

def uppertri(matin,vectin=[],eps=1e-8):
    mat=copy.deepcopy(matin)
    nr=shape(mat)[0]
    if vectin:
        vect=copy.deepcopy(vectin)
    else:
        vect=zeros((nr,1),'d')
    mat=c_[mat,vect]
#    pdb.set_trace()
    for cr in range(nr):
#    for cr in range(1):
#        print('cr='+str(cr))
        cp=argmax(abs(mat[cr:,cr]))
        cp=cp+cr
#        print('cp='+str(cp))
        if cr!=cp:
            mat=RowSwap(mat,cr,cp)
#        print('mat='+str(mat))
        curpiv=abs(mat[cr,cr])
#        print('curpiv='+str(curpiv))
        if curpiv>eps:
            mat[cr:,cr:]=mat[cr:,cr:]/mat[cr,cr]
#            print('mat='+str(mat))
            for ro in range(nr-(cr+1)):
                ci=cr+ro+1
                mat[ci,:]=mat[ci,:]-mat[cr,:]*mat[ci,cr]
#                print('mat='+str(mat))
    return mat

def backsub(mat, eps=1e-14):
#    Pdb().set_trace()
    nr=shape(mat)[0]
    nc=shape(mat)[1]
    if nc==(nr+1):
        B=mat[:,nc-1]
        A=mat[:,0:nc-1]
    else:
        B=zeros((nr,1),'D')
        A=mat
    nca=shape(A)[1]
    vectout=zeros((nr,1),'D')
#    vectout=vectout*1j
    for cr in reverse(range(nr)):
        if abs(mat[cr,cr])<eps:
            vectout[cr]=1
        else:
            if cr<(nca-1):
                dp=dot(A[cr,cr+1:],vectout[cr+1:])[0]
            else:
                dp=0
            vectout[cr]=B[cr]-dp/A[cr,cr]
    return vectout

def rrefnull(matin, eps=1e-14):
#    pdb.set_trace()
    ut=uppertri(matin,eps=eps)
    ns=backsub(ut,eps=eps)
    return ns
