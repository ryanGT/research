from scipy import *


def thresh(iterin, value, startind=0, above=1):
    if above:
        mybools = iterin >= value
    else:
        mybools = iterin <= value
    myinds = arange(len(iterin))
    keepinds = myinds[squeeze(mybools)]
    inds2 = keepinds >= startind
    keep2 = keepinds[inds2]

    if len(keep2)==0:
        return -1
    else:
        return keep2.min()


###########################
#
#  Variable Down Sampling of data for rwkbode.compress
#
###########################

def FindLogInds(f, minf, maxf, N=50):
    logf = logspace(log10(minf), log10(maxf), N)
    logfinds = [thresh(f, item) for item in logf]
    return logfinds


def FindVarLogInds(f, minf, maxf, N):
    """N is points per decade."""
    nd = log10(maxf)-log10(minf)
    NN = int(nd*N+0.5)
    return FindLogInds(f, minf, maxf, NN)


def BuildMask(inds, vect):
    mymask = [item in inds for item in range(len(vect))]
    return mymask

###########################
